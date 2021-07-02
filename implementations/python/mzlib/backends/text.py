import re
import os
import io
import logging
import warnings
import enum

from typing import List, Tuple, Union, Iterable

from mzlib.index import MemoryIndex
from mzlib.annotation import parse_annotation
from mzlib.spectrum import Spectrum
from mzlib.attributes import AttributeManager, Attributed
from mzlib.analyte import ANALYTE_MIXTURE_TERM, Analyte, Interpretation, InterpretationMember

from .base import (
    SpectralLibraryBackendBase,
    _PlainTextSpectralLibraryBackendBase,
    SpectralLibraryWriterBase,
    FORMAT_VERSION_TERM)
from .utils import try_cast

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


term_pattern = re.compile(
    r"^(?P<term>(?P<term_accession>\S+:(?:\d|X)+)\|(?P<term_name>[^=]+))")
key_value_term_pattern = re.compile(
    r"^(?P<term>(?P<term_accession>[A-Za-z0-9:.]+:(?:\d|X)+)\|(?P<term_name>[^=]+))=(?P<value>.+)")
grouped_key_value_term_pattern = re.compile(
    r"^\[(?P<group_id>\d+)\](?P<term>(?P<term_accession>\S+:(?:\d|X)+)\|(?P<term_name>[^=]+))=(?P<value>.+)")
float_number = re.compile(
    r"^\d+(.\d+)?")


class SpectrumParserStateEnum(enum.Enum):
    unknown = 0
    header = 1
    analyte = 2
    interpretation = 3
    interpretation_member = 4
    peaks = 5
    done = 6


START_OF_SPECTRUM_MARKER = re.compile(r"^<Spectrum>")
START_OF_INTERPRETATION_MARKER = re.compile(r"^<Interpretation(?:=(.+))>")
START_OF_ANALYTE_MARKER = re.compile(r"^<Analyte(?:=(.+))>")
START_OF_PEAKS_MARKER = re.compile(r"^<Peaks>")
START_OF_LIBRARY_MARKER = re.compile(r"^<mzSpecLib\s+(.+)>")
SPECTRUM_NAME_PRESENT = re.compile(r'MS:1003061\|spectrum name=')
START_OF_INTERPRETATION_MEMBER_MARKER = re.compile(r"<InterpretationMember(?:=(.+))>")


class TextSpectralLibrary(_PlainTextSpectralLibraryBackendBase):
    file_format = "mzlb.txt"
    format_name = "text"

    @classmethod
    def guess_from_header(cls, filename: str) -> bool:
        with open(filename, 'r') as stream:
            first_line = stream.readline()
            if START_OF_SPECTRUM_MARKER.match(first_line) or START_OF_LIBRARY_MARKER.match(first_line):
                return True
        return False

    def _parse_header_from_stream(self, stream: io.IOBase) -> Tuple[bool, int]:
        nbytes = 0
        first_line = stream.readline()
        nbytes += len(first_line)
        if SPECTRUM_NAME_PRESENT.match(first_line) or START_OF_SPECTRUM_MARKER.match(first_line):
            return True, 0
        elif START_OF_LIBRARY_MARKER.match(first_line):
            match = START_OF_LIBRARY_MARKER.match(first_line)
            version = match.group(1)
            attributes = AttributeManager()
            attributes.add_attribute(FORMAT_VERSION_TERM, version)
            line = stream.readline()
            while not (SPECTRUM_NAME_PRESENT.match(line) or START_OF_SPECTRUM_MARKER.match(line)):
                nbytes += len(line)
                match = key_value_term_pattern.match(line)
                if match is not None:
                    d = match.groupdict()
                    attributes.add_attribute(
                        d['term'], try_cast(d['value']))
                    line = stream.readline()
                    nbytes += len(line)
                    continue
                if line.startswith("["):
                    match = grouped_key_value_term_pattern.match(line)
                    if match is not None:
                        d = match.groupdict()
                        attributes.add_attribute(
                            d['term'], try_cast(d['value']), d['group_id'])
                        attributes.group_counter = int(d['group_id'])
                        line = stream.readline()
                        nbytes += len(line)
                        continue
                    else:
                        raise ValueError(
                            f"Malformed grouped attribute {line}")
                elif "=" in line:
                    name, value = line.split("=")
                    attributes.add_attribute(name, value)
                else:
                    raise ValueError(f"Malformed attribute line {line}")
                line = stream.readline()
            self.attributes.clear()
            self.attributes._from_iterable(attributes)
            return True, nbytes
        return False, 0

    def read_header(self) -> Tuple[bool, int]:
        with open(self.filename, 'rt') as stream:
            return self._parse_header_from_stream(stream)

    def create_index(self) -> int:
        """
        Populate the spectrum index

        Returns
        -------
        n_spectra: int
            The number of entries read
        """

        #### Check that the spectrum library filename isvalid
        filename = self.filename

        #### Determine the filesize
        file_size = os.path.getsize(filename)

        with open(filename, 'r') as infile:
            state = 'header'
            spectrum_buffer = []
            n_spectra = 0
            start_index = 0
            file_offset = 0
            line_beginning_file_offset = 0
            spectrum_file_offset = 0
            spectrum_name = ''

            # Required for counting file_offset manually (LF vs CRLF)
            infile.readline()
            file_offset_line_ending = len(infile.newlines) - 1
            infile.seek(0)

            logger.info(f"Reading {filename} ({file_size} bytes)...")
            while 1:
                line = infile.readline()
                if len(line) == 0:
                    break

                line_beginning_file_offset = file_offset

                #### tell() is twice as slow as counting it myself
                file_offset += len(line) + file_offset_line_ending

                line = line.rstrip()
                if state == 'header':
                    # if re.match(r'MS:1003061\|spectrum name=', line):
                    if START_OF_SPECTRUM_MARKER.match(line):
                        state = 'body'
                        spectrum_file_offset = line_beginning_file_offset
                    else:
                        continue
                if state == 'body':
                    if len(line) == 0:
                        continue
                    # if re.match(r'MS:1003061\|spectrum name=', line):
                    if START_OF_SPECTRUM_MARKER.match(line):
                        if len(spectrum_buffer) > 0:
                            if not spectrum_name:
                                raise ValueError("No spectrum name")
                            self.index.add(
                                number=n_spectra + start_index,
                                offset=spectrum_file_offset,
                                name=spectrum_name,
                                analyte=None)
                            n_spectra += 1
                            spectrum_buffer = []
                            #### Commit every now and then
                            if n_spectra % 1000 == 0:
                                self.index.commit()
                                percent_done = int(
                                    file_offset/file_size*100+0.5)
                                logger.info(str(percent_done)+"%...")

                        spectrum_file_offset = line_beginning_file_offset
                        spectrum_name = ''
                    if re.match(r'MS:1003061\|spectrum name', line):
                        spectrum_name = re.match(r'MS:1003061\|spectrum name=(.+)', line).group(1)

                    spectrum_buffer.append(line)


            if not spectrum_name:
                raise ValueError("No spectrum name")
            self.index.add(
                number=n_spectra + start_index,
                offset=spectrum_file_offset,
                name=spectrum_name,
                analyte=None)
            self.index.commit()
            n_spectra += 1

            #### Flush the index
            self.index.commit()

        return n_spectra

    def _buffer_from_stream(self, infile: io.IOBase) -> List:
        state = 'body'
        spectrum_buffer = []

        for line in infile:
            line = line.rstrip()
            if state == 'body':
                if len(line) == 0:
                    continue
                if START_OF_SPECTRUM_MARKER.match(line):
                    if len(spectrum_buffer) > 0:
                        return spectrum_buffer
                spectrum_buffer.append(line)
        return spectrum_buffer

    def _parse_attribute_into(self, line: str, store: Attributed, line_number_message=lambda:'') -> bool:
        match = key_value_term_pattern.match(line)
        if match is not None:
            d = match.groupdict()
            store.add_attribute(d['term'], try_cast(d['value']))
            return True
        if line.startswith("["):
            match = grouped_key_value_term_pattern.match(line)
            if match is not None:
                d = match.groupdict()
                store.add_attribute(
                    d['term'], try_cast(d['value']), d['group_id'])
                store.group_counter = int(d['group_id'])
                return True
            else:
                raise ValueError(f"Malformed grouped attribute {line}{line_number_message()}")
        elif "=" in line:
            name, value = line.split("=")
            store.add_attribute(name, try_cast(value))
            return True
        else:
            raise ValueError(f"Malformed attribute line {line}{line_number_message()}")

    def _parse(self, buffer: Iterable, spectrum_index: int = None,
               start_line_number: int=None) -> Spectrum:
        spec: Spectrum = self._new_spectrum()
        interpretation: Interpretation = None
        analyte: Analyte = None
        interpretation_member: InterpretationMember = None

        STATES = SpectrumParserStateEnum
        state = STATES.header

        peak_list = []
        line_number = -1

        def real_line_number_or_nothing():
            nonlocal start_line_number
            nonlocal line_number
            nonlocal spectrum_index

            if start_line_number is None:
                return ''
            message = f" on line {line_number + start_line_number}"
            if spectrum_index is not None:
                message += f" in spectrum {spectrum_index}"
            message += f" in state {state}"
            return message

        for line_number, line in enumerate(buffer):
            line = line.strip()
            if not line:
                break
            if state == STATES.header:
                if START_OF_SPECTRUM_MARKER.match(line):
                    continue

                elif START_OF_PEAKS_MARKER.match(line):
                    state = STATES.peaks
                    continue

                elif START_OF_INTERPRETATION_MARKER.match(line):
                    state = STATES.interpretation
                    match = START_OF_INTERPRETATION_MARKER.match(line)
                    if interpretation is not None:
                        spec.add_interpretation(interpretation)
                    interpretation = self._new_interpretation(match.group(1))
                    spec.add_interpretation(interpretation)
                    analyte = None
                    continue

                elif START_OF_ANALYTE_MARKER.match(line):
                    state = STATES.analyte
                    match = START_OF_ANALYTE_MARKER.match(line)
                    analyte = self._new_analyte(match.group(1))
                    spec.add_analyte(analyte)
                    continue
                self._parse_attribute_into(line, spec, real_line_number_or_nothing)

            elif state == STATES.interpretation:
                if START_OF_ANALYTE_MARKER.match(line):
                    warnings.warn(
                        f"An analyte found after an interpretation was encountered, {real_line_number_or_nothing()}")
                    state = STATES.analyte
                    match = START_OF_ANALYTE_MARKER.match(line)
                    if analyte is not None:
                        spec.add_analyte(analyte)
                    analyte = self._new_analyte(match.group(1))
                    spec.add_analyte(analyte)
                    continue
                elif START_OF_INTERPRETATION_MARKER.match(line):
                    state = STATES.interpretation
                    match = START_OF_INTERPRETATION_MARKER.match(line)
                    if interpretation is not None:
                        spec.add_interpretation(interpretation)
                    interpretation = self._new_interpretation(match.group(1))
                    spec.add_interpretation(interpretation)
                    analyte = None
                    continue
                elif START_OF_PEAKS_MARKER.match(line):
                    state = STATES.peaks
                    continue
                elif START_OF_INTERPRETATION_MEMBER_MARKER.match(line):
                    state = STATES.interpretation_member
                    match = START_OF_INTERPRETATION_MEMBER_MARKER.match(line)

                    if interpretation_member is not None:
                        interpretation.add_member_interpretation(interpretation_member)

                    interpretation_member = InterpretationMember(match.group(1))
                    interpretation.add_member_interpretation(interpretation_member)
                    continue

                self._parse_attribute_into(line, interpretation.attributes, real_line_number_or_nothing)
                self._analyte_interpretation_link(spec, interpretation)

            elif state == STATES.interpretation_member:
                if START_OF_PEAKS_MARKER.match(line):
                    state = STATES.peaks
                    interpretation_member = None
                    interpretation = None
                    continue
                elif START_OF_INTERPRETATION_MARKER.match(line):
                    state = STATES.interpretation
                    match = START_OF_INTERPRETATION_MARKER.match(line)
                    if interpretation is not None:
                        spec.add_interpretation(interpretation)
                    interpretation = self._new_interpretation(match.group(1))
                    spec.add_interpretation(interpretation)
                    interpretation_member = None
                    continue
                elif START_OF_INTERPRETATION_MEMBER_MARKER.match(line):
                    state = STATES.interpretation_member
                    match = START_OF_INTERPRETATION_MEMBER_MARKER.match(line)
                    if interpretation_member is not None:
                        interpretation.add_member_interpretation(interpretation_member)
                    interpretation_member = InterpretationMember(match.group(1))
                    interpretation.add_member_interpretation(interpretation_member)
                    continue

                self._parse_attribute_into(
                    line, interpretation_member, real_line_number_or_nothing)

            elif state == STATES.analyte:
                if START_OF_PEAKS_MARKER.match(line):
                    state = STATES.peaks
                    if analyte is not None:
                        spec.add_analyte(analyte)
                        analyte = None
                    continue

                elif START_OF_ANALYTE_MARKER.match(line):
                    state = STATES.analyte
                    match = START_OF_ANALYTE_MARKER.match(line)
                    if analyte is not None:
                        spec.add_analyte(analyte)
                    analyte = self._new_analyte(match.group(1))
                    continue

                elif START_OF_INTERPRETATION_MARKER.match(line):
                    state = STATES.interpretation
                    match = START_OF_INTERPRETATION_MARKER.match(line)
                    if analyte is not None:
                        spec.add_analyte(analyte)
                        analyte = None

                    # Somehow we have an in-progress Interpretation that hasn't been cleared yet.
                    # This should probably be an error strictly speaking.
                    if interpretation is not None:
                        warnings.warn(
                            f"Interleaved analytes and interpretations detected at {real_line_number_or_nothing()}")
                        spec.add_interpretation(interpretation)
                    interpretation = self._new_interpretation(match.group(1))
                    spec.add_interpretation(interpretation)
                    continue

                self._parse_attribute_into(line, analyte, real_line_number_or_nothing)

            elif state == STATES.peaks:
                match = float_number.match(line)
                if match is not None:
                    tokens = line.split("\t")
                    n_tokens = len(tokens)
                    if n_tokens == 3:
                        mz, intensity, annotation = tokens
                        annotation = parse_annotation(annotation)
                        peak_list.append([float(mz), float(intensity), annotation, ""])
                    elif n_tokens == 4:
                        mz, intensity, annotation, aggregation = tokens
                        annotation = parse_annotation(annotation)
                        peak_list.append(
                            [float(mz), float(intensity), annotation, aggregation])
                    else:
                        raise ValueError(
                            f"Malformed peak line {line} with {n_tokens} entries{real_line_number_or_nothing()}")
                else:
                    raise ValueError(f"Malformed peak line {line}{real_line_number_or_nothing()}")
            else:
                raise ValueError(f"Unknown state {state}{real_line_number_or_nothing()}")
        spec.peak_list = peak_list
        # Backfill analytes into interpretations that never explicitly listed them.
        self._default_interpretation_to_analytes(spec)
        return spec

    def get_spectrum(self, spectrum_number: int=None, spectrum_name: str=None) -> Spectrum:
        # keep the two branches separate for the possibility that this is not possible with all
        # index schemes.
        if spectrum_number is not None:
            if spectrum_name is not None:
                raise ValueError(
                    "Provide only one of spectrum_number or spectrum_name")
            offset = self.index.offset_for(spectrum_number)
        elif spectrum_name is not None:
            offset = self.index.offset_for(spectrum_name)
        buffer = self._get_lines_for(offset)
        spectrum = self._parse(buffer, spectrum_number)
        return spectrum


class TextSpectralLibraryWriter(SpectralLibraryWriterBase):
    file_format = "mzlb.txt"
    format_name = "text"
    default_version = '1.0'

    def __init__(self, filename, version=None, **kwargs):
        super(TextSpectralLibraryWriter, self).__init__(filename)
        self.version = version
        self._coerce_handle(self.filename)

    def _write_attributes(self, attributes: Attributed):
        for attribute in attributes:
            if len(attribute) == 2:
                self.handle.write(f"{attribute[0]}={attribute[1]}\n")
            elif len(attribute) == 3:
                self.handle.write(
                    f"[{attribute[2]}]{attribute[0]}={attribute[1]}\n")
            else:
                raise ValueError(
                    f"Attribute has wrong number of elements: {attribute}")

    def write_header(self, library: SpectralLibraryBackendBase):
        if self.version is None:
            version = library.attributes.get_by_name("format version")
            if version is None:
                version = self.default_version
        else:
            version = self.version
        self.handle.write("<mzSpecLib %s>\n" % (version, ))
        self._write_attributes(library.attributes)

    def write_spectrum(self, spectrum: Spectrum):
        self.handle.write("<Spectrum>\n")
        self._write_attributes(spectrum.attributes)
        for analyte in spectrum.analytes.values():
            self.handle.write(f"<Analyte={analyte.id}>\n")
            self._write_attributes(analyte.attributes)
        for interpretation in spectrum.interpretations.values():
            interpretation: Interpretation

            if len(spectrum.analytes) == 1:
                attribs_of = self._filter_attributes(interpretation, self._not_analyte_mixture_term)
            else:
                attribs_of = interpretation.attributes
            self.handle.write(f"<Interpretation={interpretation.id}>\n")
            self._write_attributes(attribs_of)

            for member in interpretation.member_interpretations.values():
                member: InterpretationMember
                self.handle.write(f"<InterpretationMember={member.id}>\n")
                self._write_attributes(member.attributes)
        self.handle.write("<Peaks>\n")
        for peak in spectrum.peak_list:
            peak_parts = [
                str(peak[0]),
                str(peak[1]),
                '?' if not peak[2] else ",".join(map(str, peak[2]))
            ]
            if peak[3]:
                peak_parts.append(str(peak[3]))
            self.handle.write("\t".join(peak_parts)+"\n")
        self.handle.write("\n")

    def close(self):
        self.handle.close()


def format_spectrum(spectrum: Spectrum, **kwargs) -> str:
    buffer = io.StringIO()
    writer = TextSpectralLibraryWriter(buffer, **kwargs)
    writer.write_spectrum(spectrum)
    return buffer.getvalue()


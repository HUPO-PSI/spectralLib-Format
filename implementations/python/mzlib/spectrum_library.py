#!/usr/bin/env python3
from __future__ import print_function
import sys
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

import re
import timeit
import os
import pathlib

from typing import Type, List, Union

from mzlib.spectrum_library_index import SpectrumLibraryIndex
from mzlib.spectrum import Spectrum
from mzlib.index import MemoryIndex, SQLIndex, IndexBase
from mzlib.backends import guess_implementation, SpectralLibraryBackendBase, SpectralLibraryWriterBase


debug = False

class SpectrumLibrary:
    """
    SpectrumLibrary - Class for a spectrum library

    Attributes
    ----------
    identifier: str
        A unique identifier string assigned to this library by a spectral library host or
        provider.
    filename: str or Path
        A location on the local file system where the spectral library is stored
    format : str
        The name of the format for the current encoding of the library.
    backend: :class:`~.SpectralLibraryBackendBase`
        The implementation used to parse the file

    Methods
    -------
    read_header - Read just the header of the whole library
    read - Read the entire library into memory
    write - Write the library to disk
    create_index - Create an index file for this library
    transform - Not quite sure what this is supposed to be
    get_spectrum - Extract a single spectrum by identifier
    find_spectra - Return a list of spectra given query constraints

    """

    backend: SpectralLibraryBackendBase
    filename: str
    identifier: str
    format: str
    index_type: Type[IndexBase]


    def __init__(self, identifier=None, filename=None, format=None, index_type=None):
        """
        __init__ - SpectrumLibrary constructor

        Parameters
        ----------
        format : string
            Name of the format for the current encoding of the library.

        """
        self.backend = None
        self.identifier = identifier
        self.index_type = index_type
        self._format = format
        self.filename = filename

    def _init_from_filename(self, filename: str, index_type: Type[IndexBase]=None):
        if index_type is None:
            index_type = self.index_type
        if self.format is None:
            self.backend = guess_implementation(self.filename, index_type)
            self._format = self.backend.format_name
        else:
            backend_type = SpectralLibraryBackendBase.type_for_format(self.format)
            if backend_type is None:
                raise ValueError(
                    f"Could not find an implementation for {self.format}")
            self.backend = backend_type(
                self.filename, index_type=index_type)
            self._format = self.backend.format_name

    def _backend_initialized(self):
        return self.backend is not None

    def _requires_backend(self):
        if not self._backend_initialized():
            raise ValueError(
                "Cannot read library data, library parser not yet initialized")

    #### Define getter/setter for attribute identifier
    @property
    def identifier(self):
        return(self._identifier)

    @identifier.setter
    def identifier(self, identifier):
        self._identifier = identifier

    #### Define getter/setter for attribute filename
    @property
    def filename(self):
        return(self._filename)

    @filename.setter
    def filename(self, filename):
        self._filename = filename
        if filename is not None:
            self._init_from_filename(filename)

    #### Define getter/setter for attribute format
    @property
    def format(self):
        return self._format

    @property
    def index(self):
        if self._backend_initialized():
            return self.backend.index
        return None

    @property
    def attributes(self):
        if self._backend_initialized():
            return self.backend.attributes
        return None

    def read_header(self) -> bool:
        """Read just the header of the whole library

        Returns
        -------
        bool:
            Whether the operation was successful
        """
        self._requires_backend()
        return self.backend.read_header()

    def read(self):
        self._requires_backend()
        return self.backend.read()

    def write(self, destination, format: str=None):
        """Write the library to disk
        """
        filename = destination
        if not isinstance(filename, (str, pathlib.Path)):
            filename = getattr(destination, "name", None)

        if format is None and filename is not None:
            basename = os.path.basename(filename)
            tokens = basename.rsplit(".", 2)
            if len(tokens) == 3:
                writer_type = SpectralLibraryWriterBase.type_for_format('.'.join(tokens[1:]))
            else:
                raise ValueError(f"Could not guess format from file name {filename}")
        else:
            writer_type = SpectralLibraryWriterBase.type_for_format(format)
        if writer_type is None:
            raise ValueError(
                f"Could not find a format writer from file name {filename} or format {format}")
        writer = writer_type(destination)
        if self._backend_initialized():
            with writer:
                writer.write_library(self.backend)
        else:
            print("Library not initialized")
            writer.close()

    def get_spectrum(self, spectrum_number: int=None, spectrum_name: str=None) -> Spectrum:
        """Retrieve a single spectrum from the library.

        Parameters
        ----------
        spectrum_number : int, optional
            The index of the specturm in the library
        spectrum_name : str, optional
            The name of the spectrum in the library

        Returns
        -------
        :class:`~.Spectrum`
        """
        self._requires_backend()
        return self.backend.get_spectrum(spectrum_number, spectrum_name)

    def find_spectra(self, specification, **query_keys) -> List[Spectrum]:
        """
        find_spectra - Return a list of spectra given query constraints
        """
        self._requires_backend()
        return self.backend.find_spectra(specification, **query_keys)

    def __getitem__(self, i: int) -> Spectrum:
        self._requires_backend()
        return self.backend[i]

    def __len__(self):
        if self._backend_initialized():
            return len(self.backend)
        return 0

    def __iter__(self):
        if self._backend_initialized():
            return iter(self.backend)
        return iter([])

    def add_attribute(self, key, value, group_identifier=None):
        """Add an attribute to the library level attributes store.

        Parameters
        ----------
        key : str
            The name of the attribute to add
        value : object
            The value of the attribute to add
        group_identifier : str, optional
            The attribute group identifier to use, if any. If not provided,
            no group is assumed.
        """
        self._requires_backend()
        return self.backend.add_attribute(key, value, group_identifier=group_identifier)

    def get_attribute(self, key, group_identifier=None):
        """Get the value or values associated with a given
        attribute key from the library level attribute store.

        Parameters
        ----------
        key : str
            The name of the attribute to retrieve
        group_identifier : str, optional
            The specific group identifier to return from.

        Returns
        -------
        attribute_value: object or list[object]
            Returns single or multiple values for the requested attribute.
        """
        self._requires_backend()
        return self.backend.get_attribute(key, group_identifier=group_identifier)

    def remove_attribute(self, key, group_identifier=None):
        """Remove the value or values associated with a given
        attribute key from the library level attribute store.

        This rebuilds the entire store, which may be expensive.

        Parameters
        ----------
        key : str
            The name of the attribute to retrieve
        group_identifier : str, optional
            The specific group identifier to return from.

        """
        self._requires_backend()
        return self.backend.remove_attribute(key, group_identifier=group_identifier)

    def has_attribute(self, key):
        """Test for the presence of a given attribute in the library
        level store.

        Parameters
        ----------
        key : str
            The attribute to test for

        Returns
        -------
        bool
        """
        self._requires_backend()
        return self.backend.has_attribute(key)



#### Example using this class
def example():

    #### Create a new RTXFeedback object
    spectrum_library = SpectrumLibrary()
    spectrum_library.filename = "../refData/sigmaups1_consensus_final_true_lib.msp"
    spectrum_library.read_header()

    spectrum_buffer = spectrum_library.get_spectrum(spectrum_index_number=2000)
    spectrum = Spectrum()
    spectrum.parse(spectrum_buffer)
    buffer = spectrum.write(format="text")
    print(buffer)
    print()

    return()


#### If this class is run from the command line, perform a short little test to see if it is working correctly
def main():

    #### Run an example
    example()
    return()


if __name__ == "__main__": main()


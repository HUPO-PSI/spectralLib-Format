from mzlib.spectrum import Spectrum
import os
import unittest
import tempfile

from mzlib.backends import MSPSpectralLibrary, TextSpectralLibrary, JSONSpectralLibrary
from mzlib.analyte import ANALYTE_MIXTURE_TERM

from .common import datafile


class LibraryBehaviorBase(object):

    def _open_library(self, file_path=None):
        if file_path is None:
            file_path = self.test_file
        lib = self.library_cls(file_path)
        return lib

    def test_sequence_behavior(self):
        lib = self._open_library()
        assert len(lib) == 7
        spec = lib[3]
        assert spec.get_attribute(
            "MS:1003061|spectrum name") == "AAAAGSTSVKPIFSR/2_0_44eV"

    # TODO: Fix clipping in _buffer_from_stream first
    # def test_iteration(self):
    #     lib = self._open_library()
    #     assert list(lib) == list(lib.read())


class TestMSPLibrary(unittest.TestCase, LibraryBehaviorBase):
    test_file = datafile("chinese_hamster_hcd_selected_head.msp")
    library_cls = MSPSpectralLibrary


class MzSpecLibLibraryBehaviorBase(LibraryBehaviorBase):
    def test_interpretation_attribute(self):
        try:
            path = self.test_interpretation_file
        except AttributeError:
            return
        lib = self._open_library(path)
        spec: Spectrum = lib[0]
        interp = spec.interpretations[1]
        assert len(interp.attributes) == 2
        assert interp.get_attribute('MS:1002357|PSM-level probability') == 0.974
        assert interp.get_attribute(ANALYTE_MIXTURE_TERM) == '1,2'
        mem = interp.get_member_interpretation(1)
        assert len(mem) == 1
        assert mem.get_attribute('MS:XXXXXXX|explained intensity fraction') == 0.287


class TestTextLibrary(unittest.TestCase, MzSpecLibLibraryBehaviorBase):
    test_file = datafile("chinese_hamster_hcd_selected_head.mzlb.txt")
    library_cls = TextSpectralLibrary
    test_interpretation_file = datafile("complex_interpretations_with_members.mzlb.txt")


class TestJSONLibrary(unittest.TestCase, MzSpecLibLibraryBehaviorBase):
    test_file = datafile("chinese_hamster_hcd_selected_head.mzlb.json")
    library_cls = JSONSpectralLibrary
    test_interpretation_file = datafile("complex_interpretations_with_members.mzlb.json")

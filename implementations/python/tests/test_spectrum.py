import os
import unittest
import tempfile

from mzlib.backends import MSPSpectralLibrary, TextSpectralLibrary

from .common import datafile


class TestSpectrum(unittest.TestCase):

    def get_library(self):
        test_file = datafile("chinese_hamster_hcd_selected_head.mzlb.txt")
        return TextSpectralLibrary(test_file)

    def get_spectrum(self, index):
        library = self.get_library()
        return library.get_spectrum(index)

    def test_write(self):
        spectrum = self.get_spectrum(0)
        buffer = spectrum.write('text')
        lines = buffer.splitlines()
        n_lines = len(lines)
        assert n_lines == 133
        assert buffer.startswith("<Spectrum>\nMS:1003061|spectrum name")

    def test_equality(self):
        spectrum = self.get_spectrum(0)
        spectrum2 = self.get_spectrum(0)
        assert spectrum == spectrum2


'''
test_globals.py - test module for globals.py in PyAR
'''
'''
Copyright (C) 2016 by Surajit Nandi, Anoop Ayyappan, and Mark P. Waller
Indian Institute of Technology Kharagpur, India and Westfaelische Wilhelms
Universitaet Muenster, Germany

This file is part of the PyAR project.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
'''

'''This test is for checking the stability of global variables inside
globals module'''
import sys
sys.dont_write_bytecode = True
import unittest
from PyAR.utils import globals

class value_check_test(unittest.TestCase):
    def test_autoang(self):
        value_get = globals.AUTOANG
        actual_value = 0.52917726
        self.assertAlmostEqual(value_get,actual_value)

    def test_autokj(self):
        value_get = globals.AUTOKJ
        actual_value = 2625.5
        self.assertAlmostEqual(value_get, actual_value)

    def test_autokcal(self):
        value_get = globals.AUTOKCAL
        actual_value = 627.509541
        self.assertAlmostEqual(value_get, actual_value)

    def test_elements_copper(self):
        all_elements = globals.ELEMENTS
        element_covalent_radii = all_elements['cu']
        actual_covalent_radii = 1.01
        self.assertAlmostEqual(element_covalent_radii, actual_covalent_radii)

    def test_elements_arsenic(self):
        all_elements = globals.ELEMENTS
        element_covalent_radii = all_elements['as']
        actual_covalent_radii = 1.15
        self.assertAlmostEqual(element_covalent_radii, actual_covalent_radii)


if __name__ == "__main__":
    unittest.main()

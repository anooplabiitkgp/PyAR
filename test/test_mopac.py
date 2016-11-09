'''
test_mopac.py - test module for mopac wrapper in PyAR
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

# This is the test file for the mopac wrapper module.

import unittest
import os
from PyAR.wrapper import mopac


class TestMopacPrepareInput(unittest.TestCase):

    def setUp(self):
        self.h2_xyz_file = 'tmp.xyz'
        with open(self.h2_xyz_file, 'w') as ftemp:
            ftemp.write("   %d\n" % 2)
            ftemp.write("  \n")
            ftemp.write("H 0.0 0.0 0.0\n")
            ftemp.write("H 0.0 0.0 0.7\n")
        mopac_inst = mopac.Mopac(self.h2_xyz_file)
        mopac_inst.prepare_input(keyword="PM7")


    def test_method(self):
        with open('tmp.mop') as finp:
            lines = finp.readlines()
        method = " ".join(lines[0].split())
        self.assertEqual('PM7', method)

    def test_no_ofcoords_lines(self):
        with open('tmp.mop') as finp:
            lines = finp.readlines()
        no_of_coords_lines = len(lines[3:])
        self.assertEqual(2, no_of_coords_lines)

    def test_return_status(self):
        mopac_inst = mopac.Mopac(self.h2_xyz_file)
        status = mopac_inst.prepare_input(keyword="PM7")
        self.assertEqual(status, 0)

    def tearDown(self):
        all_files = os.listdir("./")
        for files in all_files:
            if os.path.splitext(files)[0] == 'tmp':
                os.remove(files)

class TestMopacEnergyOptimize(unittest.TestCase):
    """ This unittest is for testing the energy and optimize function
    """
    def setUp(self):
        h2_xyz_file = 'tmp.xyz'
        with open(h2_xyz_file, 'w') as ftemp:
            ftemp.write("   %d\n" % 2)
            ftemp.write("  \n")
            ftemp.write("H 0.0 0.0 0.0\n")
            ftemp.write("H 0.0 0.0 0.7\n")
            mopac_inst = mopac.Mopac(h2_xyz_file)
            mopac_inst.prepare_input(keyword="PM7")

    def test_optimize(self):
        status = self.mopac_inst.optimize()
        self.assertEqual(status, 0)

    def test_get_energy_kcal(self):
        energy = self.mopac_inst.get_energy()
        enenergy_kcal = energy[0]
        self.assertAlmostEqual(enenergy_kcal, -32.01058)

    def test_get_energy_kj(self):
        enegy = self.mopac_inst.get_energy()
        energy_kj = enegy[1]
        self.assertAlmostEqual(energy_kj, -133.93228)


if __name__=="__main__":
    unittest.main()

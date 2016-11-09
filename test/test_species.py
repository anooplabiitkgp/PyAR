'''
test_species.py - test module for species.py in PyAR
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

import unittest
from PyAR import species

__author__ = 'surajit'


class TestExtrema(unittest.TestCase):
    coordinates = [[-4.0, 2.5, 0.0], [1.0, -0.6, 3.1], [3.0, -1.9, -3.1]]
    xcheck = species.extrema(coordinates, 'X')
    ycheck = species.extrema(coordinates, 'Y')
    zcheck = species.extrema(coordinates, 'z')

    def test_get_single_dimension_x(self):
        expected_data = [-4.0, 1.0, 3.0]
        data = self.xcheck.get_single_dimension()
        self.assertEqual(expected_data, data)

    def test_get_single_dimension_y(self):
        expected_data = [2.5, -0.6, -1.9]
        data = self.ycheck.get_single_dimension()
        self.assertEqual(expected_data, data)

    def test_get_single_dimension_z(self):
        expected_data = [0.0, 3.1, -3.1]
        data = self.zcheck.get_single_dimension()
        self.assertEqual(expected_data, data)

    def test_get_extreme_points_xaxis(self):
        expected_maxima = 3.0
        expected_minima = -4.0
        all_x_values = self.xcheck.get_single_dimension()
        maxima, minima = self.xcheck.get_extreme_points(all_x_values)
        self.assertEqual(expected_maxima, maxima)
        self.assertEqual(expected_minima, minima)

    def test_get_extreme_points_yaxis(self):
        expected_maxima = 2.5
        expected_minima = -1.9
        all_y_values = self.ycheck.get_single_dimension()
        maxima, minima = self.xcheck.get_extreme_points(all_y_values)
        self.assertEqual(expected_maxima, maxima)
        self.assertEqual(expected_minima, minima)

    def test_get_extreme_points_zaxis(self):
        expected_maxima = 3.1
        expected_minima = -3.1
        all_z_values = self.zcheck.get_single_dimension()
        maxima, minima = self.zcheck.get_extreme_points(all_z_values)
        self.assertEqual(expected_maxima, maxima)
        self.assertEqual(expected_minima, minima)


class test_get_extrema(unittest.TestCase):
    coordinates = [[1.0, -1.5, 2.7], [3.0, 0.7, -1.2], [-4.0, 0.8, -1.5]]

    def test_get_extrema_xrange(self):
        expect_maxima = 3.0
        expect_minima = -4.0
        maxima, minima = species.get_extrema(self.coordinates, 'X')
        maxima, minima = Species.get_extrema(self.coordinates, 'X')
        self.assertEqual(expect_maxima, maxima)
        self.assertEqual(expect_minima, minima)

    def test_get_extrema_yrange(self):
        expected_maxima = 0.8
        expected_minima = -1.5
        maxima, minima = species.get_extrema(self.coordinates, 'Y')
        maxima, minima = Species.get_extrema(self.coordinates, 'Y')
        self.assertEqual(expected_maxima, maxima)
        self.assertEqual(expected_minima, minima)

    def test_get_extrema_zrange(self):
        expected_maxima = 2.7
        expected_minima = -1.5
        maxima, minima = species.get_extrema(self.coordinates, 'Z')
        maxima, minima = Species.get_extrema(self.coordinates, 'Z')
        self.assertEqual(expected_maxima, maxima)
        self.assertEqual(expected_minima, minima)

class test_other_species(unittest.TestCase):
    coordinates = [[1.0, 2.5, 0.0], [3.0, 3.1, 5.9], [0.5, 4.3, 4.6]]

    def test_centroid(self):
        """ Test if the method centroid returns true center of coordinates
        """
        expected_com = [1.5, 3.3, 3.5]
        com_calculated = species.centroid(self.coordinates)
        for i in range(3):
            self.assertAlmostEqual(expected_com[i], com_calculated[i])

    def test_set_origin(self):
        """ Test if the method set_origin is actually setting the centroid
            to the point[0,0,0].
        """
        displaced_coordinates = species.set_origin(self.coordinates[:])
        com = species.centroid(displaced_coordinates)
        for i in range(3):
            self.assertAlmostEqual(com[i], 0.0)


if __name__=="__main__":
    unittest.main()

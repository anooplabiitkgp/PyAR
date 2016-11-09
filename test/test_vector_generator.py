'''
test_vector_generator.py - test module for vector_generator.py in PyAR
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
GNU General Public License for more details.#
This is a test file for vector_generator
'''

import sys
sys.dont_write_bytecode = True
from PyAR.utils.vector_gen import vector_generator
import unittest


class test_basic_inputs(unittest.TestCase):
    def test_wrong_dimension_less_number_of_vector_component_tabu_vector(self):
        """unit test for giving 'details' list which contains less elements than the dimension.
        details is an empty list. To be right, it should be a list with contents 3 lists."""
        dummy = vector_generator.vector_generator(3)
        details = []
        self.assertRaises(vector_generator.DimensionMatchError, dummy.tabu_list_vector_generator, details)
    
    def test_wrong_dimension_more_number_of_vector_component_tabu_vector(self):
        """Test for giving 'details' list which contains more elements than the dimension.
        details is an list of empty lists. to be right, it should also contains a lists of 3 lists."""
        dummy = vector_generator.vector_generator(3)
        details = [[], [], [], []]
        self.assertRaises(vector_generator.DimensionMatchError, dummy.tabu_list_vector_generator, details)

    def test_wrong_dimension_more_number_of_vector_component_normal_vector(self):
        """ Test for methods for argument that is supplied with more number of vector components
        info than dimension for in case of normal random vector generator
        """
        dimension = 3
        dummy = vector_generator.vector_generator(dimension)
        inp_list = [[(1, 2)], [(3, 4)], [(5, 6)], [(7, 8)]]
        self.assertRaises(vector_generator.DimensionMatchError, dummy.normal_vector_generator, details=inp_list)

    def test_wrong_dimension_less_number_of_vector_component_normal_vector(self):
        """ Test for methods for argument that is supplied with less number of vector components
        info than dimension for in case of normal random vector generator
        """
        dimension = 3
        dummy = vector_generator.vector_generator(dimension)
        inp_list = [[(1, 2)], [(3, 4)]]
        self.assertRaises(vector_generator.DimensionMatchError, dummy.normal_vector_generator, details=inp_list)

    def test_number_of_list_element_in_tabu_vector_method(self):
        """Test for input elements less than 2 for tabu_list_vector_generator method.
        It is expected to contains at least two variable (range, tabulist) inside details."""
        dummy = vector_generator.vector_generator(3)
        inp_list = [[1], [2], [3]]
        self.assertRaises(vector_generator.TabuVectorArgInfoError, dummy.tabu_list_vector_generator, details=inp_list)

    def test_number_of_list_elem_in_normal_vector_method(self):
        """ Test for list length less than 2"""
        dummy = vector_generator.vector_generator(3)
        inp_list = [[1, 2, 3], [2], [3]]
        self.assertRaises(vector_generator.NormalVectorArgInfoError, dummy.normal_vector_generator, details=inp_list)


class integration_test(unittest.TestCase):

    def test_tabu_dimension_1(self):
        """Integration test for the random vector generation method with
        tabu search in consideration for 1d vector"""
        inp_list = []
        prev_tabulist = []
        vrange = (-15, 15)
        prev_vector = [-111111]
        for i in range(100):
            inp_list.extend([[prev_tabulist, vrange], ])
            dummy = vector_generator.vector_generator(1)
            vector, tabu_list = dummy.tabu_list_vector_generator(details=inp_list)
            self.assertNotEqual(vector[0], prev_vector[0])
            self.assertNotEqual(prev_tabulist, tabu_list)
            prev_vector = vector[:]
            for j in tabu_list:
                prev_tabulist = j[:]
            inp_list = []

    def test_tabu_vector_dimension_3(self):
        """Integration test for the random vector generation method with
        tabu list in consideration for 3d vector"""
        inp_list = []
        vrange = (-15, 15)
        prev_tabu_list = [[], [], []]
        prev_vector = [-11111, -111111, -111111]
        for i in range(10):
            inp_list.extend([[prev_tabu_list[0], vrange], [prev_tabu_list[1], vrange], [prev_tabu_list[2], vrange]])
            dummy = vector_generator.vector_generator(3)
            vector, new_tabu_list = dummy.tabu_list_vector_generator(details=inp_list[:])
            for j in range(3):
                self.assertNotEqual(vector[j], prev_vector[j])
                self.assertNotEqual(prev_tabu_list[j], new_tabu_list[j])
                prev_tabu_list[j] = new_tabu_list[j][:]
                prev_vector[j] = vector[j]
            inp_list = []

    def test_normal_vector_dimension_1(self):
        """As we cannot perform any functional test as it is purely
        random geometry generation method, we will perform only the
        dimensionality test"""
        dimension = 1
        for i in range(10):
            dummy = vector_generator.vector_generator(dimension)
            new_vector = dummy.normal_vector_generator(details=[[(-15, 15)]])
            self.assertTrue(len(new_vector) == dimension)

    def test_normal_vector_dimension_3(self):
        """Here, we will see if the function is returning the vector
        with same dimension or not without width mension"""
        dimension = 3
        trange = (-15, 15)
        for i in range(10):
            dummy = vector_generator.vector_generator(dimension)
            new_vector = dummy.normal_vector_generator(details=[[trange], [trange], [trange]])
            self.assertEqual(dimension, len(new_vector))

    def test_normal_vector_dimension_3_with_gridsize(self):
        """Here, we will see if the function is returning the vec
        with same dimension with width mension"""
        dimension = 3
        trange = (0, 180)
        width = 6
        for i in range(10):
            dummy = vector_generator.vector_generator(dimension)
            new_vector = dummy.normal_vector_generator(details=[[trange, width], [trange, width], [trange, width]])
            self.assertEqual(dimension, len(new_vector))


if __name__ == "__main__":
    unittest.main()

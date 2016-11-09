'''
test_vector.py - test module for vector.py in PyAR
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


import sys
sys.dont_write_bytecode = True
import unittest
from PyAR.utils.vector_gen import tabu


class test_tabu(unittest.TestCase):
    test_tabu_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    test_tabu_length = 5
    instance_for_test = tabu.tabu(test_tabu_list, test_tabu_length)

    def test_test_elements_presence(self):
        """ Test for test_for_presence method of tabu class. test will
        fail if return is false
        """
        element_to_be_checked = 4
        switch = self.instance_for_test.test_elements_presence(element_to_be_checked)
        self.assertTrue(switch)

    def test_test_elements_absent(self):
        """ This function will pass if the method test_elements_presence
        return False.
        """
        element_to_be_checked = 11
        switch = self.instance_for_test.test_elements_presence(element_to_be_checked)
        self.assertFalse(switch)

    def test_update_tabu_list(self):
        """ This test will pass if a element is added at the end of a list by
        update_tabu_list method. Inside this method, a new instance of the tabu
        class will be created. Otherwise, the update_tabu_list method will
        change the original list (test_list)
        """
        element_to_be_append = 11
        new_test_list = self.test_tabu_list[:]
        new_test_instance = tabu.tabu(new_test_list, self.test_tabu_length)
        new_test_instance.update_tabu_list(element_to_be_append)
        self.assertNotEqual(self.test_tabu_list, new_test_list)
        self.assertEqual(element_to_be_append, new_test_list[-1])

    def test_pop_tabu_list_from_start(self):
        """ This test will pass if the zeroth element from the list is poped. It will
        also create a new list and a new instance because, the method pop_tabu_list
        will change the original list.
        """
        new_test_list = self.test_tabu_list[:]
        new_test_instance = tabu.tabu(new_test_list, self.test_tabu_length)
        new_test_instance.pop_tabu_list_from_start()
        # to check if the list really changed
        self.assertNotEqual(new_test_list, self.test_tabu_list)
        # to check if the length really changed
        self.assertNotEqual(len(new_test_list), len(self.test_tabu_list))

    def test_get_old_position(self):
        """ This function is for checking if the position returning by the
        method get_old_position is correct or not.
        """
        element = 7
        position = self.instance_for_test.get_old_position(element)
        actual_position = self.test_tabu_list.index(element)
        self.assertEqual(position, actual_position)

class test_vector_component(unittest.TestCase):

    def test_vector_value(self):
        ''' Test for validity of the vector returned from the class.'''
        self.own_list = [12,24,36,48,60,72,84,96]
        A = vector_component_generator(self.own_list[:], (12,37), tabulength = 5, interval = 12)
        vector,mod_tabulist = A.generator()
        self.assertTrue(vector >= 12 and vector <= 48)

    def test_return_list_tuple(self):
        ''' This object will test if the return modified tabulist
        is a tuple or a list'''
        A = vector_component_generator(self.test_list[:], (3,4), tabulength = 5)
        vector, mod_tabulist = A.generator()
        self.assertTrue(mod_tabulist == list(mod_tabulist) )

    def test_tabu_work(self):
        ''' This object will test for the validity of the principle of
        tabulist method by a loop.'''
        self.own_list = []
        prev_vector = 20
        for i in range(100):
            A = vector_component_generator(self.own_list[:], (-15, 15), tabulength = 5)
            vector, mod_tabulist = A.generator()
            self.assertNotEqual(mod_tabulist, self.own_list)
            self.assertNotEqual(vector,prev_vector)
            prev_vector = vector
            self.own_list = mod_tabulist[:]

    def test_two_instance_diff(self):
        A = vector_component_generator([], (3,4), tabulength = 5)
        B = vector_component_generator([], (3,4), tabulength = 5)
        vectorA,listA = A.generator()
        vectorB,listB = B.generator()
        self.assertEqual(len(list(listA)), len(list(listB)))

if __name__ == "__main__":
    unittest.main()

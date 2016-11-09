'''
tabu.py - to test if a solution is tabu
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

import random
from PyAR.IO import files_dirs

def generate_random_number(initial_value, final_value, interval):
    random_number = random.randrange(initial_value, final_value, interval)
    return random_number


def test_elements_presence(new_element, tabu_file):
    tabu_status = None
    if files_dirs.file_exists_check_with_return_status(tabu_file):
        tabu_status = False
        with open(tabu_file) as fp:
            all_element = fp.readlines()
        for element in all_element:
            if new_element == element:
                tabu_status = True
                break
    else:
        tabu_status = False
    return tabu_status


def update_tabu_file(new_element_trans, new_element_rot, tabu_file):
    """This object will take an element and update it to the tabu list"""
    if file is None:
        raise IOError("No file supplied for tabu element to save.")
    else:
        with open(tabu_file, 'a') as fp:
            for new_elements in new_element_trans:
                fp.write("%5s " % new_elements)
            for new_elements in new_element_rot:
                fp.write("%5s" % new_elements)
            fp.write("\n")
            fp.close()
    return


def main():
    pass

if __name__ == "__main__":
    main()

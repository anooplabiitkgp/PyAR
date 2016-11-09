'''
fragment_make.py - to make 'fragment' file from 2 xyz files in PyAR
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

class Error(Exception):pass
class IncorrectXyzFileError(Error):pass


def _get_no_of_atoms(files):
    """This function will test the xyz file format.
       if it finds it alright, then it will return
       the number of atoms.
    """
    with open(files) as xyzfile:
        filecont = xyzfile.readlines()
    no_of_atoms = len(filecont) - 2
    if no_of_atoms != int(filecont[0]):
        raise IncorrectXyzFileError
    return no_of_atoms


def get_no_of_atoms_list(*files):
    """This function is for making a list of all the atom numbers
       in all the xyz files (fragments) provided
    """
    all_file_atoms = []
    for ifile in files:
        no_of_atom = _get_no_of_atoms(ifile)
        all_file_atoms.append(no_of_atom)
    return all_file_atoms


def get_fragment_string(all_xyz_file_atom_number):
    """This function will make the lines specified by the dftd3 code
       (modified by A. Anoop & M. Waller) for fragment file.
    """
    init = 1
    final = 0
    line = []
    for i in range(len(all_xyz_file_atom_number)):
        final = int(all_xyz_file_atom_number[i]) + final
        string = str(init) + "-" + str(final)
        init = final + 1
        line.append(string)
    return line


def frag_lines_mod(lines_list):
    """This function will check if for a line the number in both side of
       the 'hiphen' is same or not. If they same, then it will replace the
       line with one number.
    """
    return_line = []
    for lines in lines_list:
        start,end = map(int,lines.split('-'))
        if start == end:
            lines = start
            return_line.append(lines)
        else:
            return_line.append(lines)
    return return_line


def fragment_make(molecule1, molecule2, fragment_file):
    """This is the functional which by help of the above functions
       create the fragment file
    """
    no_of_atoms_list = get_no_of_atoms_list(molecule1, molecule2)
    lines_raw = get_fragment_string(no_of_atoms_list)
    all_lines_for_fragment = frag_lines_mod(lines_raw)
    with open(fragment_file,'w') as f:
        for i in all_lines_for_fragment:
            f.write("%s\n" % i)
    return


def main():
    pass

if __name__ == "__main__()":
    main()

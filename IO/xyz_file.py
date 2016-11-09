''' 
xyz_file.py - handling xyz file in PyAR
'''
'''
Copyright (C) 2016 by Surajit Nandi, Anoop Ayyappan, and Mark P. Waller
Indian Institute of Technology Kharagpur, India and Westfaelische Wilhelms
Universitaet Muenster, Germany

xyz_file.py is part of the PyAR project.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
'''

from PyAR.IO import files_dirs
def flatten(inputs, input_types=(list, tuple)):
    """This function is for flattening (topic flattening with python ).
       The main purpose of this function is to make a single list from
       a list of lists (got from the internet):
       >>> flatten([1,2,[3,4]])
       >>> [1,2,3,4]
       """
    input_types = type(inputs)
    inputs = list(inputs)
    i = 0
    while i < len(inputs):
        while isinstance(inputs[i], input_types):
            if not inputs[i]:
                inputs.pop(i)
                i -= 1
                break
            else:
                inputs[i:i+1] = inputs[i]
        i += 1
    return input_types(inputs)


def make_single_xyz_from_multi_xyz(outfile, *args):
    """This function will generate a single xyz file from two or more
       xyz files. There is a use of 'flatten' function in this funciton.
       This is to make a single list from a lists of coordinates
    """

    for arg in args:
        files_dirs.important_file_check(arg)
    number_of_atoms_list = []
    all_coordinates_list = []
    all_atoms_names_list = []
    total_atoms = 0
    for arg in args:
        number_of_atoms, title, atom_names, xyz_coords = read_xyz(arg)
        number_of_atoms_list.append(number_of_atoms)
        all_coordinates_list.append(xyz_coords)
        all_atoms_names_list.append(atom_names)
    all_coordinates = flatten(all_coordinates_list)
    all_atoms = flatten(all_atoms_names_list)
    for numbers in number_of_atoms_list:
        total_atoms += int(numbers)
    outfile_object = open(outfile, 'w')
    outfile_object.write("%d\n" % total_atoms)
    outfile_object.write("XYZ file from PyReactor\n")
    for n in range(len(all_coordinates)):
        outfile_object.write("%-5s %17.10f %17.10f %17.10f \n" % (all_atoms[n], all_coords[n][0], all_coords[n][1], all_coords[n][2]))
    outfile_object.close()
    return


def read_xyz(filename):
    """ This is a function for reading xyz files. It will return the
        number of atoms, title and xyz coordinates.
    """
    f = open(filename).readlines()
    nat = int(f[0])
    title = f[1].rstrip()
    raw = [f[i].split() for i in range(2,nat+2)]
    atom_names = [(raw[i][0]) for i in range(nat)]
    xyz_coords = [(float(raw[i][1]), float(raw[i][2]), float(raw[i][3])) for i in range(nat)]
    return nat, title, atom_names, xyz_coords


def main():
    pass


if __name__ == "__main__":
    main()

'''
structure_check.py - checking structural similarity using InCHI and Smiles 
in PyAR
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

import subprocess as subp
import os


def make_smile_file(outfilename='babelcheck.dat', xyzfile='trj'):
    """This function will generate a 'babelcheck.dat' file from the xyz
    file.
    """
    with open(outfilename, 'w') as f:
        subp.check_call(["babel", "-ixyz", str(xyzfile), "-osmi"], stdout=f)
    return


def make_smile_string_from_xyz(xyzfile):
    """This function will make smile string from a xyz file.
       if more than one xyz file provide, it will return smiles
       in a list. if one xyz file supplied, it will return the
       string
    """
    if os.path.isfile(xyzfile):
        pre_smile = subp.check_output(["babel", "-ixyz", str(xyzfile), "-osmi"])
        smile = pre_smile.split()[0]
        #print "smile string inside make_smile_string_from_xyz() is: ", smile
        return smile
    else:
        raise IOError("file %s does not exists" % xyzfile)


def make_inchi_string_from_xyz(xyzfile):
    """This function will make a inchi string from a xyz file with
       babel as the tools
    """
    if os.path.isfile(xyzfile):
        inchi = subp.check_output(["babel", "-ixyz", str(xyzfile), "-oinchi"])
        return inchi
    else:
        raise IOError("file %s does not exists" % xyzfile)


def compare_string(string1, string2):
    """This function will compare between two smiles string and return
       "same" if they are same and return "different" if they are
       different
    """
    status = None
    if string1 == string2:
        status = "same"
    elif string1 != string2:
        status = "different"
    else:
        raise ValueError("unknown value string1, string2")
    return status


def prev_product_comp_status(structure, path_to_dir, cheminfo="inchi"):
    """This function actually take a xyz file and a directory where
       other xyz files are present. It then return a status of the
       similarity of those structures. Or, in another sentence, it
       compares with other structures.
    """
    status = "different"
    if os.listdir(path_to_dir):
        for files in os.listdir(path_to_dir):
            #print "file inside pdtdir,", path_to_dir, "is", files, "inside prev_pdt_comp_status"
            path_to_files = os.path.join(path_to_dir, files)
            return_status = structure_comparison(structure, path_to_files, cheminfo=cheminfo)
            if return_status == "same":
                status = "same"
    else:
        status = "different"
    return status


def structure_comparison(structure1, structure2, cheminfo="inchi"):
    """This function will check the smile string between two structures.
       if the structures are same, then it will return a string "same" 
       and if different, it will return "different" 
    """
    status = None
    if cheminfo == "inchi":
        structure1_string = make_inchi_string_from_xyz(structure1)
        structure2_string = make_inchi_string_from_xyz(structure2)
        #print "structure1, structure2 inchis are:\n ", structure1_string, structure2_string, "\ninside structure_comparison"
    else:
        structure1_string = make_smile_string_from_xyz(structure1)
        structure2_string = make_smile_string_from_xyz(structure2)
        #print "structure1, structure2 smiles are:\n ", structure1_string, structure2_string, "\ninside structure_comparison"
    status = compare_string(structure1_string, structure2_string)
    return status

def check_smi_inchi_status(*args):
    string_status = "different"
    for argument in args:
        if argument == "same":
            string_status = "same"
    return string_status


def trj_to_xyz(pdtfile, trj_file, index=None):
    """This function will give a particular snap at step n(index).
    """
    with open(trj_file) as f:
        trj_cont = f.readlines()
    natom = int(trj_cont[0])
    if index is None:
        count = 0
        for line in trj_cont:
            if "Energy" in line:
                count += 1
        index = count
    req_content = trj_cont[((index-1)*(natom + 2)):(index * (natom + 2))]
    fout = open(pdtfile, 'w')
    fout.write("  %s" % req_content[0])
    fout.write(" %s" % req_content[1])
    for i in req_content[2:]:
        j = i[:3] + "    " + i[5:16] + "    " + i[18:29] + "    " + i[31:42]
        fout.write("%s\n" % j)
    fout.close()


def main():
    pass


if __name__ == "__main__":
    main()

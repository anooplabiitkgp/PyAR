''' 
analysis.py module for analysis of the results
'''
'''
PyAR
Copyright (C) 2016 by Surajit Nandi, Anoop Ayyappan, and Mark P. Waller
Indian Institute of Technology Kharagpur, India and Westfaelische Wilhelms
Universitaet Muenster, Germany

This file is part of the PyAR project.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.'''

import sys, os
import operator
from PyAR.IO import files_dirs
from PyAR import cluster
from PyAR.wrapper import mopac
import glob
import pprint

"""This program is for analysis of the PyAR calculations. If only one path is
   given without any switch, then, it will read all the arcfiles and print the 
   lowest energy structure in this geometry. 
   If many paths are given without any switch, it will return lowest energy 
   structure from each of the paths.
   if -e switch is given along with many paths, it will search for aggregate_* 
   directory in each path and then it will report each lowest energy for each 
   run for each aggregate.
   if -c is given, along with a path, it will calculate the cooperativity and 
   report the cooperativity of each of the molecules present in there in the
   form of a table."""


class cooperativity(object):
    def __init__(self, path, monomer):
        self.path = os.path.abspath(path)
        #monomer variable is the no. of atoms in a monomeric unit in the cluster.
        self.monomer = monomer
    def get_all_xyz_files(self):
        self.all_xyz_files = files_dirs.get_dirs_files(destination = self.path, wildcard = "result*.xyz")
        return self.all_xyz_files
    def make_dirs_for_all_xyz_files(self):
        current_dir = os.getcwd()
        self.coop_dirs = []
        os.chdir(self.path)
        for i in self.all_xyz_files:
            dir_name = "c_"+i[:-4]
            self.coop_dirs.append(dir_name)
            if not os.path.exists(dir_name):
                os.makedirs(dir_name)
        os.chdir(current_dir)
        return
    def cp_xyz_files_to_dirs(self):
        for i in self.all_xyz_files:
            files_dirs.scopy(os.path.join(self.path, i), os.path.join(self.path,"c_"+i[:-4]))
        return
    def read_xyz_c(self, xyz_file):
        nat, title, xyz = cluster.read_xyz(xyz_file)
        nmer = len(xyz)/self.monomer
        monomer_coords = []
        mer_count = 0
        for i in range(nmer):
            monomer_coords.append(xyz[mer_count:mer_count+self.monomer])
            mer_count += self.monomer
        return monomer_coords
    def make_monomers(self, mers_list):
        n_of_molecules = len(mers_list)
        energy_list = []
        for i in range(n_of_molecules):
            mer_i = "mer_" + str(i) + ".xyz"
            with open(mer_i, 'w') as fp:
                fp.write("  %d\n" % len(mers_list[i]))
                fp.write("xyz for cooperativity\n")
                for j in mers_list[i]:
                    fp.write("%s %f %f %f\n" % (j[0], j[1], j[2], j[3]))
            mopac_inst = mopac.Mopac(mer_i)
            prepare_input_status = mopac_inst.prepare_input(keyword="PM7 1SCF")
            optimize = mopac_inst.optimize()
            energy = mopac_inst.get_energy()
            energy_list.append([mer_i, energy])
        return energy_list
    def make_dimers(self, mers_list):
        n_of_molecules = len(mers_list)
        no_of_atoms_in_each_mer = len(mers_list[0])
        energy_list = []
        for i in range(n_of_molecules-1):
            for j in range(i+1, n_of_molecules):
                mer_ij = "mer_" + str(i) + "_" + str(j) +".xyz"
                with open(mer_ij, 'w') as fp:
                    fp.write("  %d\n" % (2 * no_of_atoms_in_each_mer))
                    fp.write("xyz for cooperativity\n")
                    for k in mers_list[i]:
                        fp.write("%s %f %f %f\n" % (k[0], k[1], k[2], k[3]))
                    for l in mers_list[j]:
                        fp.write("%s %f %f %f\n" % (l[0], l[1], l[2], l[3]))
                mopac_inst = mopac.Mopac(mer_ij)
                prepare_input_status = mopac_inst.prepare_input(keyword="PM7 1SCF")
                optimize = mopac_inst.optimize()
                energy = mopac_inst.get_energy()
                energy_list.append([mer_ij, energy])
        return energy_list
    def single_point_nmer(self,xyz_file):
        mopac_inst = mopac.Mopac(xyz_file)
        prepare_input_status = mopac_inst.prepare_input(keyword="PM7 1SCF")
        optimize_status = mopac_inst.optimize()
        energy = mopac_inst.get_energy()
        return energy
    def calculate(self, monomers_energy_table, dimers_energy_table, nmer_energy):
        tot_mon_en = 0
        for i in monomers_energy_table:
            tot_mon_en += i[1]
        d_nmer_en = nmer_energy - tot_mon_en
        tot_di_en = 0
        for i in dimers_energy_table:
            tot_di_en += i[1]
        d_dimer_en = tot_di_en - (2 * tot_mon_en)
        dde = d_nmer_en - d_dimer_en
        return dde
    def calculate_cooperativity(self):
        cooperativity_list = []
        self.get_all_xyz_files()
        self.make_dirs_for_all_xyz_files()
        self.cp_xyz_files_to_dirs()
        for i in self.coop_dirs:
            path_name = os.path.join(self.path, i)
            working_dir = os.getcwd()
            os.chdir(path_name)
            xyz_file = glob.glob('*.xyz')
            xyz_file_with_path = os.path.join(path_name, xyz_file[0])
            mers_list = self.read_xyz_c(xyz_file_with_path)
            monomers_energy_list = self.make_monomers(mers_list)
            dimers_energy_list = self.make_dimers(mers_list)
            nmer_energy = self.single_point_nmer(xyz_file_with_path)
            os.chdir(working_dir)
            cooperativity = self.calculate(monomers_energy_list, dimers_energy_list, nmer_energy)
            cooperativity_list.append([i, cooperativity])
        return cooperativity_list


def decide_analysis(input_arguments):
    """This function is the decider function. It will decide which types of 
       calculation is asked for."""
    analysis_want = None
    if "-e" in input_arguments:
        if len(input_arguments) > 1:
            analysis_want = "comparative_lowest_en"
        else:
            print "please provide a path to calculate comparative lowest energy"
            sys.exit(0)
    elif "-c" in input_arguments:
        if len(input_arguments) > 1:
            analysis_want = "cooperativity"
        else:
            print "please provide a path to calculate cooperativity"
            sys.exit(0)
    else:
        analysis_want = "lowest_energy"
    return analysis_want


def lowest_en(paths):
    lowest_en_list = []
    for i in paths:
        abs_path_i = os.path.abspath(i)
        arc_files = files_dirs.get_files("arc", destination=abs_path_i+"/")
        energies = {j: cluster.getEnergy(os.path.join(abs_path_i,j)) for j in arc_files}
        lowest_en = sorted(energies.items(), key=operator.itemgetter(1))[0]
        lowest_en_list.append([i, lowest_en[0], lowest_en[1]])
    return lowest_en_list


def report_lowest_en(result_en_list):
    print_paths = []
    print "{:^25} {:^25} {:^10}".format('path', 'filename', 'energy')
    for i in result_en_list:
        path_last_two = i[0].split('/')[-2:]
        #print_paths.append(path_last_two[0]+"/"+path_last_two[1])
        print "{:^25} {:^25} {:^10}".format(path_last_two[0]+"/"+path_last_two[1], i[1], i[2])
    return


def  get_lowest_compare(dir_trees):
    for i in dir_trees:
        print i[0], ":"
        for dirs in range(len(i[1])):
            # this subdirs_path decision is too specific for this project.
            # the correct way would be:
            #for dirs in i[1]:
            #   subdirs_path = os.path.join(i[0],dirs)
            subdirs_path = os.path.join(i[0], "aggregate_" + str(dirs))
            lowest = lowest_en([subdirs_path])
            print "aggregate_" + str(dirs), lowest[0][1], lowest[0][2]
    return 


def print_lowest_geom(dir_trees):
    for i in dir_trees:
        xyz_file_name = "allxyz_" + os.path.basename(os.path.abspath(i[0])) + ".xyz"
        with open(xyz_file_name, 'w') as fp:
            for dirs in range(len(i[1])):
                subdirs_path = os.path.join(i[0], "aggregate_" + str(dirs))
                lowest = lowest_en([subdirs_path])
                xyz_of_lowest = cluster.getCoordinates(os.path.join(lowest[0][0], lowest[0][1]))
                fp.write(" %d\n" % len(xyz_of_lowest))
                fp.write("energy = %f\n" % lowest[0][2])
                for j in xyz_of_lowest:
                    fp.write("%s %f %f %f\n" % (j[0], j[1], j[2], j[3]))
    return


def main():
    input_argument = sys.argv[1:]
    analysis_want = decide_analysis(input_argument)
    if analysis_want == "lowest_energy":
        energy_list = lowest_en(input_argument)
        report_lowest_en(energy_list)
    if analysis_want == "comparative_lowest_en":
        paths = input_argument[1:]
        dir_trees = []
        for i in paths:
            abs_i = os.path.abspath(i)
            sub_agg_dirs = files_dirs.get_dirs_files(destination = i, wildcard = "aggregate_*")
            dir_trees.append([i, sub_agg_dirs])
        get_lowest_compare(dir_trees)
        print_lowest_geom(dir_trees)
    if analysis_want == "cooperativity":
        paths = input_argument[1:]
        for i in paths:
            abs_path = os.path.abspath(i)
            coop = cooperativity(abs_path, 3)
            print i, ":"
            cooperativity_list = coop.calculate_cooperativity()
            for i in cooperativity_list:
                print i[0],",", i[1]


if __name__ == "__main__":
    main()

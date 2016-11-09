'''
nucleation_genetic.py - interface to run aggragation program in PyAR
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

This file will first generate n random orientation. It will then optimize
each of these orientation. Then, it will choose the lowest energy orientation.
It will then again increase the number one to add one more aggregate.
'''

import time
import sys
import pprint
sys.dont_write_bytecode = True
import os
import numpy as np
from PyAR import species
from PyAR.IO import files_dirs
from PyAR.IO import xyz_file
from PyAR.utils.vector_gen import vector_generator
from PyAR.wrapper import mopac
import subprocess as subp
from PyAR import cluster

def Usage():
    print
    print "Usage: ", sys.argv[0], "molecule1.xyz molecule2.xyz for solvation."
    print "Usage: ", sys.argv[0], "molecule.xyz for clusters"
    print
    sys.exit()


def get_min_energy(dict_name):
    """
    :param dict_name: This is the dictionary which contains energy and
    file names as one to one mapping.
    :return:It will return the energy value and the string name.
    """
    min_value = min(dict_name.itervalues())
    min_keys = [k for k in dict_name if dict_name[k] == min_value]
    return min_value, min_keys


def run_cluster(fragment_1_all_set, fragment_2_all_set, tabu_file, xyz_file_adder, no_of_orientation, aggregates, atoms_name,\
                translational_vector_range, tvector_dimension, rvector_dimension):
    print "Running the Cluster calculation"
    tot_atoms = len(atoms_name)
    dict_struct_energy = {}
    for orientation in range(no_of_orientation):
        xaxis = [translational_vector_range[0], translational_vector_range[1], 1]
        yaxis = [translational_vector_range[0], translational_vector_range[1], 1]
        zaxis = [translational_vector_range[0], translational_vector_range[1], 1]
        all_trans = [xaxis, yaxis, zaxis]
        all_rot = [[0, 360, 12], [0, 180, 6], [0, 360, 12]]
        start_xyz_file_name = "start_"+str(aggregates)+"_"+xyz_file_adder+"_"+str(orientation)+".xyz"
        vectors = vector_generator.vector_generator(tvector_dimension, rvector_dimension)
        tvectors, rvectors = vectors.trans_rot_vector_generator(tabu_file, all_trans, all_rot, 0.1, 30.0)
        new_rot_coords_of_frag_2 = species.rotate_euler(fragment_2_all_set, rvectors[:])
        final_coords_of_frag_2 = species.translate(atoms_name, fragment_1_all_set, new_rot_coords_of_frag_2[:], tvectors[:])
        all_coords = np.concatenate((fragment_1_all_set, final_coords_of_frag_2), axis=0)
        fstart = open(start_xyz_file_name, 'w')
        fstart.write("%d\n" % tot_atoms)
        fstart.write("    \n")
        for n in range(len(all_coords)):
            fstart.write("%-5s %17.10f %17.10f %17.10f \n" % (atoms_name[n], all_coords[n][0], all_coords[n][1], all_coords[n][2]))
        fstart.close()
        mopac_inst = mopac.Mopac(start_xyz_file_name)
        xyz_file_result_name = "result_"+str(aggregates)+"_"+xyz_file_adder+"_"+str(orientation)
        prepare_input_status = mopac_inst.prepare_input(keyword="PM7 PRECISE SINGLET LET DDMIN=0.0 CYCLES=10000 THREADS=8")
        optimize_status = mopac_inst.optimize()
        if optimize_status != 0:
            print "optimize status: ", optimize_status
            print "The optimization for geometry:", start_xyz_file_name, "has failed."
            print "New geometry will be optimized and it will be discarded."
            continue
        #if files_dirs.file_exists_check_with_return_status("start_"+str(aggregates)+"_"+str(orientation)+".arc"):
        #    energies = mopac_inst.get_energy()
        #else:
        #    continue
        energies = mopac_inst.get_energy()
        en_kcal = energies[0]
        en_kj = energies[1]
        print "xyz_file_result_name is: ", xyz_file_result_name
        mopac_inst.extract_xyz(xyz_file_result_name+".xyz")
        dict_struct_energy[xyz_file_result_name] = en_kcal
    pprint.pprint(dict_struct_energy)
    return

def main():
    start_time = time.time()
    stopfile="STOP"
    # check for command line arguments:
    mode = None
    if len(sys.argv) == 2:
        mode = "nucleation"
    elif len(sys.argv) == 3:
        mode = "solvation"
    else:
        Usage()
    # initial variable start:

    translational_vector_range = [-20.0, 20.0]
    tvector_dimension = 3
    rvector_dimension = 3
    if mode == "nucleation":
        natom1, title1, atoms_name_fragment_1, coords1 = xyz_file.read_xyz(sys.argv[1])
        coords2 = coords1[:]
        atoms_name_fragment_2 = atoms_name_fragment_1[:]
    elif mode == "solvation":
        natom1, title1, atoms_name_fragment_1, coords1 = xyz_file.read_xyz(sys.argv[1])
        natom12, title12, atoms_name_fragment_2, coords2 = xyz_file.read_xyz(sys.argv[2])
    else:
        print "Variable 'mode' is not recognised."
    no_of_orientation = 300
    no_of_aggregation = 20
    fragment_2_all_set = species.set_origin(coords2[:])
    tabu_file_base = "trtabu"

    # the for loop starts:
    for aggregates in range(no_of_aggregation-1):
        print "Aggregate number: ", aggregates
        aggregate_home = "aggregate_"+str(aggregates)
        files_dirs.make_directories(aggregate_home)
        if aggregates == 0:
            tabu_file = tabu_file_base + "_" +str(aggregates)
            xyz_file_adder = "0"
            atoms_name = xyz_file.flatten([atoms_name_fragment_1, atoms_name_fragment_2])
            fragment_1_all_set = species.set_origin(coords1[:])
            run_cluster(fragment_1_all_set[:], fragment_2_all_set[:], tabu_file, xyz_file_adder, no_of_orientation, aggregates, atoms_name, \
                    translational_vector_range, tvector_dimension, rvector_dimension)
        else:
            best_xyz_files = files_dirs.get_files("xyz")
            #check if input file present in the best xyz file list
            # if exists, remove it from the list.
            for i in sys.argv[1:]:
                # the for loop is for multiple input xyz files
                try:
                    # to delete the input xyz file names from the list
                    input_index = best_xyz_files.index(i)
                    del best_xyz_files[input_index]
                except:
                    pass
            print "clustering for: ", best_xyz_files
            for i in range(len(best_xyz_files)):
                tabu_file = tabu_file_base + str(i) + "_" + str(aggregates)
                xyz_file_adder = str(i)
                print "best1:", best_xyz_files[i]
                natom1, title1, atoms_name_fragment_1, coords1 = xyz_file.read_xyz(best_xyz_files[i])
                atoms_name = xyz_file.flatten([atoms_name_fragment_1, atoms_name_fragment_2])
                fragment_1_all_set = species.set_origin(coords1[:])
                run_cluster(fragment_1_all_set[:], fragment_2_all_set[:], tabu_file, xyz_file_adder, no_of_orientation, aggregates, atoms_name, \
                    translational_vector_range, tvector_dimension, rvector_dimension)
        arcfiles = files_dirs.get_files("arc")
        pprint.pprint(arcfiles)
        best_clusters = cluster.quadsplit(arcfiles)
        print "The best clusters found are follows:"
        pprint.pprint(best_clusters)
        files_dirs.mmove("result", "start", "trtabu", aggregate_home)
        for files in best_clusters:
            file_name = str(os.path.basename(files))
            source_file = "./"+aggregate_home+"/"+file_name
            files_dirs.scopy(source_file, "./")
            base_file_name = os.path.splitext(file_name)[0] # this file is the xyz file corresponding to arc
            mopac.arc2xyz(file_name, base_file_name+".xyz")
            files_dirs.delete_files(base_file_name+".arc")
        files_dirs.check_stop(stopfile)
        sys.stdout.flush()
    print "\n>>>>TIME: ", time.time() - start_time, "seconds.<<<<\n"


if __name__ == "__main__":
    main()

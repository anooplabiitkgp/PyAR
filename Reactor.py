'''
Reactor.py - interface to run binary reaction in PyAR
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

import time
import sys
sys.dont_write_bytecode = True
import os
import numpy as np
from PyAR import species
from PyAR.IO import files_dirs
from PyAR.IO import xyz_file
from PyAR.utils import fragment_make
from PyAR.utils import structure_check
from PyAR.utils import gamma
from PyAR.utils.vector_gen import vector_generator
from PyAR.wrapper import turbomole


def Usage():
    print
    print "Usage: ", sys.argv[0], "molecule1.xyz molecule2.xyz"
    print
    sys.exit()


def main():
    start_time = time.time()

    # check for command line arguments:

    if len(sys.argv) < 3:
        Usage()

    # initialize variables

    N = N0 = 0
    Nmax = 30
    tabu_list_for_translation = [[], [], []]
    molecule1 = sys.argv[1]
    molecule2 = sys.argv[2]
    tvector_dimension = 3
    rvector_dimension = 3
    fragment_file = "fragment"
    define_inp_file = "define.inp"
    tabu_file = "trtabu"
    reactant_save_directory = "starting_structures"
    full_opt_geom = "full_opt"
    product_save_directory = "product_geom"
    trajectory_save_directory = "trajectory_geom"
    translational_vector_range = (-8.0, 8.0)
    jobex_file = "jobexd_sn"
    jobex_path = "/home/surajit/bin"
    jobex = os.path.join(jobex_path, jobex_file)
    stopfile = "stop_pyreactor"

    # making necessary files

    fragment_make.fragment_make(molecule1, molecule2, fragment_file)
    files_dirs.make_directories(reactant_save_directory, product_save_directory, trajectory_save_directory)

    # checking files for starting calculation

    files_dirs.important_file_check(fragment_file, "define.inp", jobex)
    files_dirs.delete_files("energy", "control", "gradient", "results_of_afir.out", tabu_file)

    # coordinate operations

    natom1, title1, atoms1, fragment_1_coordinates_unset = xyz_file.read_xyz(molecule1)
    natom2, title2, atoms2, fragment_2_coordinates_unset = xyz_file.read_xyz(molecule2)
    atoms_name = xyz_file.flatten([atoms1, atoms2])
    tot_atoms = len(atoms_name)
    print "atoms of molecule 1: \n", atoms1
    print "atoms of molecule 2: \n", atoms2
    print "All the atoms present in the two fragments are:", atoms_name
    print "total number of atoms in the system is ", tot_atoms
    fragment_1_all_set = species.set_origin(fragment_1_coordinates_unset[:])
    fragment_2_all_set = species.set_origin(fragment_2_coordinates_unset[:])
    while N-N0 < Nmax:
        print "\n\ncycle = ", N, "\n\n"
        files_dirs.delete_files('converged')
        xaxis = [translational_vector_range[0], translational_vector_range[1], 1]
        yaxis = [translational_vector_range[0], translational_vector_range[1], 1]
        zaxis = [translational_vector_range[0], translational_vector_range[1], 1]
        all_trans = [xaxis, yaxis, zaxis]
        all_rot = [[0, 360, 12], [0, 180, 6], [0, 360, 12]]
        start_xyz_file = 'start_' + str(N) + '.xyz'
        trajectory_file = 'trj_' + str(N) + '.xyz'
        product_xyz_file = 'product_' + str(N) + '.xyz'
        if len(fragment_2_all_set) < 2:
            print "The second fragment contain only one atom."
            print "No rotation will be applied on fragment 2."
            A = vector_generator.vector_generator(3, 3)
            tvectors, new_tabu_list_trans = A.trans_rot_vector_generator(tabu_file, all_trans, all_rot, 0.2, 40.0)
            for itabu in range(tvector_dimension):
                print "The tabu list generated is", new_tabu_list_trans[itabu]
                tabu_list_for_translation[itabu] = new_tabu_list_trans[itabu][:]
            print "The vector is,", tvectors
            final_coords_of_frag_2 = species.translate(atoms_name, fragment_1_all_set[:], fragment_2_all_set[:], tvectors[:])
        elif len(fragment_2_all_set) >= 2:
            print "The second fragment contain more than one atom."
            print "Rotation will be applied on second fragment."
            vectors = vector_generator.vector_generator(tvector_dimension, rvector_dimension)
            tvectors, rvectors = vectors.trans_rot_vector_generator(tabu_file, all_trans, all_rot, 0.2, 40.0)
            print "The rotational vector is:", rvectors
            print "The translational vector is:", tvectors
            new_rot_coords_of_frag_2 = species.rotate_euler(fragment_2_all_set ,rvectors[:])
            final_coords_of_frag_2 = species.translate(atoms_name, fragment_1_all_set, new_rot_coords_of_frag_2[:], tvectors[:])
        all_coords = np.concatenate((fragment_1_all_set, final_coords_of_frag_2), axis=0)
        fstart = open(start_xyz_file, 'w')
        fstart.write("%d\n" % tot_atoms)
        fstart.write("    \n")
        for n in range(len(all_coords)):
            fstart.write("%-5s %17.10f %17.10f %17.10f \n" % (atoms_name[n], all_coords[n][0], all_coords[n][1], all_coords[n][2]))
        fstart.close()
        files_dirs.scopy(start_xyz_file, reactant_save_directory)
        gamma_values = gamma.get_systematic_gamma(100, 2000, 100)
        print "the gamma list generated for random geom no. ", N, "is: ", gamma_values
        turbomole.make_coord(start_xyz_file)
        turbomole.run_define(define_inp_file)
        sys.stdout.flush()
        for igamma in gamma_values:
            optimized_geom = "tmp.xyz"
            gamma_free_opt_dir = "./tmp_afir"
            print "gamma = ", igamma
            sys.stdout.flush()
            jobex_status = turbomole.run_jobex(jobex, dftd3="-zero", afir=igamma, ri="-ri", cycle=1000)
            print "jobex with gamma,", igamma, "is completed."
            convergence_switch = turbomole.check_convergence(".")
            files_dirs.delete_files("converged", "not.converged", "dscf_problem")
            print "convergence_switch, jobex_status are:", convergence_switch, jobex_status
            sys.stdout.flush()
            if convergence_switch != "converged" or jobex_status != 0:
                print "AFIR optimization was not successful for cycle", N, "."
                print "New geometry will be generated."
                sys.stdout.flush()
                break
            turbomole.make_last_geometry(outfile=optimized_geom)
            init_geom_comp_smi = structure_check.structure_comparison(start_xyz_file, optimized_geom, cheminfo="smi")
            init_geom_comp_inchi = structure_check.structure_comparison(start_xyz_file, optimized_geom, cheminfo="inchi")
            print "The comparison status for smi and inchis are:", init_geom_comp_smi, init_geom_comp_inchi, "respectively."
            sys.stdout.flush()
            if init_geom_comp_inchi == "same" or init_geom_comp_smi == "same":
                print "The geometry optimization for gamma", igamma, "leads to same product."
                print "Optimization will be done for the next gamma."
                sys.stdout.flush()
                continue
            turbomole.prepare_for_optimization(optimized_geom, gamma_free_opt_dir, define_inp_file, fragment=fragment_file)
            jobex_free_status = turbomole.perform_optimization_in_different_location(gamma_free_opt_dir, optimized_geom, define_inp_file, jobex, dftd3="-zero", ri="-ri", cycle=50)
            free_convergence_switch = turbomole.check_convergence(gamma_free_opt_dir)
            print "jobex with zero gamma has been completed. Structures will be compared."
            if free_convergence_switch == "dscf_problem":
                print "The scf convergence for the last geometry was not satisfied within limited cycle."
                print "New Geometry will be generated for the RUN."
                sys.stdout.flush()
                break
            elif free_convergence_switch == "converged" or free_convergence_switch == "notconverged":
                print "The partial optimization was successful. The last geometry will be taken out for comparison."
                files_dirs.scopy(gamma_free_opt_dir + "/last.xyz", optimized_geom)
                files_dirs.delete_directories(gamma_free_opt_dir)
                start_last_smi = structure_check.structure_comparison(optimized_geom, start_xyz_file, cheminfo="smi")
                start_last_inchi = structure_check.structure_comparison(optimized_geom, start_xyz_file, cheminfo="inchi")
                sys.stdout.flush()
            print "The smi and inchi status after partial 0 gamma optimizaton are:", start_last_smi, start_last_inchi
            conv_test = None
            sys.stdout.flush()
            if start_last_smi == "same" or start_last_inchi == "same":
                print "The geometry does not seems to converged to a different minima than the initial one."
                print "new optimization with higher gamma will be started."
                files_dirs.delete_files(optimized_geom)
                sys.stdout.flush()
                continue
            else:
                print "The geometry converged to a local minima that is different from the initial one."
                print "Investigation on the last geometry will be done."
                turbomole.prepare_for_optimization(optimized_geom, full_opt_geom, define_inp_file, fragment=fragment_file)
                turbomole_opt_free = turbomole.perform_optimization_in_different_location(full_opt_geom, optimized_geom, define_inp_file, jobex, dftd3="-zero", ri="-ri")
                free_opt_conv = turbomole.check_convergence(full_opt_geom)
                conv_test = turbomole.converge_check(turbomole_opt_free, free_opt_conv)
                sys.stdout.flush()
            if conv_test == "converged":
                files_dirs.scopy(full_opt_geom+"/last.xyz", optimized_geom)
                same_product_inchi = structure_check.structure_comparison(optimized_geom, start_xyz_file, cheminfo="inchi")
                same_product_smi = structure_check.structure_comparison(optimized_geom, start_xyz_file, cheminfo="smi")
                previous_product_inchi = structure_check.prev_product_comp_status(optimized_geom, "./"+product_save_directory, cheminfo="inchi")
                previous_product_smi = structure_check.prev_product_comp_status(optimized_geom, "./"+product_save_directory, cheminfo="smi")
                print "The last geometry was fully converged. Structure checking will be done now."
                sys.stdout.flush()
            else:
                files_dirs.delete_files(optimized_geom)
                print "The last geometry was not converged within the limited cycles."
                print "program will take new geometry to optimize."
                sys.stdout.flush()
                break
            prod_check_status = structure_check.check_smi_inchi_status(previous_product_inchi, previous_product_smi, same_product_inchi, same_product_smi)
            sys.stdout.flush()
            if prod_check_status == "different":
                print "A new product was found. It will be saved as new product."
                N0 = N
                files_dirs.smove(optimized_geom, product_xyz_file, product_save_directory)
                print "Check the product with name: ", product_xyz_file
                sys.stdout.flush()
                break
            else:
                print "Similarity was found with other products. It will be deleted."
                files_dirs.delete_files(optimized_geom)
                sys.stdout.flush()
                break
        N += 1
        turbomole.make_trajectory_file(trajectory_save_directory+"/"+trajectory_file)
        files_dirs.delete_directories("./tmp_afir", full_opt_geom)
        files_dirs.delete_files("converged", "energy", "gradient", "coord", "control", start_xyz_file, 'last_free.xyz')
        files_dirs.check_stop(stopfile)
        sys.stdout.flush()
    print "\n>>>>TIME: ", time.time() - start_time, "seconds.<<<<\n"


if __name__ == "__main__":
    main()

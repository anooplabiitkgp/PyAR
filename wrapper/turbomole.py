'''
turbomole.py - interface to turbomole program
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
import glob
import shutil
import sys
sys.dont_write_bytecode = True
from PyAR.utils import globals
from PyAR.IO import files_dirs

def prepare_for_optimization(xyz_file, destination_directory, define_inp_file, fragment=None):
    files_dirs.make_directories(destination_directory)
    if fragment is None:
        pass
    else:
        files_dirs.scopy(fragment, destination_directory)
    files_dirs.scopy(xyz_file, destination_directory)
    files_dirs.scopy(define_inp_file, destination_directory)
    current_location = os.getcwd()
    os.chdir(destination_directory)
    make_coord(xyz_file)
    run_define(define_inp_file)
    os.chdir(current_location)
    return


def make_last_geometry(outdir="./", outfile="tmp"):
    """This function will generate last geometry from a turbomole
       optimization. It will use the command t2x -c > /path/to/<file>
    """
    output_path = os.path.join(outdir, outfile)
    with open(output_path, "w") as f:
        subp.check_call(["t2x", "-c", ">"], stdout=f)
    return


def make_coord(xyzfile, outfile='coord'):
    """this function will make coord file from turbomole x2t module
    """
    with open(outfile, 'w') as f:
        subp.check_call(["x2t", str(xyzfile), ">"], stdout=f)
    return


def make_trajectory_file(outfile="tmp"):
    """This function will generate trajectory file in a output
       destination by t2x module
    """
    with open(outfile, "w") as f:
        subp.check_call(["t2x", ">"], stdout=f)
    return


def run_define(define_input_file):
    """This function will run define module with the define_input_file
       the inputs for define
    """
    if os.path.isfile(define_input_file):
        with open('tmp_bash_command', 'w') as f:
            f.write("define<%s>&define.log" % define_input_file)
        with open('define.log', 'w') as f:
            subp.check_call(["bash", "tmp_bash_command"], stdout=f)
        os.remove('tmp_bash_command')
    return


def run_jobex(jobex_script, dftd3="-zero", afir=0.0, ri="", cycle=100):
    """ This function will run jobex script by taking the path and the 
        script name.The options are dftd3, afir, ri, cycle etc. The switch 
        will take in the same format as jobex take.
    """
    path, jobex_file = os.path.split(jobex_script)
    outfile = jobex_file + ".log"
    with open(outfile, 'w') as f:
        out = subp.Popen([jobex_script, "-a", str(afir), dftd3, ri, "-c", str(cycle)], stdout=f, stderr=f)
    output, error = out.communicate()
    poll = out.poll()
    exit_status = out.returncode
    return exit_status


def energy(location, unit="au"):
    """This function will return the energy in Hartree unit by default.
       To get energy in kcal/mol, set option in kcal
    """
    energy_file = os.path.join(location, "energy")
    with open(energy_file) as f:
        all_lines = f.readlines()
    last_energy_line = all_lines[-2]
    energy = float(last_energy_line.split()[1])
    if unit == "au":
        return energy
    else:
        return energy * globals.AUTOKCAL


def run_scf(turbo_scf_module, path_to_run="./"):
    """ This function will calculate scf of a given geometry.
        Note that, it will first go to the path where scf has
        to be calculated. After calculation, it will return to
        the current directory. Use with caution. Check if the
        necessary files are present before run this function.
        Check turbomole manual for details.
    """
    current_location = os.getcwd()
    os.chdir(path_to_run)
    out_file = turbo_scf_module + str(".log")
    with open(out_file) as f:
        subp.check_call([turbo_scf_module], stdout=f)
    os.chdir(current_location)
    return


def run_grad(turbo_grad_module, path_to_run="./"):
    """ This function will calculate rdgrad of a given geometry
        after ridft calculation. The default path where to calculate
        is set in the current location. If necessary, we can also
        calculate it in a different location. For details, see
        turbomole manual.
    """
    current_location = os.getcwd()
    os.chdir(path_to_run)
    out_file = turbo_grad_module + str(".log")
    with open(out_file) as f:
        subp.check_call([turbo_grad_module], stdout=f)
    os.chdir(current_location)
    return


def run_dftd3_with_force(dftd3_code, path_to_run="./", *options):
    """This function is for running dftd3 code with options and also for
     running afir force. Primarily is for AFIR 
    """
    current_location = os.getcwd()
    os.chdir(path_to_run)
    out_file = "dftd3.log"
    with open(out_file) as f:
        subp.check_call([dftd3_code], stdout=f)
    os.chdir(current_location)
    return


def run_statpt(statpt_code, path_to_run="./"):
    """This function is for running turbomole module statpt. 
    """
    current_location = os.getcwd()
    os.chdir(path_to_run)
    out_file = "statpt.log"
    with open(out_file) as f:
        subp.check_call([statpt_code], stdout=f)
    os.chdir(current_location)
    return


def check_convergence(path_to_check):
    """This function will check if the file 'converged' exists in a
       given directory
    """
    conv_file = os.path.join(path_to_check, 'converged')
    dscf_prob_file = os.path.join(path_to_check, 'dscf_problem')
    non_converge = os.path.join(path_to_check, 'not.converged')
    if os.path.exists(conv_file):
        conv_test =  "converged"
    elif os.path.exists(dscf_prob_file):
        conv_test = "dscf_problem"
    elif os.path.exists(non_converge):
        conv_test = "notconverged"
    else:
        conv_test = "none"
    return conv_test


def converge_check(turbomole_opt_free, free_opt_conv):
    conv_test = "noconverged"
    if turbomole_opt_free == 0 and free_opt_conv == "converged":
        conv_test = "converged"
    return conv_test


def perform_optimization_in_different_location(path_to_location, xyz_file, define_inp_file, jobex_script, dftd3 = "no", ri = "", cycle = 3000):
    """This function will optimize a geometry in a given location. It will
       also create the last geometry from t2x -c module of turbomole. It
       will first check if a xyz file is present in this location and a
       define.inp file is present in this location or not. If not present,
       it will give error and the program will break.
    """
    current_location = os.getcwd()
    os.chdir(path_to_location)
    return_code = run_jobex(jobex_script, dftd3 = dftd3, afir=0.0, ri=ri, cycle=cycle)
    if return_code == 0:
        make_last_geometry(outfile="last.xyz")
        optimization_status = 0
    else:
        optimization_status = 1
    os.chdir(current_location)
    return optimization_status


def delete_non_conv_directory():
    """This function will check for non converge directory, MPI-TEMPDIR-00*
    """
    all_MPI_TMPDIRS = glob.glob('MPI-TEMPDIR*')
    for directories in all_MPI_TMPDIRS:
        files_dirs.delete_directories(directories)
    return


def delete_non_conv_files():
    """This function will check for non converge file, dscf_problem
    """
    files_dirs.delete_files('dscf_problem', "not.converged")
    return


def main():
    pass


if __name__ == "__main__":
    main()

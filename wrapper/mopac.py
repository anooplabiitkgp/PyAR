'''
mopac.py - interface to mopac program
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

import os, sys
import subprocess as subp


class Mopac(object):

    def __init__(self, xyz_file):
        self.xyz_file = xyz_file
        self.xyz_file_without_extension = os.path.splitext(xyz_file)[0]

    def prepare_input(self, keyword=""):
        """
        :param keyword: this is the keyword for optimizations. This parameter
        should be a strings of characters which are mopac keywords
        :return: It will not return anything. It will prepare the input file for
        the purpose given in the keyword. Note that babel will be used to prepare
        the input(.mop) file.
        """
        if not keyword:
            keyword_line = "-xkPM7"
        elif keyword:
            keyword_line = "-xk" + keyword
        else:
            print "keyword:", keyword, "not recognized."
            print "program will stop"
            sys.exit()
        with open('tmp.log', 'w') as fminp:
            out = subp.Popen(["babel", "-ixyz", self.xyz_file, "-omop", self.xyz_file_without_extension+".mop", keyword_line], stdout=fminp, stderr=fminp)
        output, error = out.communicate()
        poll = out.poll()
        exit_status = out.returncode
        os.remove('tmp.log')
        return exit_status

    def optimize(self):
        """
        :return:This object will return the optimization status. It will
        optimize a structure.
        """
        with open(self.xyz_file_without_extension+".log", 'w') as fopt:
            out = subp.Popen(["mopac", self.xyz_file_without_extension+".mop"], stdout=fopt, stderr=fopt)
        output, error = out.communicate()
        poll = out.poll()
        exit_status = out.returncode
        return exit_status

    def get_energy(self):
        """
        :return:This object will return energy from a mopac calculation. It will return both the kj/mol and
        kcal/mol units.
        """
        mopac_arc_file = self.xyz_file_without_extension+".arc"
        en_kcal = 0.0
        en_kj = 0.0
        try:
            with open(mopac_arc_file,'r') as arc_out:
                arc_cont = arc_out.readlines()
            for lines in arc_cont:
                if "HEAT OF FORMATION" in lines:
                    line_cont = lines.split('=')
                    en_kcal = float(line_cont[1].split()[0])
                    en_kj = float(line_cont[2].split()[0])
                    break
        except:
            print "Warning: File ", mopac_arc_file, "was not found."
        #with open(mopac_arc_file) as arc_out:
        #    arc_cont = arc_out.readlines()
        #for lines in arc_cont:
        #    if "HEAT OF FORMATION" in lines:
        #        line_cont = lines.split('=')
        #        en_kcal = float(line_cont[1].split()[0])
        #        en_kj = float(line_cont[2].split()[0])
        #        break
        return en_kcal, en_kj

    def extract_xyz(self, out_file):
        """
        :param out_file: This is the output file in which the final xyz coordinates will be
        written
        :return: It will return nothing.
        """
        mopac_arc_file = self.xyz_file_without_extension+".arc"
        try:
            with open(mopac_arc_file) as arc_out:
                arc_cont = arc_out.readlines()
            for lines in arc_cont:
                if 'Empirical Formula:' in lines:
                    natoms = int(lines.split()[-2])
            coordinates = arc_cont[-(natoms+1):-1]
            with open(out_file, 'w') as fout:
                fout.write("%6d\n" % natoms)
                fout.write("   \n")
                for i in coordinates:
                    fout.write("%-6s %10.3f %10.3f %10.3f\n" % \
                               (str(i.split()[0]), float(i.split()[1]), float(i.split()[3]), float(i.split()[5])))
                fout.write("\n")
        except:
            print "Warning: File ", mopac_arc_file, "was not found."

        #with open(mopac_arc_file) as arc_out:
        #    arc_cont = arc_out.readlines()
        #for lines in arc_cont:
        #    if 'Empirical Formula:' in lines:
        #        natoms = int(lines.split()[-2])
        #coordinates = arc_cont[-(natoms+1):-1]
        #with open(out_file, 'w') as fout:
        #    fout.write("%6d\n" % natoms)
        #    fout.write("   \n")
        #    for i in coordinates:
        #        fout.write("%-6s %10.3f %10.3f %10.3f\n" % \
        #                   (str(i.split()[0]), float(i.split()[1]), float(i.split()[3]), float(i.split()[5])))
        #    fout.write("\n")
        return


def arc2xyz(arc_file, xyz_file):
    """This function extract xyz file from arc file"""
    print "arc_file:", arc_file
    print "xyz_file: ", xyz_file
    with open(arc_file) as fp:
        arc_cont = fp.readlines()
    for lines in arc_cont:
        if 'Empirical Formula:' in lines:
            natoms = int(lines.split()[-2])
    coordinates = arc_cont[-(natoms+1):-1]
    with open(xyz_file, 'w') as fp:
        fp.write("%d\n" % natoms)
        fp.write("%s\n" % "Geometry from arc file")
        for i in coordinates:
            fp.write(" %s \t %.8f \t %.8f \t %.8f \n" % (str(i.split()[0]),float(i.split()[1]),float(i.split()[3]),float(i.split()[5])))
    return


def main():
    pass


if __name__ == "__main__":
    main()

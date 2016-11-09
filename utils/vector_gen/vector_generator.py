'''
vector_generator.py - generating random and tabu vector
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
from PyAR.utils.vector_gen import tabu
from PyAR.utils.vector_gen import distance_tabu


class VectorGeneratorError(Exception):
    """Base class for exceptions in this module"""
    pass


class DimensionMatchError(VectorGeneratorError):
    """Exception raised for errors in vector dimension mismatch"""
    pass


#.........................................................................


class vector_generator(object):
    def __init__(self, t_dimension, r_dimension):
        self.t_dimension = t_dimension
        self.r_dimension = r_dimension

    def check_null_vector(self, vector):
        dimension = len(vector)
        null_status = False
        null_vector = []
        for i in range(dimension):
            null_vector.append(0.0)
        if vector == null_vector:
            null_status = True
        return null_status

    def get_random_vector(self, dimension, ranges):
        vector = []
        for i in range(dimension):
            vector_axis = tabu.generate_random_number(ranges[i][0], ranges[i][1], ranges[i][2])
            vector.append(vector_axis)
        return vector[:]

    def trans_rot_vector_generator(self, tabu_file_name, tranges, rranges, tcutof, rcutof):
        if (len(tranges) != self.t_dimension) or (len(rranges) != self.r_dimension):
            raise DimensionMatchError("translational or rotational vector ranges gives different dimensions.")
        distance_rotation_tabu_status = True
        orientation_count = 0
        cutoff_orientation_count = 500
        while distance_rotation_tabu_status:
            tvector = []
            rvector = []
            if orientation_count > cutoff_orientation_count:
                print "maximum limit of translational and rotational search limit crossed ", cutoff_orientation_count
                print "exiting the program"
                sys.stdout.flush()
                sys.exit(0)
            tvector = self.get_random_vector(self.t_dimension, tranges)
            if self.check_null_vector(tvector) is True:
                distance_rotation_tabu_status = True
                continue
            rvector = self.get_random_vector(self.r_dimension, rranges)
            distance_rotation_tabu_status = distance_tabu.check_rotation_translation(tvector[:], rvector[:], tabu_file_name, tcutof, rcutof)
            orientation_count += 1
        tabu.update_tabu_file(tvector[:], rvector[:], tabu_file_name)
        return tvector[:], rvector[:]

    def normal_vector_generator(self, tranges, rranges):
        """This method which generates vector randomly. It takes details about
        the vector components as a list of lists by the variable ranges. The
        variable ranges should also contain the interval of the ranges.
        example:
        [[-x1,+x1,dx],[-y1,+y1,dy],[-z1,+x1,dz]] for an 3 dimensional vector.
        """
        if (len(tranges) != self.t_dimension) or (len(rranges) != self.r_dimension):
            raise DimensionMatchError("translational or rotational vector ranges gives different dimensions.")
        tvector = self.get_random_vector(self.t_dimension, tranges)
        rvector = self.get_random_vector(self.r_dimension, rranges)
        return tvector[:], rvector[:]


def main():
    pass


if __name__ == "__main__":
    main()

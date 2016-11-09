'''
distance_tabu.py - some geometric calculations to consider as tabu 
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

import math
from PyAR.IO import files_dirs
from PyAR import species

def get_distance(point1, point2):
    distance_square = 0.0
    if len(point1) == len(point2):
        for i in range(len(point1)):
            distance_square = distance_square + (float(point1[i])-float(point2[i])) ** 2
        distance = math.sqrt(distance_square)
    else:
        raise IOError("points %s and point %s are of different dimensions" % (point1, point2))
    return distance

def get_normalized_vectors_distance(normalize_vector1, normalize_vector2):
    normalized_vector_difference = species.get_vector(normalize_vector1, normalize_vector2, 12)
    distance = species.get_norm(normalized_vector_difference)
    return distance


def get_normal_vector(point_1, point_2):
    vector_21 = species.get_vector(point_1[:], point_2[:], 21)
    normalize_vector_21 = species.get_normalize_vector(vector_21)
    return normalize_vector_21

def check_normalize_distance(vector1, vector2, cutoff_distance):
    distance = get_normalized_vectors_distance(vector1, vector2)
    if distance <= cutoff_distance:
        distance_tabu = True
    elif distance > cutoff_distance:
        distance_tabu = False
    else:
        raise IOError("distance = %s Cannot process" % distance)
    return distance_tabu

def get_all_lines(tabu_file_name):
    with open(tabu_file_name) as fp:
        all_lines = fp.readlines()
    return all_lines

def check_roational_difference(rvector, rot_point, rcutoff):
    if len(rvector) != len(rot_point):
        raise IOError("dimension fo %s, %s are different" % (rvector, rot_point))
    rotational_status = True
    for i in range(len(rvector)):
        if (math.fabs(rvector[i] - rot_point[i])) > rcutoff:
            rotational_status = False
            break
    return rotational_status


def get_tvector_rvector(line_content):
    tvector = [float(line_content[0]), float(line_content[1]), float(line_content[2])]
    rvector = [float(line_content[3]), float(line_content[4]), float(line_content[5])]
    return tvector, rvector

def check_rotation_translation(trans_point, rot_point, tabu_file_name, tcutoff, rcutoff):
    translation_rotation_tabu_status = False
    if files_dirs.file_exists_check_with_return_status(tabu_file_name) is False:
        translation_rotation_tabu_status = False
    elif files_dirs.file_exists_check_with_return_status(tabu_file_name) is True:
        all_lines = get_all_lines(tabu_file_name)
        for lines in all_lines:
            tvector, rvector = get_tvector_rvector(lines.split())
            normalize_tvector = get_normal_vector(tvector, [0., 0., 0.])
            normalize_trans_point = get_normal_vector(trans_point, [0., 0., 0.])
            distance_status = check_normalize_distance(normalize_trans_point, normalize_tvector, tcutoff)
            rotation_status = check_roational_difference(rvector, rot_point, rcutoff)
            #print "distance and rotation tabu status:", distance_status, rotation_status
            if (distance_status and rotation_status) is True:
                translation_rotation_tabu_status = True
                break
    return translation_rotation_tabu_status


def check_distance(point1, tabu_file_name, cutoff_distance):
    normalize_vector_1 = get_normal_vector(point1[:], [0.0, 0.0, 0.0])
    distance_criteria = False
    if files_dirs.file_exists_check_with_return_status(tabu_file_name):
        with open(tabu_file_name) as tabu_file:
            all_tabu_list = [[int(x) for x in line.split()] for line in tabu_file]
        for vector2 in all_tabu_list:
            normalize_vector2 = get_normal_vector(vector2[:], [0.0, 0.0, 0.0])
            distance_criteria = check_normalize_distance(normalize_vector_1[:], normalize_vector2[:], cutoff_distance)
            if distance_criteria is True:
                break
    else: distance_criteria = False
    return distance_criteria


def main():
    pass


if __name__=="__main__":
    main()

'''
species.py - geometrical operations for orientation in PyAR
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

import numpy as np
import math

#
# This module is for generating geometry orientation.
#

#
# Input is the rotational vector/translational vector and a set of coordinates.
# It also has additional methods like set_origin -- which will transform
# the coordinates to (0,0,0).
#
#    All the return elements will be in the form of a list.
#


def centroid(X):
    """This function will return the center of coordinates.
    """
    X = np.array(X)
    c = sum(X)/len(X)
    c = c.tolist()
    return c[:]


def set_origin(coordinates, displacement=None):
    """coordinates --  set of coordinates to be set at origin
       displacement -- the displacement point, i.e., (x,y,z)
                       If nothing given, then the coordinates will be
                       shifted to the center of mass,i.e.,the molecule
                       will be centered in the origin. If any point is
                       given, then the molecule will be centered in this
                       point."""
    if displacement is None:
        displacement = [0., 0., 0.]
    com = np.array(centroid(coordinates[:]))
    coordinates = np.array(coordinates)
    displacement = np.array(displacement)
    coordinates -= com
    coordinates += displacement
    displaced_coordinates = coordinates.tolist()
    return displaced_coordinates[:]


def rotate_euler(coordinates, vector):

    """This function will rotate a molecule by Euler rotation theorem.
    The Z-X-Z' rotation convention is followed. The range of phi, theta
    and psi should be (0,360),(0,180) and (0,360) respectively.

         This function will first translate the molecule to its center
    of mass(centroid). Then, it rotate the molecule and translate
    to its original position. So, to rotate a molecule around the origin,
    (0.,0.,0.), set_origin usage is necessary"""

    phi = float(vector[0]) * math.pi/180
    theta = float(vector[1]) * math.pi/180
    psi = float(vector[2]) * math.pi/180
    COM = np.array(centroid(coordinates))
    coordinates = np.array(set_origin(coordinates[:]))
    D = np.array(((math.cos(phi), math.sin(phi), 0.),
        (-math.sin(phi), math.cos(phi), 0.),
        (0., 0., 1.)))
    C = np.array(((1., 0., 0.),
        (0., math.cos(theta), math.sin(theta)),
        (0., -math.sin(theta), math.cos(theta))))
    B = np.array(((math.cos(psi), math.sin(psi), 0.),
        (-math.sin(psi), math.cos(psi), 0.),
        (0., 0., 1.)))
    A = np.dot(B, np.dot(C, D))
    coordinates = np.dot(A, np.transpose(coordinates))
    final_coordinates = np.transpose(coordinates) + COM
    final_coordinates = final_coordinates.tolist()
    return final_coordinates[:]


def translate_old(coordinates, displacement):

    """This function will translate a molecule to a given point
    The steps to translate a molecule -
        First, set the molecule at origin. Then, translate to the
    given point. Then translate w.r.t. its original com"""

    displacement = np.array(displacement[:])
    coordinates = np.array(coordinates)
    coordinates += displacement
    translated_coordinates = coordinates.tolist()
    return translated_coordinates[:]


def get_norm(vector_ij):
    square_of_norm_21 = 0
    for i in range(len(vector_ij)):
        square_of_norm_21 = square_of_norm_21 + vector_ij[i]*vector_ij[i]
    norm = math.sqrt(square_of_norm_21)
    return norm


def get_vector(point1, point2, direction):
    if len(point2) == len(point1):
        dimension = len(point1)
    else:
        raise IOError("dimension of %s %s are different" % (point1, point2))
    if direction == 12:
        vector = [point1[i] - point2[i] for i in range(dimension)]
    elif direction == 21:
        vector = [point2[i] - point1[i] for i in range(dimension)]
    else:
        raise IOError("direction %s not known please give 12 or 21 as direction." % direction)
    return vector


def get_normalize_vector(vector):
    dimension = len(vector)
    normal = get_norm(vector)
    if normal == 0.0:
        normalized_vector = vector
    else:
        normalized_vector = [vector[i]/normal for i in range(dimension)]
    return normalized_vector


def get_contact_information(frag_1_coords, frag_2_coords, fst_fragment_atoms, snd_fragment_atoms, elements):
    contact = False
    for i in range(len(frag_1_coords)):
        for j in range(len(frag_2_coords)):
            Ai = fst_fragment_atoms[i].lower()
            Bj = snd_fragment_atoms[j].lower()
            R_AB = elements[Ai] + elements[Bj]
            dxb = frag_1_coords[i][0] - frag_2_coords[j][0]
            dyb = frag_1_coords[i][1] - frag_2_coords[j][1]
            dzb = frag_1_coords[i][2] - frag_2_coords[j][2]
            r2 = dxb*dxb + dyb*dyb + dzb*dzb
            r = math.sqrt(r2)
            if r <= (R_AB + 0.5):
                contact = True
                break
            else:
                continue
    return contact


def displace_coordinates_along_vector(COM, frag_coords, vector, delta_displacement):
    new_COM = []
    for i in range(len(COM)):
        COM[i] = COM[i] + vector[i] * delta_displacement
    new_COM.append(COM[:])
    modified_frag_coords = set_origin(frag_coords[:], displacement=new_COM)
    return modified_frag_coords


#def translate(atom_name, frag_1_coords, frag_2_coords, translational_vector):
#    """This function will take two fragments coordinates as input and
#    modify the unfavorable interaction if present"""
#    from pyreactor.utils.globals import ELEMENTS as elements
#    delta_displacement = 0.1
#    fst_fragment_atoms = atom_name[:len(frag_1_coords)]
#    snd_fragment_atoms = atom_name[len(frag_1_coords):(len(frag_1_coords) + len(frag_2_coords))]
#    COM_frag_1 = centroid(frag_1_coords)
#    COM_frag_2 = centroid(frag_2_coords)
#    vector = get_vector(COM_frag_1, translational_vector, 12)
#    normalize_vector = get_normalize_vector(vector)
#    in_contact = get_contact_information(frag_1_coords, frag_2_coords, fst_fragment_atoms, snd_fragment_atoms, elements)
#    while in_contact:
#        displaced_frag_2_coordinates = displace_coordinates_along_vector(COM_frag_2, frag_2_coords, normalize_vector, delta_displacement)
#        frag_2_coords = displaced_frag_2_coordinates[:]
#        in_contact = get_contact_information(frag_1_coords, frag_2_coords, fst_fragment_atoms, snd_fragment_atoms, elements)
#    return frag_2_coords[:]


def get_molecule_com_distances(coordinates, com):
    """This function will return a list of distances between
       a coordinate set and a point (COM) """
    distances_list = []
    for i in coordinates:
        dx = i[0] - com[0]
        dy = i[1] - com[1]
        dz = i[2] - com[2]
        r = math.sqrt(dx*dx + dy*dy + dz*dz)
        distances_list.append(r)
    return distances_list


def displace_com_along_vector(com, vector, delta_displacement=0.1):
    """It will displace a given point (COM) for dr along a vector. """
    new_com = []
    for i in range(len(com)):
        com[i] += vector[i] * delta_displacement
        new_com.append(com[i])
    return new_com


def check_decrement(distance_list1, distance_list2):
    """This function will check if in distances between the COM of a fragment are
       decreasing or increasing with all the atoms of fragment 2. If with one
       small displacement of the COM of fragment 1 along a vector leads to the
       decreament of the distance with atleast one atom of fragment 2, then,
       the possibility is that the molecule is in the cage."""
    caged = False
    if len(distance_list1) == len(distance_list2):
        for i in range(len(distance_list1)):
            if distance_list1[i] > distance_list2[i]:
                caged = True
                break
    else:
        print ">>>>WARNING<<<"
        print "The two lists:", distance_list1, "and", distance_list2, "are not same."
        print "Please check carefully"
    return caged


def translate(atom_name, frag_1_coords, frag_2_coords, translational_vector):
    """This function will take two fragments coordinates as input and translate.
       This is applicable for caged molecules also."""
    from PyAR.utils.globals import ELEMENTS as elements
    delta_displacement = 0.1
    fst_fragment_atoms = atom_name[:len(frag_1_coords)]
    snd_fragment_atoms = atom_name[len(frag_1_coords):(len(frag_1_coords) + len(frag_2_coords))]
    COM_frag_1 = centroid(frag_1_coords)
    COM_frag_2 = centroid(frag_2_coords)
    vector = get_vector(COM_frag_1, translational_vector, 12)
    normalize_vector = get_normalize_vector(vector)
    while True:
        frag1_com2_distances1 = get_molecule_com_distances(frag_1_coords, COM_frag_2)
        transfered_com2 = displace_com_along_vector(COM_frag_2, vector, delta_displacement)
        frag1_com2_distances2 = get_molecule_com_distances(frag_1_coords, transfered_com2)
        caged_test = check_decrement(frag1_com2_distances1[:], frag1_com2_distances2[:])
        if caged_test:
            COM_frag_2 = transfered_com2[:]
        else:
            break
    new_frag_2_coords = set_origin(frag_2_coords[:], displacement=COM_frag_2)
    in_contact = get_contact_information(frag_1_coords, new_frag_2_coords, fst_fragment_atoms, snd_fragment_atoms, elements)
    while in_contact:
        displaced_frag_2_coordinates = displace_coordinates_along_vector(COM_frag_2, new_frag_2_coords, normalize_vector, delta_displacement)
        new_frag_2_coords = displaced_frag_2_coordinates[:]
        in_contact = get_contact_information(frag_1_coords[:], new_frag_2_coords[:], fst_fragment_atoms, snd_fragment_atoms, elements)
    return new_frag_2_coords[:]


def main():
    pass


if __name__ == "__main__":
    main()

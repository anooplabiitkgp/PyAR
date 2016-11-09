'''
restraint.py - implementation of force using Maeda and Morokuma
JCTC 2011,7,2335 in PyAR
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
from math import sqrt
from sets import Set
import operator,itertools

# covalent radii (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009, 188-197)
# values for metals decreased by 10 %
elements = { 'x ': 0.00,  \
             'h' : 0.32, 'he': 0.46, 'li': 1.20, 'be': 0.94, 'b' : 0.77, \
             'c' : 0.75, 'n' : 0.71, 'o' : 0.63, 'f' : 0.64, 'ne': 0.67, \
             'na': 1.40, 'mg': 1.25, 'al': 1.13, 'si': 1.04, 'p' : 1.10, \
             's' : 1.02, 'cl': 0.99, 'ar': 0.96, 'k' : 1.76, 'ca': 1.54, \
             'sc': 1.33, 'ti': 1.22, 'v' : 1.21, 'cr': 1.10, 'mn': 1.07, \
             'fe': 1.04, 'co': 1.00, 'ni': 0.99, 'cu': 1.01, 'zn': 1.09, \
             'ga': 1.12, 'ge': 1.09, 'as': 1.15, 'se': 1.10, 'br': 1.14, \
             'kr': 1.17, 'rb': 1.89, 'sr': 1.67, 'y' : 1.47, 'zr': 1.39, \
             'nb': 1.32, 'mo': 1.24, 'tc': 1.15, 'ru': 1.13, 'rh': 1.13, \
             'pd': 1.08, 'ag': 1.15, 'cd': 1.23, 'in': 1.28, 'sn': 1.26, \
             'sb': 1.26, 'te': 1.23, 'i' : 1.32, 'xe': 1.31, 'cs': 2.09, \
             'ba': 1.76, 'la': 1.62, 'ce': 1.47, 'pr': 1.58, 'nd': 1.57, \
             'pm': 1.56, 'sm': 1.55, 'eu': 1.51, 'gd': 1.52, 'tb': 1.51, \
             'dy': 1.50, 'ho': 1.49, 'er': 1.49, 'tm': 1.48, 'yb': 1.53, \
             'lu': 1.46, 'hf': 1.37, 'ta': 1.31, 'w' : 1.23, 're': 1.18, \
             'os': 1.16, 'ir': 1.11, 'pt': 1.12, 'au': 1.13, 'hg': 1.32, \
             'tl': 1.30, 'pb': 1.30, 'bi': 1.36, 'po': 1.31, 'at': 1.38, \
             'rn': 1.42, 'fr': 2.01, 'ra': 1.81, 'ac': 1.67, 'th': 1.58, \
             'pa': 1.52, 'u' : 1.53, 'np': 1.54, 'pu': 1.55 }


def read_xyz(filename):
    ''' This is a function for reading xyz files'''
    f = open(filename).readlines()
    nat = int(f[0])
    title = f[1].rstrip()
    raw = [f[i].split() for i in range(2,nat+2)]
    xyz = [(raw[i][0], float(raw[i][1]), float(raw[i][2]), float(raw[i][3])) for i in range(nat)]
    return nat, title, xyz


def read_turbomole_coord():
    '''returns number of atoms, and coordinates as x, y, z, sym'''
    f = open('coord').read().split('$')[1].split('\n')[1:-1]
    p = open('coord').read()
    coords = [i.split() for i in f]
    coords = [(float(i[0]),float(i[1]),float(i[2]),i[3]) for i in coords]
    return len(coords), coords


def rewrite_gradients(number_of_atoms, xyz, total_restraint_energy, total_restraint_gradients, restraint_gradient):
    '''Rewrite new gradient file in turbomole format'''
    import re
    with open('gradient','r') as f:
        read_file = f.read()
        gradients = re.split("cycle",read_file)
        prelog = "cycle".join(gradients[:-1])
        last_gradient = gradients[-1]
        list_from_last_gradient = last_gradient.split('\n')

        cycle_number = int(list_from_last_gradient[0].split()[1])
        scf_energy = float(list_from_last_gradient[0].split()[5])
        total_gradients = float(list_from_last_gradient[0].split()[8])
        new_energy = scf_energy + total_restraint_energy
        new_total_gradients = total_gradients + total_restraint_gradients
        print "  cycle =%7d    SCF energy =%20.10f   |dE/dxyz| =%10.6f " \
              %(cycle_number , new_energy, new_total_gradients)

        for line in list_from_last_gradient[1:number_of_atoms]: print line

        for line in list_from_last_gradient[number_of_atoms:number_of_atoms*2+1]:
            dx, dy,dz, Z = line.split()[0], line.split()[1],line.split()[2], line.split()[3]
            print Z,re.sub('E','D',dx)

        for line in range(number_of_atoms): print "  %.14E" % restraint_gradients[line][0]
        #print "  cycle"+last_gradient


def read_fragments(n):
    try:
        fp= open('./fragment').readlines()
    except:
        print "./fragment file should exist!"
        exit()
    lines = [lines for lines in fp if lines.strip()]
    nMol=len(lines)
    atoms_in_fragments=[]

    for line in lines:
        for delim in ',;':
            line = line.replace(delim, ' ')
        a = []
        for things in line.split():
            if "-" in things:
                start,end = map(int,things.split('-'))
                a = a + [i for i in range(start,end+1)]
            else: a = a + [int(things)]
        atoms_in_fragments.append(a)
    collected_atoms = Set([item for sublist in atoms_in_fragments for item in sublist])
    if len(collected_atoms) < n:
        all_atoms_in_molecules = Set([i for i in range(1,n+1)])
        remaining = list(all_atoms_in_molecules - collected_atoms)
        
        atoms_in_fragments.append(remaining)
    common = reduce(operator.iand, map(set,atoms_in_fragments))
    if common:
        print "Fragments contain common atoms"
        exit()
    n_fragments=len(atoms_in_fragments)
    return n_fragments,atoms_in_fragments


def usage():
    print "Usage: ", sys.argv[0], " Gamma"
    print "       Gamma is model collition energy (in kJ/mol)"


def main():
    bohr2angstrom  = 0.52917726
    atomicunit2kjoules   = 2625.5
    atomicunit2kcals = 627.509541
    
    try:
        gamma=float(sys.argv[1])/atomicunit2kjoules
    except:
        print "No Gamma value given"
        usage()
        exit()
    
    number_of_atoms, xyz = read_turbomole_coord()
    number_of_fragments,atoms_in_fragment = read_fragments(number_of_atoms)
    
    epsilon=1.0061/atomicunit2kjoules
    r_zero=3.8164/bohr2angstrom
    
    print "Gamma / a.u.: ", gamma
    
    
    #    eqn. 3 JCTC 2011,7,2335
    alpha = gamma/((2**(-1.0/6.0)-(1+sqrt(1+gamma/epsilon))**(-1.0/6.0))*r_zero)
    print "Alpha / a.u.: ", alpha
    
    # idw = inverse distance weighting factor
    idw_parameter = 6.0
    
    elements['h'] = 0.00001
    
    #    eqn. 2 JCTC 2011,7,2335
    # F(Q) = E(Q) + Alpha * factor
    
    idw_nominator=0.0
    idw_denominator=0.0
    for i in range(number_of_fragments):
        for j in range(number_of_fragments):
            if i < j:
                for k in range(len(atoms_in_fragment[i])):
                    for l in range(len(atoms_in_fragment[j])):
                        Ai = atoms_in_fragment[i][k]-1
                        Bj = atoms_in_fragment[j][l]-1
                        dxb=xyz[Ai][0]-xyz[Bj][0]
                        dyb=xyz[Ai][1]-xyz[Bj][1]
                        dzb=xyz[Ai][2]-xyz[Bj][2]
                        r2=dxb*dxb+dyb*dyb+dzb*dzb
                        r=sqrt(r2)
                        R_AB = elements[xyz[Ai][3]] + elements[xyz[Bj][3]]
                        idw_nominator=idw_nominator   + ((R_AB/r)**idw_parameter)*r
                        idw_denominator=idw_denominator + (R_AB/r)**idw_parameter
    
    fac=idw_nominator/idw_denominator
    print 'Inverse Weighting Factor =',fac
    total_restraint_energy = fac*alpha
    print '  Energy term from AFIR =',total_restraint_energy,' a.u.'
    print    '=',total_restraint_energy*atomicunit2kcals,' kcal/mol'
    
    v=idw_nominator
    w=idw_denominator
    w2=idw_denominator*idw_denominator
    
    R_AB=0.0
    
    restraint_gradients = [[0.0 for j in range(3)] for i in range(number_of_atoms)]
     
    for i in range(number_of_fragments):
        for j in range(number_of_fragments):
            if i < j:
                for k in range(len(atoms_in_fragment[i])):
                    for l in range(len(atoms_in_fragment[j])):
                        Ai = atoms_in_fragment[i][k]-1
                        Bj = atoms_in_fragment[j][l]-1
                        dxb=xyz[Ai][0]-xyz[Bj][0]
                        dyb=xyz[Ai][1]-xyz[Bj][1]
                        dzb=xyz[Ai][2]-xyz[Bj][2]
                        r2=dxb*dxb+dyb*dyb+dzb*dzb
                        r=sqrt(r2)
                        R_AB = elements[xyz[Ai][3]] + elements[xyz[Bj][3]]
                    
                        term1=alpha*((1-idw_parameter)*R_AB**idw_parameter*r2**((-1-idw_parameter)/2))/w
                        term2=alpha*(idw_parameter*v*R_AB**idw_parameter*r2**((-2-idw_parameter)/2))/w2
                    
                        restraint_gradients[Ai][0] = restraint_gradients[Ai][0] -dxb*term1 -dxb*term2
                        restraint_gradients[Ai][1] = restraint_gradients[Ai][1] -dyb*term1 -dyb*term2
                        restraint_gradients[Ai][2] = restraint_gradients[Ai][2] -dzb*term1 -dzb*term2
    
    #ger = -ger
    print 'AFIR Gradients are:'
    total_restraint_gradients = 0.0
    for i in range(number_of_atoms):
        total_restraint_gradients = total_restraint_gradients+restraint_gradients[i][0]+restraint_gradients[i][1]+restraint_gradients[i][2]
        print restraint_gradients[i]
    print total_restraint_gradients
    print 'afir all done'
    rewrite_gradients(number_of_atoms,xyz,total_restraint_energy,total_restraint_gradients,restraint_gradients)


if __name__ == '__main__':
    main()

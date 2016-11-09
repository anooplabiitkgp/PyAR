'''
clusterloop.py - to generate aggregates using PyAR
'''
'''
Copyright (C) 2016 by Surajit Nandi, Anoop Ayyappan, and Mark P. Waller
Indian Institute of Technology Kharagpur, India and Westfaelische Wilhelms
Universitaet Muenster, Germany

This file is part of the PyAR project.

PyAR is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
'''
from math import *
import operator
import sys, os.path
from numpy import linalg, ndarray
import pprint


def find_centre(coords):
    nat = len(coords)
    xcentre = sum([coords[i][1] for i in range(nat)])/nat
    ycentre = sum([coords[i][2] for i in range(nat)])/nat
    zcentre = sum([coords[i][3] for i in range(nat)])/nat
    centre = xcentre, ycentre, zcentre
    return centre


def translate(coords, centre):
    nat=len(coords)
    translated = [(coords[i][0], coords[i][1] + centre[0], coords[i][2] + centre[1], coords[i][3] + centre[2])
                  for i in range(nat)]
    return translated


def to_centre(coords):
    centre = find_centre(coords)
    centered = translate(coords, centre)
    return centered


def centre(coords1,coords2):
    nat = len(coords1)
    centre = find_centre(coords1)
    trans = (
              sum([coords1[i][1] - coords2[i][1] for i in range(nat)])/nat,
              sum([coords1[i][2] - coords2[i][2] for i in range(nat)])/nat,
              sum([coords1[i][3] - coords2[i][3] for i in range(nat)])/nat,
              )
    centered = translate(coords2, trans)
    return centered


def getrotation (ref, tgt):
    R11 = 0.0
    R12 = 0.0
    R13 = 0.0
    R21 = 0.0
    R22 = 0.0
    R23 = 0.0
    R31 = 0.0
    R32 = 0.0
    R33 = 0.0
    nat = len(ref)

    for i in range(nat):
        R11 += ref[i][1] * tgt[i][1]
        R12 += ref[i][1] * tgt[i][2]
        R13 += ref[i][1] * tgt[i][3]
        R21 += ref[i][2] * tgt[i][1]
        R22 += ref[i][2] * tgt[i][2]
        R23 += ref[i][2] * tgt[i][3]
        R31 += ref[i][3] * tgt[i][1]
        R32 += ref[i][3] * tgt[i][2]
        R33 += ref[i][3] * tgt[i][3]

    c = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]

    c[0][0] = R11 + R22 + R33
    c[0][1] = R23 - R32
    c[0][2] = R31 - R13
    c[0][3] = R12 - R21
    
    c[1][1] = R11 - R22 - R33
    c[1][2] = R12 + R21
    c[1][3] = R31 + R13
    
    c[2][2] = R22 - R33 - R11
    c[2][3] = R23 + R32
    c[3][3] = R33 - R11 - R22

    d,v = np.linalg.eigh(c)
    # d, v = jacobi(c)
    Qn = [v[0][3], v[1][3], v[2][3], v[3][3]]
    R = Q2R (Qn)
    return Qn, R


def rotate(coords, R):
    """R = rotation matrix"""
    nat = len(coords)
    rotated =[(coords[i][0],
        R[0][0] * coords[i][1] + R[0][1] * coords[i][2] + R[0][2] * coords[i][3],
        R[1][0] * coords[i][1] + R[1][1] * coords[i][2] + R[1][2] * coords[i][3],
        R[2][0] * coords[i][1] + R[2][1] * coords[i][2] + R[2][2] * coords[i][3]
        ) for i in range(nat)]
    return rotated


def Q2R(Qn):
    """ Generate a rotation matrix from a normalized quaternion"""
    R = (
    (
    (Qn[0]*Qn[0] + Qn[1]*Qn[1] - Qn[2]*Qn[2] - Qn[3]*Qn[3]),
    (2.0 * (Qn[1] * Qn[2] - Qn[0] * Qn[3])),
    (2.0 * (Qn[1] * Qn[3] + Qn[0] * Qn[2])),
    ),
    (
    (2.0 * (Qn[2] * Qn[1] + Qn[0] * Qn[3])),
    (Qn[0]*Qn[0] - Qn[1]*Qn[1] + Qn[2]*Qn[2] - Qn[3]*Qn[3]),
    (2.0 * (Qn[2] * Qn[3] - Qn[0] * Qn[1])),
    ),
    (
    (2.0 *(Qn[3] * Qn[1] - Qn[0] * Qn[2])),
    (2.0 * (Qn[3] * Qn[2] + Qn[0] * Qn[1])),
    (Qn[0]*Qn[0] - Qn[1]*Qn[1] - Qn[2]*Qn[2] + Qn[3]*Qn[3])
    )
    )
    return R


def jacobi (a):
    nrot = 1000
    v = [ [1,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,1] ]
    d = [a[j][j] for j in range(4)]

    for l in range(nrot):
        dnorm = 0.0; onorm = 0.0
        for j in range(4):
            dnorm += fabs(d[j])
            for i in range(j): onorm += fabs(a[i][j]);
        if ((onorm/dnorm) <= 1.0e-12):
            nrot = l
            for j in range(3):
                k = j
                dtemp = d[k]
                for i in range(j+1,4):
                   if(d[i] < dtemp):
                       k = i
                       dtemp = d[k]
                if(k > j):
                    d[k] = d[j]
                    d[j] = dtemp
                    for i in range(4):
                        dtemp = v[i][k]
                        v[i][k] = v[i][j]
                        v[i][j] = dtemp
            break
        else:
            for j in range(1,4):
                for i in range(j):
                    b = a[i][j]
                    if(fabs(b) > 0.0):
                         dma = d[j] - d[i]
                         if((fabs(dma) + fabs(b)) <=  fabs(dma)):
                             t = b / dma
                         else:
                             q = 0.5 * dma / b
                             t = 1.0/(fabs(q) + sqrt(1.0+q*q))
                             if(q < 0.0):
                                 t = -t
                         c = 1.0/sqrt(t * t + 1.0)
                         s = t * c
                         a[i][j] = 0.0
                         for k in range(i):
                             atemp = c * a[k][i] - s * a[k][j]
                             a[k][j] = s * a[k][i] + c * a[k][j]
                             a[k][i] = atemp
                         for k in range(i+1,j):
                             atemp = c * a[i][k] - s * a[k][j]
                             a[k][j] = s * a[i][k] + c * a[k][j]
                             a[i][k] = atemp
                         for k in range(j+1,4):
                             atemp = c * a[i][k] - s * a[j][k]
                             a[j][k] = s * a[i][k] + c * a[j][k]
                             a[i][k] = atemp
                         for k in range(4):
                             vtemp = c * v[k][i] - s * v[k][j]
                             v[k][j] = s * v[k][i] + c * v[k][j]
                             v[k][i] = vtemp
                         dtemp = c*c*d[i] + s*s*d[j] - 2.0*c*s*b
                         d[j] = s*s*d[i] + c*c*d[j] +  2.0*c*s*b
                         d[i] = dtemp
    return d, v


def fit(coords1,coords2):
    coords2 = centre(coords1,coords2)
    Qr, R = getrotation(coords2,coords1)
    coords2 = rotate(coords2, R)
    coords2 = centre(coords1,coords2)
    return coords2


def getRMSD(coords1,coords2):
    fitted = fit(coords1,coords2)
    rms = 0.0
    nat = len(coords1)
    for i in range(nat):
        dx = coords1[i][1] - fitted[i][1]
        dy = coords1[i][2] - fitted[i][2]
        dz = coords1[i][3] - fitted[i][3]
        rms += (dx*dx + dy*dy + dz*dz)
    rmsd = sqrt(rms/nat)
    return rmsd


## following five functions for kabsch rmsd
def krmsd(V, W):
    """
    Calculate Root-mean-square deviation from two sets of vectors V and W.
    This is the same as rmsd written above but with different format so,
    the function was copied from kabsch rmsd program
    """
    D = len(V[0])
    N = len(V)
    rmsd = 0.0
    for v, w in zip(V, W):
        rmsd += sum([(v[i]-w[i])**2.0 for i in range(D)])
    return np.sqrt(rmsd/N)


def kabsch_rotate(P, Q):
    """
    Rotate matrix P unto matrix Q using Kabsch algorithm
    """
    U = kabsch(P, Q)

    # Rotate P
    P = np.dot(P, U)
    return P


def kabsch(P, Q):
    """
    The optimal rotation matrix U is calculated and then used to rotate matrix
    P unto matrix Q so the minimum root-mean-square deviation (RMSD) can be
    calculated.

    Using the Kabsch algorithm with two sets of paired point P and Q,
    centered around the center-of-mass.
    Each vector set is represented as an NxD matrix, where D is the
    the dimension of the space.

    The algorithm works in three steps:
    - a translation of P and Q
    - the computation of a covariance matrix C
    - computation of the optimal rotation matrix U

    http://en.wikipedia.org/wiki/Kabsch_algorithm

    Parameters:
    P -- (N, number of points)x(D, dimension) matrix
    Q -- (N, number of points)x(D, dimension) matrix

    Returns:
    U -- Rotation matrix

    """

    # Computation of the covariance matrix
    C = np.dot(np.transpose(P), Q)

    # Computation of the optimal rotation matrix
    # This can be done using singular value decomposition (SVD)
    # Getting the sign of the det(V)*(W) to decide
    # whether we need to correct our rotation matrix to ensure a
    # right-handed coordinate system.
    # And finally calculating the optimal rotation matrix U
    # see http://en.wikipedia.org/wiki/Kabsch_algorithm
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0
    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]

    # Create Rotation matrix U
    U = np.dot(V, W)

    return U


def centroid(X):
    """
    Calculate the centroid from a vectorset X
    """
    Y = []
    for i in X:
        if len(i) == 3:
            Y = X[:]
            break
        elif len(i) == 4:
            Y.append(i[1:])
        else:
            raise "Invalid xyz format", X
    C = sum(np.array(Y))/len(np.array(Y))
    return C


def KabschRMSD(P, Q):
    """
    Rotate matrix P unto Q and calculate the RMSD
    """
    P1 = []
    Q1 = []
    P = to_centre(P)
    Q = to_centre(Q)
    for i in P:
        if len(i) == 3:
            P1 = P
            break
        elif len(i) == 4:
            P1.append(i[1:])
        else:
            raise ValueError("Wrong format P")
    for i in Q:
        if len(i) == 3:
            Q1 = Q
            break
        elif len(i) == 4:
            Q1.append(i[1:])
        else:
            raise ValueError("Wrong format Q")
    P1 = np.array(P1[:])
    Q1 = np.array(Q1[:])
    P1 = kabsch_rotate(P1, Q1)
    return krmsd(P1, Q1)


def getEnergy(mopac_arc_file):
    """
    Get heat of formation from mopac arc file
    """
    with open(mopac_arc_file) as arc_out:
        arc_contents = arc_out.readlines()
    for line in arc_contents:
        if 'Empirical Formula:' in line:
            natoms = int(line.split()[-2])
        if 'HEAT OF FORMATION' in line:
            energy = float(line.split()[4])
    return energy


def getCoordinates(mopac_arc_file):
    """
    Get coordinates from mopac arc file
    """
    with open(mopac_arc_file) as arc_out:
        arc_contents = arc_out.readlines()
    for line in arc_contents:
        if 'Empirical Formula:' in line:
            natoms = int(line.split()[-2])
    raw = [i.split() for i in arc_contents[-(natoms+1):-1]]
    xyz = [(raw[i][0], float(raw[i][1]), float(raw[i][3]), float(raw[i][5])) for i in range(natoms)]
    return xyz


def getLowest(dict_name):
    """
    :param dict_name: This is the dictionary which contains filename and value. 
    :return:[filename, value] of the lowest
    """
    return sorted(dict_name.items(), key=operator.itemgetter(1))[0]


def getWorst(dict_name):
    """
    :param dict_name: This is the dictionary which contains filename and value. 
    :return:[filename,value] of the higest
    """
    return sorted(dict_name.items(), key=operator.itemgetter(1))[-1]


def getDistance(a,b):
    dx = a[1] - b[1]
    dy = a[2] - b[2]
    dz = a[3] - b[3]
    r = sqrt(dx**2+dy**2+dz**2)
    return r


def getAverageRadii(coords):
    natoms = len(coords)
    x, y, z = find_centre(coords)
    return sum([getDistance(['x',x,y,z],coords[i]) for i in range(natoms)])/natoms


def getSortedEnergyList(fileList):
    energies = {i: getEnergy(i) for i in fileList}
    return sorted(energies.items(), key=operator.itemgetter(1))


def removeJuncCoordsFiles(fileList):
    modifiedFileList = []
    for i in fileList:
        with open(i) as file_out:
            arc_contents = file_out.readlines()
        for line in arc_contents:
            if 'Empirical Formula:' in line:
                natoms = int(line.split()[-2])
        raw = [j.split() for j in arc_contents[-(natoms+1):-1]]
        for j in range(natoms):
            try:
                x = float(raw[j][1])
                y = float(raw[j][1])
                z = float(raw[j][1])
                xyz = [(raw[j][0], x, y, z)]
            except:
                print "File:", i, "has no meaningfull coordinates"
                warnings.warn("string not convertible to float. neglecting the structure.",RuntimeWarning)
                xyz = []
                break
        if xyz:
            modifiedFileList.append(i)
        else:
            print "File:", i,"coordinate data is empty."
            print "File:", i,"will be removed from the list."
    return modifiedFileList


def makeClusterRMSD(fileList):
    sortedEnergyList = getSortedEnergyList(fileList)
    lowestEnergyMolecule = sortedEnergyList[0][0]
    RMSDwrtLowestEnergyMolecule = {i: getRMSD(getCoordinates(lowestEnergyMolecule),getCoordinates(i)) for i in fileList}
    farthestPoint = sorted(RMSDwrtLowestEnergyMolecule.items(), key=operator.itemgetter(1))[-1][0]
    
    RMSDwrtFarthestPoint = {i: getRMSD(getCoordinates(farthestPoint),getCoordinates(i)) for i in fileList}

    #Clustering into two groups, one close to lowest energy structure, and the other close to farthest point
    clusterOne = []; clusterTwo = []
    for i in fileList:
        clusterOne.append(i) if RMSDwrtLowestEnergyMolecule[i] < RMSDwrtFarthestPoint[i] else clusterTwo.append(i)
    clusterOneTable = {i: getEnergy(i) for i in clusterOne}
    clusterTwoTable = {i: getEnergy(i) for i in clusterTwo}
    return clusterOneTable, clusterTwoTable


def makeClusterKabshRMSD(fileList):
    sortedEnergyList = getSortedEnergyList(fileList)
    lowestEnergyMolecule = sortedEnergyList[0][0]
    RMSDwrtLowestEnergyMolecule = {i: KabschRMSD(getCoordinates(lowestEnergyMolecule),getCoordinates(i)) for i in fileList}
    farthestPoint = sorted(RMSDwrtLowestEnergyMolecule.items(), key=operator.itemgetter(1))[-1][0]
    RMSDwrtFarthestPoint = {i: KabschRMSD(getCoordinates(farthestPoint),getCoordinates(i)) for i in fileList}
    #Clustering into two groups, one close to lowest energy structure, and the other close to farthest point
    clusterOne = []; clusterTwo = []
    for i in fileList:
        clusterOne.append(i) if RMSDwrtLowestEnergyMolecule[i] < RMSDwrtFarthestPoint[i] else clusterTwo.append(i)
    clusterOneTable = {i: getEnergy(i) for i in clusterOne}
    clusterTwoTable = {i: getEnergy(i) for i in clusterTwo}
    return clusterOneTable, clusterTwoTable


def makeClusterAvgRad(fileList):

    averageRadii = {i: getAverageRadii(getCoordinates(i)) for i in fileList}
    smallestMolecule, smallest = sorted(averageRadii.items(), key=operator.itemgetter(1))[0]
    largestMolecule, largest = sorted(averageRadii.items(), key=operator.itemgetter(1))[-1]

    a = sorted(averageRadii.items(), key=operator.itemgetter(1))
    maximumGap = max([x[1] - a[i-1][1] for i, x in enumerate(a)][1:])

    clusterOne = []; clusterTwo = []
    for i in fileList:
        similarityWithSmallestMolecule  = abs(averageRadii[smallestMolecule]-averageRadii[i])
        similarityWithLargestMolecule = abs(averageRadii[largestMolecule]-averageRadii[i])
        clusterOne.append(i) if similarityWithSmallestMolecule <= similarityWithLargestMolecule else clusterTwo.append(i)
    clusterOneTable = {i: averageRadii[i] for i in clusterOne}
    clusterTwoTable = {i: averageRadii[i] for i in clusterTwo}
    print 'Clustering with Average Radii, Width(1) =', (getWorst(clusterOneTable)[1]-getLowest(clusterOneTable)[1])
    print 'Clustering with Average Radii, Width(2) =', (getWorst(clusterTwoTable)[1]-getLowest(clusterTwoTable)[1])
    print 'Cluster Gap       =', (getLowest(clusterTwoTable)[1]-getWorst(clusterOneTable)[1])

    return clusterOneTable, clusterTwoTable


def makeClusterEnergy(fileList):

    energyList = getSortedEnergyList(fileList)
    lowestEnergyMolecule, lowestEnergy   = energyList[0]
    highestEnergyMolecule, highestEnergy = energyList[-1]

    clusterOne = []; clusterTwo = []
    for i in fileList:
        similarityWithLowestEnergyMolecule  = abs(lowestEnergy-getEnergy(i))
        similarityWithHighestEnergyMolecule = abs(highestEnergy-getEnergy(i))
        clusterOne.append(i) if similarityWithLowestEnergyMolecule <= similarityWithHighestEnergyMolecule else clusterTwo.append(i)
    clusterOneTable = {i: getEnergy(i) for i in clusterOne}
    clusterTwoTable = {i: getEnergy(i) for i in clusterTwo}
    print 'energy width of first cluster   =', (getSortedEnergyList(clusterOne)[-1][1]-getSortedEnergyList(clusterOne)[0][1])
    print 'energy width of second cluster  =', (getSortedEnergyList(clusterTwo)[-1][1]-getSortedEnergyList(clusterTwo)[0][1])
    print 'energy gap between two clusters =', (getSortedEnergyList(clusterTwo)[0][1] -getSortedEnergyList(clusterOne)[-1][1])

    return clusterOneTable, clusterTwoTable


def quadsplit(filelist):
    A = []; B = []; a = []; b = []; c = []; d = []; clusters = []
    A, B = makeClusterRMSD(filelist)
    if len(A) > 1:
        a, b = makeClusterRMSD(A)
    if len(B) > 1:
        c, d = makeClusterRMSD(B)
    clusters = [a, b, c, d]
    print 'clustering based on similarity'
    print 'cluster size:', [len(i) for i in clusters]
    bestofeach = [getLowest(i)[0] for i in clusters if i]
    report(bestofeach)
    return bestofeach


def report(filelist):
    energyList = getSortedEnergyList(filelist)
    longestname = max([len(i[0]) for i in energyList])
    print "-"*(longestname+23)
    print 'filename'+' '*(longestname-5)+' r.e. ,avg.rad  '
    print "-"*(longestname+23)
    for i in energyList:
         relativeEnergy = i[1]-energyList[0][1]
         averageRadius = getAverageRadii(getCoordinates(i[0]))
         print "%s %s%7.2f%7.3f" %(i[0], "_"*(1+longestname-len(i[0])), relativeEnergy, averageRadius)
    print "-"*(longestname+23)


def bestCluster(arcfiles):
    fileList =  arcfiles
    radiiList = {i: getAverageRadii(getCoordinates(i)) for i in fileList}
    energyList = {i: getEnergy(i) for i in fileList}

    uniqueList = {}
    factor = 2
    radiiThreshold = 100.0
    energyThreshold = 100.0
    while len(uniqueList) <4 and energyThreshold > 0.1 and radiiThreshold > 0.001:
        radiiThreshold = (max(radiiList.values()) - min(radiiList.values()))/factor
        energyThreshold = (max(energyList.values()) - min(energyList.values()))/factor
        print factor, radiiThreshold, energyThreshold,
     
        uniqueList = fileList[:]
        for i in fileList:
            for j in fileList:
                if i != j and abs(radiiList[i]-radiiList[j]) < radiiThreshold and abs(energyList[i]-energyList[j]) < energyThreshold:
                    if energyList[i] <= energyList[j]:
                        if j in uniqueList:
                            uniqueList.remove(j)
                    else:
                        if i in uniqueList:
                            uniqueList.remove(i)
        print len(uniqueList)
        factor += 1
        if uniqueList: report(uniqueList)
    if len(uniqueList) > 4:
        finalList = [x[0] for x in getSortedEnergyList(uniqueList)][:4]
    else:
        finalList = uniqueList
    print finalList
    return finalList


# main program
def main():
    try:
        arcfiles = sys.argv[1:]
    except:
        print 'usage: ', sys.argv[0], 'mopac arc files'
        sys.exit(1)
    if len(arcfiles) < 2:
       print "Not enough files to cluster"
       exit(0)
    betCluster(arcfiles)


if __name__ == "__main__":
    main()

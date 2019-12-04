# Imports
from Angle_Calc import *
from sympy import Point3D
import scipy
import prody
import fileinput
import pandas as pd
import numpy as np
from BondDict import aa123
import sys
from itertools import islice

import shutil
shutil.copyfile('Start_peptide2.pdb', 'Start_peptide2 - Copy.pdb')

np.set_printoptions(formatter={'float': '{: 0.3f}'.format})

for line in fileinput.input("Start_peptide2 - Copy.pdb", inplace=1):
    if fileinput.filelineno() <= 4:
        line = line.replace('ASP', aa123[sys.argv[1]])
    if fileinput.filelineno() > 4:
        line = line.replace('PHE', aa123[sys.argv[2]])
    print(line.rstrip('\n'))

startpep = pd.read_table(
    "Start_peptide2 - Copy.pdb",
    delim_whitespace=True,
    header=None)

# print(getAngle(p0=np.array(startpep.iloc[1, 6:9]),
#                p1=np.array(startpep.iloc[2, 6:9]),
#                p2=np.array(startpep.iloc[3, 6:9])))

startpepPDB = prody.parsePDB("Start_peptide2 - Copy.pdb")

for p in range(1, 400):
    startpep = pd.read_table(
        "Start_peptide2 - Copy.pdb",
        delim_whitespace=True,
        header=None)
    with open("Outputfile2.pdb", "a+") as outputfile:
        with open("Start_peptide2 - Copy.pdb", 'r') as myfile:
            lines = myfile.readlines()
            outputfile.writelines(lines[0:7])
        PhiPsi_All = np.zeros((len(sys.argv) - 1, 2))
        for i in range(1, 10):
            if startpep.iloc[-1, 5] == 8:
                break
            print(PhiPsi_All)
            with open('Start_peptide2 - Copy.pdb', 'a') as myfile:
                print(startpep)
                if startpep.iloc[-1, 2] == "CA":
                    PhiPsi = PhiPsi_All[int(startpep.iloc[-1, 5])]
                    x1 = calculateCoordinates(Point3D(startpep.iloc[-4, 6:9]),
                                              Point3D(startpep.iloc[-2, 6:9]),
                                              Point3D(startpep.iloc[-1, 6:9]),
                                              1.526,
                                              cdl_spec(sys.argv[int(startpep.iloc[-5, 5])], sys.argv[int(startpep.iloc[-1, 5])],
                                                       "NAC", PhiPsi[0], PhiPsi[1]),
                                              PhiPsi[1])
                    var1 = np.around(x1, decimals=3)
                    print(var1)
                    new_row = [None]*11
                    new_row[0] = "ATOM".ljust(7)
                    new_row[1] = str(int(startpep.iloc[-1, 1]) + 1).rjust(4)
                    new_row[2] = "C".ljust(4)
                    new_row[3] = str(startpep.iloc[-1, 3]).ljust(4)
                    new_row[4] = "B".ljust(4)
                    new_row[5] = str(startpep.iloc[-1, 5]).ljust(5)
                    new_row[6] = str('%8.3f' % var1[0]).rjust(8)
                    new_row[7] = str('%8.3f' % var1[1]).rjust(8)
                    new_row[8] = str('%8.3f' % var1[2]).rjust(8)
                    new_row[9] = "1.00".rjust(6)
                    new_row[10] = "1.00".ljust(6)
                    outputfile.write("%s%s  %s%s%s%s%s%s%s%s %s\n" % (new_row[0], new_row[1], new_row[2], new_row[3],
                                                                      new_row[4], new_row[5], new_row[6], new_row[7],
                                                                      new_row[8], new_row[9], new_row[10]))
                    new_row2 = {0: "ATOM",
                               1: str(int(startpep.iloc[-1, 1]) + 1),
                               2: "C",
                               3: startpep.iloc[-1, 3],
                               4: "B",
                               5: str(startpep.iloc[-1, 5]),
                               6: var1[0],
                               7: var1[1],
                               8: var1[2],
                               9: "1.00",
                               10: "1.000"}
                    startpep.loc[len(startpep)] = new_row2
                    continue
                if startpep.iloc[-1, 2] == "C":
                    startpepPDB = prody.parsePDB("Start_peptide2 - Copy.pdb")
                    if int(startpep.iloc[-1, 5]) & int(startpep.iloc[-2, 5]) & int(startpep.iloc[-3, 5]) == 2:
                        PhiPsi_All[startpep.iloc[-1, 5], 0] = prody.calcPhi(startpepPDB[int(startpep.iloc[-1, 5]),])
                        PhiPsi_All[startpep.iloc[-1, 5], 1] = ndrd_spec_phi(sys.argv[startpep.iloc[-1, 5]],
                                                                            sys.argv[startpep.iloc[-4, 5]],
                                                                            PhiPsi_All[startpep.iloc[-1, 5], 0])
                        PhiPsi = PhiPsi_All[startpep.iloc[-1, 5]]
                    else:
                        PhiPsi = PhiPsi_All[int(startpep.iloc[-1, 5])]
                    x1 = calculateCoordinates(Point3D(startpep.iloc[-4, 6:9]),
                                              Point3D(startpep.iloc[-2, 6:9]),
                                              Point3D(startpep.iloc[-1, 6:9]),
                                              1.23,
                                              cdl_spec(sys.argv[int(startpep.iloc[-4, 5])],
                                                       sys.argv[int(startpep.iloc[-1, 5])], "ACO", PhiPsi[0], PhiPsi[1]),
                                              PhiPsi[1])  # Need metric for this
                    var1 = np.around(x1, decimals=3)
                    new_row = [None] * 11
                    new_row[0] = "ATOM".ljust(7)
                    new_row[1] = str(int(startpep.iloc[-1, 1]) + 1).rjust(4)
                    new_row[2] = "O".ljust(4)
                    new_row[3] = str(startpep.iloc[-1, 3]).ljust(4)
                    new_row[4] = "B".ljust(4)
                    new_row[5] = str(startpep.iloc[-1, 5]).ljust(5)
                    new_row[6] = str('%8.3f' % var1[0]).rjust(8)
                    new_row[7] = str('%8.3f' % var1[1]).rjust(8)
                    new_row[8] = str('%8.3f' % var1[2]).rjust(8)
                    new_row[9] = "1.00".rjust(6)
                    new_row[10] = "1.00".ljust(6)
                    outputfile.write("%s%s  %s%s%s%s%s%s%s%s %s\n" % (new_row[0], new_row[1], new_row[2], new_row[3],
                                                                      new_row[4], new_row[5], new_row[6], new_row[7],
                                                                      new_row[8], new_row[9], new_row[10]))
                    new_row2 = {0: "ATOM",
                               1: str(int(startpep.iloc[-1, 1]) + 1),
                               2: "O",
                               3: startpep.iloc[-1, 3],
                               4: "B",
                               5: str(startpep.iloc[-1, 5]),
                               6: var1[0],
                               7: var1[1],
                               8: var1[2],
                               9: "1.00",
                               10: "1.000"}
                    startpep.loc[len(startpep)] = new_row2
                    print(startpep)
                    continue
                if startpep.iloc[-1, 2] == "O":
                    PhiPsi = PhiPsi_All[int(startpep.iloc[-1, 5]) + 1]
                    x1 = calculateCoordinates(Point3D(startpep.iloc[-3, 6:9]),
                                              Point3D(startpep.iloc[-1, 6:9]),
                                              Point3D(startpep.iloc[-2, 6:9]),
                                              1.3324,
                                              cdl_spec(sys.argv[int(startpep.iloc[-4, 5])],
                                                       sys.argv[int(startpep.iloc[-1, 5]) + 1],
                                                       "ACN", PhiPsi[0], PhiPsi[1]),
                                              odl_spec(sys.argv[int(startpep.iloc[-4, 5])],
                                                       sys.argv[int(startpep.iloc[-1, 5]) + 1],
                                                       PhiPsi_All[2, 1],
                                                       PhiPsi_All[3, 0]))
                    var1 = np.around(x1, decimals=3)
                    new_row = [None] * 11
                    new_row[0] = "ATOM".ljust(7)
                    new_row[1] = str(int(startpep.iloc[-1, 1]) + 1).rjust(4)
                    new_row[2] = "N".ljust(4)
                    new_row[3] = aa123[sys.argv[int(startpep.iloc[-1, 5]) + 1]].ljust(4)
                    new_row[4] = "B".ljust(4)
                    new_row[5] = str(int(startpep.iloc[-1, 5]) + 1).ljust(5)
                    new_row[6] = str('%8.3f' % var1[0]).rjust(8)
                    new_row[7] = str('%8.3f' % var1[1]).rjust(8)
                    new_row[8] = str('%8.3f' % var1[2]).rjust(8)
                    new_row[9] = "1.00".rjust(6)
                    new_row[10] = "1.00".ljust(6)
                    outputfile.write("%s%s  %s%s%s%s%s%s%s%s %s\n" % (new_row[0], new_row[1], new_row[2], new_row[3],
                                                                      new_row[4], new_row[5], new_row[6], new_row[7],
                                                                      new_row[8], new_row[9], new_row[10]))
                    new_row2 = {0: "ATOM",
                               1: str(int(startpep.iloc[-1, 1]) + 1),
                               2: "N",
                               3: aa123[sys.argv[int(startpep.iloc[-1, 5]) + 1]],
                               4: "B",
                               5: str(int(startpep.iloc[-1, 5]) + 1),
                               6: var1[0],
                               7: var1[1],
                               8: var1[2],
                               9: "1.00",
                               10: "1.000"}
                    startpep.loc[len(startpep)] = new_row2
                    continue
                if startpep.iloc[-1, 2] == "N":
                    PhiPsi = ndrd_spec_random(
                        sys.argv[int(startpep.iloc[-1, 5]) + 1], sys.argv[int(startpep.iloc[-1, 5])])
                    PhiPsi_All[(int(startpep.iloc[-1, 5]))] = PhiPsi
                    x1 = calculateCoordinates(Point3D(startpep.iloc[-4, 6:9]),
                                              Point3D(startpep.iloc[-3, 6:9]),
                                              Point3D(startpep.iloc[-1, 6:9]),
                                              1.457,
                                              cdl_spec(sys.argv[int(startpep.iloc[-1, 5])],
                                                       sys.argv[int(startpep.iloc[-1, 5]) + 1],
                                                       "CNA", PhiPsi[0], PhiPsi[1]),
                                              odl_spec(sys.argv[int(startpep.iloc[-1, 5])],
                                                       sys.argv[int(startpep.iloc[-1, 5]) + 1],
                                                       PhiPsi_All[2, 1],
                                                       PhiPsi_All[3, 0]))
                    var1 = np.around(x1, decimals=3)
                    new_row = [None] * 11
                    new_row[0] = "ATOM".ljust(7)
                    new_row[1] = str(int(startpep.iloc[-1, 1]) + 1).rjust(4)
                    new_row[2] = "CA".ljust(4)
                    new_row[3] = startpep.iloc[-1, 3].ljust(4)
                    new_row[4] = "B".ljust(4)
                    new_row[5] = str(startpep.iloc[-1, 5]).ljust(5)
                    new_row[6] = str('%8.3f' % var1[0]).rjust(8)
                    new_row[7] = str('%8.3f' % var1[1]).rjust(8)
                    new_row[8] = str('%8.3f' % var1[2]).rjust(8)
                    new_row[9] = "1.00".rjust(6)
                    new_row[10] = "1.00".ljust(6)
                    outputfile.write("%s%s  %s%s%s%s%s%s%s%s %s\n" % (new_row[0], new_row[1], new_row[2], new_row[3],
                                                                      new_row[4], new_row[5], new_row[6], new_row[7],
                                                                      new_row[8], new_row[9], new_row[10]))
                    new_row2 = {0: "ATOM",
                               1: str(int(startpep.iloc[-1, 1]) + 1),
                               2: "CA",
                               3: startpep.iloc[-1, 3],
                               4: "B",
                               5: str(startpep.iloc[-1, 5]),
                               6: var1[0],
                               7: var1[1],
                               8: var1[2],
                               9: "1.00",
                               10: "1.000"}
                    startpep.loc[len(startpep)] = new_row2
                    continue
        outputfile.write("ENDMDL" + "\n")
        outputfile.close()
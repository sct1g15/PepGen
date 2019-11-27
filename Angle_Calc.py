import numpy as np
import math
from BondDict import *
import pandas as pd
from bisect import bisect_right
import Bio.PDB.vectors
import random

CDL = pd.read_table("KernRegr_CDL_v1.2.1_Oct14-2009.txt",
                    skiprows=87, delim_whitespace=True)
ODL = pd.read_table("omegaCDL_OmegaBetweenAsPhi1Psi0_KernRegr_v1.3.1_Aug12-2011.txt",
                    skiprows=87, delim_whitespace=True)
NDRD = pd.read_table(r"NDRD_TCBIG.txt",
                     skiprows=59, delim_whitespace=True)

def get_vector2(N) :
    return Bio.PDB.vectors.Vector(N.x, N.y, N.z)

def getAngle(p0 = np.array([0.0,0.0,0.0]), p1 = np.array([0.0,0.0,0.0]), p2 = np.array([0.0,0.0,0.0])):
    v12 = p0 - p1
    v32 = p2 - p1
    cosine_angle = np.dot(v12, v32) / (np.linalg.norm(v12) * np.linalg.norm(v32))
    angle = math.degrees(np.arccos(cosine_angle, ))
    return angle

######
from Bio.PDB import *
import warnings
def calculateCoordinates(refA, refB, refC, L, ang, di):
    AV = get_vector2(refA)
    BV = get_vector2(refB)
    CV = get_vector2(refC)

    CA = AV - CV
    CB = BV - CV

    ##CA vector
    AX = CA[0]
    AY = CA[1]
    AZ = CA[2]

    ##CB vector
    BX = CB[0]
    BY = CB[1]
    BZ = CB[2]

    ##Plane Parameters
    A = (AY * BZ) - (AZ * BY)
    B = (AZ * BX) - (AX * BZ)
    G = (AX * BY) - (AY * BX)

    ##Dot Product Constant
    F = math.sqrt(BX * BX + BY * BY + BZ * BZ) * L * math.cos(ang * (math.pi / 180.0))

    ##Constants
    const = math.sqrt(math.pow((B * BZ - BY * G), 2) * (-(F * F) * (A * A + B * B + G * G) + (
                B * B * (BX * BX + BZ * BZ) + A * A * (BY * BY + BZ * BZ) - (2 * A * BX * BZ * G) + (
                    BX * BX + BY * BY) * G * G - (2 * B * BY) * (A * BX + BZ * G)) * L * L))
    denom = (B * B) * (BX * BX + BZ * BZ) + (A * A) * (BY * BY + BZ * BZ) - (2 * A * BX * BZ * G) + (
                BX * BX + BY * BY) * (G * G) - (2 * B * BY) * (A * BX + BZ * G)

    X = ((B * B * BX * F) - (A * B * BY * F) + (F * G) * (-A * BZ + BX * G) + const) / denom

    if ((B == 0 or BZ == 0) and (BY == 0 or G == 0)):
        const1 = math.sqrt(G * G * (-A * A * X * X + (B * B + G * G) * (L - X) * (L + X)))
        Y = ((-A * B * X) + const1) / (B * B + G * G)
        Z = -(A * G * G * X + B * const1) / (G * (B * B + G * G))
    else:
        Y = ((A * A * BY * F) * (B * BZ - BY * G) + G * (-F * math.pow(B * BZ - BY * G, 2) + BX * const) - A * (
                    B * B * BX * BZ * F - B * BX * BY * F * G + BZ * const)) / ((B * BZ - BY * G) * denom)
        Z = ((A * A * BZ * F) * (B * BZ - BY * G) + (B * F) * math.pow(B * BZ - BY * G, 2) + (A * BX * F * G) * (
                    -B * BZ + BY * G) - B * BX * const + A * BY * const) / ((B * BZ - BY * G) * denom)

    # GET THE NEW VECTOR from the orgin
    D = Vector(X, Y, Z) + CV
    with warnings.catch_warnings():
        # ignore inconsequential warning
        warnings.simplefilter("ignore")
        temp = calc_dihedral(AV, BV, CV, D) * (180.0 / math.pi)

    di = di - temp
    rot = rotaxis(math.pi * (di / 180.0), CV - BV)
    D = (D - BV).left_multiply(rot) + BV

    return D.get_array()

def cdl_spec(pos1, pos2, bond, phi, psi): #Calculates bond angles within amino acids
    for key, value in Lib_class.items():
        if aa123[pos1] in value:
            x1 = key
    for key, value in x_class.items():
        if aa123[pos2] in value:
            x2 = key
    if psi == 180 or psi == 175:
        psi = -180
    if phi == 180 or phi == 175:
        phi = -180
    x3 = str(x1) + str(x2)
    x4 = CDL.loc[(CDL["ResTypeGroup"] == x3) & (CDL["Phi"] == round(phi/10.0)*10) & (CDL["Psi"] == round(psi/10.0)*10)]
    meanloc = "m" + str(bond)
    print(phi, psi)
    mean = float(x4[meanloc])
    stdloc = "s" + str(bond)
    std = float(x4[stdloc])
    return np.random.normal(mean, std)

def odl_spec(pos1, pos2, psi, phi): #Calculates Omega Bond angles
    for key, value in Lib_class.items():
        if aa123[pos1] in value:
            x1 = key
    for key, value in x_class.items():
        if aa123[pos2] in value:
            x2 = key
    if psi == 180 or phi == 175:
        psi = -180
    if phi == 180 or phi == 175:
        phi = -180
    x3 = str(x1) + str(x2)
    x4 = ODL.loc[(ODL["Phi(+1)"] == (round(phi/10)*10)) & (ODL["Psi(0)"] == (round(psi/10)*10)) & (CDL["ResTypeGroup"] == x3)]
    mean = x4["mW(+1)"]
    std = x4["sW(+1)"]
    return float(np.random.normal(mean, std))


def ndrd_spec_random(pos1, pos2): #Calculates Phi and Psi angle probabilties for adjacent amino acids, right to left (N to C)
    x1 = aa123[pos1]
    x2 = aa123[pos2]
    NDRD_subset2 = NDRD.loc[(NDRD["AminoAcid"] == x1) & (NDRD["AA2"] == x2) & (NDRD["Pos"] == "right")]
    NDRD_subset = NDRD.loc[(NDRD["AminoAcid"] == x1) & (NDRD["AA2"] == x2) & (NDRD["Pos"] == "right")].iloc[:, 7]
    NDRD_list = NDRD_subset.values.tolist()
    return NDRD_subset2.iloc[bisect_right(NDRD_list, random.uniform(0.00001, 1)), 3:5]


def ndrd_spec_phi(pos1, pos2, phi): #Calculates Phi and Psi angle probabilties for adjacent amino acids, right to left (N to C)
    x1 = aa123[pos1]
    x2 = aa123[pos2]

    NDRD_subset2 = NDRD.loc[(NDRD["AminoAcid"] == x1) & (NDRD["AA2"] == x2)
                            & (NDRD["Pos"] == "right") & (NDRD["Phi"] == (round(phi/5.0)*5))]
    NDRD_subset = NDRD.loc[(NDRD["AminoAcid"] == x1) & (NDRD["AA2"] == x2)
                           & (NDRD["Pos"] == "right") & (NDRD["Phi"] == (round(phi/5.0)*5))].iloc[:, 7]
    NDRD_list = NDRD_subset.values.tolist()
    return NDRD_subset2.iloc[bisect_right(NDRD_list, random.uniform((NDRD_subset2.iloc[1, 7] - NDRD_subset2.iloc[1, 5]),
                                                                    NDRD_subset2.iloc[-1, 7])), 3:5]["Psi"]


def ndrd_spec_psi(pos1, pos2, psi): #Calculates Phi and Psi angle probabilties for adjacent amino acids, right to left (N to C)
    x1 = aa123[pos1]
    x2 = aa123[pos2]
    NDRD_subset2 = NDRD.loc[(NDRD["AminoAcid"] == x1) & (NDRD["AA2"] == x2)
                            & (NDRD["Pos"] == "right") & (NDRD["Phi"] == (round(psi/5.0)*5))]
    NDRD_subset = NDRD.loc[(NDRD["AminoAcid"] == x1) & (NDRD["AA2"] == x2)
                           & (NDRD["Pos"] == "right") & (NDRD["Phi"] == (round(psi/5.0)*5))].iloc[:, 7]
    NDRD_list = NDRD_subset.values.tolist()
    return NDRD_subset2.iloc[bisect_right(NDRD_list, random.uniform((NDRD_subset2.iloc[1, 7] - NDRD_subset2.iloc[1, 5]),
                                                                    NDRD_subset2.iloc[-1, 7])), 3:5]["Phi"]


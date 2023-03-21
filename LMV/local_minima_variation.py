"""
Pyhton file that provides the classes neccessary to perform
local minima varitation calululations. This code is writtent to
be as agnostic to the source of data as possibel, but was developed using
nupack4, fyi.
copyright 2023 Jennifer Pearl
"""

from nupack import *
import math
import copy

import nupackAPI_Sara2_Ver1 as nupack_api
from nupackAPI_Sara2_Ver1 import Sara2SecondaryStructure, Sara2StructureList, EnsembleVariation



def test_LMV():
    print("Enter single strand RNA sequence")
    sequence = "AGGGUGGUACCGCGAUAAUCAAUCGUCCCUUCGUGUAAACGAAGGGGCG" #input()

    print("Enter Kcal delta span to look at")
    span = 7 # int(input())

    print("Enter kcal unit to plot by")
    units = .5# int(input())

    EV_test: EnsembleVariation = EnsembleVariation()

    ev_result = EV_test.process_ensemble_variation(sequence, span, units)

    print(ev_result.group_ev_list)

test_LMV()
    




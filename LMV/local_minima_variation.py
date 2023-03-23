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
from nupackAPI_Sara2_Ver1 import Sara2SecondaryStructure, Sara2StructureList, EnsembleVariation, EVResult

debug:bool = True

def test_LMV():

    sequence = ''
    target = ''
    span = 0
    units = 0

    if debug is True:
        print("using debug")
        sequence = 'AGGGUGGUACCGCGAUAAUCAAUCGUCCCUUCGUGUAAACGAAGGGGCG'
        target = '.(((.......(((((.....))))))))(((((....)))))......'
        span = 5
        units = 1
    else:
        print("Enter single strand RNA sequence")
        sequence = input()

        #print("Enter target structure")
        #target = input()

        print("Enter Kcal delta span to look at")
        span = input()

        print("Enter kcal unit to plot by")
        units = input()

        print("Enter design name")
        name = input()

    

   
    

    EV_test: EnsembleVariation = EnsembleVariation()
    ev_result:EVResult
    switch_result: EVResult

    ev_result = EV_test.process_ensemble_variation(sequence, int(span), float(units))

    #print(ev_result.group_ev_list)

    new_list_string = []
    for ev in ev_result.group_ev_list:
        ev_value = ev.ev_normalized
        new_list_string.append(ev_value)

    #new_switch_string = []
    #for ev in switch_result.group_ev_list:
    #    ev_value = ev.ev_normalized
    #    new_switch_string.append(ev_value)
   # 
    print(f'Results for name={name}, sequence={sequence} ')
    print("LMV_U")
    print(new_list_string)
    #print()
    #print("LMV_US")
    #print(new_switch_string)


test_LMV()
    




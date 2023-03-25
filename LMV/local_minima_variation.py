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

import nupackAPI_Sara2_Ver2 as nupack_api
from nupackAPI_Sara2_Ver2 import Sara2SecondaryStructure, Sara2StructureList, EnsembleVariation, EVResult

debug:bool = False

def test_LMV():

    sequence = ''
    target = ''
    folded = ''
    span = 0
    units = 0
    name = ''

    if debug is True:
        print("using debug")
        sequence = 'GCCAUCGCAUGAGGAUAUGCUCCGGUUUCCGGAGCAGAAGGCAUGUCAUAAGACAUGAGGAUCACCCAUGUAGUUAAGAUGGCA'
        target = '........(((......(((.............))).....)))........................................'
        folded = '((((((.((((......((((((((...)))))))).....))))...(((.(((((.((....)))))))..))).)))))).'
        span = 10
        units = .5
        name = "09_eli"
    else:
        print("Enter single strand RNA sequence")
        sequence = input()

        print("Enter target structure")
        target = input()

        print("Enter predicted 2nd state folded structure")
        folded = input()

        print("Enter Kcal delta span to look at")
        span = input()

        print("Enter kcal unit to plot by")
        units = input()

        print("Enter design name")
        name = input()

    

   
    

    EV_test: EnsembleVariation = EnsembleVariation()
    ev_result_mfe:EVResult
    ev_result_rel:EVResult
    switch_result_target: EVResult
    switch_result_folded: EVResult

    ev_result_mfe, ev_result_rel, switch_result_target, switch_result_folded = EV_test.process_ensemble_variation(sequence, int(span), float(units), folded, target)

    #print(ev_result.group_ev_list)

    new_list_string_mfe = []
    for ev in ev_result_mfe.group_ev_list:
        ev_value = ev.ev_normalized
        new_list_string_mfe.append(ev_value)

    new_list_string_rel = []
    for ev in ev_result_rel.group_ev_list:
        ev_value = ev.ev_normalized
        new_list_string_rel.append(ev_value)

    new_switch_string = []
    for ev in switch_result_target.group_ev_list:
        ev_value = ev.ev_normalized
        new_switch_string.append(ev_value)

    new_switch_string_folded = []
    for ev in switch_result_folded.group_ev_list:
        ev_value = ev.ev_normalized
        new_switch_string_folded.append(ev_value)
    
    print(f'Results for name={name}, sequence={sequence}, span={span}, units={units} ')
    print("LMV_U mfe")
    print(new_list_string_mfe)
    print()
    print("LMV_U rel")
    print(new_list_string_rel)
    print()
    print("LMV_US target")
    print(new_switch_string)
    print()
    print("LMV_US folded")
    print(new_switch_string_folded)
    print()


test_LMV()
    




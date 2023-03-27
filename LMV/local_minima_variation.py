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
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
import matplotlib
from typing import List


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
    designID: int= 0
    labname: str = ''

    if debug is True:
        print("using debug")
        sequence = 'GCCAUCGCAUGAGGAUAUGCUCCGGUUUCCGGAGCAGAAGGCAUGUCAUAAGACAUGAGGAUCACCCAUGUAGUUAAGAUGGCA'
        target = '........(((......(((.............))).....)))........................................'
        folded = '((((((.((((......((((((((...)))))))).....))))...(((.(((((.((....)))))))..))).)))))).'
        span = 5
        units = .5
        name = "09_eli"
        designID = 12345
        labname = "Tbox Round 1"
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

        print("Enter designID")
        designID = int(input())

        print("Enter Lab Name")
        labname = input()



    
    info_str:str = f'Eterna Lab={labname}\nDesign ID={designID}'

   
    

    EV_test: EnsembleVariation = EnsembleVariation()
    ev_result_mfe:EVResult
    ev_result_rel:EVResult
    switch_result_target: EVResult
    switch_result_folded: EVResult

    ev_result_mfe, ev_result_rel, switch_result_target, switch_result_folded = EV_test.process_ensemble_variation(sequence, int(span), float(units), folded, target)

    #print(ev_result.group_ev_list)

    time_span: List[float] = []
    tick_span = []
    index_span = range(len(ev_result_mfe.groups_list))
    

    mfe_value:float = ev_result_mfe.groups_list[0].mfe_freeEnergy
    seed_value:float = mfe_value
    tick_value:float = 0
    
    num_samples: int = len(ev_result_mfe.groups_list)
    for index in range(num_samples):
        seed_value = seed_value + float(units)
        tick_value = tick_value + float(units)
        time_span.append(seed_value)
        tick_span.append(tick_value)

    

    new_list_string_mfe: List[float] = []
    for ev in ev_result_mfe.group_ev_list:
        ev_value = ev.ev_normalized
        new_list_string_mfe.append(ev_value)

    new_list_string_rel: List[float] = []
    for ev in ev_result_rel.group_ev_list:
        ev_value = ev.ev_normalized
        new_list_string_rel.append(ev_value)

    new_switch_string: List[float] = []
    for ev in switch_result_target.group_ev_list:
        ev_value = ev.ev_normalized
        new_switch_string.append(ev_value)

    new_switch_string_folded: List[float] = []
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

    plt.title(f'LMV Switch plot for {name}')
    plt.suptitle(info_str, fontsize=14)
    #fig = plt.figure()
    #ax = fig.add_axes((0, float(span), 0, 50))
    #ax.set_xticks(tick_span)
    plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}')) # 2 decimal places
    plt.plot(time_span, new_list_string_mfe, 'b^-', label='LMV_U mfe')
    plt.plot(time_span, new_list_string_rel, 'ro-', label='LMV_U rel')
    #plt.plot(time_span, new_switch_string, 'k.-', label='LMV_US target')
    plt.plot(time_span, new_switch_string_folded, 'gs-', label='LMV_US folded')
    y_ticks = [0,5,10,15,20,25,30,35,40,45,50]
    plt.xticks(time_span)
    plt.yticks(y_ticks)
    plt.grid(True)   
    plt.legend(loc='lower right',fontsize="x-small")
    plt.subplots_adjust(top=.8, bottom=.2, left=.12, right=.95)  
    plt.tick_params(axis='x',labelrotation=90)  
    plt.ylabel("Local Minima Structure Variation (LMSV)")
    plt.xlabel("Local Kcal Energy along Ensemble")
    file_name:str = f'{name}_{designID}'
    folder_name:str = '/home/ubuntu/rna_analysis/tbox_round1'
    plt.savefig(f'{folder_name}/{file_name}.png')
    



test_LMV()
    




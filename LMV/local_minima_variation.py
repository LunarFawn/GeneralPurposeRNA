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
from datetime import datetime
import numpy as np


import nupackAPI_Sara2_Ver2 as nupack_api
from nupackAPI_Sara2_Ver2 import Sara2SecondaryStructure, Sara2StructureList, EnsembleVariation, EVResult

debug:bool = True


from bisect import bisect_left

def take_closest(myList, myNumber):
    """
    Assumes myList is sorted. Returns closest value to myNumber.

    If two numbers are equally close, return the smallest number.
    """
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return myList[0]
    if pos == len(myList):
        return myList[-1]
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
        return after
    else:
        return before
    
def test_LMV():

    sequence = ''
    target = ''
    folded = ''
    span = 0
    units = 0
    name = ''
    designID: int= 0
    labname: str = ''
    folder_name:str = ''
    ligand_oligo_energy:float = 0
    folded_energy_ligoligo: float = 0
    ligand_oligo_name:str = ''
    eterna_score:float = 0
    fold_change:float = 0

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
        folder_name:str = '/home/ubuntu/rna_analysis/tbox_round1/debug'
        ligand_oligo_energy:float = 10
        folded_energy_ligoligo: float = -29
        ligand_oligo_name:str = ''
        eterna_score:float = 100
        fold_change:float = 500
    else:
        print("Enter single strand RNA sequence")
        sequence = input()

        print("Enter target structure")
        target = input()

        print("Enter predicted 2nd state folded structure")
        folded = input()

        print("Enter Energy of folded structure with ligand/oligo bound")
        folded_energy_ligoligo = float(input())

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

        print("Enter Eterna Score")
        eterna_score = float(input())

        print("Enter Fold Change")
        fold_change = float(input())

        

        folder_name:str = '/home/ubuntu/rna_analysis/tbox_round1/poly'


    
    info_str:str = f'Eterna_Score = {eterna_score}, FoldChange = {fold_change}'

   
    

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
        #time_span is teh MFE values
        #tick_span is the index value (i.e. 0, 0.5, 1, 1.5)
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

    #get LMSV deltas
    start_value:float = 2
    #last value of span
    stop_value: float = tick_span[-1]
    #get polymorphicicyt
    mfe_fold_EV_delta = EV_test.lmsv_delta(start_value, stop_value, ev_result_mfe, switch_result_folded, tick_span)
    delta_message: str = f'Polymorphicity Level (2Kcal delta to end) = {mfe_fold_EV_delta}'

    mfe_ev_oligo:float = take_closest(time_span, folded_energy_ligoligo)
    mfe_ev_oligo_index: int = time_span.index(mfe_ev_oligo)
    ev_oligo_folded:float = new_switch_string_folded[mfe_ev_oligo_index]

    csv_log_results: List[str]=[]
    csv_log_results.append("Kcal,LMSV_U_mfe,LMSV_U_rel,LMSV_US_target,LMSV_US_folded\n")
    for index in range(len(new_list_string_mfe)):
        kcal = time_span[index]
        LMSV_U_mfe = new_list_string_mfe[index]
        LMSV_U_rel = new_list_string_rel[index]
        LMSV_US_target = new_switch_string[index]
        LMSV_US_folded = new_switch_string_folded[index]
        line:str = f'{kcal},{LMSV_U_mfe},{LMSV_U_rel},{LMSV_US_target},{LMSV_US_folded}\n'
        csv_log_results.append(line)
    
    print(f'Results for name={name}\nsequence={sequence}\nspan={span}\nunits={units}\ntarget={target}\nfolded={folded}\nDesignID={designID}\nLabName={labname}\n')
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

    #now save teh data
    import time
    timestr = time.strftime("%Y%m%d-%H%M%S")
    fig, ax = plt.subplots()
    
    plt.suptitle(f'LMV Switch plot for {name}\nEterna Lab = {labname}\nDesign ID = {designID}\n',fontsize=12)
    plt.title(info_str, fontsize=10)
    #fig = plt.figure()
    
    #ax.set_xticks(tick_span)
    plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}')) # 2 decimal places
    plt.plot(time_span, new_list_string_mfe, 'b^-', label='LMV_U mfe')
    plt.plot(time_span, new_list_string_rel, 'ro-', label='LMV_U rel')
    plt.plot(time_span, new_switch_string, 'kD-', label='LMV_US target')
    plt.plot(time_span, new_switch_string_folded, 'gs-', label='LMV_US folded')
    #y_ticks = [0,5,10,15,20,25,30,35,40,45,50]
    y_ticks = np.arange(-10,65, step=5)
    plt.xticks(time_span)
    plt.yticks(y_ticks)
    #plt.yticks()
    plt.grid(True)   
    plt.legend(loc='lower right',fontsize="x-small")
    plt.subplots_adjust(top=.8, bottom=.2, left=.12, right=.95)  
    plt.tick_params(axis='x',labelrotation=90)  
    plt.ylabel("Local Minima Structure Variation (LMSV)")
    plt.xlabel("Local Kcal Energy along Ensemble")
    plt.figtext(0.54, 0.01, delta_message, ha="center", fontsize=10, bbox={"facecolor":"orange", "alpha":.5, "pad":2})
    trans = ax.get_xaxis_transform()
    delat_energy:float = 2
    lower_range_folded_energy: float = folded_energy_ligoligo - (delat_energy/2)
    uper_range_folded_energy: float = folded_energy_ligoligo + (delat_energy/2)
    plt.axvline(x=folded_energy_ligoligo, color="green", linestyle="--")
    plt.axvline(x=lower_range_folded_energy, color="green", linestyle=":")
    plt.axvline(x=uper_range_folded_energy, color="green", linestyle=":")
    plt.text(lower_range_folded_energy, .06, '   2nd State', transform=trans, fontsize=7)
    plt.text(lower_range_folded_energy, .01, ' 2Kcal range', transform=trans, fontsize=7)
    plt.text(folded_energy_ligoligo, .06, '    2nd State', transform=trans, fontsize=7)
    plt.text(folded_energy_ligoligo, .01, '  folded Energy', transform=trans, fontsize=7)

    ev_mfe_lower:float = time_span[0]
    ev_mfe_upper:float = mfe_value + delat_energy
    
    plt.axvline(x=ev_mfe_lower, color="blue", linestyle=":")
    plt.axvline(x=ev_mfe_upper, color="blue", linestyle=":")
    plt.text(ev_mfe_lower, .06, '   1st State', transform=trans, fontsize=7)
    plt.text(ev_mfe_lower, .01, ' 2Kcal range', transform=trans, fontsize=7)
    #ax.set_ybound(lower=0, upper=70)
    

    file_name:str = f'{name}_{designID}'
    
    plt.savefig(f'{folder_name}/{file_name}_{timestr}.png')

    #now save data to csv file
    csv_record_pathstr = f'{folder_name}/{file_name}_{timestr}.csv'
    csv_lines:List[str]=[]
    with open(csv_record_pathstr, 'w') as csv_file:
        #first write teh header
        csv_lines.append(f'Local Minima Structure Variation Data\n')
        csv_lines.append(f'Creation Date={datetime.now()}\n')
        csv_lines.append("---------------------------------------\n")
        csv_lines.append("***DESIGN INFO***\n")
        csv_lines.append(f'Design Name = {name}\n')
        csv_lines.append(f'DesignID = {designID}\n')
        csv_lines.append(f'Lab Name = {labname}\n')
        csv_lines.append(f'Sequence = {sequence}\n')
        csv_lines.append(f'Eterna_Score = {eterna_score}\n')
        csv_lines.append(f'FoldChange = {fold_change}\n')
        csv_lines.append(f'2nd State Target Structure = {target}\n')
        csv_lines.append(f'2nd State Folded Structure = {folded}\n')
        csv_lines.append(f'2nd State Folded Oligo Energy = {folded_energy_ligoligo}\n')
        csv_lines.append(f'Energy Span from MFE = {span}\n')
        csv_lines.append(f'Energy span units = {units}\n')
        csv_lines.append("---------------------------------------\n")
        csv_lines.append("***RAW DATA***\n")
        csv_lines = csv_lines + csv_log_results
        csv_lines.append("---------------------------------------\n")
        csv_lines.append("***METRICS***\n")
        csv_lines.append(f'Polymorphicity Level (2kcal to end of sample) = {mfe_fold_EV_delta}\n')
        csv_lines.append(f'LMV_US_folded at folded with lignad/oligo energy = {ev_oligo_folded}\n')
        csv_lines.append("---------------------------------------\n")
        csv_lines.append("EOF\n")
        csv_file.writelines(csv_lines)

    print("done with analysis")

test_LMV()
    




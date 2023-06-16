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
#import matplotlib.pyplot as plt
#from matplotlib.ticker import StrMethodFormatter
#import matplotlib
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
    number_of_clusters:int = 0

    if debug is True:
        print("using debug")
        sequence = 'GCCAUCGCAUGAGGAUAUAAUCGGGUUUCCCGUUUAGAAGGCAUGUCAUAAGACAUGAGGAUCACCCAUGUAGUUAAGAUGGCA'
        target = '........(((......(((.............))).....)))........................................'
        folded = '((((((.((((......(((.((((...)))).))).....))))...(((.(((((.((....)))))))..))).)))))).'
        folded_energy_ligoligo: float = -21.8
        span = 6
        units = 1
        name = "09_eli"
        designID = 12345
        labname = "Tbox Round 1"
        folder_name:str = '/home/ubuntu/rna_analysis/tbox_round1/debug'
        ligand_oligo_energy:float = 10
        
        ligand_oligo_name:str = ''
        eterna_score:float = 100
        fold_change:float = 500
        number_of_clusters:int = 1000
    else:
        print("Enter single strand RNA sequence")
        sequence = input()

        print("Enter target structure")
        target = '........(((......(((.............))).....)))........................................' #input()

        print("Enter predicted 2nd state folded structure")
        folded = input()

        print("Enter Energy of folded structure with ligand/oligo bound")
        folded_energy_ligoligo = float(input())

        print("Enter Kcal delta span to look at")
        
        span = '6'#input()
        print(f'span is {span}')

        print("Enter kcal unit to plot by")
        units = '1' #input()
        print(f'units is {units}')

        print("Enter design name")
        name = 'test'#input()

        print("Enter designID")
        designID = '99'#int(input())

        print("Enter Lab Name")
        labname = 'test'#input()

        print("Enter Eterna Score")
        eterna_score = float(input())

        print("Enter Fold Change")
        fold_change = float(input())

        print("Enter Number of Clusters")
        number_of_clusters = -9 #int(input())

        

        folder_name:str = '/home/ubuntu/rna_analysis/PNAS_FMN/SSNG1'


    
    info_str:str = f'Eterna_Score = {eterna_score}, FoldChange = {fold_change}, Num Clusters = {number_of_clusters}'

   
    

    EV_test: EnsembleVariation = EnsembleVariation()
    ev_result_mfe:EVResult
    ev_result_rel:EVResult
    switch_result_target: EVResult
    switch_result_folded: EVResult

    temp_1:int = 35
    temp_2:int = 36
    temp_3:int = 37
    temp_4:int = 38
    temp_5:int = 39

    #temp_list: List[int] = [32,33,34, 35,36, 37, 38,39,40,41,42]
    temp_list: List[int] = [36, 37, 38]
    score_list:List[float] = []

    score: float = 0
    for temp in temp_list: 
        
        value = EV_test.process_ensemble_variation(sequence, int(span), float(units), folded, target, folded_energy_ligoligo, temp)
        score = score + value
        score_list.append(value)
    
    num_scores: int = len(temp_list)
    print(f'Raw score is {score} of {len(temp_list)}')
    
    num_zero_values: int = score_list.count(0)
    modified_score: float = (score - (num_zero_values)) / num_scores
    print(f'modified_score is {modified_score}')
    
    max_fold_score:float = 3    
    if modified_score > (max_fold_score):
        #its probbaly only going to be in 2nd state state so give penatly
        offset = (modified_score - max_fold_score)
        modified_score = max_fold_score - offset


    why_not:float = (27/3)*modified_score
    if why_not > 10:
        predicted_foldchange:str = f'Good Switch Predicted. Fold change predicterd to be {why_not} if greater than 10. If less then expect closer to 1'
    else:
        predicted_foldchange:str = f'Bad Switch Predicted. Fold Change Predicted to be {9*why_not}'
    print(predicted_foldchange)
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


test_LMV()
    




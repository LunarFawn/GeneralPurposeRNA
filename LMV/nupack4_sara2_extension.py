"""
File that provides an interface to standard python nupack4 install on linux
Intention is to be able to access this via a docker that will be accessable on docker hub
that has nupack4 setup and ready for this project to consume
"""

from typing import List, Dict
import struct
from nupack import *
import pandas as pd
import sys
import openpyxl
from copy import deepcopy
from dataclasses import dataclass
from datetime import datetime, timedelta
import threading
import time
import numpy as np
import collections
from enum import Enum

import nupackAPI_Sara2_Ver2 as nupack_api
from nupackAPI_Sara2_Ver2 import Sara2SecondaryStructure, Sara2StructureList, EnsembleVariation, EVResult
import serena_sara2_api as serena_api
from serena_sara2_api import SingleEnsembleGroup, MultipleEnsembleGroups

my_model = Model

rna_model='rna95-nupack3'  

class NUPACK4Interface():
    """
    Class for nupack4 interface for sara2 logic intended for serena packag
    """

    def __init__(self) -> None:
        pass

    def SetModel(self, rna_model, temp_C):
        my_model = Model(material=rna_model, celsius=int(temp_C))
        return my_model


    def GetPairProbs2DArray(self,mySequence, rna_model, temp_C ):
        my_model = self.SetModel(rna_model, temp_C)
        #convert into form the named touple is expecting
        #pairs = List[List[float]]
        nucpairs = list()
        pairsMatrix = pairs(strands=mySequence, model=my_model)
        pairsArray = pairsMatrix.to_array()
        for i in range(len(pairsArray)):
            nucpairs.append(list())
            for j in range(len(pairsArray[i])):            
                #nucpairs[i].append(list())
                pairValue = pairsArray[i][j]                    
                nucpairs[i].append(pairValue) 
        return nucpairs

    def process_ensemble_variation_make_groups(self, sequence:str, kcal_delta_span_from_mfe:int, Kcal_unit_increments: float, folded_2nd_state_structure:str='', target_2nd_state_structure:str=''):
        start_time=datetime.now()
        print(f'Starting test at {start_time}')
        print("Getting subopt\n")
        nucs_lists:List[List[str]]
        span_structures: Sara2StructureList
        span_structures = self.get_subopt_energy_gap(sequence_string=sequence, energy_delta_from_MFE=kcal_delta_span_from_mfe)       
        mfe_energy:float =  span_structures.mfe_freeEnergy
        
        
        print(f'Done with subopt gathering. {span_structures.num_structures} structures found\n')

        #this is for increments of 1 kcal need to do fraction
        num_groups: int = int(kcal_delta_span_from_mfe / Kcal_unit_increments)
        remainder: int = kcal_delta_span_from_mfe % Kcal_unit_increments

        groups_list : List[Sara2StructureList] = []
        groups_index_used: List[bool] = []
        groups_dict: Dict[int, Sara2StructureList] = {}
        group_values: List[float] = []



        #this fills up the list of energy deltas to publich EV's for
        current_energy: float = mfe_energy
        group_values.append(current_energy)
        for index in range(num_groups):
            current_energy = current_energy + Kcal_unit_increments
            group_values.append(current_energy)
        
        #if remainder > 0:
        #    current_energy = current_energy + kcal_delta_span_from_mfe
        #    group_values.append(current_energy)
        print(f'Processing group values {group_values} to \n')
        #now initialize the groups_list
        for index in range(len(group_values)-1):
            group: Sara2StructureList = Sara2StructureList()
            groups_list.append(group)
            groups_index_used.append(False)
            groups_dict[index+1] = group

        num_sara_struct: int = span_structures.num_structures
        for sara_index in range(0,num_sara_struct):
            #for sara_structure in span_structures.sara_stuctures:

            #this skips teh mfe from calulations
            sara_structure: Sara2SecondaryStructure = span_structures.sara_stuctures[sara_index]
            current_energy = sara_structure.freeEnergy

            #need to do this because there are two indexes need to look at each 
            #loop and want to avoid triggering a list index overrun
            for group_index in range(len(group_values)-1):
                #remember we are dealing with neg kcal so its you want to 
                min_energy: float = group_values[group_index]
                max_energy: float = group_values[group_index+1]
                if current_energy >= min_energy and current_energy < max_energy:
                    groups_list[group_index].add_structure(sara_structure)
                    groups_index_used[group_index] = True            
    

        for group_index in range(len(groups_list)):
            groups_dict[group_index] = groups_list[group_index]

    def get_subopt_energy_gap(self, sequence_string, energy_delta_from_MFE: int):
        #run through subopt 
        my_model = self.SetModel(rna_model, 37)
        print(f'Starting subopt at {datetime.now()}')
        kcal_group_structures_list: Sara2StructureList = Sara2StructureList()
        ensemble_kcal_group= subopt(strands=sequence_string, model=my_model, energy_gap=energy_delta_from_MFE)
        print(f'Completed subopt at {datetime.now()}')
        #get all the data out of it
        #list_of_nuc_lists: List[List[str]] = [[]]
        #num_nucs:int = len(sequence_string)
        #for index in range(num_nucs):
        #    temp_list:List[str] = []
        #    list_of_nuc_lists.append(temp_list)

        for i,kcal_group_elementInfo in enumerate(ensemble_kcal_group):
                  
            #get all the structures and energis pulled and prepped for proccessing and add them tot eh dict and the list               
            structure = str(kcal_group_elementInfo.structure)
            freeEnergy = float(kcal_group_elementInfo.energy)
            stackEnergy = float(kcal_group_elementInfo.stack_energy)
            structure_info: Sara2SecondaryStructure = Sara2SecondaryStructure(sequence=sequence_string, structure=structure, 
                                                                              freeEnergy=freeEnergy, stackEnergy=stackEnergy)
            kcal_group_structures_list.add_structure(structure_info)

            #now go throught everything
            #or index in range(num_nucs):
            #    character: str = structure[index]
            #    list_of_nuc_lists[index].append(character)


        return kcal_group_structures_list#, list_of_nuc_lists
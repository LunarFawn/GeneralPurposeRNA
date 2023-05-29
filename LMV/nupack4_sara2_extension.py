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

class MaterialParameter(Enum):
    NONE: 0

    #"Based on [Mathews99] and [Lu06] with additional parameters [Xia98,Zuker03] including coaxial stacking [Mathews99,Turner10] and dangle stacking [Serra95,Zuker03,Turner10] in 1M Na+."
    rna06_nupack4: 1
    #"Based on [Serra95] with additional parameters [Zuker03] including coaxial stacking [Mathews99,Turner10] and dangle stacking [Serra95,Zuker03,Turner10] in 1M Na+."
    rna95_nupack4: 2    
    #"Parameters from [Mathews99] with terminal mismatch free energies in exterior loops and multiloops replaced by two dangle stacking free energies. Parameters are provided only for 37 ∘C."
    rna99_nupack3: 3    
    #"Same as rna95 except that terminal mismatch free energies in exterior loops and multiloops are replaced by two dangle stacking free energies."
    rna95_nupack3: 4


class NUPACK4Interface():
    """
    Class for nupack4 interface for sara2 logic intended for serena package
    """

    def __init__(self) -> None:
        pass

    def select_material_parameters(self, parameters:MaterialParameter):
        param:str = ''
        
        match parameters :            
            case MaterialParameter.rna06_nupack4:
                param = "rna06"
            case MaterialParameter.rna95_nupack4:
                param = "rna95"
            case MaterialParameter.rna99_nupack3:
                param = "rna99-nupack3"
            case MaterialParameter.rna95_nupack3:
                param = "rna95-nupack3"
        return param                                                      

    def select_model(self, material_param:MaterialParameter, temp_C:int):
        param: str = self.select_material_parameters(material_param)
        my_model = Model(material=param, celsius=temp_C)
        return my_model

    def GetPairProbs2DArray(self, sequence, material_param:MaterialParameter, temp_C:int):
        param: str = self.select_material_parameters(material_param)
        my_model = self.select_model(param, temp_C)
        #convert into form the named touple is expecting
        #pairs = List[List[float]]
        nucpairs = list()
        pairsMatrix = pairs(strands=sequence, model=my_model)
        pairsArray = pairsMatrix.to_array()
        for i in range(len(pairsArray)):
            nucpairs.append(list())
            for j in range(len(pairsArray[i])):            
                #nucpairs[i].append(list())
                pairValue = pairsArray[i][j]                    
                nucpairs[i].append(pairValue) 
        return nucpairs

    def get_ensemble_groups(self, sequence:str, kcal_delta_span_from_mfe:int, Kcal_unit_increments: float, folded_2nd_state_structure:str='', folded_2nd_state_kcal:float=0):
        start_time=datetime.now()
        #print(f'Starting test at {start_time}')
        #print("Getting subopt\n")
        nucs_lists:List[List[str]]
        span_structures: Sara2StructureList
        span_structures = self.get_subopt_energy_gap(sequence_string=sequence, energy_delta_from_MFE=kcal_delta_span_from_mfe)       
        mfe_energy:float =  span_structures.mfe_freeEnergy
        
        
        #print(f'Done with subopt gathering. {span_structures.num_structures} structures found\n')

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
        #print(f'Processing group values {group_values} to \n')
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
    
        ensemble_groups: MultipleEnsembleGroups = MultipleEnsembleGroups()
        for group_index in range(len(groups_list)):
            this_group: SingleEnsembleGroup = SingleEnsembleGroup()
            this_group.group = groups_list[group_index]
            this_group.append_multi_state_mfe_data(span_structures.mfe_structure, span_structures.mfe_freeEnergy)
            this_group.append_multi_state_mfe_data(folded_2nd_state_structure, folded_2nd_state_kcal)
            this_group.kcal_span = Kcal_unit_increments
            this_group.kcal_start = group_values[group_index] - Kcal_unit_increments
            this_group.kcal_end = group_values[group_index]
            groups_dict[group_index] = groups_list[group_index]
            ensemble_groups.groups = groups_list
            ensemble_groups.raw_groups = groups_list
            ensemble_groups.groups_dict = groups_dict
            ensemble_groups.group_values = group_values
            
        return ensemble_groups
    
    def get_subopt_energy_gap(self, material_param:MaterialParameter, temp_C:int, sequence_string, energy_delta_from_MFE: int):
        #run through subopt
        param: str = self.select_material_parameters(material_param)
        my_model = self.select_model(param, temp_C)
        #print(f'Starting subopt at {datetime.now()}')
        kcal_group_structures_list: Sara2StructureList = Sara2StructureList()
        ensemble_kcal_group= subopt(strands=sequence_string, model=my_model, energy_gap=energy_delta_from_MFE)
        #print(f'Completed subopt at {datetime.now()}')
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
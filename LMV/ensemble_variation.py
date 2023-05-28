"""
Class that provides the ensemble variation algorithm 
used for local minima variation calculations
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


class SourceMFE(Enum):
    NONE: 0
    UNBOUND: 1
    BOUND: 2
    

@dataclass
class EV:
    ev_normalized: float = -1
    ev_ThresholdNorm: float = -1
    ev_structure: float = -1

@dataclass
class EVResult():
    groups_list : List[Sara2StructureList]
    groups_dict: Dict[int, Sara2StructureList]
    group_values: List[float]
    group_ev_list: List[EV]
    group_ev_dict: Dict[int,EV]

class LMV_Token():
    def __init__(self, num_groups: int) -> None:
        self._group_results: List[EV] = num_groups * [EV()]
        self._group_dict: Dict[int,EV] = {}
        self._group_values: List[str] = num_groups * ['']
        self._group_done_status: List[bool] = num_groups * [False]
    
    @property
    def group_dict(self):
        return self._group_dict
        
    def set_group_dict(self, index:int, value:EV):
        self._group_dict[index]=value

    @property
    def group_results(self):
        return self._group_results
        
    def set_group_result(self, index:int, value:EV):
        self._group_results[index]=value
    
    @property
    def group_values(self):
        return self._group_values
        
    def set_group_values(self, index:int, value:str):
        self._group_values[index]=value

    @property
    def group_done_status(self):
        return self._group_done_status
        
    def set_group_done_status(self, index:int, state:bool):
        self._group_done_status[index]=state
    
    @property
    def is_done(self):
        is_completed:bool = False
        if self._group_done_status.count(False) == 0:
            #its done
            is_completed = True
        return is_completed

class LMV_Shuttle():

    def __init__(self, structs_list: Sara2StructureList, mfe:Sara2SecondaryStructure, group_index:int, token:LMV_Token) -> None:
        self._kcal_group_structures_list: Sara2StructureList = structs_list
        self._sara_mfestructure:Sara2SecondaryStructure = mfe
        self._group_index:int = group_index
        self._token:LMV_Token = token
    
    @property
    def kcal_group_structures_list(self):
        return self._kcal_group_structures_list

    @kcal_group_structures_list.setter
    def kcal_group_structures_list(self, new_list: Sara2StructureList):
        self._kcal_group_structures_list = new_list
    
    @property
    def sara_mfestructure(self):
        return self._sara_mfestructure

    @sara_mfestructure.setter
    def sara_mfestructure(self, new_strucr: Sara2SecondaryStructure):
        self._sara_mfestructure = new_strucr
    
    @property
    def group_index(self):
        return self._group_index

    @group_index.setter
    def group_index(self, new_index: int):
        self._group_index = new_index
    
    @property
    def token(self):
        return self._token

    @token.setter
    def token(self, new_token: LMV_Token):
        self._token = new_token

class LMV_ThreadProcessor():
    
    def __init__(self, stuctures: List[Sara2StructureList], source_mfe: SourceMFE) -> None:
        self._sara2_groups: List[Sara2StructureList] = stuctures
        self._source_mfe: SourceMFE= source_mfe
        num_groups:int = len(stuctures)
        self._num_groups: int =  num_groups
        self._group_token: LMV_Token = LMV_Token(num_groups)
        self._LMV: EnsembleVariation = EnsembleVariation()
    
    @property
    def sara2_groups(self):
        return self._sara2_groups

    @sara2_groups.setter
    def sara2_groups(self, new_list:List[Sara2StructureList]):
        self._sara2_groups = new_list
    
    @property
    def source_mfe(self):
        return self._source_mfe

    @source_mfe.setter
    def source_mfe(self, new_struct:Sara2SecondaryStructure):
        self._source_mfe = new_struct

    @property
    def num_groups(self):
        return self._num_groups

    @num_groups.setter
    def num_groups(self, new_num:int):
        self._num_groups = new_num

    @property
    def group_token(self):
        return self._group_token

    @group_token.setter
    def group_token(self, new_token:LMV_Token):
        self._group_token = new_token
    
    @property
    def LMV(self):
        return self._LMV

    @LMV.setter
    def LMV(self, new_lmv:EnsembleVariation):
        self._LMV = new_lmv

    def run_LMV(self):
        self.start_calculations()
        self.wait_for_finish()
        #the test should be done now
        #check foor index that is -1 and if so then use prev value
        num_groups:int = len(self.group_token.group_results)
        for index in range(1, num_groups):
            if self.group_token.group_results[index].ev_normalized == -1:
                previous_EV = self.group_token.group_results[index-1]
                self.group_token.group_results[index] = previous_EV
                self.group_token.group_dict[index] = previous_EV
        return self.group_token

    def start_calculations(self):
        for thread_index in range(self.num_groups):
            sara2_structs: Sara2StructureList  = self.sara2_groups[thread_index]
            temp_source_mfe:Sara2SecondaryStructure
            if len(sara2_structs.sara_stuctures) == 0:
                #use mfe as its always there
                temp_source_mfe = self.sara2_groups[0].sara_stuctures[0]
            else:
                if self.source_mfe == SourceMFE.BOUND:
                    temp_source_mfe = sara2_structs.sara_stuctures[0]
                elif self.source_mfe == SourceMFE.UNBOUND:
                    temp_source_mfe = self.sara2_groups[0].sara_stuctures[0]
                else:
                    temp_source_mfe = self.source_mfe
            new_shuttle: LMV_Shuttle = LMV_Shuttle(structs_list=sara2_structs, mfe=temp_source_mfe, group_index=thread_index,token=self.group_token) 
            mew_thread = threading.Thread(target=self.LMV.thread_EV, args=[new_shuttle])
            mew_thread.start()

    
    def wait_for_finish(self):
                
        stop:bool = False
        while stop == False:
            print(f'Checking LMV status at {datetime.now()}')
            current_status: List[bool] = self.group_token.group_done_status
            is_done:bool = self.group_token.is_done
            
            message: str = ''
            for index in range(self.num_groups):
                goup_value:str = self.group_token.group_values[index]
                done_status: bool = self.group_token.group_done_status[index]
                message = message + f'Group_{index+1}: kcal_group={goup_value}, status={done_status}\n'
            print(message)

            if is_done == True:
                stop = True
                print(f'Its done at {datetime.now()}')
            else:
                dwell_time:int = 5
                print(f'dwelling for {dwell_time} seconds until next check')
                time.sleep(dwell_time)



class EnsembleVariation():
    """
    Class that exposes the algorithm
    """

    def __init__(self) -> None:
        pass

    def get_ensemble_variation(self, ensemble: List[Sara2StructureList], state_source:int = 1):        
        #now process all the groups
       
        source_mfe: SourceMFE = SourceMFE.UNBOUND

        
        
        #should be able to populate before hand and add to the group stuff i am working on
        #seams to push for the need a bit more
        groups_dict: Dict[int, Sara2StructureList] = {}
        group_values: List[float] = []

        if state_source == 1:    
            print(f'Begining LMV_U processing at {datetime.now()}')
        elif state_source == 2:
            source_mfe = SourceMFE.BOUND
            print(f'Begining LMV_UB processing at {datetime.now()}')
        else:
            message :str =  "Only supports 2 states right now or didi you mean 1st state if 0 entered. 1 based index please remember for this one."
            raise Exception(message)
        
        #single_ensemble_group: List[Sara2StructureList] = [ensemble.group]
        LMV_Thread: LMV_ThreadProcessor = LMV_ThreadProcessor(stuctures=ensemble,source_mfe=source_mfe)
        result_thread_LMV:LMV_Token = LMV_Thread.run_LMV()
        group_ev_list: List[EV] = result_thread_LMV.group_results
        group_ev_dict: Dict[int,EV] = result_thread_LMV.group_dict

        result_LMV_U_rel: EVResult = EVResult(groups_list=ensemble, groups_dict=groups_dict, 
                                              group_values=group_values, group_ev_list=group_ev_list, 
                                              group_ev_dict=group_ev_dict)

    def thread_EV(self, shuttle: LMV_Shuttle):
        total_EV_subscore1:int = 0
        structure_element_count = shuttle.kcal_group_structures_list.num_structures
        token:LMV_Token = shuttle.token 
        group_num:int = shuttle.group_index
        if structure_element_count != 0:
            kcal_group_structures_list: Sara2StructureList = shuttle.kcal_group_structures_list

            sara_mfestructure:Sara2SecondaryStructure = shuttle.sara_mfestructure 
           
            
            #need to do each char abd then structure
            #walk through each nucleotide but first prep containers grab what is needed
            
            #setup constants
            nuc_count = kcal_group_structures_list.nuc_count
            structure_element_count = kcal_group_structures_list.num_structures

            ensembleVariation_score_temp = 0
            nucleotide_position_variation_basescores=[0]*nuc_count
            nucleotide_position_variation_subscores=[0]*nuc_count
            energydelta_individualVariationScore_list=[]
            
            #add the step to get nuc array here
            #get all the data out of it

            #first initialize the lists
            list_of_nuc_lists: List[List[str]] = []
    
            num_nucs: int = kcal_group_structures_list.nuc_count
            for index in range(num_nucs):
                temp_list:List[str] = []
                list_of_nuc_lists.append(temp_list)
                
            
            #now go throught everything
            for sara_structure in kcal_group_structures_list.sara_stuctures:
                for index in range(num_nucs):
                    character: str = sara_structure.structure[index]
                    list_of_nuc_lists[index].append(character)

            list_of_nuc_scores_base: List[int] = [0]*nuc_count
            list_of_nuc_scores_subscores: List[int] = [0]*nuc_count
            num_structs:int = kcal_group_structures_list.num_structures

            for nucIndex in range(nuc_count):
                mfe_nuc=sara_mfestructure.structure[nucIndex]
                num_chars = list_of_nuc_lists[nucIndex].count(mfe_nuc)
                num_diff:int = num_structs - num_chars
                list_of_nuc_scores_base[nucIndex] = num_diff
                list_of_nuc_scores_subscores[nucIndex] = list_of_nuc_scores_base[nucIndex] / structure_element_count 
            
            total_EV_subscore1 = sum(list_of_nuc_scores_subscores)
        else:
            total_EV_subscore1 = -1

        result: EV =  EV(ev_normalized=total_EV_subscore1, ev_ThresholdNorm=0, ev_structure=0)  
        token.group_results[group_num]= result
        token.group_dict[group_num] = result
        token.group_done_status[group_num] = True


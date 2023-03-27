"""
Sara2 api for nupack
copyright 2023 Jennifer Pearl
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

my_model = Model

rna_model='rna95-nupack3'    
# Define physical model 

class Sara2SecondaryStructure(object):

    def __init__(self, sequence:str = '', structure: str = '', freeEnergy: float = 0, stackEnergy: float = 0) -> None:
        self._sequence: str = sequence
        self._structure: str = structure
        self._freeEnergy: float = freeEnergy
        self._stackEnergy: float = stackEnergy
        #self._nuc_count: int = len(sequence)

    @property
    def sequence(self):
        return self._sequence

    @sequence.setter
    def sequence(self, primary_struc: str):
        self._sequence = primary_struc

    @property
    def structure(self):
        return self._structure

    @structure.setter
    def structure(self, dot_parens: str):
        self._structure = dot_parens
    
    @property
    def freeEnergy(self):
        return self._freeEnergy

    @freeEnergy.setter
    def freeEnergy(self, energy: float):
        self._freeEnergy = energy
    
    @property
    def stackEnergy(self):
        return self._stackEnergy

    @stackEnergy.setter
    def stackEnergy(self, energy: float):
        self._stackEnergy = energy
    
    @property
    def nuc_count(self):
        return len(self._sequence)



class Sara2StructureList(object):
    
    def __init__(self) -> None:
        self._sara_structures_list: List[Sara2SecondaryStructure] = []
        self._structures: List[str] = []
        self._freeEnergy_list: list[float] = []
        self._stackEnergy_list: list[float] = []
        self._min_freeEnergy: float = 0
        self._max_freeEnergy: float = 0
        self._min_stackEnergy: float = 0
        self._max_stackEnergy: float = 0
        self._num_structures: int = 0
        self._nuc_count: int = 0
        self._mfe_structure: str = ''
        self._mfe_freeEnergy: float = 0
        self._mfe_stackEnergy: float = 0
        self._freeEnergy_span:float = 0
        self._stackEnergy_span:float = 0

    def process_energy(self):
            #now populate min and max
        #do free energy
        if len(self._freeEnergy_list) == 0:
            self._min_freeEnergy = 0
            self._max_freeEnergy = 0
        else:
            self._min_freeEnergy = min(self._freeEnergy_list)
            self._max_freeEnergy = max(self._freeEnergy_list)

        self._freeEnergy_span = self._max_freeEnergy - self._min_freeEnergy
        #do stack energy

        if len(self._stackEnergy_list) == 0:
            self._min_stackEnergy = 0
            self._max_stackEnergy = 0
        else:
            self._min_stackEnergy = min(self._stackEnergy_list)
            self._max_stackEnergy = max(self._stackEnergy_list)
        self._stackEnergy_span = self._max_stackEnergy - self._min_stackEnergy

        #now count
        if len(self._structures) == 0:
            self._num_structures = 0
        else:
            self._num_structures = len(self._structures)

    def add_structure(self, structure: Sara2SecondaryStructure):
        self._sara_structures_list.append(structure)
        self._structures.append(structure.structure)
        self._freeEnergy_list.append(structure.freeEnergy)
        self._stackEnergy_list.append(structure.stackEnergy)
        #self.process_energy()


    def remove_structure(self, index:int):
        del self._structures[index]
        del self._freeEnergy_list[index]
        del self._stackEnergy_list[index]
        #self.process_energy()            


    @property
    def mfe_structure(self):
        return self.sara_stuctures[0].sequence

    @property
    def mfe_freeEnergy(self):
        return self.sara_stuctures[0].freeEnergy
    
    @property
    def mfe_stackEnergy(self):
        return self.sara_stuctures[0].stackEnergy
    
    @property
    def nuc_count(self):
        return self.sara_stuctures[0].nuc_count 

    @property
    def sara_stuctures(self):
        return self._sara_structures_list

    @sara_stuctures.setter
    def sara_stuctures(self, structs_list: List[Sara2SecondaryStructure]):
        #reset list
        self._sara_structures_list=[]
        #fill it in now
        for struc in structs_list:
            self.add_structure(struc)

    
    @property
    def max_free_energy(self):
        self.process_energy()
        return self._max_freeEnergy
    
    @property
    def min_free_energy(self):
        self.process_energy()
        return self._min_freeEnergy
    
    @property
    def max_stack_energy(self):
        self.process_energy()
        return self._max_stackEnergy
    
    @property
    def min_stack_energy(self):
        self.process_energy()
        return self._min_stackEnergy
    
    @property
    def num_structures(self):
        self.process_energy()
        return self._num_structures
    
    @property
    def freeEnergy_span(self):
        self.process_energy()
        return self._freeEnergy_span

    @property
    def stackEnergy_span(self):
        self.process_energy()
        return self._stackEnergy_span 


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



def SetModel(rna_model, temp_C):
    my_model = Model(material=rna_model, celsius=int(temp_C))
    return my_model


def GetPairProbs2DArray(mySequence, rna_model, temp_C ):
    my_model = SetModel(rna_model, temp_C)
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



class EnsembleVariation:

    def __init__(self) -> None:
        pass    

    def GetEnsembleVariation(self):
        #process ensemble variation
        pass

    def process_ensemble_variation(self, sequence:str, kcal_delta_span_from_mfe:int, Kcal_unit_increments: float, folded_2nd_state_structure:str='', target_2nd_state_structure:str=''):
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



        #now process all the groups
        print(f'Begining LMV_U processing at {datetime.now()}')
        LMV_U_thread_mfe: LMV_ThreadProcessor = LMV_ThreadProcessor(stuctures=groups_list,mfe_stuct="mfe")
        result_thread_LMV_mfe:LMV_Token = LMV_U_thread_mfe.run_LMV()

        print(f'Begining LMV_U rel processing at {datetime.now()}')
        LMV_U_thread_rel: LMV_ThreadProcessor = LMV_ThreadProcessor(stuctures=groups_list,mfe_stuct="rel")
        result_thread_LMV_rel:LMV_Token = LMV_U_thread_rel.run_LMV()

        #for list_index in range(len(groups_list)-1):
        #    struct_list: Sara2StructureList = groups_list[list_index]
        #    print(f'Processing {group_values[list_index]} to {group_values[list_index+1]}')
        #    print(f'started LMV_U group EV at {datetime.now()}')
        #    ev: EV = self.advanced_EV(struct_list, struct_list.sara_stuctures[0])
        #    print(f'finished LMV_U group EV at {datetime.now()}\n')
        #    group_ev_list.append(ev)
        #    group_ev_dict[list_index+1] = ev
        
       
        switch_mfe_sub: Sara2SecondaryStructure = Sara2SecondaryStructure()
        switch_mfe_sub.sequence=span_structures.sara_stuctures[0].sequence
        switch_mfe_sub.structure=target_2nd_state_structure
        switch_mfe_sub.freeEnergy=span_structures.sara_stuctures[0].freeEnergy
        switch_mfe_sub.stackEnergy=span_structures.sara_stuctures[0].stackEnergy

        LMV_US_thread_target: LMV_ThreadProcessor = LMV_ThreadProcessor(stuctures=groups_list,mfe_stuct=switch_mfe_sub)
        

        if target_2nd_state_structure != '':
            result_thread_LMV_US:LMV_Token = LMV_US_thread_target.run_LMV()
            #now process all the groups
            #print(f'Begining LMV_US processing at {datetime.now()}')
            #for list_index in range(len(groups_list)-1):
            #    struct_list: Sara2StructureList = groups_list[list_index]
            #    print(f'Processing  {group_values[list_index]} to {group_values[list_index+1]}')
            #    print(f'started LMV_US group EV at {datetime.now()}')
            #    ev: EV = self.advanced_EV(struct_list, switch_mfe_sub)
            #    print(f'finished LMV_US group EV at {datetime.now()}\n')
            #    switch_ev_list.append(ev)
            #    switch_ev_dict[list_index+1] = ev

        
        switch_mfe_sub_folded: Sara2SecondaryStructure = Sara2SecondaryStructure()
        switch_mfe_sub_folded.sequence=span_structures.sara_stuctures[0].sequence
        switch_mfe_sub_folded.structure=folded_2nd_state_structure
        switch_mfe_sub_folded.freeEnergy=span_structures.sara_stuctures[0].freeEnergy
        switch_mfe_sub_folded.stackEnergy=span_structures.sara_stuctures[0].stackEnergy

        LMV_US_thread_folded: LMV_ThreadProcessor = LMV_ThreadProcessor(stuctures=groups_list,mfe_stuct=switch_mfe_sub_folded)

        if folded_2nd_state_structure != '':
            result_thread_LMV_US_folded:LMV_Token = LMV_US_thread_folded.run_LMV()
            
        finish_time = datetime.now()
        print(f'Started test at {start_time}')
        print(f'Completed test at {finish_time}')
        timelenght:timedelta = finish_time-start_time
        print(f'Total Time (seconds) = {timelenght.total_seconds()}')

        group_ev_list_mfe: List[EV] = result_thread_LMV_mfe.group_results
        group_ev_dict_mfe: Dict[int,EV] = result_thread_LMV_mfe.group_dict

        group_ev_list_rel: List[EV] = result_thread_LMV_rel.group_results
        group_ev_dict_rel: Dict[int,EV] = result_thread_LMV_rel.group_dict

        switch_ev_list_target: List[EV] = result_thread_LMV_US.group_results
        switch_ev_dict_target: Dict[int,EV] = result_thread_LMV_US.group_dict

        switch_ev_list_folded: List[EV] = result_thread_LMV_US_folded.group_results
        switch_ev_dict_folded: Dict[int,EV] = result_thread_LMV_US_folded.group_dict


        result_LMV_U_mfe: EVResult = EVResult(groups_list=groups_list, groups_dict=groups_dict, group_values=group_values, group_ev_list=group_ev_list_mfe, group_ev_dict=group_ev_dict_mfe)
        result_LMV_U_rel: EVResult = EVResult(groups_list=groups_list, groups_dict=groups_dict, group_values=group_values, group_ev_list=group_ev_list_rel, group_ev_dict=group_ev_dict_rel)
        result_LMV_US_target: EVResult = EVResult(groups_list=groups_list, groups_dict=groups_dict, group_values=group_values, group_ev_list=switch_ev_list_target, group_ev_dict=switch_ev_dict_target)
        result_LMV_US_folded: EVResult = EVResult(groups_list=groups_list, groups_dict=groups_dict, group_values=group_values, group_ev_list=switch_ev_list_folded, group_ev_dict=switch_ev_dict_folded)
        
        return result_LMV_U_mfe, result_LMV_U_rel, result_LMV_US_target, result_LMV_US_folded


    def get_subopt_energy_gap(self, sequence_string, energy_delta_from_MFE: int):
        #run through subopt 
        my_model = SetModel(rna_model, 37)
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

    
    def advanced_EV(self, kcal_group_structures_list: Sara2StructureList, sara_mfestructure:Sara2SecondaryStructure):
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
        result: EV =  EV(ev_normalized=total_EV_subscore1, ev_ThresholdNorm=0, ev_structure=0)   
     
        
        """

        #we have all the elements and all the stuff for the group
        for nucIndex in range(nuc_count):
            for sara_structure in kcal_group_structures_list.sara_stuctures:
                #first do the EV stuff
                mfe_nuc=sara_mfestructure.structure[nucIndex]
                delta_nuc = sara_structure.structure[nucIndex]
                if mfe_nuc != delta_nuc:
                    nucleotide_position_variation_basescores[nucIndex]=nucleotide_position_variation_basescores[nucIndex]+1
                
                structure = sara_structure.structure
                mfe_structure = sara_mfestructure.structure
                element_structure_mfe_distance=struc_distance(structure, mfe_structure)
                individualVariationScore = element_structure_mfe_distance / structure_element_count
                energydelta_individualVariationScore_list.append(individualVariationScore)

            #when done with all teh structures then do the nucScore wityh calculations for EV
            nucleotide_position_variation_subscores[nucIndex]=nucleotide_position_variation_basescores[nucIndex] / structure_element_count
        #when that is all done the calculate final score
        total_subscore=sum(nucleotide_position_variation_subscores)
        structure_ensemble_variance = sum(energydelta_individualVariationScore_list)# * 100
        
        #now do final calculations
        # original classic EV and normalized classic EV        
        nucleotide_ensemble_variance_normalized=total_subscore #/nuc_count    
        fail_threshold=0.2
        nucleotide_ensemble_variance_ThresholdNorm=(total_subscore/fail_threshold)#*100
        
        result2: EV =  EV(ev_normalized=nucleotide_ensemble_variance_normalized, 
                                                            ev_ThresholdNorm=nucleotide_ensemble_variance_ThresholdNorm, 
                                                            ev_structure=structure_ensemble_variance)   
        """   
        return  result

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

 

class LMV_ThreadProcessor():
    
    def __init__(self, stuctures: List[Sara2StructureList],mfe_stuct: Sara2SecondaryStructure) -> None:
        self._sara2_groups: List[Sara2StructureList] = stuctures
        self._mfe_stuct: Sara2SecondaryStructure= mfe_stuct
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
    def mfe_stuct(self):
        return self._mfe_stuct

    @mfe_stuct.setter
    def mfe_stuct(self, new_struct:Sara2SecondaryStructure):
        self._mfe_stuct = new_struct

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
            temp_mfe_stuct:Sara2SecondaryStructure
            if len(sara2_structs.sara_stuctures) == 0:
                #use mfe as its always there
                temp_mfe_stuct = self.sara2_groups[0].sara_stuctures[0]
            else:
                if self.mfe_stuct == "rel":
                    temp_mfe_stuct = sara2_structs.sara_stuctures[0]
                elif self.mfe_stuct == "mfe":
                    temp_mfe_stuct = self.sara2_groups[0].sara_stuctures[0]
                else:
                    temp_mfe_stuct = self.mfe_stuct
            new_shuttle: LMV_Shuttle = LMV_Shuttle(structs_list=sara2_structs, mfe=temp_mfe_stuct, group_index=thread_index,token=self.group_token) 
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
        



from typing import List, Dict
import struct
from nupack import *
import pandas as pd
import sys
import openpyxl
from copy import deepcopy
from dataclasses import dataclass

from Sara2_dev.nupackAPI_Sara2_Ver1 import Sara2SecondaryStructure, Sara2StructureList
import Sara2_dev.nupackAPI_Sara2_Ver1 as nupackAPI




@dataclass
class EV:
    ev_normalized: float = -1
    ev_ThresholdNorm: float = -1
    ev_structure: float = -1

@dataclass
class EVGroup:
    group_energy_span: float = 0
    group_starting_energy: float = 0
    ev: EV = EV()


class EVGroupResult(object):
    ev_info: EVGroup = EVGroup()
    #ev_structure_scores: Dict[Sara2StructureList, EV] = {}
    ev_group_struct_list: Sara2StructureList     

class EVResulFull(object):
    start_energy_list: List[float] = []
    #this is intended to have the float be the min mfe energy which defines teh energy group same as an index would
    result_dict: Dict[float, EVGroupResult] = {}

class EnsembleVariation():

    def __init__(self, rna_model: Model) -> None:
        self._my_model: Model = rna_model 


    def process_ensemble_variation(self, sequence:str, kcal_delta_span_from_mfe:int, Kcal_unit_increments: int):
        span_structures: Sara2StructureList = self.get_subopt_energy_gap(sequence_string=sequence, energy_delta_from_MFE=kcal_delta_span_from_mfe)       
        mfe_energy:float =  span_structures.mfe_freeEnergy

        num_groups: int = kcal_delta_span_from_mfe / Kcal_unit_increments
        remainder: int = kcal_delta_span_from_mfe % Kcal_unit_increments

        

        #groups_list : List[Sara2StructureList] = []
        #groups_dict: Dict[int, Sara2StructureList] = {}
        
        #group_ev_list: List[EV] = []
        #group_ev_dict: Dict[int,EV] = {}

        results_list: EVResulFull = EVResulFull()

        group_values: List[float] = []
        #this fills up the list of energy deltas to publich EV's for
        current_energy: float = mfe_energy
        for index in num_groups:
            num = current_energy + Kcal_unit_increments
            group_values.append(num)
        
        if remainder > 0:
            num = mfe_energy + kcal_delta_span_from_mfe
            group_values.append(num)

        #now initialize the groups_list
        for value in group_values:
            group: EVGroupResult = EVGroupResult()
            results_list.result_dict[value] = group
            results_list.start_energy_list.append(value)


        for structure in span_structures.stuctures:
            current_energy = structure.freeEnergy

            #need to do this because there are two indexes need to look at each 
            #loop and want to avoid triggering a list index overrun
            for group_index in range(len(group_values)-1):
                #first find the structure
                if current_energy in range(group_values[group_index], group_values[group_index+1]):
                    #now get the EV
                    results_list.result_dict[group_values[group_index]].ev_group_struct_list.add_structure(structure)
                    #groups_list[group_index].add_structure(structure)   
                     
                
        for group_index in range(len(group_values)):
            structure_list = results_list.result_dict[group_values[group_index]].ev_group_struct_list
            #groups_dict[group_values[group_index]].stuctures = groups_list[group_index]
            #struct_list: Sara2StructureList = groups_list[group_index]
            energy_span = group_values[group_index+1]-group_values[group_index]
            ev: EV = self.advanced_EV(structure_list, structure_list.mfe_structure)
            results_list.result_dict[group_values[group_index]].ev_info.ev = ev
            results_list.result_dict[group_values[group_index]].ev_info.group_energy_span = energy_span
            results_list.result_dict[group_values[group_index]].ev_info.group_starting_energy = group_values[group_index]
            #group_ev_list.append(ev)
            #group_ev_dict[group_values[group_index]] = ev

        #now process all the groups
        #for list_index in range(len(groups_list)):


        #result: EVResult = EVResult(groups_list=groups_list, groups_dict=groups_dict, group_values=group_values, group_ev_list=group_ev_list, group_ev_dict=group_ev_dict)
        
        results_list.start_energy_list.sort()
        
        return results_list


    def get_subopt_energy_gap(self, sequence_string, energy_delta_from_MFE: int):
        #run through subopt 
        kcal_group_structures_list: Sara2StructureList = Sara2StructureList()
        ensemble_kcal_group= subopt(strands=sequence_string, model=self._my_model, energy_gap=energy_delta_from_MFE)
        
        #get all the data out of it
        for i,kcal_group_elementInfo in enumerate(ensemble_kcal_group):        
            #get all the structures and energis pulled and prepped for proccessing and add them tot eh dict and the list               
            structure = str(kcal_group_elementInfo.structure)
            freeEnergy = float(kcal_group_elementInfo.energy)
            stackEnergy = float(kcal_group_elementInfo.stack_energy)

            structure_info: Sara2SecondaryStructure = Sara2SecondaryStructure(sequence=sequence_string, structure=structure, 
                                                                              freeEnergy=freeEnergy, stackEnergy=stackEnergy)
            kcal_group_structures_list.add_structure(structure_info)
        return kcal_group_structures_list

    
    def advanced_EV(self, kcal_group_structures_list: Sara2StructureList, mfestructure:Sara2SecondaryStructure):
        #need to do each char abd then structure
        #walk through each nucleotide but first prep containers grab what is needed
        
        #setup constants
        nuc_count = kcal_group_structures_list.nuc_count
        structure_element_count = kcal_group_structures_list.num_structures

        ensembleVariation_score_temp = 0
        nucleotide_position_variation_basescores=[0]*nuc_count
        nucleotide_position_variation_subscores=[0]*nuc_count
        energydelta_individualVariationScore_list=[]
        
        #we have all the elements and all the stuff for the group
        for nucIndex in range(nuc_count):
            for structure in kcal_group_structures_list.stuctures:
                #first do the EV stuff
                mfe_nuc=mfestructure.structure[nucIndex]
                delta_nuc = structure.structure[nucIndex]
                if mfe_nuc != delta_nuc:
                    nucleotide_position_variation_basescores[nucIndex]=nucleotide_position_variation_basescores[nucIndex]+1

                element_structure_mfe_distance=struc_distance(structure, mfestructure)
                individualVariationScore = element_structure_mfe_distance / structure_element_count
                energydelta_individualVariationScore_list.append(individualVariationScore)

            #when done with all teh structures then do the nucScore wityh calculations for EV
            nucleotide_position_variation_subscores[nucIndex]=nucleotide_position_variation_basescores[nucIndex] / structure_element_count
        #when that is all done the calculate final score
        total_subscore=sum(nucleotide_position_variation_subscores)
        structure_ensemble_variance = sum(energydelta_individualVariationScore_list)# * 100
        
        #now do final calculations
        # original classic EV and normalized classic EV        
        nucleotide_ensemble_variance_normalized=total_subscore * 100 #/nuc_count    
        fail_threshold=0.2
        nucleotide_ensemble_variance_ThresholdNorm=(total_subscore/fail_threshold)#*100
        result: EV =  EV(ev_normalized=nucleotide_ensemble_variance_normalized, 
                                                            ev_ThresholdNorm=nucleotide_ensemble_variance_ThresholdNorm, 
                                                            ev_structure=structure_ensemble_variance)      
        return  result


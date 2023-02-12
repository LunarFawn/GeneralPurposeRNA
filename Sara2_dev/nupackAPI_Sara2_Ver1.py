from typing import List, Dict
import struct
from nupack import *
import pandas as pd
import sys
import openpyxl
from copy import deepcopy
from dataclasses import dataclass

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
        self._structures: List[Sara2SecondaryStructure] = []
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
        self._min_freeEnergy = min(self._freeEnergy_list)
        self._max_freeEnergy = max(self._freeEnergy_list)
        self._freeEnergy_span = self._max_freeEnergy - self._min_freeEnergy
        #do stack energy
        self._min_stackEnergy = min(self._stackEnergy_list)
        self._max_stackEnergy = max(self._stackEnergy_list)
        self._stackEnergy_span = self._max_stackEnergy - self._min_stackEnergy

        #now count
        self._num_structures = len(self._structures)

    def add_structure(self, structure: Sara2SecondaryStructure):
        self._structures.append(structure.structure)
        self._freeEnergy_list.append(structure.freeEnergy)
        self._stackEnergy_list.append(structure.stackEnergy)
        self.process_energy()


    def remove_structure(self, index:int):
        del self._structures[index]
        del self._freeEnergy_list[index]
        del self._stackEnergy_list[index]
        self.process_energy()            


    @property
    def mfe_structure(self):
        return self._structures[0].sequence

    @property
    def mfe_freeEnergy(self):
        return self._structures[0].freeEnergy
    
    @property
    def mfe_stackEnergy(self):
        return self._structures[0].stackEnergy
    
    @property
    def nuc_count(self):
        return self._structures[0].nuc_count 

    @property
    def stuctures(self):
        return self._structures

    @stuctures.setter
    def stuctures(self, structs_list: List[Sara2SecondaryStructure]):
        self._structures =  deepcopy(structs_list)
        #now populate min and max
        #do free energy
        self.process_energy()
    
    @property
    def max_free_energy(self):
        return self._max_freeEnergy
    
    @property
    def min_free_energy(self):
        return self._min_freeEnergy
    
    @property
    def max_stack_energy(self):
        return self._max_stackEnergy
    
    @property
    def min_stack_energy(self):
        return self._min_stackEnergy
    
    @property
    def num_structures(self):
        return self._num_structures
    
    @property
    def freeEnergy_span(self):
        return self._freeEnergy_span

    @property
    def stackEnergy_span(self):
        return self._stackEnergy_span 



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


@dataclass
class EV:
    ev_normalized: float = -1
    ev_ThresholdNorm: float = -1
    ev_structure: float = -1


@dataclass
class EVResult:
    groups_list : List[Sara2StructureList] = []
    groups_dict: Dict[int, Sara2StructureList] = {}
    group_values: List[float] = []
    group_ev_list: List[EV] = []
    group_ev_dict: Dict[int,EV] = {}


class EnsembleVariation:

    def __init__(self) -> None:
        pass    

    def GetEnsembleVariation(self):
        #process ensemble variation
        pass

    def process_ensemble_variation(self, sequence:str, kcal_delta_span_from_mfe:int, Kcal_unit_increments: int):
        span_structures: Sara2StructureList = self.get_subopt_energy_gap(sequence_string=sequence, energy_delta_from_MFE=kcal_delta_span_from_mfe)       
        mfe_energy:float =  span_structures.mfe_freeEnergy

        num_groups: int = kcal_delta_span_from_mfe / Kcal_unit_increments
        remainder: int = kcal_delta_span_from_mfe % Kcal_unit_increments

        

        groups_list : List[Sara2StructureList] = []
        groups_dict: Dict[int, Sara2StructureList] = {}
        group_values: List[float] = []
        group_ev_list: List[EV] = []
        group_ev_dict: Dict[int,EV] = {}

        #this fills up the list of energy deltas to publich EV's for
        current_energy: float = mfe_energy
        for index in num_groups:
            num = current_energy + Kcal_unit_increments
            group_values.append(num)
        
        if remainder > 0:
            num = mfe_energy + kcal_delta_span_from_mfe
            group_values.append(num)

        #now initialize the groups_list
        for index in range(len(group_values)):
            group: Sara2StructureList = Sara2StructureList()
            groups_list.append(group)
            groups_dict[index+1] = group


        for structure in span_structures.stuctures:
            current_energy = structure.freeEnergy

            #need to do this because there are two indexes need to look at each 
            #loop and want to avoid triggering a list index overrun
            for group_index in range(len(group_values)-1):
                if current_energy in range(group_values[group_index], group_values[group_index+1]):
                    groups_list[group_index].add_structure(structure)        
                
        for group_index in range(len(groups_list)):
            groups_dict[group_index+1].stuctures = groups_list[group_index]

        #now process all the groups
        for list_index in range(len(groups_list)):
            struct_list: Sara2StructureList = groups_list[list_index]
            ev: EV = self.advanced_EV(struct_list, struct_list.mfe_structure)
            group_ev_list.append(ev)
            group_ev_dict[list_index+1] = ev

        result: EVResult = EVResult(groups_list=groups_list, groups_dict=groups_dict, group_values=group_values, group_ev_list=group_ev_list, group_ev_dict=group_ev_dict)
        return result


    def get_subopt_energy_gap(self, sequence_string, energy_delta_from_MFE: int):
        #run through subopt 
        kcal_group_structures_list: Sara2StructureList = Sara2StructureList()
        ensemble_kcal_group= subopt(strands=sequence_string, model=my_model, energy_gap=energy_delta_from_MFE)
        
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


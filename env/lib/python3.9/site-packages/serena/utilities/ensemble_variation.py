"""
File for the ensemble variation code to live
"""
from typing import List, Dict
from dataclasses import dataclass

from serena.utilities.ensemble_structures import Sara2SecondaryStructure, Sara2StructureList

@dataclass
class EV:
    """
    Returns the different version of enemble variation
    """
    ev_normalized: float = -1
    ev_threshold_norm: float = -1
    ev_structure: float = -1

@dataclass
class EVResult():
    """
    Class for holding the ev values for ensemble groups 
    """
    ev_values:List[EV]

class EVToken():
    """
    Class for the token that is used to pass data between thread
    and caller to record ev's and comlete flags
    """
    def __init__(self, num_groups: int) -> None:
        self._group_results: List[EV] = num_groups * [EV()]
        self._group_dict: Dict[int,EV] = {}
        self._group_values: List[str] = num_groups * ['']
        self._group_done_status: List[bool] = num_groups * [False]

    @property
    def group_dict(self)->Dict[int,EV]:
        """
        Return the group dictionary of ev's for index
        will be depreciated soon probably
        """
        return self._group_dict

    def set_group_dict(self, index:int, value:EV):
        """
        Sets the group dictionary of ev's with index
        will be depreciated soon probably
        """
        self._group_dict[index]=value

    @property
    def group_results(self)->List[EV]:
        """
        Return the ev results for the groups as a list of EV's
        """
        return self._group_results

    @property
    def ev_results(self) -> EVResult:
        """
        Return the ev results for the groups as a EVREsult
        """
        result: EVResult = EVResult(ev_values=self.group_results)
        return result

    def set_group_result(self, index:int, value:EV):
        """
        Set the ev group results for the groups
        """
        self._group_results[index]=value

    @property
    def group_values(self)->List[str]:
        """
        Return the values of the energy groups as a list of str
        """
        return self._group_values

    def set_group_values(self, index:int, value:str):
        """
        Set the values of the energy groups as a list of str by index
        """
        self._group_values[index]=value

    @property
    def group_done_status(self)->List[bool]:
        """
        Return the list of bools that denote if a group is
        done with its algorithm processing
        """
        return self._group_done_status

    def set_group_done_status(self, index:int, state:bool):
        """
        Sets and inded in the list of bools that denote if a group is
        done with its algorithm processing true or false
        """
        self._group_done_status[index]=state

    @property
    def is_done(self)->bool:
        """
        Returns the overall status of the EV processing of all groups
        """
        is_completed:bool = False
        if self._group_done_status.count(False) == 0:
            #its done
            is_completed = True
        return is_completed

class EVShuttle():
    """
    This is the controller so to speak for the EVTokens to talk back and forth
    bettween the EV threads to pass results and status
    """
    def __init__(self, structs_list: Sara2StructureList, mfe:Sara2SecondaryStructure, group_index:int, token:EVToken) -> None:#pylint: disable=line-too-long
        self._kcal_group_structures_list: Sara2StructureList = structs_list
        self._sara_mfestructure:Sara2SecondaryStructure = mfe
        self._group_index:int = group_index
        self._token:EVToken = token

    @property
    def kcal_group_structures_list(self)->Sara2StructureList:
        """
        Returns the list of structures that is being analyzed 
        """
        return self._kcal_group_structures_list

    @kcal_group_structures_list.setter
    def kcal_group_structures_list(self, new_list: Sara2StructureList):
        """
        Sets the list of structures that is being analyzed 
        """
        self._kcal_group_structures_list = new_list

    @property
    def sara_mfestructure(self)->Sara2SecondaryStructure:
        """
        Return the secondary structure used as the reference structure
        """
        return self._sara_mfestructure

    @sara_mfestructure.setter
    def sara_mfestructure(self, new_strucr: Sara2SecondaryStructure):
        """
        Sets the secondary structure used as the reference structure
        """
        self._sara_mfestructure = new_strucr

    @property
    def group_index(self)->int:
        """
        Returns the group index of this shuttle
        """
        return self._group_index

    @group_index.setter
    def group_index(self, new_index: int):
        """
        Sets the group index of this shuttle
        """
        self._group_index = new_index

    @property
    def token(self)->EVToken:
        """
        Returns the token tha tis feed between the threads
        """
        return self._token

    @token.setter
    def token(self, new_token: EVToken):
        """
        Sets the token tha tis feed between the threads
        """
        self._token = new_token

class EnsembleVariation():
    """
    Ensemble Variation algorithm that gives a estimated
    stability of the RNA controlling for nucleotide numbers
    and number of structures in the ensemble analyzed
    """

    def __init__(self) -> None:
        pass

    def thread_ev(self, shuttle: EVShuttle):
        """
        Access point for using multithreading to get
        EV quicker 
        """
        token:EVToken = shuttle.token
        group_num:int = shuttle.group_index
        structs_list:Sara2StructureList = shuttle.kcal_group_structures_list
        result: EV =  self.ensemble_variation_algorithm(kcal_group_structures_list=structs_list,
                                                        ref_structure=shuttle.sara_mfestructure )
        token.group_results[group_num]= result
        token.group_dict[group_num] = result
        token.group_done_status[group_num] = True

    def ensemble_variation_algorithm(self, kcal_group_structures_list: Sara2StructureList, ref_structure:Sara2SecondaryStructure)->EV:#pylint: disable=line-too-long, too-many-locals
        """
        This is the actual ensemble variation algorithm
        """
        total_ev_subscore1:int = 0
        structure_element_count = kcal_group_structures_list.num_structures

        if structure_element_count != 0:
            #need to do each char abd then structure
            #walk through each nucleotide but first prep containers grab what is needed

            #setup constants
            nuc_count = kcal_group_structures_list.nuc_count
            structure_element_count = kcal_group_structures_list.num_structures

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

            for nuc_index in range(nuc_count):
                mfe_nuc=ref_structure.structure[nuc_index]
                num_chars = list_of_nuc_lists[nuc_index].count(mfe_nuc)
                num_diff:int = num_structs - num_chars
                list_of_nuc_scores_base[nuc_index] = num_diff
                list_of_nuc_scores_subscores[nuc_index] = list_of_nuc_scores_base[nuc_index] / structure_element_count#pylint: disable=line-too-long

            total_ev_subscore1 = sum(list_of_nuc_scores_subscores)
        else:
            total_ev_subscore1 = -1

        result: EV =  EV(ev_normalized=total_ev_subscore1, ev_threshold_norm=0, ev_structure=0)
        return result
    
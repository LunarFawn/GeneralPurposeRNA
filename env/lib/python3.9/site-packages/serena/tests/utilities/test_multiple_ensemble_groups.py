import pytest
from typing import List, Dict, NamedTuple
#import  attrs

from serena.utilities.ensemble_structures import (Sara2SecondaryStructure, 
                                        Sara2StructureList, 
                                        KcalRanges)
from serena.utilities.ensemble_groups import SingleEnsembleGroup, MultipleEnsembleGroups,EnsembleSwitchStateMFEStructs
from serena.tests.utilities.test_sara_secondary_structure_lists import test_default_new_secondary_struct_list
from serena.tests.utilities.test_sara_secondary_structure import test_empty_secondary_struct

#@attrs.define
#class ListEnsembleGroups():
#    group_1:SingleEnsembleGroup = single_ensemble_group
#    group_2:SingleEnsembleGroup = single_ensemble_group_2

"""
test Ensemble state mfe structs
"""
def test_empty_ensemble_switch_state_mfe_strucs(empty_ensemble_state_mfe_strucs:EnsembleSwitchStateMFEStructs):
    test_empty_secondary_struct(empty_ensemble_state_mfe_strucs.non_switch_mfe_struct)
    test_empty_secondary_struct(empty_ensemble_state_mfe_strucs.switched_mfe_struct)

def test_set_non_switch_mfe_ensemble_switch_state_mfe_strucs(empty_ensemble_state_mfe_strucs:EnsembleSwitchStateMFEStructs, secondary_structure_3: Sara2SecondaryStructure):
    empty_ensemble_state_mfe_strucs.non_switch_mfe_struct = secondary_structure_3
    assert empty_ensemble_state_mfe_strucs.non_switch_mfe_struct == secondary_structure_3
    empty_ensemble_state_mfe_strucs.set_non_switch_mfe(kcal=-30,
                                                       struct="(((())))")
    assert empty_ensemble_state_mfe_strucs.non_switch_mfe_struct.structure == "(((())))"
    assert empty_ensemble_state_mfe_strucs.non_switch_mfe_struct.free_energy == -30

def test_set_switch_mfe_ensemble_switch_state_mfe_strucs(empty_ensemble_state_mfe_strucs:EnsembleSwitchStateMFEStructs, secondary_structure_3: Sara2SecondaryStructure):
    empty_ensemble_state_mfe_strucs.switched_mfe_struct = secondary_structure_3
    assert empty_ensemble_state_mfe_strucs.switched_mfe_struct == secondary_structure_3
    empty_ensemble_state_mfe_strucs.set_switch_mfe(kcal=-30,
                                                       struct="(((())))")
    assert empty_ensemble_state_mfe_strucs.switched_mfe_struct.structure == "(((())))"
    assert empty_ensemble_state_mfe_strucs.switched_mfe_struct.free_energy == -30    

"""
Now multiple ensemble groups
"""

def test_empty_multiple_ensemble_groups(empty_multiple_ensemble_groups:MultipleEnsembleGroups):
    assert empty_multiple_ensemble_groups.groups == []
    assert empty_multiple_ensemble_groups.raw_groups == []
    test_empty_secondary_struct(empty_multiple_ensemble_groups.non_switch_state_structure)
    test_empty_secondary_struct(empty_multiple_ensemble_groups.switched_state_structure)
    assert empty_multiple_ensemble_groups.groups_dict == {}
    assert empty_multiple_ensemble_groups.group_values == []
    assert empty_multiple_ensemble_groups.num_groups == 0
    assert empty_multiple_ensemble_groups.group_kcal_ranges == []

def test_initialized_multiple_ensemble_groups(multiple_ensemble_groups:MultipleEnsembleGroups):
    assert multiple_ensemble_groups.non_switch_state_structure.free_energy == -50
    assert multiple_ensemble_groups.switched_state_structure.free_energy == -40
    assert multiple_ensemble_groups.non_switch_state_structure.structure == '..().)'
    assert multiple_ensemble_groups.switched_state_structure.structure == '(...))'

"""
def test_set_num_groups(empty_multiple_ensemble_groups:MultipleEnsembleGroups):
    empty_multiple_ensemble_groups.num_groups = 10
    assert empty_multiple_ensemble_groups.num_groups == 10

def test_set_groups(empty_multiple_ensemble_groups:MultipleEnsembleGroups, single_ensemble_group: SingleEnsembleGroup, single_ensemble_group_2: SingleEnsembleGroup):
    
    #need to make a groups list with a new singel ensemble group
    
    group_list:List[SingleEnsembleGroup] = [single_ensemble_group,single_ensemble_group_2]
    empty_multiple_ensemble_groups.groups = group_list
    assert empty_multiple_ensemble_groups.groups[0] == single_ensemble_group
    assert empty_multiple_ensemble_groups.groups[1] == single_ensemble_group_2

def test_set_raw_groups(empty_multiple_ensemble_groups:MultipleEnsembleGroups,secondary_structures_list_2_item:Sara2StructureList, secondary_structures_list_2_item_alt:Sara2StructureList):
    raw_group_list:List[Sara2StructureList] = [secondary_structures_list_2_item, secondary_structures_list_2_item_alt]
    empty_multiple_ensemble_groups.raw_groups = raw_group_list
    assert empty_multiple_ensemble_groups.raw_groups[0] == secondary_structures_list_2_item
    assert empty_multiple_ensemble_groups.raw_groups[1] == secondary_structures_list_2_item_alt
    assert empty_multiple_ensemble_groups.total_structures == 4

def test_set_groups_dict(empty_multiple_ensemble_groups:MultipleEnsembleGroups,secondary_structures_list_2_item:Sara2StructureList, secondary_structures_list_2_item_alt:Sara2StructureList):
    group_dict: Dict[Sara2StructureList] = {}
    group_dict[2] = secondary_structures_list_2_item
    group_dict[4] = secondary_structures_list_2_item_alt
    empty_multiple_ensemble_groups.groups_dict = group_dict
    assert empty_multiple_ensemble_groups.groups_dict[2] == secondary_structures_list_2_item
    assert empty_multiple_ensemble_groups.groups_dict[4] == secondary_structures_list_2_item_alt

def test_set_group_values(empty_multiple_ensemble_groups:MultipleEnsembleGroups):
    values:List[float] = [1.1,2.2,3.3]
    empty_multiple_ensemble_groups.group_values = values
    assert empty_multiple_ensemble_groups.group_values[0] == 1.1
    assert empty_multiple_ensemble_groups.group_values[1] == 2.2
    assert empty_multiple_ensemble_groups.group_values[2] == 3.3

def test_set_kcal_ranges(empty_multiple_ensemble_groups:MultipleEnsembleGroups,kcal_range:KcalRanges, kcal_range_2:KcalRanges):
    kcal_list:List[KcalRanges]= [kcal_range, kcal_range_2]
    empty_multiple_ensemble_groups.group_kcal_ranges = kcal_list
    assert empty_multiple_ensemble_groups.group_kcal_ranges[0] == kcal_range
    assert empty_multiple_ensemble_groups.group_kcal_ranges[1] == kcal_range_2 
"""
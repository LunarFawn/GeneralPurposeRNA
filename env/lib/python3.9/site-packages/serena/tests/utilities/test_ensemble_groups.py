#from re import S
import pytest
from typing import List, Dict, NamedTuple

from serena.utilities.ensemble_structures import (Sara2SecondaryStructure, 
                                        Sara2StructureList, 
                                        KcalRanges)
from serena.utilities.ensemble_groups import SingleEnsembleGroup, MultipleEnsembleGroups, EnsembleSwitchStateMFEStructs
from serena.tests.utilities.test_sara_secondary_structure_lists import test_default_new_secondary_struct_list
from serena.tests.utilities.test_sara_secondary_structure import test_empty_secondary_struct

def test_empty_single_ensemble_group(empty_single_ensemble_group:SingleEnsembleGroup):
    #test that the ensemble list is empty default
    test_default_new_secondary_struct_list(empty_secondary_structure_list=empty_single_ensemble_group.group)
    assert empty_single_ensemble_group.multi_state_mfe_kcal == []
    assert empty_single_ensemble_group.multi_state_mfe_struct == []
    assert empty_single_ensemble_group.kcal_end == 0
    assert empty_single_ensemble_group.kcal_span == 0
    assert empty_single_ensemble_group.kcal_start == 0
    assert empty_single_ensemble_group.switch_state_structures == None

def test_set_single_ensemble_group_properties(single_ensemble_group:SingleEnsembleGroup, secondary_structures_list_2_item:Sara2StructureList):
    assert single_ensemble_group.group.sara_stuctures[0].structure == '((.)))'
    assert single_ensemble_group.group.sara_stuctures[1].structure == '..().)'
    assert single_ensemble_group.multi_state_mfe_struct[0] == '((..))'
    assert single_ensemble_group.multi_state_mfe_struct[1] == '(...))'
    assert single_ensemble_group.multi_state_mfe_kcal[0] == -10
    assert single_ensemble_group.multi_state_mfe_kcal[1] == -20
    assert single_ensemble_group.kcal_end == 10
    assert single_ensemble_group.kcal_span == 20
    assert single_ensemble_group.kcal_start == 30
    assert single_ensemble_group.switch_state_structures.non_switch_mfe_struct.structure == '(...))'
    assert single_ensemble_group.switch_state_structures.switched_mfe_struct.structure == '..().)'

def test_fancy_single_ensemble_group_properties(empty_single_ensemble_group:SingleEnsembleGroup):
    empty_single_ensemble_group.append_multi_state_mfe_data('((..))',-10)
    empty_single_ensemble_group.append_multi_state_mfe_data('(...))', -20)
    assert empty_single_ensemble_group.multi_state_mfe_struct[0] == '((..))'
    assert empty_single_ensemble_group.multi_state_mfe_struct[1] == '(...))'
    assert empty_single_ensemble_group.multi_state_mfe_kcal[0] == -10
    assert empty_single_ensemble_group.multi_state_mfe_kcal[1] == -20
    
    empty_single_ensemble_group.update_kcals(start=30,
                                            stop=10,
                                            span=20)
    assert empty_single_ensemble_group.kcal_end == 10
    assert empty_single_ensemble_group.kcal_span == 20
    assert empty_single_ensemble_group.kcal_start == 30


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

def test_initialized_multiple_ensemble_groups(initialized_multiple_ensemble_groups:MultipleEnsembleGroups):
    assert initialized_multiple_ensemble_groups.non_switch_state_structure.free_energy == -50
    assert initialized_multiple_ensemble_groups.switched_state_structure.free_energy == -40
    assert initialized_multiple_ensemble_groups.non_switch_state_structure.structure == '..().)'
    assert initialized_multiple_ensemble_groups.switched_state_structure.structure == '(...))'

def test_add_group_multiple_ensemble_groups(single_ensemble_group:SingleEnsembleGroup, empty_multiple_ensemble_groups:MultipleEnsembleGroups):
    #group_value:float = -31.5
    empty_multiple_ensemble_groups.add_group(group=single_ensemble_group)
    assert empty_multiple_ensemble_groups.num_groups == 1
    assert empty_multiple_ensemble_groups.groups == [single_ensemble_group]
    assert empty_multiple_ensemble_groups.raw_groups == [single_ensemble_group.group]
    assert empty_multiple_ensemble_groups.groups_dict == {0:single_ensemble_group.group}
    assert empty_multiple_ensemble_groups.group_values == [single_ensemble_group.kcal_start]
    assert empty_multiple_ensemble_groups.group_kcal_ranges == [KcalRanges(start=single_ensemble_group.kcal_start, stop=single_ensemble_group.kcal_end)]


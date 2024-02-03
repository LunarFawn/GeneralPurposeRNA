
import pytest

from serena.utilities. ensemble_structures import (Sara2SecondaryStructure, 
                                        Sara2StructureList, 
                                        )


"""
Test sara2 secondary structure lists
"""


def test_default_new_secondary_struct_list(empty_secondary_structure_list:Sara2StructureList):
    assert empty_secondary_structure_list.mfe_structure == ''
    assert empty_secondary_structure_list.mfe_free_energy == 0
    assert empty_secondary_structure_list.mfe_stack_energy == 0
    assert empty_secondary_structure_list.nuc_count == 0
    assert empty_secondary_structure_list.sara_stuctures == []
    assert empty_secondary_structure_list.max_free_energy == 0
    assert empty_secondary_structure_list.min_free_energy == 0
    assert empty_secondary_structure_list.max_stack_energy == 0
    assert empty_secondary_structure_list.min_stack_energy == 0
    assert empty_secondary_structure_list.num_structures == 0 
    assert empty_secondary_structure_list.free_energy_span == 0
    assert empty_secondary_structure_list.stack_energy_span == 0
    assert empty_secondary_structure_list.weighted_structure == ''

def test_add_sara_struct_sara_list(empty_secondary_structure_list:Sara2StructureList, secondary_structure_1: Sara2SecondaryStructure, secondary_structure_2: Sara2SecondaryStructure):
    empty_secondary_structure_list.add_structure(secondary_structure_1)
    empty_secondary_structure_list.add_structure(secondary_structure_2)
    assert empty_secondary_structure_list.sara_stuctures[0] == secondary_structure_1
    assert empty_secondary_structure_list.sara_stuctures[1] == secondary_structure_2
    assert empty_secondary_structure_list.mfe_structure == '((((((.((((......((((((((...)))))))).....))))((.....(((((.((....))))))).))...)))))).'
    assert empty_secondary_structure_list.mfe_free_energy == -30
    assert empty_secondary_structure_list.mfe_stack_energy == -10
    assert empty_secondary_structure_list.nuc_count == 84
    assert empty_secondary_structure_list.max_free_energy == -30
    assert empty_secondary_structure_list.min_free_energy == -50
    assert empty_secondary_structure_list.max_stack_energy == -10
    assert empty_secondary_structure_list.min_stack_energy == -20
    assert empty_secondary_structure_list.num_structures == 2 
    assert empty_secondary_structure_list.free_energy_span == 20
    assert empty_secondary_structure_list.stack_energy_span == 10

def test_set_weighted_struct_sara_list(empty_secondary_structure_list:Sara2StructureList):
    empty_secondary_structure_list.weighted_structure = '.....'
    assert empty_secondary_structure_list.weighted_structure == '.....'

def test_secondary_structure_list_2_item(secondary_structures_list_2_item:Sara2StructureList):
    assert len(secondary_structures_list_2_item.sara_stuctures) == 2
    #test structures
    assert secondary_structures_list_2_item.sara_stuctures[0].sequence == 'GCCAUA'
    assert secondary_structures_list_2_item.sara_stuctures[0].structure == '((.)))'
    assert secondary_structures_list_2_item.sara_stuctures[0].free_energy == -30
    assert secondary_structures_list_2_item.sara_stuctures[0].stack_energy == -10
    assert secondary_structures_list_2_item.sara_stuctures[1].sequence == 'GCCAUA'
    assert secondary_structures_list_2_item.sara_stuctures[1].structure == '..().)'
    assert secondary_structures_list_2_item.sara_stuctures[1].free_energy == -50
    assert secondary_structures_list_2_item.sara_stuctures[1].stack_energy == -20
    #now test the meta data stuff
    assert secondary_structures_list_2_item.mfe_free_energy == -30
    assert secondary_structures_list_2_item.mfe_stack_energy == -10
    assert secondary_structures_list_2_item.nuc_count == 6
    assert secondary_structures_list_2_item.max_free_energy == -30
    assert secondary_structures_list_2_item.min_free_energy == -50
    assert secondary_structures_list_2_item.max_stack_energy == -10
    assert secondary_structures_list_2_item.min_stack_energy == -20
    assert secondary_structures_list_2_item.num_structures == 2 
    assert secondary_structures_list_2_item.free_energy_span == 20
    assert secondary_structures_list_2_item.stack_energy_span == 10

def test_secondary_structure_list_2_item_alt(secondary_structures_list_2_item_alt:Sara2StructureList):
    assert len(secondary_structures_list_2_item_alt.sara_stuctures) == 2
    #test structures
    assert secondary_structures_list_2_item_alt.sara_stuctures[0].sequence == 'GCCAUA'
    assert secondary_structures_list_2_item_alt.sara_stuctures[0].structure == '((..))'
    assert secondary_structures_list_2_item_alt.sara_stuctures[0].free_energy == -30
    assert secondary_structures_list_2_item_alt.sara_stuctures[0].stack_energy == -10
    assert secondary_structures_list_2_item_alt.sara_stuctures[1].sequence == 'GCCAUA'
    assert secondary_structures_list_2_item_alt.sara_stuctures[1].structure == '(...))'
    assert secondary_structures_list_2_item_alt.sara_stuctures[1].free_energy == -40
    assert secondary_structures_list_2_item_alt.sara_stuctures[1].stack_energy == -30
    #now test the meta data stuff
    assert secondary_structures_list_2_item_alt.mfe_free_energy == -30
    assert secondary_structures_list_2_item_alt.mfe_stack_energy == -10
    assert secondary_structures_list_2_item_alt.nuc_count == 6
    assert secondary_structures_list_2_item_alt.max_free_energy == -30
    assert secondary_structures_list_2_item_alt.min_free_energy == -40
    assert secondary_structures_list_2_item_alt.max_stack_energy == -10
    assert secondary_structures_list_2_item_alt.min_stack_energy == -30
    assert secondary_structures_list_2_item_alt.num_structures == 2 
    assert secondary_structures_list_2_item_alt.free_energy_span == 10
    assert secondary_structures_list_2_item_alt.stack_energy_span == 20



import pytest
from typing import List

from serena.utilities.ensemble_structures import Sara2StructureList, Sara2SecondaryStructure
from serena.utilities.weighted_structures import WeightedStructure, WeightedEnsembleResult, WeightedNucCounts, WeightedComparisonResult

def test_weighted_ensemble_result(secondary_structure_3:Sara2SecondaryStructure, secondary_structure_4:Sara2SecondaryStructure):
    struct_list:List[Sara2SecondaryStructure] = []
    struct_list.append(secondary_structure_3)
    struct_list.append(secondary_structure_4)
    result:WeightedEnsembleResult = WeightedEnsembleResult(structs=struct_list)
    assert result.structs[0] == secondary_structure_3
    assert result.structs[1] == secondary_structure_4

def test_empty_weighted_nuc_counts(empty_weighted_nuc_count:WeightedNucCounts):
    assert empty_weighted_nuc_count.num_both == -1
    assert empty_weighted_nuc_count.num_bound == -1
    assert empty_weighted_nuc_count.num_dot == -1
    assert empty_weighted_nuc_count.num_nucs == -1
    assert empty_weighted_nuc_count.num_unbound == -1

def test_initialized_nuc_counts(weighted_nuc_count:WeightedNucCounts):
    assert weighted_nuc_count.num_unbound == 1
    assert weighted_nuc_count.num_both == 2
    assert weighted_nuc_count.num_bound == 3
    assert weighted_nuc_count.num_dot == 4
    assert weighted_nuc_count.num_nucs == 5

def test_setting_weighted_nuc_counts(empty_weighted_nuc_count:WeightedNucCounts):
    empty_weighted_nuc_count.num_both = 1
    empty_weighted_nuc_count.num_bound = 2
    empty_weighted_nuc_count.num_dot = 3
    empty_weighted_nuc_count.num_nucs = 4
    empty_weighted_nuc_count.num_unbound = 5
    assert empty_weighted_nuc_count.num_both == 1
    assert empty_weighted_nuc_count.num_bound == 2
    assert empty_weighted_nuc_count.num_dot == 3
    assert empty_weighted_nuc_count.num_nucs == 4
    assert empty_weighted_nuc_count.num_unbound == 5

def test_empty_weighted_result(empty_weighted_comparison_result:WeightedComparisonResult):
    assert empty_weighted_comparison_result.comp_struct == ''
    
    #do unbound struct first
    assert empty_weighted_comparison_result.unbound_mfe_struct.sequence == ''
    assert empty_weighted_comparison_result.unbound_mfe_struct.structure == ''
    assert empty_weighted_comparison_result.unbound_mfe_struct.free_energy == 0
    assert empty_weighted_comparison_result.unbound_mfe_struct.stack_energy == 0
    assert empty_weighted_comparison_result.unbound_mfe_struct.nuc_count == 0

    #now do bound struct
    assert empty_weighted_comparison_result.bound_mfe_struct.sequence == ''
    assert empty_weighted_comparison_result.bound_mfe_struct.structure == ''
    assert empty_weighted_comparison_result.bound_mfe_struct.free_energy == 0
    assert empty_weighted_comparison_result.bound_mfe_struct.stack_energy == 0
    assert empty_weighted_comparison_result.bound_mfe_struct.nuc_count == 0

    #now weighted nuc counts
    assert empty_weighted_comparison_result.weighted_nuc_counts.num_both == -1
    assert empty_weighted_comparison_result.weighted_nuc_counts.num_bound == -1
    assert empty_weighted_comparison_result.weighted_nuc_counts.num_dot == -1
    assert empty_weighted_comparison_result.weighted_nuc_counts.num_nucs == -1
    assert empty_weighted_comparison_result.weighted_nuc_counts.num_unbound == -1

def test_make_weighted_structure(weighted_struct_class:WeightedStructure, secondary_structures_list:Sara2StructureList):
    weighted_struct:Sara2SecondaryStructure = weighted_struct_class.make_weighted_struct(structure_list=secondary_structures_list)
    assert weighted_struct.structure == '(...))'

def test_make_weighted_comparison_struct(weighted_struct_class:WeightedStructure, secondary_structure_3:Sara2SecondaryStructure, secondary_structure_4:Sara2SecondaryStructure, secondary_structure_5:Sara2SecondaryStructure):
    nuc_count: int = len(secondary_structure_3.sequence)
    result: WeightedComparisonResult = weighted_struct_class.compair_weighted_structure(unbound_mfe_struct=secondary_structure_3,
                                                            bound_mfe_struct=secondary_structure_4,
                                                            weighted_result=secondary_structure_5,
                                                            nuc_count=nuc_count)
    comp_nuc_counts:WeightedNucCounts = result.weighted_nuc_counts
    assert result.comp_struct == '|-|.|+'
    assert comp_nuc_counts.num_both == 1
    assert comp_nuc_counts.num_bound == 1
    assert comp_nuc_counts.num_dot == 1
    assert comp_nuc_counts.num_unbound == 3
    assert comp_nuc_counts.num_nucs == 6

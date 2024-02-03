import pytest

from serena.utilities.comparison_structures import (
        ComparisonStructures, 
        ComparisonResult, 
        ComparisonNucCounts, 
        ComparisonNucResults
)
from serena.utilities.ensemble_structures import Sara2SecondaryStructure

def test_empty_comparison_nuc_count(empty_comparison_nuc_count:ComparisonNucCounts):
    assert empty_comparison_nuc_count.bound_count == -1
    assert empty_comparison_nuc_count.unbound_count == -1
    assert empty_comparison_nuc_count.both_count == -1
    assert empty_comparison_nuc_count.dot_count == -1
    assert empty_comparison_nuc_count.num_nucs == -1

def test_populated_comparison_nuc_count(comparison_nuc_count:ComparisonNucCounts):
    assert comparison_nuc_count.bound_count == 1
    assert comparison_nuc_count.unbound_count == 2
    assert comparison_nuc_count.both_count == 3
    assert comparison_nuc_count.dot_count == 4
    assert comparison_nuc_count.num_nucs == 5

def test_set_populated_comparison_nuc_count(empty_comparison_nuc_count:ComparisonNucCounts):
    empty_comparison_nuc_count.bound_count = 2
    empty_comparison_nuc_count.unbound_count = 4
    empty_comparison_nuc_count.both_count = 6
    empty_comparison_nuc_count.dot_count = 8
    empty_comparison_nuc_count.num_nucs  = 10
    assert empty_comparison_nuc_count.bound_count == 2
    assert empty_comparison_nuc_count.unbound_count == 4
    assert empty_comparison_nuc_count.both_count == 6
    assert empty_comparison_nuc_count.dot_count == 8
    assert empty_comparison_nuc_count.num_nucs == 10

def test_comparison_nuc_results(comparison_nuc_result:ComparisonNucResults, comparison_nuc_count:ComparisonNucCounts, comparison_nuc_count_2:ComparisonNucCounts):
    assert comparison_nuc_result.comparison_nuc_counts[0] == comparison_nuc_count
    assert comparison_nuc_result.comparison_nuc_counts[1] == comparison_nuc_count_2

def test_comparison_result(comparison_result:ComparisonResult,secondary_structure_3:Sara2SecondaryStructure,comparison_nuc_count:ComparisonNucCounts):
    assert comparison_result.comp_counts == comparison_nuc_count
    assert comparison_result.comp_struct == secondary_structure_3

def test_make_comparison_struct(secondary_structure_3:Sara2SecondaryStructure, secondary_structure_4:Sara2SecondaryStructure, secondary_structure_5:Sara2SecondaryStructure):
    nuc_count: int = len(secondary_structure_3.sequence)
    comp:ComparisonStructures = ComparisonStructures()
    result: ComparisonResult = comp.compair_structures(unbound_struct=secondary_structure_3,
                                                            bound_struct=secondary_structure_4,
                                                            reference_struct=secondary_structure_5,
                                                            nuc_count=nuc_count)
    comp_struct:Sara2SecondaryStructure = result.comp_struct
    comp_nuc_counts:ComparisonNucCounts = result.comp_counts
    assert comp_struct.structure == '|-|.|+'
    assert comp_nuc_counts.both_count == 1
    assert comp_nuc_counts.bound_count == 1
    assert comp_nuc_counts.dot_count == 1
    assert comp_nuc_counts.unbound_count == 3
    assert comp_nuc_counts.num_nucs == 6
    


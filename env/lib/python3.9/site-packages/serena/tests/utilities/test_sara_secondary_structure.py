import pytest
from typing import List

from serena.utilities.ensemble_structures import Sara2SecondaryStructure, KcalRanges


"""
Test sara2 secondary structure
"""

def test_empty_secondary_struct(empty_secondary_structure:Sara2SecondaryStructure):
    assert empty_secondary_structure.sequence == ''
    assert empty_secondary_structure.structure == ''
    assert empty_secondary_structure.free_energy == 0
    assert empty_secondary_structure.stack_energy == 0
    assert empty_secondary_structure.nuc_count == 0

def test_set_secondary_struct(secondary_structure_1:Sara2SecondaryStructure):
    assert secondary_structure_1.sequence == 'GCCAUCGCAUGAGGAUAUGCUCCCGUUUCGGGAGCAGAAGGCAUGUCACAAGACAUGAGGAUCACCCAUGUAGAUAAGAUGGCA'
    assert secondary_structure_1.structure == '((((((.((((......((((((((...)))))))).....))))((.....(((((.((....))))))).))...)))))).'
    assert secondary_structure_1.free_energy == -30
    assert secondary_structure_1.stack_energy == -10
    assert secondary_structure_1.nuc_count == 84

def test_setting_secondary_stuct_sequence(empty_secondary_structure:Sara2SecondaryStructure):
    empty_secondary_structure.sequence = 'GCCAUCGCAUGAGGAUAUGCUCCCGUUUCGGGAGCAGAAGGCAUGUCACAAGACAUGAGGAUCACCCAUGUAGAUAAGAUGGCA'
    assert empty_secondary_structure.sequence == 'GCCAUCGCAUGAGGAUAUGCUCCCGUUUCGGGAGCAGAAGGCAUGUCACAAGACAUGAGGAUCACCCAUGUAGAUAAGAUGGCA'

def test_setting_secondary_stuct_structure(empty_secondary_structure:Sara2SecondaryStructure):
    empty_secondary_structure.structure = '((((((.((((......((((((((...)))))))).....))))((.....(((((.((....))))))).))...)))))).'
    assert empty_secondary_structure.structure == '((((((.((((......((((((((...)))))))).....))))((.....(((((.((....))))))).))...)))))).'

def test_setting_secondary_stuct_free_energy(empty_secondary_structure:Sara2SecondaryStructure):
    empty_secondary_structure.free_energy = -10
    assert empty_secondary_structure.free_energy == -10

def test_setting_secondary_stuct_stack_energy(empty_secondary_structure:Sara2SecondaryStructure):
    empty_secondary_structure.stack_energy = -20
    assert empty_secondary_structure.stack_energy == -20

def test_setting_secondary_stuct_nuc_count(empty_secondary_structure:Sara2SecondaryStructure):
    empty_secondary_structure.sequence = 'GCC'
    assert empty_secondary_structure.nuc_count == 3

def test_kcal_range(empty_kcal_range:KcalRanges):
    assert empty_kcal_range.start == 0
    assert empty_kcal_range.stop == 0

def test_set_kcal_range(empty_kcal_range:KcalRanges):
    empty_kcal_range.start = 5
    empty_kcal_range.stop = 3
    assert empty_kcal_range.start == 5
    assert empty_kcal_range.stop == 3

def test_secondary_structure_5(secondary_structure_5:Sara2SecondaryStructure):
    assert secondary_structure_5.sequence == 'GCCAUA'
    assert secondary_structure_5.structure == '(...))'
    assert secondary_structure_5.free_energy == -40
    assert secondary_structure_5.stack_energy == -30
    assert secondary_structure_5.nuc_count == 6

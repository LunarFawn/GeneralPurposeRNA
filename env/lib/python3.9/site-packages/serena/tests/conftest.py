from typing import List
import pytest

from serena.utilities.ensemble_structures import (Sara2SecondaryStructure, 
                                        Sara2StructureList, 
                                        KcalRanges)

from serena.utilities.comparison_structures import ComparisonNucCounts, ComparisonNucResults, ComparisonResult
from serena.utilities.weighted_structures import WeightedNucCounts,WeightedComparisonResult, WeightedStructure, WeightedEnsembleResult
from serena.utilities.ensemble_groups import SingleEnsembleGroup, MultipleEnsembleGroups, EnsembleSwitchStateMFEStructs
from serena.utilities.ensemble_variation import EV, EVResult, EVShuttle, EVToken, EnsembleVariation
from serena.utilities.local_minima_variation import ComparisonLMV, ComparisonLMVResponse
from serena.utilities.thread_manager import EV_ThreadProcessor
from serena.interfaces.nupack4_0_28_wsl2_interface import NUPACK4Interface, NupackSettings


"""
Secondary Structure Fixtures
"""
@pytest.fixture
def kcal_range():
    return KcalRanges(start=1, stop=3)

@pytest.fixture
def kcal_range_2():
    return KcalRanges(start=2, stop=5)

@pytest.fixture()
def empty_secondary_structure():
    """
       Returns an empty sara2 secondary structure
    """
    return Sara2SecondaryStructure()


@pytest.fixture()
def secondary_structure_1():
    """
    Returns a initialed secondary structure object
    """
    sequence:str = 'GCCAUCGCAUGAGGAUAUGCUCCCGUUUCGGGAGCAGAAGGCAUGUCACAAGACAUGAGGAUCACCCAUGUAGAUAAGAUGGCA'
    structure: str = '((((((.((((......((((((((...)))))))).....))))((.....(((((.((....))))))).))...)))))).'
    free_energy:float = -30
    stack_energgy:float = -10
    return Sara2SecondaryStructure(sequence=sequence,
                                   structure=structure,
                                   free_energy=free_energy,
                                   stack_energy=stack_energgy)

@pytest.fixture()
def secondary_structure_2():
    """
    Returns a initialed secondary structure object
    """
    sequence:str = 'GCCAUCGCAUGAGGAUAUGCUCCCGUUUCGGGAGCAGAAGGCAUGUCACAAGACAUGAGGAUCACCCAUGUAGAUAAGAUGGCG'
    structure: str = '((((((.((((......((((((((...)))))))).....))))((.....(((((.((....))))))).))...)))))))'
    free_energy:float = -50
    stack_energgy:float = -20
    return Sara2SecondaryStructure(sequence=sequence,
                                   structure=structure,
                                   free_energy=free_energy,
                                   stack_energy=stack_energgy)

@pytest.fixture()
def secondary_structure_3():
    """
    Returns a initialed secondary structure object
    """
    sequence:str = 'GCCAUA'
    structure: str = '((.)))'
    free_energy:float = -30
    stack_energgy:float = -10
    return Sara2SecondaryStructure(sequence=sequence,
                                   structure=structure,
                                   free_energy=free_energy,
                                   stack_energy=stack_energgy)

@pytest.fixture()
def secondary_structure_3_1():
    """
    Returns a initialed secondary structure object
    """
    sequence:str = 'GCCAUA'
    structure: str = '((..))'
    free_energy:float = -30
    stack_energgy:float = -10
    return Sara2SecondaryStructure(sequence=sequence,
                                   structure=structure,
                                   free_energy=free_energy,
                                   stack_energy=stack_energgy)

@pytest.fixture()
def secondary_structure_4():
    """
    Returns a initialed secondary structure object
    """
    sequence:str = 'GCCAUA'
    structure: str = '..().)'
    free_energy:float = -50
    stack_energgy:float = -20
    return Sara2SecondaryStructure(sequence=sequence,
                                   structure=structure,
                                   free_energy=free_energy,
                                   stack_energy=stack_energgy)

@pytest.fixture()
def secondary_structure_5():
    """
    Returns a initialed secondary structure object
    """
    sequence:str = 'GCCAUA'
    structure: str = '(...))'
    free_energy:float = -40
    stack_energgy:float = -30
    return Sara2SecondaryStructure(sequence=sequence,
                                   structure=structure,
                                   free_energy=free_energy,
                                   stack_energy=stack_energgy)

@pytest.fixture
def empty_secondary_structure_list():
    """
       Returns an empty sara2 secondary structure
    """
    return Sara2StructureList()

@pytest.fixture
def secondary_structures_list_2_item_2(empty_secondary_structure_list:Sara2StructureList, secondary_structure_1: Sara2SecondaryStructure, secondary_structure_2: Sara2SecondaryStructure):
    empty_secondary_structure_list.add_structure(secondary_structure_1)
    empty_secondary_structure_list.add_structure(secondary_structure_2)
    return empty_secondary_structure_list

@pytest.fixture
def secondary_structures_list_2_item(secondary_structure_3, secondary_structure_3_1, secondary_structure_4, secondary_structure_5):
    structure_list:Sara2StructureList = Sara2StructureList()
    structure_list.add_structure(secondary_structure_3)
    structure_list.add_structure(secondary_structure_4)
    return structure_list

@pytest.fixture
def secondary_structures_list_2_item_alt(secondary_structure_3, secondary_structure_3_1, secondary_structure_4, secondary_structure_5):
    structure_list:Sara2StructureList = Sara2StructureList()
    structure_list.add_structure(secondary_structure_3_1)
    structure_list.add_structure(secondary_structure_5)
    return structure_list


@pytest.fixture
def secondary_structures_list(secondary_structure_3, secondary_structure_3_1, secondary_structure_4, secondary_structure_5):
    structure_list:Sara2StructureList = Sara2StructureList()
    structure_list.add_structure(secondary_structure_3)
    structure_list.add_structure(secondary_structure_3_1)
    structure_list.add_structure(secondary_structure_4)
    structure_list.add_structure(secondary_structure_5)
    return structure_list

@pytest.fixture
def empty_kcal_range():
    """
    Returns an initialized kcal range with default values
    """
    return KcalRanges()



"""
Comparison Structure Fixtures
"""

@pytest.fixture
def empty_comparison_nuc_count():
    """
    Returns an empty comparions nuc pair dataclass
    """
    return ComparisonNucCounts()

@pytest.fixture
def comparison_nuc_count():
    """
    Returns populated comparions nuc pair dataclass
    """
    return ComparisonNucCounts(bound_count = 1,
                                unbound_count=2,
                                both_count=3,
                                dot_count=4,
                                num_nucs=5)

@pytest.fixture
def comparison_nuc_count_2():
    """
    Returns populated comparions nuc pair dataclass
    """
    return ComparisonNucCounts(bound_count = 2,
                                unbound_count=3,
                                both_count=4,
                                dot_count=5,
                                num_nucs=5) 

@pytest.fixture
def comparison_nuc_count_3():
    """
    Returns populated comparions nuc pair dataclass
    """
    return ComparisonNucCounts(bound_count = 4,
                                unbound_count=3,
                                both_count=2,
                                dot_count=1,
                                num_nucs=5)                            

@pytest.fixture
def comparison_nuc_result(comparison_nuc_count:ComparisonNucCounts, comparison_nuc_count_2:ComparisonNucCounts):
    """
    Returns a comparison nuc result
    """
    comp_list:List[ComparisonNucResults] = []
    comp_list.append(comparison_nuc_count)
    comp_list.append(comparison_nuc_count_2)
    return ComparisonNucResults(comparison_nuc_counts=comp_list)

@pytest.fixture
def comparison_nuc_result_2(comparison_nuc_count:ComparisonNucCounts, comparison_nuc_count_2:ComparisonNucCounts, comparison_nuc_count_3:ComparisonNucCounts):
    """
    Returns a comparison nuc result
    """
    comp_list:List[ComparisonNucResults] = []
    comp_list.append(comparison_nuc_count)
    comp_list.append(comparison_nuc_count_2)
    comp_list.append(comparison_nuc_count_3)
    return ComparisonNucResults(comparison_nuc_counts=comp_list)

@pytest.fixture
def comparison_result(secondary_structure_3:Sara2SecondaryStructure, comparison_nuc_count:ComparisonNucCounts):
    """
    Returns a comparison result (not the list one)
    """
    return ComparisonResult(comp_struct=secondary_structure_3,
                            comp_counts=comparison_nuc_count)


"""
Weighted Structure Fixtures
"""

@pytest.fixture
def weighted_ensemble_result(secondary_structure_4:Sara2SecondaryStructure, secondary_structure_5:Sara2SecondaryStructure):
    struct_list:List[Sara2SecondaryStructure] = []
    struct_list.append(secondary_structure_4)
    struct_list.append(secondary_structure_5)
    return WeightedEnsembleResult(structs=struct_list)


@pytest.fixture
def empty_weighted_nuc_count():
    """
    Returns a empty weighted nuc count
    """
    return WeightedNucCounts()

@pytest.fixture
def weighted_nuc_count():
    """
    Returns a weighted nuc count populated at initialization
    """
    return WeightedNucCounts(num_unbound=1,
                            num_both=2,
                            num_bound=3,
                            num_dot=4,
                            num_nucs=5)

@pytest.fixture
def empty_weighted_comparison_result():
    """
    Return an empty comparison result
    """
    return WeightedComparisonResult()

@pytest.fixture
def weighted_struct_class():
    return WeightedStructure()

"""
Fixtures for ensemble groups
"""

@pytest.fixture
def empty_single_ensemble_group():
    """
    Return a empty single ensemble group class
    """
    return SingleEnsembleGroup()

@pytest.fixture
def single_ensemble_group(secondary_structures_list_2_item:Sara2StructureList, secondary_structure_4:Sara2SecondaryStructure, secondary_structure_5:Sara2SecondaryStructure):
    """
    Return a empty single ensemble group class
    """
    ensemble_group:SingleEnsembleGroup = SingleEnsembleGroup()
    ensemble_group.group = secondary_structures_list_2_item
    
    mfe_structs_list:List[str] = ['((..))','(...))']
    ensemble_group.multi_state_mfe_struct = mfe_structs_list
    
    mfe_kcal_list:List[float] = [-10,-20]
    ensemble_group.multi_state_mfe_kcal = mfe_kcal_list

    ensemble_group.switch_state_structures = EnsembleSwitchStateMFEStructs(switched_mfe_struct=secondary_structure_4,
                                                                           non_switch_mfe_struct=secondary_structure_5)
    
    ensemble_group.kcal_end = 10
    ensemble_group.kcal_span = 20
    ensemble_group.kcal_start = 30
    return ensemble_group

@pytest.fixture
def single_ensemble_group_2(secondary_structures_list_2_item_alt:Sara2StructureList):
    """
    Return a empty single ensemble group class
    """
    ensemble_group:SingleEnsembleGroup = SingleEnsembleGroup()
    ensemble_group.group = secondary_structures_list_2_item_alt
    
    mfe_structs_list:List[str] = ['(....)','..()..']
    ensemble_group.multi_state_mfe_struct = mfe_structs_list
    
    mfe_kcal_list:List[float] = [-30,-40]
    ensemble_group.multi_state_mfe_kcal = mfe_kcal_list
    
    ensemble_group.kcal_end = 40
    ensemble_group.kcal_span = 50
    ensemble_group.kcal_start = 60
    return ensemble_group


@pytest.fixture
def empty_ensemble_state_mfe_strucs():
    return EnsembleSwitchStateMFEStructs()

@pytest.fixture
def ensemble_state_mfe_structs(empty_ensemble_state_mfe_strucs:EnsembleSwitchStateMFEStructs, secondary_structure_4:Sara2SecondaryStructure, secondary_structure_5:Sara2SecondaryStructure ):
    empty_ensemble_state_mfe_strucs.non_switch_mfe_struct = secondary_structure_4
    empty_ensemble_state_mfe_strucs.switched_mfe_struct = secondary_structure_5
    return empty_ensemble_state_mfe_strucs    


@pytest.fixture
def empty_multiple_ensemble_groups(empty_ensemble_state_mfe_strucs:EnsembleSwitchStateMFEStructs):
    """
    Return a empty multiple ensemble group class
    """
    return MultipleEnsembleGroups(switch_state_structures=empty_ensemble_state_mfe_strucs)

@pytest.fixture
def initialized_multiple_ensemble_groups(empty_ensemble_state_mfe_strucs:EnsembleSwitchStateMFEStructs, secondary_structure_4:Sara2SecondaryStructure, secondary_structure_5:Sara2SecondaryStructure):
    """
    Returns a multiple ensemble groups class with
    values provided at instantiation
    """
    empty_ensemble_state_mfe_strucs.non_switch_mfe_struct = secondary_structure_4
    empty_ensemble_state_mfe_strucs.switched_mfe_struct = secondary_structure_5
    return MultipleEnsembleGroups(switch_state_structures=empty_ensemble_state_mfe_strucs)  

@pytest.fixture
def multiple_ensemble_groups(initialized_multiple_ensemble_groups:MultipleEnsembleGroups, single_ensemble_group:SingleEnsembleGroup, single_ensemble_group_2:SingleEnsembleGroup):
    """
    Returns a multiple ensemble groups class with
    values provided at instantiation
    """
    initialized_multiple_ensemble_groups.add_group(group=single_ensemble_group)
    initialized_multiple_ensemble_groups.add_group(group=single_ensemble_group_2)
    return initialized_multiple_ensemble_groups

"""
Ensemble variation fixtures
"""

@pytest.fixture
def empty_ensemble_variation():
    return EnsembleVariation()

@pytest.fixture
def empty_ev():
    """
    Return empty ev
    """
    return EV()

@pytest.fixture
def initialized_ev():
    """
    Return an initialized ev
    """
    return EV(ev_normalized=1.1,
              ev_structure=2.2,
              ev_threshold_norm=3.3)

@pytest.fixture
def initialzed_ev_2():
    """
    Return an initialized ev
    """
    return EV(ev_normalized=4.4,
              ev_structure=5.5,
              ev_threshold_norm=6.6)

@pytest.fixture
def ev_result(initialized_ev:EV, initialzed_ev_2:EV):
    """
    Return a initialized ev result
    """
    ev_list: List[EV] = [initialized_ev, initialzed_ev_2]
    return EVResult(ev_values=ev_list)

@pytest.fixture
def empty_ev_token_3_groups():
    """
    Return empty EV token initialized with 3 groups
    """
    return EVToken(num_groups=3)

@pytest.fixture
def ev_token_3_groups(empty_ev_token_3_groups: EVToken, initialized_ev:EV, initialzed_ev_2:EV):
    empty_ev_token_3_groups.set_group_dict(0,initialzed_ev_2)
    empty_ev_token_3_groups.set_group_dict(2,initialized_ev)
    empty_ev_token_3_groups.set_group_result(index=0,
                                             value=initialzed_ev_2)
    empty_ev_token_3_groups.set_group_result(index=2,
                                             value=initialized_ev)
    empty_ev_token_3_groups.set_group_done_status(0, False)
    empty_ev_token_3_groups.set_group_done_status(1, False)
    return empty_ev_token_3_groups

"""
EV shuttle
"""

@pytest.fixture
def empty_ev_shuttle_num_3():
    return EVShuttle(structs_list=Sara2StructureList(),
                        mfe=Sara2SecondaryStructure(),
                        group_index=2,
                        token=EVToken(num_groups=3))

@pytest.fixture
def ev_shuttle_group_num_3(secondary_structures_list_2_item:Sara2StructureList, secondary_structure_5:Sara2SecondaryStructure, ev_token_3_groups:EVToken):
    return EVShuttle(structs_list=secondary_structures_list_2_item,
                      mfe=secondary_structure_5,
                      group_index=1,
                      token=ev_token_3_groups)

"""
EV threadprocessor
"""
#@pytest.fixture
#def initialized_ev_thread_processor(secondary_structures_list_2_item:Sara2StructureList, secondary_structures_list_2_item_alt:Sara2StructureList, secondary_structures_list_2_item_2:Sara2StructureList ):
#    structs_list:List[Sara2StructureList] = [secondary_structures_list_2_item, secondary_structures_list_2_item_alt, secondary_structures_list_2_item_2]
#    return EV_ThreadProcessor(stuctures=structs_list)

@pytest.fixture
def empty_ev_thread_processor():
    return EV_ThreadProcessor(stuctures=[], 
                              comp_structure=Sara2SecondaryStructure(),
                              comp_struct_list_option=[])

@pytest.fixture
def ev_thread_proc_struc_list(secondary_structures_list_2_item:Sara2StructureList, secondary_structures_list_2_item_alt:Sara2StructureList, secondary_structures_list_2_item_2:Sara2StructureList):
    return [secondary_structures_list_2_item, secondary_structures_list_2_item_alt, secondary_structures_list_2_item_2]


"""
local minima variation
"""

@pytest.fixture
def empty_comparison_lmv():
    return ComparisonLMV()

@pytest.fixture
def initiailized_comparison_lmv():
    return ComparisonLMV(lmv_comp=EV(ev_normalized=1,
                                     ev_structure=2,
                                     ev_threshold_norm=3),
                        lmv_mfe=EV(ev_normalized=4,
                                   ev_structure=5,
                                   ev_threshold_norm=6),
                        lmv_rel=EV(ev_normalized=7,
                                   ev_structure=8,
                                   ev_threshold_norm=9))

@pytest.fixture
def comparison_lmv_1():
    return ComparisonLMV(lmv_comp=EV(ev_normalized=1),
                        lmv_mfe=EV(ev_normalized=4),
                        lmv_rel=EV(ev_normalized=7))

@pytest.fixture
def comparison_lmv_2():
    return ComparisonLMV(lmv_comp=EV(ev_normalized=4),
                        lmv_mfe=EV(ev_normalized=1),
                        lmv_rel=EV(ev_normalized=7))

@pytest.fixture
def comparison_lmv_3():
    return ComparisonLMV(lmv_comp=EV(ev_normalized=3),
                        lmv_mfe=EV(ev_normalized=5),
                        lmv_rel=EV(ev_normalized=7))

@pytest.fixture
def comparison_lmv_4():
    return ComparisonLMV(lmv_comp=EV(ev_normalized=5),
                        lmv_mfe=EV(ev_normalized=3),
                        lmv_rel=EV(ev_normalized=7))




@pytest.fixture
def empty_comparison_lmv_response():
    return ComparisonLMVResponse()

@pytest.fixture
def initialized_comparison_lmv_response(initiailized_comparison_lmv:ComparisonLMV):
    return ComparisonLMVResponse(lmv_comps=[initiailized_comparison_lmv])

@pytest.fixture
def comparison_lmv_response_no_switch(comparison_lmv_2, comparison_lmv_1, comparison_lmv_3):
    return ComparisonLMVResponse(lmv_comps=[comparison_lmv_2,comparison_lmv_1,comparison_lmv_3])

@pytest.fixture
def comparison_lmv_response_yes_switch(comparison_lmv_1, comparison_lmv_2, comparison_lmv_4):
    return ComparisonLMVResponse(lmv_comps=[comparison_lmv_1,comparison_lmv_2,comparison_lmv_4])

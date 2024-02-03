import pytest
from typing import List

from serena.utilities.ensemble_variation import EVShuttle, EVResult, EV, EnsembleVariation, EVToken
from serena.utilities.ensemble_structures import Sara2SecondaryStructure, Sara2StructureList
from serena.utilities.thread_manager import EV_ThreadProcessor
from serena.tests.utilities.test_sara_secondary_structure_lists import test_secondary_structure_list_2_item
from serena.tests.utilities.test_sara_secondary_structure import test_secondary_structure_5, test_empty_secondary_struct


def test_empty_ev(empty_ev:EV):
    assert empty_ev.ev_normalized == -1
    assert empty_ev.ev_structure == -1
    assert empty_ev.ev_threshold_norm == -1

def test_initialized_ev(initialized_ev:EV):
    assert initialized_ev.ev_normalized == 1.1
    assert initialized_ev.ev_structure == 2.2
    assert initialized_ev.ev_threshold_norm == 3.3

def test_set_ev_values(empty_ev:EV):
    empty_ev.ev_normalized = 2
    empty_ev.ev_structure = 4
    empty_ev.ev_threshold_norm = 6
    assert empty_ev.ev_normalized == 2
    assert empty_ev.ev_structure == 4
    assert empty_ev.ev_threshold_norm == 6   

def test_ensemble_variation_algorithm(empty_ensemble_variation:EnsembleVariation, secondary_structures_list_2_item: Sara2StructureList, secondary_structure_5:Sara2SecondaryStructure):
    result:EV = empty_ensemble_variation.ensemble_variation_algorithm(kcal_group_structures_list=secondary_structures_list_2_item,
                                                    ref_structure=secondary_structure_5)
    assert result.ev_normalized == 3.0
    assert result.ev_structure == 0
    assert result.ev_threshold_norm == 0

#def test_ensemble_variation_algorithm(secondary_structures_list_2_item: Sara2StructureList, secondary_structure_5:Sara2SecondaryStructure):
#    ensemble_variation:EnsembleVariation = EnsembleVariation()
#    result:EV = ensemble_variation.ensemble_variation_algorithm(kcal_group_structures_list=secondary_structures_list_2_item,
#                                                    ref_structure=secondary_structure_5)
#    assert result.ev_normalized == 3.0
#    assert result.ev_structure == 0
#    assert result.ev_threshold_norm == 0

def test_ev_result(ev_result:EVResult):
    assert ev_result.ev_values[0].ev_normalized == 1.1
    assert ev_result.ev_values[0].ev_structure == 2.2
    assert ev_result.ev_values[0].ev_threshold_norm == 3.3
    assert ev_result.ev_values[1].ev_normalized == 4.4
    assert ev_result.ev_values[1].ev_structure == 5.5
    assert ev_result.ev_values[1].ev_threshold_norm == 6.6


def test_empty_3_group_ev_token_group_results(empty_ev_token_3_groups:EVToken):
    #first test group results initialization
    assert len(empty_ev_token_3_groups.group_results) == 3
    assert empty_ev_token_3_groups.group_results[0].ev_normalized == -1
    assert empty_ev_token_3_groups.group_results[0].ev_structure == -1
    assert empty_ev_token_3_groups.group_results[0].ev_threshold_norm == -1
    assert empty_ev_token_3_groups.group_results[1].ev_normalized == -1
    assert empty_ev_token_3_groups.group_results[1].ev_structure == -1
    assert empty_ev_token_3_groups.group_results[1].ev_threshold_norm == -1
    assert empty_ev_token_3_groups.group_results[2].ev_normalized == -1
    assert empty_ev_token_3_groups.group_results[2].ev_structure == -1
    assert empty_ev_token_3_groups.group_results[2].ev_threshold_norm == -1

def test_empty_3_group_ev_token_group_dict(empty_ev_token_3_groups:EVToken):
    #first test group results initialization
    assert empty_ev_token_3_groups.group_dict == {}
    #assert len(empty_ev_token_3_groups.group_dict) == 3
    #assert empty_ev_token_3_groups.group_dict[0].ev_normalized == -1
    #assert empty_ev_token_3_groups.group_dict[0].ev_structure == -1
    #assert empty_ev_token_3_groups.group_dict[0].ev_threshold_norm == -1
    #assert empty_ev_token_3_groups.group_dict[1].ev_normalized == -1
    #assert empty_ev_token_3_groups.group_dict[1].ev_structure == -1
    #assert empty_ev_token_3_groups.group_dict[1].ev_threshold_norm == -1
    #assert empty_ev_token_3_groups.group_dict[2].ev_normalized == -1
    #assert empty_ev_token_3_groups.group_dict[2].ev_structure == -1
    #assert empty_ev_token_3_groups.group_dict[2].ev_threshold_norm == -1

def test_empty_3_group_ev_token_group_values(empty_ev_token_3_groups:EVToken):
    #first test group results initialization
    assert len(empty_ev_token_3_groups.group_values) == 3
    assert empty_ev_token_3_groups.group_values[0] == ''   
    assert empty_ev_token_3_groups.group_values[1] == ''
    assert empty_ev_token_3_groups.group_values[2] == ''
    

def test_empty_3_group_ev_token_group_done_status(empty_ev_token_3_groups:EVToken):
    #first test group results initialization
    assert len(empty_ev_token_3_groups.group_done_status) == 3
    assert empty_ev_token_3_groups.group_done_status[0] == False  
    assert empty_ev_token_3_groups.group_done_status[1] == False
    assert empty_ev_token_3_groups.group_done_status[2] == False
    assert empty_ev_token_3_groups.is_done == False

def test_empty_3_group_ev_token(empty_ev_token_3_groups:EVToken):
    assert len(empty_ev_token_3_groups.group_results) == 3
    assert empty_ev_token_3_groups.group_results[0].ev_normalized == -1
    assert empty_ev_token_3_groups.group_results[0].ev_structure == -1
    assert empty_ev_token_3_groups.group_results[0].ev_threshold_norm == -1
    assert empty_ev_token_3_groups.group_results[1].ev_normalized == -1
    assert empty_ev_token_3_groups.group_results[1].ev_structure == -1
    assert empty_ev_token_3_groups.group_results[1].ev_threshold_norm == -1
    assert empty_ev_token_3_groups.group_results[2].ev_normalized == -1
    assert empty_ev_token_3_groups.group_results[2].ev_structure == -1
    assert empty_ev_token_3_groups.group_results[2].ev_threshold_norm == -1
    assert empty_ev_token_3_groups.group_dict == {}
    assert len(empty_ev_token_3_groups.group_values) == 3
    assert empty_ev_token_3_groups.group_values[0] == ''   
    assert empty_ev_token_3_groups.group_values[1] == ''
    assert empty_ev_token_3_groups.group_values[2] == ''
    assert len(empty_ev_token_3_groups.group_done_status) == 3
    assert empty_ev_token_3_groups.group_done_status[0] == False  
    assert empty_ev_token_3_groups.group_done_status[1] == False
    assert empty_ev_token_3_groups.group_done_status[2] == False
    assert empty_ev_token_3_groups.is_done == False


def test_set_ev_token_group_dict(empty_ev_token_3_groups:EVToken, initialized_ev:EV, initialzed_ev_2:EV):
    empty_ev_token_3_groups.set_group_dict(2,initialized_ev)
    assert empty_ev_token_3_groups.group_dict[2] == initialized_ev

def test_set_ev_token_group_results(empty_ev_token_3_groups:EVToken, initialized_ev:EV, initialzed_ev_2:EV):
    empty_ev_token_3_groups.set_group_result(index=0,
                                             value=initialzed_ev_2)
    empty_ev_token_3_groups.set_group_result(index=2,
                                             value=initialized_ev)
    assert len(empty_ev_token_3_groups.group_results) == 3
    assert empty_ev_token_3_groups.group_results[0] == initialzed_ev_2
    assert empty_ev_token_3_groups.group_results[1].ev_normalized == -1
    assert empty_ev_token_3_groups.group_results[1].ev_structure == -1
    assert empty_ev_token_3_groups.group_results[1].ev_threshold_norm == -1
    assert empty_ev_token_3_groups.group_results[2] == initialized_ev
    #now test the EVREsult part
    ev_results:EVResult = empty_ev_token_3_groups.ev_results
    assert ev_results.ev_values[0] == initialzed_ev_2
    assert ev_results.ev_values[1].ev_normalized == -1
    assert ev_results.ev_values[1].ev_structure == -1
    assert ev_results.ev_values[1].ev_threshold_norm == -1
    assert ev_results.ev_values[2] == initialized_ev

def test_set_ev_token_group_done_status(empty_ev_token_3_groups:EVToken):
    empty_ev_token_3_groups.set_group_done_status(0, True)
    empty_ev_token_3_groups.set_group_done_status(1, True)
    assert len(empty_ev_token_3_groups.group_done_status) == 3
    assert empty_ev_token_3_groups.group_done_status[0] == True  
    assert empty_ev_token_3_groups.group_done_status[1] == True
    assert empty_ev_token_3_groups.group_done_status[2] == False
    assert empty_ev_token_3_groups.is_done == False
    empty_ev_token_3_groups.set_group_done_status(2, True)
    assert empty_ev_token_3_groups.group_done_status[0] == True  
    assert empty_ev_token_3_groups.group_done_status[1] == True
    assert empty_ev_token_3_groups.group_done_status[2] == True
    assert empty_ev_token_3_groups.is_done == True

def test_ev_token_3_groups(ev_token_3_groups:EVToken):
    #first test the groups dict
    assert ev_token_3_groups.group_dict[2].ev_normalized == 1.1
    assert ev_token_3_groups.group_dict[2].ev_structure == 2.2
    assert ev_token_3_groups.group_dict[2].ev_threshold_norm == 3.3
    #then test the group results
    assert len(ev_token_3_groups.group_results) == 3
    assert ev_token_3_groups.group_results[0].ev_normalized == 4.4
    assert ev_token_3_groups.group_results[0].ev_structure == 5.5
    assert ev_token_3_groups.group_results[0].ev_threshold_norm == 6.6
    assert ev_token_3_groups.group_results[1].ev_normalized == -1
    assert ev_token_3_groups.group_results[1].ev_structure == -1
    assert ev_token_3_groups.group_results[1].ev_threshold_norm == -1
    assert ev_token_3_groups.group_results[2].ev_normalized == 1.1
    assert ev_token_3_groups.group_results[2].ev_structure == 2.2
    assert ev_token_3_groups.group_results[2].ev_threshold_norm == 3.3
    #then test group done status
    assert len(ev_token_3_groups.group_done_status) == 3
    assert ev_token_3_groups.group_done_status[0] == False  
    assert ev_token_3_groups.group_done_status[1] == False
    assert ev_token_3_groups.group_done_status[2] == False
    assert ev_token_3_groups.is_done == False


"""
EV shuttle
"""

def test_empty_ev_shuttle_group_num_3(empty_ev_shuttle_num_3:EVShuttle):
    #test structs list first
    assert empty_ev_shuttle_num_3.kcal_group_structures_list.mfe_structure == ''
    assert empty_ev_shuttle_num_3.kcal_group_structures_list.mfe_free_energy == 0
    assert empty_ev_shuttle_num_3.kcal_group_structures_list.mfe_stack_energy == 0
    assert empty_ev_shuttle_num_3.kcal_group_structures_list.nuc_count == 0
    assert empty_ev_shuttle_num_3.kcal_group_structures_list.sara_stuctures == []
    assert empty_ev_shuttle_num_3.kcal_group_structures_list.max_free_energy == 0
    assert empty_ev_shuttle_num_3.kcal_group_structures_list.min_free_energy == 0
    assert empty_ev_shuttle_num_3.kcal_group_structures_list.max_stack_energy == 0
    assert empty_ev_shuttle_num_3.kcal_group_structures_list.min_stack_energy == 0
    assert empty_ev_shuttle_num_3.kcal_group_structures_list.num_structures == 0 
    assert empty_ev_shuttle_num_3.kcal_group_structures_list.free_energy_span == 0
    assert empty_ev_shuttle_num_3.kcal_group_structures_list.stack_energy_span == 0
    assert empty_ev_shuttle_num_3.kcal_group_structures_list.weighted_structure == ''
    
    #now test mfe structure
    assert empty_ev_shuttle_num_3.sara_mfestructure.sequence == ''
    assert empty_ev_shuttle_num_3.sara_mfestructure.structure == ''
    assert empty_ev_shuttle_num_3.sara_mfestructure.free_energy == 0
    assert empty_ev_shuttle_num_3.sara_mfestructure.stack_energy == 0
    assert empty_ev_shuttle_num_3.sara_mfestructure.nuc_count == 0  

    #now group index
    assert empty_ev_shuttle_num_3.group_index == 2

    #now test the token
    assert len(empty_ev_shuttle_num_3.token.group_results) == 3
    assert empty_ev_shuttle_num_3.token.group_results[0].ev_normalized == -1
    assert empty_ev_shuttle_num_3.token.group_results[0].ev_structure == -1
    assert empty_ev_shuttle_num_3.token.group_results[0].ev_threshold_norm == -1
    assert empty_ev_shuttle_num_3.token.group_results[1].ev_normalized == -1
    assert empty_ev_shuttle_num_3.token.group_results[1].ev_structure == -1
    assert empty_ev_shuttle_num_3.token.group_results[1].ev_threshold_norm == -1
    assert empty_ev_shuttle_num_3.token.group_results[2].ev_normalized == -1
    assert empty_ev_shuttle_num_3.token.group_results[2].ev_structure == -1
    assert empty_ev_shuttle_num_3.token.group_results[2].ev_threshold_norm == -1
    assert empty_ev_shuttle_num_3.token.group_dict == {}
    assert len(empty_ev_shuttle_num_3.token.group_values) == 3
    assert empty_ev_shuttle_num_3.token.group_values[0] == ''   
    assert empty_ev_shuttle_num_3.token.group_values[1] == ''
    assert empty_ev_shuttle_num_3.token.group_values[2] == ''
    assert len(empty_ev_shuttle_num_3.token.group_done_status) == 3
    assert empty_ev_shuttle_num_3.token.group_done_status[0] == False  
    assert empty_ev_shuttle_num_3.token.group_done_status[1] == False
    assert empty_ev_shuttle_num_3.token.group_done_status[2] == False
    assert empty_ev_shuttle_num_3.token.is_done == False


def test_ev_shuttle_group_num_3(ev_shuttle_group_num_3:EVShuttle):
    test_secondary_structure_list_2_item(ev_shuttle_group_num_3.kcal_group_structures_list)
    test_secondary_structure_5(ev_shuttle_group_num_3.sara_mfestructure)
    assert ev_shuttle_group_num_3.group_index == 1
    test_ev_token_3_groups(ev_shuttle_group_num_3.token)

def test_set_ev_shuttle_kcal_group_structures(empty_ev_shuttle_num_3:EVShuttle, secondary_structures_list_2_item:Sara2StructureList):
    empty_ev_shuttle_num_3.kcal_group_structures_list = secondary_structures_list_2_item
    assert len(empty_ev_shuttle_num_3.kcal_group_structures_list.sara_stuctures) == 2
    #test structures
    assert empty_ev_shuttle_num_3.kcal_group_structures_list.sara_stuctures[0].sequence == 'GCCAUA'
    assert empty_ev_shuttle_num_3.kcal_group_structures_list.sara_stuctures[0].structure == '((.)))'
    assert empty_ev_shuttle_num_3.kcal_group_structures_list.sara_stuctures[0].free_energy == -30
    assert empty_ev_shuttle_num_3.kcal_group_structures_list.sara_stuctures[0].stack_energy == -10
    assert empty_ev_shuttle_num_3.kcal_group_structures_list.sara_stuctures[1].sequence == 'GCCAUA'
    assert empty_ev_shuttle_num_3.kcal_group_structures_list.sara_stuctures[1].structure == '..().)'
    assert empty_ev_shuttle_num_3.kcal_group_structures_list.sara_stuctures[1].free_energy == -50
    assert empty_ev_shuttle_num_3.kcal_group_structures_list.sara_stuctures[1].stack_energy == -20
    #now test the meta data stuff
    assert empty_ev_shuttle_num_3.kcal_group_structures_list.mfe_free_energy == -30
    assert empty_ev_shuttle_num_3.kcal_group_structures_list.mfe_stack_energy == -10
    assert empty_ev_shuttle_num_3.kcal_group_structures_list.nuc_count == 6
    assert empty_ev_shuttle_num_3.kcal_group_structures_list.max_free_energy == -30
    assert empty_ev_shuttle_num_3.kcal_group_structures_list.min_free_energy == -50
    assert empty_ev_shuttle_num_3.kcal_group_structures_list.max_stack_energy == -10
    assert empty_ev_shuttle_num_3.kcal_group_structures_list.min_stack_energy == -20
    assert empty_ev_shuttle_num_3.kcal_group_structures_list.num_structures == 2 
    assert empty_ev_shuttle_num_3.kcal_group_structures_list.free_energy_span == 20
    assert empty_ev_shuttle_num_3.kcal_group_structures_list.stack_energy_span == 10

def test_set_ev_shuttle_mfe_structure(empty_ev_shuttle_num_3:EVShuttle, secondary_structure_5:Sara2SecondaryStructure):
    empty_ev_shuttle_num_3.sara_mfestructure = secondary_structure_5
    assert empty_ev_shuttle_num_3.sara_mfestructure.sequence == 'GCCAUA'
    assert empty_ev_shuttle_num_3.sara_mfestructure.structure == '(...))'
    assert empty_ev_shuttle_num_3.sara_mfestructure.free_energy == -40
    assert empty_ev_shuttle_num_3.sara_mfestructure.stack_energy == -30
    assert empty_ev_shuttle_num_3.sara_mfestructure.nuc_count == 6

def test_set_ev_shuttle_group_index(empty_ev_shuttle_num_3:EVShuttle):
    empty_ev_shuttle_num_3.group_index = 3
    assert empty_ev_shuttle_num_3.group_index == 3

def test_set_ev_shuttle_ev_token(empty_ev_shuttle_num_3:EVShuttle, ev_token_3_groups:EVToken):
    empty_ev_shuttle_num_3.token = ev_token_3_groups
    #first test the groups dict
    assert empty_ev_shuttle_num_3.token.group_dict[2].ev_normalized == 1.1
    assert empty_ev_shuttle_num_3.token.group_dict[2].ev_structure == 2.2
    assert empty_ev_shuttle_num_3.token.group_dict[2].ev_threshold_norm == 3.3
    #then test the group results
    assert len(empty_ev_shuttle_num_3.token.group_results) == 3
    assert empty_ev_shuttle_num_3.token.group_results[0].ev_normalized == 4.4
    assert empty_ev_shuttle_num_3.token.group_results[0].ev_structure == 5.5
    assert empty_ev_shuttle_num_3.token.group_results[0].ev_threshold_norm == 6.6
    assert empty_ev_shuttle_num_3.token.group_results[1].ev_normalized == -1
    assert empty_ev_shuttle_num_3.token.group_results[1].ev_structure == -1
    assert empty_ev_shuttle_num_3.token.group_results[1].ev_threshold_norm == -1
    assert empty_ev_shuttle_num_3.token.group_results[2].ev_normalized == 1.1
    assert empty_ev_shuttle_num_3.token.group_results[2].ev_structure == 2.2
    assert empty_ev_shuttle_num_3.token.group_results[2].ev_threshold_norm == 3.3
    #then test group done status
    assert len(empty_ev_shuttle_num_3.token.group_done_status) == 3
    assert empty_ev_shuttle_num_3.token.group_done_status[0] == False  
    assert empty_ev_shuttle_num_3.token.group_done_status[1] == False
    assert empty_ev_shuttle_num_3.token.group_done_status[2] == False
    assert empty_ev_shuttle_num_3.token.is_done == False

def test_thread_ev(ev_shuttle_group_num_3:EVShuttle):
    ensemble_variation:EnsembleVariation = EnsembleVariation()
    ensemble_variation.thread_ev(ev_shuttle_group_num_3)
    assert ev_shuttle_group_num_3.token.group_results[2].ev_normalized == 1.1
    assert ev_shuttle_group_num_3.token.group_results[2].ev_structure == 2.2
    assert ev_shuttle_group_num_3.token.group_results[2].ev_threshold_norm == 3.3
    assert ev_shuttle_group_num_3.token.group_results[1].ev_normalized == 3.0
    assert ev_shuttle_group_num_3.token.group_results[1].ev_structure == 0
    assert ev_shuttle_group_num_3.token.group_results[1].ev_threshold_norm == 0
    assert ev_shuttle_group_num_3.token.group_results[0].ev_normalized == 4.4
    assert ev_shuttle_group_num_3.token.group_results[0].ev_structure == 5.5
    assert ev_shuttle_group_num_3.token.group_results[0].ev_threshold_norm == 6.6

    assert ev_shuttle_group_num_3.token.group_dict[0].ev_normalized == 4.4
    assert ev_shuttle_group_num_3.token.group_dict[0].ev_structure == 5.5
    assert ev_shuttle_group_num_3.token.group_dict[0].ev_threshold_norm == 6.6
    assert ev_shuttle_group_num_3.token.group_dict[1].ev_normalized == 3.0
    assert ev_shuttle_group_num_3.token.group_dict[1].ev_structure == 0
    assert ev_shuttle_group_num_3.token.group_dict[1].ev_threshold_norm == 0 
    assert ev_shuttle_group_num_3.token.group_dict[2].ev_normalized == 1.1
    assert ev_shuttle_group_num_3.token.group_dict[2].ev_structure == 2.2
    assert ev_shuttle_group_num_3.token.group_dict[2].ev_threshold_norm == 3.3 
    
    assert ev_shuttle_group_num_3.token.group_done_status == [False,True, False] 

def test_empty_ev_thread_processor(empty_ev_thread_processor:EV_ThreadProcessor):
    assert empty_ev_thread_processor.sara2_groups == []
    test_empty_secondary_struct(empty_ev_thread_processor.comparison_structure)
    assert empty_ev_thread_processor.comp_struct_list_option == []
    assert empty_ev_thread_processor.num_groups == 0
    assert len(empty_ev_thread_processor.group_token.group_values) == 0
    assert type(empty_ev_thread_processor.EV) is EnsembleVariation

def test_initialized_empty_ev_thread_processor(secondary_structures_list_2_item:Sara2StructureList, secondary_structures_list_2_item_alt:Sara2StructureList, secondary_structures_list_2_item_2:Sara2StructureList ):
    structs_list:List[Sara2StructureList] = [secondary_structures_list_2_item, secondary_structures_list_2_item_alt, secondary_structures_list_2_item_2]
    comp_struct:Sara2SecondaryStructure = Sara2SecondaryStructure()
    initialized_ev_thread_processor:EV_ThreadProcessor = EV_ThreadProcessor(stuctures=structs_list,
                                                                            comp_structure=comp_struct)
    #test groups
    assert initialized_ev_thread_processor.sara2_groups[0] == secondary_structures_list_2_item
    assert initialized_ev_thread_processor.sara2_groups[1] == secondary_structures_list_2_item_alt
    assert initialized_ev_thread_processor.sara2_groups[2] == secondary_structures_list_2_item_2
    assert initialized_ev_thread_processor.comparison_structure == comp_struct
    assert initialized_ev_thread_processor.comp_struct_list_option == []
    assert initialized_ev_thread_processor.num_groups == 3
    test_empty_3_group_ev_token(empty_ev_token_3_groups=initialized_ev_thread_processor.group_token)
    assert type(initialized_ev_thread_processor.EV) is EnsembleVariation
    #now test by assertion
    
def test_set_ev_thread_proc_group_structs_list(empty_ev_thread_processor:EV_ThreadProcessor, secondary_structures_list_2_item:Sara2StructureList):
    new_list:List[Sara2StructureList] = [secondary_structures_list_2_item]
    empty_ev_thread_processor.sara2_groups = new_list
    assert empty_ev_thread_processor.sara2_groups == new_list

def test_set_comp_struct_ev_thread_proc(empty_ev_thread_processor:EV_ThreadProcessor, secondary_structure_1:Sara2SecondaryStructure):
    empty_ev_thread_processor.comparison_structure = secondary_structure_1
    assert empty_ev_thread_processor.comparison_structure == secondary_structure_1

def test_set_compt_struct_option_ev_thread_proc(empty_ev_thread_processor:EV_ThreadProcessor, secondary_structures_list_2_item:Sara2StructureList):
    new_list:List[Sara2StructureList] = [secondary_structures_list_2_item]
    empty_ev_thread_processor.comp_struct_list_option = new_list
    assert empty_ev_thread_processor.comp_struct_list_option == new_list

def test_set_num_groups_ev_thread_proc(empty_ev_thread_processor:EV_ThreadProcessor):
    empty_ev_thread_processor.num_groups = 5
    assert empty_ev_thread_processor.num_groups == 5

def test_set_ev_token_ev_thread_proc(empty_ev_thread_processor:EV_ThreadProcessor, ev_token_3_groups:EVToken):
    empty_ev_thread_processor.group_token = ev_token_3_groups
    assert empty_ev_thread_processor.group_token == ev_token_3_groups

def test_ev_thread_processor_comp_struct(secondary_structures_list_2_item:Sara2StructureList, secondary_structures_list_2_item_alt:Sara2StructureList, secondary_structure_5:Sara2SecondaryStructure):
    struct_list:List[Sara2StructureList]= [secondary_structures_list_2_item, secondary_structures_list_2_item_alt]
    thread_proc:EV_ThreadProcessor = EV_ThreadProcessor(stuctures=struct_list,
                                                        comp_structure=secondary_structure_5)
    result_token:EVToken = thread_proc.run_EV()
    result_token.ev_results.ev_values[0]
    assert result_token.ev_results.ev_values[0].ev_normalized == 3.0
    assert result_token.ev_results.ev_values[0].ev_structure == 0
    assert result_token.ev_results.ev_values[0].ev_threshold_norm == 0
    assert result_token.ev_results.ev_values[1].ev_normalized == 0.5
    assert result_token.ev_results.ev_values[1].ev_structure == 0
    assert result_token.ev_results.ev_values[1].ev_threshold_norm == 0





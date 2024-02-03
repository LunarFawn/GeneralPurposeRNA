import pytest



from serena.interfaces.nupack4_0_28_wsl2_interface import NUPACK4Interface, NupackSettings, MaterialParameter
from serena.utilities.ensemble_structures import  Sara2SecondaryStructure, Sara2StructureList
from serena.utilities.ensemble_groups import MultipleEnsembleGroups, EnsembleSwitchStateMFEStructs
from serena.local_minima_variation import RunLocalMinimaVariation
from serena.ensemble_variation import RunEnsembleVariation, EV
from serena.utilities.ensemble_variation import EVResult
from serena.utilities.local_minima_variation import LocalMinimaVariation
from serena.switchyness_score import RunInvestigateEnsemble, InvestigateEnsembleResults
from serena.interfaces.vienna2_fmn_hack_interface import Vienna2FMNInterface

@pytest.fixture
def nupack_4_settings():
    return NupackSettings()

@pytest.fixture
def initialized_nupack_4_settings():
    return NupackSettings(material_param=MaterialParameter.rna95_nupack4,
                          temp_C=37,
                          kcal_span_from_mfe=5,
                          Kcal_unit_increments=1,
                          sequence='ACGUACAUGAC')

@pytest.fixture
def real_world_nupack_4_settings():
    return NupackSettings(material_param=MaterialParameter.rna95_nupack4,
                          temp_C=37,
                          kcal_span_from_mfe=5,
                          Kcal_unit_increments=1,
                          sequence='GCCAUCGCAUGAGGAUAUGCUCCCGUUUCGGGAGCAGAAGGCAUGUCACAAGACAUGAGGAUCACCCAUGUAGAUAAGAUGGCA')


@pytest.fixture
def nupack_switch_states():
    unbound_struct:Sara2SecondaryStructure = Sara2SecondaryStructure(sequence='ACGUACAUGAC',
                                                                     structure='(((....).))')
    bound_struct:Sara2SecondaryStructure = Sara2SecondaryStructure(sequence='ACGUACAUGAC',
                                                                     structure='((.......))')
    
    switch_state:EnsembleSwitchStateMFEStructs = EnsembleSwitchStateMFEStructs(non_switch_mfe_struct=unbound_struct,
                                                                               switched_mfe_struct=bound_struct)
    return switch_state

@pytest.fixture
def simple_nupack_multi_group_ensemble(initialized_nupack_4_settings:NupackSettings, nupack_switch_states:EnsembleSwitchStateMFEStructs):    
    nupack4: NUPACK4Interface = NUPACK4Interface()   
    structs:Sara2StructureList = nupack4.get_subopt_energy_gap(material_param=initialized_nupack_4_settings.material_param,
                                  temp_C=initialized_nupack_4_settings.temp_C,
                                  sequence_string=initialized_nupack_4_settings.sequence,
                                  energy_delta_from_MFE=initialized_nupack_4_settings.kcal_span_from_mfe,
                                  )
    ensemble:MultipleEnsembleGroups = nupack4.load_nupack_subopt_as_ensemble(span_structures=structs,
                                                                             settings=initialized_nupack_4_settings,
                                                                             switch_state=nupack_switch_states
                                                                             )
    return ensemble


def test_empty_nupack_4_settings(nupack_4_settings:NupackSettings):
    assert nupack_4_settings.material_param == MaterialParameter.NONE
    assert nupack_4_settings.temp_C == 0
    assert nupack_4_settings.kcal_span_from_mfe == 0
    assert nupack_4_settings.Kcal_unit_increments == 0
    assert nupack_4_settings.sequence == ''

def test_initialized_nupack_4_settings(initialized_nupack_4_settings:NupackSettings):
    assert initialized_nupack_4_settings.material_param == MaterialParameter.rna95_nupack4
    assert initialized_nupack_4_settings.temp_C == 37
    assert initialized_nupack_4_settings.kcal_span_from_mfe == 5
    assert initialized_nupack_4_settings.Kcal_unit_increments == 1
    assert initialized_nupack_4_settings.sequence == 'ACGUACAUGAC'
    

def test_get_subopt_groups(initialized_nupack_4_settings:NupackSettings):
    nupack4: NUPACK4Interface = NUPACK4Interface()
    structs:Sara2StructureList = nupack4.get_subopt_energy_gap(material_param=initialized_nupack_4_settings.material_param,
                                  temp_C=initialized_nupack_4_settings.temp_C,
                                  sequence_string=initialized_nupack_4_settings.sequence,
                                  energy_delta_from_MFE=initialized_nupack_4_settings.kcal_span_from_mfe,
                                  )
    assert structs.num_structures == 10
    assert structs.sara_stuctures[0].structure == '...........'
    assert structs.sara_stuctures[0].free_energy == 0.0
    assert structs.sara_stuctures[0].stack_energy == 0.0
    assert structs.sara_stuctures[0].nuc_count == 11
    assert structs.sara_stuctures[0].sequence == 'ACGUACAUGAC'    
    
    assert structs.sara_stuctures[1].structure == '...(....)..'
    assert structs.sara_stuctures[1].free_energy == 2.0991673469543457
    assert structs.sara_stuctures[1].stack_energy == 0.0
    assert structs.sara_stuctures[1].nuc_count == 11
    assert structs.sara_stuctures[1].sequence == 'ACGUACAUGAC'
     
    assert structs.sara_stuctures[2].structure == '...(.....).'
    assert structs.sara_stuctures[2].free_energy == 1.338149070739746
    assert structs.sara_stuctures[2].stack_energy == 0.0
    assert structs.sara_stuctures[2].nuc_count == 11
    assert structs.sara_stuctures[2].sequence == 'ACGUACAUGAC'
    
    assert structs.sara_stuctures[3].structure == '..(....)...'
    assert structs.sara_stuctures[3].free_energy == 2.453690528869629
    assert structs.sara_stuctures[3].stack_energy == 0.0
    assert structs.sara_stuctures[3].nuc_count == 11
    assert structs.sara_stuctures[3].sequence == 'ACGUACAUGAC'
    
    assert structs.sara_stuctures[4].structure == '..(.......)'
    assert structs.sara_stuctures[4].free_energy == 3.604722738265991
    assert structs.sara_stuctures[4].stack_energy == 0.0
    assert structs.sara_stuctures[4].nuc_count == 11
    assert structs.sara_stuctures[4].sequence == 'ACGUACAUGAC'
    
    assert structs.sara_stuctures[5].structure == '..((....).)'
    assert structs.sara_stuctures[5].free_energy == 4.704723358154297
    assert structs.sara_stuctures[5].stack_energy == 0.0
    assert structs.sara_stuctures[5].nuc_count == 11
    assert structs.sara_stuctures[5].sequence == 'ACGUACAUGAC'

    assert structs.sara_stuctures[6].structure == '..((.....))'
    assert structs.sara_stuctures[6].free_energy == -0.0952768325805664
    assert structs.sara_stuctures[6].stack_energy == 0.0
    assert structs.sara_stuctures[6].nuc_count == 11
    assert structs.sara_stuctures[6].sequence == 'ACGUACAUGAC'
    
    assert structs.sara_stuctures[7].structure == '.(......)..'
    assert structs.sara_stuctures[7].free_energy == 1.1241339445114136
    assert structs.sara_stuctures[7].stack_energy == 0.0
    assert structs.sara_stuctures[7].nuc_count == 11
    assert structs.sara_stuctures[7].sequence == 'ACGUACAUGAC'    
    
    assert structs.sara_stuctures[8].structure == '.((....))..'
    assert structs.sara_stuctures[8].free_energy == 1.2241339683532715
    assert structs.sara_stuctures[8].stack_energy == 0.0
    assert structs.sara_stuctures[8].nuc_count == 11
    assert structs.sara_stuctures[8].sequence == 'ACGUACAUGAC'
    
    assert structs.sara_stuctures[9].structure == '(......)...'
    assert structs.sara_stuctures[9].free_energy == 3.1283416748046875
    assert structs.sara_stuctures[9].stack_energy == 0.0
    assert structs.sara_stuctures[9].nuc_count == 11
    assert structs.sara_stuctures[9].sequence == 'ACGUACAUGAC'

def test_load_nupack_subopt_as_ensemble(initialized_nupack_4_settings:NupackSettings, nupack_switch_states:EnsembleSwitchStateMFEStructs):
    nupack4: NUPACK4Interface = NUPACK4Interface()   
    structs:Sara2StructureList = nupack4.get_subopt_energy_gap(material_param=initialized_nupack_4_settings.material_param,
                                  temp_C=initialized_nupack_4_settings.temp_C,
                                  sequence_string=initialized_nupack_4_settings.sequence,
                                  energy_delta_from_MFE=initialized_nupack_4_settings.kcal_span_from_mfe,
                                  )
    ensemble:MultipleEnsembleGroups = nupack4.load_nupack_subopt_as_ensemble(span_structures=structs,
                                                                                kcal_span_from_mfe=initialized_nupack_4_settings.kcal_span_from_mfe,
                                                                                Kcal_unit_increments=initialized_nupack_4_settings.Kcal_unit_increments,
                                                                                switch_state=nupack_switch_states
                                                                                )
    assert ensemble.num_groups == 5

def test_get_lmv_nupack_ensemble(initialized_nupack_4_settings:NupackSettings, nupack_switch_states:EnsembleSwitchStateMFEStructs):
    nupack4: NUPACK4Interface = NUPACK4Interface()   
    structs:Sara2StructureList = nupack4.get_subopt_energy_gap(material_param=initialized_nupack_4_settings.material_param,
                                  temp_C=initialized_nupack_4_settings.temp_C,
                                  sequence_string=initialized_nupack_4_settings.sequence,
                                  energy_delta_from_MFE=initialized_nupack_4_settings.kcal_span_from_mfe,
                                  )
    ensemble:MultipleEnsembleGroups = nupack4.load_nupack_subopt_as_ensemble(span_structures=structs,
                                                                                kcal_span_from_mfe=initialized_nupack_4_settings.kcal_span_from_mfe,
                                                                                Kcal_unit_increments=initialized_nupack_4_settings.Kcal_unit_increments,
                                                                                switch_state=nupack_switch_states
                                                                                )
    lmv:LocalMinimaVariation = LocalMinimaVariation()
    groups_results: EVResult = lmv.get_multi_group_lmv_single_ref(ensemble=ensemble,
                                                       reference_structure=nupack_switch_states.non_switch_mfe_struct)
    assert groups_results.ev_values[0].ev_normalized == 6.0
    assert groups_results.ev_values[1].ev_normalized == 5.333333333333333
    assert groups_results.ev_values[2].ev_normalized == 6.0
    assert groups_results.ev_values[3].ev_normalized == 4.0
    assert groups_results.ev_values[4].ev_normalized == 6.0
    
def test_get_real_ev_nupack_struct_list(real_world_nupack_4_settings:NupackSettings):
    nupack4: NUPACK4Interface = NUPACK4Interface()   
    structs:Sara2StructureList = nupack4.get_subopt_energy_gap(material_param=real_world_nupack_4_settings.material_param,
                                  temp_C=real_world_nupack_4_settings.temp_C,
                                  sequence_string=real_world_nupack_4_settings.sequence,
                                  energy_delta_from_MFE=real_world_nupack_4_settings.kcal_span_from_mfe,
                                  )
    run_ev:RunEnsembleVariation = RunEnsembleVariation()
    ensemble_variation:float = run_ev.ev_from_structures_list(structures_list=structs, mfe_structure=structs.sara_stuctures[0])
    assert ensemble_variation == 10.441805225653205
    
def test_get_mfe_lmv_nupack(real_world_nupack_4_settings:NupackSettings):
    material_param=real_world_nupack_4_settings.material_param
    temp_C=real_world_nupack_4_settings.temp_C
    kcal_span_from_mfe=real_world_nupack_4_settings.kcal_span_from_mfe
    Kcal_unit_increments=real_world_nupack_4_settings.Kcal_unit_increments
    sequence=real_world_nupack_4_settings.sequence
    lmv:RunLocalMinimaVariation = RunLocalMinimaVariation()
    result:EVResult = lmv.get_mfe_multi_group_lmv_nupack(sequence=sequence,
                                                         material_param=material_param,
                                                         temp_c=temp_C,
                                                         kcal_span_from_mfe=kcal_span_from_mfe,
                                                         kcal_unit_increments=Kcal_unit_increments)
    
    assert result.ev_values[0].ev_normalized == 1
    assert result.ev_values[1].ev_normalized == 4.555555555555555
    assert result.ev_values[2].ev_normalized == 6.428571428571427
    assert result.ev_values[3].ev_normalized == 8.480392156862745
    assert result.ev_values[4].ev_normalized == 10.924444444444438
    

def test_get_comp_lmv_nupack(real_world_nupack_4_settings:NupackSettings):
    material_param=real_world_nupack_4_settings.material_param
    temp_C=real_world_nupack_4_settings.temp_C
    kcal_span_from_mfe=real_world_nupack_4_settings.kcal_span_from_mfe
    Kcal_unit_increments=real_world_nupack_4_settings.Kcal_unit_increments
    sequence=real_world_nupack_4_settings.sequence
    lmv:RunLocalMinimaVariation = RunLocalMinimaVariation()
    result:EVResult = lmv.get_comp_multi_group_lmv_nupack(sequence=sequence,
                                                        material_param=material_param,
                                                        temp_c=temp_C,
                                                        kcal_span_from_mfe=kcal_span_from_mfe,
                                                        kcal_unit_increments=Kcal_unit_increments)
    assert result.ev_values[0].ev_normalized == 1
    assert result.ev_values[1].ev_normalized == 6.222222222222221
    assert result.ev_values[2].ev_normalized == 6.19047619047619
    assert result.ev_values[3].ev_normalized == 8.872549019607842
    assert result.ev_values[4].ev_normalized == 11.826666666666659

def test_get_relative_lmv_nupack(real_world_nupack_4_settings:NupackSettings):
    material_param=real_world_nupack_4_settings.material_param
    temp_C=real_world_nupack_4_settings.temp_C
    kcal_span_from_mfe=real_world_nupack_4_settings.kcal_span_from_mfe
    Kcal_unit_increments=real_world_nupack_4_settings.Kcal_unit_increments
    sequence=real_world_nupack_4_settings.sequence
    lmv:RunLocalMinimaVariation = RunLocalMinimaVariation()
    result:EVResult = lmv.get_relative_multi_group_lmv_nupack(sequence=sequence,
                                                        material_param=material_param,
                                                        temp_c=temp_C,
                                                        kcal_span_from_mfe=kcal_span_from_mfe,
                                                        kcal_unit_increments=Kcal_unit_increments)
    assert result.ev_values[0].ev_normalized == 1
    assert result.ev_values[1].ev_normalized == 5.666666666666666
    assert result.ev_values[2].ev_normalized == 6.523809523809522
    assert result.ev_values[3].ev_normalized == 12.852941176470587
    assert result.ev_values[4].ev_normalized == 19.68444444444445
    
def test_get_folded_lmv_nupack(real_world_nupack_4_settings:NupackSettings):
    material_param=real_world_nupack_4_settings.material_param
    temp_C=real_world_nupack_4_settings.temp_C
    kcal_span_from_mfe=real_world_nupack_4_settings.kcal_span_from_mfe
    Kcal_unit_increments=real_world_nupack_4_settings.Kcal_unit_increments
    sequence=real_world_nupack_4_settings.sequence
    folded_structure:Sara2SecondaryStructure = Sara2SecondaryStructure(sequence=real_world_nupack_4_settings.sequence,
                                                                        structure='(((((...........((((((((...)))))))).....)(()))((.....((((.((....))))))).))...)))))))')
    lmv:RunLocalMinimaVariation = RunLocalMinimaVariation()
    result:EVResult = lmv.get_folded_multi_group_lmv_nupack(ensemble_sequence=sequence,
                                                        material_param=material_param,
                                                        temp_c=temp_C,
                                                        kcal_span_from_mfe=kcal_span_from_mfe,
                                                        kcal_unit_increments=Kcal_unit_increments,
                                                        folded_structure=folded_structure)
    assert result.ev_values[0].ev_normalized == 32.0
    assert result.ev_values[1].ev_normalized == 33.22222222222222
    assert result.ev_values[2].ev_normalized == 30.333333333333336
    assert result.ev_values[3].ev_normalized == 30.068627450980394
    assert result.ev_values[4].ev_normalized == 28.786666666666676
    
#def test_investigate_and_score_ensemble_nupack(real_world_nupack_4_settings:NupackSettings):
#    material_param=MaterialParameter.rna06_nupack4
#    temp_C=37
#    kcal_span_from_mfe=7
#    Kcal_unit_increments=1
#    sequence='CAGUACUGAGGAUAUACAUGAGGAUCACCCAUGUGUUCCUUGUAGGAUAUGCGCAAAAGUAGAAGGCAGAAGGACAUGUGCCUA'
#    aggressive:bool = False
#    vienna2_fmn_hack: Vienna2FMNInterface = Vienna2FMNInterface()
#    folded_structure = vienna2_fmn_hack.rnafold_fmn(input_sequence=sequence,
#                                                            do_fmn=True)
#    #folded_structure:Sara2SecondaryStructure = Sara2SecondaryStructure(sequence=real_world_nupack_4_settings.sequence,
#    #                                                                    structure='(((((...........((((((((...)))))))).....)(()))((.....((((.((....))))))).))...)))))))')
#    investigation:RunInvestigateEnsemble = RunInvestigateEnsemble()
#    result: InvestigateEnsembleResults = investigation.investigate_and_score_ensemble_nupack(sequence=sequence,
#                                                                                             folded_referenec_struct=folded_structure.structure,
#                                                                                            material_param=material_param,
#                                                                                             temp_c=temp_C,
#                                                                                             kcal_span_from_mfe=kcal_span_from_mfe,
#                                                                                             kcal_unit_increments=Kcal_unit_increments,
#                                                                                             aggressive=aggressive)
#    assert result.basic_scores.total_score == 15
#    assert result.advanced_scores.total_score == -5.526
#    assert result.number_structures == 0
    
    

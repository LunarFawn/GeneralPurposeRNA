#pylint: disable=line-too-long, too-many-arguments, too-many-locals
"""
This file is the entry point to get lmv results for an ensemble
"""

from typing import List

from serena.utilities.ensemble_structures import Sara2SecondaryStructure, Sara2StructureList
from serena.utilities.ensemble_variation import EV, EVResult, EVToken
from serena.utilities.ensemble_groups import MultipleEnsembleGroups, EnsembleSwitchStateMFEStructs
from serena.utilities.thread_manager import EV_ThreadProcessor
from serena.utilities.weighted_structures import WeightedEnsembleResult, WeightedStructure
from serena.utilities.local_minima_variation import LocalMinimaVariation
from serena.interfaces.nupack4_0_28_wsl2_interface import MaterialParameter, NUPACK4Interface

class RunLocalMinimaVariation(LocalMinimaVariation):
    """
    Class that acts as the entry point into LMV algorithm
    """

    def get_relative_mutli_group_lmv(self, ensemble: MultipleEnsembleGroups)->EVResult:
        """
        Function for getting the lmv_r for a RNA sequence for provided MultipleEnsembleGroups object
        """
        ev_values:List[EV] = []
        ref_list:List[Sara2SecondaryStructure] = []
        #for group in ensemble.groups:
        #    ref_structure:Sara2SecondaryStructure = group.group.sara_stuctures[0]
        #    ev_result:EVResult = self.get_single_group_lmv(ensemble_group=group,
        #                                                reference_structure=ref_structure)
        #    ev_values.append(ev_result.ev_values[0])
        
        for group in ensemble.groups:
            ref_structure:Sara2SecondaryStructure = group.group.sara_stuctures[0]
            ref_list.append(ref_structure)
        
        result: EVResult = self.get_multi_group_lmv_list_ref(ensemble=ensemble,
                                          reference_list=ref_list)
        return result

    def get_relative_multi_group_lmv_nupack(self, sequence:str, material_param:MaterialParameter, temp_c: int, kcal_span_from_mfe:int, kcal_unit_increments: float = 1)->EVResult:
        """
        Function for getting the lmv_r for a RNA sequence folded in silico by NUPACK4
        """
        nupack4:NUPACK4Interface = NUPACK4Interface()
        structs:Sara2StructureList = nupack4.get_subopt_energy_gap(material_param=material_param,
                                                                    temp_C=temp_c,
                                                                    sequence_string=sequence,
                                                                    energy_delta_from_MFE=kcal_span_from_mfe,
                                                                    )
        switch_states:EnsembleSwitchStateMFEStructs = EnsembleSwitchStateMFEStructs(non_switch_mfe_struct=structs.sara_stuctures[0])
        ensemble:MultipleEnsembleGroups = nupack4.load_nupack_subopt_as_ensemble(span_structures=structs,
                                                                                    kcal_span_from_mfe=kcal_span_from_mfe,
                                                                                    Kcal_unit_increments=kcal_unit_increments,
                                                                                    switch_state=switch_states
                                                                                    )
        result:EVResult = self.get_relative_mutli_group_lmv(ensemble=ensemble)
        return result

    def get_comp_multi_group_lmv(self, ensemble: MultipleEnsembleGroups, weighted_structures: WeightedEnsembleResult)->EVResult:
        """
        Function for getting the lmv_c for a RNA sequence for provided MultipleEnsembleGroups object and WeightedEnembleResult
        """
        ev_values:List[EV] = []
        ref_list:List[Sara2SecondaryStructure] = []
        #for group_index in range(len(ensemble.groups)):#pylint: disable=consider-using-enumerate
        #
        #    ref_structure:Sara2SecondaryStructure = weighted_structures.structs[group_index]
        #    ev_result:EVResult = self.get_single_group_lmv(ensemble_group=ensemble.groups[group_index],
        #                                                reference_structure=ref_structure)
        #    ev_values.append(ev_result.ev_values[0])
        
        for group_index in range(len(ensemble.groups)):#pylint: disable=consider-using-enumerate

            ref_structure:Sara2SecondaryStructure = weighted_structures.structs[group_index]
            ref_list.append(ref_structure)
        
        result: EVResult = self.get_multi_group_lmv_list_ref(ensemble=ensemble,
                                                             reference_list=ref_list)
        #result: EVResult = EVResult(ev_values=ev_values)
        return result

    def get_comp_multi_group_lmv_nupack(self, sequence:str, material_param:MaterialParameter, temp_c: int, kcal_span_from_mfe:int, kcal_unit_increments: float = 1)->EVResult:
        """
        Function for getting the lmv_c for a RNA sequence folded in silico by NUPACK4
        """
        nupack4:NUPACK4Interface = NUPACK4Interface()
        nupack_structs:Sara2StructureList = nupack4.get_subopt_energy_gap(material_param=material_param,
                                                                    temp_C=temp_c,
                                                                    sequence_string=sequence,
                                                                    energy_delta_from_MFE=kcal_span_from_mfe,
                                                                    )
        switch_states:EnsembleSwitchStateMFEStructs = EnsembleSwitchStateMFEStructs(non_switch_mfe_struct=nupack_structs.sara_stuctures[0])
        ensemble:MultipleEnsembleGroups = nupack4.load_nupack_subopt_as_ensemble(span_structures=nupack_structs,
                                                                                    kcal_span_from_mfe=kcal_span_from_mfe,
                                                                                    Kcal_unit_increments=kcal_unit_increments,
                                                                                    switch_state=switch_states
                                                                                    )
        structs:List[Sara2SecondaryStructure] = []
        #make weighted structures
        for group in ensemble.groups:
            weighted:WeightedStructure = WeightedStructure()
            struct:Sara2SecondaryStructure = weighted.make_weighted_struct(structure_list=group.group)
            structs.append(struct)

        weighted_structures:WeightedEnsembleResult = WeightedEnsembleResult(structs=structs)
        result:EVResult = self.get_comp_multi_group_lmv(ensemble=ensemble,
                                                        weighted_structures=weighted_structures)
        return result

    def get_mfe_mult_group_lmv(self, ensemble: MultipleEnsembleGroups)->EVResult:
        """
        Function for getting the lmv_m for a RNA sequence for provided MultipleEnsembleGroups object
        """
        lmv_thread: EV_ThreadProcessor = EV_ThreadProcessor(stuctures=ensemble.raw_groups,
                                                            comp_structure=ensemble.non_switch_state_structure)
        result_thread_lmv:EVToken = lmv_thread.run_EV()
        lmv_results: EVResult = result_thread_lmv.ev_results
        return lmv_results

    def get_mfe_multi_group_lmv_nupack(self, sequence:str, material_param:MaterialParameter, temp_c: int, kcal_span_from_mfe:int, kcal_unit_increments: float = 1)->EVResult:
        """
        Function for getting the lmv_m for a RNA sequence folded in silico by NUPACK4
        """
        nupack4:NUPACK4Interface = NUPACK4Interface()
        structs:Sara2StructureList = nupack4.get_subopt_energy_gap(material_param=material_param,
                                                                    temp_C=temp_c,
                                                                    sequence_string=sequence,
                                                                    energy_delta_from_MFE=kcal_span_from_mfe,
                                                                    )
        switch_states:EnsembleSwitchStateMFEStructs = EnsembleSwitchStateMFEStructs(non_switch_mfe_struct=structs.sara_stuctures[0])
        ensemble:MultipleEnsembleGroups = nupack4.load_nupack_subopt_as_ensemble(span_structures=structs,
                                                                                    kcal_span_from_mfe=kcal_span_from_mfe,
                                                                                    Kcal_unit_increments=kcal_unit_increments,
                                                                                    switch_state=switch_states
                                                                                    )
        results:EVResult = self.get_mfe_mult_group_lmv(ensemble=ensemble)
        return results

    def get_folded_group_lmv(self, ensemble: MultipleEnsembleGroups, folded_structure:Sara2SecondaryStructure)->EVResult:
        """
        Function for getting the lmv_f for a RNA sequence for provided MultipleEnsembleGroups object 
        and folded structure as Sara2SecondaryStructure
        """
        lmv_thread: EV_ThreadProcessor = EV_ThreadProcessor(stuctures=ensemble.raw_groups,
                                                            comp_structure=folded_structure)
        result_thread_lmv:EVToken = lmv_thread.run_EV()
        lmv_results: EVResult = result_thread_lmv.ev_results
        return lmv_results

    def get_folded_multi_group_lmv_nupack(self, ensemble_sequence:str, material_param:MaterialParameter, temp_c: int, kcal_span_from_mfe:int, folded_structure:Sara2SecondaryStructure, kcal_unit_increments: float = 1)->EVResult:
        """
        Function for getting the lmv_f for a RNA sequence folded in silico by NUPACK4
        """
        nupack4:NUPACK4Interface = NUPACK4Interface()
        structs:Sara2StructureList = nupack4.get_subopt_energy_gap(material_param=material_param,
                                                                    temp_C=temp_c,
                                                                    sequence_string=ensemble_sequence,
                                                                    energy_delta_from_MFE=kcal_span_from_mfe,
                                                                    )
        switch_states:EnsembleSwitchStateMFEStructs = EnsembleSwitchStateMFEStructs(non_switch_mfe_struct=structs.sara_stuctures[0],
                                                                                    switched_mfe_struct=folded_structure)
        ensemble:MultipleEnsembleGroups = nupack4.load_nupack_subopt_as_ensemble(span_structures=structs,
                                                                                    kcal_span_from_mfe=kcal_span_from_mfe,
                                                                                    Kcal_unit_increments=kcal_unit_increments,
                                                                                    switch_state=switch_states
                                                                                    )
        result:EVResult = self.get_folded_group_lmv(ensemble=ensemble,
                                                    folded_structure=folded_structure)
        return result

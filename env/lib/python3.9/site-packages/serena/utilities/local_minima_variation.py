"""
File to hold the local minima variation code
"""
import attrs
from typing import List
from dataclasses import dataclass

from serena.utilities.ensemble_groups import MultipleEnsembleGroups, SingleEnsembleGroup
from serena.utilities.ensemble_structures import Sara2SecondaryStructure
from serena.utilities.ensemble_variation import EVToken, EV, EVResult
from serena.utilities.thread_manager import EV_ThreadProcessor

@dataclass
class ComparisonLMV():
    """
    Holds the LMV's for the comparison structure
    algorithm
    """
    lmv_comp:EV = EV()
    lmv_mfe:EV = EV()
    lmv_rel:EV = EV()

@attrs.define
class ComparisonLMVResponse():
    """
    Holds the LMV's responses for the comparison structure
    algorithm
    """
    lmv_comps:List[ComparisonLMV] = []

class LocalMinimaVariation():
    """
    Local Minima Variation main algorithm access point. 
    This is teh base for finding the many flavors of LMV
    """
    def __init__(self) -> None:
        pass

    def get_multi_group_lmv_single_ref(self, ensemble: MultipleEnsembleGroups, reference_structure:Sara2SecondaryStructure):#pylint: disable=line-too-long
        """
        Return the lmv for a provided MultipleEnsemble Groups with seconfary structure as reference
        """
        lmv_thread: EV_ThreadProcessor = EV_ThreadProcessor(stuctures=ensemble.raw_groups,
                                                            comp_structure=reference_structure)
        result_thread_lmv:EVToken = lmv_thread.run_EV()
        lmv_results: EVResult = result_thread_lmv.ev_results
        return lmv_results

    def get_single_group_lmv_single_ref(self, ensemble_group: SingleEnsembleGroup, reference_structure:Sara2SecondaryStructure):#pylint: disable=line-too-long
        """
        Return the lmv for a provided SingleEnsembleGroup with seconfary structure as reference
        """
        lmv_thread: EV_ThreadProcessor = EV_ThreadProcessor(stuctures=[ensemble_group.group],
                                                              comp_structure=reference_structure)
        result_thread_lmv:EVToken = lmv_thread.run_EV()
        lmv_results: EVResult = result_thread_lmv.ev_results
        return lmv_results
    
    def get_multi_group_lmv_list_ref(self, ensemble: MultipleEnsembleGroups, reference_list:List[Sara2SecondaryStructure]):#pylint: disable=line-too-long
        """
        Return the lmv for a provided MultipleEnsemble Groups with seconfary structure as reference
        """
        lmv_thread: EV_ThreadProcessor = EV_ThreadProcessor(stuctures=ensemble.raw_groups,
                                                            comp_struct_list_option=reference_list)
        result_thread_lmv:EVToken = lmv_thread.run_EV()
        lmv_results: EVResult = result_thread_lmv.ev_results
        return lmv_results

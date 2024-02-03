"""
This is the file that holds the code that serena calls
that calls all the lower level utilities
"""

from typing import List
from dataclasses import dataclass

from serena.utilities.ensemble_structures import Sara2SecondaryStructure, Sara2StructureList
from serena.utilities.ensemble_groups import MultipleEnsembleGroups
from serena.utilities.weighted_structures import WeightedEnsembleResult, WeightedStructure
from serena.utilities.ensemble_variation import EVResult
from serena.utilities.local_minima_variation import ComparisonLMV, ComparisonLMVResponse, LocalMinimaVariation
from serena.utilities.comparison_structures import ComparisonNucCounts, ComparisonResult, ComparisonNucResults, ComparisonStructures

from serena.analysis.investigator import ComparisonInvestigator, ComparisonEvalResults, LocalMinimaVariationInvestigator, LMVAssertionResult, SettingsAssertionLMV,InvestigatorResults
from serena.analysis.judge_pool import AnalysisJudgePool, JudgesResults
from serena.analysis.scoring import SerenaScoring, BasicScoreResults, AdvancedScoreResults
from serena.local_minima_variation import RunLocalMinimaVariation

@dataclass
class ReferenceStructures():
    """
    Reference structures to use throughout the analysis
    """
    mfe_structure:Sara2SecondaryStructure
    weighted_structures: WeightedEnsembleResult


class ProcessEnsemble():
    """
    Process a MulitEnsembleGroup representation of the ensemble for various features
    This information is then given to the investigator to collect the data and hand off
    to the judges and later on scoring
    """
    def __init__(self) -> None:
        pass

    def process_ensemble_for_weighted_structures(self, ensemble:MultipleEnsembleGroups) -> WeightedEnsembleResult:
        """
        Finds and returns the weighted structures for each ensemble group and outputs it as a list of 
        weighted structures in sara2secondarystructure form.
        """
        ensemble_weighted_structures: List[Sara2SecondaryStructure] = []
        for singel_group in ensemble.groups:
            structs_list: Sara2StructureList = singel_group.group
            weighted:WeightedStructure =  WeightedStructure()
            ensemble_weighted_structures.append(weighted.make_weighted_struct(structure_list=structs_list))
        
        ensemble_result:WeightedEnsembleResult = WeightedEnsembleResult(structs=ensemble_weighted_structures)
        return ensemble_result
    
    def process_ensemble_for_lmv(self, ensemble: MultipleEnsembleGroups, ref_structures:ReferenceStructures)->ComparisonLMVResponse:
        """
        Finds and returns the lmv flavors for the ensemble groups and outputs it as a ComparisonLMVResonse
        """
        #first get mfe lmv then weighted for groups
        lmv:RunLocalMinimaVariation = RunLocalMinimaVariation()
        mfe_result:EVResult = lmv.get_mfe_mult_group_lmv(ensemble=ensemble)
        #now get ref ev
        rel_result:EVResult = lmv.get_relative_mutli_group_lmv(ensemble=ensemble)

        #now get weightedEV
        weight_result:EVResult = lmv.get_comp_multi_group_lmv(ensemble=ensemble, 
                                                                weighted_structures=ref_structures.weighted_structures)
        comparisons_lmv_response: List[ComparisonLMV] = []
        for group_index in range(len(ensemble.groups)):
            lmv_data:ComparisonLMV = ComparisonLMV()
            lmv_data.lmv_comp = weight_result.ev_values[group_index]
            lmv_data.lmv_mfe = mfe_result.ev_values[group_index]
            lmv_data.lmv_rel = rel_result.ev_values[group_index]
            comparisons_lmv_response.append(lmv_data)
        
        serena_lmv_respone: ComparisonLMVResponse = ComparisonLMVResponse(lmv_comps=comparisons_lmv_response)
        return serena_lmv_respone

    # need to feed the ensemble to this and then process it
    #compaire each weighted struct against the unbound mfe and bound structs
    def process_ensemble_for_comparison_structures(self, raw_ensemble:MultipleEnsembleGroups, weighted_ensemble:WeightedEnsembleResult)->ComparisonNucResults:
        """
        Find and return the comparison structures for the ensemble groups. This compares the
        weighted structure for each group with the unbound mfe and the folded mfe to get a compariosn
        structure for each group and returns ComparisonNucResults
        """
        nuc_count:int = raw_ensemble.groups[0].group.nuc_count

        comparison_nucs_list:List[ComparisonNucCounts]= []
        for group_index in range(raw_ensemble.num_groups):
            unbound_mfe_struct:Sara2SecondaryStructure = raw_ensemble.non_switch_state_structure
            bound_mfe_struct: Sara2SecondaryStructure = raw_ensemble.switched_state_structure
            weighted_struct:Sara2SecondaryStructure = weighted_ensemble.structs[group_index]
            comp:ComparisonStructures = ComparisonStructures()
            comparison_data:ComparisonResult = comp.compair_structures(unbound_struct=unbound_mfe_struct,
                                                                       bound_struct=bound_mfe_struct,
                                                                       reference_struct=weighted_struct,
                                                                       nuc_count=nuc_count)
            
            comparison_nuc_counts: ComparisonNucCounts = comparison_data.comp_counts
            comparison_nucs_list.append(comparison_nuc_counts)
        
        result: ComparisonNucResults = ComparisonNucResults(comparison_nuc_counts=comparison_nucs_list)
        return result
    
@dataclass
class InvestigateEnsembleResults():
    """
    Container for all the scores that are returned by the 
    scoring algorithms
    """
    basic_scores:BasicScoreResults
    advanced_scores:AdvancedScoreResults
    number_structures:int
    investigator_results:InvestigatorResults
    lmv_references:ReferenceStructures

class InvestigateEnsemble():
    """
    Entry point for automated analysis of the ensemble for
    switchyness scores
    """
    def __init__(self) -> None:
        pass

    def investigate_and_score_ensemble(self, ensemble:MultipleEnsembleGroups, is_aggressive:bool = False)->InvestigateEnsembleResults:
        """
        Does what it says. Process and investigate the MultipleEnsembleGroup
        for switchyness and report the score after judging.
        """
        process_ensemble: ProcessEnsemble = ProcessEnsemble()

        #first get weighted structs
        weighted_result: WeightedEnsembleResult = process_ensemble.process_ensemble_for_weighted_structures(ensemble=ensemble)

        #then get lmv
        lmv_references:ReferenceStructures = ReferenceStructures(mfe_structure=ensemble.non_switch_state_structure,
                                                                 weighted_structures=weighted_result)
        lmv_results:ComparisonLMVResponse = process_ensemble.process_ensemble_for_lmv(ensemble=ensemble,
                                                                                      ref_structures=lmv_references)

        #now get comparison structures
        comparison_result:ComparisonNucResults = process_ensemble.process_ensemble_for_comparison_structures(raw_ensemble=ensemble,
                                                                                                             weighted_ensemble=weighted_result)
        
        #now do the investigation
        comparison_investigator:ComparisonInvestigator = ComparisonInvestigator()

        comparison_eval_result: ComparisonEvalResults = comparison_investigator.evalulate_comparison_nucs(comparison_nucs=comparison_result)

        #use default values
        lmv_eval_settings:SettingsAssertionLMV = SettingsAssertionLMV()

        lmv_investigator:LocalMinimaVariationInvestigator = LocalMinimaVariationInvestigator()
        lmv_eval_results:LMVAssertionResult =  lmv_investigator.evaluate_lmv_for_structure_presence(lmv_data=lmv_results,
                                                             setting=lmv_eval_settings)

        investigation_results: InvestigatorResults = InvestigatorResults(comparison_eval_results=comparison_eval_result,
                                                                         comp_nuc_counts=comparison_result,
                                                                         lmv_values=lmv_results,
                                                                         lmv_assertions=lmv_eval_results,
                                                                         num_groups=ensemble.num_groups,
                                                                         total_structures_ensemble=ensemble.total_structures)
        
        #now judge the investigation
        judges:AnalysisJudgePool = AnalysisJudgePool()
        judges_decisions: JudgesResults = judges.run_all_judges(investigator=investigation_results,
                                                                is_aggressive=is_aggressive)

        #now apply scoreing to the decisions
        scoring:SerenaScoring = SerenaScoring()

        basic_scores:BasicScoreResults = scoring.basic_score_groups(judge_results=judges_decisions,
                                                                       investigator=investigation_results)
        
        advanced_scores:AdvancedScoreResults = scoring.advanced_score_groups(judge_results=judges_decisions,
                                                                             investigator=investigation_results)
        
        #add a tweak
        #if advanced_scores.total_score < 0 and ensemble.total_structures > 1000:
        #    #this is a suspect weak switch so half basic score
        #    half = basic_scores.total_score / 2
        #    basic_scores.penalties = basic_scores.penalties + half
        #    basic_scores.total_score = basic_scores.total_score - half
        
        analysis_results:InvestigateEnsembleResults = InvestigateEnsembleResults(basic_scores=basic_scores,
                                                                                 advanced_scores=advanced_scores,
                                                                                 number_structures=ensemble.total_structures,
                                                                                 investigator_results=investigation_results,
                                                                                 lmv_references=lmv_references)

        return analysis_results
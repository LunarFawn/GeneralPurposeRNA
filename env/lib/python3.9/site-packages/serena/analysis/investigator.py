#pylint: disable=line-too-long, too-few-public-methods, too-many-instance-attributes, too-many-locals,bare-except,invalid-name
"""
File for class for analysis stuff
"""

from typing import List
from dataclasses import dataclass
import attrs

from serena.utilities.comparison_structures import  ComparisonNucResults
from serena.utilities.local_minima_variation import  ComparisonLMVResponse

@dataclass
class SettingsAssertionLMV():
    """
    Settings use by algortihm that determines
    what side of switch is asserted based on lmv_m
    and lmv_c comparisons
    """
    diff_limit_mfe:float = 0
    """ The limt that the delta of lmv_c minus lmnv_m must be greater than """
    diff_limit_comp:float = 1
    """ The limt that the delta of lmv_m minus lmnv_c must be greater than """

#@dataclass
@attrs.define
class LMVAssertionResult():
    """
    Results from LVM comparisons
    """
    comp_compare_to_mfe:List[str] = []
    """ the lmv_c to lmv_m comparision as a <, >, or = symbol for each energy group"""
    unbouund_pronounced:List[bool] = []
    """ A bool that indicates if that energy group has the unbound state pronounced via lmv comparisons """
    bound_pronounced: List[bool] = []
    """ A bool that indicates if that energy group has the bound state pronounced via lmv comparisons """
    is_on_off_switch:List[bool] = []
    """ List of Bools that indicates if that energy group has indications that it is a on/off switch based on lmv comparisons in enemble
        A list is sued so tath the amount of on/off can be judged
    """

#@dataclass
#class SwitchabilitySettings():
#    limit: float = 1.5

#@dataclass
#class SwitchynessResult():
#    is_switchable_group:List[bool]
#    switchable_groups_list:List[int]
#    is_powerfull_switch_group:List[bool]
#    powerfull_groups_list:List[int]

@dataclass
class RatioResults():
    """
    Ratios for the comparison structure under analysis
    """
    unbound_to_total_ratio:float = -1
    bound_ratio: float = -1
    last_unbound_ratio: float = -1
    last_bound_ratio: float = -1
    last_both_ratio: float = -1
    bound_to_both_ratio: float = -1
    bound_to_total_ratio:float = -1
    both_nuc_total:float = -1
    dot_to_total_ratio: float = -1
    unbound_to_both:float = -1

@attrs.define
class ComparisonEvalResults():
    """
    Results from comparsion structure eval for ensemble by group
    """
    #last_count_unbound:float=0
    #last_count_bound:float=0
    #last_count_both: float = 0
    ratios:List[RatioResults] = []
    BRaise_list:List[float] =[]
    BUratio_list:List[float] = []
    bound_total_list: List[float] = []
    unbound_total_list: List[float] = []
    nuc_penatly_count:int = 0
    first_BUratio:float = 0

@dataclass
class InvestigatorResults():
    """
    Hold the results from the investigation by the investigator
    this includes all the evidence so to speak that back up the
    claims
    """
    comparison_eval_results: ComparisonEvalResults
    comp_nuc_counts: ComparisonNucResults
    lmv_values: ComparisonLMVResponse
    lmv_assertions: LMVAssertionResult
    num_groups:int = 0
    total_structures_ensemble:int = 0

class ComparisonInvestigator():
    """
    The investigator for comparison algorithm to determine if
    the comparison structures indicate that a RNA sequence is 
    capabale of acting as a switch
    """
    def __init__(self) -> None:
        pass

    def evalulate_comparison_nucs(self, comparison_nucs:ComparisonNucResults)->ComparisonEvalResults:#pylint: disable=too-many-branches,too-many-statements
        """
        Evaluate the results from the comparison nucs steps and return
        the findings
        """
        BRaise_list:List[float] = []
        BUratio_list:List[float] = []
        bound_total_list: List[float] = []
        unbound_total_list: List[float] = []
        ratios:List[RatioResults] = []
        nuc_penatly_count:int = 0
        bound_hold:int = -1
        for group_index in range(len(comparison_nucs.comparison_nuc_counts)):#pylint: disable=consider-using-enumerate
            last_index:int = 0
            if group_index > 0:
                last_index = group_index -1
            unbound:float = comparison_nucs.comparison_nuc_counts[group_index].unbound_count
            last_unbound:float = comparison_nucs.comparison_nuc_counts[last_index].unbound_count

            bound:float = comparison_nucs.comparison_nuc_counts[group_index].bound_count
            last_bound:float = comparison_nucs.comparison_nuc_counts[last_index].bound_count

            both_nuc:float = comparison_nucs.comparison_nuc_counts[group_index].both_count
            last_both:float = comparison_nucs.comparison_nuc_counts[last_index].both_count

            dot_nuc:float = comparison_nucs.comparison_nuc_counts[group_index].dot_count

            nuc_count:int = comparison_nucs.comparison_nuc_counts[last_index].num_nucs

            unbound_to_total_ratio:float = 0
            bound_to_total_ratio:float = 0
            both_nuc_total:float = 0
            bound_ratio: float = 0
            last_unbound_ratio = 0
            last_bound_ratio = 0
            last_both_ratio = 0
            bound_to_both_ratio = 0
            unbound_to_both_ratio = 0
            try:
                last_unbound_ratio = last_unbound/unbound
            except:
                pass

            try:
                bound_ratio = bound/unbound
            except:
                pass

            try:

                if bound_hold != -1:
                    #do normal
                    if bound_hold < last_bound:
                        if bound_hold == 0:
                            bound_hold = 1
                        last_bound_ratio = bound/bound_hold
                    else:
                        last_bound_ratio = bound/last_bound
                else:
                    last_bound_ratio = bound/last_bound

                if bound > last_bound:
                    #its getting bigger so record that
                    bound_hold = last_bound
                else:
                    bound_hold = -1
            except:
                pass

            #added to address the ones with 0 in the first group
            if group_index > 0:
                if BRaise_list[group_index-1] == 0 and bound > 0:
                    last_bound_ratio = bound


            try:
                last_both_ratio = both_nuc/last_both
            except:
                pass
            
            try:
                unbound_to_both_ratio = unbound/both_nuc
            except:
                pass
            
            try:
                bound_to_both_ratio = bound/(both_nuc - unbound)
            except:
                pass

            unbound_to_total_ratio = unbound/nuc_count
            bound_to_total_ratio = bound/nuc_count
            both_nuc_total= both_nuc/nuc_count
            dot_nuc_total= dot_nuc/nuc_count

            bound_total_list.append(bound_to_total_ratio)
            unbound_total_list.append(unbound_to_total_ratio)

            #now round teh data to make it more managable
            last_unbound_ratio = round(last_unbound_ratio,2)
            last_bound_ratio = round(last_bound_ratio,2)
            unbound_to_total_ratio = round(unbound_to_total_ratio,2)
            bound_ratio = round(bound_ratio,2)
            #bound_stats: str = f'BURatio:{round(bound_ratio,2)},both_Raise:{round(last_both_ratio,2)} BRaise:{round(last_bound_ratio,2)}, UDrop:{round(last_unbound_ratio,2)},BothTotal:{round(both_nuc_total,2)}, BoundTotal:{round(bound_to_total_ratio,2)}, UTotal:{round(unbound_to_total_ratio,2)}, bound_both:{round(bound_to_both_ratio,2)} B:{bound}, U:{unbound}. both:{both_nuc}'


            #this is only for the fist kcal group
            if group_index == 0:
                nuc_penatly_count = bound
                first_BUratio = float(round(bound_ratio,2))

            BUratio_list.append(float(round(bound_ratio,2)))
            BRaise_list.append(float(round(bound,2)))

            ratio_results:RatioResults = RatioResults(unbound_to_total_ratio=unbound_to_total_ratio,
                                                      bound_ratio=bound_ratio,
                                                      last_unbound_ratio=last_unbound_ratio,
                                                      last_bound_ratio=last_bound_ratio,
                                                      last_both_ratio=last_both_ratio,
                                                      bound_to_both_ratio=bound_to_both_ratio,
                                                      bound_to_total_ratio=bound_to_total_ratio,
                                                      both_nuc_total=both_nuc_total,
                                                      dot_to_total_ratio=dot_nuc_total,
                                                      unbound_to_both=unbound_to_both_ratio
                                                      )
            ratios.append(ratio_results)

        comparison_eval_results: ComparisonEvalResults = ComparisonEvalResults(ratios=ratios,
                                                                               BRaise_list=BRaise_list,
                                                                               BUratio_list=BUratio_list,
                                                                               bound_total_list=bound_total_list,
                                                                               unbound_total_list=unbound_total_list,
                                                                               nuc_penatly_count=nuc_penatly_count,
                                                                               first_BUratio=first_BUratio)
        return comparison_eval_results

class LocalMinimaVariationInvestigator():
    """
    Investigators for local minima variation results and whether 
    the lmv's of the ensemble groups indicate that the rna sequence is
    capable of performing as a switch
    """
    def __init__(self) -> None:
        pass

    def evaluate_lmv_for_structure_presence(self, lmv_data:ComparisonLMVResponse, setting:SettingsAssertionLMV)->LMVAssertionResult:
        """
        Evalute the comparison structures lmv values and determine
        if the enembled groups indicate a on/off switch. return the proof for
        this determination as well for the judges to review
        """
        ev_comp_limit: float = 25#25
        ev_mfe_limit:float = 30

        diff_limit_mfe:float = setting.diff_limit_mfe
        diff_limit_comp:float = setting.diff_limit_comp
        ev_min_limit:float = 30#15

        comp_pronounced:List[bool] = []
        is_on_off_switch:List[bool] = []
        mfe_pronounced:List[bool] = []

        for group_index in range(len(lmv_data.lmv_comps)):#pylint: disable=consider-using-enumerate
            ev_comp:float = lmv_data.lmv_comps[group_index].lmv_comp.ev_normalized
            ev_mfe:float = lmv_data.lmv_comps[group_index].lmv_mfe.ev_normalized

            comp_asserted:bool = False

            mfe_asserted:bool = False

            """"""
            diff_comp:float = round(ev_mfe,2) - round(ev_comp,2)
            if round(ev_comp,1) < round(ev_mfe,1):# and diff_comp >= diff_limit_mfe:
                #new stuff
                if group_index > 0:
                    if round(ev_comp,2) < ev_comp_limit and round(ev_mfe,2) < ev_mfe_limit:
                        mfe_asserted = True  
                else:
                    if round(ev_comp,2) < 15 and round(ev_mfe,2) < 15:
                        mfe_asserted = True
                        
            elif round(ev_comp,1) == round(ev_mfe,1):# and diff_comp >= diff_limit_mfe:
                #new stuff
                if round(ev_comp,2) < ev_comp_limit and round(ev_mfe,2) < ev_mfe_limit:
                    mfe_asserted = True
                    comp_asserted = True

            diff_mfe = round(ev_comp,2) - round(ev_mfe,2)
            if round(ev_mfe,2) < round(ev_comp,2):# and (diff_mfe >= diff_limit_comp):
                #new stuff
                if round(ev_comp,2) < ev_comp_limit and round(ev_mfe,2) < ev_mfe_limit:
                    comp_asserted = True
                
            if group_index > 0  and comp_asserted is True:
                if mfe_pronounced[0] is True:
                    is_on_off_switch.append(True)
            else:
                is_on_off_switch.append(False)

            comp_pronounced.append(comp_asserted)
            mfe_pronounced.append(mfe_asserted)

        ev_comp_to_mfe_list:List[str] = self.comp_compared_mfe_lmv(lmv_data=lmv_data)

        lmv_presence_result: LMVAssertionResult = LMVAssertionResult(comp_compare_to_mfe=ev_comp_to_mfe_list,
                                                                        unbouund_pronounced=mfe_pronounced,
                                                                        bound_pronounced=comp_pronounced,
                                                                        is_on_off_switch=is_on_off_switch)

        return lmv_presence_result

    def comp_compared_mfe_lmv(self, lmv_data:ComparisonLMVResponse)->List[str]:
        """
        Determine if the lmv_c or lmv_m is asserted per group
        """
        ev_comp_to_mfe_list:List[str] = []

        for group_index in range(len(lmv_data.lmv_comps)):#pylint: disable=consider-using-enumerate
            ev_comp:float = lmv_data.lmv_comps[group_index].lmv_comp.ev_normalized
            ev_mfe:float = lmv_data.lmv_comps[group_index].lmv_mfe.ev_normalized
            if ev_comp < ev_mfe:
                ev_comp_to_mfe_list.append('<')
            elif ev_comp == ev_mfe:
                ev_comp_to_mfe_list.append('=')
            elif ev_comp > ev_mfe:
                ev_comp_to_mfe_list.append('>')

        return ev_comp_to_mfe_list

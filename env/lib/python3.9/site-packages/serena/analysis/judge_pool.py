"""
File to hold the code to judge if RNA is a switch
"""

from dataclasses import dataclass
from typing import List
import attrs

from serena.utilities.comparison_structures import ComparisonNucCounts, ComparisonResult
from serena.utilities.ensemble_variation import EV, EVResult
from serena.utilities.local_minima_variation import ComparisonLMV
from serena.analysis.investigator import (InvestigatorResults,
                                              LMVAssertionResult)

#@dataclass
#class SwitchabilitySettings():
#    """
#    Class for holding values for limits when
#    deterining switchability
#    """
#    limit: float = 1.5 

@attrs.define
class CompSwitchJudgeResult():
    is_good_switch:bool = False
    switchable_groups_list:List[int] = []
    is_good_count:int = 0
    
    is_powerful_switch:bool  =False 
    powerfull_groups_list:List[int] = []
    is_powerful_count:int = 0
    
    is_unbound_asserted:bool = False
    unbound_asserted_groups_list: List[int]=[]
    unbound_asserted_count:int = 0

@attrs.define
class LMVSwitchJudgeResult():
    is_on_off_count: int = 0
    is_on_off_switch:bool = False
    on_off_groups_list:List[int] = []

@attrs.define
class JudgesResults():
    """
    Class to hold the results from
    the judge decisions
    """
    comp_switch_judge:CompSwitchJudgeResult = CompSwitchJudgeResult()
    lmv_switch_judge:LMVSwitchJudgeResult = LMVSwitchJudgeResult()


class AnalysisJudgePool():
    """
    Class for all the different specialized judges
    """
    def __init__(self) -> None:
        pass
        #is_powerful_switch:bool = False
        #is_good_switch:bool = False
        #is_good_count:int = 0
        #is_excelent_count:int = 0
        #current_group_index:int = -1

    def run_all_judges(self, investigator:InvestigatorResults, is_aggressive:bool = False):
        """
        Main entry point to judges to run all the judges currently
        """
        comp_judge_results: CompSwitchJudgeResult = self.is_comp_switch_judge(investigator=investigator,
                                                                              is_aggressive=is_aggressive)
        lmv_judge_results: LMVSwitchJudgeResult = self.is_lmv_switch_judge(investigator=investigator)
        results:JudgesResults = JudgesResults(comp_switch_judge=comp_judge_results,
                                              lmv_switch_judge=lmv_judge_results)
        return results 
    
    def is_comp_switch_judge(self, investigator:InvestigatorResults, is_aggressive:bool = False)->CompSwitchJudgeResult:
        """
        Judge the comp nuc comparison results from the investigator
        and return a judgment on its switchyness based on comp nuc ratios
        per energy group in the ensemble
        """        
        num_groups: int = investigator.num_groups
        
        limit: float = 1.5
        is_switchable_group:List[bool] = []
        switchable_groups_list:List[int] = []
        is_powerfull_switch_group:List[bool] = []
        powerfull_groups_list:List[int] = []
        is_good_count:int = 0
        is_excelent_count:int = 0
        is_powerful_switch:bool = False
        is_good_switch:bool = False
        is_unbound_asserted:bool = False
        unbound_asserted_list:List[int] = []
        unbound_asserted_count:int = 0

        for current_group_index in range(num_groups):
            #last_index:int = 0
            #if current_group_index>0:
            #    last_index = current_group_index-1

            last_unbound_ratio:float = investigator.comparison_eval_results.ratios[current_group_index].last_unbound_ratio
            last_unbound_ratio = round(last_unbound_ratio,2)
            last_bound_ratio: float = investigator.comparison_eval_results.ratios[current_group_index].last_bound_ratio
            last_bound_ratio = round(last_bound_ratio,2)        
            
            unbound_to_total_ratio:float = investigator.comparison_eval_results.ratios[current_group_index].unbound_to_total_ratio
            unbound_to_total_ratio = round(unbound_to_total_ratio,2)   
            bound_to_total_ratio:float = investigator.comparison_eval_results.ratios[current_group_index].bound_to_total_ratio
            bound_to_total_ratio = round(bound_to_total_ratio,2) 
              
            bound_ratio:float = investigator.comparison_eval_results.ratios[current_group_index].bound_ratio
            bound_ratio = round(bound_ratio,2)
            
            bound_to_both:float = investigator.comparison_eval_results.ratios[current_group_index].bound_to_both_ratio
            bound_to_both = round(bound_to_both,2)
            both_to_total:float = investigator.comparison_eval_results.ratios[current_group_index].both_nuc_total
            both_to_total = round(both_to_total,2)

            bound: int = investigator.comp_nuc_counts.comparison_nuc_counts[current_group_index].bound_count
            unbound:int = investigator.comp_nuc_counts.comparison_nuc_counts[current_group_index].unbound_count

            lmv_data:List[ComparisonLMV] = investigator.lmv_values.lmv_comps
            ev_weight_asserted:bool = investigator.lmv_assertions.bound_pronounced[current_group_index]#   [current_group_index].comp_pronounced
            ev_weigth_under_limit:bool = False
            ev_weight_limit_upper:int = 30
            ev_weight_limit_lower:int = 2
            if lmv_data[current_group_index].lmv_comp.ev_normalized < ev_weight_limit_upper and lmv_data[current_group_index].lmv_comp.ev_normalized > ev_weight_limit_lower:
                ev_weigth_under_limit = True 

            #sweet spot for unbound to total ratio
            unbound_to_total_limit_max:float = .5
            if is_aggressive is True:
                unbound_to_total_limit_max = .3
            in_unbound_to_total_sweet:bool = False
            if unbound_to_total_ratio <=unbound_to_total_limit_max and unbound_to_total_ratio >= .15:
                in_unbound_to_total_sweet = True
            
            in_unbound_to_total_strong:bool = False
            if unbound_to_total_ratio <=.2 and unbound_to_total_ratio >= .15:
                in_unbound_to_total_strong = True
            
            last_unbound_ratio_sweet:bool = False
            #if (last_unbound_ratio >= limit and last_unbound_ratio < 3) or (last_bound_ratio >= limit and last_bound_ratio < 3):
            if (last_unbound_ratio >= limit or last_bound_ratio >= limit):# and bound_to_both > 0:
                last_unbound_ratio_sweet = True
            
            
            bound_is_sweet:bool = False
            if bound > 2:# and bound_to_both > 0:
                bound_is_sweet = True
            
            if last_unbound_ratio_sweet is True and in_unbound_to_total_sweet is True and ev_weigth_under_limit is True and bound_is_sweet is True:#bound > 2:# and bound_to_both >= .4:
                is_good_switch = True
                switchable_groups_list.append(current_group_index)
                is_good_count = is_good_count+1

            if last_bound_ratio >=2 and unbound_to_total_ratio < .5 and ev_weigth_under_limit is True and bound_is_sweet is True:
                is_powerful_switch = True
                powerfull_groups_list.append(current_group_index)
                is_excelent_count = is_excelent_count +1
                
            #this is new stuff
            if bound_to_both > 1.5 and in_unbound_to_total_sweet is True and ev_weigth_under_limit is True and bound >= 2:# and ev_weight_asserted is True:
                is_powerful_switch = True
                powerfull_groups_list.append(current_group_index)
                is_excelent_count = is_excelent_count +1
                #is_good_switch = True
                #switchable_groups_list.append(current_group_index)
                #is_good_count = is_good_count+1

            if bound_ratio >= .5 and bound_ratio <.9 and bound_to_total_ratio > 0.09 and unbound_to_total_ratio < .4 and bound_is_sweet is True:#bound > 2:
                is_good_switch = True
                switchable_groups_list.append(current_group_index)
                is_good_count = is_good_count+1
            
            if both_to_total >= .85 and bound == 0 and unbound == 0:
                is_good_switch = True
                switchable_groups_list.append(current_group_index)
                is_good_count = is_good_count+1
            
            #if both_to_total < .75 and both_to_total > .5 and bound_ratio > .1 and bound > 2:
            #    is_good_switch = True
            #    switchable_groups_list.append(current_group_index)
            #    is_good_count = is_good_count+1

            if last_unbound_ratio >= limit and last_bound_ratio >= limit and in_unbound_to_total_strong is True and ev_weigth_under_limit is True :# and ev_weight_asserted is True:
                is_powerful_switch = True
                powerfull_groups_list.append(current_group_index)
                is_excelent_count = is_excelent_count +1
           
           
           
           #same as above basically 
           # if (last_unbound_ratio >= limit or last_bound_ratio >= limit) and in_unbound_to_total_strong is True and ev_weight_asserted is True:
           #     is_powerful_switch = True
           #     powerfull_groups_list.append(current_group_index)
           #     is_excelent_count = is_excelent_count +1
            
            #this is dangerous I think as it will hit on way to strong of switches
            #maybe so nee to use not 
            if bound_ratio >=  limit and unbound_to_total_ratio <=.15 and ev_weight_asserted is True:
                is_powerful_switch = True
                powerfull_groups_list.append(current_group_index)
                is_excelent_count = is_excelent_count +1

            if last_bound_ratio >=  2 and in_unbound_to_total_strong is True:
                is_powerful_switch = True
                powerfull_groups_list.append(current_group_index)
                is_excelent_count = is_excelent_count +1

            #if last_bound_ratio > 3 and ev_weight_asserted is True and in_unbound_to_total_sweet is True:
            #    is_good_switch = True
            #    is_powerful_switch = True
            #    is_good_count = is_good_count + 1
            #    is_excelent_count = is_excelent_count + 1
            #    switchable_groups_list.append(current_group_index)
            #    powerfull_groups_list.append(current_group_index)
            
            
            #do unbound assertions now
            #unbound_to_both:float = round(investigator.comparison_eval_results.ratios[current_group_index].unbound_to_both,2)
            #if unbound_to_both > 1.5 and both_to_total < .4:
                #is_unbound_asserted = True
                #unbound_asserted_list.append(current_group_index)
                #unbound_asserted_count += 1


        results: JudgesResults = CompSwitchJudgeResult(is_powerful_count=is_excelent_count,
                                               is_good_count=is_good_count,
                                               is_good_switch=is_good_switch,
                                               is_powerful_switch=is_powerful_switch,
                                               switchable_groups_list=switchable_groups_list,
                                               powerfull_groups_list=powerfull_groups_list,
                                               unbound_asserted_count=unbound_asserted_count,
                                               is_unbound_asserted=is_unbound_asserted,
                                               unbound_asserted_groups_list=unbound_asserted_list)

        return results

    def is_lmv_switch_judge(self, investigator:InvestigatorResults)->LMVSwitchJudgeResult:
        """
        Judge the lmv comparison results from the investigator
        and return a judgment on its switchyness based on lmv assertions
        per energy group in the ensemble
        """
        lmv_data:LMVAssertionResult = investigator.lmv_assertions
        result:LMVSwitchJudgeResult = LMVSwitchJudgeResult()
        #determine if on/off switch
        on_off_switch_list:List[bool] = lmv_data.is_on_off_switch
        
        if True in on_off_switch_list:
            result.is_on_off_switch = True
        
        #now the count
        result.is_on_off_count = lmv_data.is_on_off_switch.count(True)
    
        #decide which groups are on-off switch groups that meet the criteria
        #the investigator used
        for group_index,value in enumerate(lmv_data.is_on_off_switch):
            if value == True:
                result.on_off_groups_list.append(group_index)
        return result
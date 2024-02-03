"""
File to handles teh calsses for dealing with scores
"""

from dataclasses import dataclass
from typing import List

from serena.analysis.judge_pool import JudgesResults
from serena.analysis.investigator import InvestigatorResults

@dataclass
class BasicScoreResults():
    """
    Basic scores for switchyness
    """
    total_score:float = 0
    functional_switch_score:float = 0
    powerful_switch_score:float = 0
    on_off_switch_score:float = 0
    bonuses:float = 0
    penalties:float = 0

@dataclass
class AdvancedScoreResults():
    """
    Bonuses that amplify ability to decern switchyness
    """
    lmv_bonus:float =0
    lmv_penalty:float = 0    
    comp_bonus:float = 0
    comp_penalty:float = 0
    excess_struct_penalty:float=0
    total_score:float = 0
    


class SerenaScoring():
    """
    Scoring the results from the judges decisions
    """
    def __init__(self) -> None:
        pass

    def basic_score_groups(self, judge_results:JudgesResults, investigator: InvestigatorResults)->BasicScoreResults:
        """
        Perform basic scoring functions that determine switchyness of rna sequence
        """
        #inititalization data
        found_functional_switch: List[int] = judge_results.comp_switch_judge.switchable_groups_list 
        found_powerful_switch: List[int] = judge_results.comp_switch_judge.powerfull_groups_list
        found_on_off_switch: List[int] = judge_results.lmv_switch_judge.on_off_groups_list


        bound_range_index_plus_one:List[int] = judge_results.comp_switch_judge.switchable_groups_list
        is_powerful_switch:bool = judge_results.comp_switch_judge.is_powerful_switch
        is_functional_switch:bool = judge_results.comp_switch_judge.is_good_switch
        is_off_on_switch:bool = judge_results.lmv_switch_judge.is_on_off_switch

        #SetupScores
        total_score:float = 0
        functional_switch_score:float = 0
        powerful_switch_score:float = 0
        on_off_switch_score:float = 0
        bonuses:float = 0
        penalties:float = 0

        #main scores
        if is_powerful_switch is True:
            multiplier:int = 1
            message:str = 'Potential High Fold Change'
            #result_messages = self.log_message(message, result_messages) 
            powerful_switch_score = powerful_switch_score + (len(found_powerful_switch) * multiplier)
        
        if is_functional_switch is True: 
            multiplier:int = 2
            message:str = "Potential  Functional Switch"
            #result_messages = self.log_message(message, result_messages)
            functional_switch_score = functional_switch_score + (len(found_functional_switch) * multiplier)
        
        if is_off_on_switch is True:
            multiplier:int = 1
            message:str = "Potential  off/on leaning design via LMV"
            #result_messages = self.log_message(message, result_messages)
            #on_off_switch_score= on_off_switch_score + (len(found_on_off_switch) * multiplier)
            on_off_switch_score = on_off_switch_score + 1

        #now do penalties for assertion of unbound
        if judge_results.comp_switch_judge.is_unbound_asserted is True:
            multiplier:int = 1
            penalties = penalties + (len(judge_results.comp_switch_judge.unbound_asserted_groups_list)* multiplier)

        #now bonuses
        #for value in found_functional_switch:
        #    if value >= 0 and value <= 1 and value != -1:
        #        message:str = "Confirmned good. Add bonus point for point for functional being in first two groups"
        #        #result_messages = self.log_message(message, result_messages)
        #        #functional_switch_score += 1
        #        bonuses += 1

            #if value in found_on_off_switch:
            #    message:str = "Add bonus for functional being in range of on/off prediction"
            #    #result_messages = self.log_message(message, result_messages)
            #    #functional_switch_score += 1
            #    bonuses += .5

        #for value in found_powerful_switch:
        #    if value >= 0 and value <= 1 and value != -1:
        #        message:str = "Confirmned good. Add bonus point for high performing being in first two groups"
        #        #result_messages = self.log_message(message, result_messages)
        #        #powerful_switch_score += 1
        #        bonuses += 1

            #if value in found_on_off_switch:
            #    message:str = "Add bonus for high performing being in range of on/off prediction"
            #    #result_messages = self.log_message(message, result_messages)
            #    #powerful_switch_score += 1
            #    bonuses += .5
        
        #only count a design sas good if functionial has at leat 1 point if not
        #it is probbaly flase powerfull
        #if functional_switch_score > 0:        
        total_score = powerful_switch_score + functional_switch_score + on_off_switch_score + bonuses - penalties
        #else:
        #    total_score = 0

        basic_score_results:BasicScoreResults = BasicScoreResults(total_score=total_score,
                                                                  functional_switch_score=functional_switch_score,
                                                                  powerful_switch_score=powerful_switch_score,
                                                                  on_off_switch_score=on_off_switch_score,
                                                                  bonuses=bonuses,
                                                                  penalties=penalties)

        return basic_score_results
      
    def excessive_structures_penalties(self, num_structures: int, excess_divisor:float,excess_limit:float):
        """
        Algorithm for determining the penalty for excessive number of secondary structures in the 
        whole ensemble
        """
        #excess_divisor:float = 2000#2500
        penalty:float = 0
        if num_structures > excess_limit:
            #factor:float = ((float(num_structures) - excess_limit) / excess_divisor ) * .5
            factor:float = (float(num_structures) / excess_divisor ) * .5
            message:str = f'Exsessive structs. Found:{num_structures} penalizing {factor} points'
            #result_messages = self.log_message(message, result_messages)
            sixty_range_num:float = 50000#15000
            #penalize for too many structs
            penalty += factor
            if num_structures > sixty_range_num:
                message:str = f'Significant excess structures found: found {num_structures - sixty_range_num} structures over limit of {sixty_range_num}'
                #result_messages = self.log_message(message, result_messages)
                message:str = f'Eterna_score should be ~60 for temp group and could be good design currently has high penalty for excess structures and now yet one more penalty'
                #result_messages = self.log_message(message, result_messages)
                penalty += .5
        
        return penalty
    
    def advanced_score_groups(self, judge_results:JudgesResults, investigator: InvestigatorResults):
        """
        Bonuses and penalties that affect the fine tunning of swithyness determinations
        """
        lmv_bonus:float =0
        lmv_penalty:float = 0    
        comp_bonus:float = 0
        comp_penalty:float = 0
        total_score:float = 0

        comp_less_ratio: float = investigator.lmv_assertions.comp_compare_to_mfe.count('<') / investigator.num_groups
        com_great_ratio: float = investigator.lmv_assertions.comp_compare_to_mfe.count('>')  / investigator.num_groups
        message:str = f'ev comp great:{com_great_ratio}, ev comp less:{comp_less_ratio}'
        #result_messages = self.log_message(message, result_messages)
        if com_great_ratio > comp_less_ratio and com_great_ratio >= .7:
            message:str = "EV for comparison struct is Greater MORE OFTEN than unbound mfe so add bonus"
            #result_messages = self.log_message(message, result_messages)
            #lmv_bonus += 1
        elif comp_less_ratio > com_great_ratio and comp_less_ratio >= .5:
            message:str = "EV for mfe struct is GREATER MORE OFTEN than unbound mfe so penatly"
            #result_messages = self.log_message(message, result_messages)
            #lmv_penalty += 1
            if comp_less_ratio >= .8:
                message:str = "EV for mfe is GREATER EXTRA MORE OFTEN then mfe so minus penalty point"
                #result_messages = self.log_message(message, result_messages)
            #    lmv_penalty += 1
        
        if investigator.comparison_eval_results.nuc_penatly_count > 0:
            if investigator.comparison_eval_results.BUratio_list[0] >= .6:
                new_penalty: float = 1#investigator.comparison_eval_results.nuc_penatly_count * 1
                message:str = f'Bound unbound ratio higher than 75% so it will most likely just fold into what should have been a switch so minus {new_penalty} points'
                #result_messages = self.log_message(message, result_messages)
                comp_penalty += new_penalty
            #elif BUratio_list[0] > .60 and BUratio_list[1] < .3:
            #    new_penalty: float = nuc_penatly_count * 1
            #    message:str = f'Bound unbound ratio higher than 50% and then the 2nd energy group less than 20% so it will likely be blocked from switching so minus {new_penalty} points'
            #    result_messages = self.log_message(message, result_messages)
            #    score = score - new_penalty            
            else:
                bonus: float =  1           
                message:str = f'Bound nucs found in first energy group. Design is primed to switch so add bonus of {bonus} points'
                #result_messages = self.log_message(message, result_messages)
                comp_bonus += bonus

        #penalize for being too strong of a switch and only forms in the 2nd
        #bound:float = investigator.comp_nuc_counts.comparison_nuc_counts[index].bound_count
        #unbound: float = investigator.comp_nuc_counts.comparison_nuc_counts[index].unbound_count
        #for index, ratios in enumerate(investigator.comparison_eval_results.ratios):
        #    if ratios.unbound_to_total_ratio <= .15 and :
        #        comp_penalty += 1
        #        #its probably too strong of a switch
        is_good_actually:bool = False
       
        
        #for index, ratios in enumerate(investigator.lmv_values.lmv_comps):
        #    if ratios.lmv_rel.ev_normalized >= 40:
        #        lmv_penalty +=1 
        
        if investigator.lmv_values.lmv_comps[0].lmv_mfe.ev_normalized > 13 and investigator.lmv_values.lmv_comps[0].lmv_comp.ev_normalized > 13:
            #it is a bit too unstable in the first group and will probably fall apart
            lmv_penalty +=1
        
        if investigator.lmv_values.lmv_comps[0].lmv_mfe.ev_normalized < 13 and investigator.lmv_values.lmv_comps[0].lmv_comp.ev_normalized < 13 and investigator.lmv_values.lmv_comps[0].lmv_mfe.ev_normalized > 5:
            if investigator.lmv_values.lmv_comps[1].lmv_comp.ev_normalized < 20:# and investigator.lmv_values.lmv_comps[1].lmv_comp.ev_normalized < 20:
                lmv_bonus +=1
            else:
                #it is probably a bit too unstable to act like a proper switch
                lmv_penalty +=1
        
        #if investigator.lmv_values.lmv_comps[0].lmv_mfe.ev_normalized < 5 and investigator.lmv_values.lmv_comps[0].lmv_comp.ev_normalized < 5:
        #    #not likely to pull out of the single state as the first group is too strong for the mfe
        #    lmv_penalty +=1
        index_0_mfe:float = 0
        for index, lmv_comps in enumerate(investigator.lmv_values.lmv_comps):
            if index == 0:
                index_0_mfe = lmv_comps.lmv_mfe.ev_normalized
            if index_0_mfe == 1:
                #it probably has a very low struct count as this indicates group 1 has no variation from mfe
                if lmv_comps.lmv_mfe.ev_normalized < 3 and index > 0:
                    lmv_penalty +=1
                if lmv_comps.lmv_mfe.ev_normalized < 5 and index > 3:
                    lmv_penalty +=1
                if investigator.total_structures_ensemble > 1000:
                    #too high of struct count for mfe to be so low so penalty I think
                    lmv_penalty +=1
            elif index_0_mfe > 1:
                #it probably has a very low struct count
                if lmv_comps.lmv_mfe.ev_normalized < 2 and index == 0:
                    lmv_penalty +=1
                if lmv_comps.lmv_mfe.ev_normalized < 8 and index > 0:
                    lmv_penalty +=1
                if lmv_comps.lmv_mfe.ev_normalized < 12 and index > 2:
                    lmv_penalty +=1
            #else:
            #    if lmv_comps.lmv_comp.ev_normalized < 10 and lmv_comps.lmv_mfe.ev_normalized < 10:
            #        lmv_penalty +=1
        
        for index, ratios in enumerate(investigator.comparison_eval_results.ratios):
            if ratios.bound_ratio > .3 and ratios.bound_ratio < 1:
                comp_penalty +=1
            #if ratios.bound_to_both_ratio < 0:
                #comp_penalty +=1
        #for index, ratios in enumerate(investigator.comparison_eval_results.ratios):
            if ratios.last_bound_ratio > 5.0 or ratios.last_unbound_ratio > 5.0:
                if investigator.lmv_values.lmv_comps[0].lmv_comp.ev_normalized > 15 and investigator.lmv_values.lmv_comps[index].lmv_comp.ev_normalized > 25:
                    comp_penalty +=1   
                elif investigator.lmv_values.lmv_comps[0].lmv_comp.ev_normalized < 15 and investigator.lmv_values.lmv_comps[index].lmv_comp.ev_normalized < 25:
                    is_good_actually = True
        #for index, ratios in enumerate(investigator.comparison_eval_results.ratios):
            if ratios.both_nuc_total >= .9:
                comp_penalty +=1 
                
        #not sure if I want to use... i was ify about before and it seams not fully baked in implementatiuon. 
        # need to make a ticket for this funciton
        #if is_good_switch is True and bound_to_both_ratio >= 0.08:
        #    message:str = "Low number of both and mfe nucs in relation to bound. Add bonus point"
        #    result_messages = self.log_message(message, result_messages)
        #    score= score + 1
        excess_limit:float = 35000#this is based on new data 7500
        excess_divisor:float = 2000#2500
        excess_struct_penalty:float = self.excessive_structures_penalties(num_structures=investigator.total_structures_ensemble,
                                                                    excess_limit=excess_limit,
                                                                    excess_divisor=excess_divisor)
        if excess_struct_penalty > 5 and is_good_actually is True:
            #disregard the penatly
            total_score = lmv_bonus - lmv_penalty + comp_bonus - comp_penalty - (excess_struct_penalty/2)
        else:            
            total_score = lmv_bonus - lmv_penalty + comp_bonus - comp_penalty - excess_struct_penalty
            
        advanced_score_response: AdvancedScoreResults = AdvancedScoreResults(lmv_bonus=lmv_bonus,
                                                                             lmv_penalty=lmv_penalty,
                                                                             comp_bonus=comp_bonus,
                                                                             comp_penalty=comp_penalty,
                                                                             excess_struct_penalty=excess_struct_penalty,
                                                                             total_score=total_score)
        return advanced_score_response

    
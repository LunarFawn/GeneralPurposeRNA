
from dataclasses import dataclass
from typing import List
import collections

from serena.utilities.ensemble_structures import Sara2StructureList, Sara2SecondaryStructure
from serena.utilities.ensemble_groups import MultipleEnsembleGroups, SingleEnsembleGroup

@dataclass
class WeightedEnsembleResult():
    """
    Class that holds the resuls from weighted structurs as sara2SecondaryStructures
    """
    structs: List[Sara2SecondaryStructure]

@dataclass
class WeightedNucCounts():
    """
    Class for the weighted structure counts
    """
    num_bound:float = -1
    num_unbound:float = -1
    num_both:float = -1
    num_dot:float = -1
    num_nucs:int = -1

@dataclass
class WeightedComparisonResult():
    """
    Holds the results from weighting the structures
    """
    comp_struct: str = ''
    unbound_mfe_struct:Sara2SecondaryStructure = Sara2SecondaryStructure()
    bound_mfe_struct: Sara2SecondaryStructure = Sara2SecondaryStructure()
    weighted_nuc_counts:WeightedNucCounts = WeightedNucCounts()

class WeightedStructure():

    def __init__(self) -> None:
        pass

    def make_weighted_struct(self, structure_list: Sara2StructureList)->Sara2SecondaryStructure:
        is_bond_value: int = 2
        not_bond_value: int = -1

        nuc_poistion_values: List[int] = []
        nuc_pairs_comp_list: List[List[str]] = []
        good_nucs_each_pos: List[bool] = []

        struct_count: int = structure_list.num_structures

        for nucIndex in range(structure_list.nuc_count):
            nuc_poistion_values.append(0)
            pairs_list: List[str] = []            
            nuc_pairs_comp_list.append(pairs_list)
            #good_nucs_each_pos.append(False)

        for struct in structure_list.sara_stuctures:
            for nucIndex in range(structure_list.nuc_count):
                nuc_bond_type:str = struct.structure[nucIndex]
                nuc_pairs_comp_list[nucIndex].append(nuc_bond_type)
                adder: int = 0
                if nuc_bond_type == '.':
                    adder = not_bond_value
                else:
                    adder = is_bond_value
                nuc_poistion_values[nucIndex] = nuc_poistion_values[nucIndex] + adder

        #now record if the nuc position has a weghted bond
        for nucIndex in range(structure_list.nuc_count):
            is_weighted_bond=False
            if nuc_poistion_values[nucIndex] > struct_count:
                is_weighted_bond = True
            good_nucs_each_pos.append(is_weighted_bond)

        weighted_structure:str = ''
        for nucIndex in range(structure_list.nuc_count):
            is_bonded = good_nucs_each_pos[nucIndex]
            new_counter: collections.Counter = collections.Counter(nuc_pairs_comp_list[nucIndex])
            most_common_char: str= '.'
            if is_bonded is True:
                #most_common_char = '|'
                new_char:str = new_counter.most_common(2)[0][0]
                length = len(new_counter.most_common(2))
                if new_char == '.' and length > 1:
                    #then get second most common
                    new_char = new_counter.most_common(2)[1][0]
                most_common_char = new_char
            weighted_structure = weighted_structure + most_common_char

        weighted_structure: Sara2SecondaryStructure = Sara2SecondaryStructure(sequence=structure_list.sara_stuctures[0].sequence,
                                                                                structure=weighted_structure,)

        return weighted_structure

    def compair_weighted_structure(self, unbound_mfe_struct:Sara2SecondaryStructure, bound_mfe_struct:Sara2SecondaryStructure, weighted_result:Sara2SecondaryStructure, nuc_count:int):
        """
        Compaire the weighted structure against the folded and not-folded mfe's.
        If a element is present in the folded mfe then it gets a '-'
        if element is in unbound only then it gets a '|'.
        The idea is that if you have a straight line in the list then it is very close to the
        folded mfe and if it is not straight then it is more like the unbound mfe.
        """
        unbound:str = '|'
        num_unbound:int = 0
        bound:str = '-'
        num_bound:int = 0
        both:str = '+'
        num_both:int = 0
        dot:str = '.'
        num_dot:int = 0
        compared_struct:str = ''            

        for nuc_index in range(nuc_count):
            weighted_nuc:str = weighted_result.structure[nuc_index]
            unbound_nuc:str = unbound_mfe_struct.structure[nuc_index]
            bound_nuc: str = bound_mfe_struct.structure[nuc_index]

            comp_nuc_symbol:str = ''

            if weighted_nuc == bound_nuc and weighted_nuc != unbound_nuc:
                comp_nuc_symbol = bound
                num_bound += 1
            elif weighted_nuc != bound_nuc and weighted_nuc == unbound_nuc:
                comp_nuc_symbol = unbound
                num_unbound += 1
            elif weighted_nuc == bound_nuc and weighted_nuc == unbound_nuc:
                comp_nuc_symbol = both
                num_both += 1
            else:
                comp_nuc_symbol = dot
                num_dot += 1

            weighted_nuc_counts:WeightedNucCounts = WeightedNucCounts(num_unbound=num_unbound,
                                                                        num_bound=num_bound,
                                                                        num_both=num_both,
                                                                        num_dot=num_dot,
                                                                        num_nucs=nuc_count
                                                                        )
            compared_struct = compared_struct + comp_nuc_symbol

            weighted_nuc_counts.num_nucs = nuc_count

        compared_data: WeightedComparisonResult = WeightedComparisonResult(comp_struct=compared_struct,
                                                                           unbound_mfe_struct=unbound_mfe_struct,
                                                                           bound_mfe_struct=bound_mfe_struct,
                                                                           weighted_nuc_counts=weighted_nuc_counts)    
        return compared_data

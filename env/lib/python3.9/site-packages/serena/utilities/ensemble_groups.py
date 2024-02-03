"""
File for the classes associated with the ensemble groups
"""
from typing import List, Dict
import attrs

from serena.utilities.ensemble_structures import (Sara2SecondaryStructure,
                                                  Sara2StructureList,
                                                  KcalRanges)

@attrs.define
class EnsembleSwitchStateMFEStructs():
    """
    Class that holds the secondary structures for representing
    the mfe and folded structures uses for LMV
    """
    non_switch_mfe_struct:Sara2SecondaryStructure= Sara2SecondaryStructure()
    switched_mfe_struct:Sara2SecondaryStructure = Sara2SecondaryStructure()

    def set_non_switch_mfe(self, kcal:float, struct:str):
        """
        Sets the unbound mfe secondary structure
        """
        self.non_switch_mfe_struct = Sara2SecondaryStructure(structure=struct,
                                                             free_energy=kcal)

    def set_switch_mfe(self, kcal:float, struct:str):
        """
        Sets the folded mfe secondary structure.
        """
        self.switched_mfe_struct = Sara2SecondaryStructure(structure=struct,
                                                             free_energy=kcal)

#need to turn into a dataclass or attrs.define
class SingleEnsembleGroup():#pylint: disable=too-many-instance-attributes
    """
    Datatype that represents and hold the info for a single
    energy span group in the ensemble
    """
    def __init__(self) -> None:
        self._group: Sara2StructureList = Sara2StructureList()
        self._switch_state_structures: EnsembleSwitchStateMFEStructs = None
        self._multi_state_mfe_struct: List[str] = []
        """
        0 is mfe for unbound and 1 is mfe for bound
        """
        self._multi_state_mfe_kcal: List[float] = []
        self._kcal_start: float = 0
        self._kcal_end: float = 0
        self._kcal_span: float = 0

    @property
    def group(self)->Sara2StructureList:
        """
        Return the ensemble group as a sara2structurelist
        """
        return self._group

    @group.setter
    def group(self, the_group:Sara2StructureList):
        """
        Set the ensemble group using a Sara2SecondaryStructure
        """
        self._group = the_group

    @property
    def multi_state_mfe_struct(self)->List[str]:
        """
        Return the multi state mfe struct that holds
        the unbound and bound structs
        """
        return self._multi_state_mfe_struct

    @multi_state_mfe_struct.setter
    def multi_state_mfe_struct(self, structs:List[str]):
        """
        Sets teh multi state struct
        """
        self._multi_state_mfe_struct = structs

    def append_multi_state_mfe_data(self, structure: str, kcal: float):
        """
        Appends structures to the multo state mfe data
        """
        self._multi_state_mfe_struct.append(structure)
        self._multi_state_mfe_kcal.append(kcal)

    @property
    def multi_state_mfe_kcal(self)->List[float]:
        """
        Return the multi state mfe kcal list
        """
        return self._multi_state_mfe_kcal

    @multi_state_mfe_kcal.setter
    def multi_state_mfe_kcal(self, kcals:List[float]):
        """
        Set the multi state mfe kcal list
        """
        self._multi_state_mfe_kcal = kcals

    @property
    def kcal_span(self)->float:
        """
        Return the kcal span of the group
        """
        return self._kcal_span

    @kcal_span.setter
    def kcal_span(self, kcal:float):
        """
        Set the kcal span of the group
        """
        self._kcal_span = kcal

    @property
    def kcal_start(self)->float:
        """
        Return the kcal start of the group
        """
        return self._kcal_start

    @kcal_start.setter
    def kcal_start(self, kcal:float):
        """
        Set the kcal start of the group
        """
        self._kcal_start = kcal

    @property
    def kcal_end(self)->float:
        """
        Return the kcal end of the group
        """
        return self._kcal_end

    @kcal_end.setter
    def kcal_end(self, kcal:float):
        """
        Set the kcal end of the group
        """
        self._kcal_end = kcal

    def update_kcals(self, start:float, stop:float, span:float):
        "Update all the kcal properties at once"
        self._kcal_start = start
        self._kcal_end = stop
        self._kcal_span = span

    @property
    def switch_state_structures(self)->EnsembleSwitchStateMFEStructs:
        """
        Return the switch state structures that hold the 
        unbound and bound mfe structures
        """
        return self._switch_state_structures

    @switch_state_structures.setter
    def switch_state_structures(self, structs:EnsembleSwitchStateMFEStructs):
        """
        Set the switch state structures that hold the 
        unbound and bound mfe structures
        """
        self._switch_state_structures = structs

class MultipleEnsembleGroups():
    """"
    Multiple Ensemble Groups class
    """
    def __init__(self, switch_state_structures: EnsembleSwitchStateMFEStructs) -> None:
        self._groups: List[SingleEnsembleGroup] = []
        self._raw_groups: List[Sara2StructureList] = []
        self._switch_state_structures: EnsembleSwitchStateMFEStructs = switch_state_structures
        self._groups_dict: Dict[int, Sara2StructureList] = {}
        self._group_values: List[float] = []
        self._num_groups: int = 0
        self._group_kcal_ranges: List[KcalRanges] =  []

    @property
    def switch_state_structures(self)->EnsembleSwitchStateMFEStructs:
        """
        Return the switch state structures for holdng the unbound and 
        bound mfe structs for analysis
        """
        return self._switch_state_structures

    @property
    def num_groups(self)->int:
        """
        Return the number of energy groups in the ensemble
        """
        return self._num_groups

    def add_group(self, group:SingleEnsembleGroup):
        """
        Prefered way to add a ensemble group to the multi group class
        """
        self._groups.append(group)
        self._raw_groups.append(group.group)
        self._groups_dict[self._num_groups]= group.group
        self._group_values.append(group.kcal_start)
        kcal_range: KcalRanges = KcalRanges(start=group.kcal_start, stop=group.kcal_end)
        self._group_kcal_ranges.append(kcal_range)
        self._num_groups = self._num_groups + 1

    @property
    def groups(self)->List[SingleEnsembleGroup]:
        """
        Return the ensemble groupds as a list of 
        SingleEnsembleGroups
        """
        return self._groups

    @property
    def raw_groups(self)->List[Sara2StructureList]:
        """
        Return the raw ensembled groups as a list of 
        Sara2StructureLists
        """
        return self._raw_groups

    @property
    def non_switch_state_structure(self)->Sara2SecondaryStructure:
        """
        Return the unbound mfe secondary structure
        """
        return self._switch_state_structures.non_switch_mfe_struct

    @property
    def switched_state_structure(self)->Sara2SecondaryStructure:
        """
        Return the folded mfe secondary structure
        """
        return self._switch_state_structures.switched_mfe_struct

    @property
    def groups_dict(self)->Dict[int, Sara2StructureList]:
        """
        Return the groups dict for values and structures
        """
        return self._groups_dict

    @property
    def group_values(self)->List[float]:
        """
        Return the group values as a list of floats
        """
        return self._group_values

    @property
    def group_kcal_ranges(self)->List[KcalRanges]:
        """
        Return the list of kcal ranges
        """
        return self._group_kcal_ranges

    @property
    def total_structures(self)->int:
        """
        Return the total number of structures
        """
        total:int = 0
        for group in self.raw_groups:
            total += group.num_structures
        return total

class MakeEnsembleGroups():
    """
    Class for generating the ensemble groups consumed by serena and sara
    """

    def make_switch_mfe_states_from_secondary_strucures(self, switched_mfe_struc:Sara2SecondaryStructure, non_switch_mfe_struct:Sara2SecondaryStructure):#pylint: disable=line-too-long
        """
        Make switch states
        """
        return EnsembleSwitchStateMFEStructs(non_switch_mfe_struct=non_switch_mfe_struct,
                                             switched_mfe_struct=switched_mfe_struc)

    def make_singel_ensemble_group(self, ensemble_structures:Sara2StructureList, mfe_switch_structures:EnsembleSwitchStateMFEStructs, kcal_start:float, kcal_end:float):#pylint: disable=line-too-long
        """
        Function to make a single ensemble group from a sara2structure list
        """
        single_ensemble_group:SingleEnsembleGroup = SingleEnsembleGroup()
        single_ensemble_group.group = ensemble_structures
        single_ensemble_group.switch_state_structures = mfe_switch_structures
        single_ensemble_group.kcal_start = kcal_start
        single_ensemble_group.kcal_end = kcal_end
        single_ensemble_group.kcal_span = kcal_end - kcal_start
        return single_ensemble_group

    def make_multiple_ensemple_groups(self, ensemble_groups:List[SingleEnsembleGroup], mfe_switch_structures:EnsembleSwitchStateMFEStructs):#pylint: disable=line-too-long
        """
        Function to make a multi ensemble group from a list of single ensemble goups
        """
        multi_group:MultipleEnsembleGroups = MultipleEnsembleGroups(switch_state_structures=mfe_switch_structures)#pylint: disable=line-too-long
        for group in ensemble_groups:
            multi_group.add_group(group=group)
        return multi_group

"""
Sara2 api for accessing and manipulating secondary structures 
in dot parenthisis form
copyright 2023 GrizzlyEngineer
"""
from typing import List
from dataclasses import dataclass

@dataclass
class KcalRanges():
    """
    Class to hold kcal ranges
    """
    start: float = 0
    stop: float = 0

class Sara2SecondaryStructure():
    """
    Sara 2 Secondary Structure that is used to hold all the info for
    each secondary structure in ensemble
    """

    def __init__(self, sequence:str = '', structure: str = '', free_energy: float = 0, stack_energy: float = 0) -> None:#pylint: disable=line-too-long
        self._sequence: str = sequence
        self._structure: str = structure
        self._free_energy: float = free_energy
        self._stack_energy: float = stack_energy

    @property
    def sequence(self):
        """
        Returns the sequence as a string
        """
        return self._sequence

    @sequence.setter
    def sequence(self, primary_struc: str):
        """
        Sets the sequence using string
        """
        self._sequence = primary_struc

    @property
    def structure(self):
        """
        Returns the secondary strucuture in dot parens notation
        """
        return self._structure

    @structure.setter
    def structure(self, dot_parens: str):
        """
        Sets the secondary structure using dot parense notation string
        """
        self._structure = dot_parens

    @property
    def free_energy(self):
        """
        Returns the total free energy as float
        """
        return self._free_energy

    @free_energy.setter
    def free_energy(self, energy: float):
        """
        Sets the total free energy with float
        """
        self._free_energy = energy

    @property
    def stack_energy(self):
        """
        Returns the stack energy as float
        """
        return self._stack_energy

    @stack_energy.setter
    def stack_energy(self, energy: float):
        """
        Sets the stack energy with float
        """
        self._stack_energy = energy

    @property
    def nuc_count(self):
        """
        Returns the number of nucleotides as a int
        """
        return len(self._sequence)


class Sara2StructureList():#pylint: disable=too-many-instance-attributes
    """
    Sara2 Structure List that holds all the Sar2SecondaryStructurs
    that represent the ensemble in raw form
    """
    def __init__(self) -> None:
        self._sara_structures_list: List[Sara2SecondaryStructure] = []
        self._structures: List[str] = []
        self._free_energy_list: list[float] = []
        self._stack_energy_list: list[float] = []
        self._min_free_energy: float = 0
        self._max_free_energy: float = 0
        self._min_stack_energy: float = 0
        self._max_stack_energy: float = 0
        self._num_structures: int = 0
        self._free_energy_span:float = 0
        self._stack_energy_span:float = 0
        self._weighted_structure:str = ''

    def process_energy(self):
        """
        Process min and max energies in list as well
        as populate counts. It always ran after adding 
        structure
        """
            #now populate min and max
        #do free energy
        if len(self._free_energy_list) == 0:
            self._min_free_energy = 0
            self._max_free_energy = 0
        else:
            self._min_free_energy = min(self._free_energy_list)
            self._max_free_energy = max(self._free_energy_list)

        self._free_energy_span = self._max_free_energy - self._min_free_energy
        #do stack energy

        if len(self._stack_energy_list) == 0:
            self._min_stack_energy = 0
            self._max_stack_energy = 0
        else:
            self._min_stack_energy = min(self._stack_energy_list)
            self._max_stack_energy = max(self._stack_energy_list)
        self._stack_energy_span = self._max_stack_energy - self._min_stack_energy

        #now count
        self._num_structures = len(self._sara_structures_list)

    def add_structure(self, structure: Sara2SecondaryStructure):
        """
        main way to add a structure to the list
        """
        self._sara_structures_list.append(structure)
        #self._structures.append(structure.structure)
        self._free_energy_list.append(structure.free_energy)
        self._stack_energy_list.append(structure.stack_energy)
        #self.process_energy()

    def remove_structure(self, index:int):
        """
        remove a structure from memory
        """
        del self._structures[index]
        del self._free_energy_list[index]
        del self._stack_energy_list[index]
        #self.process_energy()

    @property
    def mfe_structure(self):
        """
        Returns the mfe secibdary structure as a string
        """
        structure:str = ''
        if len(self.sara_stuctures) > 0:
            structure = self.sara_stuctures[0].structure
        return structure

    @property
    def mfe_free_energy(self):
        """
        Returns the mfe total free energy as float
        """
        self.process_energy()
        energy: float = 0
        if len(self.sara_stuctures) > 0:
            energy = self.sara_stuctures[0].free_energy
        return energy

    @property
    def mfe_stack_energy(self):
        """
        Returns the mfe stack energy as float
        """
        self.process_energy()
        energy: float = 0
        if len(self.sara_stuctures) > 0:
            energy = self.sara_stuctures[0].stack_energy
        return energy

    @property
    def nuc_count(self):
        """
        Returns the total number of nucleotides as int
        """
        count: int = 0
        if len(self.sara_stuctures) > 0:
            count = self.sara_stuctures[0].nuc_count
        return count

    @property
    def sara_stuctures(self):
        """
        Returns the sara structures that make up list
        """
        return self._sara_structures_list

    @sara_stuctures.setter
    def sara_stuctures(self, structs_list: List[Sara2SecondaryStructure]):
        """
        Sets the sara structures list using a List of Sara2Structures
        """
        #reset list
        self._sara_structures_list=[]
        #fill it in now
        for struc in structs_list:
            self.add_structure(struc)

    @property
    def max_free_energy(self):
        """
        Returns the maximum free energy of the structures in the list
        """
        self.process_energy()
        return self._max_free_energy

    @property
    def min_free_energy(self):
        """
        Returns the minimum free energy of the structures in the list
        """
        self.process_energy()
        return self._min_free_energy

    @property
    def max_stack_energy(self):
        """
        Returns the maximum stack energy of the structures in the list
        """
        self.process_energy()
        return self._max_stack_energy

    @property
    def min_stack_energy(self):
        """
        Returns the minimum stack energy of the structures in the list
        """
        self.process_energy()
        return self._min_stack_energy

    @property
    def num_structures(self):
        """
        Returns the number of structures in the list
        """
        self.process_energy()
        return self._num_structures

    @property
    def free_energy_span(self):
        """
        Returns the span of the free energy of the structures in the list
        """
        self.process_energy()
        return self._free_energy_span

    @property
    def stack_energy_span(self):
        """
        Returns the span of the stack energy of the structures in the list
        """
        self.process_energy()
        return self._stack_energy_span

    @property
    def weighted_structure(self):
        """
        Returns the weighted structure as a string
        """
        return self._weighted_structure

    @weighted_structure.setter
    def weighted_structure(self, structure: str):
        """
        sets the weigthed structure
        """
        self._weighted_structure = structure

class MakeSecondaryStructures():
    """
    Class to genereate the secondary structure
    framework used by serena and sara
    """
    def make_secondary_structure(self, primary_structure:str, secondary_structure:str, free_energy:float, stack_free_energy:float)->Sara2SecondaryStructure:#pylint: disable=line-too-long
        """
        Function to make a secondary structue
        """
        return Sara2SecondaryStructure(sequence=primary_structure,
                                       structure=secondary_structure,
                                       free_energy=free_energy,
                                       stack_energy=stack_free_energy
                                       )

    def make_secondary_strucuture_list(self, secondary_structures_list: List[Sara2SecondaryStructure])->Sara2StructureList:#pylint: disable=line-too-long
        """
        Function to make a secondary structure list
        """
        structure_list:Sara2StructureList = Sara2StructureList()
        for structure in secondary_structures_list:
            structure_list.add_structure(structure)
        return structure_list

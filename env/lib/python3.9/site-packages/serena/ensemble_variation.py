"""
This is the file to just get plain old ensemble variation for a list
of structures that make up an ensemble
"""

from serena.utilities.ensemble_structures import Sara2SecondaryStructure, Sara2StructureList
from serena.utilities.ensemble_variation import EV, EnsembleVariation
from serena.interfaces.nupack4_0_28_wsl2_interface import NUPACK4Interface, MaterialParameter

class RunEnsembleVariation(EnsembleVariation):
    """
    Class to get the ensemble variation
    """
    def ev_from_structures_list(self, structures_list:Sara2StructureList, mfe_structure:Sara2SecondaryStructure):#pylint: disable=line-too-long
        """
        Get the Ensemble Variation for a sara2Structure List amd reference mfe structure
        """
        e_v:EV = self.ensemble_variation_algorithm(kcal_group_structures_list=structures_list,
                                                        ref_structure=mfe_structure)
        return e_v.ev_normalized

    def ev_from_nupack4(self, material:MaterialParameter, temp_c:int, span_from_mfe:int, sequence:str):#pylint: disable=line-too-long
        """
        Get the Ensemble Variation for an ensemble after folding in silico with nupack
        """
        nupack4: NUPACK4Interface = NUPACK4Interface()
        structs:Sara2StructureList = nupack4.get_subopt_energy_gap(material_param=material,
                                  temp_C=temp_c,
                                  sequence_string=sequence,
                                  energy_delta_from_MFE=span_from_mfe,
                                  )
        ensemble_variation:float = self.ev_from_structures_list(structures_list=structs,
                                                                mfe_structure=structs.sara_stuctures[0])#pylint: disable=line-too-long
        return ensemble_variation

"""
Main entry for the algorithm to determine the swithchyness 
of a sequence baed on analysis of the ensemble
"""

from pathlib import Path
from typing import List
import os
import sys

from serena.analysis.ensemble_analysis import InvestigateEnsemble, InvestigateEnsembleResults
from serena.interfaces.nupack4_0_28_wsl2_interface import MaterialParameter, NUPACK4Interface
from serena.utilities.ensemble_structures import Sara2SecondaryStructure, Sara2StructureList
from serena.utilities.ensemble_groups import MultipleEnsembleGroups, EnsembleSwitchStateMFEStructs

# from serena.bin.ba import ArchiveSecondaryStructureList
# import serena.bin.backup_serena

class RunInvestigateEnsemble(InvestigateEnsemble):
    """
    Class that is the main entry point for ensemble investigation
    """
    def get_and_store_ensemble_data(self,sequence:str, material_param:MaterialParameter, temp_c: int, kcal_span_from_mfe:int, backup_folder:Path, record_name:str):
        nupack4:NUPACK4Interface = NUPACK4Interface()
        structs:Sara2StructureList = nupack4.get_subopt_energy_gap(material_param=material_param,
                                                                    temp_C=temp_c,
                                                                    sequence_string=sequence,
                                                                    energy_delta_from_MFE=kcal_span_from_mfe,
                                                                    )
       
        # backup_records:ArchiveSecondaryStructureList = ArchiveSecondaryStructureList(working_folder=backup_folder,
        #                                      var_name=record_name,
        #                                      use_db=True)
        
        # # backup_records.structs.mfe_free_energy = structs.mfe_free_energy
        # # backup_records.structs.mfe_structure = structs.mfe_structure
        # # backup_records.structs.mfe_stack_energy = structs.mfe_stack_energy
        # # backup_records.structs.nuc_count = structs.nuc_count                    
        # backup_records.structs.sara_stuctures = structs.sara_stuctures
        # # backup_records.structs.max_free_energy = structs.max_free_energy
        # # backup_records.structs.min_free_energy = structs.min_free_energy
        # # backup_records.structs.max_stack_energy = structs.max_stack_energy
        # # backup_records.structs.min_stack_energy = structs.min_stack_energy
        # # backup_records.structs.num_structures = structs.num_structures
        # # backup_records.structs.free_energy_span = structs.free_energy_span
        # # backup_records.structs.stack_energy_span = structs.stack_energy_span
        # backup_records.structs.weighted_structure = structs.weighted_structure

    def pull_archived_data(self, data_folder:str, design_id:str)->Sara2StructureList:
        
        
        
        # directories = [f for f in os.listdir(path=data_folder) if os.path.isdir(f)]
        
        target_folder:str = f'{data_folder}/{design_id}'
        
        if os.path.isdir(target_folder) == False:
            raise Exception(f'Path {target_folder} does not exist')
        
        structs:Sara2StructureList = Sara2StructureList()
        
        # backup_records:ArchiveSecondaryStructureList = ArchiveSecondaryStructureList(working_folder=data_folder,
        #                                     var_name=design_id,
        #                                     use_db=True)
        
        # # structs.mfe_free_energy = backup_records.structs.mfe_free_energy
        # #structs.mfe_structure = backup_records.structs.mfe_structure
        # #structs.mfe_stack_energy = backup_records.structs.mfe_stack_energy
        # #structs.nuc_count = backup_records.structs.nuc_count             
        # structs.sara_stuctures = backup_records.structs.sara_stuctures
        # #structs.max_free_energy = backup_records.structs.max_free_energy
        # #structs.min_free_energy = backup_records.structs.min_free_energy
        # #structs.max_stack_energy = backup_records.structs.max_stack_energy
        # #structs.min_stack_energy = backup_records.structs.min_stack_energy
        # #structs.num_structures = backup_records.structs.num_structures
        # #structs.free_energy_span = backup_records.structs.free_energy_span
        # #structs.stack_energy_span = backup_records.structs.stack_energy_span
        # structs.weighted_structure = backup_records.structs.weighted_structure
        
        return structs

        
        
    def investigate_and_score_ensemble_nupack(self,sequence:str, folded_referenec_struct:str, material_param:MaterialParameter, temp_c: int, kcal_span_from_mfe:int, kcal_unit_increments: float = 1, aggressive:bool= False)->InvestigateEnsembleResults:
        """
        Use the nupack folding enginer to generate a MultipleEnsembleGroups from a sequence 
        and refence folded structure (folded mfe maybe?) and analyze it for switchyness score
        """
        nupack4:NUPACK4Interface = NUPACK4Interface()
        structs:Sara2StructureList = nupack4.get_subopt_energy_gap(material_param=material_param,
                                                                    temp_C=temp_c,
                                                                    sequence_string=sequence,
                                                                    energy_delta_from_MFE=kcal_span_from_mfe,
                                                                    )
        switch_states:EnsembleSwitchStateMFEStructs = EnsembleSwitchStateMFEStructs(non_switch_mfe_struct=structs.sara_stuctures[0],
                                                                                    switched_mfe_struct=Sara2SecondaryStructure(sequence=sequence,
                                                                                                                                structure=folded_referenec_struct))
        ensemble:MultipleEnsembleGroups = nupack4.load_nupack_subopt_as_ensemble(span_structures=structs,
                                                                                    kcal_span_from_mfe=kcal_span_from_mfe,
                                                                                    Kcal_unit_increments=kcal_unit_increments,
                                                                                    switch_state=switch_states
                                                                                    )
        results:InvestigateEnsembleResults = self.investigate_and_score_ensemble(ensemble=ensemble,
                                                                                 is_aggressive=aggressive)
        return results
"""
File to hold the main entry point for determining if a rna
sequence's ensemble is accepting of a switch
"""

from dataclasses import dataclass
from pathlib import Path
import os
import statistics
import openpyxl
import pandas as pd
from pandas import DataFrame
from typing import List, Dict
import time

from serena.interfaces.Sara2_API_Python3 import Sara2API, puzzleData
from serena.interfaces.vienna2_fmn_hack_interface import Vienna2FMNInterface
from serena.analysis.ensemble_analysis import InvestigateEnsemble, InvestigateEnsembleResults
from serena.interfaces.nupack4_0_28_wsl2_interface import NUPACK4Interface, MaterialParameter, NupackSettings, EnsembleSwitchStateMFEStructs
from serena.utilities.weighted_structures import WeightedStructure
from serena.utilities.ensemble_groups import SingleEnsembleGroup, MultipleEnsembleGroups
from serena.utilities.logging_serena import PNASAnalysisLogging
from serena.utilities.ensemble_structures import Sara2SecondaryStructure, Sara2StructureList
from serena.switchyness_score import RunInvestigateEnsemble



class SwitchAccetance():

    def __init__(self) -> None:
        pass

    def run_eterna_pnas(self):

        vienna2_fmn_hack: Vienna2FMNInterface = Vienna2FMNInterface()

        details:str= 'all'#f'20k_filtered_weighted_100K_gtrequal2_nucpenalty_run_1ish'
        pnas_round101_sheet:str = 'R101 Filtered good bad'
        #Round 7 (R101)
        same_state:str='3'
        sublab_name:str = "good"#f'Same State NG {same_state}'
        is_aggressive:bool = False
        save_title:str = sublab_name + "_open"
        run_name:str = "data_nut_test"#f'SSNG{same_state}_{details}'


        pnas_path:str = '/home/rnauser/test_data/pnas_testing_tweak.xlsx'
        timestr = time.strftime("%Y%m%d-%H%M%S")
        save_path:str = f'/home/rnauser/test_data/run_data/{run_name}/pnas_eternacon_{timestr}.xlsx'
        
        new_sara:Sara2API = Sara2API()
        puzzle_data: puzzleData
        pandas_sheet: DataFrame
        puzzle_data, pandas_sheet = new_sara.ProcessLab(path=pnas_path,
                                                      designRound_sheet=pnas_round101_sheet,
                                                      sublab_name=sublab_name
                                                      )

       

        predicted_foldchange_list: List[float] = []
        avg_raw_score_list: List[float] = []
        avg_num_structures_list: List[float] = []

        raw_score_36_list: List[float] = []
        raw_score_37_list: List[float] = []
        raw_score_38_list: List[float] = []

        num_structures_36_list: List[float] = []
        num_structures_37_list: List[float] = []
        num_structures_38_list: List[float] = []
        nupack_4:NUPACK4Interface = NUPACK4Interface()
        scoreing:InvestigateEnsemble = InvestigateEnsemble()
        
        flag:int =0
        for design in puzzle_data.designsList:
            design_id= str(design.design_info.DesignID)
            sequence = design.design_info.Sequence
            fold_change = design.wetlab_results.FoldChange
            eterna_score = design.wetlab_results.Eterna_Score
            folding_subscore = design.wetlab_results.Folding_Subscore
            switch_subscore = design.wetlab_results.Switch_Subscore
            baseline_subscore = design.wetlab_results.Baseline_Subscore

            #make a new line of just this designs row
            
            
            do_weighted:bool = False
            struct_to_use:Sara2SecondaryStructure = ''
            if do_weighted is True:
            #this is the fmn bound mfe struct, subopt list and weighted struck
                fmn_subopt = vienna2_fmn_hack.rnasubopt_fmn(input_sequence=sequence,
                                                            do_fmn=True)
                fmn_weighted_struct: WeightedStructure = WeightedStructure()
                struct_to_use= fmn_weighted_struct.make_weighted_struct(fmn_subopt)
            else:
                struct_to_use = vienna2_fmn_hack.rnafold_fmn(input_sequence=sequence,
                                                            do_fmn=True)

            material:MaterialParameter = MaterialParameter.rna06_nupack4
            temp:int = 37
            kcal_span:int = 7
            kcal_unit_increments:float = 1

            #nupack_settings:NupackSettings = NupackSettings(material_param=material,
            #                                                temp_C=temp,
            #                                                kcal_delta_span_from_mfe=kcal_span,
            #                                                Kcal_unit_increments=kcal_unit_increments,
            #                                                sequence=sequence,
            #                                                folded_2nd_state_structure=struct_to_use,
            #                                                folded_2nd_state_kcal=0
            #                                                )
            
            #do 37 deg first
            archive_folder:Path = Path('/mnt/s/Sync/RNA/serena_data_nut/pnas_test_good/')
            investigate:RunInvestigateEnsemble = RunInvestigateEnsemble()
            # investigate.get_and_store_ensemble_data(sequence=sequence,
            #                                         temp_c=temp,
            #                                         material_param=material,
            #                                         kcal_span_from_mfe=kcal_span,
            #                                         backup_folder=archive_folder,
            #                                         record_name=design_id)
            
            sara2_list:Sara2StructureList = investigate.pull_archived_data(data_folder=archive_folder, design_id=design_id)
            
            # sara2_list:Sara2StructureList = nupack_4.get_subopt_energy_gap(material_param=material,
            #                                                                temp_C=temp,
            #                                                                sequence_string=sequence,
            #                                                                energy_delta_from_MFE=kcal_span)
         
            if sara2_list.num_structures < 500000:
                reference_structures:EnsembleSwitchStateMFEStructs = EnsembleSwitchStateMFEStructs(switched_mfe_struct=struct_to_use,
                                                                                                non_switch_mfe_struct=sara2_list.sara_stuctures[0])
                
                ensemble_groups: MultipleEnsembleGroups = nupack_4.load_nupack_subopt_as_ensemble(span_structures=sara2_list,
                                                                                                kcal_span_from_mfe=kcal_span,
                                                                                                Kcal_unit_increments=kcal_unit_increments,
                                                                                                switch_state=reference_structures
                                                                                                )
                
                investigation_results:InvestigateEnsembleResults = scoreing.investigate_and_score_ensemble(ensemble=ensemble_groups,
                                                                                                           is_aggressive=is_aggressive)
                
                total_scores: float = 0
                #if investigation_results.basic_scores.total_score > 0:
                total_scores: float = investigation_results.basic_scores.total_score + investigation_results.advanced_scores.total_score
                #else:
                #    total_scores = 0 - investigation_results.advanced_scores.excess_struct_penalty
                
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'SerenaTotalScore'] = total_scores
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'SerenaTotalScore_NoExcessStructs'] = total_scores + investigation_results.advanced_scores.excess_struct_penalty
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Basic_score'] = investigation_results.basic_scores.total_score
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Functional_score'] = investigation_results.basic_scores.functional_switch_score
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Powerfull_score'] = investigation_results.basic_scores.powerful_switch_score
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Bonuses'] = investigation_results.basic_scores.bonuses
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'On_off_score'] = investigation_results.basic_scores.on_off_switch_score
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Penalties'] = investigation_results.basic_scores.penalties
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Advanced_Score'] = investigation_results.advanced_scores.total_score
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Advanced_comp_bonus'] = investigation_results.advanced_scores.comp_bonus
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Advanced_comp_penalty'] = investigation_results.advanced_scores.comp_penalty
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Advanced_lmv_bonus'] = investigation_results.advanced_scores.lmv_bonus
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Advanced_lmv_penalty'] = investigation_results.advanced_scores.lmv_penalty
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Structure_Penalty'] = investigation_results.advanced_scores.excess_struct_penalty
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'NumStructs'] =  investigation_results.number_structures
                
                
            else:
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'SerenaTotalScore'] = -100
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'SerenaTotalScore_NoExcessStructs'] = -100
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Basic_score'] = -100
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Functional_score'] = -100
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Powerfull_score'] = -100
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Bonuses'] = -100
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'On_off_score'] = -100
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Penalties'] = -100
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Advanced_Score'] = -100
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Advanced_comp_bonus'] = -100
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Advanced_comp_penalty'] = -100
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Advanced_lmv_bonus'] = -100
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Advanced_lmv_penalty'] = -100
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Structure_Penalty'] = -100
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Advanced_Score'] = -100
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'Structure_Penalty'] = 100
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'NumStructs'] =  sara2_list.num_structures
            
            design_data_df:DataFrame = pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID]
            logging: PNASAnalysisLogging = PNASAnalysisLogging()
            if os.path.isdir(f'/home/rnauser/test_data/run_data/{run_name}') == False:
                os.makedirs(f'/home/rnauser/test_data/run_data/{run_name}')
            logging.save_excel_sheet(design_data_df, save_path, save_title)

switch_acceptance:SwitchAccetance = SwitchAccetance()
switch_acceptance.run_eterna_pnas()
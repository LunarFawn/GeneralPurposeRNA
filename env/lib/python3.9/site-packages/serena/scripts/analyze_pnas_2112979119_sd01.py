"""
File for managing analysis of dataset pnas.2112979119.sd01 from https://www.pnas.org/doi/full/10.1073/pnas.2112979119
"""

from enum import Enum
from typing import List, Union
from pathlib import Path
from pandas import DataFrame
import os
from dataclasses import dataclass
import argparse

from serena.utilities.ensemble_structures import Sara2SecondaryStructure, Sara2StructureList
from serena.utilities.weighted_structures import WeightedStructure

from serena.interfaces.Sara2_API_Python3 import Sara2API, puzzleData, DesignPerformanceData, DesignInformation
from serena.interfaces.vienna2_fmn_hack_interface import Vienna2FMNInterface
from serena.interfaces.nupack4_0_28_wsl2_interface import (
    MaterialParameter, 
    NUPACK4Interface,
    NupackSettings
)

from serena.bin.backup_serena_v2 import ArchiveSecondaryStructureList
from serena.bin.backup_investigator_v1 import ArchiveInvestigator

from serena.interfaces.nupack4_0_28_wsl2_interface import NUPACK4Interface, MaterialParameter, NupackSettings, EnsembleSwitchStateMFEStructs
from serena.analysis.ensemble_analysis import InvestigateEnsemble, InvestigateEnsembleResults, InvestigatorResults, ReferenceStructures
from serena.utilities.ensemble_groups import SingleEnsembleGroup, MultipleEnsembleGroups
from serena.utilities.logging_serena import PNASAnalysisLogging
from serena.analysis.investigator import (
    ComparisonEvalResults,
    ComparisonNucResults,
    ComparisonLMVResponse,
    LMVAssertionResult
)


class ArchiveFlow(Enum):
    GET="GET"
    PUT="PUT"

@dataclass
class ArchiveData():
    design_info:DesignPerformanceData = None
    nupack_settings:NupackSettings = None
    structs:Sara2StructureList = None
    fmn_folded_mfe:Sara2SecondaryStructure = None
    fmn_folded_weighted:Sara2SecondaryStructure = None
    
@dataclass
class ArchiveInvestigatorData():
    investigator: InvestigateEnsembleResults = None
    design_info:DesignPerformanceData = None
    
class ProcessPNAS():
    
    def __init__(self) -> None:
        self._vienna2_fmn_hack: Vienna2FMNInterface = Vienna2FMNInterface()
        self._nupack4:NUPACK4Interface = NUPACK4Interface()
        self._scoreing:InvestigateEnsemble = InvestigateEnsemble()
    
    @property
    def vienna2_fmn_hack(self)->Vienna2FMNInterface:
        return self._vienna2_fmn_hack
    
    @property
    def nupack4(self)->NUPACK4Interface:
        return self._nupack4

    @property
    def scoring(self)->InvestigateEnsemble:
        return self._scoreing
    
    def get_fmn_state_fold(self, sequence:str, do_weighted:bool)->Sara2SecondaryStructure:
        
        #  = False
        struct_to_use:Sara2SecondaryStructure = ''
        if do_weighted is True:
        #this is the fmn bound mfe struct, subopt list and weighted struck
            fmn_subopt = self.vienna2_fmn_hack.rnasubopt_fmn(input_sequence=sequence,
                                                        do_fmn=True)
            fmn_weighted_struct: WeightedStructure = WeightedStructure()
            struct_to_use= fmn_weighted_struct.make_weighted_struct(fmn_subopt)
        else:
            struct_to_use = self.vienna2_fmn_hack.rnafold_fmn(input_sequence=sequence,
                                                        do_fmn=True)
        return struct_to_use

        
    def record_nupack_ensemble_structs(self, pnas_dataset_path:Path, round:str, sublab:str, nupack_settings:NupackSettings, archive_path:str=None)->List[Sara2StructureList]:
        new_sara:Sara2API = Sara2API()
        puzzle_data: puzzleData
        pandas_sheet: DataFrame
        puzzle_data, pandas_sheet = new_sara.ProcessLab(path=pnas_dataset_path.as_posix(),
                                                      designRound_sheet=round,
                                                      sublab_name=sublab
                                                      )
        pnas_data:List[Sara2StructureList] = []
        for design in puzzle_data.designsList:


            # design_id= str(design.design_info.DesignID)
            sequence = design.design_info.Sequence
            # fold_change = design.wetlab_results.FoldChange
            # eterna_score = design.wetlab_results.Eterna_Score
            # folding_subscore = design.wetlab_results.Folding_Subscore
            # switch_subscore = design.wetlab_results.Switch_Subscore
            # baseline_subscore = design.wetlab_results.Baseline_Subscor    
            
            
            structs:Sara2StructureList = self.nupack4.get_subopt_energy_gap(material_param=nupack_settings.material_param,
                                                                    temp_C=nupack_settings.temp_C,
                                                                    sequence_string=sequence,
                                                                    energy_delta_from_MFE=nupack_settings.kcal_span_from_mfe,
                                                                    )
            pnas_data.append(structs)
            fmn_mfe:Sara2SecondaryStructure = self.get_fmn_state_fold(sequence=sequence,
                                                                    do_weighted=False)
            fmn_weighted_struct:Sara2SecondaryStructure = self.get_fmn_state_fold(sequence=sequence,
                                                                                do_weighted=True)
            
            # if do_archive is True:
            if os.path.isdir(archive_path) is False:
                raise FileExistsError(f'File {archive_path} is not a valid path')
            archive_data:ArchiveData = ArchiveData(design_info=design,
                                                    nupack_settings=nupack_settings,
                                                    structs=structs,
                                                    fmn_folded_mfe=fmn_mfe,
                                                    fmn_folded_weighted=fmn_weighted_struct)
            
            self.archive_ensemble_data(dest_folder=archive_path,
                                        flow=ArchiveFlow.PUT,data=archive_data)
            
        print("Done!")

    def perform_investigation_computations(self, pnas_dataset_path:Path, round:str, sublab:str, source_archive_path:str, target_archive_path:str, do_weighted:bool, max_num_structs:int = 500000, is_agressive:bool = False):
        new_sara:Sara2API = Sara2API()
        puzzle_data: puzzleData
        pandas_sheet: DataFrame
        puzzle_data, pandas_sheet = new_sara.ProcessLab(path=pnas_dataset_path.as_posix(),
                                                      designRound_sheet=round,
                                                      sublab_name=sublab
                                                      )
        # pnas_data:List[Sara2StructureList] = []
        for design in puzzle_data.designsList:
            if os.path.isdir(source_archive_path) is False:
                raise FileExistsError(f'File {source_archive_path} is not a valid path')
            
            temp_archive:ArchiveData = ArchiveData(design_info=design)
            found_data:ArchiveData = self.archive_ensemble_data(dest_folder=source_archive_path,
                                        flow=ArchiveFlow.GET,
                                        data=temp_archive)
            
            if found_data.structs.num_structures < max_num_structs:
                struct_to_use:Sara2SecondaryStructure = found_data.fmn_folded_mfe
                if do_weighted is True:
                    struct_to_use = found_data.fmn_folded_weighted
                    
                reference_structures:EnsembleSwitchStateMFEStructs = EnsembleSwitchStateMFEStructs(switched_mfe_struct=struct_to_use,
                                                                                                non_switch_mfe_struct=found_data.structs.sara_stuctures[0])
                
                ensemble_groups: MultipleEnsembleGroups = self.nupack4.load_nupack_subopt_as_ensemble(span_structures=found_data.structs,
                                                                                                kcal_span_from_mfe=found_data.nupack_settings.kcal_span_from_mfe,
                                                                                                Kcal_unit_increments=found_data.nupack_settings.Kcal_unit_increments,
                                                                                                switch_state=reference_structures
                                                                                                )
                
                investigation_results:InvestigateEnsembleResults = self.scoring.investigate_and_score_ensemble(ensemble=ensemble_groups,
                                                                                                           is_aggressive=is_agressive)
                
                data_to_archive:ArchiveInvestigatorData = ArchiveInvestigatorData(investigator=investigation_results,
                                                                          design_info=design)
                
                self.archive_investigation_computations(dest_folder=target_archive_path,
                                                        flow=ArchiveFlow.PUT,
                                                        data=data_to_archive)
        
        print("Done!!!")
    
    def archive_investigation_computations(self, dest_folder:Path, flow:ArchiveFlow, data:ArchiveInvestigatorData)->Union[None, ArchiveInvestigatorData]:
        
        backup_investigator:ArchiveInvestigator = ArchiveInvestigator(working_folder=dest_folder,
                                             var_name=str(data.design_info.design_info.DesignID),
                                             use_db=True)
        
           
        if flow == ArchiveFlow.PUT:
            backup_investigator.design_info.design_info = data.design_info.design_info
            backup_investigator.design_info.wetlab_data = data.design_info.wetlab_results
            
            backup_investigator.investigator.comparison_eval_result.ratios = data.investigator.investigator_results.comparison_eval_results.ratios
            # print(type(data.investigator.investigator_results.comparison_eval_results.BRaise_list[0]))
            # input()
            backup_investigator.investigator.comparison_eval_result.BRaise_list = data.investigator.investigator_results.comparison_eval_results.BRaise_list
            backup_investigator.investigator.comparison_eval_result.BUratio_list = data.investigator.investigator_results.comparison_eval_results.BUratio_list
            backup_investigator.investigator.comparison_eval_result.bound_total_list = data.investigator.investigator_results.comparison_eval_results.bound_total_list
            backup_investigator.investigator.comparison_eval_result.unbound_total_list = data.investigator.investigator_results.comparison_eval_results.unbound_total_list
            backup_investigator.investigator.comparison_eval_result.nuc_penatly_count = data.investigator.investigator_results.comparison_eval_results.nuc_penatly_count
            backup_investigator.investigator.comparison_eval_result.first_BUratio = data.investigator.investigator_results.comparison_eval_results.first_BUratio
                       
            backup_investigator.investigator.lmv_values.lmv_comps = data.investigator.investigator_results.lmv_values.lmv_comps
            
            backup_investigator.investigator.lmv_assertions.comp_compare_to_mfe = data.investigator.investigator_results.lmv_assertions.comp_compare_to_mfe
            backup_investigator.investigator.lmv_assertions.unbouund_pronounced = data.investigator.investigator_results.lmv_assertions.unbouund_pronounced
            backup_investigator.investigator.lmv_assertions.bound_pronounced = data.investigator.investigator_results.lmv_assertions.bound_pronounced
            backup_investigator.investigator.lmv_assertions.is_on_off_switch = data.investigator.investigator_results.lmv_assertions.is_on_off_switch

            backup_investigator.investigator.comp_nuc_counts = data.investigator.investigator_results.comp_nuc_counts
            backup_investigator.investigator.num_groups = data.investigator.investigator_results.num_groups
            backup_investigator.investigator.total_structures_ensemble = data.investigator.investigator_results.total_structures_ensemble
            
            backup_investigator.investigator.lmv_references = data.investigator.lmv_references
            
            backup_investigator.scores.basic_scores = data.investigator.basic_scores
            backup_investigator.scores.advanced_scores = data.investigator.advanced_scores
            backup_investigator.scores.number_structures = data.investigator.number_structures
            
            
  
            return None
        
        elif flow == ArchiveFlow.GET:
            
            
            retreived_comparison_eval_results:ComparisonEvalResults = ComparisonEvalResults(ratios=backup_investigator.investigator.comparison_eval_result.ratios,
                                                                                            BRaise_list=backup_investigator.investigator.comparison_eval_result.BRaise_list,
                                                                                            BUratio_list=backup_investigator.investigator.comparison_eval_result.BUratio_list,
                                                                                            bound_total_list=backup_investigator.investigator.comparison_eval_result.bound_total_list,
                                                                                            unbound_total_list=backup_investigator.investigator.comparison_eval_result.unbound_total_list,
                                                                                            nuc_penatly_count=backup_investigator.investigator.comparison_eval_result.nuc_penatly_count,
                                                                                            first_BUratio=backup_investigator.investigator.comparison_eval_result.first_BUratio)
            
            retrieved_comparison_nucs:ComparisonNucResults = backup_investigator.investigator.comp_nuc_counts
            
            retrieved_lmv_values:ComparisonLMVResponse = ComparisonLMVResponse(lmv_comps=backup_investigator.investigator.lmv_values.lmv_comps)
            
            retreived_lmv_assertions:LMVAssertionResult = LMVAssertionResult(comp_compare_to_mfe=backup_investigator.investigator.lmv_assertions.comp_compare_to_mfe,
                                                                             unbouund_pronounced=backup_investigator.investigator.lmv_assertions.unbouund_pronounced,
                                                                             bound_pronounced=backup_investigator.investigator.lmv_assertions.bound_pronounced,
                                                                             is_on_off_switch=backup_investigator.investigator.lmv_assertions.is_on_off_switch)
            
            retreived_investigator:InvestigatorResults = InvestigatorResults(comparison_eval_results=retreived_comparison_eval_results,
                                                                             comp_nuc_counts=retrieved_comparison_nucs,
                                                                             lmv_values=retrieved_lmv_values,
                                                                             lmv_assertions=retreived_lmv_assertions,
                                                                             num_groups=backup_investigator.investigator.num_groups,
                                                                             total_structures_ensemble=backup_investigator.investigator.total_structures_ensemble)
            retrieved_reference_strucs:ReferenceStructures = backup_investigator.investigator.lmv_references
            
            retrieved_results:InvestigateEnsembleResults = InvestigateEnsembleResults(investigator_results=retreived_investigator,
                                                                                      basic_scores=backup_investigator.scores.basic_scores,
                                                                                      advanced_scores=backup_investigator.scores.advanced_scores,
                                                                                      number_structures=backup_investigator.scores.number_structures,
                                                                                      lmv_references=retrieved_reference_strucs)
            
            retrieved_design_info:DesignInformation = backup_investigator.design_info.design_info
            
            retrieved_archive: ArchiveInvestigatorData = ArchiveInvestigatorData(investigator=retrieved_results,
                                                                                 disign_info=retrieved_design_info)
            
            return retrieved_archive
                
                

    def retrieve_archive_and_run_analysis(self, save_folder:Path, run_name:str, pnas_dataset_path:Path, round:str, sublab:str, do_weighted:bool, is_agressive:bool = False, max_num_structs:int= 500000, archive_path:str=None):
        new_sara:Sara2API = Sara2API()
        puzzle_data: puzzleData
        pandas_sheet: DataFrame
        puzzle_data, pandas_sheet = new_sara.ProcessLab(path=pnas_dataset_path.as_posix(),
                                                      designRound_sheet=round,
                                                      sublab_name=sublab
                                                      )
        # pnas_data:List[Sara2StructureList] = []
        for design in puzzle_data.designsList:
            if os.path.isdir(archive_path) is False:
                raise FileExistsError(f'File {archive_path} is not a valid path')

            
            found_data:ArchiveData = self.archive_ensemble_data(dest_folder=archive_path,
                                        flow=ArchiveFlow.GET)
            
            if found_data.structs.num_structures < max_num_structs:
                struct_to_use:Sara2SecondaryStructure = found_data.fmn_folded_mfe
                if do_weighted is True:
                    struct_to_use = found_data.fmn_folded_weighted
                    
                reference_structures:EnsembleSwitchStateMFEStructs = EnsembleSwitchStateMFEStructs(switched_mfe_struct=struct_to_use,
                                                                                                non_switch_mfe_struct=found_data.structs.sara_stuctures[0])
                
                ensemble_groups: MultipleEnsembleGroups = self.nupack4.load_nupack_subopt_as_ensemble(span_structures=found_data.structs,
                                                                                                kcal_span_from_mfe=found_data.nupack_settings.kcal_span_from_mfe,
                                                                                                Kcal_unit_increments=found_data.nupack_settings.Kcal_unit_increments,
                                                                                                switch_state=reference_structures
                                                                                                )
                
                investigation_results:InvestigateEnsembleResults = self.scoring.investigate_and_score_ensemble(ensemble=ensemble_groups,
                                                                                                           is_aggressive=is_agressive)
                
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
                pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID, 'NumStructs'] =  found_data.structs.num_structures
            
            design_data_df:DataFrame = pandas_sheet.loc[pandas_sheet['DesignID']==design.design_info.DesignID]
            logging: PNASAnalysisLogging = PNASAnalysisLogging()
            if os.path.isdir(save_folder.as_posix()) == False:
                os.makedirs(save_folder.as_posix())
            logging.save_excel_sheet(design_data_df, save_folder, run_name)
            
    
    def archive_ensemble_data(self, dest_folder:Path, flow:ArchiveFlow, data: ArchiveData=None)->Union[None, ArchiveData]:
        
        backup_records:ArchiveSecondaryStructureList = ArchiveSecondaryStructureList(working_folder=dest_folder,
                                             var_name=str(data.design_info.design_info.DesignID),
                                             use_db=True)
        
           
        if flow == ArchiveFlow.PUT:
            print(data.design_info)
            backup_records.data.design_info = data.design_info
            backup_records.data.nupack_settings = data.nupack_settings
            backup_records.data.structs = data.structs
            backup_records.data.fmn_folded_mfe = data.fmn_folded_mfe
            backup_records.data.fmn_folded_weighted = data.fmn_folded_weighted
          
            return None
        elif flow == ArchiveFlow.GET:
            retrieved_archive:ArchiveData = ArchiveData(design_info=backup_records.data.design_info,
                                                        nupack_settings=backup_records.data.nupack_settings,
                                                        structs=backup_records.data.structs,
                                                        fmn_folded_mfe=backup_records.data.fmn_folded_mfe,
                                                        fmn_folded_weighted=backup_records.data.fmn_folded_weighted)
            return retrieved_archive
                                 

    
    # def switchyness_analysis(self, designs_structures:List[Sara2StructureList]):
    #     pass
    

def get_nupack_ensemble_structs_and_archive_r101():
    parser = argparse.ArgumentParser(description='Get and process R101 PNAS data')
    
    parser.add_argument('--pnas', 
                        type=Path,
                        required=True,
                        help='Path to the pnas file for analysis')
    
    parser.add_argument('--sublab', 
                        type=str,
                        default='',
                        required=True,
                        help='sublab in R101 to generate and archive enemble data for')
    
    parser.add_argument('--round', 
                        type=str,
                        default='Round 7 (R101)',
                        required=False,
                        help='Round to run')
    
    parser.add_argument('--archive', 
                        type=Path,
                        required=True,
                        help='Path to the data nut squirrel archive folder')
    
    parser.add_argument('--material',
                        type=str,
                        choices=[MaterialParameter.rna06_nupack4.name,
                                 MaterialParameter.rna95_nupack4.name,
                                 MaterialParameter.rna99_nupack3.name,
                                 MaterialParameter.rna95_nupack3.name],
                        required=True,
                        help='Material parameters that you can use'
                        )
    
    parser.add_argument('--temp',
                        type=int,
                        required=True,
                        help='Temperature to perform the fold at in degrees C'
                        )
    
    parser.add_argument('--span',
                        type=int,
                        required=True,
                        help='Span from MFE to grab ensemble data for.'
                        )
    
    parser.add_argument('--unit',
                        type=float,
                        required=True,
                        help='Increment size in Kcal of each unit of the ensemble when it is divided up'
                        )
    
    
    args = parser.parse_args()
    
    material: MaterialParameter = MaterialParameter[args.material]
    
    nupack_settings:NupackSettings = NupackSettings(material_param=material,
                                                    temp_C=args.temp,
                                                    kcal_span_from_mfe=args.span,
                                                    Kcal_unit_increments=args.unit,
                                                    sequence='')
    
    # parser.add_argument('--r',
    #                     '--round', 
    #                     type=str,
    #                     default='Round 7 (R101)',
    #                     help='Name to use for the output f')
    
    
    details:str= 'all'#f'20k_filtered_weighted_100K_gtrequal2_nucpenalty_run_1ish'

    
    process_pnas:ProcessPNAS = ProcessPNAS()
    process_pnas.record_nupack_ensemble_structs(pnas_dataset_path=args.pnas,
                                                round=args.round,
                                                sublab=args.sublab,
                                                nupack_settings=nupack_settings,
                                                archive_path=args.archive)
    
    # #Round 7 (R101)
    # same_state:str='3'
    # sublab_name:str = "good"#f'Same State NG {same_state}'
    # is_aggressive:bool = False
    # save_title:str = sublab_name + "_open"
    # run_name:str = "data_nut_test"#f'SSNG{same_state}_{details}'


    # pnas_path:str = '/home/rnauser/test_data/pnas_testing_tweak.xlsx'
    # timestr = time.strftime("%Y%m%d-%H%M%S")
    # save_path:str = f'/home/rnauser/test_data/run_data/{run_name}/pnas_eternacon_{timestr}.xlsx'
        

def perform_serena_ensemble_computations():
    parser = argparse.ArgumentParser(description='Get and process R101 PNAS data')
    
    parser.add_argument('--pnas', 
                        type=Path,
                        required=True,
                        help='Path to the pnas file for analysis')
    
    parser.add_argument('--round', 
                        type=str,
                        default='Round 7 (R101)',
                        required=False,
                        help='Round to run')
    
    parser.add_argument('--sublab', 
                        type=str,
                        default='',
                        required=True,
                        help='sublab in R101 to generate and archive enemble data for')
    
    
    parser.add_argument('--source', 
                        type=Path,
                        required=True,
                        help='Path to the data nut squirrel archive folder')
    
    parser.add_argument('--target', 
                        type=Path,
                        required=True,
                        help='Path to save the new data nut squirrel archive folder')
    
    parser.add_argument('--do-weighted', 
                        action="store_true",
                        required=False,
                        help='Use weighted struct from Vienna2 enemble')
    
    parser.add_argument('--do-agressive', 
                        action="store_true",
                        required=False,
                        help='use agressive limits')
    
    parser.add_argument('--max-structs', 
                        type=int,
                        default=500000,
                        required=False,
                        help='sublab in R101 to generate and archive enemble data for')
    
    
    args = parser.parse_args()
    

    process_pnas:ProcessPNAS = ProcessPNAS()
    
    # print(args.do_agressive)
    
    process_pnas.perform_investigation_computations(pnas_dataset_path=args.pnas,
                                                    sublab=args.sublab,
                                                    source_archive_path=args.source,
                                                    target_archive_path=args.target,
                                                    do_weighted=args.do_weighted,
                                                    max_num_structs=args.max_structs,
                                                    is_agressive=args.do_agressive,
                                                    round=args.round)
    

def switchyness_analysis():
    pass
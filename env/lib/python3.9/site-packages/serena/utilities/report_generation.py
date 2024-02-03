"""
Class for generating a report file for each analysis type
"""
from pathlib import Path
import time
from typing import List, Dict
from dataclasses import dataclass
from datetime import datetime
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
import matplotlib

from serena.utilities.ensemble_variation import EnsembleVariation, EVResult

@dataclass
class SequenceResults():
    result_data: EVResult
    result_name: str
    
@dataclass
class SequenceInfo():
    lab_name:str = ''
    sequence: str = ''
    sequence_name: str = ''
    sequence_ID: int = -1
    folded_energy: float = 0
    ligand_oligo_name:str = ''
    eterna_score:float = 0
    fold_change:float = 0
    number_of_clusters:int = 0
    temp_C: int = 0
    span = 0
    units = 0
    folded_struct:str


class FullRunInfo():

    def __init__(self, run_name:str, run_ID: int, sequence_info: SequenceInfo) -> None:
        self._run_results: List[SequenceResults] = []
        self._sequence_info: SequenceInfo  = sequence_info
        self._run_name: str = run_name
        self._run_ID: int = run_ID
    
    def add_sequence_result(self, result: SequenceResults):
        self._run_results.append(result)
    
    @property
    def run_results(self):
        return self._run_results
    
    @run_results.setter
    def run_results(self, results:List[SequenceResults] ):
        self._run_results = results

    @property
    def sequence_info(self):
        return self._sequence_info
    
    @sequence_info.setter
    def sequence_info(self, info:SequenceInfo):
        self._sequence_info = info
    
    @property
    def run_name(self):
        return self._run_name
    
    @run_name.setter
    def run_name(self, name:str):
        self._run_name = name

    @property
    def run_ID(self):
        return self._run_ID
    
    @run_ID.setter
    def run_ID(self, ID:int):
        self._run_ID = ID

@dataclass
class PrepedReportData():
    new_list_string_mfe: List[float]
    new_list_string_rel: List[float]
    new_switch_string_folded: List[float]
    time_span: List[float]
    #tick_span = []
    mfe_value:float 
    #seed_value:float
    #tick_value:float
    #units:int


class LocalMinimaVariationReport():
    """
    Class for making the local minima variation (ensemble variation) report
    """
    
    def __init__(self, working_folder: Path) -> None:
        self._working_folder: Path = working_folder
    
    @property
    def working_folder(self):
        return self._working_folder
    
    @working_folder.setter
    def working_folder(self, folder_path:Path):
        self._working_folder = folder_path

    

    def prep_data_for_report(self, run_info: FullRunInfo):
        #now save data to csv file        
        timestr = time.strftime("%Y%m%d-%H%M%S")

        ev_result_mfe: EVResult = run_info.run_results[0].result_data
        ev_result_rel: EVResult = run_info.run_results[1].result_data
        switch_result_folded: EVResult = run_info.run_results[2].result_data

        time_span: List[float] = []
        tick_span = []
        index_span = range(len(ev_result_mfe.groups_list))
        

        mfe_value:float = ev_result_mfe.groups_list[0].mfe_freeEnergy
        seed_value:float = mfe_value
        tick_value:float = 0
        
        units:int = run_info.sequence_info.units

        num_samples: int = len(ev_result_mfe.groups_list)
        for index in range(num_samples):
            seed_value = seed_value + float(units)
            tick_value = tick_value + float(units)
            #time_span is teh MFE values
            #tick_span is the index value (i.e. 0, 0.5, 1, 1.5)
            time_span.append(seed_value)
            tick_span.append(tick_value)

        

        new_list_string_mfe: List[float] = []
        for ev in ev_result_mfe.group_ev_list:
            ev_value = ev.ev_normalized
            new_list_string_mfe.append(ev_value)

        new_list_string_rel: List[float] = []
        for ev in ev_result_rel.group_ev_list:
            ev_value = ev.ev_normalized
            new_list_string_rel.append(ev_value)

        new_switch_string_folded: List[float] = []
        for ev in switch_result_folded.group_ev_list:
            ev_value = ev.ev_normalized
            new_switch_string_folded.append(ev_value)

    

    def generate_text_report(self, run_info: FullRunInfo):
        #now save data to csv file        
        timestr = time.strftime("%Y%m%d-%H%M%S")

        ev_result_mfe: EVResult = run_info.run_results[0].result_data
        ev_result_rel: EVResult = run_info.run_results[1].result_data
        switch_result_folded: EVResult = run_info.run_results[2].result_data

        time_span: List[float] = []
        tick_span = []
        index_span = range(len(ev_result_mfe.groups_list))
        

        mfe_value:float = ev_result_mfe.groups_list[0].mfe_freeEnergy
        seed_value:float = mfe_value
        tick_value:float = 0
        
        units:int = run_info.sequence_info.units

        num_samples: int = len(ev_result_mfe.groups_list)
        for index in range(num_samples):
            seed_value = seed_value + float(units)
            tick_value = tick_value + float(units)
            #time_span is teh MFE values
            #tick_span is the index value (i.e. 0, 0.5, 1, 1.5)
            time_span.append(seed_value)
            tick_span.append(tick_value)

        

        new_list_string_mfe: List[float] = []
        for ev in ev_result_mfe.group_ev_list:
            ev_value = ev.ev_normalized
            new_list_string_mfe.append(ev_value)

        new_list_string_rel: List[float] = []
        for ev in ev_result_rel.group_ev_list:
            ev_value = ev.ev_normalized
            new_list_string_rel.append(ev_value)

        new_switch_string_folded: List[float] = []
        for ev in switch_result_folded.group_ev_list:
            ev_value = ev.ev_normalized
            new_switch_string_folded.append(ev_value)



        csv_log_results: List[str]=[]
        csv_log_results.append("Kcal,LMSV_U_mfe,LMSV_U_rel,LMSV_US_target,LMSV_US_folded\n")
        for index in range(len(new_list_string_mfe)):
            kcal = time_span[index]
            LMSV_U_mfe = new_list_string_mfe[index]
            LMSV_U_rel = new_list_string_rel[index]
            LMSV_US_folded = new_switch_string_folded[index]
            line:str = f'{kcal},{LMSV_U_mfe},{LMSV_U_rel},{LMSV_US_folded}\n'
            csv_log_results.append(line)


        csv_record_pathstr = f'{self._working_folder}/{run_info.run_name}_{timestr}.csv'
        csv_lines:List[str]=[]
        with open(csv_record_pathstr, 'w') as csv_file:
            #first write teh header
            csv_lines.append(f'Local Minima Structure Variation Data\n')
            csv_lines.append(f'Creation Date={datetime.now()}\n')
            csv_lines.append("---------------------------------------\n")
            csv_lines.append("***DESIGN INFO***\n")
            csv_lines.append(f'Design Name = {run_info.sequence_info.sequence_name}\n')
            csv_lines.append(f'DesignID = {run_info.sequence_info.sequence_ID}\n')
            csv_lines.append(f'Lab Name = {run_info.sequence_info.lab_name}\n')
            csv_lines.append(f'Sequence = {run_info.sequence_info.sequence}\n')
            csv_lines.append(f'Eterna_Score = {run_info.sequence_info.eterna_score}\n')
            csv_lines.append(f'FoldChange = {run_info.sequence_info.fold_change}\n')
            csv_lines.append(f'2nd State Folded Structure = {run_info.sequence_info.folded_struct}\n')
            csv_lines.append(f'2nd State Folded Oligo Energy = {run_info.sequence_info.folded_energy}\n')
            csv_lines.append(f'Energy Span from MFE = {run_info.sequence_info.span}\n')
            csv_lines.append(f'Energy span units = {units}\n')
            csv_lines.append("---------------------------------------\n")
            csv_lines.append("***RAW DATA***\n")
            csv_lines = csv_lines + csv_log_results
            csv_lines.append("---------------------------------------\n")
            csv_lines.append("EOF\n")
            csv_file.writelines(csv_lines)
    
    def generate_plot(self, run_info: FullRunInfo, report_data:PrepedReportData, save_folder_name:Path):
        #now save teh data

        timestr = time.strftime("%Y%m%d-%H%M%S")
        fig, ax = plt.subplots()
        info_str:str = f'Eterna_Score = {run_info.sequence_info.eterna_score}, FoldChange = {run_info.sequence_info.fold_change}, Num Clusters = {run_info.sequence_info.number_of_clusters}'
        plt.suptitle(f'LMV Switch plot for {run_info.sequence_info.sequence_name}\nEterna Lab = {run_info.sequence_info.lab_name}\nDesign ID = {run_info.sequence_info.sequence_ID}\n',fontsize=12)
        plt.title(info_str, fontsize=10)
        #fig = plt.figure()
        
        #ax.set_xticks(tick_span)
        plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}')) # 2 decimal places
        plt.plot(report_data.time_span, report_data.new_list_string_mfe, 'b^-', label='LMV_U mfe')
        plt.plot(report_data.time_span, report_data.new_list_string_rel, 'ro-', label='LMV_U rel')
        plt.plot(report_data.time_span, report_data.new_switch_string_folded, 'gs-', label='LMV_US folded')
        #y_ticks = [0,5,10,15,20,25,30,35,40,45,50]
        y_ticks = np.arange(-10,65, step=5)
        plt.xticks(report_data.time_span)
        plt.yticks(y_ticks)
        #plt.yticks()
        plt.grid(True)   
        plt.legend(loc='lower right',fontsize="x-small")
        plt.subplots_adjust(top=.8, bottom=.2, left=.12, right=.95)  
        plt.tick_params(axis='x',labelrotation=90)  
        plt.ylabel("Local Minima Structure Variation (LMSV)")
        plt.xlabel("Local Kcal Energy along Ensemble")
        #plt.figtext(0.54, 0.01, delta_message, ha="center", fontsize=10, bbox={"facecolor":"orange", "alpha":.5, "pad":2})
        trans = ax.get_xaxis_transform()
        delat_energy:float = 2
        lower_range_folded_energy: float = run_info.sequence_info.folded_energy - (delat_energy/2)
        uper_range_folded_energy: float = run_info.sequence_info.folded_energy + (delat_energy/2)
        plt.axvline(x=run_info.sequence_info.folded_energy, color="green", linestyle="--")
        plt.axvline(x=lower_range_folded_energy, color="green", linestyle=":")
        plt.axvline(x=uper_range_folded_energy, color="green", linestyle=":")
        plt.text(lower_range_folded_energy, .06, '   2nd State', transform=trans, fontsize=7)
        plt.text(lower_range_folded_energy, .01, ' 2Kcal range', transform=trans, fontsize=7)
        plt.text(run_info.sequence_info.folded_energy, .06, '    2nd State', transform=trans, fontsize=7)
        plt.text(run_info.sequence_info.folded_energy, .01, '  folded Energy', transform=trans, fontsize=7)

        ev_mfe_lower:float = report_data.time_span[0]
        ev_mfe_upper:float = report_data.mfe_value + delat_energy
        
        plt.axvline(x=ev_mfe_lower, color="blue", linestyle=":")
        plt.axvline(x=ev_mfe_upper, color="blue", linestyle=":")
        plt.text(ev_mfe_lower, .06, '   1st State', transform=trans, fontsize=7)
        plt.text(ev_mfe_lower, .01, ' 2Kcal range', transform=trans, fontsize=7)
        #ax.set_ybound(lower=0, upper=70)
        

        file_name:str = f'{run_info.sequence_info.sequence_name}_{run_info.sequence_info.sequence_ID}'
        
        plt.savefig(f'{save_folder_name.as_posix()}/{file_name}_{timestr}.png')

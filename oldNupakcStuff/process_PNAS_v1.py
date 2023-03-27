#!/usr/bin/env python
# coding: utf-8

# In[25]:


# Import NUPACK Python module
from nupack import *
import NUPACK4_module1_V5 as nup4
import pandas as pd
import sys
import openpyxl


rna06_nupack4_parameter=["rna06", "Based on [Mathews99] and [Lu06] with additional parameters [Xia98,Zuker03] including coaxial stacking [Mathews99,Turner10] and dangle stacking [Serra95,Zuker03,Turner10] in 1M Na+."]
rna95_nupack4_parameters=["rna95", "Based on [Serra95] with additional parameters [Zuker03] including coaxial stacking [Mathews99,Turner10] and dangle stacking [Serra95,Zuker03,Turner10] in 1M Na+."]
rna99_nupack3_parameters=["rna99-nupack3","Parameters from [Mathews99] with terminal mismatch free energies in exterior loops and multiloops replaced by two dangle stacking free energies. Parameters are provided only for 37 ∘C."]
rna95_nupack3_parameters=["rna95-nupack3", "Same as rna95 except that terminal mismatch free energies in exterior loops and multiloops are replaced by two dangle stacking free energies."]
rna_parameterSet=[rna06_nupack4_parameter, rna95_nupack4_parameters, rna99_nupack3_parameters, rna95_nupack3_parameters]
#use this param
paramNum=0

deltasToInspect_list = [1,2,3,4]

temperature = 37

pnasPath = r'pnas.2112979119.sd01.xlsx'
designRound = "Round 7 (R101)"
roundData = ['DesignID', 'Design', 'Player', 'Puzzle_Name', 'Eterna_Score', 'FoldChange', 'Sequence']
sheet = pd.read_excel(pnasPath, sheet_name=designRound)
dataframe = sheet[roundData]


print("opened file and processing")

design_header_defaults=['Index','DesignID','Design','Sequence','Player','Puzzle_Name','Eterna_Score','FoldChange','mfe_freeEnergy','mfe_stackEnergy','ensemble_defect']
ev_header=[]
for deltaIndex in range(len(deltasToInspect_list)):
    delta = str(deltasToInspect_list[deltaIndex])
    ev_header.append("EV_classic_"+delta)
    ev_header.append("EV_classic_Norm_"+delta)
    ev_header.append("EV_new_"+delta)
desing_header = design_header_defaults+ev_header

wb = openpyxl.Workbook()
ws_write = wb.create_sheet(0)
ws_write.append(desing_header)

index=1
for lab_design_index in dataframe.index:
    sequence = dataframe['Sequence'][lab_design_index]
    DesignID = dataframe['DesignID'][lab_design_index]
    DesignName = dataframe['Design'][lab_design_index]
    Player = dataframe['Player'][lab_design_index]
    Eterna_Score=dataframe['Eterna_Score'][lab_design_index]
    FoldChange=dataframe['FoldChange'][lab_design_index]
    Puzzle_Name=dataframe['Puzzle_Name'][lab_design_index]
    
    print("Running design "+ str(DesignID))
    mfe_structure, mfe_freeEnergy, mfe_stackEnergy, ensemble_defect, EV_Delta_dict  = nup4.processSequence_Nupack4(sequence, temperature, paramNum, deltasToInspect_list)
    print("Completed running "+ str(DesignID))
    
    desing_info_list_defauls = [index,DesignID,DesignName,sequence,Player,Puzzle_Name,Eterna_Score,FoldChange,mfe_freeEnergy,mfe_stackEnergy,ensemble_defect]
    desing_info_EV_list=[]
    for deltaIndex in range(len(deltasToInspect_list)):
        delta = deltasToInspect_list[deltaIndex]
        desing_info_EV_list.append(EV_Delta_dict[delta]["EV_old_raw"])
        desing_info_EV_list.append(EV_Delta_dict[delta]["EV_old_norm"])
        desing_info_EV_list.append(EV_Delta_dict[delta]["EV_new"])
    desing_info_List = desing_info_list_defauls + desing_info_EV_list
    ws_write.append(desing_info_List)
    print("Wrote data to excel for "+ str(DesignID))
    
   
import time
timestr = time.strftime("%Y%m%d-%H%M%S")
filename_withtime="PNAS_anlaysis{0}.xlsx".format(timestr)
wb.save(filename=filename_withtime)   
print("its done. did it work")


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





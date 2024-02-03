"""
pandas baesd log reporting
"""

from typing import List
import pandas as pd
from pandas import DataFrame
import os

class PNASAnalysisLogging():

    def __init__(self) -> None:
        pass

    def open_sublab_from_excel(self,path:str, sheet_name:str, sublab:str):
        sheet:DataFrame = pd.DataFrame()
        
        try:            
            sheet:DataFrame = pd.read_excel(path, sheet_name=sheet_name)
        except:
            raise FileExistsError(f'Unable to open {path} to grab sublab info')
        
        sublab_sheet:DataFrame = sheet
        if sublab != "":
            sublab_sheet:DataFrame  = sheet[sheet['Puzzle_Name'] == sublab]
        return sublab_sheet
    
    def save_excel_sheet(self, df:DataFrame, excel_path:str, sheet_name:str):
        if not os.path.exists(excel_path):
            df.to_excel(excel_path, sheet_name=sheet_name, index=False)
        else:
            with pd.ExcelWriter(excel_path, engine='openpyxl', if_sheet_exists='overlay', mode='a') as writer:
                df.to_excel(writer, sheet_name=sheet_name, startrow=writer.sheets[sheet_name].max_row, header=None, index=False)

    def save_dataframe_to_excel(self, src_dataframe:DataFrame, dst_path:str, dst_sheet_name:str):
        with pd.ExcelWriter(dst_path) as writer:
            src_dataframe.to_excel(writer, sheet_name=dst_sheet_name)

    def add_new_named_collumn(self,src_datafrm:DataFrame, column_name:str, column_list:List[str]):
        src_datafrm[column_name]=column_list
        return src_datafrm


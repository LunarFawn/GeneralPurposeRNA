# create the class for holding all the data associated with each design. 
# use namedtuple for memory for for what is parsed from spreedsheets and 
# fold data generated by nupack as should not change
from decimal import Decimal
import string
from tokenize import Double
from turtle import st
from typing import Dict, NamedTuple, List, Optional
from unicodedata import decimal
import pandas as pd
import sys
import openpyxl
from dataclasses import dataclass
import mysql.connector
from mysql.connector import Error
import pandas as pd
from collections import OrderedDict
from enum import Enum

import currentStuff.nupackAPI_Sara2_Ver1 as nupackAPI

#this is the stuff that is generated by Rhijus wet lab
@dataclass
class WetlabData:
    Sequence: str = ""
    Eterna_Score: float = -1
    Baseline_Subscore: float = -1
    Folding_Subscore: float = -1
    Switch_Subscore: float = -1
    NumberOfClusters1: int = -1
    FoldChange: float = -1
    FoldChange_err_factor: float = -1
    KDOFF: float = -1
    KDON: float = -1
    ddG: float = -1
    ddG_err: float = -1

#data commonly associated with secondary structures. written with alternate structures in mind
#class structureData(NamedTuple):
#    secondaryStructure: str
#    freeEnergy: float
#    deltaMfeEnergy: float  

#this is the data the generated by folding software and in this case it is nupack
class NupackFoldData(NamedTuple):
    temperature: int
    doPknot: bool
    isPknot: bool
    #mfeInfo: structureData
    #alternateStructureList: List[structureData]
    pairprobsList: List[List[float]]
    #primaryPairsList:list
    #primaryPairsSortedList: list
    #secondaryPairsList: list
    #secondaryPairsSortedList: list

#class Oligos(NamedTuple):
#    oligoSequence: str
#    oligoSequenceLenght:str
#    OligoConcentration: str

class DesignInformation(NamedTuple):
    Sequence: str
    Sequence_Length: int
    #oligosList: List[Oligos]
    DesignID: int
    Design: str
    Player: str
    Puzzle_Name: str

#entry point for each design in a puzzle/lab
class DesignPerformanceData(NamedTuple):    
    DesignInfo: DesignInformation
    wetlabResults: WetlabData
    nupackFoldResults: NupackFoldData

#entry point to puzzle or lab as we call them
class puzzleData(NamedTuple):
    Puzzle_Name: str
    designsList: List[DesignPerformanceData]
    designsDict: Dict[str, DesignPerformanceData]


#class RoundData(NamedTuple):
#    roundName: str
#    roundID: int
#    puzzleLabList: List[puzzleData]

       

# for each design entry we make a DesignData object and add that to a Lst once that is done 
# we add it to a puzzleData object and once that is done for all the puzzles we then add it to a round data and wea re finished.
# I think round data should be changable

def GetNamesClass(className):
    variables = [i for i in vars(className).keys() if not callable(i) and not i.startswith('__') ]
    #variables = list(vars(className).keys())
    return variables

def openExcelWetlab(path, designRoundSheet):
    #$designRound = "Round 7 (R101)"
    pnasPath = r'pnas.2112979119.sd01.xlsx'
    designRound = "Round 7 (R101)"
    roundData = ['DesignID', 'Design', 'Player', 'Puzzle_Name', 'Eterna_Score', 'FoldChange', 'Sequence']
    sheet = pd.read_excel(path, sheet_name=designRoundSheet).itertuples()
    return sheet

def GenerateNupackEntry(wetlabDataObject: WetlabData, temperature, doPknot):
    #now need to run each design through nupack
    #pairProbsList = List[List[float]]
    pairProbsList = nupackAPI.GetPairProbs2DArray(wetlabDataObject.Sequence, 'rna95-nupack3', temperature)
    
    isPknot=False
    if doPknot is False:
        isPknot=False
    else:
        doPknot=True
    
    nupackEntry = NupackFoldData(temperature=temperature, doPknot=doPknot, isPknot=isPknot, pairprobsList=pairProbsList)
    return nupackEntry

def GetSetValue(getObject, variableSeek, setObject):
    value = getattr(getObject, variableSeek)
    setattr(setObject, variableSeek, value)
    return setObject

def GenerateWetlabEntry(row: NamedTuple):
    wetlab = WetlabData
    #this is a list of the variable names from the wetlab class
    variables = GetNamesClass(WetlabData)
    for variable in variables:
        #sourceValue = getattr(row, variable)
        #setattr(wetlab, variable, sourceValue)
        wetlab =  GetSetValue(row, variable, wetlab)
    #for name, value in WetlabData.__dict__().items():
    #    wetlab =  GetSetValue(row, name, wetlab)
    #all the contents of the wetlab excel should now be loaded into wetlabData object wetlab
    #return it then
    return wetlab


def GenerateDesignInfo(row: NamedTuple):
    designInfo = DesignInformation
    #this is a list of the variable names from the wetlab class
    members = GetNamesClass(DesignInformation)
    for variable in members:
        GetSetValue(row, variable, designInfo)
    return designInfo
    

def ProcessLab(path, designRound_sheet):
    sheet = openExcelWetlab(path, designRound_sheet)
    #first do the Design entry stuff  

    designs: List[DesignPerformanceData] = []
    designsDict: Dict[str, DesignInformation] = {}
     
    for row in sheet:
        wetlabResults: WetlabData
        nupackRestuls: NupackFoldData
        desingInfo: DesignInformation
        # this is a single line from teh file and representas a single design. load this into wetlabdata
        # and then do nupack. each row is in a namedtouple formate
        wetlabResults = GenerateWetlabEntry(row)
        nupackRestuls = GenerateNupackEntry(wetlabResults, 37, False)
        desingInfo = GenerateDesignInfo(row)
        DesingData = DesignPerformanceData(DesignInfo=desingInfo, wetlabResults=wetlabResults, nupackFoldResults=nupackRestuls)
        designs.append(DesingData)
        designsDict[DesingData.DesignInfo.DesignID]=DesingData
    puzzlename = designs[0].DesignInfo.Puzzle_Name
    #lets stop at puzzle data until this is fully gigure out and tested a bit
    puzzleInfo = puzzleData(Puzzle_Name=puzzlename, designsList=designs, designsDict=designsDict)    
    return puzzleInfo



 #this will give a list of pairs based on prob level
def NucProbSearch(self, foldData: NupackFoldData, probThresh):
    pairsList = foldData.pairprobsList
    snuppPairsDict= dict(string, float)
    snuppList = []

    for i in range(len(pairsList)):        
        for j in range(len(pairsList[i])):
            pairValue = pairsList[i][j]
            if pairValue >= probThresh:
                #pair = NucPair(i,j)
                pairName = "{i}:{j}".format(i=i, j=j)
                snuppPairsDict[pairName]=pairValue
                snuppList.append(pairName)
    return snuppPairsDict, snuppList     
#this is new stuff to do snupp pairs like it should have been in the C# code.

#need a collection of nupackfold data

#need a dict 
#pairsDict= dict(string, list)






#find if a list of elements are in a molecule based on in silico calulations

#make rainbow of pairs in stuct


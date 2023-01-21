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
import numpy as np

import currentStuff.nupackAPI_Sara2_Ver1 as nupackAPI
import Sara2_dev.Sara2_API_Python3_V1 as sara2Root

from pathlib import Path

from draw_rna import draw, draw_all, draw_from_rdat, draw_utils




class SnuppPair(object):

    def __init__(self, nuc1: int, nuc2: int) -> None:
        self._i = nuc1
        self._j = nuc2
        self._pair = f'{self._i}:{self._j}'
    
    def __str__(self) -> str:
        return self._pair
    
    def __repr__(self):
        return f'SnuppPair(pair={self._pair}, i={self._i}, j={self._j})'

    @property
    def snuppPair(self):
        return self._pair

    @property
    def i(self):
        return self._i

    @property
    def j(self):
        return self._j

#returns teh ranges and designs in each group for a rainbow color map
#this is the entrance
class NucPair(object):

    def __init__(self, nuc1: int, nuc2: int, pairProb:float, designID_str: str, _foldchange: Optional[float]=None)-> None:
        #convert to comon i, j pair
        self._pair:SnuppPair = SnuppPair(nuc1, nuc2)
        self._probability=pairProb
        #self._designId = designID_str
        self._designIdList: List[str]
        self._designIdList.append(designID_str)
        self._foldChangeList: List[float]
        if _foldchange is not None:
            self._foldChangeList.append(_foldchange)
    
    def __str__(self) -> str:
        return self._pair
    
    def __repr__(self):
        return f'NucPair(pair={self._pair}, prob={self._probability}, deisgnIDList={self._designIdList})'

    @property
    def snuppPair(self):
        return self._pair

    @property
    def i(self):
        return self.snuppPair.i

    @property
    def j(self):
        return self.snuppPair.j

    @property
    def probability(self):
        return self._probability

    @property
    def designIDs(self):
        return self._designIdList
    
    @property
    def foldChangeList(self):
        return self._foldChangeList

    #@property
    #def designIDList(self):
    #    return self._designIdList
        
    def appendDesignId(self, value: List[str]):
        self._designIdList = self._designIdList + value
    
    def appendFoldChange(self, value: float):
        self._foldChangeList.append(value)       

class NucPairList(object):
    
    def __init__(self) -> None:
        #self._nucPairListFull : List[NucPair] = _nucs
        self._nucPairListFull : List[NucPair]
        self._sortedPairListFull: List[NucPair]
        self._snuppOnly: List[SnuppPair]
        self._pairsDict: OrderedDict[SnuppPair, NucPair]
        self._isEmpty:bool = True

        #self.refresh()

    def makeNewPairList(self, _nucs:List[NucPair], _foldchange: Optional[float]=None):
        self._nucPairListFull : List[NucPair] = _nucs
        self._isEmpty:bool = False
        if _foldchange is not None:
            pass
        self.refresh()

    @property
    def empty(self):
        return self._isEmpty

    def refresh(self):
        self.sort_probability()
        self.makeJustPairs()
        self.makePairsDict()

    def sort_probability(self):
        #sort list by probability
        sortedPairList: List[NucPair] = sorted(self._nucPairListFull, key=lambda x: x.probability)
        self._sortedPairListFull = sortedPairList

    def makeJustPairs(self):
        justPairs = []
        for pair in self._nucPairListFull:
            justPairs.append(pair.snuppPair)
        self._snuppOnly = justPairs

    def makePairsDict(self):
        pairsDict = OrderedDict[SnuppPair, NucPair]
        for pair in self.fullPairsList:
            pairsDict[pair.snuppPair] = pair
        self._pairsDict = pairsDict

    #this will add two lists of NucPairs and only keep the pairs that are common
    def appendNucPairList(self, value, _foldchange: Optional[float]=None) -> bool:
        if (self._isEmpty==False):
            #first need see what nuc pairs and probs are common
            secondList: NucPairList = NucPairList()
            if isinstance(value, NucPairList):
                secondList = value
            elif isinstance(value, List[NucPair]):
                secondList.makeNewPairList(value)    

            tempNucList: List[NucPair] = []
            for secondSnupp, secondNucPair in secondList.pairsDictFull.items():
                snuppPair: SnuppPair = secondSnupp.snuppPair
                if snuppPair in self.snuppOnlyList:
                    
                    #do fold change first so thta it is easier to handle probabilities dict
                    if  _foldchange is not None and _foldchange in secondList.pairsDictFull[snuppPair].foldChangeList:
                        #if they are the same then just need to append the design ID
                        self.pairsDictFull[snuppPair].appendFoldChange(_foldchange)

                    if secondList.pairsDictFull[snuppPair].probability == self.pairsDictFull[snuppPair].probability:
                        #if they are the same then just need to append the design ID
                        self.pairsDictFull[snuppPair].appendDesignId(secondNucPair.designIDs)
                        tempNucList.append(self.pairsDictFull[snuppPair])
                    else:
                        #if not the same pair prob then need to add both orignal and new to list as they are different
                        tempNucList.append(self.pairsDictFull[snuppPair])
                        tempNucList.append(secondNucPair)            
        else:
            newNucList: List[NucPair] 
            if isinstance(value, NucPairList):
                newNucList = value.fullPairsList
            elif isinstance(value, List[NucPair]):
                newNucList = value 
            self.makeNewPairList(newNucList)
   

    @property
    def pairsDictFull(self):
        return self._pairsDict

    @property
    def snuppOnlyList(self):
        return self._snuppOnly

    @property
    def fullPairsList(self):
        return self._sortedPairListFull

class SearchProtocol(Enum):
    FOLDCHANGE = 1
    PAIRPROB = 2
    NUCPAIR = 3

class SearchResult(object):
    #its really a dict

    def __init__(self, numSections: Optional[int] = 0) -> None:
        #currently only does FoldChange and PairProb
        self._resultsDict_float: Dict[float, NucPairList] = {}
        self._parameterValuesList: List[float] = None
        self._minParameterValue: float = None
        self._maxParameterValue: float = None
        self._rangeParameterValue: float = None
        self._numSections: int = numSections
        self._listSectionValues: List[float] = None
        #self.appendResult(searchParameter, searchResult)

    
    def appendResult(self, searchParameter: float, searchResult: NucPairList, numSections: Optional[int] = None):
        if searchParameter not in self._resultsDict_float:
            self._resultsDict_float[searchParameter] = NucPairList           
        else:
            primaryNucPairList: NucPairList = self._resultsDict_float[searchParameter]
            primaryNucPairList.addNucPairList(searchResult, searchParameter)
            self._resultsDict_float[searchParameter]=primaryNucPairList

        self.createParameterList()
        self.get_min_max()
        if(numSections is not None):
            self._numSections = numSections
        #only defne sections is neccessary
        if numSections > 0:
            self.defineSections()
        successString = f"Parameter {searchParameter} successfully added to SearchResult"
        return successString

    def defineSections(self):
        #first get the max and min
        numSections = self._numSections
        rangeOfValues = self._maxParameterValue - self._minParameterValue
        pointsPerSection:float = rangeOfValues/float(numSections)
        listSectionValues: List[float] =[]
        
        value = self._minParameterValue
        listSectionValues.append(value)
        #stop at the last one since that need to be the max to prevent weird math things from creating errors
        for num in range(numSections-1):            
            value = value + pointsPerSection
            listSectionValues.append(value)
        listSectionValues.append(self._maxParameterValue)
        self._listSectionValues = listSectionValues

    @property
    def sections(self):
        return self._numSections, self._listSectionValues

    @sections.setter
    def sections(self, numSections: int):
        self._numSections = numSections
        self.defineSections(numSections)   

    def createParameterList(self):
        self._parameterValuesList = self._resultsDict_float.keys()
    
    def get_min_max(self):
        self._minParameterValue = min(self._parameterValuesList)
        self._maxParameterValue = max(self._parameterValuesList)
        self._rangeParameterValue =  self._maxParameterValue - self._minParameterValue
    
    @property
    def parameterList(self):
        return self._parameterValuesList

    @property
    def keys(self):
        return self._parameterValuesList

    @property
    def min(self):
        return self._minParameterValue
    
    @property
    def max(self):
        return self._maxParameterValue
    
    @property
    def fullResults(self):
        return self._resultsDict_float

    @property
    def dictResults(self):
        return self._resultsDict_float

    @property
    def parameterResult(self, searchParameter:float):
        if (searchParameter not in self._resultsDict_float):
            return None
        else:     
            return self._resultsDict_float[searchParameter]

    @property
    def keySearch(self, searchKey:float):
        if (searchKey not in self._resultsDict_float):
            return None
        else:     
            return self._resultsDict_float[searchKey]
    
    @property
    def clear(self):
        self._resultsDict_float = {}
        self._parameterValuesList: List[float] = None
        self._minParameterValue: float = None
        self._maxParameterValue: float = None
        self._rangeParameterValue: float = None
        self._listSectionValues: List[float] = None
        return self._resultsDict_float



    


class NupackFoldDataEnum(Enum):
    PAIRPROBS = 1
    #NOT IMPLEMENTED YET
    PAIRSLIST =2

class GenerateRainbowStructurePlot:

    def __init__(self, _searchType: SearchProtocol, _sourceType: NupackFoldDataEnum, lab: sara2Root.puzzleData, minProbForPair: Optional[float] = 0) -> None:
        self.searchType: SearchProtocol = _searchType
        self.sourceType: NupackFoldDataEnum = _sourceType
        self.foldChangeDict: SearchResult = SearchResult(numSections=5)
        self.rawLabData: sara2Root.puzzleData = lab
        self._minProbForPair: float = minProbForPair
        self.initialize()
        #now a dictionary is loaded with all the fold changes and goodies

    def initialize(self):
        if self.searchType==SearchProtocol.FOLDCHANGE:
            if self.sourceType == NupackFoldDataEnum.PAIRPROBS: 
                #load fold change for lab into memory
                self.LoadLab()
        pass


    #this will give a list of pairs based on prob level
    #remember that nupackfolddata is for a single design for a single lab
    def PairsSearch_SingleDesign(self, foldData: sara2Root, probThresh:float, doInit:Optional[bool]=False, foldChange:Optional[float]=None):
        pairsList = foldData.NupackFoldData.pairprobsList
        designID = foldData.DesignInformation.DesignID
        
        tempSnuppPairsDict: Dict[NucPair, float] = {}  
        tempSnuppList: List[NucPair]= []
        #snuppPairsResults = self.SearchResultsGrouping

        for i in range(len(pairsList)):        
            for j in range(len(pairsList[i])):
                pairValue = pairsList[i][j]
                if pairValue >= probThresh:
                    pair = NucPair(i,j, pairValue, designID, foldChange)
                    if doInit==True:
                        #pairName = "{i}:{j}".format(i=i, j=j)
                        tempSnuppPairsDict[str(pair)]=pair
                        tempSnuppList.append(pair)
                    else:
                        #pairName = "{i}:{j}".format(i=i, j=j)
                        tempSnuppPairsDict[str(pair)]=pair
                        tempSnuppList.append(pair)
        return tempSnuppPairsDict, tempSnuppList     

    #this will load a datatype for furute searches to be faster
    #IMPORTANT need to remember to assume that currently teh data loaded needs to be 
    #ensured to be common to each fold change level. that means use NucPairList + override
    def LoadLab(self):
        #this will be orded Dict[fold_change, Dict[pairing_prob, List[NucPair]]]
        self.foldChangeDict.clear()
        if self.searchType==SearchProtocol.FOLDCHANGE:
            for design in self.rawLabData.designsList:
                nupackData = design.nupackFoldResults
                wetlabData = design.wetlabResults                
                currentFoldChange = wetlabData.FoldChange
                currentPairDict, currentPairsList = self.PairsSearch_SingleDesign(foldData=nupackData, probThresh=self._minProbForPair, doInit=True, foldChange=currentFoldChange)
                newNucPairList: NucPairList = NucPairList()
                newNucPairList.appendNucPairList(currentPairsList, currentFoldChange)
                self.foldChangeDict.appendResult(currentFoldChange, newNucPairList)    
            # now you should have a dictionary of fold changes with ech fold change having a list of nucpairs sorted by pairing prob
            # and the foldchange(s) the nuc pair is found in is(are) retained as well 

    #now need to write functions that filter through the loaded dic. the benifiet of the dict is that hopefully
    #we only need to do a few searches as they are all insorted order

    #need to write a function that is used as a base class for checking a value againsta a rna stucture vs finding commonalityin structures


    def new_foldChangeSearchPairs(self, foldChange_min: float, foldChange_max: float):
        #need to work off of fold change dict
        #self.foldChangeDict
        resultNucPairList: NucPairList = NucPairList()
        for foldChangeValue in self.foldChangeDict.keys:
            if (foldChangeValue >= foldChange_min and foldChangeValue <= foldChange_max):
                #then this fold change needs to be recorded
                listFromDict: NucPairList = self.foldChangeDict.parameterResult(foldChangeValue)
                resultNucPairList.appendNucPairList(listFromDict)
        return resultNucPairList

    def generate_Marker_Grouping_For_Analysis(self, minFoldChangePair:float):
        #first get the list of sections
        results: SearchResult = SearchResult()
        sectionNums: int
        sectionList: list[float]
        sectionNums, sectionList= self.foldChangeDict.sections

        for sectionIndex in range(len(sectionList)-1):
            lowerSectionValue = sectionList[sectionIndex]
            upperSectionValue = sectionList[sectionIndex+1]
            commonPairsList: NucPairList = self.new_foldChangeSearchPairs(lowerSectionValue, upperSectionValue)
            if commonPairsList.empty == False:
                results.appendResult(upperSectionValue, commonPairsList)
        return results


    def foldChangeSearch(self, roundInfo: sara2Root.puzzleData, foldChange_min: float, foldChange_max: float, probMin: float):
        
        commonPairsList:NucPairList = NucPairList()
        
        #i dont care about the disgn name or anything right now just the pairs
        commonFoldChangeList=[]
        for designData in roundInfo.designsList:
            #designData is everything to do with each individual design in the lab
            #get the foldchange
            if((designData.wetlabResults.FoldChange >=  foldChange_min) and (designData.wetlabResults.FoldChange<=foldChange_max) ):
                #get nuc pairs
                pairDict, pairList = self.PairsSearch_SingleDesign(designData.nupackFoldResults, probMin)
                commonPairsList.appendNucPairList(pairList)
                #now write to the list
        return commonPairsList
        


    # make this have a argument passed that define sget/set
    def FindPairs(self, roundInfo: sara2Root.puzzleData, wetLabElement: string, wetLabElement_min: float, wetLabElemt_max: float, probMin: float):
        #do a search 
        #initially do fold change, no. of clusters, and EternScores
        if wetLabElement is "foldchange":
            pairsList = self.foldChangeSearch(roundInfo, wetLabElement_min, wetLabElemt_max, probMin)
        return pairsList

   
   
   
   
   #old stuff and maybe some still used code... idk
   
    #this is agrouping and what defines a grouping is done by the function that calls this
    @dataclass
    class SearchResultsGrouping:
        _RawResults_Dict: Dict = {}
        _SortedResults_Dict: Dict = {}
        Prob_Pair: OrderedDict[float, List[NucPair]] = {}      


    def getRainbowColorMap(self):
        pass

    #find all pairs that have a parameter in common

    def writePair(self, i: int, j: int):
        pairName = "{i}:{j}".format(i=i, j=j)
        return pairName
    
  
    
    
    

    
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
#import mysql.connector
#from mysql.connector import Error
import pandas as pd
from collections import OrderedDict
from enum import Enum
import numpy as np
import copy

import nupackAPI_Sara2_Ver1 as nupackAPI
import Sara2_API_Python3_V1 as sara2Root
#from Sara2_API_Python3_V1 import Sara2 as sara2Root

from pathlib import Path

#from draw_rna import draw, draw_all, draw_from_rdat, draw_utils




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
    
    def __hash__(self):
        return hash((self._i, self._j, self._pair))

    def __eq__(self, other):
        return (self._i, self._j, self._pair) == (other._i, other._j, other._pair)

    def __ne__(self, other):
        # Not strictly necessary, but to avoid having both x==y and x!=y
        # True at the same time
        return not(self == other)

#returns teh ranges and designs in each group for a rainbow color map
#this is the entrance
class NucPair(object):
    round_value: int = 0
    def __init__(self, nuc1: int, nuc2: int, pairProb:float, designID_str: str, _foldchange: Optional[float]=None, fcl_round_digits:int = round_value)-> None:
        #convert to comon i, j pair
        self._pair:SnuppPair = SnuppPair(nuc1, nuc2)
        self._probability=pairProb
        #self._designId = designID_str
        self._designIdList: List[str] = []
        self._designIdList.append(designID_str)
        self._foldChangeList: List[float] = []
        self._foldChangeList_rounded: List[float] = []
        if _foldchange is not None:
            self._foldChangeList.append(_foldchange)
            self._foldChangeList_rounded.append(round(_foldchange, fcl_round_digits))
    
    def __hash__(self):
        return hash((self._pair, self._probability, self._designIdList, self._foldChangeList))

    def __eq__(self, other):
        return (self._pair, self._probability, self._designIdList, self._foldChangeList) == (other._pair, other._probability, other._designIdList, other._foldChangeList)

    def __ne__(self, other):
        # Not strictly necessary, but to avoid having both x==y and x!=y
        # True at the same time
        return not(self == other)  

    def __str__(self) -> str:
        return self._pair._pair
    
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
    
    @property
    def foldChangeList_rounded(self):
        return self._foldChangeList_rounded
    

    #@property
    #def designIDList(self):
    #    return self._designIdList
        
    def appendDesignId(self, value: List[str]):
        self._designIdList = self._designIdList + value
    
    def appendFoldChange(self, value: float, fcl_round_digits:int=round_value):
        self._foldChangeList.append(value)
        self._foldChangeList_rounded.append(round(value, fcl_round_digits))

     

class NucPairList(object):
    
    def __init__(self) -> None:
        #self._nucPairListFull : List[NucPair] = _nucs
        self._nucPairListFull : List[NucPair]
        self._sortedPairListFull: List[NucPair]
        self._snuppOnly: List[SnuppPair]
        self._snuppOnly_str: List[str]
        self._pairsDict: OrderedDict[SnuppPair, NucPair]
        self._isEmpty:bool = True

        #self.refresh()

    def makeNewPairList(self, _nucs:List[NucPair], _foldchange: Optional[float]=None):
        self._nucPairListFull : List[NucPair] = _nucs
        self._isEmpty:bool = False
        if _foldchange is not None:
            pass
        self.refresh()
        self._isEmpty=False

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
        just_pairs_str=[]
        for pair in self._nucPairListFull:
            justPairs.append(pair.snuppPair)
            just_pairs_str.append(pair.snuppPair.snuppPair)
        self._snuppOnly = justPairs
        self._snuppOnly_str = just_pairs_str

    def makePairsDict(self):
        pairsDict: OrderedDict[SnuppPair, NucPair] = {}
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
            else: # isinstance(value, List[NucPair]):
                secondList.makeNewPairList(value)    

            tempNucList: List[NucPair] = []
            for secondSnupp, secondNucPair in secondList.pairsDictFull.items():
                snuppPair: SnuppPair = secondSnupp.snuppPair
                #print(f'Appending Nuc Pair list for snuppPair {str(snuppPair)}')
                #num_pairs = self.snuppOnlyList.count(snuppPair)
                if snuppPair in self.snuppOnlyList_str:                    
                    #do fold change first so thta it is easier to handle probabilities dict
                    if  _foldchange is not None:
                        rounded_foldchange = round(_foldchange, 0)
                        search_list = secondList.pairsDictFull[secondSnupp].foldChangeList_rounded
                        num_times = search_list.count(rounded_foldchange)
                        if num_times == 0:
                            #if they are the same then just need to append the design ID
                            self.pairsDictFull[secondSnupp].appendFoldChange(_foldchange)

                    #the numbers are never close due to precision
                    #round to 6 decimal places?
                    secondList_prob = round(secondList.pairsDictFull[secondSnupp].probability, 10)
                    current_list_prob = round(self.pairsDictFull[secondSnupp].probability,10)

                    if secondList_prob == current_list_prob:
                        #if they are the same then just need to append the design ID
                        self.pairsDictFull[secondSnupp].appendDesignId(secondNucPair.designIDs)
                        tempNucList.append(self.pairsDictFull[secondSnupp])
                    else:
                        #if not the same pair prob then need to add both orignal and new to list as they are different
                        tempNucList.append(self.pairsDictFull[secondSnupp])
                        tempNucList.append(secondNucPair)            
        else:
            newNucList: List[NucPair] 
            if isinstance(value, NucPairList):
                newNucList = value.fullPairsList
            else:
                #it must be a List[NUCPAIR]
                newNucList = value 
            self.makeNewPairList(newNucList)

    def filter_out_NucPairList(self, filter_outList):

        if self._isEmpty == False:
            temp_filter_outList: NucPairList = NucPairList()
            if isinstance(filter_outList, NucPairList):
                temp_filter_outList = filter_outList
            else: # isinstance(value, List[NucPair]):
                temp_filter_outList.makeNewPairList(filter_outList)

            #so now you should have a new nucpair list from which you an work from
            #need to filter out nuc pairs that are similar
            #remember that this is only filtering out based on nuc pairs
            #do this based on snupp pairs list

            current_list_snupps_set: set = set(self.snuppOnlyList_str)
            filter_out_list_snupps_set: set =set(temp_filter_outList.snuppOnlyList_str)
            unique_snupps_current_list = list(current_list_snupps_set-filter_out_list_snupps_set)
            #unique_snupps_filter_out_list =  list(filter_out_list_snupps_set-current_list_snupps_set)

            #now remove anythng not only in current
            #going to need to iterate through

            temp_nucPairList: List[NucPair] = []

            #for nucpair_index in range(len(self._nucPairListFull)):
            #    nucPair: NucPair = self._nucPairListFull[nucpair_index]
            #    if nucPair.snuppPair.snuppPair in unique_snupps_current_list:
            #        #remove from current list then
            #        temp_nucPairList.append()

            for nucPair in self._nucPairListFull:
                #if the nucPair is not in the current only list then remove it
                if nucPair.snuppPair.snuppPair in unique_snupps_current_list:
                    temp_nucPairList.append(nucPair)
                    #remove from current list then
                    #self._nucPairListFull.remove(nucPair)

                # snupp_str:str = nucPair.snuppPair.snuppPair
                #num_count = unique_snupps_current_list.count(snupp_str)
                #if num_count > 0:
                ##if snupp_str in unique_snupps_filter_out_list:
                #    #remove from current list then
                #    self._nucPairListFull.remove(nucPair)
                #    pass
            
            #now you should have a list of just unique nuc pairs
            #so now do nucpair list housekeeping
            self._nucPairListFull = temp_nucPairList
            self.refresh()


    @property
    def pairsDictFull(self):
        return self._pairsDict

    @property
    def snuppOnlyList(self):
        return self._snuppOnly
    
    @property
    def snuppOnlyList_str(self):
        return self._snuppOnly_str

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
        self._unique_results_dict_float:  Dict[float, NucPairList] = {}
        self._parameterValuesList: List[float] = None
        self._minParameterValue: float = None
        self._maxParameterValue: float = None
        self._rangeParameterValue: float = None
        self._numSections: int = numSections
        self._listSectionValues: List[float] = None
        #self.appendResult(searchParameter, searchResult)

    def define_unique_section_nucs(self):
        #first get the list of section values and iterate through them
        #start at first group and then filter nucs based on the rest
        #and then go on from there

        #make a deap copy of the current dict
        self._unique_results_dict_float = copy.deepcopy(self._resultsDict_float)

        for value, current_keep_NucPairList in self._unique_results_dict_float.items():
            #now itereat through the dict
            for value_key, filter_out_nucpair_value in self._resultsDict_float.items():
                #dont do anything if it is the current one
                if value_key != value:
                    print(f'working on group {value} and filterin out nucs that appear in group {value_key}')
                    current_keep_NucPairList.filter_out_NucPairList(filter_out_nucpair_value)
            #now save teh temp list as the current list for the unique dict
            self._unique_results_dict_float[value]=current_keep_NucPairList



    def appendResult(self, searchParameter: float, searchResult: NucPairList, numSections: Optional[int] = None):
        print(f'Appending result to result object')
        if searchParameter not in self._resultsDict_float:
            self._resultsDict_float[searchParameter] = searchResult           
        else:
            primaryNucPairList: NucPairList = self._resultsDict_float[searchParameter]
            primaryNucPairList.appendNucPairList(searchResult, searchParameter)
            self._resultsDict_float[searchParameter]=primaryNucPairList

        self.createParameterList()
        self.get_min_max()
        if(numSections is not None):
            self._numSections = numSections
        #only defne sections is neccessary
        if self._numSections > 0:
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
    def unique_dict_results(self):
        return self._unique_results_dict_float

    
    def parameterResult(self, searchParameter:float):
        if (searchParameter not in self._resultsDict_float):
            return None
        else:     
            return self._resultsDict_float[searchParameter]

   
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

    def __init__(self, _searchType: SearchProtocol, _sourceType: NupackFoldDataEnum, lab: sara2Root.sublabData, minProbForPair: Optional[float] = 0, numberSections:Optional[int]=5) -> None:
        self.searchType: SearchProtocol = _searchType
        self.sourceType: NupackFoldDataEnum = _sourceType
        self.foldChangeDict: SearchResult = SearchResult(numSections=numberSections)
        self.rawLabData: sara2Root.sublabData = lab
        self._minProbForPair: float = minProbForPair
        self.initialize()
        #now a dictionary is loaded with all the fold changes and goodies

    def initialize(self):
        if self.searchType==SearchProtocol.FOLDCHANGE:
            if self.sourceType == NupackFoldDataEnum.PAIRPROBS: 
                #load fold change for lab into memory
                self.LoadLab()
        pass

    def get_dict(self):
        return self.foldChangeDict.dictResults

    #this will give a list of pairs based on prob level
    #remember that nupackfolddata is for a single design for a single lab
    def PairsSearch_SingleDesign(self, foldData: sara2Root.DesignPerformanceData, probThresh:float, doInit:Optional[bool]=False, foldChange:Optional[float]=None):
        print(f'Performing pairs search for {foldData.DesignInfo.DesignID}')
        pairsList = foldData.nupackFoldResults.pairprobsList
        designID = foldData.DesignInfo.DesignID
        
        tempSnuppPairsDict: Dict[NucPair, float] = {}  
        tempSnuppList: List[NucPair]= []
        #snuppPairsResults = self.SearchResultsGrouping

        for i in range(len(pairsList)):        
            for j in range(len(pairsList[i])):
                #only process pairs with first digit lower than first so you dont do them over
                #if i is not less than j then it has already been done
                if i < j:
                    pairValue = pairsList[i][j]
                    if pairValue >= probThresh and pairValue != 0.0:
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
        self.foldChangeDict.clear
        if self.searchType==SearchProtocol.FOLDCHANGE:
            for design in self.rawLabData.designsList:
                print(f'processing desing ID {design.DesignInfo.DesignID} for rainbow analysis')
                nupackData = design.nupackFoldResults
                wetlabData = design.wetlabResults                
                currentFoldChange = wetlabData.FoldChange
                currentPairDict, currentPairsList = self.PairsSearch_SingleDesign(foldData=design, probThresh=self._minProbForPair, doInit=True, foldChange=currentFoldChange)
                newNucPairList: NucPairList = NucPairList()
                newNucPairList.appendNucPairList(currentPairsList, currentFoldChange)
                self.foldChangeDict.appendResult(currentFoldChange, newNucPairList)    
            # now you should have a dictionary of fold changes with ech fold change having a list of nucpairs sorted by pairing prob
            # and the foldchange(s) the nuc pair is found in is(are) retained as well
        print("lab is now loaded into memory for analysis")

    #now need to write functions that filter through the loaded dic. the benifiet of the dict is that hopefully
    #we only need to do a few searches as they are all insorted order

    #need to write a function that is used as a base class for checking a value againsta a rna stucture vs finding commonalityin structures


    def new_foldChangeSearchPairs(self, foldChange_min: float, foldChange_max: float):
        #need to work off of fold change dict
        #self.foldChangeDict
        print("processing foldchang pairs")
        resultNucPairList: NucPairList = NucPairList()
        for foldChangeValue in self.foldChangeDict.keys:
            print(f'processing value {foldChangeValue}')
            if (foldChangeValue >= foldChange_min and foldChangeValue <= foldChange_max):
                #then this fold change needs to be recorded
                listFromDict: NucPairList = self.foldChangeDict.parameterResult(foldChangeValue)
                resultNucPairList.appendNucPairList(listFromDict, foldChangeValue)
                print(f'finished appending nuc pairs for value {foldChangeValue}')
        return resultNucPairList

    def generate_Marker_Grouping_For_Analysis(self, minFoldChangePair: Optional[float]=0):
        #first get the list of sections
        results: SearchResult = SearchResult()
        sectionNums: int
        sectionList: list[float]
        sectionNums, sectionList= self.foldChangeDict.sections

        for sectionIndex in range(len(sectionList)-1):
            print(f'processing section index {sectionIndex}')
            lowerSectionValue = sectionList[sectionIndex]
            upperSectionValue = sectionList[sectionIndex+1]
            commonPairsList: NucPairList = self.new_foldChangeSearchPairs(lowerSectionValue, upperSectionValue)
            if commonPairsList.empty == False:
                results.appendResult(upperSectionValue, commonPairsList)
        print(f'Defining unique nuc lists for search result')
        results.define_unique_section_nucs()
        return results


    def foldChangeSearch(self, roundInfo: sara2Root.sublabData, foldChange_min: float, foldChange_max: float, probMin: float):
        
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
    def FindPairs(self, roundInfo: sara2Root.sublabData, wetLabElement: string, wetLabElement_min: float, wetLabElemt_max: float, probMin: float):
        #do a search 
        #initially do fold change, no. of clusters, and EternScores
        if wetLabElement is "foldchange":
            pairsList = self.foldChangeSearch(roundInfo, wetLabElement_min, wetLabElemt_max, probMin)
        return pairsList

   
   
   
   
   #old stuff and maybe some still used code... idk
   
    #this is agrouping and what defines a grouping is done by the function that calls this
    #@dataclass
    #class SearchResultsGrouping:
    #    _RawResults_Dict: Dict = {}
    #    _SortedResults_Dict: Dict = {}
    #    Prob_Pair: OrderedDict[float, List[NucPair]] = {}      


    #def getRainbowColorMap(self):
    #    pass

    #find all pairs that have a parameter in common

    #def writePair(self, i: int, j: int):
    #    pairName = "{i}:{j}".format(i=i, j=j)
    #    return pairName
    
  
    
    
    

    
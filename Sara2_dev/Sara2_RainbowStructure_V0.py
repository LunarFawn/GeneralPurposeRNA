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
import Sara2_dev.Sara2_API_Python3_V1 as sara2Root



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

    def __init__(self, nuc1: int, nuc2: int, pairProb:float, designID_str: str)-> None:
        #convert to comon i, j pair
        self._pair:SnuppPair = SnuppPair(nuc1, nuc2)
        self._probability=pairProb
        #self._designId = designID_str
        self._designIdList: List[str]
        self._designIdList.append(designID_str)
    
    def __str__(self) -> str:
        return self._pair
    
    def __repr__(self):
        return f'NucPair(pair={self._pair}, prob={self._probability})'

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

    #@property
    #def designIDList(self):
    #    return self._designIdList
        
    def appendDesignId(self, value: List[str]):
        self._designIdList = self._designIdList + value
        



class NucPairList(object):
    
    def __init__(self, _nucs:List[NucPair]) -> None:
        self._nucPairListFull : List[NucPair] = _nucs
        self._sortedPairListFull: List[NucPair]
        self._snuppOnly: List[SnuppPair]
        self._pairsDict: OrderedDict[SnuppPair, NucPair]
        self.refresh()

    def refresh(self):
        self.sort_probability()
        self.makeJustPairs()
        self.makePairsDict()

    def makePairsDict(self):
        pairsDict = OrderedDict[SnuppPair, NucPair]
        for pair in self.fullPairsList:
            pairsDict[pair.snuppPair] = pair
        self._pairsDict = pairsDict


    def makeJustPairs(self):
        justPairs = []
        for pair in self._nucPairListFull:
            justPairs.append(pair.snuppPair)
        self._snuppOnly = justPairs

    def sort_probability(self):
        #sort list by probability
        sortedPairList: List[NucPair] = sorted(self._nucPairListFull, key=lambda x: x.probability)
        self._sortedPairListFull = sortedPairList

    #this will add two lists of NucPairs and only keep the pairs that are common
    def __add__(self, value):
        #first need see what nuc pairs and probs are common
        secondList: NucPairList
        if isinstance(value, NucPairList):
            secondList = value
        elif isinstance(value, List[NucPair]):
            secondList = NucPairList(value)        

        tempNucList: List[NucPair] = []
        for secondSnupp, secondNucPair in secondList.pairsDictFull.items():
            snuppPair: SnuppPair = secondSnupp.snuppPair
            if snuppPair in self.snuppOnlyList:
                if secondList.pairsDictFull[snuppPair].probability == self.pairsDictFull[snuppPair].probability:
                    #if they are the same then just need to append the design ID
                    self.pairsDictFull[snuppPair].appendDesignId(secondNucPair.designIDs)
                    tempNucList.append(self.pairsDictFull[snuppPair])
                else:
                    #if not the same pair prob then need to add both orignal and new to list as they are different
                    tempNucList.append(self.pairsDictFull[snuppPair])
                    tempNucList.append(secondNucPair)
        #now there shoudl be a new list so use that to init this one
        self.__init__(tempNucList)
   

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

class NupackFoldDataEnum(Enum):
    PAIRPROBS = 1
    #NOT IMPLEMENTED YET
    PAIRSLIST =2

class GenerateRainbowStructurePlot:

    def __init__(self, _searchType: SearchProtocol, _sourceType: NupackFoldDataEnum) -> None:
        self.searchType: SearchProtocol = _searchType
        self.sourceType: NupackFoldDataEnum = _sourceType
        self.initialize()

    def initialize(self):
        if self.searchType==SearchProtocol.FOLDCHANGE:
            if self.sourceType == NupackFoldDataEnum.PAIRPROBS: 
                #load fold change for lab into memory
                pass
        pass

    #this will load a datatype for furute searches to be faster
    def LoadLab(self, _lab: sara2Root.puzzleData):
        self.rawLabData=_lab
        self.foldChangeDict: Dict[float, Dict[float, List[NucPair]]]={}
        if self.searchType==SearchProtocol.FOLDCHANGE:
            for design in _lab.designsList:
                nupackData = design.nupackFoldResults
                wetlabData = design.wetlabResults
                currentPairDict, currentPairsList = self.PairsSearch_SingleDesign(nupackData, 0)
                currentFoldChange = wetlabData.FoldChange


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
    
  
    
    
    #this will give a list of pairs based on prob level
    #remember that nupackfolddata is for a single design for a single lab
    def PairsSearch_SingleDesign(self, foldData: sara2Root, probThresh:float, doInit:Optional[bool]=False):
        pairsList = foldData.NupackFoldData.pairprobsList
        designID = foldData.DesignInformation.DesignID
        
        tempSnuppPairsDict= dict(NucPair, float)        
        tempSnuppList = []
        #snuppPairsResults = self.SearchResultsGrouping

        for i in range(len(pairsList)):        
            for j in range(len(pairsList[i])):
                pairValue = pairsList[i][j]
                if pairValue >= probThresh:
                    pairName = NucPair(i,j, pairValue, designID)
                    if doInit==True:
                        pass
                    else:
                        #pairName = "{i}:{j}".format(i=i, j=j)
                        tempSnuppPairsDict[str(pairName)]=pairValue
                        tempSnuppList.append(pairName)
        return tempSnuppPairsDict, tempSnuppList     

    #need to write a function that is used as a base class for checking a value againsta a rna stucture vs finding commonalityin structures

    def foldChangeSearch(self, roundInfo: sara2Root.puzzleData, foldChange_min: float, foldChange_max: float, probMin: float):
        commonPairsList=[]
        
        #i dont care about the disgn name or anything right now just the pairs
        commonFoldChangeList=[]
        for designData in roundInfo.designsList:
            #designData is everything to do with each individual design in the lab
            #get the foldchange
            if((designData.wetlabResults.FoldChange >=  foldChange_min) and (designData.wetlabResults.FoldChange<=foldChange_max) ):
                #get nuc pairs
                pairDict, pairList = self.NucProbSearch(designData.nupackFoldResults, probMin)
                #now write to the list
                commonFoldChangeList.append(pairList)
        
        #now should have all the lists of pairs so do intersection and return just the pairs that are in common
        commonPairsList = list(set.intersection(*map(set, commonFoldChangeList)))
        return commonPairsList
        


    # make this have a argument passed that define sget/set
    def FindPairs(self, roundInfo: sara2Root.puzzleData, wetLabElement: string, wetLabElement_min: float, wetLabElemt_max: float, probMin: float):
        #do a search 
        #initially do fold change, no. of clusters, and EternScores
        if wetLabElement is "foldchange":
            pairsList = self.foldChangeSearch(roundInfo, wetLabElement_min, wetLabElemt_max, probMin)
        return pairsList

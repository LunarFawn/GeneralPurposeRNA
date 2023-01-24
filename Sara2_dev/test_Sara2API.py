import pytest

import currentStuff.nupackAPI_Sara2_Ver1 as nupackAPI
import currentStuff.nupackAPI_Sara2_Ver1 as nupy

import Sara2_dev.Sara2_API_Python3_V1 as saraApi
#from Sara2_API_Python3_V1 import Sara2 as saraApi
import Sara2_dev.Sara2_RainbowStructure_V0 as rainbow

import os

absolute_path = os.path.dirname(__file__)
relative_path = "srcData"
full_path = os.path.join(absolute_path, relative_path)

pnanFilename = 'pnas.2112979119.sd01.xlsx'
roundName = 'Round 7 (R101)'
path = os.path.join(full_path, pnanFilename)



def TestSara():

    results: saraApi.puzzleData = saraApi.Sara2.ProcessLab(path, roundName)      
    design: saraApi.DesignPerformanceData  = results.designsList[0]
    test = results.designsList[0].nupackFoldResults.pairprobsList
    assert design, "Design entry not built right"
    assert design.DesignInfo.Sequence, "sequence not defined"

    #this creates a search result with all teh nucs associated with the groupings defined by the number of sections
    newSearch: rainbow.GenerateRainbowStructurePlot = rainbow.GenerateRainbowStructurePlot(rainbow.SearchProtocol.FOLDCHANGE, rainbow.NupackFoldDataEnum.PAIRPROBS, results, 0, 10)
    markerSearchResult: rainbow.SearchResult = newSearch.generate_Marker_Grouping_For_Analysis()


TestSara()
import pytest
import time
from datetime import datetime, time, timedelta
#import currentStuff.nupackAPI_Sara2_Ver1 as nupackAPI
#import currentStuff.nupackAPI_Sara2_Ver1 as nupy

import Sara2_API_Python3_V1 as saraApi
#from Sara2_API_Python3_V1 import Sara2 as saraApi
import Sara2_RainbowStructure_V0 as rainbow

import os

absolute_path = os.path.dirname(__file__)
relative_path = "srcData"
full_path = os.path.join(absolute_path, relative_path)

train_Filename = "pnas_training_1-25-23.xlsx"
main_research = "pnas.2112979119.sd01.xlsx"
roundName = 'Round 7 (R101)'
path = os.path.join(full_path, train_Filename)



def TestSara():
    start_time = datetime.now()
    print(f'Start Time = {datetime.now()}')
    sara2_fuck = saraApi.Sara2()
    results: saraApi.puzzleData = sara2_fuck.ProcessLab(path, roundName)      
    design: saraApi.DesignPerformanceData  = results.designsList[0]
    test = results.designsList[0].nupackFoldResults.pairprobsList
    assert design, "Design entry not built right"
    assert design.DesignInfo.Sequence, "sequence not defined"

    #this creates a search result with all teh nucs associated with the groupings defined by the number of sections
    newSearch: rainbow.GenerateRainbowStructurePlot = rainbow.GenerateRainbowStructurePlot(rainbow.SearchProtocol.FOLDCHANGE, rainbow.NupackFoldDataEnum.PAIRPROBS, results, 0, 10)
    newDict = newSearch.get_dict()
    markerSearchResult: rainbow.SearchResult = newSearch.generate_Marker_Grouping_For_Analysis()    
    dict = markerSearchResult.dictResults
    unique_dict = markerSearchResult.unique_dict_results
    end_time = datetime.now()
    time_span: float = (end_time - start_time).total_seconds()
    thing = time_span / 60
    print("done")
    print(f'End time = {datetime.now()}')
    print(f'elapsed time is {thing} min total')
    print(f'elapsed time is {time_span} sec total')
    print("hi")

TestSara()
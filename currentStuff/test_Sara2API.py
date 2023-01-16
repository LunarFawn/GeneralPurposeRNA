import pytest

import nupackAPI_Sara2_Ver1 as nupackAPI
import Sara2_API_Python3 as sara2

import os

absolute_path = os.path.dirname(__file__)
relative_path = "srcData"
full_path = os.path.join(absolute_path, relative_path)

pnanFilename = 'pnas.2112979119.sd01.xlsx'
roundName = 'Round 7 (R101)'
path = os.path.join(full_path, pnanFilename)



def TestSara():

    results = sara2.ProcessLab(path, roundName)      
    design = results.designsList[0]
    test = results.designsList[0].nupackFoldResults.pairprobsList
    assert design, "Design entry not built right"
    assert design.Sequence, "sequence no defined"

TestSara()
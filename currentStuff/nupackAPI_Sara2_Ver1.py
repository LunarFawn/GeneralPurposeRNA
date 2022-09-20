from typing import List
import struct
from nupack import *
import pandas as pd
import sys
import openpyxl

rna_model='rna95-nupack3'    
# Define physical model 

temperature = 37
my_model = Model

def SetModel(rna_model, temp_C):
    temperature=temp_C
    my_model = Model(material=rna_model, celsius=int(temperature))

def GetPairProbs2DArray(mySequence):
    #convert into form the named touple is expecting
    pairs = List[List[float]]
    pairsMatrix = pairs(strands=mySequence, model=my_model)
    pairsArray = pairsMatrix.to_array()
    for i in range(len(pairsArray)):
        for j in range(len(pairsArray[i])):
            pairs[i][j]= pairsArray[i][j]

    return pairs


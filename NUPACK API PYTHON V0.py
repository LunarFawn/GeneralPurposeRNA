#NUPACK API base code at least initially in python
from nupack import *
from enum import Enum


import Sara2_BaseClass_0

class NUPACK_API:
    def run_single_sequence(sequence):
        #run single sequence
        new_dict = Sara2_BaseClass_0.DataCollection("test")
        new_dict.init_dict()
        new_dict.add_entry('sequence', sequence)
        test = new_dict.get_entry('sequence')
        result = test.get_value()
        result2 = 2


        #print("entry added")

        #now get teh entry info
        #test = new_dict.get_entry("sequence")
        #test_Dict = {1: Sara2_BaseClass_0.DataCollection, 2: Sara2_BaseClass_0.DataCollection}






    run_single_sequence("CUACGCUACUCUACGUAAGUAGUGUAGUGUAGUUUCGGAACUUAGUUGUGCAGAUUGUAUUGAUCGCUUAGGUAACUAUGUGAUC")







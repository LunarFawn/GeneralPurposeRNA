#NUPACK API base code at least initially in python
from nupack import *
from enum import Enum


import Sara2_BaseClass_0

class NUPACK_API:

    def run_single_sequence(sequence):
        # run single sequence
        new_dict = Sara2_BaseClass_0.DataCollection("test")
        new_dict.init_dict()
        new_dict.add_entry("sequence", sequence)
        print("DONE")

    run_single_sequence("CUACGCUACUCUACGUAAGUAGUGUAGUGUAGUUUCGGAACUUAGUUGUGCAGAUUGUAUUGAUCGCUUAGGUAACUAUGUGAUC")







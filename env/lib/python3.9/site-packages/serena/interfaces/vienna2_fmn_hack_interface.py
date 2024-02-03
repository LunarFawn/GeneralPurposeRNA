
import subprocess
from subprocess import CompletedProcess
from dataclasses import dataclass

from serena.utilities.ensemble_structures import Sara2SecondaryStructure, Sara2StructureList



class Vienna2FMNInterface():
    """
    class to interface with vienna2_fmn hack by Elnando888
    """
    @dataclass
    class CommandResponse():
        stdout:str
        stderr:str
        return_code:int

    def __init__(self) -> None:
        pass

    def rnafold_fmn(self,input_sequence:str,do_fmn:bool = True, fmn_amount:int = 200, fmn_delta:float = 0)-> Sara2SecondaryStructure:
        command = []
        if do_fmn is True:
            command= ["RNAfold", "--ligand", f'FMN:{fmn_amount}']
        else:
            command= ["RNAfold"]
        sequence: str = ''
        structure:str = ''
        energy:float = 0
        process = subprocess.run(command, input=input_sequence, encoding="utf-8", capture_output=True)
        raw_result = process.stdout.split('\n')
        if len(raw_result) > 0:
            sequence = raw_result[0]
            raw_struct = raw_result[1].split(' ')
            if len(raw_struct) > 0:
                structure = raw_struct[0]
                energy_str:str = raw_struct[1]
                energy_str = energy_str.strip('(')
                energy_str = energy_str.strip(')')
                try:
                    energy = float(energy_str)
                    energy = energy - fmn_delta
                except ValueError as error:
                        energy = 0
        
        sara_struct: Sara2SecondaryStructure = Sara2SecondaryStructure(sequence=sequence,
                                                                       structure=structure,
                                                                       free_energy=energy)

        return sara_struct

    def rnasubopt_fmn(self, input_sequence:str, do_fmn:bool = True, fmn_amount:int = 200):
        struct_list_response:Sara2StructureList = Sara2StructureList()
        sequence: str = ''
        command = []
        if do_fmn is True:
            command= ["RNAsubopt", "--ligand", f'FMN:{fmn_amount}']
        else:
            command= ["RNAsubopt"]
        subopt_response = self.run_command_locally(command=command,
                                            extra_input=input_sequence)
        raw_result = subopt_response.stdout.split('\n')
        found_sequence:bool = False
        if len(raw_result) > 1:
            sequence_raw = raw_result[0].split(' ')
            if len(sequence_raw) > 1:
                sequence = sequence_raw[0]
                found_sequence = True
        
        if found_sequence is True:
            #skip first one as it is the sequence
            for index in range(1, len(raw_result)):
                sara_struct:Sara2SecondaryStructure = Sara2SecondaryStructure()
                raw_struct = raw_result[index].split(' ')
                if len(raw_struct) > 1:
                    sara_struct.sequence = sequence
                    sara_struct.structure = raw_struct[0]
                    try:
                        energy = float(raw_struct[1])
                        sara_struct.free_energy = energy
                    except ValueError as error:
                        energy = 0
                    struct_list_response.add_structure(sara_struct)
        
        return struct_list_response
                


    def run_command_locally(self, command:str, extra_input:str ='/n')-> CommandResponse:
        process = subprocess.run(command, input=extra_input, encoding="utf-8", capture_output=True)
        response: self.CommandResponse = self.CommandResponse(stdout=process.stdout,
                                                    stderr=process.stderr,
                                                    return_code=process.returncode)
        return response


#new_viena: Vienna2FMNInterface = Vienna2FMNInterface()
#seq = 'GCCAUCGCAUGAGGAUAUGCUCCCGUUUCGGGAGCAGAAGGCAUGUCAUAAGACAUGAGGAUCACCCAUGUGGUUAAGAUGGCA'
#result = new_viena.rnasubopt_fmn(seq)
#print(result)
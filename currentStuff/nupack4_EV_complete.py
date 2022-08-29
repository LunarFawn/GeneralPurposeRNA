from nupack import *
import math
import copy


rna06_nupack4_parameter=["rna06", "Based on [Mathews99] and [Lu06] with additional parameters [Xia98,Zuker03] including coaxial stacking [Mathews99,Turner10] and dangle stacking [Serra95,Zuker03,Turner10] in 1M Na+."]
rna95_nupack4_parameters=["rna95", "Based on [Serra95] with additional parameters [Zuker03] including coaxial stacking [Mathews99,Turner10] and dangle stacking [Serra95,Zuker03,Turner10] in 1M Na+."]
rna99_nupack3_parameters=["rna99-nupack3","Parameters from [Mathews99] with terminal mismatch free energies in exterior loops and multiloops replaced by two dangle stacking free energies. Parameters are provided only for 37 ∘C."]
rna95_nupack3_parameters=["rna95-nupack3", "Same as rna95 except that terminal mismatch free energies in exterior loops and multiloops are replaced by two dangle stacking free energies."]

rna_parameterSet=[rna06_nupack4_parameter, rna95_nupack4_parameters, rna99_nupack3_parameters, rna95_nupack3_parameters]

print("Enter single strand RNA sequence")
sequence = "GGGAACGACUCGAGUAGAGUCGAAAAGUUGAAACGACUGAACAUGGUAACAUGAGUGGUUUGAACACAUACGAACAGGGUUCUUUCGAGGAUCCAAAAGAAACAACAACAACAAC" #input()

print("Enter the temperature to simulate at")
temperature = 37 #input()

#print("Select the energy parameter set you wish to use. enter the number that corresponds")
#for parameterSet_index in range(len(rna_parameterSet)):
##    parameterSet="#={0}, name={1}, description={2}".format(parameterSet_index, rna_parameterSet[parameterSet_index][0], rna_parameterSet[parameterSet_index][1])
 #   print(parameterSet)

#paramNum=int(input())

                                                           
#rna_model= rna_parameterSet[paramNum][0]                                                         

rna_model='rna95-nupack3'    
# Define physical model 

def energy_model(rna_model, temperature):
    model = Model(material=rna_model, celsius=int(temperature))
    return model


my_model = energy_model(rna_model, temperature)


# Partition function and complex free energy
partionFunction, freeEnergy_complex = pfunc(strands=sequence, model=my_model)

# Equilibrium pair probability matrix
probabilityMatrix_Raw = pairs(strands=sequence, model=my_model)

kcal_group_elements_list_global=[]
def get_subOpt_kcal_delta_elements(sequence_string, energymodel, energy_delta_MFE):
    #run through subopt 
    ensemble_kcal_group= subopt(strands=sequence_string, model=energymodel, energy_gap=energy_delta_MFE)
    
    #get all the data out of it
    kcal_group_elements_list=[]
    for i,kcal_group_elementInfo in enumerate(ensemble_kcal_group):        
        #get all the structures and energis pulled and prepped for proccessing and add them tot eh dict and the list
        element_info_Dict = {}        
        element_structure = str(kcal_group_elementInfo.structure)
        element_info_Dict["structure"] = element_structure
        element_freeEnergy = kcal_group_elementInfo.energy
        element_info_Dict["freeEnergy"] = element_freeEnergy
        element_stackEnergy = kcal_group_elementInfo.stack_energy
        element_info_Dict["stackEnergy"] = element_stackEnergy
        kcal_group_elements_list.append(element_info_Dict)
        
    
    kcal_group_elements_list_global = kcal_group_elements_list
    return kcal_group_elements_list



def get_EV(ensemble_kcal_group_elements, sequence):
    #need to do each char abd then structure
    #walk through each nucleotide but first prep containers grab what is needed
    
    #setup constants
    nuc_count = len(sequence)
    structure_element_count = len(ensemble_kcal_group_elements)
    
    
    ensembleVariation_score_temp = 0
    nucleotide_position_variation_basescores=[0]*nucleotide_count
    nucleotide_position_variation_subscores=[0]*nucleotide_count
    
    #we have all the elements and all the stuff for the group
    for nucIndex in range(nuc_count):
        for structureIndex in range(structure_element_count):
            mfe_nuc=ensemble_kcal_group_elements[0]["structure"][nucIndex]
            delta_nuc=ensemble_kcal_group_elements[structureIndex]["structure"][nucIndex]
            if mfe_nuc != delta_nuc:
                #then add a point since there is variation if nuc dotbracket not the same
                #this is a subscore for that nuc position only across the ensemble
                nucleotide_position_variation_basescores[nucIndex]=nucleotide_position_variation_basescores[nucIndex]+1
        #when done with all teh structures then do the nucScore wityh calculations
        nucleotide_position_variation_subscores[nucIndex]=nucleotide_position_variation_basescores[nucIndex] / structure_element_count
    #when that is all done the calculate final score
    total_subscore=sum(nucleotide_position_variation_subscores)
    
    #now do final calculations
    # original classic EV and normalized classic EV
    nucleotide_ensemble_variance_classic=total_subscore/nuc_count    
    fail_threshold=0.2
    nucleotide_ensemble_variance_normalized=(nucleotide_ensemble_variance_classic/fail_threshold)*100
    
    
    #new method that uses structure distance in nupack for points... distance is points
    
    energydelta_individualVariationScore_list=[]
    for element in ensemble_kcal_group_elements:
        element_structure = element["structure"]
        mfe_structure = ensemble_kcal_group_elements[0]["structure"]
        element_structure_mfe_distance=struc_distance(element_structure, mfe_structure)
        individualVariationScore = element_structure_mfe_distance / structure_element_count
        energydelta_individualVariationScore_list.append(individualVariationScore)
    
    structure_ensemble_variance = sum(energydelta_individualVariationScore_list) * 100
    return nucleotide_ensemble_variance_classic, nucleotide_ensemble_variance_normalized, structure_ensemble_variance
    
mfe_info = get_subOpt_kcal_delta_elements(sequence, my_model,0)[0]
mfe_structure=mfe_info["structure"]
mfe_freeEnergy=mfe_info["freeEnergy"]
mfe_stackEnergy=mfe_info["stackEnergy"]

#now calulate ensemble variance
print("enter a list of kcal deltas to look at seperated by a single space")
energyDeltas_raw= input()
deltasToInspect_list=energyDeltas_raw.split(" ")
deltacount=len(deltasToInspect_list)

deltaMessage="Getting EV score for {0} deltas: {1}".format(deltacount, deltasToInspect_list)

#Get ensemble variance data
EV_Delta_dict={}
nucleotide_count=len(sequence)
for delta in deltasToInspect_list:
    delta = int(delta)
    delta_elements=get_subOpt_kcal_delta_elements(sequence, my_model,delta) 
    nucleotide_ensemble_variance_classic, nucleotide_ensemble_variance_normalized, structure_ensemble_variance=get_EV(delta_elements, sequence)
    EV_Delta_dict[delta]={}
    structures_count=len(delta_elements)
    
    EV_Delta_dict[delta]["count"]=structures_count 
    EV_Delta_dict[delta]["struct_elements"]=delta_elements    
    EV_Delta_dict[delta]["EV_old_raw"]=nucleotide_ensemble_variance_classic
    EV_Delta_dict[delta]["EV_old_norm"]=nucleotide_ensemble_variance_normalized   
    EV_Delta_dict[delta]["EV_new"]=structure_ensemble_variance
    #now we have everything in easy access for the group of ensemble

#Get ensemble defect
ensemble_defect=defect(strands=sequence, structure=mfe_structure, model=my_model)

#print out the data now
#print MFE styff
mfe_string = "MFE Secondary Structure = {0}\nMFE Free Energy={1}\nMFE Stack Energy={2}".format(mfe_structure, mfe_freeEnergy,mfe_stackEnergy)
print(mfe_string)

#print out NUPACK ED data
ensemble_defect_string = "Ensemble Defect (ED) = {0}".format(ensemble_defect)
print(ensemble_defect_string)

#print out JenniferPearl EV data and ensemble strucutre info stuctures
for deltaIndex in range(len(EV_Delta_dict)):
    #get eh list of keys
    thisDelta=list(EV_Delta_dict.keys())[deltaIndex]    
    kcalDelta=thisDelta
    EV_old=EV_Delta_dict[thisDelta]["EV_old_raw"]
    EV_old_norm=EV_Delta_dict[thisDelta]["EV_old_norm"]
    EV_new=EV_Delta_dict[thisDelta]["EV_new"]
    strutureCount=EV_Delta_dict[thisDelta]["count"]
    deltaMessage = "Kcal_Delta={0}, EV_old={1}, EV_old_norm={2}, EV_new={3}, structure_count={4}".format(kcalDelta,EV_old,EV_old_norm,EV_new,strutureCount)
    print(deltaMessage)
    print("Ensemble kcal delta secondary structures and energy's")
    index = 1
    for element in EV_Delta_dict[thisDelta]["struct_elements"]:
        element_structure=element["structure"]
        element_freeEnergy=element["freeEnergy"]
        element_stackEnergy=element["stackEnergy"]
        struct_message="{0}: {1} FreeEnergy: {2:+.2f}: Stack Energy: {3:+.2f}".format(index, element_structure ,element_freeEnergy, element_stackEnergy)
        print(struct_message)
        index += 1
    print("\n")
print(probabilityMatrix_Raw.to_array())  
    
    
#need to dop some math on this stuff

def get_AdvancedEV(ensemble_kcal_group_elements, mfestructure):
    #need to do each char abd then structure
    #walk through each nucleotide but first prep containers grab what is needed
    
    #setup constants
    nuc_count = len(mfestructure)
    structure_element_count = len(ensemble_kcal_group_elements)
    
    
    ensembleVariation_score_temp = 0
    nucleotide_position_variation_basescores=[0]*nucleotide_count
    nucleotide_position_variation_subscores=[0]*nucleotide_count
    
    #we have all the elements and all the stuff for the group
    for nucIndex in range(nuc_count):
        for structureIndex in range(structure_element_count):
            mfe_nuc=mfestructure[nucIndex]
            delta_nuc=ensemble_kcal_group_elements[structureIndex]["structure"][nucIndex]
            if mfe_nuc != delta_nuc:
                #then add a point since there is variation if nuc dotbracket not the same
                #this is a subscore for that nuc position only across the ensemble
                nucleotide_position_variation_basescores[nucIndex]=nucleotide_position_variation_basescores[nucIndex]+1
        #when done with all teh structures then do the nucScore wityh calculations
        nucleotide_position_variation_subscores[nucIndex]=nucleotide_position_variation_basescores[nucIndex] / structure_element_count
    #when that is all done the calculate final score
    total_subscore=sum(nucleotide_position_variation_subscores)
    
    #now do final calculations
    # original classic EV and normalized classic EV
    nucleotide_ensemble_variance_classic=total_subscore#/nuc_count    
    fail_threshold=0.2
    nucleotide_ensemble_variance_normalized=(nucleotide_ensemble_variance_classic/fail_threshold)#*100
    
    
    #new method that uses structure distance in nupack for points... distance is points
    
    energydelta_individualVariationScore_list=[]
    for element in ensemble_kcal_group_elements:
        element_structure = element["structure"]
        mfe_structure = ensemble_kcal_group_elements[0]["structure"]
        element_structure_mfe_distance=struc_distance(element_structure, mfe_structure)
        individualVariationScore = element_structure_mfe_distance / structure_element_count
        energydelta_individualVariationScore_list.append(individualVariationScore)
    
    structure_ensemble_variance = sum(energydelta_individualVariationScore_list)# * 100
    return nucleotide_ensemble_variance_classic, nucleotide_ensemble_variance_normalized, structure_ensemble_variance

if deltasToInspect_list[0] == "2":
    #get each 1 kcal EV for entire range of value
    #to start it will be 5 for development purposes
    
    #initialize dictionary for lists of energy deltas
    energyDict={}
    energyRangeDict={}
    energyRangeDictCount={}
    #make +1 so that you can also get the EV at 0 energy delta so basically the index is teh actuall energy delta
    for eneryIndex in range(int(deltasToInspect_list[0])+2):
        energyDict[eneryIndex]=[]
        energyRangeDict[eneryIndex] = [] 
        energyRangeDictCount[eneryIndex] = [] 
           
    #now add to the dictionary
    mfeEnergy=0
    mfeStruct=""
    index=0    
    for deltaIndex in range(len(EV_Delta_dict)):
        #get eh list of keys        
        thisDelta=list(EV_Delta_dict.keys())[deltaIndex]    
        kcalDelta=thisDelta      
        for element in EV_Delta_dict[thisDelta]["struct_elements"]:            
            element_structure=element["structure"]
            element_freeEnergy=element["freeEnergy"]
            element_stackEnergy=element["stackEnergy"]
            tempList=[]
            if index==0:
                #this is the mfe and the base for everything record this and then move on
                mfeEnergy=element_stackEnergy
                mfeStruct=element_structure
                tempList=energyDict[0]
                tempList.append(element)
                energyDict[0]=tempList
            else:
                #this vreates teh individual energy list
                #energyDeltaNow=math.ceil(abs(mfeEnergy)-abs(element_freeEnergy))
                currentStackEnergy_Preped = abs(element_stackEnergy)
                mfeStackEnergy_Preped = abs(mfeEnergy); mfeEnergy
                energyDeltaNow=math.ceil(mfeStackEnergy_Preped - currentStackEnergy_Preped)

                
                
                tempList=energyDict[energyDeltaNow]
                tempList.append(element)
                energyDict[energyDeltaNow]=tempList
                
                #now work on the range energy list and use index for that                ``
            index=index+1
            
    #now that the energy dict is built lets build a energy range dict
    currentEnergyList = []
   
    for delta in range(int(deltasToInspect_list[0])+2):
        #at each pass the element list increases with each new delta. do the delta rtange stuff before each new one.. yay cryptic messageds
        count=0
        for element in energyDict[delta]:
            if count == 0:
                if delta == 0:
                    #do not record this one
                    pass
                else:
                    currentEnergyList.append(element)
            else:
                currentEnergyList.append(element)                               
        count= count+1;
        energyRangeDict[delta]=copy.deepcopy(currentEnergyList)
         #minus 1 from len since thre is the mfe in there)   
       
              
    #now should have a fully built energy dictionary need to do caclulations on each group
    for energyIndex in sorted(energyDict.keys()):
        energyDeltaCeil = energyIndex
        energyDeltaFloor=energyDeltaCeil-1
        #this is now a list of data structures so this next thing should work
        currentElementList=energyDict[energyIndex]
        rangeElementList = energyRangeDict[energyIndex]
        #make sure there is something there
        if len(currentElementList) > 0:                    
            #get EV for this chunk
            AdvancedEV_classic, AdvancedEV_normalized, AdvancedEV_structure=get_AdvancedEV(currentElementList, mfeStruct)
            if len(rangeElementList) > 0:                
                RangeAdvancedEV_classic, RangeAdvancedEV_normalized, RangeAdvancedEV_structure=get_AdvancedEV(rangeElementList, mfeStruct)
            structureCountEv=len(currentElementList)
            structurecountRange=len(rangeElementList)
            instantaniousEV_message="Energy Delta {0}: EV_Classic: {1}: EV_Normalized {2}: EV_Structure {3}: Number Elements {4}".format(energyDeltaCeil, AdvancedEV_classic ,AdvancedEV_normalized, AdvancedEV_structure, structureCountEv)
            print(instantaniousEV_message)      
            if energyIndex is not 0:
                rangeEV_message = "Energy Delta 1-{0}: EV_Classic: {1}: EV_Normalized {2}: EV_Structure {3}: Number Elements {4}".format(energyDeltaCeil, RangeAdvancedEV_classic ,RangeAdvancedEV_normalized, RangeAdvancedEV_structure, structurecountRange)
                print(rangeEV_message)
            else:
                #do nothing
                 pass
            print("\n")       
        else:
            AdvancedEV_classic=0
            AdvancedEV_normalized=0
            AdvancedEV_structure=0
            RangeAdvancedEV_classic=0
            RangeAdvancedEV_normalized=0
            RangeAdvancedEV_structure=0
            structureCountEv=1
            structurecountRange=1
            rangeEV_message="none"
            instantaniousEV_message="none"          
        
         
       
       
       
import nupack4_research_api as nupppy

def mechpaperRun():
    sequence=input("What is the sequence?")
    temp = input("what temp?")
    cutoff = input("what threshold cuttoff for pairs?")
    
    if temp == "d":
        temp = "37"
    
    if cutoff == "d":
        cutoff="0.0001"
    else:
        cutoff = float(cutoff)
    pairsdict, pairs, snuppPairs = nupppy.getPairProbs(sequence,temp,cutoff)
    
    pairNum=1
    for pair in snuppPairs:
        pairdata = "pair_{0}  {1}".format(pairNum, pair)
        print(pairdata)
        pairNum += 1
    

mechpaperRun()


import rungeofold as gf
import numpy as np

if __name__ == "__main__": main()

def main():
    pass

def getScore(pdb,itr):
    '''after GeoFold runs, calculates the score for this iteration'''
    try:
        with open("stf/%s.pair"%(pdb)) as fin:
            pairs = [line.split() for line in fin]
    except:
        hbfile = "tmp/%s_%i/%s_%i.hb"%(pdb,itr,pdb,itr)
        hxfile = "stf/%s.hx"%(pdb)
        output = "stf/%s.pair"%(pdb)
        pairs = hxToPair(hxfile,hbfile,output)

def getTotalScore(scores):
    '''given the scores for each entry, calculate the total score of the loss
    function'''
    return (np.array(scores).sum())/len(scores)



def run(breakcut,pivotcut,hingecut,itr=0):
    '''runs everything and calculates loss function for a given set of breakcut,
    pivotcut, and hingecut'''
    scores = []
    pdbs = [] #The test set of pdbs we're using and the appropriate temperatures to run them
    for (pdb,temp) in pdbs:
        writeParamFile(pdb,temp,breakcut,pivotcut,hingecut,itr)
        gf.main(["","parameters/%s_%i"%(pdb,itr),"walcob.conf"])
        scores.append(getScore(pdb,temp,breakcut,pivotcut,hingecut))
    return getTotalScore(scores)


def hxToPair(hxfile,hbfile,output):
    '''given the hx file from start2fold and the h-bond file, determine the contact
    pairs that we have data for as well as their experimental ages.  Note: This
    only needs to occur once and should be read from an appropriate output file'''
    pass

def writeParamFile(pdb,temp,breakcut,pivotcut,hingecut,itr=0):
    '''writes the parameters file for the given PDB file, temperature, and other
    parameters'''
    fout = open("parameters/%s_%i"%(pdb,itr),"w+")

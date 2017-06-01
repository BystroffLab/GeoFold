import rungeofold as gf
import numpy as np

if __name__ == "__main__": main()

def main():
    pass

def getScore(pdb,itr):
    '''after GeoFold runs, calculates the score for this iteration'''
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

def parseHBFile(hbfile):
    '''Takes the HBFile and parses it into a list of pairs of backbone hbonds'''
    output = {} # {donor:[list of residues]}
    with fin as open(hbfile):
        for line in fin:
            if line[0] != "!":
                line = [x.strip() for x in line.split()]
                if line[1] == 'N' and line[3] == 'O' and line[4] == 'H':
                    try:
                        output[int(line[0])].append(int(line[2]))
                    except:
                        output[int(line[0])] = [int(line[2])]
    return output

def parseHXFile(hxfile):
    '''Takes the HXFile and parses it into a list of residues and their relative age'''
    rgbToNum = {'R':1,'G':2,'B':3}
    with open(hxfile) as fin:
        output = [line.split() for line in fin]
    output = [[int(resn),rgbToNum(rgb)] for [resn,rgb] in output]
    return output


def hxToPair(hxfile,hbfile,output):
    '''given the hx file from start2fold and the h-bond file, determine the contact
    pairs that we have data for as well as their experimental ages.  Note: This
    only needs to occur once and should be read from an appropriate output file'''
    #HX FILE FORMAT
    #Resn [R || G || B]
    try:
        pairs = [line.split() for line in open(output)]
    except:
        hx = parseHXFile(hxfile)
        hb = parseHBFile(hbfile)
        pairs = []
        for [resn,rgb] in hx:
            for pair in hb[resn]:
                pairs.append([resn,pair,rgb])
    return pairs

def writeParamFile(pdb,temp,breakcut,pivotcut,hingecut,itr=0):
    '''writes the parameters file for the given PDB file, temperature, and other
    parameters'''
    fout = open("parameters/%s_%i"%(pdb,itr),"w+")
    fout.write("LNAME %s_%i\nEMAIL walcob@rpi.edu\nPDBCODE %s\nOMEGA 1.\n"%(pdb,itr,pdb))
    fout.write("INTERMEDIATES 0\nBARRELMOVES 0\nORANGE 1.0\nRUNGEOFOLD 1\n")
    fout.write("MOLSCRIPT 0\nBREAKCUT %f\nPIVOTCUT %f\nHINGECUT %f\n"%(breakcut,pivotcut,hingecut))
    fout.write("SEAMCUT 10\nBREAKPOINTENTROPY 90.\nHINGEPOINTENTROPY 30.\n")
    fout.write("TEMPERATURE %f\nCONCENTRATION 1.\nVOIDENTROPY 0.\nSOLIDITY 1000\n"%(temp))
    fout.write("HBONDENERGY 100.\nHAMMONDSCALE 1000.\nSIDECHAINENTROPY 1.\nHINGEBARRIER 0.\n")
    fout.write("PIVOTBARRIER 0.\nWATER 1.0\nMAXSPLIT 4\nMAXTIME 10.\nMINSEG 4\n")
    fout.write("CAVITATION 0.000001\nFLORY 0\nFLORYW 8.\nREDUCING 0\nHLFE 0\n")
    fout.write("FING 0\nCHAIN .\n")
    fout.close()

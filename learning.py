import rungeofold as gf
import numpy as np
import os


def purge(pdb,itr):
    from commands import getstatusoutput as run
    run("rm -rfv tmp/%s_%i output/%s_%i"%(pdb,itr,pdb,itr))

def main():
    scores = {}
    itr = 0
    breakcut = 0.
    while breakcut <= 1.:
        scores[breakcut] = run(breakcut,0.25,0.5,itr)
        breakcut += 0.1
        itr += 1
    output = open("breakcut.dat","w+")
    for bc in scores.keys().sort():
        output.write("%f, %f\n"%(bc,scores[bc]))
        print "%f, %f"%(bc,scores[bc])
    output.close()

def getScore(pdb,itr):
    '''after GeoFold runs, calculates the score for this iteration'''
    os.chdir("/Users/walcob/GeoFold")
    hxfile = "%s.hx"%(pdb)
    hbfile = "tmp/%s_%i/%s_%i.hb"%(pdb,itr,pdb,itr)
    output = "%s.pair"%(pdb)
    pairs = hxToPair(hxfile,hbfile,output)
    agefile = "tmp/%s_%i/%s_%i_1.dag.age"%(pdb,itr,pdb,itr)
    ages = parseAgeFile(agefile)
    age_pairs = []
    for [res1,res2,xage] in pairs:
        try:
            age_pairs.append([res1,res2,xage,ages[(res1,res2)]])
        except KeyError:
            age_pairs.append([res1,res2,xage,10.])
    age_pairs.sort(key=lambda x: x[3])
    late,early,intermediate = countLEI(age_pairs)
    last_age = (-1,0)
    for i in range(early):
        if age_pairs[i][3] == last_age[0]:
            age_pairs[i][3] = last_age[1]
        else:
            last_age = (age_pairs[i][3],1)
            age_pairs[i][3] = 1
    for i in range(early,intermediate+early):
        if age_pairs[i][3] == last_age[0]:
            age_pairs[i][3] = last_age[1]
        else:
            last_age = (age_pairs[i][3],2)
            age_pairs[i][3] = 2
    for i in range(intermediate+early,intermediate+early+late):
        if age_pairs[i][3] == last_age[0]:
            age_pairs[i][3] = last_age[1]
        else:
            last_age = (age_pairs[i][3],3)
            age_pairs[i][3] = 3
    score = 0.
    for [res1,res2,xage,gage] in age_pairs:
        score += float((xage-gage)**2)
    score = score/float(len(age_pairs))
    print score
    return score

def countLEI(age_pairs):
    late = 0
    early = 0
    intermediate = 0
    for [x,y,z,a] in age_pairs:
        if z == 1: early += 1
        if z == 2: intermediate += 1
        if z == 3: late += 1
    return late,early,intermediate

def parseAgeFile(agefile):
    age_dict = {}
    with open(agefile) as fin:
        ages = [line.split() for line in fin]
    for [x,y,age] in ages:
        age_dict[(int(x),int(y))] = float(age)
        age_dict[(int(y),int(x))] = float(age)
    return age_dict

def getTotalScore(scores):
    '''given the scores for each entry, calculate the total score of the loss
    function'''
    return (np.array(scores).sum())/float(len(scores))

def readPDBList(pdblist):
    with open(pdblist) as fin:
        return [tuple(line.split()) for line in fin]


def run(breakcut,pivotcut,hingecut,itr=0):
    '''runs everything and calculates loss function for a given set of breakcut,
    pivotcut, and hingecut'''
    scores = []
    pdbs = readPDBList("pdblist.txt") #The test set of pdbs we're using and the appropriate temperatures to run them
    for (pdb,temp) in pdbs:
        print pdb
        writeParamFile(pdb,temp,breakcut,pivotcut,hingecut,itr)
        gf.main(["","parameters/%s_%i"%(pdb,itr),"walcob.conf"])
        os.chdir("/Users/walcob/GeoFold")
        scores.append(getScore(pdb,itr))
        purge(pdb,itr)
    return getTotalScore(scores)

def parseHBFile(hbfile):
    '''Takes the HBFile and parses it into a list of pairs of backbone hbonds'''
    output = {} # {donor:[list of residues]}
    with open(hbfile) as fin:
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
    output = [[int(resn),rgbToNum[rgb]] for [resn,rgb] in output]
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
            try:
                for pair in hb[resn]:
                    pairs.append([resn,pair,rgb])
            except:
                pass
    return pairs

def writeParamFile(pdb,temp,breakcut,pivotcut,hingecut,itr=0):
    '''writes the parameters file for the given PDB file, temperature, and other
    parameters'''
    fout = open("parameters/%s_%i"%(pdb,itr),"w+")
    fout.write("LNAME %s_%i\nEMAIL walcob@rpi.edu\nPDBCODE %s\nOMEGA 1.\n"%(pdb,itr,pdb))
    fout.write("INTERMEDIATES 0\nBARRELMOVES 0\nORANGE 1.0\nRUNGEOFOLD 1\n")
    fout.write("MOLSCRIPT 0\nBREAKCUT %f\nPIVOTCUT %f\nHINGECUT %f\n"%(breakcut,pivotcut,hingecut))
    fout.write("SEAMCUT 10\nBREAKPOINTENTROPY 90.\nHINGEPOINTENTROPY 30.\n")
    fout.write("TEMPERATURE %s\nCONCENTRATION 1.\nVOIDENTROPY 0.\nSOLIDITY 1000\n"%(temp))
    fout.write("HBONDENERGY 100.\nHAMMONDSCALE 1000.\nSIDECHAINENTROPY 1.\nHINGEBARRIER 0.\n")
    fout.write("PIVOTBARRIER 0.\nWATER 1.0\nMAXSPLIT 4\nMAXTIME 10.\nMINSEG 4\n")
    fout.write("CAVITATION 0.000001\nFLORY 0\nFLORYW 8.\nREDUCING 0\nHLFE 0\n")
    fout.write("FING 0\nCHAIN .\n")
    fout.close()


if __name__ == "__main__": main()

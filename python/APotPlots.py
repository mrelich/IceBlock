
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#
# Make vector potential plots and the electric field histograms #
# from the output of Geant4. The goal is to be dynamic for N    #
# antennas in the output file...                                #
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#

import sys
from tools import *
from ROOT import *

rootfiledir = "efieldroot/"

#---------------------------------------------#
# Unfortunately I need some human input here.
# Add the number of bins to save in figure
# including the step size.
#---------------------------------------------#

nbins = 2000
#tStep = 0.015 # ns
#tStep = 0.02 # ns
tStep = 0.005 # ns

#---------------------------------------------#
# Specify input variables necessary here
#---------------------------------------------#

args = sys.argv
if len(args) != 3:
    print "Need to pass two arguments for"
    print "the input file and output file name"
    print "Number of args given: ", len(args) - 1
    sys.exit()

inFileName = ""
outFileName = ""
for arg in args:
    if ".dat" in arg:  inFileName = arg
    if ".root" in arg: outFileName = arg
    
# Check we have input and output name
if inFileName == "": 
    print "Input file name not specified"
    sys.exit()
if outFileName == "":
    print "Output file name not specified"
    sys.exit()

# Open files
infile  = open(inFileName,"r")
outfile = TFile(rootfiledir + outFileName, "recreate")

#---------------------------------------------#
# Define some user methods to make
#---------------------------------------------#

def profName(antNum, x, y, z):
    name = "AntNum_"+str(antNum)
    name += "_pos_"
    name += str(x) + "_" 
    name += str(y) + "_" 
    name += str(z)
    return name

#---------------------------------------------#
# Loop and read in data
#---------------------------------------------#

prof_VPot_X = []
prof_VPot_Y = []
prof_VPot_Z = []
prof_VPot   = []

prof_VPot_perEvent = []

currentEvent = -1
curAnt = 0
curAntName = ""
for line in infile:

    # Skip lines in event
    if "Event" in line:
        currentEvent += 1
        if currentEvent % 5 == 0: 
            print "Current Event: ", currentEvent
        eventVector = []
        prof_VPot_perEvent.append( eventVector )
        continue



    # Specify what profile to fill
    # based on the current antenna
    if "Antenna" in line:
        sl = line.split()
        curAnt     = int(sl[1])
        curAntName = profName(curAnt, sl[3], sl[4], sl[5])
        continue
    
    # Parse the input line
    sl = line.split()
    t  = float(sl[0])
    Ax = float(sl[1])
    Ay = float(sl[2])
    Az = float(sl[3])

    # Now check if antenna exists. If 
    # not create it with antenna name
    if len(prof_VPot) <= curAnt:
        VPName = "A_" + curAntName
        prof_VPot_X.append( makeProfile(VPName+"_X", nbins, t,
                                        t+nbins*tStep,
                                        kBlue,
                                        25,
                                        "time [ns]",
                                        "A_{x} [Vs/m]") )
        
        prof_VPot_Y.append( makeProfile(VPName+"_Y", nbins, t,
                                        t+nbins*tStep,
                                        kBlue,
                                        25,
                                        "time [ns]",
                                        "A_{y} [Vs/m]") )

        prof_VPot_Z.append( makeProfile(VPName+"_Z", nbins, t,
                                        t+nbins*tStep,
                                        kBlue,
                                        25,
                                        "time [ns]",
                                        "A_{z} [Vs/m]") )

        prof_VPot.append( makeProfile(VPName, nbins, t,
                                      t+nbins*tStep,
                                      kBlue,
                                      25,
                                      "time [ns]",
                                      "A [Vs/m]") )

    # Now handle the per event level results
    if len(prof_VPot_perEvent[currentEvent]) <= curAnt:
        VPEvtName = "A_" + curAntName + "_Event" + str(currentEvent)
        #print "Making: ", VPEvtName
        prof_VPot_perEvent[currentEvent].append( makeProfile(VPEvtName,nbins,t,
                                                             t+nbins*tStep,
                                                             kBlue,
                                                             25,
                                                             "time [ns]",
                                                             "A [Vs/m]") )
                                                 
                                                             
    
    # Fill Antenna vector potential
    prof_VPot_X[curAnt].Fill(t,Ax)
    prof_VPot_Y[curAnt].Fill(t,Ay)
    prof_VPot_Z[curAnt].Fill(t,Az)
    prof_VPot[curAnt].Fill(t,sqrt(Ax*Ax+Ay*Ay+Az*Az))
    prof_VPot_perEvent[currentEvent][curAnt].Fill(t,sqrt(Ax*Ax+Ay*Ay+Az*Az))

# Clean up
infile.close()
outfile.Write()
outfile.Close()

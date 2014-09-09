
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

nbins = 1000
tStep = 0.02 # ns

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
prof_EField = []

currentEvent = -1
curAnt = 0
curAntName = ""
for line in infile:

    # Skip lines in event
    if "Event" in line:
        currentEvent += 1
        if currentEvent % 5 == 0: 
            print "Current Event: ", currentEvent
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
    
    # Fill Antenna vector potential
    prof_VPot_X[curAnt].Fill(t,Ax)
    prof_VPot_Y[curAnt].Fill(t,Ay)
    prof_VPot_Z[curAnt].Fill(t,Az)
    prof_VPot[curAnt].Fill(t,sqrt(Ax*Ax+Ay*Ay+Az*Az))


#----------------------------------------------------#
# Now calculate the electric field from the profiles
# Whereas this will represent that average electric
# field from the average vector potential
#----------------------------------------------------#
def makeEField(profList):
    
    outfile.cd()
    for prof in profList:
    
        # Specify the name
        name  = prof.GetName()
        sl = name.split("A_")
        EName = "E_" + sl[1]

    
        # Make the new histogram
        t = prof.GetXaxis().GetBinLowEdge(1)
        Efield = makeProfile(EName, 
                             nbins,
                             t,
                             t+nbins*tStep,
                             kRed,
                             20,
                             "time [ns]",
                             "E [V/m]")

        # Loop over bins
        for i in range(nbins):
        
            binWidth = prof.GetBinWidth(i)
            time     = prof.GetBinCenter(i)
            A0       = prof.GetBinContent(i)
            A1       = prof.GetBinContent(i+1)
            #print "Bin width: ", binWidth, " time: ", time, " A0: ", A0, " A1: ", A1
            Efield.Fill(time, (A1-A0)/(binWidth*1e-9))

        # end loop over bins
        Efield.Write()
        Efield.Delete()
        
    # end loop over profiles
    return

makeEField(prof_VPot_X)
makeEField(prof_VPot_Y)
makeEField(prof_VPot_Z)
makeEField(prof_VPot)

# Clean up
infile.close()
outfile.Write()
outfile.Close()

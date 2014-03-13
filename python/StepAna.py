
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=#
# Analysis of the step output. This will include:
#    * Average energy of a Electron 
#    * Moliere Radius plots
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#

import sys
import os
pwd = os.environ['PWD']
sys.path.insert(0, pwd+'/python')

from tools import * 
from parseStepFile import *
from math import sqrt
from ROOT import *

#
## Specify input file
#

if len(sys.argv) < 2:
    print "Must specify input file name"
    sys.exit()

# Parse input file
fname = sys.argv[1]
f_energy = fname.split("_")[2]
f_nEvent = fname.split("_")[1]
f_mat    = fname.split("_")[3]

infile = open(fname,"r")
output = TFile("rootfiles/StepAna_"+f_nEvent+"_"+f_energy+"_"+f_mat+".root","recreate")

#
## Constants
#

threshold = 0.611 #100             # MeV
beamE     = float(f_energy) # MeV
print "Beam energy: ", beamE

# Specify binning for moliere radius
# Ice is default
moliereBins  = 100
moliereStep  = 0.5
if f_mat == "iron": 
    moliereBins = 100
    moliereStep = 0.1
elif f_mat == "lead":
    moliereBins = 100
    moliereStep = 0.1
    

#
## Plots and constants
#

p_dVsE = makeProfile("dVsE",50,0,200,kBlack,20,"Distance [cm]", "Energy [GeV]")
p_rVsShowerE = makeProfile("rVsShowerE",moliereBins, 0, moliereStep*moliereBins,
                           kBlack, 20, "Radius [cm]", "U/E_{0}")

#
## Useful Methods and placeholders
#

def getMoliereBins():
    return [0] * moliereBins

def getRad(x,y):
    return sqrt(x*x + y*y)

MoliereE = []
prevTrkID = -999

#
## Loop
#
for line in infile:
    
    # Initialize stuff at beginning of evnet
    if "Event" in line:
        MoliereE = getMoliereBins()
        prevTrkID = -999
        continue

    # End of event. Fill profiles for event
    if "End" in line:
        for i in range(1, len(MoliereE)+1):
            p_rVsShowerE.Fill(i*moliereStep, MoliereE[i-1])
        continue

    # Load some constants
    PDG    = getStepPDG(line)
    trkID  = getStepTrkID(line)
    energy = getStepE(line)
    dE     = getStepdE(line)

    # Plots right now are for electron only
    if abs(PDG) != 11: continue

    # Save Moliere Information for the
    # charged particles
    if prevTrkID != trkID:
        prevTrkID = trkID

    # Get the Radius for this point
    if trkID != 1:
        radius  = getRad(getStepX(line),getStepY(line))
        for i in range(1,moliereBins+1):
            if not radius > i*moliereStep: continue
            MoliereE[i-1] += (dE)/beamE

    # Place some energy threshold cut
    if energy < threshold: continue

    # Fill for leading particle only
    if trkID == 1:
        p_dVsE.Fill(getStepZ(line),energy/1000.0)        

# end loop over file

# Clean up
output.Write()
output.Close()
infile.close()
        

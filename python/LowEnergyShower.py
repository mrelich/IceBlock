
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
# So I now realize that 40 MeV showers don't really create an    #
# excess profile on their own.  Their showers are actually       #
# quite weird, where the initial particle ranges out before man  #
# secondaries are created and range out further. I want to study #
# those properties.                                              #
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#

from tools import * 
from parseTrackFile import *
from ROOT import *
import sys
import os
from array import array

pwd = os.environ['PWD']
sys.path.insert(0,pwd+'/python')

#
## Specify number of desired events
#

maxEvents = 5000

#
## Specify the input file
#

if len(sys.argv) < 2:
    print "Must specify input file name"
    sys.exit()

# Parse input file
fname = sys.argv[1]
f_energy = fname.split("_")[2]
f_nEvent = fname.split("_")[1]
f_mat    = fname.split("_")[3]
f_part   = fname.split("_")[4]
f_npart  = (fname.split("_")[5]).split(".")[0]

infile = open(fname,"r")
output = TFile("rootfiles/LowEnergyStudy_"+f_energy+"_"+f_mat+"_"+f_part+"_"+f_npart+".root","recreate")

#
## Plots to make
#

# Shower depths
stepsize = 0.5 # cm
zmin     = 0   #cm
zmax     = 200
nzbins   = int(zmax/stepsize)

h_initialZ  = makeHist("initialZ",nzbins,zmin,zmax,kBlack,20,"z_{i} [cm]","Entries")
h_finalZ    = makeHist("finalZ",nzbins,zmin,zmax,kBlack,20,"z_{f} [cm]","Entries")
h_dZ        = makeHist("dZ",100,0,25,kBlack,20,"#Deltaz [cm]","Entries")
h_E         = makeHist("E",100,0,50,kBlack,20,"E [MeV]","Entries")
h_EvsZI     = TH2F("EvsZI","",nzbins,zmin,zmax,100,0,50)

p_NPartSum  = makeProfile("NPartSum",nzbins,zmin,zmax,kBlack,20,"z_{i} [cm]","Entries")
p_NPartDiff = makeProfile("NPartDiff",nzbins,zmin,zmax,kBlue,25,"z_{i} [cm]","Entries")

#
## Useful for profiles
#

def getZbins():
    return [0] * nzbins

NPartSum  = []
NPartDiff = []

#
## Loop 
#
counter = 0
for line in infile:

    # Limit events for now
    if counter >= maxEvents: break

    # Skip event line
    if "Event" in line:
        NPartSum  = getZbins()
        NPartDiff = getZbins()
        continue
    
    # If end of event, skip
    if "End" in line:
        counter +=1
        
        # Fill profiles
        for i in range(len(NPartSum)):
            p_NPartSum.Fill(i*stepsize, NPartSum[i])
            p_NPartDiff.Fill(i*stepsize, NPartDiff[i])

        continue
    
    # Only look at charged particles
    pdg = getTrkPDG(line)
    if not abs(pdg) == 11: 
        continue

    # Get Variables
    zi = getTrkZ(line)
    zf = getTrkZFinal(line)
    E  = getTrkE(line)

    # Fill plots
    h_initialZ.Fill(zi,1)
    h_finalZ.Fill(zf,1)
    h_dZ.Fill(zf-zi)
    h_E.Fill(E)
    h_EvsZI.Fill(zi,E)

    # Save variables for profiles
    int_zi = int(zi/0.5)
    int_zf = int(zf/0.5)

    if not int_zi == int_zf:
        if int_zf > nzbins: int_zf = nzbins
        for i in range(int_zi, int_zf):
            NPartSum[i] += 1
            if pdg == 11:    NPartDiff[i] += 1
            elif pdg == -11: NPartDiff[i] -= 1
#
## Let user know events processed
#    
print "Events processed: ", counter

#
## Save file
#

output.Write()
output.Close()

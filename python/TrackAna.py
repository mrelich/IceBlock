
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
# Analysis of the track output. So far this will include #
# looking at the total track length as a function of     #
# the threshold cut and the counting of NParticles as    #
# function of radiation lengths.                         #
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#

from tools import * 
from parseTrackFile import *
from ROOT import *
import sys
import os
from array import array

pwd = os.environ['PWD']
sys.path.insert(0, pwd+'/python')

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
f_part   = (fname.split("_")[4]).split(".")[0]

infile = open(fname,"r")
#output = TFile("rootfiles/TrkAna_"+f_nEvent+"_"+f_energy+"_"+f_mat+"_"+f_part+"_5MeV.root","recreate")
output = TFile("rootfiles/TrkAna_"+f_nEvent+"_"+f_energy+"_"+f_mat+"_"+f_part+".root","recreate")

#
## Constants
#

radLength = 0
threshold = 0

if "ice" == f_mat:
    radLength = 39.0522 # cm
    threshold = 0.611   # MeV
elif "iron" == f_mat:
    radLength = 1.75749 # cm
    threshold = 100     # MeV
elif "lead" == f_mat:
    radLength = 5.61253 # cm
    threshold = 100     # MeV
else:
    print "Material not recognized"
    sys.exit()

#
## Plots and constants
#

# Nparticles
stepsize = 0.5 # Determines bin size for NPart plot
npMin    = 0
npMax    = 20
npBins   = int(npMax/stepsize)
npXTitle = "Radiation Lengths"
npYTitle = "<NParticles>"
p_NPartSum  = makeProfile("NPartSum",npBins, npMin, npMax, kBlack, 20, npXTitle, npYTitle)
p_NPartDiff = makeProfile("NPartDiff",npBins, npMin, npMax, kBlue, 25, npXTitle, npYTitle)

# Track length vs. Kinetic energy
kinBins = array('f',[0,0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                     1,1.5,2,3,4,5,6,7,8,9,10])
kinXTitle = "Kinetic Thershold [MeV]"
kinYTitle = "Track Length [m]"
p_TrkLVsE_tot  = makeProfileVar("TrkLVsE_tot",kinBins, kBlack, 20, kinXTitle, kinYTitle)
p_TrkLVsE_sum  = makeProfileVar("TrkLVsE_sum",kinBins, kBlack, 20, kinXTitle, kinYTitle)
p_TrkLVsE_diff = makeProfileVar("TrkLVsE_diff",kinBins, kBlack, 20, kinXTitle, kinYTitle)

#
## Useful Methods and placeholders
#

def getNParticles():
    return [0] * npBins

def getTrkBins():
    return [0] * len(kinBins)

nPartSum  = []
nPartDiff = []
trkL_tot  = []
trkL_sum  = []
trkL_diff = []

#
## Loop
#

for line in infile:
    
    # Initialize stuff at beginning of evnet
    if "Event" in line:
        nPartSum  = getNParticles()
        nPartDiff = getNParticles()
        trkL_sum  = getTrkBins()
        trkL_diff = getTrkBins()
        trkL_tot  = getTrkBins()
        continue

    # End of event. Fill profiles for event
    if "End" in line:
        for np in range(len(nPartSum)):
            p_NPartSum.Fill(stepsize*np, nPartSum[np])
            p_NPartDiff.Fill(stepsize*np, nPartDiff[np])
        for nt in range(len(trkL_sum)):
            bc = kinBins[nt]
            #p_TrkLVsE.Fill(bc, trkLength[nt])
            p_TrkLVsE_tot.Fill(bc, trkL_tot[nt]/100.)
            p_TrkLVsE_sum.Fill(bc, trkL_sum[nt]/100.)
            p_TrkLVsE_diff.Fill(bc, trkL_diff[nt]/100.)            

        continue

    # Plots right now are for electron only
    PDG = getTrkPDG(line)
    energy = getTrkE(line)
    if abs(PDG) != 11: continue

    # Fill Trk Length vs. Energy threshold
    trkL = getTrkZFinal(line) - getTrkZ(line)
    initial = TVector3(getTrkX(line), getTrkY(line), getTrkZ(line))
    final = TVector3(getTrkXFinal(line), getTrkYFinal(line), getTrkZFinal(line))
    for i in range(len(kinBins)):
        if energy > kinBins[i]:
            trkL_tot[i] += (final - initial).Mag()
            trkL_sum[i] += trkL
            if PDG == 11: trkL_diff[i] += trkL
            else:         trkL_diff[i] -= trkL
        else: break

    # Place some energy threshold cut
    if energy < threshold: continue

    # Fill for NPart plot
    rad_int = int(getTrkZ(line) / (stepsize*radLength))
    rad_fin = int(getTrkZFinal(line) / (stepsize*radLength))

    if rad_fin > npBins: rad_fin = npBins
    for i in range(rad_int, rad_fin):
        nPartSum[i] += 1
        if PDG == 11: nPartDiff[i] += 1
        else:         nPartDiff[i] -= 1

# end loop over file


# Clean up
output.Write()
output.Close()
infile.close()
        

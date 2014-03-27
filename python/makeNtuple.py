
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
# For the event by event analysis to calculate the E-field #
# it is a good idea to have an ntuple handy.  Otherwise    #
# I will have to make a ton of histograms (1 per event)    #
# which seems like it will be more work to deal with..     #
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#

import os
import sys
pwd = os.environ['PWD']
sys.path.insert(0, pwd+'/python')

from tools import *
from parseTrackFile import *
from ROOT import *
import numpy as n

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
f_part   = fname.split("_")[4]
f_npart  = (fname.split("_")[5]).split(".")[0]

infile = open(fname,"r")
output = TFile("rootfiles/TrkTree_"+f_nEvent+"_"+f_energy+"_"+f_mat+"_"+f_part+"_"+f_npart+".root","recreate")

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
## Setup the tree
#

event = n.zeros(1,dtype=int)
#x     = n.zeros(1,dtype=float)
#y     = n.zeros(1,dtype=float)
z     = n.zeros(1,dtype=float)
#xf    = n.zeros(1,dtype=float)
#yf    = n.zeros(1,dtype=float)
zf    = n.zeros(1,dtype=float)
E     = n.zeros(1,dtype=float)
pdg   = n.zeros(1,dtype=int)

tree = TTree("ntuple","Simple Tree")
tree.Branch('evtNum', event,'event/I')
tree.Branch('zi', z,'normal/D')
tree.Branch('zf', zf,'normal/D')
tree.Branch('E', E,'energy/D')
tree.Branch('pdg', pdg, 'pdg/I')

for line in infile:

    # Initialize stuff at beginning of event
    if "Event" in line:
        if event[0] % 50 ==0: 
            print "*** Processing Event ", event[0]
        continue

    # End of event. Fill profiles for event
    if "End" in line:
        event[0] += 1
        continue

    # Set Variables
    #x[0]  = getTrkX(line)
    #y[0]  = getTrkY(line)
    z[0]  = getTrkZ(line)
    #xf[0] = getTrkXFinal(line)
    #yf[0] = getTrkYFinal(line)
    zf[0] = getTrkZFinal(line)
    E[0]  = getTrkE(line)
    pdg[0] = getTrkPDG(line)

    tree.Fill()


# Clean up
infile.close()
output.Write("",TObject.kOverwrite)
output.Close()









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

#                                                                                                                               
## Specify input file                                                                                                           
#                                                                                                                               

if len(sys.argv) < 2:
    print "Must specify input file name"
    sys.exit()

# Parse input file                                                                                                             \
                                                                                                                                
fname = sys.argv[1]
f_energy = fname.split("_")[2]
f_nEvent = fname.split("_")[1]
f_mat    = fname.split("_")[3]
f_part   = (fname.split("_")[4]).split(".")[0]

infile = open(fname,"r")
output = TFile("rootfiles/TrkTree_"+f_nEvent+"_"+f_energy+"_"+f_mat+"_"+f_part+".root","recreate")

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

event = int(0)
x     = float(0)
y     = float(0)
z     = float(0)
xf    = float(0)
yf    = float(0)
zf    = float(0)
E     = float(0)


tree = TTree("ntuple","Simple Tree")
tree.Branch('evtNum', event,'normal/I')
tree.Branch('xi', hex(id(x)),'normal/F')
tree.Branch('yi', hex(id(y)),'normal/F')
tree.Branch('zi', hex(id(z)),'normal/F')
tree.Branch('xf', hex(id(xf)),'normal/F')
tree.Branch('yf', hex(id(yf)),'normal/F')
tree.Branch('zf', hex(id(zf)),'normal/F')
tree.Branch('E', hex(id(E)),'normal/F')

for line in infile:

    # Initialize stuff at beginning of evnet                                                                                    
    if "Event" in line:
        if event % 50 ==0: print "Event: ", event
        continue

    # End of event. Fill profiles for event                                                                                     
    if "End" in line:
        event += 1
        continue

    # Set Variables
    x  = getTrkX(line)
    y  = getTrkY(line)
    z  = getTrkZ(line)
    xf = getTrkXFinal(line)
    yf = getTrkYFinal(line)
    zf = getTrkZFinal(line)
    E  = getTrkE(line)
    
    tree.Fill()

# Clean up
infile.close()
output.Write()
output.Close()








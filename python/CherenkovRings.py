
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
# This script will read in the ouput from APotPlots for now and try to #
# show the cherenkov ring from each detector.  Right now it will have  #
# to be fed in grid parameters.                                        #
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#

from ROOT import *
import sys

#-------------------------------------#
# Read some input here
# Specify the input file as well
# as the antenna setup file
#-------------------------------------#

argv = sys.argv
if len(argv) != 3:
    print "Need to specify two arguments to run script"
    print "1.) input root file"
    print "2.) the antenna coordinates"
    sys.exit()

infname = argv[1]
antname = argv[2]
outname = "plots/" + (infname.split("/")[1]).split(".")[0] + "_XYPlot.root"

#-------------------------------------#
# Open file and read in antenna coords
#-------------------------------------#

#
## Input file
#

infile = TFile(infname)

#
## Get antenna coordinates
#

# Open file and some place holders
antfile = open(antname,"r")
antNum  = 0
antList = []

# Specify the min and max x-y coordinate
xmin   = 999
xmax   = -999
nxbins = 0
ymin   = 999
ymax   = -999
nybins = 0

# Method to create ant name
def makeAntName(num,x,y,z):
    name = "A_AntNum_"+str(num)
    name += "_pos_" + x + "_" + y + "_" + z
    #print name
    return name

# Some placeholders
prevX = 0

# Loop over input file
for line in antfile:

    # The Split Line
    sline = line.split()
    
    # Set coords
    x = sline[0]
    y = sline[1]
    z = sline[2]

    # Deal with floating point buisness
    if len(x.split(".")) > 1 and x.split(".")[1] == "0": x = x.split(".")[0]
    if len(y.split(".")) > 1 and y.split(".")[1] == "0": y = y.split(".")[0]
    if len(z.split(".")) > 1 and z.split(".")[1] == "0": z = z.split(".")[0]

    # The Antenna List
    antList.append(makeAntName(antNum,x,y,z))

    # Increment antenna number
    antNum += 1

    # Update the min and max info
    if float(x) < xmin: xmin = float(x)
    if float(x) > xmax: xmax = float(x)
    if float(y) < ymin: ymin = float(y)
    if float(y) > ymax: ymax = float(y)

    # Count nbins
    if prevX != float(x):
        nxbins += 1
        nybins += 1
        prevX = float(x)

# end loop over antenna file

#-------------------------------------#
# Now make the histogram and fill
#-------------------------------------#

# simple method to re-extract x-y coordinate
def getX(pname):
    return float(pname.split("_")[4])

def getY(pname):
    return float(pname.split("_")[5])

hist = TH2D("XYFixedZ","",nxbins,xmin,xmax,nybins,ymin,ymax)
hist.GetXaxis().SetTitle("x [m]")
hist.GetYaxis().SetTitle("y [m]")
hist.GetZaxis().SetTitle("A [Vm/s]")

for profName in antList:
    prof = infile.Get(profName)
    #print profName
    maximum = prof.GetMaximum()
    
    hist.Fill(getX(profName),getY(profName),maximum)

# end loop over profiles

outfile = TFile(outname,"recreate")
hist.Write()
outfile.Write()
outfile.Close()

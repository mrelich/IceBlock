
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
# This script is to generate a txt file to be read in by Geant4  #
# simulation for electric field. This particular script will     #
# generate a grid in x and y for a fixed z position. The goal is #
# be able to see the Chrenkov angle.                             #
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#

from math import sqrt, atan, pi

#-------------------------#
# Output file name
#-------------------------#

#outfile = open("antennaConfig/xz_40MeV.txt","w")
outfile = open("antennaConfig/xz_730cmAway_40MeV.txt","w")

#-------------------------#
# Fix x-y position
#-------------------------#

xpos = 7.3 # [m]
ypos = 0 # [m]

#-------------------------#
# Specify max for X and Y
#-------------------------#

zmin = 0
zmax = 10

#-------------------------#
# Specify step size
#-------------------------#

stepsize = 0.1 # [m]

#-------------------------#
# Now fill file
#-------------------------#

for z in range(int((zmax-zmin)/stepsize)+1):

    # Need to add some dummy information
    zpos = zmin + z*stepsize
    R = sqrt(xpos*xpos+ypos*ypos+zpos*zpos)
    theta = pi/2
    if zpos != 0:
        theta = atan(xpos/zpos)
    refAngle = -1
    phi = -1
    Rprime = R
    zprime = zpos

    # Write x
    outfile.write(str(xpos)+"\t")
    
    # Write y
    outfile.write(str(ypos)+"\t")
    
    # Write z
    outfile.write(str(zpos)+"\t")
    
    # dummy angles
    outfile.write(str(theta*180/pi)+"\t")
    outfile.write(str(refAngle*180/pi)+"\t")
    outfile.write(str(zprime)+"\n")
        

outfile.close()
    

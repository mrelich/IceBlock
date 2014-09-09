
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
# This script is to generate a txt file to be read in by Geant4  #
# simulation for electric field. This particular script will     #
# generate a grid in x and y for a fixed z position. The goal is #
# be able to see the Chrenkov angle.                             #
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#


#-------------------------#
# Output file name
#-------------------------#

outfile = open("antennaConfig/xz_40MeV.txt","w")

#-------------------------#
# Fix x-y position
#-------------------------#

xpos = 6 # [m]
ypos = 0 # [m]

#-------------------------#
# Specify max for X and Y
#-------------------------#

zmin = 0
zmax = 15

#-------------------------#
# Specify step size
#-------------------------#

stepsize = 0.2 # [m]

#-------------------------#
# Now fill file
#-------------------------#

for z in range(int((zmax-zmin)/stepsize)+1):
    
    # Write x
    outfile.write(str(xpos)+"\t")
    
    # Write y
    outfile.write(str(ypos)+"\t")
    
    # Write z
    outfile.write(str(zmin+z*stepsize)+"\n")
        

outfile.close()
    

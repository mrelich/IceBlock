
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
# Fix Z position
#-------------------------#

#zpos = 100 # [m]
#zpos = 5 # [m]
ypos = 0 # [m]

#-------------------------#
# Specify max for X and Y
#-------------------------#

maxi = 15
xmin = 0 #-maxi # [m]
xmax = maxi  # [m]
#ymin = -maxi # [m]
#ymax = maxi  # [m]
zmin = 0
zmax = 15

#-------------------------#
# Specify step size
#-------------------------#

stepsize = 0.5 # [m]

#-------------------------#
# Now fill file
#-------------------------#

for x in range(int((xmax-xmin)/stepsize)+1):
    #for y in range(int((ymax-ymin)/stepsize)+1):
    for z in range(int((zmax-zmin)/stepsize)+1):
        
        # Write x
        outfile.write(str(xmin+x*stepsize)+"\t")
        
        # Write y
        #outfile.write(str(ymin+y*stepsize)+"\t")
        outfile.write(str(ypos)+"\t")

        # Write z
        #outfile.write(str(zpos)+"\n")
        outfile.write(str(zmin+z*stepsize)+"\n")
        

outfile.close()
    

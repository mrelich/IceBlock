
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
# This script is to generate a txt file to be read in by Geant4  #
# simulation for electric field. This particular script will     #
# generate a grid in x and y for a fixed z position. The goal is #
# be able to see the Chrenkov angle.                             #
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#

import sys
from math import *

#-------------------------#
# Output file name
#-------------------------#

#outfile = open("antennaConfig/xzRefracted_40MeV.txt","w")
outfile = open("antennaConfig/temp_Rot30.txt","w")
#outfile = open("antennaConfig/temp_Rot60.txt","w")
#outfile = open("antennaConfig/xzUnRefracted_40MeV.txt","w")

#-------------------------#
# Fix x-y position
#-------------------------#

xpos = 6 # [m]
ypos = 0 # [m]

#-------------------------#
# Specify rotation angle
#-------------------------#

rotAngle = 30 * pi/180.
#rotAngle = 60 * pi/180.

#-------------------------#
# Specify angular steps
#-------------------------#

n_ice = 1.78
n_air = 1.0

criticalAng = asin(n_air/n_ice) * 180/pi # degrees
angleStep   = 1                          # degrees
#angleStep   = 0.5                          # degrees

#-------------------------#
# Define a few methods
#-------------------------#

# Get the incident angle given 
# the angle relative to beam
def getIncident(theta, rotAngle):
    return theta - rotAngle

# Get the refracted angle
def getRefracted(incident):
    return asin(n_ice/n_air * sin(incident))
    
# This is the angle between the un-refracted ray
# and the refracted ray
def getPhi(theta, rotAngle, refAngle):
    return refAngle - theta + rotAngle

# Get Rprime
def getRprime(z,R,x,phi):
    
    # Using quadratic method to solve
    #a = z*z - R*R*cos(phi)*cos(phi)
    #b = 2*R*R*R*cos(phi) - 2*R*z*z*cos(phi)
    #c = z*z*R*R - R*R*R*R
    a = z*z - R*R*cos(phi)*cos(phi)
    b = 2*R*x*x*cos(phi)
    c = -z*z*x*x - x*x*x*x

    # Deal with inside square brackets
    inSqrt = b*b-4*a*c
    if abs(inSqrt) < 1e-8: inSqrt = 0

    # Now get two solutions    
    R1 = (-b - sqrt(inSqrt)) / (2*a)
    R2 = (-b + sqrt(inSqrt)) / (2*a)

    # Take the smaller of the two angles!
    if R2 > 0 and R2 < R1: return R2    
    elif R1 > 0 and R1 <= R2: return R1

    else:
        print "Something is wrong, Rprime is"
        print "negative", R1, R2
        return 1000.

#-------------------------#
# Now fill file
#-------------------------#

# Start at whatever the rotation angle is.
# this will be normal to the surface!
theta    = rotAngle
incident = getIncident(theta,rotAngle)

# Hack to handle the case of when 
# the heights become negative
prevZ = 999999

while incident*180/pi < criticalAng:

    refAng = getRefracted(incident)
    phi = getPhi(theta,rotAngle,refAng)
    R = xpos / sin(theta)
    z = xpos / tan(theta)
    Rprime = getRprime(z,R,xpos, phi)
    zprime = sqrt(Rprime*Rprime - xpos*xpos)
    
    if prevZ - zprime < 0: 
        prevZ = zprime
        zprime = -zprime
    else:
        prevZ = zprime

    print theta * 180/pi, refAng*180/pi, Rprime/R
    
    # Write x
    outfile.write(str(xpos)+"\t")
    
    # Write y
    outfile.write(str(ypos)+"\t")
    
    # Write z
    outfile.write(str(z)+"\t")
        
    # Write angles
    outfile.write(str(theta*180/pi)+"\t")
    outfile.write(str(refAng*180/pi)+"\t")
    outfile.write(str(zprime)+"\n")
    

    # Increment angle
    theta += angleStep*pi/180
    incident = getIncident(theta,rotAngle)

outfile.close()
    

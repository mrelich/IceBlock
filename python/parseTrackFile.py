
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#
# The current output from Geant4 is for tracks and for steps. This file #
# will return the useful quantities related to the steps file when      #
# passed an input line from a file                                      #
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#

#-----------------------------#
# Get XPos
#-----------------------------#
def getTrkX(line):
    return float((line.split(" "))[0])

#-----------------------------#
# Get YPos
#-----------------------------#
def getTrkY(line):
    return float((line.split(" "))[1])

#-----------------------------#
# Get ZPos
#-----------------------------#
def getTrkZ(line):
    return float((line.split(" "))[2])

#-----------------------------#
# Get Kinetic Energy
#-----------------------------#
def getTrkE(line):
    return float((line.split(" "))[3])

#-----------------------------#
# Get TotalE
#-----------------------------#
def getTrkTotalE(line):
    return float((line.split(" "))[4])

#-----------------------------#
# Get XPos
#-----------------------------#
def getTrkXFinal(line):
    return float((line.split(" "))[5])

#-----------------------------#
# Get YPos
#-----------------------------#
def getTrkYFinal(line):
    return float((line.split(" "))[6])

#-----------------------------#
# Get ZPos
#-----------------------------#
def getTrkZFinal(line):
    return float((line.split(" "))[7])

#-----------------------------#
# Get TrkID
#-----------------------------#
def getTrkID(line):
    return int((line.split(" "))[8])

#-----------------------------#
# Get ParentID
#-----------------------------#
def getTrkParentID(line):
    return int((line.split(" "))[9])

#-----------------------------#
# Get PDG
#-----------------------------#
def getTrkPDG(line):
    return int((line.split(" "))[10])


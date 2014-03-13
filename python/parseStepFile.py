
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#
# The current output from Geant4 is for tracks and for steps. This file #
# will return the useful quantities related to the steps file when      #
# passed an input line from a file                                      #
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#

#-----------------------------#
# Get XPos
#-----------------------------#
def getStepX(line):
    return float((line.split(" "))[0])

#-----------------------------#
# Get YPos
#-----------------------------#
def getStepY(line):
    return float((line.split(" "))[1])

#-----------------------------#
# Get ZPos
#-----------------------------#
def getStepZ(line):
    return float((line.split(" "))[2])

#-----------------------------#
# Get dX
#-----------------------------#
def getStepdX(line):
    return float((line.split(" "))[3])

#-----------------------------#
# Get Kinetic Energy
#-----------------------------#
def getStepE(line):
    return float((line.split(" "))[4])

#-----------------------------#
# Get dE
#-----------------------------#
def getStepdE(line):
    return float((line.split(" "))[5])

#-----------------------------#
# Get PDG
#-----------------------------#
def getStepPDG(line):
    return int((line.split(" "))[6])

#-----------------------------#
# Get trkID
#-----------------------------#
def getStepTrkID(line):
    return int((line.split(" "))[7])

#-----------------------------#
# Get ParentID
#-----------------------------#
def getStepParentID(line):
    return int((line.split(" "))[8])

#-----------------------------#
# Get ELoss ProcID
#-----------------------------#
def getStepProcID(line):
    return int((line.split(" "))[9])


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#
# Script to analyze the photon distributions in the ice #
# in order to estimate the effects of radiation.        #
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#

from parseTrackFile import *
from tools import *
from ROOT import *

#-----------------------------------------#
# Input file
#-----------------------------------------#

#f_in = open("tracks/output_50_40_ice_eBeam_np10000_AngleScane_100m_singlePos.dat","r")
#f_in = open("tracks/output_50_40_ice_eBeam_np1_antFileUsed.dat","r")
f_in = open("tracks/output_10000_40_ice_eBeam_np1_HardCodedAntenna_R7m_singlePos.dat","r")
#nparticles = 10000.
#nbunches   = 5.
#scale = 1.e9 / nbunches / nparticles
scale = 1

#-----------------------------------------#
# Open output file
#-----------------------------------------#

f_out = TFile("rootfiles/tracks/test.root","recreate")

#-----------------------------------------#
# Define histograms
#-----------------------------------------#

nbins = 100
emin  = 0
emax  = 40
xtitle = "Energy [MeV]"
ytitle = "Entries/bin"

h_E = makeHist("g_E",nbins,emin,emax,kBlack,20,xtitle,ytitle)

# Total events staying inside ice
h_Ein = makeHist("tempIn",nbins,emin,emax,kBlack,20,xtitle,ytitle)
p_Ein = makeProfile("p_Ein",nbins,emin,emax,kBlack,20,xtitle,"Photons / primary")

# Total number of events leaving
h_Eleave = makeHist("tempLeave",nbins,emin,emax,kBlack,20,xtitle,ytitle)
p_Eleave = makeProfile("p_Eleave",nbins,emin,emax,kBlack,20,xtitle,"Photons / primary")

# Events leaving through top
h_Etop = makeHist("tempTop",nbins,emin,emax,kBlack,20,xtitle,ytitle)
p_Etop = makeProfile("p_Etop",nbins,emin,emax,kBlack,20,xtitle,"Photons / primary")

# Events leaving through top
h_Eside = makeHist("tempSide",nbins,emin,emax,kBlack,20,xtitle,ytitle)
p_Eside = makeProfile("p_Eside",nbins,emin,emax,kBlack,20,xtitle,"Photons / primary")


#-----------------------------------------#
# Method to determine if the particle 
# starts in the ice
#-----------------------------------------#
def inIce(x,y,z,side):
    
    # Simplicity, assume it is cube with
    # sides 30 cm. This will be conservative
    # estimate.  The base of the box is at (0,0,0)
    if not (-side/2 < x and x < side/2):
        return False
    if not (-side/2 < y and y < side/2):
        return False
    if not (0 < z and z < side):
        return False
    
    return True



#-----------------------------------------#
# Loop and fill histos
#-----------------------------------------#
counter = 0
for line in f_in:

    # Break out of the loop if counter greater than number
    #if counter > nbunches: break

    # Skip lines that have event number in them
    if "Event" in line:
        counter += 1
        continue

    # Skip Lines that have end in them
    if "End" in line:
        for ibin in range(1,nbins+1):
            con = h_Ein.GetBinContent(ibin)
            cen = h_Ein.GetBinCenter(ibin)
            p_Ein.Fill(cen,con)

            con = h_Eleave.GetBinContent(ibin)
            cen = h_Eleave.GetBinCenter(ibin)
            p_Eleave.Fill(cen,con)

            con = h_Etop.GetBinContent(ibin)
            cen = h_Etop.GetBinCenter(ibin)
            p_Etop.Fill(cen,con)

            con = h_Eside.GetBinContent(ibin)
            cen = h_Eside.GetBinCenter(ibin)
            p_Eside.Fill(cen,con)


        h_Ein.Reset()
        h_Eleave.Reset()
        h_Etop.Reset()
        h_Eside.Reset()

        continue


    # Parse the input
    pdg    = getTrkPDG(line)
    energy = getTrkE(line)

    # Only look at photons
    if pdg != 22: continue

    # Get the initial points
    xi = getTrkX(line)
    yi = getTrkY(line)
    zi = getTrkZ(line)

    # Get final points
    xf = getTrkXFinal(line)
    yf = getTrkYFinal(line)
    zf = getTrkZFinal(line)


    # Fill if the final points are inside the ice
    if inIce(xf,yf,zf,50):
        h_E.Fill(energy,scale)
        h_Ein.Fill(energy,scale)

    # Fill if the initial point is iniside, but
    # the final point is outside, implying particle
    # leaves
    if inIce(xi,yi,zi,30) and not inIce(xf,yf,zf,30):
        h_Eleave.Fill(energy,scale)

        # See if the event leaves through the top of the ice
        if inIce(xi,yi,zi,30) and not inIce(xi,yi,zf,30):
            h_Etop.Fill(energy,scale)

    
        elif inIce(xi,yi,zi,30) and not inIce(xf,yf,zi,30):
            h_Eside.Fill(energy,scale)



h_Ein.Delete()
h_Eleave.Delete()
h_Etop.Delete()
h_Eside.Delete()
f_out.Write()
f_out.Close()

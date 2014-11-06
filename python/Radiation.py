
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#
# Script to analyze the photon distributions in the ice #
# in order to estimate the effects of radiation.        #
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-#

from parseTrackFile import *
from tools import *
from ROOT import *
from math import pi

#-----------------------------------------#
# Input file
#-----------------------------------------#

#f_in = open("tracks/output_50_40_ice_eBeam_np10000_AngleScane_100m_singlePos.dat","r")
#f_in = open("tracks/output_50_40_ice_eBeam_np1_antFileUsed.dat","r")
#f_in = open("tracks/output_10000_40_ice_eBeam_np1_HardCodedAntenna_R7m_singlePos.dat","r")
f_in = open("tracks/output_10000_40_ice_eBeam_np1_HardCodedAntenna_R100m_singlePos_shortenedZ.dat","r")
#f_in = open("tracks/output_100000_40_ice_eBeam_np1_HardCodedAntenna_R100m_singlePos_shortenedZ.dat","r")
#nparticles = 10000.
#nbunches   = 5.
#scale = 1.e9 / nbunches / nparticles
scale = 1

#-----------------------------------------#
# Open output file
#-----------------------------------------#

f_out = TFile("rootfiles/tracks/test.root","recreate")
#f_out = TFile("rootfiles/tracks/testOld.root","recreate")
#f_out = TFile("rootfiles/tracks/Radiation_40MeV_10kEvents.root","recreate")
#f_out = TFile("rootfiles/tracks/Radiation_40MeV_100kEvents.root","recreate")

#-----------------------------------------#
# Define histograms
#-----------------------------------------#

generic = "Photons / primary"

# Energy bins
nbins = 100
emin  = 0
emax  = 40
xtitle = "Energy [MeV]"
ytitle = "Entries/bin"

h_E = makeHist("g_E",nbins,emin,emax,kBlack,20,xtitle,ytitle)

# Total events staying inside ice
h_Ein = makeHist("tempIn",nbins,emin,emax,kBlack,20,xtitle,ytitle)
p_Ein = makeProfile("p_Ein",nbins,emin,emax,kBlack,20,xtitle,generic)

# Total number of events leaving ice
h_Eleave = makeHist("tempLeave",nbins,emin,emax,kBlack,20,xtitle,ytitle)
p_Eleave = makeProfile("p_Eleave",nbins,emin,emax,kBlack,20,xtitle,generic)

# Events leaving through top ice
h_Etop = makeHist("tempTop",nbins,emin,emax,kBlack,20,xtitle,ytitle)
p_Etop = makeProfile("p_Etop",nbins,emin,emax,kBlack,20,xtitle,generic)

# Events leaving through top ice
h_Eside = makeHist("tempSide",nbins,emin,emax,kBlack,20,xtitle,ytitle)
p_Eside = makeProfile("p_Eside",nbins,emin,emax,kBlack,20,xtitle,generic)

# All Photons outside the ice
h_Eout = makeHist("tempOut",nbins,emin,emax,kBlack,20,xtitle,ytitle)
p_Eout = makeProfile("p_Eout",nbins,emin,emax,kBlack,20,xtitle,generic)

# Photons made outside the ice
h_EMadeOut = makeHist("tempProducedOut",nbins,emin,emax,kBlack,20,xtitle,ytitle)
p_EMadeOut = makeProfile("p_EMadeOut",nbins,emin,emax,kBlack,20,xtitle,generic)

# Angular bins
nabins = 90
amin   = 0
amax   = 180
xtitle = "Angle Relative to beam [deg]"
ytitle = "Entries/bin"

# Photons angular distribution that leave ice
h_AngleOut = makeHist("tempAngleOut",nabins,amin,amax,kBlack,20,xtitle,ytitle)
p_AngleOut = makeProfile("p_AngleOut",nabins,amin,amax,kBlack,20,xtitle,generic)

# Plot Energy vs angle
h_EvsAngOut = makeHist2("tempEvsAngOut",nbins,emin,emax,
                        nabins,amin,amax,"","","")
p_EvsAngOut = makeProf2("p_EvsAngOut",nbins,emin,emax,
                        nabins,amin,amax,"Energy [MeV]",
                        "Angle Relative to beam [deg]",
                        generic)

# Height bins
nzbins = 100
zmin   = -5  # [m]
zmax   = 20  # [m]

# Final height
h_zfOut = makeHist("tempzfOut",nzbins,zmin,zmax,kBlack,20,"","")
p_zfOut = makeProfile("p_zfOut",nzbins,zmin,zmax,kBlack,20,"z [m]",generic)

# Number of paritlces to pass through a cross-section
# of a human
ndbins = 18
dmin   = 2
dmax   = 20
h_cross = makeHist("tempCross",ndbins,dmin,dmax,kBlack,20,"","")
p_cross = makeProfile("p_cross",ndbins,dmin,dmax,kBlack,20,"distance [m]", generic)

h_EvsCross = makeHist2("tempEvsCross",nbins,emin,emax,ndbins,dmin,dmax,"","","")
p_EvsCross = makeProf2("p_EvsCross",nbins,emin,emax,
                       ndbins,dmin,dmax,
                       "Energy [Mev]",
                       "Distance [m]",
                       generic)


#-----------------------------------------#
# Method to determine if the particle 
# starts in the ice
#-----------------------------------------#
def inIce(x,y,z):
    
    # The box is at (0,0,0)
    sidex = 100.0 # cm
    sidey = 30.0 # cm
    sidez = 30.0 # cm
    if not (-sidex/2 < x and x < sidex/2):
        return False
    if not (-sidey/2 < y and y < sidey/2):
        return False
    if not (0 < z and z < sidez):
        return False
    
    return True


#-----------------------------------------#
# Get angle
#-----------------------------------------#
def getAngle(x,y,z):
    unit = TVector3(0,0,1)
    vect = TVector3(x,y,z)
    return unit.Angle(vect) * 180./pi

#-----------------------------------------#
# Method to fill the profile
#-----------------------------------------#
def fillProf(hist,prof):
    nbins = hist.GetNbinsX()
    for ibin in range(1,nbins+1):
        con = hist.GetBinContent(ibin)
        cen = hist.GetBinCenter(ibin)
        prof.Fill(cen,con)
    hist.Reset()

#-----------------------------------------#
# Method to fill 2D profile
#-----------------------------------------#
def fillProf2(hist,prof):
    nbinsx = hist.GetXaxis().GetNbins()
    nbinsy = hist.GetYaxis().GetNbins()
    for xbin in range(1,nbinsx+1):
        for ybin in range(1,nbinsy+1):
            con  = hist.GetBinContent(xbin,ybin)
            cenx = hist.GetXaxis().GetBinCenter(xbin)
            ceny = hist.GetYaxis().GetBinCenter(ybin)
            prof.Fill(cenx,ceny,con)
    hist.Reset()

#-----------------------------------------#
# Define information relevant for figuring
# out how many photons pass through some 
# cross-section
#-----------------------------------------#
def passThrough(d_obs, xi, xf):

    # Specify height and width
    # of generic person
    height = 1.8 # m
    width  = 0.5 # m

    # Define plane parameters
    normal = TVector3(1,0,0)
    p0     = TVector3(d_obs,0,0)
    
    # Setup particle params. Add 5m
    # to the z position to account for
    # the Ice being raised up
    x0     = TVector3(xi[0],xi[1],xi[2])
    v      = TVector3(xf[0]-xi[0],
                      xf[1]-xi[1],
                      xf[2]-xi[2] + 5)
    
    # If v.Dot(normal) = 0 then lines parallel
    # Handle cases
    if v.Dot(normal) == 0:
        if normal.Dot(p0-x0) == 0: 
            print "Line inside plane"
            return True # line inside plane
        else:
            #print "Lines never cross"
            return False # Never cross
    
    # Solve for t in vector equation of points
    t = normal.Dot(p0-x0)/normal.Dot(v)
    #print "t: ", t

    # Now that we have t, find if the z and y 
    # positions are within the cross-section
    # of the human
    y_obs = x0.Y() + t * v.Y()
    z_obs = x0.Z() + t * v.Z()
    if not (-width/2. < y_obs and y_obs < width/2.):
        #print "Not in width: ", y_obs
        return False
    if not (0 < z_obs and z_obs < height):
        #print "Not in height", z_obs
        return False

    # Otherwise we are good
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
        fillProf(h_Ein, p_Ein)
        fillProf(h_Eleave, p_Eleave)
        fillProf(h_Etop, p_Etop)
        fillProf(h_Eside, p_Eside)
        fillProf(h_Eout, p_Eout)
        fillProf(h_EMadeOut, p_EMadeOut)
        fillProf(h_AngleOut, p_AngleOut)
        fillProf2(h_EvsAngOut, p_EvsAngOut)
        fillProf(h_zfOut, p_zfOut)
        fillProf(h_cross, p_cross)
        fillProf2(h_EvsCross, p_EvsCross)
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

    # Get angle relative to beam
    angle = getAngle(xf-xi,yf-yi,zf-zi)

    # Fill if the final points are inside the ice
    if inIce(xf,yf,zf):
        h_E.Fill(energy,scale)
        h_Ein.Fill(energy,scale)

    # Fill if the initial point is iniside, but
    # the final point is outside, implying particle
    # leaves
    if inIce(xi,yi,zi) and not inIce(xf,yf,zf):
        h_Eleave.Fill(energy,scale)

        # See if the event leaves through the top of the ice
        if inIce(xi,yi,zi) and not inIce(xi,yi,zf):
            h_Etop.Fill(energy,scale)

    
        elif inIce(xi,yi,zi) and not inIce(xf,yf,zi):
            h_Eside.Fill(energy,scale)
    
    if not inIce(xf,yf,zf):
        h_Eout.Fill(energy,scale)
        h_AngleOut.Fill(angle,scale)
        h_EvsAngOut.Fill(energy,angle,scale)
        h_zfOut.Fill(zf/100.,scale)
        
        # Determine if the point crosses for several
        # distances.
        v_xi = TVector3(xi/100,yi/100,zi/100)
        v_xf = TVector3(xf/100,yf/100,zf/100)
        for ibin in range(1,h_cross.GetNbinsX()+1):
            bc = h_cross.GetXaxis().GetBinLowEdge(ibin)
            if passThrough(bc, v_xi, v_xf):
                #print "found one: ", bc
                h_cross.Fill(bc,scale)
                h_EvsCross.Fill(energy,bc,scale)
                
    if not inIce(xi,yi,zi) and not inIce(xf,yf,zf):
        h_EMadeOut.Fill(energy,scale)

h_Ein.Delete()
h_Eleave.Delete()
h_Etop.Delete()
h_Eside.Delete()
h_Eout.Delete()
h_EMadeOut.Delete()
h_AngleOut.Delete()
h_EvsAngOut.Delete()
h_zfOut.Delete()
h_cross.Delete()
h_EvsCross.Delete()
f_out.Write()
f_out.Close()


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=#
# Helper class to make histograms and things. #
# Should make plotting easier for future      #
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=#

from ROOT import *

#---------------------------#
# Make Profile
#---------------------------#
def makeProfile(name, nxbins, xmin, xmax, color, marker,
                xtitle, ytitle):
    prof = TProfile(name,name,nxbins,xmin,xmax)
    prof.SetLineColor(color)
    prof.SetMarkerColor(color)
    prof.SetMarkerStyle(marker)
    prof.GetXaxis().SetTitle(xtitle)
    prof.GetYaxis().SetTitle(ytitle)
    prof.SetMarkerSize(0.5);
    prof.Sumw2()
    return prof

def makeProfileVar(name, xbins, color, marker,
                   xtitle, ytitle):
    prof = TProfile(name,name,len(xbins)-1, xbins)
    prof.SetLineColor(color)
    prof.SetMarkerColor(color)
    prof.SetMarkerStyle(marker)
    prof.GetXaxis().SetTitle(xtitle)
    prof.GetYaxis().SetTitle(ytitle)
    prof.Sumw2()
    return prof



#---------------------------#
# Make Histogram
#---------------------------#
def makeHist(name, nxbins, xmin, xmax, color, marker,
             xtitle, ytitle):
    hist = TH1F(name,name,nxbins,xmin,xmax)
    hist.SetLineColor(color)
    hist.SetMarkerColor(color)
    hist.SetMarkerStyle(marker)
    hist.GetXaxis().SetTitle(xtitle)
    hist.GetYaxis().SetTitle(ytitle)
    hist.Sumw2()
    return hist


#---------------------------#
# Make Canvas
#---------------------------#
def makeCanvas(name):
    c = TCanvas(name,name,700,600)
    c.SetTopMargin(0.05)
    c.SetBottomMargin(0.12)
    c.SetRightMargin(0.05)
    c.SetLeftMargin(0.12)
    return c

#---------------------------#
# Make Legend
#---------------------------#
def makeLegend(x0,x1,y0,y1):
    leg = TLegend(x0, y0, x1, y1);
    leg.SetFillColor(0);
    leg.SetFillStyle(0);
    leg.SetLineColor(0);
    leg.SetBorderSize(0);
    leg.SetTextSize(0.04);
    return leg;

#---------------------------#
# Set attributes
#---------------------------#
def setAtt(hist, color, marker):
    hist.SetLineColor(color)
    hist.SetMarkerColor(color)
    hist.SetMarkerStyle(marker)
    hist.SetMarkerSize(0.75)
    return hist


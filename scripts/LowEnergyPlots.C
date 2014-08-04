
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
// So far my focus has been on confirming my understanding and simulation //
// of Askaryan effect.  However for ELS test beam we will be using 40 MeV //
// electrons which have a much different shower profile. This script will //
// format the figures from LowEnergyShower.py                             //
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//

#include "myHist.C"

TString pEnergy = "40MeV";
TString nPrimary = "1000";
TString savedir = "plots/LowEnergyShower/";

//-------------------------------------------//
// Main
//-------------------------------------------//
void LowEnergyPlots(int opt)
{

  // Input File
  //TFile* infile = new TFile("rootfiles/LowEnergyStudy_40_ice_eBeam_np1.root");
  //TFile* infile = new TFile("rootfiles/LowEnergyStudy_40_ice_eBeam_np100.root");
  TFile* infile = new TFile("rootfiles/LowEnergyStudy_40_ice_eBeam_np1000.root");
  
  // Options to run
  if(opt == 0) plotNPart(infile);
  if(opt == 1) plotKin(infile);
  if(opt == 2){
    plotEvsZI(infile,false);
    plotEvsZI(infile,true);
  }
}

//-------------------------------------------//
// Plot NPart
//-------------------------------------------//
void plotNPart(TFile* file)
{

  // Make Canvas
  TCanvas* c = makeCanvas("c");
  c->SetLogy();

  // Make Legend
  TLegend* leg = makeLegend(0.55,0.75,0.7,0.9);
  leg->SetHeader(("E^{primary}_{e^{-}} = "+pEnergy+", N^{primary}="+nPrimary).Data());

  // Load histograms
  TString xtit = "Depth [cm]";
  TString ytit = "Number of Particles";
  TH1F* h_NPartSum  = getHist(file,"NPartSum",xtit,ytit,kBlack,20);
  TH1F* h_NPartDiff = getHist(file,"NPartDiff",xtit,ytit,kBlue,25);

  // Add to Legend
  leg->AddEntry(h_NPartSum,"N(e^{-}+e^{+})","lep");
  leg->AddEntry(h_NPartDiff,"N(e^{-}-e^{+})","lep");

  // Now Draw
  h_NPartSum->SetMaximum(1.2*h_NPartSum->GetMaximum());
  h_NPartSum->Draw();
  h_NPartDiff->Draw("same");
  leg->Draw("same");

  // Save
  c->SaveAs((savedir+"NPart_NPrimary"+nPrimary+".png").Data());
}

//-------------------------------------------//
// Plot kinematic distributions
//-------------------------------------------//
void plotKin(TFile* file)
{

  // Make Canvas
  TCanvas* c = makeCanvas("c");
  c->SetLogy();

  // Make Legend
  TLegend* leg = makeLegend(0.55,0.75,0.83,0.9);
  leg->SetHeader(("E^{primary}_{e^{-}} = "+pEnergy+", N^{primary}="+nPrimary).Data());

  // Specify the plots
  vector<TString> pnames;
  pnames.push_back("initialZ");
  pnames.push_back("finalZ");
  pnames.push_back("dZ");
  pnames.push_back("E");

  // Specify xtitles
  vector<TString> xtitles;
  xtitles.push_back("z_{initial} [cm]");
  xtitles.push_back("z_{final} [cm]");
  xtitles.push_back("z_{final} - z_{initial} [cm]");
  xtitles.push_back("E_{initial} [MeV]");

  // Specify ytitle
  TString ytitle = "Entries/bin";
  
  // Loop and get histograms
  TH1F* h = NULL;
  for(unsigned int i=0; i<pnames.size(); ++i){
    TString pname = pnames.at(i);
    TString xtitle = xtitles.at(i);
    
    h = getHist(file,pname,xtitle,ytitle,kBlack,20);
    h->Draw("hist");
    leg->Draw("same");
    c->SaveAs((savedir+pname+"_NPrimary"+nPrimary+"_kin.png").Data());
  }

  //delete c;
  //delete h;


}

//-------------------------------------------//
// Plot E vs. ZI, which is 2D histogram
//-------------------------------------------//
void plotEvsZI(TFile* file, bool xylog)
{

  // Legend for show
  TLegend* leg = makeLegend(0.15,0.25,0.87,0.9);
  leg->SetHeader(("E^{primary}_{e^{-}} = "+pEnergy+", N^{primary}="+nPrimary).Data());

  // Make Canvas
  TCanvas* c = makeCanvas("c");
  c->SetLogz();
  c->SetRightMargin(0.12);
  if(xylog){
    c->SetLogy();
    c->SetLogx();
  }

  // Now get figure
  TH2F* h = getHist2(file,"EvsZI","z_{initial} [cm]",
		     "E_{initial} [MeV]", "Entries");

  // Draw
  h->Draw("colz");
  leg->Draw("same");

  // Save
  c->SaveAs((savedir+"EvsZI_NPrimary" + nPrimary + 
	     (xylog ? "_logxy" : "") + ".png").Data());

}

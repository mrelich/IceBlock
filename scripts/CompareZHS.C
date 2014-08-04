
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-//
// This script will compare the results for ZHS and G4 from  //
// as many shower energies as I can construct.               //
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-//

#include "myHist.C"

typedef unsigned int uint;

TString savedir = "plots/ZHSComp/";

//---------------------------------------//
// Main
//---------------------------------------//
void CompareZHS()
{

  // Specify the directories for the input files
  //TString zhsdir = "../ZHS_lowEnergy/rootfiles/";
  //TString zhsdir = "../ZHS_100GeV/rootfiles/";
  TString zhsdir = "../ZHS_1TeV/rootfiles/";
  TString g4dir  = "rootfiles/";

  // Now specify the paths to each file
  vector<TString> g4fnames;
  //g4fnames.push_back(g4dir+"TrkAna_100_100000_ice_eBeam_thresh1MeV_1MeV_larger.root");
  g4fnames.push_back(g4dir+"TrkAna_50_1000000_ice_eBeam_thresh1MeV_1MeV_larger.root");
  //g4fnames.push_back(g4dir+"TrkAna_20_10000000_ice_eBeam_thresh1MeV_1MeV_larger.root");

  // Specify paths to ZHS files
  vector<TString> zhsfnames;
  //zhsfnames.push_back(zhsdir+"beam100e3MeV_NPart.root");
  zhsfnames.push_back(zhsdir+"beam1e6MeV_NPart.root");
  //zhsfnames.push_back(zhsdir+"beam10e6MeV_NPart.root");


  // Labels for Energy
  vector<TString> showerEs;
  //showerEs.push_back("100GeV");
  showerEs.push_back("1TeV");
  //showerEs.push_back("10TeV");

  // Loop and plot
  for(uint i=0; i<g4fnames.size(); ++i)
    plot(g4fnames.at(i),zhsfnames.at(i),showerEs.at(i));
  
}

//---------------------------------------//
// Plot
//---------------------------------------//
void plot(TString g4fname, TString zhsfname, TString showerE)
{

  // Make Canvas
  TCanvas* c = makeCanvas("c");

  // Put files in vector
  vector<TString> fnames;
  fnames.push_back(g4fname);
  fnames.push_back(zhsfname);

  // Leg names
  vector<TString> names;
  names.push_back("GEANT4");
  names.push_back("ZHS");

  // Make legend
  TLegend* leg = makeLegend(0.7,0.9,0.75,0.92);
  leg->SetHeader(("E_{e^{-}} = " + showerE).Data());

  // Decorations
  int colors[]  = {kBlack, kBlue};
  int markers[] = {20, 25};
  TString xtitle = "Shower Depth [#chi_{0}]";
  TString ytitle = "N(e^{-}-e^{+})";
  
  // Loop and load histograms
  TH1* hists[2] = {NULL, NULL};
  for(uint f=0; f<fnames.size(); ++f){
    TFile* file = new TFile(fnames.at(f).Data());
    hists[f] = getHist(file,"NPartDiff",xtitle,ytitle,
		       colors[f],markers[f]);
    hists[f]->SetMarkerSize(1.);
    hists[f]->SetDirectory(0);
    leg->AddEntry(hists[f],names.at(f).Data(),"lep");
    delete file;
  }

  // Get maximum
  double maximum = getMax(hists[0],hists[1]);
  hists[0]->SetMaximum(maximum*1.2);
  
  // Now plot
  hists[0]->Draw();
  hists[1]->Draw("same");
  leg->Draw("same");

  c->SaveAs((savedir+showerE+".png").Data());

}

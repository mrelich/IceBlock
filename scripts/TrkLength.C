
#include "myHist.C"

int colors[] = {kBlack, kRed, kBlue};
int markers[] = {20,25,23};

//------------------------------------//
// Main
//------------------------------------//
void TrkLength()
{

  // Track length vs E Threshold
  //trkLVsE();

  // Track length vs. beam energy
  trkLVsBeamE();

}

//------------------------------------//
// Plot Track length vs. Energy 
// threshold
//------------------------------------//
void trkLVsE()
{
  // Leave the input file as only one file now. 
  // Later I can make it an option if needed
  TFile* file = new TFile("rootfiles/TrkAna_50_100000_ice.root");

  // It is less code to loop
  vector<TString> pnames;
  pnames.push_back("TrkLVsE_tot");
  pnames.push_back("TrkLVsE_sum");
  pnames.push_back("TrkLVsE_diff");

  // Legend Names
  vector<TString> names;
  names.push_back("Absolute"); 
  names.push_back("Projected (e+p)"); 
  names.push_back("Projected (e-p)");

  // Make Canvas
  TCanvas* c = makeCanvas("c");
  c->SetLogy();
  c->SetLogx();

  // Make legend
  TLegend* leg = makeLegend(0.2,0.4,0.3,0.5);
  
  // Get Dummy histogram
  TH1F* dummy = getDummy("Energy Threshold [MeV]",
			 "Total Track Length [m]",
			 0.1, 10, 1, 1000);
  dummy->Draw();

  // Loop 
  TProfile* prof = NULL;
  TGraph* graph = NULL;
  for(unsigned int i=0; i<pnames.size(); ++i){
    
    // Get Prof
    TString pname = pnames.at(i);
    prof = getProfile(file,pname,"","",colors[i],markers[i]);
    
    // Get graph and draw
    graph = makeGraph(prof);
    setAtt(graph,colors[i],markers[i],"","");
    graph->Draw("epsame");
    
    // Add to legend
    leg->AddEntry(graph, names.at(i).Data(), "ep");
    
  }
  
  leg->Draw("same");
  
  c->SaveAs("plots/TrkLength/energy_100000_ice.png");

}

//------------------------------------//
// Plot Track length vs. Energy 
// threshold
//------------------------------------//
void trkLVsBeamE()
{

  float cutval = 1.0; // MeV

  // Leave the input file as only one file now. 
  // Later I can make it an option if needed
  TString fbase = "rootfiles/TrkAna_";
  const int nFiles = 4;
  TFile* files[nFiles];
  files[0] = new TFile((fbase+"100_100000_ice_eBeam.dat.root").Data());
  files[1] = new TFile((fbase+"100_200000_ice_eBeam.dat.root").Data());
  files[2] = new TFile((fbase+"100_400000_ice_eBeam.dat.root").Data());
  //files[3] = new TFile((fbase+"100_700000_ice_eBeam.dat.root").Data());
  files[3] = new TFile((fbase+"20_1000000_ice_eBeam.dat.root").Data());

  // It is less code to loop
  vector<TString> pnames;
  pnames.push_back("TrkLVsE_tot");
  pnames.push_back("TrkLVsE_sum");
  pnames.push_back("TrkLVsE_diff");

  // Legend Names
  vector<TString> names;
  names.push_back("Absolute"); 
  names.push_back("Projected (e+p)"); 
  names.push_back("Projected (e-p)");

  // Make Canvas
  TCanvas* c = makeCanvas("c");
  c->SetLogy();
  c->SetLogx();

  // Make legend
  TLegend* leg = makeLegend(0.6,0.8,0.2,0.4);
  
  // Get Dummy histogram
  TH1F* dummy = getDummy("Shower Energy [GeV]",
			 "Total Track Length [m]",
			 100, 1000, 10, 10000);

  dummy->GetXaxis()->SetMoreLogLabels();
  dummy->Draw();

  // Loop 
  TProfile* prof = NULL;
  TGraph* graph = NULL;
  TFile* file = NULL;

  for(unsigned int i=0; i<pnames.size(); ++i){
    TString pname = pnames.at(i);
    //float x[nFiles] = {100,200,400,700,1000};
    float x[nFiles] = {100,200,400,1000};
    float y[nFiles];
    int np = 0;
    
    for(int f=0; f<nFiles; ++f){
      file = files[f];
      
      // Get Prof
      prof = getProfile(file,pname,"","",colors[i],markers[i]);
      int bin = prof->FindBin(cutval);
      y[f] = prof->GetBinContent(bin);      
      np++;
      
    }

    graph = new TGraph(np, x, y);
    setAtt(graph,colors[i],markers[i],"","");
    graph->Draw("epsame");
    leg->AddEntry(graph, names.at(i).Data(), "ep");
  }

  leg->Draw("same");
  
  c->SaveAs("plots/TrkLength/energy_vsTrkLength_ice.png");
  
}

//------------------------------------------------//
// Make graph from Profile
//------------------------------------------------//
TGraph* makeGraph(TProfile* prof)
{
  
  float x[100];
  float y[100];
  
  int nbins = prof->GetNbinsX();
  for(int bin=0; bin<nbins; ++bin){
    x[bin] = prof->GetBinLowEdge(bin+1);
    y[bin] = prof->GetBinContent(bin+1);
  }
  
  TGraph* gr = new TGraph(nbins, x, y);
  return gr;

}

//------------------------------------------------//
// Set TGraph properties
//------------------------------------------------//
void setAtt(TGraph* &gr, int color, int marker, TString xtitle, 
            TString ytitle){
  gr->SetLineColor(color);
  gr->SetMarkerColor(color);
  gr->SetMarkerStyle(marker);
  gr->GetXaxis()->SetTitle(xtitle.Data());
  gr->GetYaxis()->SetTitle(ytitle.Data());
  
}

//------------------------------------------------//
// Get Dummy histogram
//------------------------------------------------//
TH1F* getDummy(TString xtitle, TString ytitle,
	       float xmin, float xmax,
	       float minimum, float maximum)
{
  
  TH1F* h = new TH1F("h","h",10,xmin,xmax);
  h->SetMaximum(maximum);
  h->SetMinimum(minimum);
  h->GetXaxis()->SetTitle(xtitle.Data());
  h->GetYaxis()->SetTitle(ytitle.Data());
  h->SetStats(0);
  h->SetTitle("");
   
   return h;
   
}

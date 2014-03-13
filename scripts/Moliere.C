
#include "scripts/myHist.C"
#include "TGraph.h"
#include "TString.h"
#include <sstream>

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
// Simple script to draw the results for ice, //
// iron, and lead for the Moliere Radius. In  //
// each case the Moliere Radius will be drawn //
// on the figure as well.                     //
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//

const int nFiles = 3;
TString materials[nFiles] = {"ice", "iron", "lead"};
int colors[nFiles]        = {kBlack, kBlue, kRed};
int markers[nFiles]       = {20, 25, 22};
TString energy = "";
TString eLbl   = "";

//------------------------------------------------//
// Main
//------------------------------------------------//
void Moliere(int opt)
{

  // Plot name
  TString pname = "rVsShowerE";

  // Base file name and path
  TString fbase = "rootfiles/StepAna_100_";

  // Energies supported
  if(opt == 0){
    energy = "1000";
    eLbl = "1 GeV";
  }
  else if(opt == 1){
    energy = "10000";
    eLbl = "10 GeV";
  }
  else{
    cout<<"Energy option is not supported"<<endl;
    return;
  }
  
  fbase += energy;
  
  // Make canvas
  TCanvas* c = makeCanvas("C");
  c->SetLogy();
  c->SetLogx();

  // Constant line
  TF1* f1 = new TF1("f1","0.1",0,100);
  f1->SetLineColor(kBlack);
  f1->SetLineStyle(2);

  // Make Legend
  TLegend* leg = makeLegend(0.6,0.8,0.7,0.9);

  // Get Dummy histogram 
  getDummy()->Draw();

  // Open the files
  // Load plots
  // Save
  for(int f=0; f<nFiles; ++f){

    // Open file
    TString fname = fbase;
    fname += "_" + materials[f] + ".root";
    TFile* file = new TFile(fname.Data());
    
    // Get Profile
    TProfile* prof = file->Get(pname.Data());
    printRadius(prof, materials[f]);

    
    // Make the graph from the profile
    TGraph* graph  = makeGraph(prof);
    setAtt(graph,colors[f],markers[f],"Radius [cm]","U/E_{0}");
    graph->Draw("sameep");
    //if(f == 0) graph->Draw("ep");
    //else       graph->Draw("epsame");
    
    // Get Radius
    float RM =  getRM(graph);
    stringstream ss;
    ss << materials[f] << ": R_{M} = " << RM;
    

    // Add to legend
    leg->AddEntry(graph, ss.str().c_str(), "ep");
    leg->SetHeader(("E_{beam} = " + eLbl).Data());

  }// end loop over files
  
  leg->Draw("same");
  f1->Draw("same");

  // Save in plots/Moliere Dir
  TString save = "plots/Moliere/beam_"+energy+".png";
  c->SaveAs(save.Data());

}

//------------------------------------------------//
// Get Dummy histogram
//------------------------------------------------//
TH1F* getDummy()
{

  TH1F* h = new TH1F("h","h",100,1,50);
  h->SetMaximum(1);
  h->SetMinimum(0.01);
  h->GetXaxis()->SetTitle("Radius [cm]");
  h->GetYaxis()->SetTitle("U/E_{0}");
  h->SetStats(0);
  h->SetTitle("");

  return h;

}

//------------------------------------------------//
// Make a graph from Profile
//------------------------------------------------//
TGraph* makeGraph(TProfile* prof)
{
  
  float x[100];
  float y[100];
  
  int nbins = prof->GetNbinsX();
  for(int bin=0; bin<nbins; ++bin){
    x[bin] = prof->GetBinLowEdge(bin+1);// + prof->GetBinWidth(bin+1);
    //x[bin] = prof->GetBinCenter(bin+1);
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
// Get moliere Radius
//------------------------------------------------//
float getRM(TGraph* gr)
{

  double x, y;
  int nP = gr->GetN();
  for(int p=1; p<nP; ++p){
    gr->GetPoint(p,x,y);
    if( y < 0.1 ){
      gr->GetPoint(p-1,x,y);
      return x;
    }
  }
  
  return -100;

}

//------------------------------------------------//
// Print radius
//------------------------------------------------//
void printRadius(TProfile* prof, TString lbl)
{

  // Threshold for Moliere Radius is 0.1
  float threshold = 0.1;
  
  // Find out when U/E_0 < 0.1
  int binsx = prof->GetNbinsX();
  for(int bin=1; bin<=binsx; ++bin){
    float bc = prof->GetBinContent(bin);
    if( bc == 0 ) continue;
    if( bc < threshold ){
      float R  = prof->GetBinLowEdge(bin);
      cout<<lbl<<" R = "<<R<<endl;
      //cout<<lbl<<" R = "<<prof->GetBinCenter(bin)<<endl;
      break;
    }// end if crosses threshold
  }// end loop over bins

}


#include "myHist.C"

TString m_savedir = "plots/EField/CherenkovPlots/";
bool m_save       = true; //false;


//----------------------------------//
// Main
//----------------------------------//
void FormatCherenkov()
{

  // Specify the files here
  vector<TString> fnames;
  //fnames.push_back("plots/TestAug19_50Evt_100GeV_AntArray_XYPlot.root");
  fnames.push_back("plots/Output_100Evt_40MeV_100Prim_xy_XYPlot.root");
  
  // Put the names to save the plot as
  vector<TString> snames;
  //snames.push_back("XYPlot_100GeV_50Evt");
  snames.push_back("XYPlot_40MeV_100Evt_100Prim");

  // Plot name
  TString pname = "XYFixedZ";

  // Make canvas
  TCanvas* c = makeCanvas("c");
  c->SetRightMargin(0.15);
  
  // Get plots draw and save
  for(unsigned int i=0; i<fnames.size(); ++i){
    TString fname = fnames.at(i);
    TFile* file = new TFile(fname.Data());
    
    TH2D* hist = (TH2D*) file->Get(pname.Data());
    hist->SetStats(0);
    hist->Draw("colz");

    if(m_save){
      c->SaveAs((m_savedir+snames.at(i)+".png").Data());
    }
    
  }
  

}

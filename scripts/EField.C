
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-//
// Script to plot the electric field from the vector potential.  //
// Currently this is work in progress for single pulse that will //
// later be extedned to handle the multipulse case.              //
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-//

#include "myHist.C"

// Need to specify some scaling factors depending
// on how you want to mix the events.
double m_scale  = 1;  // scaling factor for 

//-------------------------------------------------//
// Main
//-------------------------------------------------//
void EField()
{

  // Specify the input file name
  TString dir   = "efieldroot/";
  TString fname = dir + "Output_5Evt_40MeV_1000Prim_HardCodedAntenna_R10m.root";
  
  // Plot single
  //plotSingle(fname);


  // Plot multiple
  plotMultiple(fname);

}

//-------------------------------------------------//
// Method to plot single bunch
//-------------------------------------------------//
void plotSingle(TString infile)
{

  // Specify the input file
  TFile* file = new TFile(infile.Data());

  // Currently plot at the Cherenkov Angle
  TString pname = "A_AntNum_0_pos_8.27371_0_5.61656";
  
  // Get profile
  TH1D* A = ((TProfile*) file->Get(pname.Data()))->ProjectionX("single");

  // Make canvas
  TCanvas* c = makeCanvas("c");

  // Graph labels and stuff
  TString xtitle = "time [ns]";
  TString ytitle = "E [V/m]";
  int color = kBlue;

  // Get the graph and set Attributes
  TGraph* gr = getEField(A);
  setAtt(gr,xtitle,ytitle,color);

  // Draw
  gr->Draw();


}

//-------------------------------------------------//
// Plot multiple events
//-------------------------------------------------//
void plotMultiple(TString infile)
{

  // Make canvas
  TCanvas* c = makeCanvas("c");

  // Specify the antenna position
  TString antPos = "A_AntNum_0_pos_8.27371_0_5.61656";

  // Open input file
  TFile* file = new TFile(infile.Data());

  // Load the vector potential
  TH1D* A = makeSumVPot(file, 5, antPos, 0.350); 

  // Now get the E-field
  TGraph* E = getEField(A);
  
  // Set some attributes
  TString xtitle = "time [ns]";
  TString ytitle = "E [V/m]";
  int color = kBlue;
  setAtt(E, xtitle, ytitle, color);

  // Draw
  E->Draw();

}

//-------------------------------------------------//
// Method to retrive efield from vector potential
//-------------------------------------------------//
TGraph* getEField(TH1D* A)
{

  // Graph points
  int npoints = 0;
  float t[100000];
  float E[100000];

  // Calculat the E-field
  int nbins = A->GetNbinsX();
  for(int bin=1; bin<nbins; ++bin){
    float A0 = A->GetBinContent(bin);
    float A1 = A->GetBinContent(bin+1);
    float bw = A->GetXaxis()->GetBinWidth(bin);
    
    float time = A->GetXaxis()->GetBinCenter(bin);
    
    t[bin-1] = time;
    E[bin-1] = (A1-A0)/(bw*1e-9);
    npoints ++;
  }
  
  // Make and return graph
  TGraph* gr = new TGraph(npoints, t, E);
  return gr;
  
}

//-------------------------------------------------//
// Make vector potential for N Events
//-------------------------------------------------//
TH1D* makeSumVPot(TFile* file, int nEvents,
		      TString antPos, float dt)
{
  
  // Get the first profile
  TString pname = makePName(0, antPos);
  TH1D* A = ((TProfile*) file->Get(pname.Data()))->ProjectionX("sum");
  int nbins = A->GetNbinsX();
  
  // Now loop over the events, get other profiles, and reset
  // the bin content of the histogram including the dt shift
  TProfile* temp = NULL;
  for(int iEvt=1; iEvt<nEvents; ++iEvt){
    // Get next profile
    pname = makePName(iEvt,antPos);
    temp = (TProfile*) file->Get(pname.Data());
    
    // Now add the bin content, including the offset
    for(int bin = 1; bin<=nbins; ++bin){
      float t = temp->GetBinCenter(bin) + dt*iEvt;
      float bc = temp->GetBinContent(bin);
      int newbin = A->FindBin(t);
      
      // Check we are not out of range
      if( newbin > nbins ) continue;

      // Add content for this bin
      float oldbc = A->GetBinContent(newbin);
      A->SetBinContent(newbin, oldbc + bc);
      
    }// end loop over bins

    // Delete temporary histogram
    delete temp;

  }// end loop over events

  // Return hist
  return A;

}

//-------------------------------------------------//
// Pname builder
//-------------------------------------------------//
TString makePName(int Evt, TString antPos)
{

  stringstream ss;
  ss << antPos << "_Event" << Evt;

  return TString(ss.str().c_str());

}
    

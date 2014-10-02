
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
// Script to read the output and plot the beam profiles //
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//

#include "myHist.C"

// Saving options
TString m_savedir = "plots/beamProfile/";

//---------------------------------------//
// Main
//---------------------------------------//
void PlotBeamProfile()
{

  // Generic things
  TCanvas* c = makeCanvas("c");  
  c->SetTopMargin(0.1);
  c->SetRightMargin(0.15);

  // Files to include
  vector<TString> fnames;
  fnames.push_back("beamProfile/Gauss.txt");
  fnames.push_back("beamProfile/Flat.txt");

  // Titles to add
  vector<TString> titles;
  titles.push_back("Gaussian #sigma = 3.5");
  titles.push_back("Uniform radius = 3.5");

  // Save names
  vector<TString> savenames;
  savenames.push_back("gauss");
  savenames.push_back("flat");

  // Histogram object
  TH2F* h = NULL;

  // Loop over
  for(unsigned int i=0; i<fnames.size(); ++i){
    
    // Make hist
    h = makeHist2("hist",100,-10,10,100,-10,10,"x [mm]","y [mm]","Entries/bin");
  
    double x = 0;
    double y = 0;
    ifstream input (fnames.at(i).Data());
  
    double prev = -9999;
    while( !input.eof() ){
    
      input >> x >> y;
    
      if( x == prev ) continue;
      prev = x;
    
      h->Fill(x,y);

    }// end loop over file

    
    // Draw
    h->SetTitle(titles.at(i).Data());
    h->Draw("colz");

    // Save
    c->SaveAs((m_savedir + savenames.at(i) + ".png").Data());
    delete h;
    input.close();

  }// end loop over files

}


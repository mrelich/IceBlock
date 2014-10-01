
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
// This script is meant to help visualizing the impact of including //
// the refractive index of ice.                                     //
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//

#include "myHist.C"

bool m_save       = true;
TString m_savedir = "plots/EField/RefractedEField/";

//---------------------------------------------------//
// Main
//---------------------------------------------------//
void refractedAngles()
{

  // Specify the input file containing 
  // the antenna information
  ifstream m_file ("antennaConfig/xzRefracted_40MeV.txt");

  // Plot how the z-direction changes for different angles
  plotZChange(m_file);

}

//---------------------------------------------------//
// Plot change in Z
//---------------------------------------------------//
void plotZChange(ifstream &file)
{

  // Setup holders for z positions and angles
  int npoints = 0;
  double Z[1000];
  double Zprime[1000];
  double Angles[1000];
  
  // Loop over file and load results
  double dummy;
  double z=0, zprime=0, angle=0;
  while( !file.eof() ){

    getLine(file,dummy,dummy,z,angle,dummy,zprime);
    
    // protect against last line
    if(dummy == z) continue;
    
    // Store values
    Z[npoints]      = z;
    Zprime[npoints] = zprime;
    Angles[npoints] = angle;
    npoints++;

  }// end loop over file

  // Define graph attributes
  TString ytitle = "z [m]";
  TString xtitle = "Angle [deg]";
  int colors[]   = {kBlue, kRed};

  // Make graphs
  TGraph* gr_nom = new TGraph(npoints,Angles,Z);
  TGraph* gr_ref = new TGraph(npoints,Angles,Zprime);
  setAtt(gr_nom,xtitle,ytitle,colors[0]);
  setAtt(gr_ref,xtitle,ytitle,colors[1]);

  // Differ from generic and add markers
  gr_nom->SetMarkerStyle(20);
  //gr_nom->SetMarkerSize(1);
  gr_ref->SetMarkerStyle(20);
  //gr_ref->SetMarkerSize(1);

  // Make canvas
  TCanvas* c =makeCanvas("c");
  
  // make Legend
  TLegend* leg = makeLegend(0.7,0.8,0.8,0.93);
  leg->AddEntry(gr_nom, "Nominal","l");
  leg->AddEntry(gr_ref, "Refracted","l");

  // Make a frame
  TH1F* frame = makeFrame(gr_nom, gr_ref);
  
  // Draw
  frame->Draw();
  gr_nom->Draw("samelp");
  gr_ref->Draw("samelp");
  leg->Draw("same");

  // Save
  if( m_save )
    c->SaveAs((m_savedir+"Angle_Vs_Zpos.png").Data());

}

//---------------------------------------------------//
// Method to extract info from file
//---------------------------------------------------//
void getLine(ifstream &file, double &x, double &y, double &z,
	     double &angle, double &refAngle, double &zprime)
{

  file >> x >> y >> z >> angle >> refAngle >> zprime;

}

//-------------------------------------------------//
// Make a frame histogram
//-------------------------------------------------//
TH1F* makeFrame(TGraph* gr1, TGraph* gr2)
{

  // Get minimum and maximum for x and y
  Double_t xmin = 9999; 
  Double_t xmax = -9999;
  Double_t ymin = 9999;
  Double_t ymax = -9999;
  Double_t x = 0, y=0;
  // Loop over points
  int np = gr1->GetN();
  for(int ip=0; ip<np; ++ip){
    gr1->GetPoint(ip,x,y);
    if( x > xmax ) xmax = x;
    if( x < xmin ) xmin = x;
    if( y > ymax ) ymax = y;
    if( y < ymin ) ymin = y;

    gr2->GetPoint(ip,x,y);
    if( x > xmax ) xmax = x;
    if( x < xmin ) xmin = x;
    if( y > ymax ) ymax = y;
    if( y < ymin ) ymin = y;
  }

  // Make histogram
  TString xtitle = gr1->GetXaxis()->GetTitle();
  TString ytitle = gr1->GetYaxis()->GetTitle();
  TH1F* h = makeHist("frame",1,xmin,xmax,xtitle,ytitle,kWhite,0);
  h->SetMinimum(1.1*ymin);
  h->SetMaximum(1.1*ymax);
  
  return h;

}

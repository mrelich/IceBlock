
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-//
// Script to plot the electric field from the vector potential.  //
// Currently this is work in progress for single pulse that will //
// later be extedned to handle the multipulse case.              //
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-//

// TODO: Clean up plotMaxEVsZ and make sure plotRefractedCompare still works

#include "myHist.C"

// Need to specify some scaling factors depending
// on how you want to mix the events.
double m_scale  = 1;  // scaling factor for 

// Save information
bool m_save = false;
//TString m_savedir = "plots/EField/RefractedEField/";
TString m_savedir = "plots/EField/RefractedEField_BinningCheck/";

//-------------------------------------------------//
// Main
//-------------------------------------------------//
void EField(int opt = 0, bool save=false)
{
  
  m_save = save;

  // Specify the input file name
  TString dir   = "efieldroot/";
  //TString fname = dir + "Output_5Evt_40MeV_1000Prim_HardCodedAntenna_R10m.root";
  //TString fname = dir + "Output_5Evt_40MeV_10000Prim_HardCodedAntenna_R10m.root";
  //TString fname = dir + "Output_57Evt_40MeV_1000Prim_HardCodedAntenna_ROldConfig.root";
  //TString fname = dir + "Output_5Evt_40MeV_10000Prim_xzUnRefracted_40MeV.root";
  TString fname = dir + "Output_5Evt_40MeV_10000Prim_xzRefracted_40MeV.root";

  // Calculate the scale factor
  m_scale = 1.e9 / 5 / 10000; 
  cout<<"Scaling by : "<<m_scale<<endl;
 
  // Plot single
  if(opt == 0) 
    plotSingle(fname);

  // Plot multiple
  if(opt == 1)
    plotMultiple(fname);

  // Plot refracted result
  if(opt == 2){
    cout<<"remember to set the files in this method!!!"<<endl;
    // REDO THIS METHOD. I should just read in from one file. 
    // include the scale factor for Rprime/R and I should add the time
    // delay which corresponds to n*(R-Rprime)/c (in seconds)
    plotRefractedCompare(true);
  }

  // Plot |E| vs z position
  if(opt == 3)
    plotMaxEVsZ();
  
  // Investigate vector potential 
  if(opt == 4)
    plotVPComp();

}

//-------------------------------------------------//
// Method to plot single bunch
//-------------------------------------------------//
void plotSingle(TString infile)
{

  // Specify the input file
  TFile* file = new TFile(infile.Data());

  // Currently plot at the Cherenkov Angle
  //TString pname = "A_AntNum_0_pos_8.27371_0_5.61656";
  //TString pname = "A_AntNum_0_pos_6_0_4.04";
  //TString pname = "A_AntNum_26_pos_6_0_4.04705";
  TString pname = "A_AntNum_26_pos_6_0_0.919403";
  
  // Get profile
  TH1D* A = ((TProfile*) file->Get(pname.Data()))->ProjectionX("single");
  A->Scale(m_scale);

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
  //TString antPos = "A_AntNum_0_pos_8.27371_0_5.61656";
  TString antPos = "A_AntNum_0_pos_6_0_4.04";

  // Open input file
  TFile* file = new TFile(infile.Data());

  // Load the vector potential
  TH1D* A = makeSumVPot(file, 57, antPos, 0.350); 
  A->Scale(m_scale);

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
TGraph* getEField(TH1D* A, double dt = 0)
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
    
    t[bin-1] = time + dt;
    E[bin-1] = (A1-A0)/(bw*1e-9);
    npoints ++;
  }
  
  // Make and return graph
  TGraph* gr = new TGraph(npoints, t, E);
  return gr;
  
}

//-------------------------------------------------//
// Plot E-field comparing nominal vs. refracted ang
//-------------------------------------------------//
void plotRefractedCompare(bool doSingle)
{

  // Specify the antenna files to consider
  // Important: It is assumed file lengths are the 
  // same, so there is a one-to-one coreespondence
  ifstream m_nom ("antennaConfig/xzRefracted_40MeV.txt");

  // Specify the root files that corespond to the 
  // two input files above
  TString indir = "efieldroot/";
  //TFile* f_nom = new TFile((indir+"Output_5Evt_40MeV_10000Prim_xzRefracted_40MeV.root").Data());
  TFile* f_nom = new TFile((indir+"Output_50Evt_40MeV_10000Prim_xzRefracted_40MeV_BinningCheck.root").Data());

  // Make canvas
  TCanvas* c = makeCanvas("c");
  
  // Make Legend
  TLegend* leg = makeLegend(0.7,0.8,0.8,0.93);

  // Graph labels and stuff
  TString xtitle = "time [ns]";
  TString ytitle = "E [V/m]";
  int colors[] = {kBlue,kRed};

  // Loop over files and get the antenna positions
  int counter = 0;
  TString antPos   = "";
  double angle     = 0;
  double refAngle  = 0;
  double  R        = 0;
  double  Rprime   = 0;
  double  x=0, y=0, z=0;
  double zprime = 0;
  double prevZ = 0;

  while( !m_nom.eof() ){
  //while( !m_nom.eof() && counter < 1){
    
    getAntPos(m_nom, counter, antPos, 
	      angle, refAngle, 
	      x,y,z,
	      R, Rprime,
	      zprime);
    counter++;

    // Protect against last step going 
    // beyond the end of file
    if( z == prevZ ) continue;
    prevZ = z;


    cout<<antPos<<" "<<angle<<" "<<refAngle<<" "
	<<R<<" "<<Rprime<<endl;
    
    // So now I have the antenna names from the file
    // we can plot single or multiple plots. Go with 
    // single first
    TH1D* A_nom = NULL;
    TH1D* A_ref = NULL;
    if(doSingle)
      A_nom = ((TProfile*) f_nom->Get(antPos.Data()))->ProjectionX("singleNom");
    else
      A_nom = makeSumVPot(f_nom, 5, antPos, 0.350); 

    // Copy nominal to reference
    A_ref = (TH1D*) A_nom->Clone("refracted");
    
    // Scale
    A_nom->Scale(m_scale);
    A_ref->Scale(m_scale * R/Rprime);

    // Get nominal graph and set Attributes
    TGraph* gr_nom = getEField(A_nom);
    setAtt(gr_nom,xtitle,ytitle,colors[0]);

    // Get the refracted result with timing offset
    double dt = 1.78 * (Rprime - R)/ 2.99792458e8 / 1e-9;
    TGraph* gr_ref = getEField(A_ref, dt);
    setAtt(gr_ref,xtitle,ytitle,colors[1]);

    // Make a frame histogram
    TH1F* frame = makeFrame(gr_nom, gr_ref);

    // Draw
    frame->Draw();
    gr_nom->Draw("same");
    gr_ref->Draw("same");

    // Add legend info
    leg->Clear();
    leg->SetHeader(Form("Angle: %.0f#circ",angle));
    leg->AddEntry(gr_nom, "Nominal","l");
    leg->AddEntry(gr_ref, "Refracted", "l");
    leg->Draw("same");
    
    // Need to delete the histograms and graphs to avoid 
    // the memory leaks.  Only clean up for case where we
    // save the plots to png.
    if( m_save ){
      TString savename = m_savedir + antPos;
      if(doSingle) savename += "_singlePulse";
      else         savename += "_multiPulse";
      savename += ".png";
      c->SaveAs(savename.Data());

      delete frame;
      delete gr_nom;
      delete gr_ref;
      delete A_nom;
      delete A_ref;
    }// end if save

  }// end loop over antennas

}

//-------------------------------------------------//
// Plot E-field comparing nominal vs. refracted ang
//-------------------------------------------------//
void plotMaxEVsZ()
{

  // Specify the antenna files to consider
  // Important: It is assumed file lengths are the 
  // same, so there is a one-to-one coreespondence
  ifstream m_nom ("antennaConfig/xzRefracted_40MeV.txt");

  // Specify the root files that corespond to the 
  // two input files above
  TString indir = "efieldroot/";
  //TFile* f_nom = new TFile((indir+"Output_5Evt_40MeV_10000Prim_xzRefracted_40MeV.root").Data());
  //TFile* f_nom = new TFile((indir+"Output_50Evt_40MeV_10000Prim_xzRefracted_40MeV.root").Data());
  TFile* f_nom = new TFile((indir+"Output_50Evt_40MeV_10000Prim_xzRefracted_40MeV_BinningCheck.root").Data());

  // Make canvas
  TCanvas* c = makeCanvas("c");
  
  // Make Legend
  TLegend* leg = makeLegend(0.17,0.33,0.75,0.93);

  // Graph labels and stuff
  TString legTitle = "Ice Tilt = 30#circ";
  TString xtitle = "Angle Relative to Beam [deg]";
  TString ytitle = "Maximum |E| [Vm]";
  int colors[] = {kBlue,kRed};

  // Holders for graph
  int counter = 0;
  double Z[1000];
  double Zprime[1000];
  double maxE[1000];
    
  // Holders
  TString antPos   = "";
  double angle    = 0;
  double refAngle = 0;
  double x=0,y=0,z=0;
  double R=0, Rprime = 0;
  double zprime = 0;
  double prevZ  = 0;
  int npoints = 0;

  // Loop over files and get the antenna positions
  while( !m_nom.eof() ){
    
    getAntPos(m_nom, counter, antPos, angle, refAngle, 
	      x,y,z,
	      R, Rprime,
	      zprime);
    counter++;
    cout<<"Working on : "<<counter-1<<" "<<antPos<<endl;

    // Protect against last step going 
    // beyond the end of file
    if( z == prevZ ) continue;
    prevZ = z;

    // Get Vector potential
    TH1D* A = ((TProfile*) f_nom->Get(antPos.Data()))->ProjectionX("singleNom");

    // Scale
    A->Scale(m_scale * R/Rprime);
    
    // Get the refracted result with timing offset
    double dt = 1.78 * (Rprime - R)/ 2.99792458e8 / 1e-9;
    TGraph* gr_E = getEField(A, dt);

    // Now get the maximum E field
    int np = gr_E->GetN();
    double t =0, E=0;
    double absE = -999;
    for(int ip=0; ip<np; ++ip){
      gr_E->GetPoint(ip,t,E);
      if( fabs(E) > absE ) absE = fabs(E);
    }

    // Don't need to count the z points
    

    // Now store results
    double conv = 180/TMath::Pi();

    if( atan(x/zprime) < 0 ) continue;
    Z[npoints]      = atan(x/z)*conv; ////z;
    Zprime[npoints] = atan(x/zprime)*conv; //zprime;
    maxE[npoints]   = absE;
    npoints++;

    // Delete objects
    delete A;
    delete gr_E;
    

  }// end loop over file
  
  // Make graphs
  TGraph* gr_nom = new TGraph(npoints, Z, maxE);
  setAtt(gr_nom,xtitle,ytitle,colors[0]);
  TGraph* gr_ref = new TGraph(npoints, Zprime, maxE);
  setAtt(gr_ref,xtitle,ytitle,colors[1]);

  // Make a frame histogram
  TH1F* frame = makeFrame(gr_nom, gr_ref);
  
  // Draw
  frame->Draw();
  gr_nom->Draw("same");
  gr_ref->Draw("same");

  // Add legend info
  leg->Clear();
  leg->SetHeader("Ice Tilt = 30#circ");
  leg->AddEntry(gr_nom, "Nominal","l");
  leg->AddEntry(gr_ref, "Refracted", "l");
  leg->Draw("same");
    
  // Need to delete the histograms and graphs to avoid 
  // the memory leaks.  Only clean up for case where we
  // save the plots to png.
  if( m_save ){
    TString savename = m_savedir + "EFieldVsZ.png";
    c->SaveAs(savename.Data());

  }// end if save
  
  
}

//-------------------------------------------------//
// Get antenna position for a given fil
//-------------------------------------------------//
void getAntPos(ifstream &infile, 
	       int &counter,
	       TString &antPos,
	       double &angle, double &refAngle,
	       double &x, double &y, double &z,
	       double &R, double &Rprime, double &zprime)
{
  
  // Specify variables to read in the lines in the pos
  infile >> x >> y >> z >> angle >> refAngle >> zprime;
  //cout<<"\t"<<x<<" "<<y<<" "<<z<<endl;

  // make the strings
  stringstream ss;
  ss << "A_AntNum_" << counter << "_pos_"
     << x << "_" << y << "_" << z;
  antPos = TString(ss.str().c_str());
  
  // Save R and Rprime
  R      = sqrt(x*x+y*y+z*z);
  Rprime = sqrt(x*x+y*y+zprime*zprime);

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

//-------------------------------------------------//
// Vector potential comparison for studies
void plotVPComp()
{

  // Specify the input file
  //TFile* file = new TFile("efieldroot/Output_50Evt_40MeV_10000Prim_xzRefracted_40MeV.root");
  TFile* file = new TFile("efieldroot/Output_50Evt_40MeV_10000Prim_xzRefracted_40MeV_BinningCheck.root");

  // Currently plot at the Cherenkov Angle
  vector<TString> pnames;
  pnames.push_back("A_AntNum_27_pos_6_0_3.89645");
  pnames.push_back("A_AntNum_28_pos_6_0_3.74922");
  pnames.push_back("A_AntNum_29_pos_6_0_3.60516");
  
  // Names for legend
  vector<TString> names;
  names.push_back("57#circ");
  names.push_back("58#circ");
  names.push_back("59#circ");

  // Make canvas
  TCanvas* c = makeCanvas("c");
  c->SetLogy();

  // Make legend
  TLegend* leg = makeLegend(0.15,0.3,0.75,0.93);

  // Titles and shiz for hists
  TString xtitle = "t [ns]";
  TString ytitle = "|A| [Vs/m]";
  int color[] = {kBlack, kRed, kBlue,kMagenta};

  // Loop and get the profile
  TProfile* profs[3];
  float maximum = -99999;
  for(unsigned int i=0; i<pnames.size(); ++i){
    TString pname = pnames.at(i);
    TString name  = names.at(i);

    profs[i] = getProfile(file,pname,xtitle,ytitle,color[i],20);
    leg->AddEntry(profs[i], name.Data(), "le");

    if(maximum < profs[i]->GetMaximum())
      maximum = profs[i]->GetMaximum();

  }

  // Set max
  profs[0]->SetMaximum(5*maximum);
  profs[0]->SetMinimum(1e-3*maximum);

  // Set min and max for x range
  float xmin = 41; float xmax = 44;
  float xminB = profs[0]->GetXaxis()->FindBin(xmin);
  float xmaxB = profs[0]->GetXaxis()->FindBin(xmax);
  profs[0]->GetXaxis()->SetRange(xminB,xmaxB);

  // Draw
  profs[0]->Draw();
  for(unsigned int i=1; i<pnames.size(); ++i)
    profs[i]->Draw("same");
  
  leg->Draw("same");

  // Save
  if( m_save )
    c->SaveAs((m_savedir+"VPCheckForSpikes_LargerBins.png").Data());

}

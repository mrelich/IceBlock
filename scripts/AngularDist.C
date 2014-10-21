
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-//
// Plotting script needed to make the angular distribution for //
// several frequencies to compare with Keiichi.  The goal is   //
// to determine where the peak actually is and what angle we   //
// should be using for the iceblock.                           //
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-//

#include "AraSees.C"

// Save information
TString m_savedir = "plots/EField/Angular/";

//-----------------------------------------------//
// Main
//-----------------------------------------------//
void AngularDist()
{
  
  // Specify the input file
  TString rootdir = "rootfiles/EField/";
  //TFile* f_in = new TFile(rootdir+"efield_50Evt_40MeV_10000Prim_xzRef.root");
  //TFile* f_in = new TFile(rootdir+"efield_50Evt_40MeV_10000Prim_angleScan_R100m.root");
  //TFile* f_in = new TFile(rootdir+"efield_50Evt_40MeV_10000Prim_angleScan_R100m_singlePulse.root");
  //TFile* f_in = new TFile(rootdir+"efield_50Evt_40MeV_100Prim_angleScan_R100m_singlePulse.root");
  //TFile* f_in = new TFile(rootdir+"efield_50Evt_40MeV_100Prim_angleScan_R100m_RandFlat3.5.root");
  //TFile* f_in = new TFile(rootdir+"Testing_ShiftedSinglePos.root");
  //TFile* f_in = new TFile(rootdir+"efield_50Evt_40MeV_10000Prim_angleScan_100M_singlePos.root");
  //TFile* f_in = new TFile(rootdir+"efield_50Evt_40MeV_10000Prim_angleScan_100M_singlePos_thresh1MeV.root");
  //TFile* f_in = new TFile(rootdir+"efield_50Evt_40MeV_10000Prim_angleScan_100M_RandFlat3.5.root");
  TFile* f_in = new TFile(rootdir+"efield_50Evt_40MeV_10000Prim_angleScan_100M_singlePos_TestingEndpointIssue.root");  

  // Specify the antenna file
  //ifstream f_ant ("antennaConfig/xzRefracted_40MeV.txt");
  //ifstream f_ant ("antennaConfig/AngleScane_100m.txt");
  ifstream f_ant ("antennaConfig/AngleScan_100m.txt");

  // Plot
  plot(f_in,f_ant);
  //plotFT(f_in,f_ant);
  //compareAnalytic(f_in,f_ant);
  //compareMultipleFiles();
}

//-----------------------------------------------//
// Plot peak for several frequencies
//-----------------------------------------------//
void plot(TFile* f_in, ifstream &f_ant)
{

  // For now, just do this for one frequency
  //float freq = 800.e6;
  vector<float> freqs;
  freqs.push_back(800.e6);
  freqs.push_back(600.e6);
  freqs.push_back(400.e6);
  freqs.push_back(200.e6);

  vector<TString> fnames;
  fnames.push_back("800 MHz");
  fnames.push_back("600 MHz");
  fnames.push_back("400 MHz");
  fnames.push_back("200 MHz");

  int colors[] = {kBlack,kViolet,kRed,kBlue};

  // Specify some stuff for the graph
  TString xtitle = "Angle Relative to Beam [deg]";
  TString ytitle = "|RE(f)| [V]";
  //TString ytitle = "Normalized";
  //TString ytitle = "|E(f)| [V]";

  // Place holders for result
  int npoints = 0;
  float angles[10000];
  vector<vector<float> > v_E_freqs;
  for(unsigned int i=0; i<freqs.size(); ++i){
    vector<float> temp;
    v_E_freqs.push_back( temp );
  }
  
  // Loop over the antenna positions to 
  // get the input histograms.
  //TString x, y, z;
  float x=0, y=0, z=0;
  float angle=0;
  float dummy=0;
  //TString prevZ = "";
  float prevZ = 0;
  int counter = 0;
  //stringstream ss;
  float maximum = -999;
  while( !f_ant.eof() ){
    
    // Get antenna positions
    f_ant >> x >> y >> z >> angle >> dummy >> dummy;
    if( z == prevZ ) continue;
    prevZ = z;

    float R = sqrt(x*x+y*y+z*z);

    // Make antenna name
    stringstream ss; //ss.str("");
    ss << "A_AntNum_" << counter << "_pos_" 
       << x << "_"
       << y << "_"
       << z;
    
    string antPos = ss.str();
    counter++;


    // Get the histogram from input file
    TH1F* h_Et = (TH1F*) f_in->Get(antPos.c_str());
    cout<<antPos<<" "<<h_Et<<endl;
    if( !h_Et ) continue;
    // Get fourier transform
    // First obtain the sum of the electric field
    int nbins = h_Et->GetNbinsX();
    float xmin = h_Et->GetBinLowEdge(1);
    float xmax = h_Et->GetBinLowEdge(nbins) + h_Et->GetBinWidth(nbins);
    //cout<<"Min: "<<xmin<<" max: "<<xmax<<" nbins: "<<nbins<<" Step: "<<((xmax-xmin)/nbins)*1e9<<endl;

    // Now fourier transform that shit
    TH1F* h_Ef = h_Et->FFT(NULL,"MAG");
    TH1F* h_Ef_format = remake(h_Ef,nbins,xmin,xmax);

    // Now get the point we are interested in
    angles[npoints] = angle;
    //int bin = h_Ef_format->FindBin(freq);
    //E_freq[npoints] = h_Ef_format->GetBinContent(bin);
    //npoints++;

    for(unsigned int i=0; i<freqs.size(); ++i){
      int bin = h_Ef_format->FindBin(freqs.at(i));
      v_E_freqs.at(i).push_back( R*h_Ef_format->GetBinContent(bin) );
      //v_E_freqs.at(i).push_back( h_Ef_format->GetBinContent(bin) );
      if( maximum < v_E_freqs.at(i).back() )
	maximum = v_E_freqs.at(i).back();
    }
    npoints++;

    delete h_Et;
    delete h_Ef;
    delete h_Ef_format;

  }// end loop over antenna file
  
  // Canvas
  TCanvas* c = makeCanvas("c");

  // Legend
  TLegend* leg = makeLegend(0.15,0.3,0.6,0.93);

  // Make a frame object
  TH1F* h = makeHist("frame",1,0,90,xtitle,ytitle,0,0);
  h->SetMinimum(0);
  h->SetMaximum(1.2*maximum);
  
  h->Draw();

  float scale = 0;
  for(unsigned int i=0; i<freqs.size(); ++i){
    float E_freqs[10000];
    for(int p=0; p<npoints; ++p){
      E_freqs[p] = v_E_freqs[i][p];
      scale += E_freqs[p];
    }
    TGraph* gr = new TGraph(npoints,angles,E_freqs);
    setAtt(gr,xtitle,ytitle,colors[i]);
    //scale = gr->Integral(0,-1);
    gr->Draw("same");

    leg->AddEntry(gr,fnames.at(i),"l");

  }
  leg->Draw("same");

  //c->SaveAs(m_savedir+"AngularScan_fixedX.png");
  //c->SaveAs(m_savedir+"AngularScan_singleBunch.png");

}


//-----------------------------------------------//
// Compare to analytic formula
//-----------------------------------------------//
void compareAnalytic(TFile* f_in, ifstream &f_ant)
{

  // For now, just do this for one frequency
  float freq = 800.e6;

  TString fname = "800 MHz";

  int colors[] = {kBlack,kViolet,kRed,kBlue};

  // Specify some stuff for the graph
  TString xtitle = "Angle Relative to Beam [deg]";
  TString ytitle = "Normalized";

  // Make histogram to contain the results
  TH1F* h_angle = makeHist("angle",90,0,90,xtitle,ytitle,kBlack,20);
  
  // Loop over the antenna positions to 
  // get the input histograms.
  //TString x, y, z;
  float x=0, y=0, z=0;
  float angle=0;
  float dummy=0;
  //TString prevZ = "";
  float prevZ = 0;
  int counter = 0;
  //stringstream ss;
  float maximum = -999;
  while( !f_ant.eof() ){
    
    // Get antenna positions
    f_ant >> x >> y >> z >> angle >> dummy >> dummy;
    if( z == prevZ ) continue;
    prevZ = z;

    float R = sqrt(x*x+y*y+z*z);

    // Make antenna name
    stringstream ss; //ss.str("");
    ss << "A_AntNum_" << counter << "_pos_" 
       << x << "_"
       << y << "_"
       << z;
    
    string antPos = ss.str();
    counter++;


    // Get the histogram from input file
    TH1F* h_Et = (TH1F*) f_in->Get(antPos.c_str());
    cout<<antPos<<" "<<h_Et<<endl;
    if( !h_Et ) continue;
    // Get fourier transform
    // First obtain the sum of the electric field
    int nbins = h_Et->GetNbinsX();
    float xmin = h_Et->GetBinLowEdge(1);
    float xmax = h_Et->GetBinLowEdge(nbins) + h_Et->GetBinWidth(nbins);
    //cout<<"Min: "<<xmin<<" max: "<<xmax<<" nbins: "<<nbins<<" Step: "<<((xmax-xmin)/nbins)*1e9<<endl;

    // Now fourier transform that shit
    TH1F* h_Ef = h_Et->FFT(NULL,"MAG");
    TH1F* h_Ef_format = remake(h_Ef,nbins,xmin,xmax);

    // Now get the point we are interested in
    int bin  = h_Ef_format->FindBin(freq);
    float Ef = R*h_Ef_format->GetBinContent(bin);
    bin = h_angle->FindBin(angle);
    h_angle->SetBinContent(bin, Ef);
    h_angle->SetBinError(bin,0);

    delete h_Et;
    delete h_Ef;
    delete h_Ef_format;
    
  }// end loop over antenna file
  
  // Canvas
  TCanvas* c = makeCanvas("c");

  // Legend
  TLegend* leg = makeLegend(0.15,0.3,0.6,0.93);

  // Draw result
  h_angle->Draw("p");

  //TF1* approx = getApprox(800.e6,0,90,0.068);
  //TF1* approx = getApprox(800.e6,0,90,0.068);
  TF1* approx = getApprox(800.e6,0,90,0.06);
  TH1F* htemp = (TH1F*) approx->GetHistogram();
  htemp->Scale(h_angle->Integral()/htemp->Integral());
  //htemp->DrawNormalized("psame");
  setHistAtt(htemp,"","",kRed,25);
  htemp->Draw("psame");



  //c->SaveAs(m_savedir+"AngularScan_fixedX.png");
  //c->SaveAs(m_savedir+"AngularScan_singleBunch.png");

}

//-----------------------------------------------//
// Plot peak for several frequencies
//-----------------------------------------------//
void plotFT(TFile* f_in, ifstream &f_ant)
{

  // Plot for folowing angles
  vector<int> angles;
  /*angles.push_back(0);
  angles.push_back(10);
  angles.push_back(20);
  angles.push_back(30);
  angles.push_back(40);
  angles.push_back(50);
  angles.push_back(60);
  angles.push_back(70);
  angles.push_back(80);
  angles.push_back(90);
  */
  angles.push_back(50);
  angles.push_back(52);
  angles.push_back(54);
  angles.push_back(56);
  angles.push_back(58);
  angles.push_back(60);

  int colors[] = {kBlack,kViolet,kRed,kBlue,kGreen,
                  kYellow-6, kGreen+2, kCyan+1,kBlue+2,kPink};

  // Specify some stuff for the graph
  TString xtitle = "f [Hz]";
  TString ytitle = "|RE(f)| [V]";

  // Holder for histograms
  TH1F* hists[10];
  int curHist = 0;

  // Legend
  TLegend* leg = makeLegend(0.15,0.3,0.6,0.93);
 
  // Loop over the antenna positions to 
  // get the input histograms.
  //TString x, y, z;
  float x=0, y=0, z=0;
  float angle=0;
  float dummy=0;
  //TString prevZ = "";
  float prevZ = 0;
  int counter = 0;
  //stringstream ss;
  float maximum = -999;
  while( !f_ant.eof() ){
    
    // Get antenna positions
    f_ant >> x >> y >> z >> angle >> dummy >> dummy;
    if( z == prevZ ) continue;
    prevZ = z;


    float R = sqrt(x*x+y*y+z*z);

    // Make antenna name
    stringstream ss; //ss.str("");
    ss << "A_AntNum_" << counter << "_pos_" 
       << x << "_"
       << y << "_"
       << z;
    
    string antPos = ss.str();
    counter++;

    // Check if this one is in our list
    bool found = false;
    for(unsigned int ia=0; ia<angles.size(); ++ia)
      if( angle == angles.at(ia) ) 
	found = true;
    if(!found) continue;

    // Get the histogram from input file
    TH1F* h_Et = (TH1F*) f_in->Get(antPos.c_str());
    cout<<antPos<<" "<<h_Et<<endl;
    if( !h_Et ) continue;
    // Get fourier transform
    // First obtain the sum of the electric field
    int nbins = h_Et->GetNbinsX();
    float xmin = h_Et->GetBinLowEdge(1);
    float xmax = h_Et->GetBinLowEdge(nbins) + h_Et->GetBinWidth(nbins);
    //cout<<"Min: "<<xmin<<" max: "<<xmax<<" nbins: "<<nbins<<" Step: "<<((xmax-xmin)/nbins)*1e9<<endl;

    // Now fourier transform that shit
    TH1F* h_Ef = h_Et->FFT(NULL,"MAG");
    TH1F* h_Ef_format = remake(h_Ef,nbins,xmin,xmax);
    setHistAtt(h_Ef_format,xtitle,ytitle,colors[curHist],20);
    int low = h_Ef_format->FindBin(1e6);
    int high = h_Ef_format->FindBin(4000e6);
    h_Ef_format->GetXaxis()->SetRange(low,high);
    h_Ef_format->SetName((antPos+"_FT").c_str());

    // Save
    hists[curHist] = h_Ef_format;
    curHist++;
    
    // Add to legend
    leg->AddEntry(h_Ef_format, Form("%.0f#circ",angle),"l");

    // Get maximum
    if( maximum < h_Ef_format->GetMaximum() )
      maximum = h_Ef_format->GetMaximum();

    delete h_Et;
    delete h_Ef;
    //delete h_Ef_format;

  }// end loop over antenna file
  
  // Canvas
  TCanvas* c = makeCanvas("c");
  c->SetLogy();
  c->SetLogx();
  c->SetGridx();
  c->SetGridy();

  
  // Draw
  hists[0]->SetMaximum(5*maximum);
  hists[0]->Draw();
  for(unsigned int i=1; i<angles.size(); ++i)
    hists[i]->Draw("same");
  leg->Draw("same");

  //c->SaveAs(m_savedir+"FT_MultipleAngles_5Bunch.png");

}

//---------------------------------------------//
// The Gaussian approximation
//---------------------------------------------//
TF1* getApprox(float freq, float xmin, float xmax, float length)
{

  TF1* f_approx = new TF1("Approx",GausApprox,xmin,xmax,3);
  f_approx->SetParameter(0,length);
  f_approx->SetParameter(1,freq);
  f_approx->SetParameter(2,1e9);
  
  return f_approx;
}

//---------------------------------------------//
Double_t GausApprox(Double_t* x, Double_t* par)
{

  // Constants
  double pi  = TMath::Pi();
  double mu0 = 4 * pi * 1e-7;
  double n   = 1.78;
  double c   = 3e8;

  // Variable is angle
  Double_t theta = x[0] * pi / 180.;

  // Parameters will be:
  // 0 -- L
  // 1 -- frequency
  // 2 -- Max charge
  double L = par[0];
  double w = 2 * pi * par[1];
  double Q = par[2];

  double k = w * n / c;

  double C = mu0*w*Q*L*sin(theta) / sqrt(2*pi);
  double F = TMath::Exp(-(k*k*L*L)*pow(cos(theta) - 1/n,2)/2.);

  // Return function value
  return C * F;

}

//-----------------------------------------------//
// Compare various results for different N(Primaries)
//-----------------------------------------------//
void compareMultipleFiles()
{

  // Antenna file
  ifstream f_ant ("antennaConfig/AngleScan_100m.txt");

  // Specify the files to compare
  vector<TString> f_ins;
  TString indir = "rootfiles/EField/";
  if(true){
    f_ins.push_back(indir+"efield_50Evt_40MeV_10000Prim_angleScan_100M_singlePos.root");
    f_ins.push_back(indir+"efield_50Evt_40MeV_1000Prim_AngleScan_100m_singlePos.root");
    f_ins.push_back(indir+"efield_50Evt_40MeV_100Prim_AngleScan_100m_singlePos.root");
  }
  else{
    f_ins.push_back(indir+"efield_50Evt_40MeV_10000Prim_angleScan_100M_RandFlat3.5.root");
    f_ins.push_back(indir+"efield_50Evt_40MeV_1000Prim_AngleScan_100m_RandFlat3.5.root");
    f_ins.push_back(indir+"efield_50Evt_40MeV_100Prim_AngleScan_100m_RandFlat3.5.root");
  }

  // For now, just do this for one frequency
  float freq = 800.e6;

  vector<TString> fnames;
  fnames.push_back("10000");
  fnames.push_back("1000");
  fnames.push_back("100");

  int colors[] = {kBlack,kViolet,kRed,kBlue};

  // Specify some stuff for the graph
  TString xtitle = "Angle Relative to Beam [deg]";
  TString ytitle = "|RE(f)| [V]";

  // Place holders for result
  int npoints = 0;
  float angles[10000];
  vector<vector<float> > v_E_freqs;
  for(unsigned int i=0; i<f_ins.size(); ++i){
    vector<float> temp;
    v_E_freqs.push_back( temp );
  }
  
  // Loop over the antenna positions to 
  // get the input histograms.
  //TString x, y, z;
  float x=0, y=0, z=0;
  float angle=0;
  float dummy=0;
  //TString prevZ = "";
  float prevZ = 0;
  int counter = 0;
  //stringstream ss;
  float maximum = -999;
  TFile* f_in = NULL;
  while( !f_ant.eof() ){
    
    // Get antenna positions
    f_ant >> x >> y >> z >> angle >> dummy >> dummy;
    if( z == prevZ ) continue;
    prevZ = z;
    
    float R = sqrt(x*x+y*y+z*z);
    
    // Make antenna name
    stringstream ss; //ss.str("");
    ss << "A_AntNum_" << counter << "_pos_" 
       << x << "_"
       << y << "_"
       << z;
    
    string antPos = ss.str();
    counter++;

    // Now get the point we are interested in
    angles[npoints] = angle;
    
    // Loop over the input files
    for(unsigned int f=0; f<f_ins.size(); ++f){
      f_in = new TFile(f_ins.at(f));
    
      // Get the histogram from input file
      TH1F* h_Et = (TH1F*) f_in->Get(antPos.c_str());
      cout<<antPos<<" "<<h_Et<<endl;
      if( !h_Et ) continue;
      // Get fourier transform
      // First obtain the sum of the electric field
      int nbins = h_Et->GetNbinsX();
      float xmin = h_Et->GetBinLowEdge(1);
      float xmax = h_Et->GetBinLowEdge(nbins) + h_Et->GetBinWidth(nbins);
      //cout<<"Min: "<<xmin<<" max: "<<xmax<<" nbins: "<<nbins<<" Step: "<<((xmax-xmin)/nbins)*1e9<<endl;

      // Now fourier transform that shit
      TH1F* h_Ef = h_Et->FFT(NULL,"MAG");
      TH1F* h_Ef_format = remake(h_Ef,nbins,xmin,xmax);


      int bin = h_Ef_format->FindBin(freq);
      v_E_freqs.at(f).push_back( R*h_Ef_format->GetBinContent(bin) );
      //v_E_freqs.at(f).push_back( h_Ef_format->GetBinContent(bin) );
      if( maximum < v_E_freqs.at(f).back() )
	maximum = v_E_freqs.at(f).back();

      //delete h_Et;
      delete h_Ef;
      delete h_Ef_format;
      f_in->Close();
      delete f_in;
    }
    npoints++;

    
  }// end loop over antenna file
  
  // Canvas
  TCanvas* c = makeCanvas("c");

  // Legend
  TLegend* leg = makeLegend(0.15,0.3,0.6,0.93);

  // Make a frame object
  TH1F* h = makeHist("frame",1,0,90,xtitle,ytitle,0,0);
  h->SetMinimum(0);
  h->SetMaximum(1.2*maximum);
  
  h->Draw();

  float scale = 0;
  for(unsigned int i=0; i<f_ins.size(); ++i){
    float E_freqs[10000];
    for(int p=0; p<npoints; ++p){
      E_freqs[p] = v_E_freqs[i][p];
      scale += E_freqs[p];
    }
    TGraph* gr = new TGraph(npoints,angles,E_freqs);
    setAtt(gr,xtitle,ytitle,colors[i]);
    //scale = gr->Integral(0,-1);
    gr->Draw("same");

    leg->AddEntry(gr,fnames.at(i),"l");

  }
  leg->Draw("same");

  //c->SaveAs(m_savedir+"AngularScan_fixedX.png");
  //c->SaveAs(m_savedir+"AngularScan_singleBunch.png");

}

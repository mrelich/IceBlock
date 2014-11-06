
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-//
// The purpose of this script is to estimate what ARA detector //
// will see at the ELS test site. We will extract the voltage  //
// result by:                                                  //
//    1.) Transform E(t) --> E(f)                              //
//    2.) Get V(f) = E(f) * g(f), where g(f) is the freq       //
//        dependent gain in units of m.                        //
//    3.) Transform V(f) --> V(t) and plot                     //
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-//

#include "myHist.C"
#include "BuildSum.C"
#include "TVirtualFFT.h"

TString savedir = "plots/EField/FourierPlots/";

//--------------------------------------------//
// Main
//--------------------------------------------//
void AraSees(TString efname,    // File containing electric field
	     TString gfname,    // File containing the gain
	     TString ename,     // Specify histogram name for e-field
	     TString gname,     // The histogram name for the gain
	     int opt)           // Input options
{

  // Open the files
  TFile* efile = new TFile(efname.Data());
  TFile* gfile = new TFile(gfname.Data());

  // Plot sin
  if( opt == 0 ){
    TString freq = "5";
    plotSin(freq);
  }
  else if( opt <= 3 ){
    plotSumE(efile, gfile, ename, gname, opt);
  }
  else return;

}

//--------------------------------------------//
// Plot Sin
//--------------------------------------------//
void plotSin(TString freq)
{

  // Make a canvas
  TCanvas* c = makeCanvas("c");

  // for saving
  stringstream save;

  // Build a sin function with the input frequency
  TF1* f = new TF1("f",("sin(2*TMath::Pi()*"+freq+"*x)").Data(),0,10);

  // Now fill histogram with that information
  int nbins = 10000;
  float xmin  = 0;
  float xmax  = 5;
  TH1F* h_sinx = makeHist("sinx",nbins,xmin,xmax,"time [s]","Entries",kBlack,20);
  for(int bin=1; bin<=h_sinx->GetNbinsX(); ++bin){
    float bc = h_sinx->GetBinCenter(bin);
    h_sinx->SetBinContent(bin, f->Eval(bc));
  }

  // Plot and save this guy
  h_sinx->Draw();
  save << savedir << "sinx_freq" << freq << ".png";
  //c->SaveAs(save.str().c_str());
  save.str("");

  // Now get the fourier transform of this guy
  // in order to save it to a file
  TH1F* h_sinf = h_sinx->FFT(NULL,"MAG");

  // Scale this histogram to normal units... stupid root
  // puts it in weird units related to input histogram
  TH1F* h_sinf_scaled = remake(h_sinf,nbins,xmin,xmax);
  setHistAtt(h_sinf_scaled,"f [Hz]","Entries",kBlack,20);

  // Save
  c->SetLogx();
  h_sinf_scaled->Draw();
  save << savedir << "sinf_freq" << freq << ".png";
  //c->SaveAs(save.str().c_str());
  c->SetLogx(false);
  save.str("");

  // Now get the FFT Object for 
  // this partticular input
  TVirtualFFT* fft = TVirtualFFT::GetCurrentTransform();
  
  // Get the re and im components
  Double_t *re_full = new Double_t[nbins];
  Double_t *im_full = new Double_t[nbins];
  fft->GetPointsComplex(re_full,im_full);

  // Now lets make a backward transform
  TVirtualFFT* fft_back = TVirtualFFT::FFT(1, &nbins, "C2R M K");
  fft_back->SetPointsComplex(re_full,im_full);
  fft_back->Transform();
  TH1F* h_back = TH1::TransformHisto(fft_back,NULL,"Re");
  
  // Format reverse transform
  TH1F* h_back_format = remakeTF(h_back, nbins, xmin, xmax);
  setHistAtt(h_back_format, "time [s]", "Entries", kBlack, 20);
  
  // Save this result
  save << savedir << "sinx_back_freq" << freq << ".png";
  h_back_format->Draw();
  ////c->SaveAs(save.str().c_str());
  

}

//--------------------------------------------//
// Plot Sin
//--------------------------------------------//
void plotE(TFile* efile, TFile* gfile, 
	   TString ename, TString gname, 
	   int opt)
{
  
  // Make Canvas
  TCanvas* c = makeCanvas("c");
  c->SetRightMargin(0.10);
  
  // sstring to save
  stringstream save;

  // Load the gain histogram
  TH1F* h_gain = (TH1F*) gfile->Get(gname.Data());

  // First obtain the sum of the electric field
  TGraph* Et = (TGraph*) efile->Get(ename.Data());
  int nbins = 3000;
  float xmin = 0;
  float xmax = 100e-9;
  TH1F* h_Et = convertToSeconds(Et,nbins,xmin,xmax);
  setHistAtt(h_Et, "time [s]", "E(t) [V m^{-1}]", kBlack, 20);
  
  // Draw and save
  h_Et->Draw();
  save << savedir << ename << "_input.png";
  //c->SaveAs(save.str().c_str());
  save.str("");

  // Now fourier transform that shit
  TH1F* h_Ef = h_Et->FFT(NULL,"MAG");
  TH1F* h_Ef_format = remake(h_Ef,nbins,xmin,xmax);
  setHistAtt(h_Ef_format,"freq [Hz]", "E(f) [V m^{-1} Hz^{-1}]",kBlack,20);
  int low = h_Ef_format->FindBin(1e6);
  int high = h_Ef_format->FindBin(4000e6);
  h_Ef_format->GetXaxis()->SetRange(low,high);

  // Draw and save
  save << savedir << ename << "_FT.png";
  c->SetLogx();
  c->SetLogy();
  h_Ef_format->Draw();
  //c->SaveAs(save.str().c_str());
  save.str("");
  c->SetLogx(false);
  c->SetLogy(false);

  // Now get the FFT Object:
  TVirtualFFT* fft = TVirtualFFT::GetCurrentTransform();
  
  // Get the re and im components
  Double_t *re_full = new Double_t[nbins];
  Double_t *im_full = new Double_t[nbins];
  fft->GetPointsComplex(re_full,im_full);

  // This is where the options come into play.
  // opt == 1 -- Fourier transform the original back
  // opt == 2 -- Fourier transform specified range with flat gain
  // opt == 3 -- Fourier transform specified range with gain
  
  TString append = "_min0Hz";
  if( opt != 1 ){
    
    // Reset bins
    //float fmin = 100e6;
    //float fmax = 800e6;
    float fmin = 100e6;
    float fmax = 300e6;
    resetbins(re_full, im_full, nbins, fmin, fmax, 1/(xmax-xmin), h_Ef, h_gain,opt);

  }


  // Now lets make a backward transofmr
  TVirtualFFT* fft_back = TVirtualFFT::FFT(1, &nbins, "C2R M K");
  fft_back->SetPointsComplex(re_full,im_full);
  fft_back->Transform();
  TH1F* h_back = TH1::TransformHisto(fft_back,NULL,"Re");
  TH1F* h_back_format = remakeTF(h_back, nbins, xmin, xmax);
  if( opt == 1 )
    setHistAtt(h_back_format, "time [s]", "E(t) [V m^{1}]", kBlack, 20);
  else
    setHistAtt(h_back_format, "time [s]", "V(t) [V]", kBlack, 20);

  // Now save
  save << savedir << ename << "_IFT";
  if( opt == 1 )    save << "_backToOriginal.png";
  else if(opt == 2) save << "_flatGain" << append << ".png";
  else if(opt == 3) save << "_" << gname << ".png";
  h_back_format->SetMaximum(0.02);
  h_back_format->SetMinimum(-0.02);
  h_back_format->Draw();
  //c->SaveAs(save.str().c_str());

}

//--------------------------------------------//
// Plot SumE
//--------------------------------------------//
void plotSumE(TFile* efile, TFile* gfile, 
	      TString ename, TString gname, 
	      int opt)
{

  // Make Canvas
  TCanvas* c = makeCanvas("c");
  c->SetRightMargin(0.10);

  // sstring to save
  stringstream save;

  // Load the gain histogram
  TH1F* h_gain = (TH1F*) gfile->Get(gname.Data());

  // First obtain the sum of the electric field
  TH1F* h_Et = (TH1F*) efile->Get(ename.Data());
  int nbins = h_Et->GetNbinsX();
  float xmin = h_Et->GetBinLowEdge(1);
  float xmax = h_Et->GetBinLowEdge(nbins) + h_Et->GetBinWidth(nbins);
  cout<<"Min: "<<xmin<<" max: "<<xmax<<" nbins: "<<nbins<<" Step: "<<((xmax-xmin)/nbins)*1e9<<endl;
  setHistAtt(h_Et, "time [s]", "E(t) [V m^{-1}]", kBlack, 20);
  
  // Draw and save
  h_Et->Draw();
  save << savedir << ename << "_input.png";
  c->SaveAs(save.str().c_str());
  save.str("");

  // Now fourier transform that shit
  TH1F* h_Ef = h_Et->FFT(NULL,"MAG");
  TH1F* h_Ef_format = remake(h_Ef,nbins,xmin,xmax);
  setHistAtt(h_Ef_format,"freq [Hz]", "E(f) [V m^{-1} Hz^{-1}]",kBlack,20);
  int low = h_Ef_format->FindBin(1e6);
  int high = h_Ef_format->FindBin(4000e6);
  h_Ef_format->GetXaxis()->SetRange(low,high);

  // Draw and save
  save << savedir << ename << "_FT.png";
  c->SetLogx();
  c->SetLogy();
  h_Ef_format->Draw();
  c->SaveAs(save.str().c_str());
  save.str("");
  c->SetLogx(false);
  c->SetLogy(false);

  // Now get the FFT Object:
  TVirtualFFT* fft = TVirtualFFT::GetCurrentTransform();
  
  // Get the re and im components
  Double_t *re_full = new Double_t[nbins];
  Double_t *im_full = new Double_t[nbins];
  fft->GetPointsComplex(re_full,im_full);

  // This is where the options come into play.
  // opt == 1 -- Fourier transform the original back
  // opt == 2 -- Fourier transform specified range with flat gain
  // opt == 3 -- Fourier transform specified range with gain
  
  //float fmaxs[] = {800e6, 700e6, 600e6, 500e6, 400e6, 300e6, 200e6};
  //int nfmax = 7;

  if( opt != 1 ){

    // Reset bins
    float fmin = 100e6;
    float fmax = 800e6;
    resetbins(re_full, im_full, nbins, fmin, fmax, 1/(xmax-xmin), h_Ef, h_gain,opt);

  }


  // Now lets make a backward transofmr
  TVirtualFFT* fft_back = TVirtualFFT::FFT(1, &nbins, "C2R M K");
  fft_back->SetPointsComplex(re_full,im_full);
  fft_back->Transform();
  TH1F* h_back = TH1::TransformHisto(fft_back,NULL,"Re");
  TH1F* h_back_format = remakeTF(h_back, nbins, xmin, xmax);
  if( opt == 1 )
    setHistAtt(h_back_format, "time [s]", "E(t) [V m^{1}]", kBlack, 20);
  else
    setHistAtt(h_back_format, "time [s]", "V(t) [V]", kBlack, 20);

  // Now save
  save << savedir << ename << "_IFT";
  if( opt == 1 )    save << "_backToOriginal.png";
  else if(opt == 2) save << "_flatGain_fmax"<<fmax<<".png";
  else if(opt == 3) save << "_" << gname << ".png";
  h_back_format->Draw();
  c->SaveAs(save.str().c_str());

  TFile* fout = new TFile(("outroot/"+ename+".root").Data(),"recreate");
  h_back_format->Write();
  fout->Write();
  fout->Close();

}

//--------------------------------------------//
// reset bins
//--------------------------------------------//
void resetbins(Double_t* re, Double_t* im,
	       int nbins, float fmin, float fmax,
	       float scale, TH1F* h, TH1F* gain = NULL,
	       int opt = 2)
{

  // Loop over the bins
  for(int n=0; n<nbins; ++n){
    
    // Get the frequency with appropriate scale
    // to take it back to Hz
    double f = h->GetBinCenter(n+1) * scale;

    // Check if it is in range we are after
    // otherwise set it to zero
    if( fmin < f && f < fmax){
      
      // if it is in range, see if we 
      // should use the gain
      if( opt == 3 ){
	int bin = gain->FindBin( f/1.e6 );
	double g = gain->GetBinContent(bin);
	re[n] = g * re[n];
	im[n] = g * im[n];
      }
      else; // do nothing
    }
    else{
      re[n] = 0;
      im[n] = 0;
    }
  }// end loop over points
  
}

//--------------------------------------------//
// Remake histogram
//--------------------------------------------//
TH1F* remakeTF(TH1F* hin, int nbins, float xmin, float xmax)
{

  int nbins = hin->GetNbinsX();
  
  TH1F* h = new TH1F("scaledTF","",nbins, xmin, xmax);
  double shift = (xmax-xmin)/nbins;
  for(int bin=1; bin<=nbins; ++bin){
    float bc = hin->GetBinContent(bin)/nbins;
    float bcent = hin->GetBinCenter(bin);
    
    h->Fill(bcent*shift,bc);
  }

  return h;
}

//--------------------------------------------//
// Remake histogram
//--------------------------------------------//
TH1F* remake(TH1F* hin, int nbins, float xmin, float xmax)
{

  float xminin = hin->GetXaxis()->GetXmin();
  float xmaxin = hin->GetXaxis()->GetXmax();
  int nbins = hin->GetNbinsX();

  cout<<"Using Nbins: "<<nbins<<" min: "<<xmin<<" max: "<<xmax<<" diff: "<<1/(xmax-xmin)<<endl;

  TH1F* h = new TH1F("scaled","",nbins, xminin/(xmax-xmin), xmaxin/(xmax-xmin));
  for(int bin=1; bin<=nbins; ++bin){
    float bc = hin->GetBinContent(bin)/sqrt(nbins);
    float bcent = hin->GetBinCenter(bin);

    h->Fill(bcent/(xmax-xmin),bc);
  }

  return h;
}

//--------------------------------------------//
// Convert to bins in seconds
//--------------------------------------------//
TH1F* convertToSeconds(TGraph* &gr, int nbins,
		       float xmin, float xmax)
{

  TH1F* h = new TH1F("h","",nbins,xmin,xmax);

  int np = gr->GetN();
  double E=0, t=0;
  for(int n=0; n<np; ++n){
    gr->GetPoint(n,t,E);
    h->Fill(t*1e-9,E);
  }

  return h;
}

//--------------------------------------------//
// Get real part
//--------------------------------------------//
TGraph* getReal(TGraph* Mag, TGraph* Phase)
{

  int npoints = 0;
  int np = Mag->GetN();
  double mag=0, phase=0, x=0;
  double Real[10000];
  double var[10000];
  for(int n=0; n<np; ++n){
    Mag->GetPoint(n,x,mag);
    Phase->GetPoint(n,x,phase);
    Real[n] = mag * cos(phase);
    var[n] = x;
    npoints++;
  }

  return new TGraph(npoints,var,Real);

}

//--------------------------------------------//
// Get the potential
//--------------------------------------------//
TGraph* getVF(TGraph* gr_E, TH1F* h_gain)
{
  
  // Placeholders for potential
  double Voltage[10000];
  double freq[10000];
  int npoints = 0;

  // Get each point in electric field
  // and multiply by the gain
  int np = gr_E->GetN();
  double E=0, f=0;
  for(int n=0; n<np; ++n){

    // Get points and gain that
    // is associated with frequency
    gr_E->GetPoint(n,f,E);
    int bin = h_gain->FindBin(f);

    // Save
    freq[n] = f;
    //Voltage[n] = E * h_gain->GetBinCenter(bin);
    Voltage[n] = E;
    npoints++;

  }// end loop over bins

  // Make graph and return
  return new TGraph(npoints, freq, Voltage);

}

 

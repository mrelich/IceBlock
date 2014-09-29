
#include "TH1.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TFile.h"

#include <vector>

//------------------------------------------------------------//
// Make Histogram
//------------------------------------------------------------//
TH1F* makeHist(TString name, int nbins, float xmin, float xmax,
               TString xtitle, TString ytitle, int color, int marker)
{

  TH1F* h = new TH1F(name.Data(),name.Data(),nbins,xmin,xmax);
  setHistAtt(h,xtitle,ytitle,color,marker);
  return h;

}

//------------------------------------------------------------//
// Set histogram attributes
//------------------------------------------------------------//
void setHistAtt(TH1F* &h, TString xtitle, TString ytitle, 
                int color, int marker)
{

  h->GetXaxis()->SetTitle(xtitle.Data());
  h->GetYaxis()->SetTitle(ytitle.Data());
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetMarkerStyle(marker);
  h->SetStats(0);
  h->SetTitle("");
  h->SetLineWidth(2);
  h->GetYaxis()->SetTitleOffset(1.5);

}


//------------------------------------------------------------//
// Make pads, needed for ratio plotting
//------------------------------------------------------------//
void makePads(TCanvas* &c, TPad* &top, TPad* &bot, float ysep=0.3)
{

  c->cd();
  top = new TPad("top","",0,ysep,1,1);
  top->SetTopMargin(0.05);
  top->SetBottomMargin(0.02);
  top->SetRightMargin(0.05);
  top->SetLeftMargin(0.12);
  top->SetNumber(1);

  bot =new TPad("bot","",0,0,1,ysep);
  bot->SetTopMargin(0.015);
  bot->SetBottomMargin(0.3);
  bot->SetRightMargin(0.05);
  bot->SetLeftMargin(0.12);
  bot->SetNumber(2);

}



//------------------------------------------------------------//
// Fit for radiation length
//------------------------------------------------------------//
pair<float,float> getRadLength(TProfile* prof, TF1* &fit, int color)
{
  
  // Must pass the distance travelled vs the average energy.
  // Fit for the radiation length following:
  // <E> = E0 e^{-x/chi_0)

  //TF1 func = TF1("radLength","[0]*TMath::Exp(-x/[1]) + [2]",0,200);
  TF1 func = TF1("radLength","[0]*TMath::Exp(-x/[1])",0,200);
  func.SetParameter(0,0.5);
  func.SetParameter(1,4.0);

  
  // Fit
  prof->Fit("radLength","RQ");

  fit = prof->GetFunction("radLength");
  fit->SetLineColor(color);

  // Print info
  float E0 = fit->GetParameter(0);
  float radLength = fit->GetParameter(1);
  float radErr = fit->GetParError(1);

  cout<<"E0: "<<E0<<" chi: "<<radLength<<" +/- "<<radErr<<endl;

  //delete fit;

  pair<float,float> radPair (radLength, radErr);

  return radPair;

}

//------------------------------------------------------------//
// Make TLatex
//------------------------------------------------------------//
TLatex* makeLatex()
{

  TLatex* lat = new TLatex();
  lat->SetTextSize(0.04);
  lat->SetTextFont(42);
  lat->SetTextColor(kBlack);
  lat->SetNDC();
  return lat;
}

//------------------------------------------------------------//
// Get maximum of two histograms
//------------------------------------------------------------//
double getMax(TH1* h0, TH1* h1)
{
  
  if(h0->GetMaximum() < h1->GetMaximum())
    return h1->GetMaximum();
  
  return h0->GetMaximum();

}


//------------------------------------------------------------//
// Make legend
//------------------------------------------------------------//
TLegend* makeLegend(float x0, float x1, float y0, float y1)
{

  TLegend* leg = new TLegend(x0, y0, x1, y1);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetLineColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  return leg;

}

//------------------------------------------------------------//
// Make canvas
//------------------------------------------------------------//
TCanvas* makeCanvas(TString name)
{

  TCanvas* c = new TCanvas(name.Data(), name.Data(), 700, 500);
  c->SetTicks(1,1);
  c->SetTopMargin(0.05);
  c->SetRightMargin(0.05);
  c->SetBottomMargin(0.12);
  c->SetLeftMargin(0.12);
  
  return c;
}

//------------------------------------------------------------//
// Get TProfile
//------------------------------------------------------------//
TProfile* getProfile(TFile* file, TString pname, TString xtitle,
		     TString ytitle, int color, int marker)
{
  
  TProfile* prof = (TProfile*) (file->Get(pname.Data())->Clone(Form("%s_%i",pname.Data(),color)));
  prof->GetXaxis()->SetTitle(xtitle.Data());
  prof->GetYaxis()->SetTitle(ytitle.Data());
  prof->SetMarkerStyle(marker);
  prof->SetMarkerColor(color);
  //prof->SetMarkerSize(0.5);
  prof->SetLineColor(color);
  prof->SetTitle("");
  prof->SetStats(0);
  prof->GetYaxis()->SetTitleOffset(1.5);
  
  return prof;

}

//------------------------------------------------------------//
// Get histogram
//------------------------------------------------------------//
TH1F* getHist(TFile* file, TString pname, TString xtitle,
		  TString ytitle, int color, int marker)
{
  
  TH1F* hist = (TH1F*) (file->Get(pname.Data())->Clone(Form("%s_%i",pname.Data(),color)));
  hist->GetXaxis()->SetTitle(xtitle.Data());
  hist->GetYaxis()->SetTitle(ytitle.Data());
  hist->SetMarkerStyle(marker);
  hist->SetMarkerColor(color);
  hist->SetMarkerSize(0.5);
  hist->SetLineColor(color);
  hist->SetTitle("");
  hist->SetStats(0);
  hist->GetYaxis()->SetTitleOffset(1.5);
  hist->SetLineWidth(2);
  return hist;

}

//------------------------------------------------------------//
// Get 2D histogram
//------------------------------------------------------------//
TH2F* getHist2(TFile* file, TString pname, TString xtitle,
	       TString ytitle, TString ztitle)
{
  
  TH2F* hist = (TH2F*) (file->Get(pname.Data())->Clone());
  hist->GetXaxis()->SetTitle(xtitle.Data());
  hist->GetYaxis()->SetTitle(ytitle.Data());
  hist->GetZaxis()->SetTitle(ztitle.Data());
  hist->SetTitle("");
  hist->SetStats(0);
  hist->GetYaxis()->SetTitleOffset(1.5);
  hist->GetZaxis()->SetTitleOffset(1.);
  hist->SetLineWidth(2);

  return hist;

}

//------------------------------------------------------------//
// Set graph attributes
//------------------------------------------------------------//
void setAtt(TGraph* &gr, TString xtitle, TString ytitle, int color)
{

  gr->GetXaxis()->SetTitle(xtitle.Data());
  gr->GetYaxis()->SetTitle(ytitle.Data());
  gr->SetLineColor(color);
  gr->SetMarkerColor(color);
  gr->SetTitle("");

}

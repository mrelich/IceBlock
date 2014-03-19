
#include "TH1.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TFile.h"

#include <vector>

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
// Plot to fine Moliere Radius
//------------------------------------------------------------//
void plotMoliereRadius()
{

  // The Moliere radius can be estimated by looking at the
  // ratio of the shower energy outside some radius R 
  // of a hypothetical cylinder centered on the shower max.
  // In practice, we look at the average showere energy
  // outside this radius in steps.  Then look when U/E(0) = 0.1
  // That is Moliere radius.

  // The plot name
  TString pname = "rVsShowerE";

  // Constant line
  TF1* f1 = new TF1("f1","0.1",0,100);
  f1->SetLineColor(kBlack);
  f1->SetLineStyle(2);

  // Make canvas
  TCanvas* c = makeCanvas("c");
  c->SetLogx();
  c->SetLogy();

  // Make Legend
  TLegend* leg = makeLegend(0.6, 0.8, 0.7, 0.9);

  // Loop and get profiles
  TProfile* profs[nFiles];
  for(int f=0; f<nFiles; ++f){
    profs[f] = getProfile(m_files[f], pname, "Radius [cm]", "U/E_{0}",
			  m_colors[f], m_markers[f]);
    leg->AddEntry(profs[f], m_labels[f].Data(),"l");
    
    // Just have some uniform min and max
    profs[f]->SetMinimum(1e-2);
    profs[f]->SetMaximum(1);
    int xmin = profs[f]->FindBin(1);
    int xmax = profs[f]->FindBin(50);
    profs[f]->GetXaxis()->SetRange(xmin,xmax);

    // draw
    if( f==0 ) profs[f]->Draw();
    else       profs[f]->Draw("same");
  }

  // draw legend
  leg->Draw("same");

  // Draw function
  f1->Draw("same");

  // Save
  c->SaveAs("/home/mrelich/workarea/ara/ELS_ana/plots/MoliereRadius.png");

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
  prof->SetMarkerSize(0.5);
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


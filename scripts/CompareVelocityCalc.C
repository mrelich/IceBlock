
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-//
// This script will compare the vector potential from ZHS to //
// what is calculated using Geant4.                          //
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-//

#include "myHist.C"

// Scale factor is needed.  Since ZHS gives R*A, and G4 gives 
// just A, we need to scale by the distance to antenna.  So
// in this case it is 100m
double m_R = 1000; //100.;
//double m_R = 100.; //100.;
//double m_R = 400.;
//double m_R = 1000.;

bool m_save = false;

TString m_savedir = "plots/EField/VelocityCalcComp/";

//----------------------------------------//
// Main
//----------------------------------------//
void CompareVelocityCalc(bool save = false)
{

  m_save = save;

  // Specify the two input files
  TString indir = "efieldroot/";
  TString f_oldV = indir + "Output_50Evt_100GeV_1Prim_HardCodedAntenna_R1000m_oldV.root";
  TString f_newV = indir + "Output_50Evt_100GeV_1Prim_HardCodedAntenna_R1000m_newV.root";

  // Specify pname
  TString p_g4name  = "A_AntNum_0_pos_827.371_0_561.656";


  // Now plot
  plotWithRatio(f_oldV, f_newV, p_g4name);
  
}

//----------------------------------------//
// Plot ZHS and G4 on same figure
//----------------------------------------//
void plotWithRatio(TString f_oldV,
		   TString f_newV,
		   TString pname)
{

  // Make canvas
  TCanvas* c = makeCanvas("c");

  // Make a legend
  TLegend* leg = makeLegend(0.15,0.3,0.79,0.92);  
  leg->SetTextSize(0.055);

  // Titles
  TString xtitle = "time [ns]";
  TString ytitle = "A [Vs/m]";    

  // Loop and get plots
  TProfile* p_oldV = NULL;
  TProfile* p_newV = NULL;
  TH1D* h_resetZ   = NULL;

  // Load ZHS profile
  TFile* file_old = new TFile(f_oldV);
  p_oldV = getProfile(file_old, pname, xtitle, ytitle, kBlue, 20);
  p_oldV->SetDirectory(0);
  file_old->Close();
  leg->AddEntry(p_oldV,"Old","lep");

  // Load Geant profile
  TFile* file_new = new TFile(f_newV);
  p_newV = getProfile(file_new, pname, xtitle, ytitle, kRed, 20);
  p_newV->SetDirectory(0);
  file_new->Close();
  leg->AddEntry(p_newV,"New","lep");
  
  // Scale Geant4 profile up by factor of R
  //p_newV->Scale(m_R);
  
  // Get maximum
  float maximum = p_newV->GetMaximum();
  if( maximum < p_oldV->GetMaximum() )
    maximum = p_oldV->GetMaximum();
  p_newV->SetMaximum(maximum*5);
  p_newV->SetMinimum(maximum*1e-5);
  
  // Make two pads
  TPad* top = NULL;
  TPad* bot = NULL;
  makePads(c, top, bot);
  
  // Set some attributes
  p_newV->GetYaxis()->SetLabelSize(0.05);
  p_newV->GetYaxis()->SetTitleSize(0.05);
  p_newV->GetYaxis()->SetTitleOffset(1.0);
  p_newV->GetXaxis()->SetLabelSize(0);
  
  // Draw top plots
  c->cd();
  top->Draw();
  top->cd();
  top->SetLogy();
  p_newV->Draw();
  p_oldV->Draw("sameep");
  leg->Draw("same");
  top->Update();
  
  // Draw bot ratio
  c->cd();
  bot->Draw();
  bot->cd();
  TH1D* prof = p_newV->ProjectionX("ratio");
  prof->SetLineColor(p_newV->GetLineColor());
  prof->SetMarkerColor(p_newV->GetMarkerColor());
  prof->SetStats(0);
  prof->Divide(p_oldV->ProjectionX("old_px"));
  prof->GetYaxis()->SetTitle("New/Old");
  prof->SetMinimum(0);
  prof->SetMaximum(2);
  prof->GetXaxis()->SetLabelSize(0.1);
  prof->GetYaxis()->SetLabelSize(0.1);
  prof->GetXaxis()->SetTitleSize(0.12);
  prof->GetYaxis()->SetTitleSize(0.12);
  prof->GetYaxis()->SetTitleOffset(0.4);
  prof->GetYaxis()->SetNdivisions(405);
  
  // Fit a line to ratio
  TF1* line = new TF1("line","[0]",0,1000000);
  line->SetLineColor(kBlack);
  prof->Fit("line");

  // Draw bot
  prof->Draw();
  line->Draw("same");

  // Plot fit results
  TLatex* lat = makeLatex();
  lat->SetTextSize(0.15);
  lat->DrawLatex(0.15,0.8,Form("Fit = %1.3f",line->GetParameter(0)));
  //lat->DrawLatex(0.2,0.8,"test");

  // Update bot
  bot->Update();


  
  if(m_save){
    //TString save = m_savedir + p_g4 + "_ZHSCherAngle.png";
    TString save = m_savedir + "comparing_old_new_Vcalc.png";
    c->SaveAs(save.Data());
  }
  

}


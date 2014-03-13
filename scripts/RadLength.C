
#include "myHist.C"

TString energy = "";
TString eLbl   = "";
TString material = "";

//----------------------------------------//
// Main
//----------------------------------------//
void RadLength(int eopt, int mopt)
{

  // Plot name
  TString pname = "dVsE";

  // Specify the energy to fit for
  if(eopt == 0){
    energy = "1000000";
    eLbl   = "1 TeV";
  }
  else if(eopt == 1){
    energy = "100000";
    eLbl   = "100 GeV";
  }
  else if(eopt == 2){
    energy = "10000";
    eLbl   = "10 GeV";
  }
  else if(eopt == 3){
    energy = "1000";
    eLbl   = "1 GeV";
  }
  else{
    cout<<"Energy option is not supported"<<endl;
    return;
  }
  
  // Specify material
  if(mopt == 0)      material = "ice";
  else{
    cout<<"Material not supported"<<endl;
    return;
  }
    
  // Open file
  TString fname = "rootfiles/StepAna_100_";
  fname += energy + "_"+material+".root";
  TFile* m_file = new TFile(fname.Data());
  m_file->Print();

  // Make Canvas
  TCanvas* c = makeCanvas("c");
  
  // Get Legned
  TLegend* leg = makeLegend(0.7,0.9,0.7,0.9);
  leg->SetHeader(("E_{beam} = "+eLbl).Data());

  // Get the profile
  TProfile* prof = getProfile(m_file, pname, "Distance [cm]",
			      "Energy [GeV]", kBlack, 20);
  leg->AddEntry(prof, "Data", "lep");
  
  // Get Fit
  TF1* fit = getFit(prof);
  leg->AddEntry(fit, "Fit", "l");
  
  // Draw
  prof->Draw();
  fit->Draw("same");
  leg->Draw("same");

  // Add Fit results
  TLatex* lat = makeLatex();
  lat->SetTextSize(0.06);
  lat->DrawLatex(0.4,0.8,
		 Form("#chi_{0} = %4.2f",fit->GetParameter(1))
		 );

  c->SaveAs(("plots/RadLength/energy_"+energy+"_"+material+".png").Data());

}

//----------------------------------------//
// Fit function
//----------------------------------------//
TF1* getFit(TProfile* prof)
{

  // Set fit function
  TF1 func = TF1("radLength","[0]*TMath::Exp(-x/[1])",0,200);
  func.SetParameter(0,0.5);
  func.SetParameter(1,4.0);

  // Fit
  prof->Fit("radLength","RQ");
  TF1* fit = prof->GetFunction("radLength");
  fit->SetLineColor(kBlue);

  return fit;
}

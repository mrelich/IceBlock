
#include "myHist.C"
#include <string>
#include <sstream>

int colors[] = {kBlack, kBlue, kRed};
int markers[] = {20, 25, 23};

float E_critical = 79.0255; // From Geant4 (MeV)
//float E_critical = 68.8; // From paper
//float E_critical = 50.; // guess

//-------------------------------------------//
// Main
//-------------------------------------------//
void NParticle()
{

  basic();
  //singleFit();
  //multipleEnergies();
}

//-------------------------------------------//
// Basic N(e+p) and N(e-p) plot
//-------------------------------------------//
void basic()
{
  
  // Specify the beam energy
  TString energy = "100000"; // MeV
  TString eLbl   = "100 GeV";

  // Get the file
  TString fname = "rootfiles/TrkAna_100_";
  //TString fname = "rootfiles/TrkAna_50_";
  fname += energy;
  fname += "_ice_eBeam.root";
  TFile* file = new TFile(fname.Data());

  // Plots
  vector<TString> pnames;
  pnames.push_back("NPartSum");
  pnames.push_back("NPartDiff");

  // Names
  vector<TString> names;
  names.push_back("N(e+p)");
  names.push_back("N(e-p)");

  // Make canvas
  TCanvas* c = makeCanvas("c");
  
  // Make legend
  TLegend* leg = makeLegend(0.6,0.7,0.9,0.7);
  leg->SetHeader(("E_{beam} = " + eLbl).Data());

  // Loop and get plots
  TProfile* profs[2];
  float maximum = -999;
  for(unsigned int i=0; i<pnames.size(); ++i){
    TString pname = pnames.at(i);
    profs[i] = getProfile(file,pname,"Radiation Length",
			  "Number of Particles", colors[i],
			  markers[i]);
    
    if( maximum < profs[i]->GetMaximum())
      maximum = profs[i]->GetMaximum();
    
    leg->AddEntry(profs[i], names.at(i).Data(), "lep");
  }

  // Now draw
  profs[0]->SetMaximum(1.2*maximum);
  profs[0]->Draw();
  for(unsigned int i=0; i<pnames.size(); ++i)
    profs[i]->Draw("same");
  leg->Draw("same");

  c->SaveAs("plots/NParticles/basic_100_100000_eBeam.png");

}

//-------------------------------------------//
// Basic N(e+p) w/Fit
//-------------------------------------------//
void singleFit()
{

  //TString beam = "eBeam";
  TString beam = "gBeam";
  
  // Energies
  vector<TString> energies;
  energies.push_back("100000");
  energies.push_back("500000");
  energies.push_back("1000000");
  float E[] = {100000,500000,1000000};
  
  // E Labels
  vector<TString> eLbls;
  eLbls.push_back("100 GeV");
  eLbls.push_back("500 GeV");
  eLbls.push_back("1 TeV");

  // NEvents
  vector<TString> nEvts;
  nEvts.push_back("100");
  nEvts.push_back("50");
  nEvts.push_back("20");


  // Plots
  TString pname = "NPartSum";
  TString name  = "N(e+p)";
  
  // Make canvas
  TCanvas* c = makeCanvas("c");
  
  // Make legend
  TLegend* leg = makeLegend(0.7,0.8,0.9,0.7);

  // Latex object
  //TLatex* lat = makeLatex();
  TLatex* lat = new TLatex();
  lat->SetNDC();
  lat->SetTextSize(0.04);

  for(unsigned int i=0; i<energies.size(); ++i){
    TString energy = energies.at(i);
    TString eLbl   = eLbls.at(i);
    TString nEvt   = nEvts.at(i);
    float e        = E[i];

    TString fname = "rootfiles/TrkAna_"+nEvt+"_"+energy+"_ice_"+beam+"_5MeV.root";
    TFile* file = new TFile(fname.Data());

    leg->Clear();
    leg->SetHeader(("E_{beam} = " + eLbl).Data());

    // Get Plot
    TProfile* prof = getProfile(file,pname,"Radiation Length",
				"Number of Particles",kBlack,20);
    //prof->Scale(0.65);
    TF1* fit = fitGreisen(prof,e,kBlue,2);
    leg->AddEntry(prof, name.Data(), "lep");
    leg->AddEntry(fit, "Fit", "l");
        
    // Now draw
    //prof->SetMaximum(1.2*prof->GetMaximum());
    prof->Draw();
    fit->Draw("same");
    leg->Draw("same");
    
    // Add some fit info
    lat->DrawLatex(0.73,0.66,Form("A(E) = %1.2f",fit->GetParameter(0)));
    lat->DrawLatex(0.73,0.60,Form("a(E) = %1.2f",fit->GetParameter(1)));

    cout<<"Energy: "<<eLbl
	<<" A(E): "<<fit->GetParameter(0)
	<<" a(E): "<<fit->GetParameter(1)
	<<endl;

    TString save = "plots/NParticles/singleFits_"+nEvt+"_"+energy+"_"+beam+".png";
    c->SaveAs(save.Data());

  }// end loop over energies

}

//-------------------------------------------//
// Draw N(e+p) for many showers with
// Greison fits
//-------------------------------------------//
void multipleEnergies()
{
  
  // Specify the files
  const int nFiles = 3;
  TFile* files[nFiles];
  TString fbase = "rootfiles/TrkAna";
  TString beam  = "eBeam";
  //TString beam  = "gBeam";
  files[0] = new TFile((fbase+"_100_100000_ice_"+beam+"_5MeV.root"));
  files[1] = new TFile((fbase+"_50_500000_ice_"+beam+"_5MeV.root"));
  files[2] = new TFile((fbase+"_20_1000000_ice_"+beam+"_5MeV.root"));


  float energy[] = {100000,500000,1000000};

  // Plots
  TString pname = "NPartSum";

  // Names
  vector<TString> names;
  names.push_back("100 GeV");
  names.push_back("500 GeV");
  names.push_back("1 TeV");

  // Make canvas
  TCanvas* c = makeCanvas("c");
  c->SetLogy();
  
  // Make legend
  TLegend* leg   = makeLegend(0.2,0.4,0.2,0.4);
  TLegend* lFits = makeLegend(0.4,0.6,0.2,0.4);
  
  // Loop and get plots
  TProfile* profs[nFiles];
  TF1* fits[nFiles];
  for(int i=0; i<nFiles; ++i){
    profs[i] = getProfile(files[i],pname,"Shower Depth [#chi_{0}]",
			  "N(e+p)", colors[i],
			  markers[i]);
    leg->AddEntry(profs[i], names.at(i).Data(), "lep");

    // Get Fit
    fits[i] = fitGreisen(profs[i], energy[i], colors[i], 2);
    float A = fits[i]->GetParameter(0);
    float a = fits[i]->GetParameter(1);
    lFits->AddEntry(fits[i], 
		    Form("Fit: A(E)=%1.2f, a(E)=%1.2f",A,a),
		    "lep");
  }
  
  // Now draw
  profs[0]->SetMinimum(0.1);
  profs[0]->SetMaximum(2000);
  profs[0]->Draw();
  for(int i=1; i<nFiles; ++i){
    profs[i]->Draw("same");
    //fits[i]->Draw("same");
  }
  
  leg->Draw("same");
  lFits->Draw("same");

  // Save
  TString save = "plots/NParticles/ShowerFits_"+beam+".png";
  c->SaveAs(save.Data());
}

//-------------------------------------------//
// Fit Greisen parameterization for shower
//-------------------------------------------//
TF1* fitGreisen(TProfile* prof, float E0, int color, int style)
{

  // Function 7 from this paper:
  //    * http://prd.aps.org/pdf/PRD/v65/i10/e103002
  // Function of shower energy, energy, and shower 
  // depth.

  // [0] -- A(E)
  // [1] -- a(E)

  // y term
  stringstream y;
  y << "TMath::Log("
    << E0 << "/" << E_critical
    << ")";

  // Constant Term
  stringstream C;
  C << "(0.31*[0])/(TMath::Sqrt("
    << y.str() << "))";


  // Define t1
  string t1 = "(x+[1])";  
    
  // Ln(s1)
  stringstream lns1;
  lns1 << "TMath::Log(3*"
       << t1
       << "/(" 
       << t1 << "+2*" << y.str()
       <<"))";

  // expone
  stringstream exp;
  exp << "TMath::Exp("
      << t1 <<"*(1-1.5*"
      << lns1.str()
      <<"))";
  
  //cout <<"Function: "<<C.str()+"*"+exp.str()<<endl;

  // Function
  TF1 func = TF1("greisen",(C.str()+"*"+exp.str()).c_str(),0,20);
  func.SetParameter(0,0.5);
  func.SetParameter(1,1.);
  //func.FixParameter(1,0.76);

  // Fit
  prof->Fit("greisen","RQ");
  TF1* fit = prof->GetFunction("greisen");
  fit->SetLineColor(color);
  fit->SetLineStyle(style);

  return fit;
  
}

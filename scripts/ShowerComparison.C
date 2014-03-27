
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
// If we want to use Geant we will need to simulate //
// many particles in a given shower to effectively  //
// replicate the high energy showers. This is both  //
// useful for ARA and the ELS test stand.           //
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//

#include "myHist.C"

TString savedir = "plots/ShowerComp/";

//------------------------------------------//
// Main
//------------------------------------------//
void ShowerComparison(int opt)
{

  // Specify the files to compare
  TString indir = "rootfiles/";
  vector<TString> fnames;
  vector<TString> names;
  TString saveAppend = "";

  // 100 GeV option
  if(opt == 0){
    fnames.push_back(indir+"TrkAna_100_100000_ice_eBeam_1MeV.root");
    fnames.push_back(indir+"TrkAna_100_10000_ice_eBeam_np10_1MeV.root");
  
    names.push_back("1 Particle @ 100GeV");
    names.push_back("10 Particles @ 10GeV");

    saveAppend = "100GeV";
  }
  else if(opt == 1){
    fnames.push_back(indir+"TrkAna_50_1000000_ice_eBeam_1MeV.root");
    fnames.push_back(indir+"TrkAna_50_100000_ice_eBeam_np10_1MeV.root");
    fnames.push_back(indir+"TrkAna_50_10000_ice_eBeam_np100_1MeV.root");
  
    names.push_back("1 Particle @ 1TeV");
    names.push_back("10 Particles @ 100 GeV");
    names.push_back("100 Particles @ 10 GeV");

    saveAppend = "1TeV";
  }
  else{
    cout<<"Option not supported"<<endl;
    return;
  }

  // Plot
  plotComp(fnames, names, saveAppend);

}

//------------------------------------------//
// Plot showers for given set of parameters
//------------------------------------------//
void plotComp(vector<TString> fnames, 
	      vector<TString> names,
	      TString saveAppend)
{
  
  // There are two plots we are after from 
  // each file.
  vector<TString> plots;
  plots.push_back("NPartSum");
  plots.push_back("NPartDiff");
  
  vector<TString> pnames;
  pnames.push_back("N(e+p)");
  pnames.push_back("N(e-p)");

  // set colors and markers
  int colors[] = {kBlack, kBlue, kRed};
  int markers[] = {20,25};
  
  // Make Canvas
  TCanvas* c = makeCanvas("c");

  // Make Legend
  TLegend* leg = makeLegend(0.5,0.8,0.7,0.9);

  // Set xandy title
  TString xtitle = "Radiation Lengths";
  TString ytitle = "Number of Particles";

  // Loop and Plot
  TProfile* prof[4][4];
  float maximum = -999;
  for(unsigned int f=0; f<fnames.size(); ++f){
    TString fname = fnames.at(f);
    TString name  = names.at(f);
    
    TFile* file = new TFile(fname.Data());
    int color = colors[f];

    // loop over plots and get profiles
    for(unsigned int i=0; i<plots.size(); ++i){
      TString plot = plots.at(i);
      TString pname = pnames.at(i);
      int marker = markers[i];

      prof[f][i] = getProfile(file,plot,xtitle,ytitle,color,marker);
      prof[f][i]->SetMarkerSize(1);
      leg->AddEntry(prof[f][i],(pname+" "+name).Data(),"lep");
      
      if( maximum < prof[f][i]->GetMaximum() )
	maximum = prof[f][i]->GetMaximum();
    }

  }
  
  prof[0][0]->SetMaximum(1.2*maximum);
  prof[0][0]->Draw();
  for(unsigned int f=0; f<fnames.size(); ++f)
    for(unsigned int i=0; i<plots.size(); ++i)
      if( !(i==0&&f==0) ) prof[f][i]->Draw("same");
  
  leg->Draw("same");
  
  c->SaveAs((savedir+"ChargeProfile_"+saveAppend+".png").Data());

}

//------------------------------------------//
// Calculate the shower integral
//------------------------------------------//


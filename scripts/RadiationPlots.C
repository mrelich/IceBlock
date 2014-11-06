
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-//
// This is a quick script trying to estimate the gammas that //
// are escaping the ice and staying in the ice.              //
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-//

#include "myHist.C"

TString m_savedir = "plots/Radiation/";
//bool m_save = false;
bool m_save = true;

//---------------------------------------------//
// Main
//---------------------------------------------//
void RadiationPlots(int opt)
{
  
  //TFile* file = new TFile("rootfiles/tracks/test.root");
  //TFile* file = new TFile("rootfiles/tracks/testOld.root");
  //TFile* file = new TFile("rootfiles/tracks/Radiation_40MeV_10kEvents.root");
  TFile* file = new TFile("rootfiles/tracks/Radiation_40MeV_100kEvents.root");
  
  if(opt == 0) plotVsEnergy(file);
  if(opt == 1) plotVsDistance(file);
  if(opt == 2) plotEvsDistance(file);
  if(opt == 3) plotEVariousDistance(file);
  if(opt == 4) plotTotEvsDistance(file);

}

//---------------------------------------------//
// Plot various results
//---------------------------------------------//
void plotVsEnergy(TFile* file)
{


  TCanvas* c = makeCanvas("c");
  c->SetLogy();

  vector<TString> pnames;
  pnames.push_back("p_Ein");
  pnames.push_back("p_Eout");
  pnames.push_back("p_Eleave");
  //pnames.push_back("p_Etop");
  //pnames.push_back("p_Eside");
  pnames.push_back("p_EMadeOut");
  

  vector<TString> legnames;
  legnames.push_back("Total In Ice");
  legnames.push_back("Total Outside");
  legnames.push_back("Total Leaving");
  //legnames.push_back("Leaving Top");
  //legnames.push_back("Leaving Side");
  legnames.push_back("Made Outside");
  

  TLegend* leg = makeLegend(0.6,0.75,0.65,0.93);

  int colors[] = {kBlack, kBlue, kMagenta, kRed, kGreen, kOrange};
  int markers[] = {20,25,22,23,21,23};

  TString xtitle = "Energy [MeV]";
  TString ytitle = "Photons / Primary";

  for(unsigned int i=0; i<pnames.size(); ++i){
  
    TH1F* h = getHist(file,pnames.at(i),xtitle,ytitle,colors[i],markers[i]);
    h->SetMarkerSize(1);
    if(i==0) h->Draw();
    else h->Draw("same");

    leg->AddEntry(h,legnames.at(i),"lep");
    
    int b15MeV = h->GetXaxis()->FindBin(15);
    cout<<legnames.at(i)
	<<"\tGreater than 15 MeV: "<<h->Integral(b15MeV,-1)<<" Photons / Primary"<<endl;
    
    float total = 0;
    for(int bin=1; bin<=h->GetNbinsX(); ++bin)
      total += h->GetBinContent(bin) * h->GetBinCenter(bin);
    cout<<"\tTotal Energy / primary: "<<total<<" MeV or "
    	<<Form("%.2f",(100*total/40))<<"%"<<endl;
    
    
  }
  leg->Draw("same");

  if( m_save )
    c->SaveAs((m_savedir+"RadiationEstimate.png"));


}

//---------------------------------------------//
// Plot vs distance
//---------------------------------------------//
void plotVsDistance(TFile* file)
{


  TCanvas* c = makeCanvas("c");
  c->SetLogy();

  TString pname = "p_cross";
  int color = kBlack;

  TString xtitle = "Distance to human [MeV]";
  TString ytitle = "Photons / Primary";

  TH1F* h = getHist(file,pname,xtitle,ytitle,color,20);
  h->SetMarkerSize(1);
  h->Draw();

  if( m_save )
    c->SaveAs((m_savedir+"DistanceToHuman.png"));

}



//---------------------------------------------//
// Plot E vs distance
//---------------------------------------------//
void plotEvsDistance(TFile* file)
{

  TCanvas* c = makeCanvas("c");
  c->SetLogz();
  c->SetRightMargin(0.17);

  TString pname = "p_EvsCross";

  TString xtitle = "Energy [MeV]";
  TString ytitle = "Distance to human [MeV]";
  TString ztitle = "Photons / Primary";
  
  TProfile2D* prof = (TProfile2D*) file->Get(pname);
  prof->SetStats(0);
  prof->SetTitle("");
  prof->GetZaxis()->SetTitleOffset(1.5);

  prof->SetMinimum(1e-6);
  prof->Draw("colz");
  

  if( m_save )
    c->SaveAs((m_savedir+"DistanceToHumanVsEnergy.png"));

}

//---------------------------------------------//
// Plot E for different distance
//---------------------------------------------//
void plotEVariousDistance(TFile* file)
{

  // Make canvas
  TCanvas* c = makeCanvas("c");
  c->SetLogy();

  // Plot name
  TString pname = "p_EvsCross";

  // Set some plot info
  TString xtitle = "Energy [MeV]";
  TString ytitle = "Photons / Primary";
  int colors[] = {kBlack, kBlue, kRed, kGreen, kMagenta};
  int markers[] = {20,25,21,24,23};

  // I only want to plot 4, 8, 12, 16, 20 m
  // So pick out those bins.
  const int pbins = 5;
  int bins[pbins] = {3,7,11,15,18};
  vector<TString> lnames;
  lnames.push_back("4m");
  lnames.push_back("8m");
  lnames.push_back("12m");
  lnames.push_back("16m");
  lnames.push_back("20m");  

  // Make legend
  TLegend* leg = makeLegend(0.6,0.7,0.7,0.93);

  // Load the profile
  TProfile2D* prof = (TProfile2D*) file->Get(pname);

  // Now we will make 1D projections from this
  TH1F* hists[pbins];
  float maximum = -999;
  stringstream ss;
  for(int i=0; i<pbins; ++i){
    ss.str("");
    ss << "proj" << i;
    hists[i] = (TH1F*)  prof->ProjectionX(ss.str().c_str(),bins[i],bins[i]);
    setHistAtt(hists[i],xtitle,ytitle,colors[i],markers[i]);
    leg->AddEntry(hists[i], lnames.at(i), "lep");
    
    if( maximum < hists[i]->GetMaximum() )
      maximum = hists[i]->GetMaximum();
    
  }

  // Now plot
  hists[0]->SetMaximum(5*maximum);
  hists[0]->SetMinimum(1e-4*maximum);
  hists[0]->Draw();
  for(int i=1; i<pbins; ++i)
    hists[i]->Draw("same");
  leg->Draw("same");

  if(m_save )
    c->SaveAs((m_savedir+"EnergyVariousDistances.png"));

}


//---------------------------------------------//
// Plot E for different distance
//---------------------------------------------//
void plotTotEvsDistance(TFile* file)
{

  // Make canvas
  TCanvas* c = makeCanvas("c");
  //c->SetLogy();

  // Plot name
  TString pname = "p_EvsCross";

  // Set some plot info
  TString xtitle = "Distance [m]";
  TString ytitle = "Total Energy [MeV] / Primary";
  int colors[] = {kBlack, kBlue, kRed, kGreen, kMagenta};
  int markers[] = {20,25,21,24,23};

  // Make legend
  TLegend* leg = makeLegend(0.6,0.7,0.7,0.93);

  // Load the profile
  TProfile2D* prof = (TProfile2D*) file->Get(pname);
  TH1F* h_dVsE = (TH1F*) prof->ProjectionY("proj");
  h_dVsE->Reset();
  setHistAtt(h_dVsE,xtitle,ytitle,kBlack,20);

  // Now loop over the profile and sum up the 
  // energy per distance.
  int nybins = prof->GetNbinsY();
  int nxbins = prof->GetNbinsX();
  for(int ybin = 1; ybin<=nybins; ++ybin){

    float totalE = 0;
    for(int xbin = 1; xbin<=nxbins; ++xbin){
      
      float bcen = prof->GetXaxis()->GetBinCenter(xbin);
      float bcon = prof->GetBinContent(xbin,ybin);
      
      totalE += bcen * bcon;

    } // end loop over energy bins

    h_dVsE->SetBinContent(ybin,totalE);

  }// end loop over ybins
    
  h_dVsE->Draw();
      
  if(m_save )
    c->SaveAs((m_savedir+"TotalEnergyVsDistance.png"));

}


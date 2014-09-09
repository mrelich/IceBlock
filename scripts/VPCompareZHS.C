
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

TString m_savedir = "plots/EField/ZHSVPComp/";

//----------------------------------------//
// Main
//----------------------------------------//
void VPCompareZHS(int opt = -1, bool save = false)
{

  m_save = save;

  // Specify ZHS and G4 files
  TString zhsdir = "../ZHS_ELSEnergy/rootfiles/";
  vector<TString> f_zhs;
  TString g4dir = "efieldroot/";
  vector<TString> f_g4;
  vector<TString> savenames;

  // 100 GeV Shower
  if(opt == 0){
    f_zhs.push_back(zhsdir+"beam100e3MeV_1Prim_50NEvt_Angular.root");
    f_g4.push_back(g4dir+"Output_50Evt_100GeV_1Prim_HardCodedAntenna_R1000m_newV.root");
    //f_g4.push_back(g4dir+"Output_50Evt_100GeV_1Prim_HardCodedAntenna_R1000m_oldV.root");
    savenames.push_back("VP_50Evt_100GeV_1Prim");
  }
  else if(opt == 1){
    // 1 TeV Shower
    f_zhs.push_back(zhsdir+"beam1e6MeV_1Prim_20NEvt_Angular.root");
    //f_g4.push_back(g4dir+"Output_20Evt_1TeV_1Prim_HardCoded.root");
    //f_g4.push_back(g4dir+"Output_20Evt_1TeV_1Prim_HardCoded_R1000.root");
    //f_g4.push_back(g4dir+"TestAug19_20Evt_1TeV.root");
    //f_g4.push_back(g4dir+"Output_19Evt_1TeV_1Prim_HardCoded_R1000m.root");
    f_g4.push_back(g4dir+"Output_20Evt_1TeV_1Prim_HardCodedAntenna_R1000m_newV.root");
    savenames.push_back("VP_20Evt_1TeV_1Prim");
  }
  else if(opt ==2 ){
    // 40 MeV Shower -- 1 primary
    //f_zhs.push_back(zhsdir+"beam40MeV_100Prim_100NEvt_Angular.root");
    //f_g4.push_back(g4dir+"Output_100Evt_40MeV_100Prim_HardCoded.root");
    //savenames.push_back("VP_100Evt_40GeV_1Prim");
    //f_g4.push_back(g4dir+"Output_20Evt_40MeV_100Prim_HardCodedAntenna_R100m.root");
    //f_g4.push_back(g4dir+"Output_50Evt_40MeV_100Prim_HardCodedAntenna_R1000m.root");
    //savenames.push_back("VP_100Evt_40GeV_1Prim");
  }
  else if(opt == 3){
    plotPeakComp();
    return;
  }

  // Specify the ZHS histogram name
  TString p_zhsname = "VP_avg_55.829616";
  //TString p_g4name  = "A_AntNum_0_pos_82.7371_0_56.1656";
  //TString p_g4name  = "A_AntNum_1_pos_82.7371_0_56.1656";
  //TString p_g4name  = "A_AntNum_1_pos_827.371_0_561.656";
  TString p_g4name  = "A_AntNum_0_pos_827.371_0_561.656";


  // Now plot
  //plot(f_zhs, f_g4, p_zhsname, p_g4name,savenames);
  plotWithRatio(f_zhs, f_g4, p_zhsname, p_g4name,savenames);
  
}

//----------------------------------------//
// Plot ZHS and G4 on same figure
//----------------------------------------//
void plot(vector<TString> f_zhs,
	  vector<TString> f_g4,
	  TString p_zhs,
	  TString p_g4,
	  vector<TString> savenames)
{

  // Make canvas
  TCanvas* c = makeCanvas("c");
  c->SetLogy();

  // Make a legend
  TLegend* leg = makeLegend(0.15,0.3,0.8,0.92);  

  // Titles
  TString xtitle = "time [ns]";
  TString ytitle = "A [Vs/m]";    

  // Loop and get plots
  TProfile* prof_Z = NULL;
  TProfile* prof_G = NULL;
  TH1D* h_resestZ  = NULL;
  for(unsigned int i=0; i<f_zhs.size(); ++i){

    // Load ZHS profile
    TFile* file_Z = new TFile(f_zhs.at(i).Data());
    prof_Z = getProfile(file_Z, p_zhs, xtitle, ytitle, kBlue, 20);
    prof_Z->SetDirectory(0);
    file_Z->Close();
    leg->AddEntry(prof_Z,"ZHS","lep");

    // Load Geant profile
    TFile* file_G = new TFile(f_g4.at(i).Data());
    prof_G = getProfile(file_G, p_g4, xtitle, ytitle, kRed, 20);
    prof_G->SetDirectory(0);
    file_G->Close();
    leg->AddEntry(prof_G,"Geant4","lep");
    
    // Scale Geant4 profile up by factor of R
    //prof_G->Scale(m_R);
    
    // reset the zhs timing info
    //prof_Z = resetZHS(prof_Z, prof_G, kBlue, m_R);
    h_resetZ = resetZHS(prof_Z, prof_G, kBlue, m_R);

    // Get maximum
    float maximum = prof_G->GetMaximum();
    if( maximum < h_resetZ->GetMaximum() )
      maximum = h_resetZ->GetMaximum();
    prof_G->SetMaximum(maximum*5);
    prof_G->SetMinimum(maximum*1e-5);

    // Now plot
    prof_G->Draw();
    //prof_Z->Draw("sameep");
    h_resetZ->Draw("sameep");
    leg->Draw("same");
    //prof_G->Draw("same");

    if(m_save){
      //TString save = m_savedir + p_g4 + "_ZHSCherAngle.png";
      TString save = m_savedir + savenames.at(i) + ".png";
      c->SaveAs(save.Data());
    }

  }// end loop over files


}

//----------------------------------------//
// Plot ZHS and G4 on same figure
//----------------------------------------//
void plotWithRatio(vector<TString> f_zhs,
		   vector<TString> f_g4,
		   TString p_zhs,
		   TString p_g4,
		   vector<TString> savenames)
{

  // Make canvas
  TCanvas* c = makeCanvas("c");
  //c->SetLogy();

  // Make a legend
  TLegend* leg = makeLegend(0.15,0.3,0.79,0.92);  
  leg->SetTextSize(0.055);

  // Titles
  TString xtitle = "time [ns]";
  TString ytitle = "A [Vs/m]";    

  // Loop and get plots
  TProfile* prof_Z = NULL;
  TProfile* prof_G = NULL;
  TH1D* h_resetZ   = NULL;
  for(unsigned int i=0; i<f_zhs.size(); ++i){

    // Load ZHS profile
    TFile* file_Z = new TFile(f_zhs.at(i).Data());
    prof_Z = getProfile(file_Z, p_zhs, xtitle, ytitle, kBlue, 20);
    prof_Z->SetDirectory(0);
    file_Z->Close();
    leg->AddEntry(prof_Z,"ZHS","lep");

    // Load Geant profile
    TFile* file_G = new TFile(f_g4.at(i).Data());
    prof_G = getProfile(file_G, p_g4, xtitle, ytitle, kRed, 20);
    prof_G->SetDirectory(0);
    file_G->Close();
    leg->AddEntry(prof_G,"Geant4","lep");
    
    // Scale Geant4 profile up by factor of R
    //prof_G->Scale(m_R);
    
    // reset the zhs timing info
    //prof_Z = resetZHS(prof_Z, prof_G, kBlue, m_R);
    h_resetZ = resetZHS(prof_Z, prof_G, kBlue, m_R);
    //prof_Z->Scale(1/m_R);

    // Get maximum
    float maximum = prof_G->GetMaximum();
    if( maximum < h_resetZ->GetMaximum() )
      maximum = h_resetZ->GetMaximum();
    prof_G->SetMaximum(maximum*5);
    prof_G->SetMinimum(maximum*1e-5);

    // Make two pads
    TPad* top = NULL;
    TPad* bot = NULL;
    makePads(c, top, bot);

    // Set some attributes
    prof_G->GetYaxis()->SetLabelSize(0.05);
    prof_G->GetYaxis()->SetTitleSize(0.05);
    prof_G->GetYaxis()->SetTitleOffset(1.0);
    prof_G->GetXaxis()->SetLabelSize(0);

    // Draw top plots
    c->cd();
    top->Draw();
    top->cd();
    top->SetLogy();
    prof_G->Draw();
    //prof_Z->Draw("sameep");
    h_resetZ->Draw("sameep");
    leg->Draw("same");
    top->Update();

    // Draw bot ratio
    c->cd();
    bot->Draw();
    bot->cd();
    //TProfile* prof = prof_G->Clone("ratio");
    TH1D* prof = prof_G->ProjectionX("ratio");
    prof->SetLineColor(prof_G->GetLineColor());
    prof->SetMarkerColor(prof_G->GetMarkerColor());
    prof->SetStats(0);
    //prof->Divide(prof_Z);
    prof->Divide(h_resetZ);
    prof->GetYaxis()->SetTitle("G4/ZHS");
    prof->SetMinimum(0);
    //prof->SetMaximum(50);
    prof->SetMaximum(2);
    prof->GetXaxis()->SetLabelSize(0.1);
    prof->GetYaxis()->SetLabelSize(0.1);
    prof->GetXaxis()->SetTitleSize(0.12);
    prof->GetYaxis()->SetTitleSize(0.12);
    prof->GetYaxis()->SetTitleOffset(0.4);
    prof->GetYaxis()->SetNdivisions(405);
    prof->Draw();
    bot->Update();
    
    if(m_save){
      //TString save = m_savedir + p_g4 + "_ZHSCherAngle.png";
      TString save = m_savedir + savenames.at(i) + "_ratio.png";
      c->SaveAs(save.Data());
    }

  }// end loop over files


}

//----------------------------------------//
// Reset ZHS to proper time
//----------------------------------------//
TH1D* resetZHS(TProfile* zhs, TProfile* g4,
	       int color, double scale)
{

  //TProfile* zhs_reset = g4->Clone("zhs_reset");
  TH1D* zhs_reset = g4->ProjectionX("zhs_reset");
  zhs_reset->Reset();
  zhs_reset->SetLineColor(color);
  zhs_reset->SetMarkerColor(color);

  // Find bin containing maximum value of geant4 profile
  double maximum = g4->GetXaxis()->GetBinCenter(g4->GetMaximumBin());

  // Offset ZHS by this maximum bin
  int nbins = zhs->GetNbinsX(); 
  for(int bin = 1; bin<=nbins; ++bin){
    double bc = zhs->GetBinContent(bin);
    double be = zhs->GetBinError(bin);
    double time = zhs->GetXaxis()->GetBinCenter(bin);
    int newbin = zhs_reset->FindBin(time+maximum);
    double prev = zhs_reset->GetBinContent(bin);
    zhs_reset->SetBinContent(newbin,bc/scale);
    zhs_reset->SetBinError(newbin,be/scale);
    //zhs_reset->SetBinEntries(newbin,1);

    //cout<<"Filling: "<<time<<" "<<time+maximum<<" "<<bc<<endl;
  }
  
  //cout<<"Maximum : "<<maximum<<endl;

  return zhs_reset;

}  

//----------------------------------------//
// Plot the peak of a few values
//----------------------------------------//
void plotPeakComp()
{

  // Make Canvas
  TCanvas* c = makeCanvas("c");
  c->SetLogy();
  //c->SetLogx();

  // G4 and ZHS plots
  TString g4pname = "A_AntNum_0_pos_827.371_0_561.656";
  TString zhspname = "VP_avg_55.829616";

  // Specify the G4 plots
  TString g4dir = "efieldroot/";
  vector<TString> g4files;
  g4files.push_back(g4dir+"Output_50Evt_1GeV_1Prim_HardCodedAntenna_R1000m_newV.root");
  g4files.push_back(g4dir+"Output_50Evt_10GeV_1Prim_HardCodedAntenna_R1000m_newV.root");
  g4files.push_back(g4dir+"Output_50Evt_100GeV_1Prim_HardCodedAntenna_R1000m_newV.root");
  g4files.push_back(g4dir+"Output_20Evt_1TeV_1Prim_HardCodedAntenna_R1000m_newV.root");
  //g4files.push_back(g4dir+"Output_20Evt_10TeV_1Prim_HardCodedAntenna_R1000m_newV.root");

  // Specify the G4 plots
  TString zhsdir = "../ZHS_ELSEnergy/rootfiles/";
  vector<TString> zhsfiles;
  zhsfiles.push_back(zhsdir+"beam1e3MeV_1Prim_50NEvt_Angular.root");
  zhsfiles.push_back(zhsdir+"beam10e3MeV_1Prim_50NEvt_Angular.root");
  zhsfiles.push_back(zhsdir+"beam100e3MeV_1Prim_50NEvt_Angular.root");
  zhsfiles.push_back(zhsdir+"beam1e6MeV_1Prim_20NEvt_Angular.root");


  // What energies correspond to file [GeV]
  vector<double> NRG;
  NRG.push_back(1);
  NRG.push_back(10);
  NRG.push_back(100);
  NRG.push_back(1000);
  NRG.push_back(10000);

  // Specify some histogram shiz
  TString xtitle = "log(Energy/GeV)";
  TString ytitle = "Max(A) [Vm/s]";
  TH1F* frame = makeHist("frame",1,-1,5,xtitle,ytitle,kBlack,20);

  // Now loop over files and get points
  int npoints = 0;
  double E[1000];
  double g4_max[1000];
  double zhs_max[1000];
  double maximum = -999;
  for(unsigned int i=0; i<g4files.size(); ++i){

    // Store energy 
    E[i] = log10(NRG.at(i));

    // Get G4 result
    TFile* f_g4 = new TFile(g4files.at(i).Data());
    TProfile* prof = getProfile(f_g4, g4pname,"","",kBlack,20);
    //g4_max[i] = prof->GetMaximum();
    g4_max[i] = getMax(prof);
    if( maximum < g4_max[i] ) maximum = g4_max[i];

    // Get G4 result
    TFile* f_zhs = new TFile(zhsfiles.at(i).Data());
    TProfile* prof = getProfile(f_zhs, zhspname,"","",kBlack,20);
    //zhs_max[i] = prof->GetMaximum() / 1000.;
    //zhs_max[i] = prof->GetMaximum() /1000. ;
    zhs_max[i] = getMax(prof) / 1000.;
    if( maximum < zhs_max[i] ) maximum = zhs_max[i];

    // Increment points
    npoints++;

  }

  // Make TGraph
  TGraph* gr_g4 = new TGraph(npoints, E, g4_max);
  gr_g4->SetLineWidth(1);
  gr_g4->SetMarkerSize(1);
  gr_g4->SetMarkerStyle(20);
  gr_g4->SetLineColor(kBlue);
  gr_g4->SetMarkerColor(kBlue);

  // Make TGraph
  TGraph* gr_zhs = new TGraph(npoints, E, zhs_max);
  gr_zhs->SetLineWidth(1);
  gr_zhs->SetMarkerSize(1);
  gr_zhs->SetMarkerStyle(20);
  gr_zhs->SetLineColor(kRed);
  gr_zhs->SetMarkerColor(kRed);

  // Draw graph
  frame->SetMaximum(5*maximum);
  frame->SetMinimum(1e-5*maximum);
  frame->Draw();
  gr_g4->Draw("samelp");
  gr_zhs->Draw("samelp");

  // Add some legend
  TLegend* leg = makeLegend(0.2,0.4,0.8,0.9);
  leg->AddEntry(gr_g4,"Geant4","lp");  
  leg->AddEntry(gr_zhs,"ZHS","lp");
  leg->Draw("same");

}

double getMax(TProfile* p)
{

  int nbins = p->GetNbinsX();
  double maxi = 0;
  for(int bin=1; bin<=nbins; ++bin){
    double bc = p->GetBinContent(bin);
    if( bc > maxi ) maxi = bc;
  }

  return maxi;
}

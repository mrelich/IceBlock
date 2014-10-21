
#include "myHist.C"
#include "EField.C"

// Where to save
TString m_savedir = "plots/beamProfile/";

//bool m_save = false;
bool m_save = true;

void VPCompBeamProfile()
{
  
  // Update -- loop over antennas from file
  // Plot name
  //TString pname = "A_AntNum_0_pos_827.371_0_561.656";
  //TString pname = "A_AntNum_0_pos_5.7916_0_3.93159";
  //TString pname = "A_AntNum_1_pos_6_0_9.98568";
  ifstream f_in ("antennaConfig/xzRefracted_40MeV.txt");
  
  // Fnames
  TString indir = "efieldroot/";
  vector<TString> fnames;
  //fnames.push_back(indir+"Output_50Evt_40meV_10000Prim_HCAnt_R1000m_singlePos.root");
  //fnames.push_back(indir+"Output_50Evt_40meV_10000Prim_HCAnt_R1000m_RandFlat3.5.root");
  //fnames.push_back(indir+"Output_50Evt_40meV_10000Prim_HCAnt_R1000m_RandGauss3.5.root");
  //fnames.push_back(indir+"Output_50Evt_40MeV_10000Prim_HCAnt_R7m_singlePos.root");
  //fnames.push_back(indir+"Output_50Evt_40MeV_10000Prim_HCAnt_R7m_RandFlat3.5.root");
  //fnames.push_back(indir+"Output_50Evt_40MeV_10000Prim_HCAnt_R7m_RandGauss3.5.root");
  fnames.push_back(indir+"Output_50Evt_40MeV_10000Prim_xzRefracted_40MeV_singlePos.root");
  fnames.push_back(indir+"Output_50Evt_40MeV_10000Prim_xzRefracted_40MeV_RandFlat3.5.root");
  fnames.push_back(indir+"Output_50Evt_40MeV_10000Prim_xzRefracted_40MeV_RandGauss3.5.root");
  
  // Leg names
  vector<TString> legnames;
  legnames.push_back("Single Position");
  legnames.push_back("Flat R = 3.5mm");
  legnames.push_back("Gaussian #sigma = 3.5mm");

  // Legend
  TLegend* leg = makeLegend(0.6,0.8,0.65,0.93);
  leg->SetTextSize(0.05);

  // Stuff for plots
  TString xtit = "time [ns]";
  TString ytit = "|A| [Vs/m]";
  int colors[] = {kBlack, kRed, kBlue};

  // Holders
  TString antPos   = "";
  double angle    = 0;
  double refAngle = 0;
  double x=0,y=0,z=0;
  double R=0, Rprime = 0;
  double zprime = 0;
  double prevZ  = 0;
  int npoints = 0;
  int counter = 0;

  // Make canvas
  TCanvas* c = makeCanvas("c");
  //c->SetLogy();
  TPad* top = NULL; TPad* bot = NULL;
  makePads(c,top,bot);
  
  // Loop over evennts in file
  while( !f_in.eof() ){

    getAntPos(f_in, counter, antPos, angle, refAngle, 
              x,y,z,
              R, Rprime,
              zprime);
    
    counter ++;
    cout<<"Working on: "<<antPos<<endl;

    if( z == prevZ ) continue;
    prevZ = z;

    // Set legend header
    leg->Clear();
    leg->SetHeader(Form("Angle: %.0f#circ",angle));
    
    
    // Loop and load profs
    float maximum = -9999;
    TProfile* profs[3];
    for(unsigned int i=0; i<fnames.size(); ++i){

      TFile* file = new TFile(fnames.at(i).Data());
      profs[i] = getProfile(file,antPos,xtit,ytit,colors[i],20);
      profs[i]->SetDirectory(0);
      file->Close();
      
      leg->AddEntry(profs[i],legnames.at(i).Data(),"l");
      
      if( maximum < profs[i]->GetMaximum() )
	maximum = profs[i]->GetMaximum();
    }

    // Set min max
    profs[0]->SetMaximum(5*maximum);
    profs[0]->SetMinimum(1e-4*maximum);
    
    //
    profs[0]->GetYaxis()->SetLabelSize(0.05);
    profs[0]->GetYaxis()->SetTitleSize(0.05);
    profs[0]->GetYaxis()->SetTitleOffset(1.0);
    profs[0]->GetXaxis()->SetLabelSize(0);
    
    
    // Set a range
    int xmin = profs[0]->GetXaxis()->FindBin(40);
    int xmax = profs[0]->GetXaxis()->FindBin(45);
    //profs[0]->GetXaxis()->SetRange(xmin,xmax);   
    
    // Now draw top stuff
    c->cd();
    top->Draw();
    top->cd();
    top->SetLogy();
    profs[0]->Draw();
    profs[1]->Draw("same");
    profs[2]->Draw("same");
    leg->Draw("same");
    top->Update();
    
    
    // Get bottom objects
    TH1D* ratio1 = profs[0]->ProjectionX("ratio1");
    TH1D* ratio2 = profs[0]->ProjectionX("ratio2");
    
    ratio1->Divide(profs[1]);
    ratio2->Divide(profs[2]);
    
    ratio1->SetLineColor(profs[1]->GetLineColor());
    ratio2->SetLineColor(profs[2]->GetLineColor());
    ratio1->SetMarkerColor(profs[1]->GetLineColor());
    ratio2->SetMarkerColor(profs[2]->GetLineColor());
    
    ratio1->GetYaxis()->SetTitle("Nom/Dist");
    ratio1->GetXaxis()->SetLabelSize(0.1);
    ratio1->GetYaxis()->SetLabelSize(0.1);
    ratio1->GetXaxis()->SetTitleSize(0.12);
    ratio1->GetYaxis()->SetTitleSize(0.12);
    ratio1->GetYaxis()->SetTitleOffset(0.4);
    ratio1->GetYaxis()->SetNdivisions(405);
    ratio1->SetStats(0);
    
    ratio1->SetMinimum(0);
    ratio1->SetMaximum(2);
    //ratio1->GetXaxis()->SetRange(xmin,xmax);
    
    // Now Draw
    c->cd();
    bot->Draw();
    bot->cd();
    ratio1->Draw();
    ratio2->Draw("same");
    
    // Save
    if( m_save ){
      TString savename = m_savedir + "ComparisonPointFlatGaus_"+antPos+".png";
      //c->SaveAs((m_savedir+"ComparisonPointFlatGaussian_R7m.png").Data());
      c->SaveAs(savename.Data());
      delete ratio1;
      delete ratio2;
      for(unsigned int i=0; i<fnames.size(); ++i)
	delete profs[i];
    }
    
  }// end while loop over antenna file
}

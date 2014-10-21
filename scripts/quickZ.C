
#include "myHist.C"

TString m_savedir = "plots/EField/Angular/";
//bool m_save = false;
bool m_save = true;

void quickZ()
{

  float xmin = 30;
  float xmax = 90;

  // Make Canvas
  TCanvas* c = makeCanvas("c");
  
  // Legend
  //TLegend* leg = makeLegend(0.7,0.85,0.63,0.93);
  TLegend* leg = makeLegend(0.7,0.85,0.73,0.93);
  leg->SetHeader("x = 6m");

  // Make frame
  TH1F* frame = makeHist("frame", 1, xmin, xmax, 
			 "Angle Relative to Beam [deg]",
			 "Antenna Height [m]", 0,0);

  //frame->SetMinimum(-2);
  frame->SetMinimum(-1);
  frame->SetMaximum(20);
  frame->Draw();

  // Unfrefracted angle
  TF1* unref = new TF1("unref","6./tan(x*TMath::Pi()/180)",xmin,xmax);
  unref->SetLineColor(kBlack);
  unref->Draw("same");
  leg->AddEntry(unref,"Nominal","l");

  // Now handle refracted angles:
  vector<float> angles;
  //angles.push_back(20);
  //angles.push_back(30);
  //angles.push_back(40);
  //angles.push_back(50);
  //angles.push_back(60);
  angles.push_back(45);

  int colors[] = {kRed,kBlue,kViolet+2,kGreen,kMagenta};

  TF1* func = NULL;
  for(unsigned int i=0; i<angles.size(); ++i){
    float angle = angles.at(i);
    func = new TF1(TString(Form("ref_%i",i)),refracted,xmin,xmax,1);
    func->SetParameter(0,angle);
    func->SetLineColor(colors[i]);
    func->Draw("same");
    leg->AddEntry(func,Form("#alpha = %.0f#circ",angle),"l");
  }

  // Refracted w/ rotation 30 deg
  //TF1* ref30 = new TF1("ref30",refracted,xmin,xmax,1);
  //ref30->SetParameter(0,30); // 30 deg
  //ref30->SetLineColor(kRed);
  //ref30->Draw("same");
  //leg->AddEntry(ref30,"Tilt 30#circ", "l");
  //leg->AddEntry(ref30,"#alpha = 60#circ", "l");

  // Refracted w/ rotation 60 deg
  //TF1* ref60 = new TF1("ref60",refracted,xmin,xmax,1);
  //ref60->SetParameter(0,60); // 60 deg
  //ref60->SetLineColor(kViolet);
  //ref60->Draw("same");
  //leg->AddEntry(ref60,"Tilt 60#circ", "l");
  //leg->AddEntry(ref60,"#alpha = 30#circ", "l");
  

  // Draw legend
  leg->Draw("same");

  // Get a line for the height of the wall
  //float height = 10 - 3.6576;
  float h_wall = 5;
  vector<float> heights;
  heights.push_back(17-h_wall);
  heights.push_back(14-h_wall);
  heights.push_back(10-h_wall);

  int wcolors[] = {kBlack, kBlue, kRed};

  TLegend* wleg = makeLegend(0.4,0.55,0.75,0.92);
  wleg->SetHeader("Tower Height");
  for(unsigned int i=0; i<heights.size(); ++i){
    float height = heights.at(i);
    TLine* l_ant = makeLine(xmin,xmax,height,height,wcolors[i],2);
    l_ant->SetLineWidth(2);
    l_ant->Draw("same");
    wleg->AddEntry(l_ant,Form("%.0fm",height+h_wall),"l");
  }
  wleg->Draw("same");

  // Save Result
  if( m_save )
    //c->SaveAs((m_savedir+"IceTilt_AntHeightVsRelativeAngle.png"));
    c->SaveAs((m_savedir+"IceTilt_AntHeightVsRelativeAngle_45Only.png"));

}

Double_t refracted(Double_t *x, Double_t *par)
{

  // Define pi for convenience
  double pi = TMath::Pi();

  // Angle is given in degrees
  double angle = x[0] * pi / 180.;
  //cout<<"**********************"<<endl;
  //cout<<"Angle: "<<angle*180/pi<<endl;

  // Load the rotation
  double rot = par[0] * pi / 180.;
  ////cout<<"Rot: "<<rot*180/pi<<endl;

  // Determine the incident angle
  double inc = fabs(angle - rot);
  ////cout<<"Inc: "<<inc*180/pi<<endl;  

  // Make sure we don't go past criticle angle
  if( inc > asin(1/1.78) ) return -9999;

  // Get refracted angle
  double ref = asin(1.78*sin(inc));
  //cout<<"Ref: "<<ref*180/pi<<endl;

  // Take refracted angle back to 
  // the ref frame.  Be careful to
  // handle the boundaries right
  double ref_prime = pi/2. - ref - rot;
  if( angle < rot ) ref_prime = pi/2. + ref - rot;
  //cout<<"Angle in frame: "<<ref_prime*180/pi<<endl;


  // Now return the new z
  return 6 * tan(ref_prime);
  

}

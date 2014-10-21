
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-//
// This is a quick script trying to estimate the gammas that //
// are escaping the ice and staying in the ice.              //
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-//

#include "myHist.C"

void RadiationPlots()
{

  TCanvas* c = makeCanvas("c");
  c->SetLogy();

  TFile* file = new TFile("rootfiles/tracks/test.root");

  vector<TString> pnames;
  pnames.push_back("p_Ein");
  pnames.push_back("p_Eleave");
  pnames.push_back("p_Etop");
  pnames.push_back("p_Eside");
  
  vector<TString> legnames;
  legnames.push_back("Total In Ice");
  legnames.push_back("Total Leaving");
  legnames.push_back("Leaving Top");
  legnames.push_back("Leaving Side");

  TLegend* leg = makeLegend(0.6,0.75,0.7,0.93);

  int colors[] = {kBlack, kBlue, kGreen, kRed};

  TString xtitle = "Energy [MeV]";
  TString ytitle = "Photons / Primary";

  for(unsigned int i=0; i<pnames.size(); ++i){
  
    TH1F* h = getHist(file,pnames.at(i),xtitle,ytitle,colors[i],20);
    if(i==0) h->Draw();
    else h->Draw("same");

    leg->AddEntry(h,legnames.at(i),"l");
    
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

  c->SaveAs("RadiationEstimate.png");


}

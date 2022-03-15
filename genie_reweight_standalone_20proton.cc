#include <iostream>

#include <math.h>

#include "BDTReweighter.h"
#include "MissingProtonFakeData_BDTRW_FHC.h"
#include "MissingProtonFakeData_BDTRW_RHC.h"



#include "TChain.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TRatioPlot.h"
#include "TFile.h"
#include "TString.h"


using namespace std;

int mode_func(int j){
  if (j==0) return 1;
  else return 0;
}

int main(int argc, char ** argv){

  // Location of CAF files
  std::string pathToCAFfile = "~/Physics/PhD/root_analysis_files/CAFv7_0mgsimple_subset131.root";
  
  //std::string pathToCAFfile = argv[1];

  TChain * caf = new TChain("cafTree");
  caf->Add(pathToCAFfile.c_str());

  TString pic = argv[2];
  cout << " check 1 \n";

  // Variables we need from the CAF files.
  int nuPDG, isCC, nP, nN, nipip, nipim, nipi0, niem, isRES, isMEC, isDIS, isCOH, isQES, mode;

  double Ev, LepE, LepNuAngle, Q2, W, X, Y, eP, eN, ePip, ePim, ePi0;
  double eRecoP , eRecoN , eRecoPip , eRecoPim , eRecoPi0 , eRecoOther;
  double eRec_FromDep, eRec_FromDep20;

  caf->SetBranchAddress("nuPDG", &nuPDG);
  caf->SetBranchAddress("isCC", &isCC);
  caf->SetBranchAddress("mode", &mode);

  /*
  caf->SetBranchAddress("isDIS", &isDIS);
  caf->SetBranchAddress("isQES", &isQES);
  caf->SetBranchAddress("isMEC", &isMEC);
  caf->SetBranchAddress("isRES", &isRES);
  caf->SetBranchAddress("isCOH", &isCOH);
  */

  caf->SetBranchAddress("nP", &nP);
  caf->SetBranchAddress("nN", &nN);
  caf->SetBranchAddress("nipip", &nipip);
  caf->SetBranchAddress("nipim", &nipim);
  caf->SetBranchAddress("nipi0", &nipi0);
  caf->SetBranchAddress("niem", &niem);
  caf->SetBranchAddress("eP", &eP);
  caf->SetBranchAddress("eN", &eN);
  caf->SetBranchAddress("ePip", &ePip);
  caf->SetBranchAddress("ePim", &ePim);
  caf->SetBranchAddress("ePi0", &ePi0);
  caf->SetBranchAddress("Ev", &Ev);
  caf->SetBranchAddress("LepE", &LepE);
  caf->SetBranchAddress("LepNuAngle", &LepNuAngle);
  caf->SetBranchAddress("Q2", &Q2);
  caf->SetBranchAddress("W", &W);
  caf->SetBranchAddress("X", &X);
  caf->SetBranchAddress("Y", &Y);

  caf->SetBranchAddress("eRecoP", &eRecoP);
  caf->SetBranchAddress("eRecoN", &eRecoN);
  caf->SetBranchAddress("eRecoPip", &eRecoPip);
  caf->SetBranchAddress("eRecoPim", &eRecoPim);
  caf->SetBranchAddress("eRecoPi0", &eRecoPi0);
  caf->SetBranchAddress("eRecoOther", &eRecoOther);
  caf->SetBranchAddress("LepE", &LepE);


  cout << " check 2\n";

  // Get BDT reweighter and set ovarall normalization factors.
  std::vector<BDTReweighter*> bdt_reweighter;

  bdt_reweighter.push_back(new MissingProtonFakeData_BDTRW_FHC());
  bdt_reweighter.push_back(new MissingProtonFakeData_BDTRW_RHC());


  // Set up array of features. This is the input to the BDT reweighter
  union BDTReweighter::BDTReweighterFeature features[8];

  // Loop through events in CAF file
  unsigned int n_entries = caf->GetEntries();

   // Set up a histogram
  TH1D * heRec = new TH1D("heRec", "eRec", 100, 0, 10);
  TH1D * heRecReweighted = new TH1D("heRecReweighted", "eRec", 100, 0, 10);

  for (unsigned int i_entry = 0; i_entry < n_entries; i_entry++){
    caf->GetEntry(i_entry);

    //std::cout << " ----> event " << i_entry << "\n";


    // Fill features array
    features[0].fvalue = mode_func(mode-1);   // QES
    features[1].fvalue = mode_func(mode-3);   // DIS
    features[2].fvalue = mode_func(mode-4);   // RES
    features[3].fvalue = mode_func(mode-5);   // COH
    features[4].fvalue = mode_func(mode-10);  // MEC

    features[5].fvalue = Ev;  // compare the values
    features[6].fvalue = eP;
    features[7].fvalue = 1.0 - LepE/Ev; 


    // Evaluate the BDT and convert to weight.
    float wght_val = bdt_reweighter[nuPDG > 0 ? 0 : 1]->GetWeight(features, 1);
    eRec_FromDep = eRecoP + eRecoN + eRecoPip + eRecoPim + eRecoPi0 + eRecoOther + LepE;
    eRec_FromDep20 = eRecoP + eRecoN + eRecoPip + eRecoPim + eRecoPi0 + eRecoOther + LepE;

    heRec -> Fill(eRec_FromDep); heRecReweighted->Fill(eRec_FromDep20,wght_val);

    //std::cout << "Event " << i_entry << " weight " << wght_val << std::endl;
  }
  
  TCanvas * c1 = new TCanvas();
  double factor = heRecReweighted->Integral()*1.0/heRec->GetEntries();
  cout << " factor -->  " << factor << "\n";
  //heRec->Scale(factor); 
  //heRec->Scale(1.1); 


  heRec->SetLineColor(kRed); heRecReweighted->SetLineColor(kGreen);

  heRecReweighted->Draw("HIST");
  heRec->Draw("SAMEHIST");

  heRec->SetXTitle("recoNuErg (GeV)");

  auto legend = new TLegend(0.5,0.5,0.7,0.75);
  legend->SetHeader("recoNuErg","C"); // option "C" allows to center the header
  legend->AddEntry(heRecReweighted,"reweighted events","f");
  legend->AddEntry(heRec,"events","f");
  legend->Draw("same");

  c1->SaveAs("Test_reweight_eP20.png");


};

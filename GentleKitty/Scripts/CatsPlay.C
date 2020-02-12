#include "PlayWithCats.h"
#include "TidyCats.cxx"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TApplication.h"

int main(int argc, char *argv[]) {
  TApplication theApp("App",&argc, argv);
  PlayWithCats *catsPlay = new PlayWithCats();
  catsPlay->GenerateSourceDistpxi();
  theApp.Run();
//  catsPlay->CloseFile();

  //const char* Data = (argv[1]) ? argv[1] : "";
  //const char* Fit = (argv[2]) ? argv[2] : "";
  //if (Fit != "")catsPlay->ExtractUncertaintyFit(Fit);
  //if (Data != "")catsPlay->ExtractUncertaintyData(Data);
  //catsPlay->GenerateDefault();
//  catsPlay->ShiftBinning();
  //catsPlay->GenerateCoulombOnly();
  //catsPlay->PlotPotentials();
  //catsPlay->PlotPotentialSum();

//   TidyCats* tidy = new TidyCats();
//
//  TFile* out = TFile::Open("out.root","recreate");
//  for (int irad = 0 ; irad < 7; ++irad) {
//
//    CATS AB_ppup;
//    CATS AB_pplow;
//
//    tidy->GetCatsProtonProton(&AB_ppup, 125, 0, 250, TidyCats::sGaussian);
//    tidy->GetCatsProtonProton(&AB_pplow, 125, 0, 250, TidyCats::sGaussian);
//
//
//    TGraphErrors* graph = new TGraphErrors();
//    AB_ppup.SetAnaSource(0,1.48-0.1*irad);
//    AB_ppup.KillTheCat();
//    AB_pplow.SetAnaSource(0,1.52-0.1*irad);
//    AB_pplow.KillTheCat();
//
//    for (int ikStar = 0; ikStar < 125; ++ikStar ) {
//      double mean  = (AB_ppup.GetCorrFun(ikStar)+AB_pplow.GetCorrFun(ikStar))/2.;
//      graph->SetPoint(ikStar,AB_ppup.GetMomentum(ikStar),mean);
//      double Err = TMath::Abs(AB_ppup.GetCorrFun(ikStar)-mean);
//      graph->SetPointError(ikStar,0,Err);
//    }
//    graph->SetName(TString::Format("mTBin_%u",irad).Data());
//    out->cd();
//    graph->Write();
//  }
//  std::cout << "Done looping \n";
//  out->Close();
  return 0;
}

#include "TROOT.h"
#include "TSystem.h"
#include "ReadDreamFile.h"
#include "TCanvas.h"
#include <iostream>
#include <vector>
void AddFakeReweighting(DreamPair* Pair, DreamPair* aPair) {
  TFile* inFake = TFile::Open(
      "~/cernbox/pPb/v0offlineFix/woDetadPhi/CFOutput_pXi.root");
  TList* PairList = (TList*) inFake->Get("PairDist");
  if (!PairList) {
    std::cout << "No PairDist List\n";
    return;
  }
  TList* ThePair = (TList*) PairList->FindObject("Pair");
  if (!ThePair) {
    std::cout << "No Pair List\n";
    return;
  }
  TH2F* SEMultPart =  (TH2F*) ThePair->FindObject("SEMultDist_Particle0_Particle4_clone");
  if (!SEMultPart) {
    std::cout << "SEMultDist_Particle0_Particle4_clone missing\n ";
    return;
  } else {
    Pair->GetPair()->SetSEMultDist(SEMultPart,"");
  }
  TH2F* MEMultPart = (TH2F*) ThePair->FindObject(
      "MEMultDist_Particle0_Particle4_clone");
  if (!MEMultPart) {
    std::cout << "MEMultDist_Particle0_Particle4_clone missing\n ";
    return;
  } else {
    Pair->GetPair()->SetMEMultDist(MEMultPart, "");
  }

  TList* aPairList = (TList*) inFake->Get("AntiPairDist");
  if (!aPairList) {
    std::cout << "No aPairDist List\n";
    return;
  }
  TList* aThePair = (TList*) aPairList->FindObject("Pair");
  if (!aThePair) {
    std::cout << "No Anti Pair List\n";
    return;
  }
  TH2F* SEMultaPart =  (TH2F*) aThePair->FindObject("SEMultDist_Particle1_Particle5_clone");
  if (!SEMultaPart) {
    std::cout << "SEMultDist_Particle1_Particle5_clone missing\n ";
    return;
  } else {
    aPair->GetPair()->SetSEMultDist(SEMultaPart,"");
  }
  TH2F* MEMultaPart = (TH2F*) aThePair->FindObject(
      "MEMultDist_Particle1_Particle5_clone");
  if (!MEMultaPart) {
    std::cout << "MEMultDist_Particle1_Particle5_clone missing\n ";
    return;
  } else {
    aPair->GetPair()->SetMEMultDist(MEMultaPart, "");
  }
}

void ReweightManually(TH1F* CF) {
  TFile* inFake = TFile::Open(
      "~/cernbox/pPb/v0offlineFix/woDetadPhi/CFOutput_pXi.root");
  TH1F* noRew = (TH1F*)inFake->Get("hCk_RebinnedMeV_1");
  TH1F* Rew = (TH1F*)inFake->Get("hCk_ReweightedMeV_1");
  Rew->Divide(noRew);
  CF->Multiply(Rew);
}

void ProcessVariation(ReadDreamFile* DreamFile, TCanvas* c1) {
  static int outNum = 1;
  DreamCF* CF_pXi = new DreamCF();
  DreamPair* pXi = new DreamPair("Part", 0.24, 0.34);
  DreamPair* ApAXi = new DreamPair("AntiPart", 0.24, 0.34);
  pXi->SetPair(DreamFile->GetPairDistributions(0, 4, ""));
  ApAXi->SetPair(DreamFile->GetPairDistributions(1, 5, ""));
  AddFakeReweighting(pXi, ApAXi);

  pXi->FixShift(pXi->GetPair(),nullptr,0.008,true);
  ApAXi->FixShift(ApAXi->GetPair(),nullptr,0.008,true);

  pXi->Rebin(pXi->GetPairFixShifted(0), 5);
  ApAXi->Rebin(ApAXi->GetPairFixShifted(0), 5);

  pXi->Rebin(pXi->GetPairFixShifted(0), 5);
  ApAXi->Rebin(ApAXi->GetPairFixShifted(0), 5);

//  pXi->ReweightMixedEvent(pXi->GetPairRebinned(0), 0.2, 0.9);
//  ApAXi->ReweightMixedEvent(ApAXi->GetPairRebinned(0), 0.2, 0.9);

  CF_pXi->SetPairs(pXi,ApAXi);
  CF_pXi->GetCorrelations();
  c1->cd();
  TH1F* CFDraw = CF_pXi->FindCorrelationFunction("hCk_RebinnedMeV_1");
  ReweightManually(CFDraw);
  CFDraw->SetStats(false);
  CFDraw->GetXaxis()->SetRangeUser(0,800);
  CFDraw->GetYaxis()->SetRangeUser(0.8,3.);
  if (outNum == 1) CFDraw->DrawCopy();
  else CFDraw->DrawCopy("SAME");
  CF_pXi->WriteOutput(TString::Format("%s/CFpXiVariations_%u.root",gSystem->pwd(),outNum++));
}

int main(int argc, char* argv[]) {
  gROOT->ProcessLine("gErrorIgnoreLevel = 2001;");
  const char* filename = argv[1];
  const char* prefix = argv[2];
//  std::vector<int> varsNumbers = { 1 };
  std::vector<int> varsNumbers =
      { 1, 2, 3, 4, 5, 6, 7, 8, 9, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
          30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41 };
  auto c1 = new TCanvas();
  for (auto it : varsNumbers) {
    ReadDreamFile* DreamFile = new ReadDreamFile(6, 6);
    DreamFile->SetQuite();
    DreamFile->SetAnalysisFile(filename, prefix, TString::Format("%u", it));
    std::cout << "All fine, dont worry! \n";
    ProcessVariation(DreamFile,c1);
  }
  TFile* inFake = TFile::Open(
        "~/cernbox/pPb/v0offlineFix/woDetadPhi/HEP/CFOutput_pXi.root");
  TH1F* cfDef = (TH1F*)(inFake->Get("hCk_ReweightedMeV_1"))->Clone("hCk_RebinnedMeV_1");
  cfDef->SetLineWidth(3);
  cfDef->SetLineColor(2);
  c1->cd();
  cfDef->DrawCopy("SAME");
  TFile* defFile = TFile::Open(TString::Format("%s/CFpXiVariations_33.root",gSystem->pwd()), "RECREATE");
  defFile->cd();
  cfDef->Write();
  defFile->Close();
  c1->SaveAs("CFVars.pdf");

}

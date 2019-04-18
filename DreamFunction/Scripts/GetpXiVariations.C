#include "TROOT.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "ReadDreamFile.h"
#include "TCanvas.h"
#include "TF1.h"
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
  TH2F* SEMultPart = (TH2F*) ThePair->FindObject(
      "SEMultDist_Particle0_Particle4_clone");
  if (!SEMultPart) {
    std::cout << "SEMultDist_Particle0_Particle4_clone missing\n ";
    return;
  } else {
    Pair->GetPair()->SetSEMultDist(SEMultPart, "");
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
  TH2F* SEMultaPart = (TH2F*) aThePair->FindObject(
      "SEMultDist_Particle1_Particle5_clone");
  if (!SEMultaPart) {
    std::cout << "SEMultDist_Particle1_Particle5_clone missing\n ";
    return;
  } else {
    aPair->GetPair()->SetSEMultDist(SEMultaPart, "");
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
  TH1F* noRew = (TH1F*) inFake->Get("hCk_RebinnedMeV_1");
  TH1F* Rew = (TH1F*) inFake->Get("hCk_ReweightedMeV_1");
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
//  AddFakeReweighting(pXi, ApAXi);
  pXi->ShiftForEmpty(pXi->GetPair());
  ApAXi->ShiftForEmpty(ApAXi->GetPair());

  pXi->FixShift(pXi->GetPairShiftedEmpty(0), ApAXi->GetPairShiftedEmpty(0), ApAXi->GetFirstBin());
  ApAXi->FixShift(ApAXi->GetPairShiftedEmpty(0), pXi->GetPairShiftedEmpty(0), pXi->GetFirstBin());

//  pXi->FixShift(pXi->GetPair(), nullptr, 0.008, true);
//  ApAXi->FixShift(ApAXi->GetPair(), nullptr, 0.008, true);

  pXi->Rebin(pXi->GetPairFixShifted(0), 5);
  ApAXi->Rebin(ApAXi->GetPairFixShifted(0), 5);

//  pXi->Rebin(pXi->GetPairFixShifted(0), 5);
//  ApAXi->Rebin(ApAXi->GetPairFixShifted(0), 5);

  pXi->ReweightMixedEvent(pXi->GetPairRebinned(0), 0.2, 0.9);
  ApAXi->ReweightMixedEvent(ApAXi->GetPairRebinned(0), 0.2, 0.9);

  CF_pXi->SetPairs(pXi, ApAXi);
  CF_pXi->GetCorrelations();
  c1->cd();
  TH1F* CFDraw = CF_pXi->FindCorrelationFunction("hCk_ReweightedMeV_0");
//  ReweightManually(CFDraw);
  CFDraw->SetStats(false);
  CFDraw->GetXaxis()->SetRangeUser(0, 800);
  CFDraw->GetYaxis()->SetRangeUser(0.8, 3.);
  if (outNum == 1)
    CFDraw->DrawCopy();
  else
    CFDraw->DrawCopy("SAME");
  CF_pXi->WriteOutput(
      TString::Format("%s/CFpXiVariations_%u.root", gSystem->pwd(), outNum++));
}

void SystematicLimits(TH1F* lower, TH1F* upper, TH1F* def, TCanvas *c1) {
  lower = (TH1F*) def->Clone("upper");
  upper = (TH1F*) def->Clone("lower");
  lower->Reset();
  upper->Reset();
  TFile* SystErrFile = TFile::Open(
      "~/cernbox/SystematicsAndCalib/pPbRun2_MB_StreamFIX/C2totalsysPXi.root",
      "READ");
  TH1F* outputParam = (TH1F*) SystErrFile->Get(Form("SysParamPXi"));
  TF1 *RelSyst = new TF1("sys", "pol2", 0, 3);
  if (outputParam) {
    RelSyst->SetParameter(0, outputParam->GetBinContent(1));
    RelSyst->SetParameter(1, outputParam->GetBinContent(2));
    RelSyst->SetParameter(2, outputParam->GetBinContent(3));
  } else {
    std::cout << "Sytematics not available, exiting \n";
    return;
  }
  for (int ikS = 1; ikS <= def->FindBin(500); ++ikS) {
    float cK = def->GetBinContent(ikS);
    float err = cK
                * RelSyst->Eval(def->GetBinCenter(ikS) / 1000.);
    lower->SetBinContent(ikS, cK - err);
    lower->SetBinError(ikS, 1e-5);
    upper->SetBinContent(ikS, cK + err);
    upper->SetBinError(ikS, 1e-5);
  }
  c1->cd();
  def->SetLineWidth(3);
  upper->SetLineWidth(3);
  lower->SetLineWidth(3);
  def->SetLineColor(2);
  upper->SetLineColor(2);
  lower->SetLineColor(2);
  upper->SetLineStyle(5);
  lower->SetLineStyle(5);
  def->DrawCopy("SAME");
  upper->DrawCopy("SAME");
  lower->DrawCopy("SAME");
}

void SampleVariations(TH1F* inSample, TCanvas* c1) {
  static int outNumVars = 1;
  TF1 *myGauss = new TF1("myGauss", "gausn", 0, 5);
  TH1F* outHist = (TH1F*) inSample->Clone("hCk_RebinnedMeV_1");
  TFile* SystErrFile = TFile::Open(
      "~/cernbox/SystematicsAndCalib/pPbRun2_MB_StreamFIX/C2totalsysPXi.root",
      "READ");
  TH1F* outputParam = (TH1F*) SystErrFile->Get(Form("SysParamPXi"));
  TF1 *RelSyst = new TF1("sys", "pol2", 0, 3);
  TFile *outFile = TFile::Open(
      TString::Format("%s/CFpXiVariations_%u.root", gSystem->pwd(),
                      outNumVars++),
      "RECREATE");
  if (outputParam) {
    RelSyst->SetParameter(0, outputParam->GetBinContent(1));
    RelSyst->SetParameter(1, outputParam->GetBinContent(2));
    RelSyst->SetParameter(2, outputParam->GetBinContent(3));
  } else {
    std::cout << "Sytematics not available, exiting \n";
    return;
  }
  for (int ikS = 1; ikS <= inSample->FindBin(500); ++ikS) {
    myGauss->SetParameter(0, 1);
    myGauss->SetParameter(1, inSample->GetBinContent(ikS));
    myGauss->SetParameter(
        2,
        inSample->GetBinContent(ikS)
            * RelSyst->Eval(inSample->GetBinCenter(ikS) / 1000.));
    outHist->SetBinContent(ikS, myGauss->GetRandom());
  }
  c1->cd();
  outHist->GetXaxis()->SetRangeUser(0, 500);
  outHist->GetYaxis()->SetRangeUser(0.8, 3.);
  if (outNumVars == 2)
    outHist->DrawCopy();
  else
    outHist->DrawCopy("SAME");

  outHist->GetXaxis()->SetRangeUser(0, 2500);
  outFile->cd();
  outHist->Write();
  outFile->Close();
}

int main(int argc, char* argv[]) {
  gRandom = new TRandom3(0);
  gRandom->SetSeed(0);
  gROOT->ProcessLine("gErrorIgnoreLevel = 2001;");
  const char* filename = argv[1];
  const char* prefix = argv[2];
//  std::vector<int> varsNumbers = { 1 };
  std::vector<int> varsNumbers =
      { 1, 2, 3, 4, 5, 6, 7, 8, 9, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
          30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42 };
  TH1F* upper;
  TH1F* lower;
  auto c1 = new TCanvas();
  TFile* inFake = TFile::Open(
//      "~/cernbox/pPb/v0offlineFix/woDetadPhi/HEP/CFOutput_pXi.root","read");
      "~/cernbox/HM13TeV/AnalysisData/ClosePairRej/SelectedPairs/CFOutput_pXi.root","read");
  TH1F* cfDef = (TH1F*) (inFake->Get("hCk_ReweightedMeV_1")->Clone("hCk_ReweightedMeV_0"));
  cfDef->SetName("hCk_ReweightedMeV_0");
  cfDef->SetTitle("hCk_ReweightedMeV_0");
  cfDef->SetLineColor(kPink);
  for (auto it : varsNumbers) {
//    SampleVariations(cfDef, c1);
    ReadDreamFile* DreamFile = new ReadDreamFile(6, 6);
    DreamFile->SetQuite();
    DreamFile->SetAnalysisFile(filename, prefix, TString::Format("%u", it));
    ProcessVariation(DreamFile,c1);
  }

//  SystematicLimits(lower,upper,cfDef,c1);
  c1->cd();
  cfDef->Draw("SAME");
  TFile* defFile = TFile::Open(
      TString::Format("%s/CFpXiVariations_34.root", gSystem->pwd()),
      "RECREATE");
  defFile->cd();
  cfDef->Write();
  defFile->Close();

  c1->SaveAs("CFVars.pdf");
}


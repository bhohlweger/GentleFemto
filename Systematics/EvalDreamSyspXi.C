#include "CATSInput.h"
#include "DreamPlot.h"
#include "DreamSystematics.h"
#include "ForgivingReader.h"
#include "CandidateCounter.h"
#include "TCanvas.h"
#include <iostream>
#include "TSystem.h"
#include "DecayQA.h"

void EvalDreamSystematics(TString InputFile, TString prefix,
                          float upperFitRange) {
  gROOT->ProcessLine("gErrorIgnoreLevel = 3001");
  std::cout << InputFile.Data() << std::endl;
  DreamPlot::SetStyle(false, true);
  auto CATSinput = new CATSInput();
  CATSinput->SetNormalization(0.500, 0.700);
  CATSinput->SetFixedkStarMinBin(true, 0.);
  const int rebin = 5;
  auto counter = new CandidateCounter();

  ReadDreamFile* DreamFile = new ReadDreamFile(6, 6);
  DreamFile->SetAnalysisFile(InputFile.Data(), prefix, "0");

  ForgivingReader* ForgivingFileDef = new ForgivingReader(InputFile.Data(),
                                                          prefix, "0");
  counter->SetNumberOfCandidates(ForgivingFileDef);
  const int nTracks = counter->GetNumberOfTracks();
  const int nCascades = counter->GetNumberOfCascades();
  counter->ResetCounter();
  DreamDist* pXi = DreamFile->GetPairDistributions(0, 4, "");
  DreamDist* ApAXi = DreamFile->GetPairDistributions(1, 5, "");
  DreamCF* CFpXiDef = CATSinput->ObtainCFSyst(rebin, "pXiVar0", pXi, ApAXi);
  const int pairCountsDefault = CFpXiDef->GetFemtoPairs(0, 0.2);
  DreamSystematics protonXi(DreamSystematics::pXi);
  protonXi.SetUpperFitRange(upperFitRange);
  protonXi.SetBarlowUpperRange(400);
  protonXi.SetDefaultHist(CFpXiDef, "hCk_ReweightedpXiVar0MeV_1");

  TH1F* PurityXi = new TH1F("PurityXiVar", "PurityXiVar", 45, 0.5, 44.5);
  TH1F* PurityAXi = new TH1F("PurityAXiVar", "PurityAXiVar", 45, 0.5, 44.5);

  const TString currentWordDir = gSystem->pwd();
  std::cout << "Current working directory: " << currentWordDir << std::endl;
  int outCounter = 1;
  for (int i = 1; i <= 44; ++i) {
    ReadDreamFile* DreamVarFile = new ReadDreamFile(6, 6);
    DreamVarFile->SetAnalysisFile(InputFile.Data(), prefix, Form("%u", i));
    TString VarName = TString::Format("ppVar%u", outCounter);
    DreamCF* CFpXiVar = CATSinput->ObtainCFSyst(
        rebin, VarName.Data(), DreamVarFile->GetPairDistributions(0, 4, ""),
        DreamVarFile->GetPairDistributions(1, 5, ""));
    DreamCF* CFpXiOut = CATSinput->ObtainCFSyst(
        rebin, VarName.Data(), DreamVarFile->GetPairDistributions(0, 4, ""),
        DreamVarFile->GetPairDistributions(1, 5, ""));
    int femtoPairVar = CFpXiVar->GetFemtoPairs(0, 0.2);
    float relDiff = (femtoPairVar - pairCountsDefault)
        / (float) pairCountsDefault;
    if (TMath::Abs(relDiff) > 0.2) {
      continue;
    }

    TString QADirectory = TString::Format("%s/Var_%u/", currentWordDir.Data(), i);
    gSystem->MakeDirectory(QADirectory.Data());
    gSystem->cd(QADirectory.Data());

    protonXi.SetVarHist(CFpXiVar,
                        TString::Format("Reweighted%sMeV_1", VarName.Data()));
    TString VarString = TString::Format("%u", i);
    ForgivingReader* ForgivingFile = new ForgivingReader(InputFile.Data(),
                                                         prefix,
                                                         VarString.Data());
    DecayQA* cascQA = new DecayQA("#Xi^{-}", "#pi#Lambda");
    cascQA->SetCanvasDivisions(4, 4);
    cascQA->SetInvMasspTStartBin(2);
    cascQA->SetIMHistoScale(2.5, 0.8, 0.45);
    cascQA->SetDecayCuts(ForgivingFile->GetCascadeCuts());
    cascQA->SetAntiDecayCuts(ForgivingFile->GetAntiCascadeCuts());
    cascQA->SetRangesFitting(1.31, 1.33, 1.285, 1.365);
    cascQA->InvariantMassXiMinBooking(1.317, 1.327);
    PurityXi->SetBinContent(i, cascQA->GetIntegratedPurity(0));
    PurityXi->SetBinError(i, cascQA->GetIntegratedPurity(0)*0.001);
    PurityAXi->SetBinContent(i, cascQA->GetIntegratedPurity(1));
    PurityAXi->SetBinError(i, cascQA->GetIntegratedPurity(1)*0.001);
    counter->SetNumberOfCandidates(ForgivingFile);
    protonXi.SetPair(pairCountsDefault, CFpXiVar->GetFemtoPairs(0, 0.2));
    protonXi.SetParticles(nTracks, nCascades, counter->GetNumberOfTracks(),
                          counter->GetNumberOfCascades());
    TString OutputDirCF = TString::Format("%s/CF_pXi_Var%u.root",
                                          currentWordDir.Data(), outCounter++);
    CFpXiOut->WriteOutput(OutputDirCF.Data());
    counter->ResetCounter();
    gSystem->cd(currentWordDir.Data());

  }
  protonXi.EvalSystematics();
  protonXi.EvalDifferenceInPairs();
  protonXi.EvalDifferenceInParticles();

  protonXi.GetDefault()->SetLineColor(kBlack);
  protonXi.GetDefault()->SetLineWidth(3);
  protonXi.GetDefault()->SetMarkerColor(kBlack);
  protonXi.GetDefault()->SetMarkerStyle(24);
  protonXi.GetDefault()->SetMarkerSize(1.2);

  protonXi.WriteOutput();
  CFpXiDef->WriteOutput(
      TString::Format("%s/CF_pXi_Var0.root", gSystem->pwd()).Data());
  TCanvas* purity = new TCanvas("puritiesVar","puritiesVar",0,0,1200,1000);
  purity->Divide(2,1);
  purity->cd(1);
  PurityXi->GetXaxis()->SetTitle("Variation");
  PurityXi->GetYaxis()->SetTitle("Purity");
  PurityXi->GetYaxis()->SetRangeUser(0.9,1.);
  PurityXi->SetMarkerColor(kBlue);
  PurityXi->SetLineColor(kBlue);
  PurityXi->SetMarkerStyle(21);
  PurityXi->SetMarkerSize(1.2);
  PurityXi->SetLineWidth(3);
  PurityXi->Draw();
  purity->cd(2);
  PurityAXi->GetXaxis()->SetTitle("Variation");
  PurityAXi->GetYaxis()->SetTitle("Purity");
  PurityAXi->GetYaxis()->SetRangeUser(0.9,1.);
  PurityAXi->SetMarkerStyle(22);
  PurityAXi->SetMarkerSize(1.2);
  PurityAXi->SetLineWidth(3);
  PurityAXi->Draw();
  purity->SaveAs("PurityVar.pdf");
  TFile* purityFile = TFile::Open("Purityies.root","recreate");
  purityFile->cd();
  PurityXi->Write();
  PurityAXi->Write();
  purityFile->Write();
  purityFile->Close();


}

int main(int argc, char* argv[]) {
  EvalDreamSystematics(argv[1], argv[2], atof(argv[3]));

  return 1;
}

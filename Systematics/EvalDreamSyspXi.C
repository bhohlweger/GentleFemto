#include "CATSInput.h"
#include "DreamPlot.h"
#include "DreamSystematics.h"
#include "ForgivingReader.h"
#include "CandidateCounter.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "DecayQA.h"
#include "AnalyseProXi.h"
#include <iostream>
#include <string>

void EvalDreamSystematics(TString InputFile, TString prefix,
                          float upperFitRange) {
  //to do: get # of femto pairs out of analysis pxi and fix the suffix ...
  gROOT->ProcessLine("gErrorIgnoreLevel = 3001");
  std::cout << InputFile.Data() << std::endl;
  DreamPlot::SetStyle();
  AnalyseProXi* ana = new AnalyseProXi(1000, 0.95);
  ana->SetAnalysisFile(InputFile, prefix);
  ana->Default();

  ReadDreamFile* DreamFile = new ReadDreamFile(6, 6);
  DreamFile->SetAnalysisFile(InputFile.Data(), prefix, "0");

  DreamDist* pxi = DreamFile->GetPairDistributions(0, 4, "");
  DreamDist* ApAxi = DreamFile->GetPairDistributions(1, 5, "");
  const int pairCountsDefault = pxi->GetFemtoPairs(0, 0.200)
      + ApAxi->GetFemtoPairs(0, 0.200);

  DreamSystematics protonXi(DreamSystematics::pXi);
  protonXi.SetUpperFitRange(800);
  protonXi.SetBarlowUpperRange(400);
  protonXi.SetDefaultHist(ana->GetVariation(0, false));  //set the default histogram from the ana

  TH1F* PurityXi = new TH1F("PurityXiVar", "PurityXiVar", 45, 0.5, 44.5);
  TH1F* PurityAXi = new TH1F("PurityAXiVar", "PurityAXiVar", 45, 0.5, 44.5);

  const TString currentWordDir = gSystem->pwd();
  std::cout << "Current working directory: " << currentWordDir << std::endl;
  int outCounter = 1;

  for (int i = 1; i <= 44; ++i) {
    ReadDreamFile* DreamFileVar = new ReadDreamFile(6, 6);
    std::cout << "TString::Format(, i).Data(): " << TString::Format("%i", i).Data() << std::endl;
    DreamFileVar->SetAnalysisFile(InputFile.Data(), prefix,
                                  TString::Format("%i", i).Data());
    DreamDist* pxiVar = DreamFileVar->GetPairDistributions(0, 4, "");
    DreamDist* ApAxiVar = DreamFileVar->GetPairDistributions(1, 5, "");
    int femtoPairVar = pxiVar->GetFemtoPairs(0, 0.200)
        + ApAxiVar->GetFemtoPairs(0, 0.200);

    delete DreamFileVar;
    delete pxiVar;
    delete ApAxiVar;

    float relDiff = (femtoPairVar - pairCountsDefault)
        / (float) pairCountsDefault;
    if (TMath::Abs(relDiff) > 0.2) {
      continue;
    }
    std::string s = std::to_string(i);
    char const *pchar = s.c_str();
    ana->SetAnalysisFile(InputFile, prefix, pchar);
    protonXi.SetVarHist(ana->GetVariation(outCounter, false));

    TString QADirectory = TString::Format("%s/Var_%u/", currentWordDir.Data(),
                                          i);
    gSystem->MakeDirectory(QADirectory.Data());
    gSystem->cd(QADirectory.Data());
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
    PurityXi->SetBinError(i, cascQA->GetIntegratedPurity(0) * 0.001);
    PurityAXi->SetBinContent(i, cascQA->GetIntegratedPurity(1));
    PurityAXi->SetBinError(i, cascQA->GetIntegratedPurity(1) * 0.001);
    protonXi.SetPair(pairCountsDefault, femtoPairVar);
    TString OutputDirCF = TString::Format("%s/CF_pXi_Var%u.root",
                                          currentWordDir.Data(), outCounter++);
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

  TCanvas* purity = new TCanvas("puritiesVar", "puritiesVar", 0, 0, 1200, 1000);
  purity->Divide(2, 1);
  purity->cd(1);
  PurityXi->GetXaxis()->SetTitle("Variation");
  PurityXi->GetYaxis()->SetTitle("Purity");
  PurityXi->GetYaxis()->SetRangeUser(0.9, 1.);
  PurityXi->SetMarkerColor(kBlue);
  PurityXi->SetLineColor(kBlue);
  PurityXi->SetMarkerStyle(21);
  PurityXi->SetMarkerSize(1.2);
  PurityXi->SetLineWidth(3);
  PurityXi->Draw();
  purity->cd(2);
  PurityAXi->GetXaxis()->SetTitle("Variation");
  PurityAXi->GetYaxis()->SetTitle("Purity");
  PurityAXi->GetYaxis()->SetRangeUser(0.9, 1.);
  PurityAXi->SetMarkerStyle(22);
  PurityAXi->SetMarkerSize(1.2);
  PurityAXi->SetLineWidth(3);
  PurityAXi->Draw();
  purity->SaveAs("PurityVar.pdf");
  TFile* purityFile = TFile::Open("Purityies.root", "recreate");
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

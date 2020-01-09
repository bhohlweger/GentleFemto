#include "AnalyseProXi.h"
#include "CATSInput.h"
#include "CATSLambdaParam.h"
#include "DreamSystematics.h"
#include "LambdaGami.h"
#include "SideBandFit.h"
#include "TidyCats.h"

#include "TApplication.h"
#include "TF1.h"
#include "TSystem.h"

#include "DLM_CkDecomposition.h"

int main(int argc, char *argv[]) {
//  TApplication *app = new TApplication("myapp",argc, argv);
  TFile* inFileDist = TFile::Open(
      "~/cernbox/HM13TeV/AnalysisData/InjectedMC.root", "read");
  inFileDist->ls();
  TH1F* cloneme = (TH1F*) inFileDist->Get("SETruthNorm");
  TH1F* SEDist = (TH1F*) cloneme->Clone("SERaw");
  TH1F* SEDistRebin = (TH1F*) cloneme->Clone("SERawRebin");
  cloneme = (TH1F*) inFileDist->Get("METruthNorm");
  TH1F* MEDist = (TH1F*) cloneme->Clone("MERaw");
  TH1F* MEDistRebin = (TH1F*) cloneme->Clone("MERawRebin");
  cloneme = (TH1F*) inFileDist->Get("cf");
  TH1F* CFDist = (TH1F*) cloneme->Clone("cfRaw");
  cloneme = (TH1F*) inFileDist->Get("ProjYNorm");
  TH1F* MEReco = (TH1F*) cloneme->Clone("MEReco");

  SEDistRebin->Rebin(4);
  MEDistRebin->Rebin(4);

  SEDistRebin->Scale(
      1. / SEDistRebin->Integral(SEDistRebin->FindBin(0.2), SEDistRebin->FindBin(0.4)));
  MEDistRebin->Scale(
      1. / MEDistRebin->Integral(MEDistRebin->FindBin(0.2), MEDistRebin->FindBin(0.4)));
  TH1F* CFRebin = (TH1F*) SEDistRebin->Clone("cfRebin");
  CFRebin->Divide(MEDistRebin);


  MomentumGami* momGami = new MomentumGami(0.997);
  DreamPlot::SetStyle();
  auto* canOne = new TCanvas("c1", "c1");
  TString CalibBaseDir =
      "/home/schmollweger/cernbox/HM13TeV/AnalysisData/1436_AODXioton/ResolutionpXi.root";
  TFile* inFile = TFile::Open(CalibBaseDir.Data(), "read");
  if (!inFile) {
    std::cout << "No Infile set, no Momentum resolution set, RIP \n";
  } else {
    TH2F* momReso = (TH2F*) inFile->Get("FiveMeV");
    momGami->SetResolution(momReso, 1);
  }

  //Smear
  TH1F* SEFolded = momGami->Fold(SEDist);
  TH1F* MEFolded = momGami->Fold(MEDist);

  TH1F* RatioFolding = (TH1F*)MEFolded->Clone("MERatioFolded");
  RatioFolding->Divide(MEReco);
//  momGami->Fold(SEDist);
  SEFolded->Scale(
      1. / SEFolded->Integral(SEFolded->FindBin(0.2), SEFolded->FindBin(0.4)));
  MEFolded->Scale(
      1. / MEFolded->Integral(MEFolded->FindBin(0.2), MEFolded->FindBin(0.4)));
  TH1F* CFFolded = (TH1F*) SEFolded->Clone("cffolded");
  CFFolded->Divide(MEFolded);


  TH1F* SEFoldedRebin = (TH1F*)SEFolded->Clone("SEFoldedRebin");
  SEFoldedRebin->Rebin(4);
  TH1F* MEFoldedRebin = (TH1F*)MEFolded->Clone("MEFoldedRebin");
  MEFoldedRebin->Rebin(4);

  SEFoldedRebin->Scale(
      1. / SEFoldedRebin->Integral(SEFoldedRebin->FindBin(0.2), SEFoldedRebin->FindBin(0.4)));
  MEFoldedRebin->Scale(
      1. / MEFoldedRebin->Integral(MEFoldedRebin->FindBin(0.2), MEFoldedRebin->FindBin(0.4)));
  TH1F* CFFoldedRebin = (TH1F*)SEFoldedRebin->Clone("CFFoldedRebin");
  CFFoldedRebin->Divide(MEFoldedRebin);

  //Unfold
  TH1F* SEUnfolded = momGami->UnfoldviaRooResp(SEFolded);
  TH1F* MEUnfolded = momGami->UnfoldviaRooResp(MEFolded);

  TH1F* RatioME = (TH1F*)MEUnfolded->Clone("MERatio");
  RatioME->Divide(MEDist);

  SEUnfolded->Rebin(4);
  MEUnfolded->Rebin(4);
  SEUnfolded->Scale(
      1.
          / SEUnfolded->Integral(SEUnfolded->FindBin(0.2),
                                 SEUnfolded->FindBin(0.4)));
  MEUnfolded->Scale(
      1.
          / MEUnfolded->Integral(MEUnfolded->FindBin(0.2),
                                 MEUnfolded->FindBin(0.4)));

  TH1F* CFUnfolded = (TH1F*) SEUnfolded->Clone("cfunfolded");
  CFUnfolded->Divide(MEUnfolded);

  TH1F* CFRatioRawUnfolded = (TH1F*)CFUnfolded->Clone("RatioRawUnfoldedCF");
  CFRatioRawUnfolded->Divide(CFRebin);

  TH1F* CFRatioRawFolded = (TH1F*)CFFoldedRebin->Clone("RatioRawFoldedCF");
  CFRatioRawFolded->Divide(CFRebin);

  TFile* output = TFile::Open("UnfoldererBayes_5.root", "recreate");
  output->cd();

  SEDist->Write();
  MEDist->Write();
  CFDist->Write();

  SEDistRebin->Write();
  MEDistRebin->Write();
  CFRebin->Write();

  SEFolded->Write();
  MEFolded->Write();
  CFFolded->Write();

  SEFoldedRebin->Write();
  MEFoldedRebin->Write();
  CFFoldedRebin->Write();

  SEUnfolded->Write();
  MEUnfolded->Write();
  CFUnfolded->Write();
  RatioME->Write();
  RatioFolding->Write();
  CFRatioRawUnfolded->Write();
  CFRatioRawFolded->Write();

  momGami->GetQAList()->Write("QAList", 1);
  output->Close();
//  app->Run();
  return 0;
}


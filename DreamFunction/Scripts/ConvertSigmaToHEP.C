#include "DreamHEP.h"
#include "TFile.h"

void ConvertProtonProton(TString path) {
  DreamHEP hep;
  hep.SetRootS(13000);
  hep.SetMaxkStar(225);
  auto CFFile = TFile::Open(Form("%s/CFOutput_pp.root", path.Data()));
  auto hepstogram = (TH1F*) CFFile->Get("hCk_ReweightedMeV_1");

  TFile* CFFileSyst = TFile::Open("~/cernbox/SystematicsAndCalib/ppRun2_HM/Systematics_pp.root");
  TF1* sysParam = (TF1*) CFFileSyst->Get("SystError");

  auto heptographsystErr = hep.GetSystErrHist(hepstogram, sysParam);
  hep.printTH1HEPdata(hepstogram, heptographsystErr, Form("%s/PP", path.Data()));
}

void ConvertProtonSigma(TString path) {
  DreamHEP hep;
  hep.SetRootS(13000);
  hep.SetMaxkStar(365);
  auto CFFile = TFile::Open(Form("%s/data/CFOutput_pSigma.root", path.Data()));
  auto hepstograph = (TGraphAsymmErrors*) CFFile->Get("Graph_from_hCk_Reweighted_3MeV");

  auto CFFileSyst = TFile::Open(
      Form("%s/systematics/Systematics_pSigma0.root", path.Data()));
  TF1* sysParam = (TF1*) CFFileSyst->Get("SystError");

  auto heptographsystErr = hep.GetSystErrHist(hepstograph, sysParam);
  hep.printTGAsymmHEPdata(hepstograph, heptographsystErr, Form("%s/PSigma0", path.Data()));
}

int main(int argc, char *argv[]) {
  ConvertProtonProton(argv[1]);
  ConvertProtonSigma(argv[2]);
  return 0;
}

#include "DreamHEP.h"
#include "TFile.h"
#include <iostream>
#include "TF1.h"
void ConvertTH1ToHEP(TString inHistFile, TString HistName, TString SysFile,
                     TString outname);
int main(int argc, char *argv[]) {
  ConvertTH1ToHEP(argv[1], argv[2], argv[3], argv[4]);
  return 0;
}

void ConvertTH1ToHEP(TString inHistFile, TString HistName, TString SysFile,
                     TString outname) {
  DreamHEP* hep;
  TFile* inputHist = TFile::Open(inHistFile.Data(), "read");
  if (!inputHist) {
    std::cout << "No input file for Histograms \n";
    return;
  }
  TH1F* hepstogram = nullptr;
  hepstogram = (TH1F*) inputHist->Get(HistName.Data());
  if (!hepstogram) {
    std::cout << "Histogram not found, what is in this file? \n";
    inputHist->ls();
    return;
  }
  TFile* inputSys = TFile::Open(SysFile.Data(), "read");
  if (!inputSys) {
    std::cout << "No input file for Systematic errors \n";
    return;
  }
  TString PairName;
  if (SysFile.Contains("PP")) {
    PairName = "PP";
  } else if (SysFile.Contains("PL")) {
    PairName = "PL";
  } else if (SysFile.Contains("LL")) {
    PairName = "LL";
  } else if (SysFile.Contains("PXi")) {
    PairName = "PXi";
  } else {
    std::cout << "Pair Name not detected! \n";
    return;
  }
  std::cout << "PairName correct? " << PairName.Data() << std::endl;
  TH1F* outputParam = nullptr;
  outputParam = (TH1F*) inputSys->Get(Form("SysParam%s", PairName.Data()));
  if (!outputParam) {
    std::cout << "No output Parameter, Pair name correct? \n";
    inputSys->ls();
    return;
  }
  TGraphErrors* HistWithSysErr = new TGraphErrors();
  float convFac = HistName.Contains("MeV") ? 1000. : 1.;
  std::cout << "Conversion factor correct? " << convFac << std::endl;
  TF1 *RelSyst = new TF1("sys", "pol2", 0, 3);
  RelSyst->SetParameter(0, outputParam->GetBinContent(1));
  RelSyst->SetParameter(1, outputParam->GetBinContent(2));
  RelSyst->SetParameter(2, outputParam->GetBinContent(3));
  int NumSEB = hepstogram->FindBin(0.5 * convFac);
  for (int iBin = 0; iBin < NumSEB; iBin++) {
    const float x = hepstogram->GetBinCenter(iBin + 1);
    const float y = hepstogram->GetBinContent(iBin + 1);
    HistWithSysErr->SetPoint(iBin, x, y);
    HistWithSysErr->SetPointError(iBin, 0, y * RelSyst->Eval(x / convFac));
  }
  hep->printTH1HEPdata(hepstogram, HistWithSysErr, outname.Data());
}


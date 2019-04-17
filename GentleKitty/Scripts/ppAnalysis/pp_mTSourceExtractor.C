#include "TFile.h"
#include "TGraphErrors.h"
#include "TNtuple.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <iostream>
void mTRadiusExtractor(const char *inputFile, float &radius, float &radErr) {
  TFile* inFile = TFile::Open(inputFile, "READ");
  if(!inFile) {
    std::cout << "No input found for " << inputFile << ", exiting \n";
    return;
  }
  TNtuple* mTSysVar = (TNtuple*) inFile->Get("ntResult");
  auto histRad = new TH1D("hRad", "hRad", 200, 0.4, 1.2);
  mTSysVar->Draw("Radius_pp>>hRad");
  radius = histRad->GetMean();
  radErr = histRad->GetRMS();
  delete histRad;
}

void mTRadiusPlot(const char* inputFolder, const char* avgmTFile) {
  TFile* mTFile = TFile::Open(avgmTFile, "READ");
  TGraphErrors* avgmT = (TGraphErrors*) mTFile->Get("AveragemT");
  TGraphErrors* mTRadius = new TGraphErrors();
  for (int imT = 0; imT < 7; ++imT) {
    double mT, dummy;
    float radius, radErr;
    TString inputFile = TString::Format("%s/OutFileVarpp_%u.root",
                                        inputFolder, imT);
    avgmT->GetPoint(imT, dummy, mT);
    std::cout << avgmT->GetErrorY(imT) << std::endl;
    mTRadiusExtractor(inputFile.Data(), radius, radErr);
    std::cout << "radius: " << radius << " radErr: " << radErr << std::endl;
    mTRadius->SetPoint(imT, mT, radius);
    mTRadius->SetPointError(imT, avgmT->GetErrorY(imT), radErr);
  }
  TFile* output = TFile::Open("mTRad.root", "RECREATE");
  output->cd();
  mTRadius->SetName("mTRadius");
  mTRadius->Write();
  output->Close();
}

int main(int argc, char *argv[]) {
  mTRadiusPlot(argv[1], argv[2]);
  return 0;
}

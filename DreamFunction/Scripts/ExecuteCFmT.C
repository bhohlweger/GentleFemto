#include "DreamKayTee.h"
#include "ReadDreamFile.h"
#include "TROOT.h"
#include "TSystem.h"

int main(int argc, char* argv[]) {
  const char* filename = argv[1];
  const char* CalibPath = argv[2];
  const char* prefix = argv[3];
  const char* addon = (argv[4]) ? argv[4] : "";

  ReadDreamFile* DreamFile = new ReadDreamFile(6, 6);
  DreamFile->SetQuite();
  TString foldername = filename;
  foldername.ReplaceAll("AnalysisResults.root", "");
  DreamKayTee* mTppDists;
  DreamFile->ReadmTHistos(filename, prefix, addon);
  mTppDists = DreamFile->GetmTPairDistributions(0, 0, 1, 1);

  std::vector<float> mTppBins = { 1.02, 1.14, 1.2, 1.26, 1.38, 1.56, 1.86, 4.5};
//  std::vector<float> mTBins = { 1.01,1.13,1.25,1.55,4.5 };
  TString CalibPP = Form("%s/CFOutput_pp.root", CalibPath);
  mTppDists->SetSEMEReweightingRatio(CalibPP, "1", false);
  mTppDists->SetKayTeeBins(mTppBins);
  mTppDists->SetNormalization(0.24, 0.34);
  mTppDists->SetRebin( { 2 });
  mTppDists->FixShift({true,true,true,true,true,true,true},{0.004,0.004,0.008,0.004,0.008,0.008,0.012});
  mTppDists->ObtainTheCorrelationFunction(gSystem->pwd(), prefix, "pp");

  DreamKayTee* mTpLDists;
  DreamFile->ReadmTHistos(filename, prefix, addon);
  mTpLDists = DreamFile->GetmTPairDistributions(0, 2, 1, 3);

//  std::vector<float> mTpLBins = { 1.02, 1.26, 1.32, 1.44, 1.62, 1.68, 4.5 };
//
//  TString CalibPL = Form("%s/CFOutput_pL.root",CalibPath);
//  mTpLDists->SetSEMEReweightingRatio(CalibPL, "1", true);
//  mTpLDists->SetKayTeeBins(mTpLBins);
//  mTpLDists->SetNormalization(0.24, 0.34);
//  mTpLDists->ObtainTheCorrelationFunction(gSystem->pwd(), prefix, "pL");
//
  DreamKayTee* mTpXiDists;
  DreamFile->ReadmTHistos(filename, prefix, addon);
  mTpXiDists = DreamFile->GetmTPairDistributions(0, 4, 1, 5);

//  std::vector<float> mTpXiBins = { 1.08, 4.5 };
  std::vector<float> mTpXiBins = { 1.08, 1.74, 4.5 };

  TString CalibpXi = Form("%s/CFOutput_pXi.root",CalibPath);
  mTpXiDists->SetSEMEReweightingRatio(CalibpXi, "1", true);
  mTpXiDists->SetKayTeeBins(mTpXiBins);
  mTpXiDists->SetNormalization(0.24, 0.34);
  mTpXiDists->SetRebin({5});
  mTpXiDists->FixShift({false,true},{0.000,0.004});
  mTpXiDists->ObtainTheCorrelationFunction(gSystem->pwd(), prefix, "pXi");
//
  return 0;
}

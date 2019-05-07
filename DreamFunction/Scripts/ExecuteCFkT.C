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
  DreamKayTee* kTppDists;
  DreamFile->ReadkTHistos(filename, prefix, addon);
  kTppDists = DreamFile->GetkTPairDistributions(0, 0, 1, 1);

  std::vector<float> kTppBins = { 0.13,0.3,0.5,0.7,0.9,1.2};
//  std::vector<float> kTBins = { 1.01,1.13,1.25,1.55,4.5 };
  TString CalibPP = Form("%s/CFOutput_pp.root", CalibPath);
  kTppDists->SetSEMEReweightingRatio(CalibPP, "1", false);
  kTppDists->SetKayTeeBins(kTppBins);
  kTppDists->SetNormalization(0.35,0.4);
//  kTppDists->SetRebin( { 2 });
//  kTppDists->FixShift({true,true,true,true,true,true,true},{0.004,0.004,0.008,0.004,0.008,0.008,0.012});
  kTppDists->ObtainTheCorrelationFunction(gSystem->pwd(), prefix, "pp");

//  DreamKayTee* kTpLDists;
//  DreamFile->ReadkTHistos(filename, prefix, addon);
//  kTpLDists = DreamFile->GetkTPairDistributions(0, 2, 1, 3);
//
//  std::vector<float> kTpLBins = { 1.02, 1.26, 1.32, 1.44, 1.62, 1.68, 4.5 };
//
//  TString CalibPL = Form("%s/CFOutput_pL.root",CalibPath);
//  kTpLDists->SetSEMEReweightingRatio(CalibPL, "1", true);
//  kTpLDists->SetKayTeeBins(kTpLBins);
//  kTpLDists->SetNormalization(0.24, 0.34);
//  kTpLDists->ObtainTheCorrelationFunction(gSystem->pwd(), prefix, "pL");
//
//  DreamKayTee* kTpXiDists;
//  DreamFile->ReadkTHistos(filename, prefix, addon);
//  kTpXiDists = DreamFile->GetkTPairDistributions(0, 4, 1, 5);
//
////  std::vector<float> kTpXiBins = { 1.08, 4.5 };
//  std::vector<float> kTpXiBins = { 1.08, 1.74, 4.5 };
//
//  TString CalibpXi = Form("%s/CFOutput_pXi.root",CalibPath);
//  kTpXiDists->SetSEMEReweightingRatio(CalibpXi, "1", true);
//  kTpXiDists->SetKayTeeBins(kTpXiBins);
//  kTpXiDists->SetNormalization(0.24, 0.34);
//  kTpXiDists->SetRebin({5});
//  kTpXiDists->FixShift({false,true},{0.000,0.004});
//  kTpXiDists->ObtainTheCorrelationFunction(gSystem->pwd(), prefix, "pXi");
//
  return 0;
}

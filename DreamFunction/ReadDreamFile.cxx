/*
 * ReadDreamFile.cxx
 *
 *  Created on: Aug 21, 2018
 *      Author: hohlweger
 */
#include "ReadDreamFile.h"
#include <iostream>
ReadDreamFile::ReadDreamFile(int nPart1, int nPart2)
    : fNPart1(nPart1),
      fNPart2(nPart2),
      fSE(nullptr),
      fSEMult(nullptr),
      fSEkT(nullptr),
      fSEmT(nullptr),
      fME(nullptr),
      fMEMult(nullptr),
      fMEkT(nullptr),
      fMEmT(nullptr) {
  TH1::AddDirectory(kFALSE);
  TH2::AddDirectory(kFALSE);
}

ReadDreamFile::~ReadDreamFile() {
  for (int iPart1 = 0; iPart1 < fNPart1; ++iPart1) {
    for (int iPart2 = iPart1; iPart2 < fNPart2; ++iPart2) {
      if (fSE && fSE[iPart1][iPart2])
        delete fSE[iPart1][iPart2];
      if (fSEMult && fSEMult[iPart1][iPart2])
        delete fSEMult[iPart1][iPart2];
      if (fSEkT && fSEkT[iPart1][iPart2])
        delete fSEkT[iPart1][iPart2];
      if (fSEmT && fSEmT[iPart1][iPart2])
        delete fSEmT[iPart1][iPart2];
      if (fME && fME[iPart1][iPart2])
        delete fME[iPart1][iPart2];
      if (fMEMult && fMEMult[iPart1][iPart2])
        delete fMEMult[iPart1][iPart2];
      if (fMEkT && fMEkT[iPart1][iPart2])
        delete fMEkT[iPart1][iPart2];
      if (fMEmT && fMEmT[iPart1][iPart2])
        delete fMEmT[iPart1][iPart2];
    }
  }
}

void ReadDreamFile::SetAnalysisFile(const char* PathAnalysisFile,
                                    const char* Prefix, const char* Addon) {
  TFile* _file0 = TFile::Open(PathAnalysisFile, "READ");

  TDirectoryFile *dirResults = (TDirectoryFile*) (_file0->FindObjectAny(
      Form("%sResults%s", Prefix, Addon)));
  TList *Results;
  dirResults->GetObject(Form("%sResults%s", Prefix, Addon), Results);
  ExtractResults(Results);
  _file0->Close();
}

void ReadDreamFile::SetSigmaAnalysisFile(const char* PathAnalysisFile,
                                         const char* suffix) {
  auto file = TFile::Open(PathAnalysisFile);
  TString name = "Sigma0_Femto_";
  name += suffix;
  TDirectory *dir = file->GetDirectory(name);
  name = "femto_";
  name += suffix;
  auto histoList = (TList *) dir->Get(name);
  auto Results = (TList*) histoList->FindObject("Results");
  ExtractResults(Results);
  TIter next(Results);
  TObject *obj = nullptr;
  while (obj = next()) {
    TList *list = dynamic_cast<TList*>(obj);
    if (list)
      list->Delete();
  }
  Results->Delete();
  histoList->Delete();
  dir->Close();
  file->Close();
  delete file;
}

void ReadDreamFile::ExtractResults(const TList *Results) {
  TList *PartList;

  fSE = new TH1F**[fNPart1];
  fSEMult = new TH2F**[fNPart1];
  fME = new TH1F**[fNPart1];
  fMEMult = new TH2F**[fNPart1];

  for (int iPart1 = 0; iPart1 < fNPart1; ++iPart1) {
    fSE[iPart1] = new TH1F*[fNPart2];
    fSEMult[iPart1] = new TH2F*[fNPart2];
    fME[iPart1] = new TH1F*[fNPart2];
    fMEMult[iPart1] = new TH2F*[fNPart2];

    for (int iPart2 = iPart1; iPart2 < fNPart2; ++iPart2) {
      TString FolderName = Form("Particle%i_Particle%i", iPart1, iPart2);
      PartList = (TList*) Results->FindObject(FolderName.Data());
      fSE[iPart1][iPart2] = nullptr;
      auto hist1D = (TH1F*) PartList->FindObject(
          Form("SEDist_%s", FolderName.Data()));
      fSE[iPart1][iPart2] = (TH1F*) hist1D->Clone(
          Form("%s_clone", hist1D->GetName()));
      if (!fSE[iPart1][iPart2]) {
        std::cout << "SE Histogramm missing from " << FolderName.Data()
                  << std::endl;
      }
      fSE[iPart1][iPart2]->Sumw2();

      fSEMult[iPart1][iPart2] = nullptr;
      auto hist2D = (TH2F*) PartList->FindObject(
          Form("SEMultDist_%s", FolderName.Data()));
      fSEMult[iPart1][iPart2] = (TH2F*) hist2D->Clone(
          Form("%s_clone", hist2D->GetName()));
      if (!fSEMult[iPart1][iPart2]) {
        std::cout << "SEMult Histogramm missing from " << FolderName.Data()
                  << std::endl;
      }
      fSEMult[iPart1][iPart2]->Sumw2();
      if (iPart1 == 1 && iPart2 == 5) {
        fSE[iPart1][iPart2]->SetBinContent(1, 0);
        fSEMult[iPart1][iPart2]->SetBinContent(
            fSEMult[iPart1][iPart2]->GetXaxis()->FindBin(0.1),
            fSEMult[iPart1][iPart2]->GetYaxis()->FindBin(12.1), 0);
      }
      fME[iPart1][iPart2] = nullptr;
      hist1D = (TH1F*) PartList->FindObject(
          Form("MEDist_%s", FolderName.Data()));
      fME[iPart1][iPart2] = (TH1F*) hist1D->Clone(
          Form("%s_clone", hist1D->GetName()));
      if (!fME[iPart1][iPart2]) {
        std::cout << "ME Histogramm missing from " << FolderName.Data()
                  << std::endl;
      }
      fME[iPart1][iPart2]->Sumw2();

      fMEMult[iPart1][iPart2] = nullptr;
      hist2D = (TH2F*) PartList->FindObject(
          Form("MEMultDist_%s", FolderName.Data()));
      fMEMult[iPart1][iPart2] = (TH2F*) hist2D->Clone(
          Form("%s_clone", hist2D->GetName()));
      if (!fMEMult[iPart1][iPart2]) {
        std::cout << "ME Mult Histogramm missing from " << FolderName.Data()
                  << std::endl;
      }
      fMEMult[iPart1][iPart2]->Sumw2();
    }
  }
}

void ReadDreamFile::ReadkTHistos(const char* AnalysisFile, const char* prefix,
                                 const char* addon) {
  fSEkT = new TH2F**[fNPart1];
  fMEkT = new TH2F**[fNPart1];

  TFile* _file0 = TFile::Open(AnalysisFile, "READ");
  TDirectoryFile *dirResults = (TDirectoryFile*) (_file0->FindObjectAny(
      Form("%sResults%s", prefix, addon)));
  TList *Results;
  dirResults->GetObject(Form("%sResults%s", prefix, addon), Results);
  TList *PartList;
  for (int iPart1 = 0; iPart1 < fNPart1; ++iPart1) {
    fSEkT[iPart1] = new TH2F*[fNPart2];
    fMEkT[iPart1] = new TH2F*[fNPart2];

    for (int iPart2 = iPart1; iPart2 < fNPart2; ++iPart2) {
      TString FolderName = Form("Particle%i_Particle%i", iPart1, iPart2);
      PartList = (TList*) Results->FindObject(FolderName.Data());

      fSEkT[iPart1][iPart2] = nullptr;
      fSEkT[iPart1][iPart2] = (TH2F*) PartList->FindObject(
          Form("SEkTDist_%s", FolderName.Data()));
      if (!fSEkT[iPart1][iPart2]) {
        std::cout << "SEkT Histogramm missing from " << FolderName.Data()
                  << std::endl;
      }

      fMEkT[iPart1][iPart2] = nullptr;
      fMEkT[iPart1][iPart2] = (TH2F*) PartList->FindObject(
          Form("MEkTDist_%s", FolderName.Data()));
      if (!fMEkT[iPart1][iPart2]) {
        std::cout << "MEkT Histogramm missing from " << FolderName.Data()
                  << std::endl;
      }
    }
  }
  return;
}

void ReadDreamFile::ReadmTHistos(const char* AnalysisFile, const char* prefix,
                                 const char* addon) {
  fSEmT = new TH2F**[fNPart1];
  fMEmT = new TH2F**[fNPart1];

  TFile* _file0 = TFile::Open(AnalysisFile, "READ");
  TDirectoryFile *dirResults = (TDirectoryFile*) (_file0->FindObjectAny(
      Form("%sResults%s", prefix, addon)));
  TList *Results;
  dirResults->GetObject(Form("%sResults%s", prefix, addon), Results);
  TList *PartList;
  for (int iPart1 = 0; iPart1 < fNPart1; ++iPart1) {
    fSEmT[iPart1] = new TH2F*[fNPart2];
    fMEmT[iPart1] = new TH2F*[fNPart2];

    for (int iPart2 = iPart1; iPart2 < fNPart2; ++iPart2) {
      TString FolderName = Form("Particle%i_Particle%i", iPart1, iPart2);
      PartList = (TList*) Results->FindObject(FolderName.Data());

      fSEmT[iPart1][iPart2] = nullptr;
      fSEmT[iPart1][iPart2] = (TH2F*) PartList->FindObject(
          Form("SEmTDist_%s", FolderName.Data()));
      if (!fSEmT[iPart1][iPart2]) {
        std::cout << "SEmT Histogramm missing from " << FolderName.Data()
                  << std::endl;
      }

      fMEmT[iPart1][iPart2] = nullptr;
      fMEmT[iPart1][iPart2] = (TH2F*) PartList->FindObject(
          Form("MEmTDist_%s", FolderName.Data()));
      if (!fMEmT[iPart1][iPart2]) {
        std::cout << "MEmT Histogramm missing from " << FolderName.Data()
                  << std::endl;
      }
    }
  }
  return;
}

DreamDist* ReadDreamFile::GetPairDistributions(int iPart1, int iPart2,
                                               const char* name) {
  //user needs to ensure deletion
  if (iPart2 < iPart1) {
    std::cout << "Particle Combination does not exist \n";
    return nullptr;
  }
  DreamDist* pair = new DreamDist();
  pair->SetSEDist(fSE[iPart1][iPart2], name);
  pair->SetSEMultDist(fSEMult[iPart1][iPart2], name);
  pair->SetMEDist(fME[iPart1][iPart2], name);
  pair->SetMEMultDist(fMEMult[iPart1][iPart2], name);
  return pair;
}

DreamKayTee* ReadDreamFile::GetkTPairDistributions(int iPart1, int iPart2,
                                                   int iAPart1, int iAPart2) {
  //user needs to ensure deletion
  if (iPart2 < iPart1) {
    std::cout << "Particle Combination does not exist \n";
    return nullptr;
  }
  DreamKayTee* pair = new DreamKayTee();
  pair->SetSEkTDist(0, fSEkT[iPart1][iPart2]);
  pair->SetMEkTDist(0, fMEkT[iPart1][iPart2]);

  pair->SetSEkTDist(1, fSEkT[iPart1][iPart2]);
  pair->SetMEkTDist(1, fMEkT[iPart1][iPart2]);

  return pair;
}

DreamKayTee* ReadDreamFile::GetmTPairDistributions(int iPart1, int iPart2,
                                                   int iAPart1, int iAPart2) {
  //user needs to ensure deletion
  if (iPart2 < iPart1) {
    std::cout << "Particle Combination does not exist \n";
    return nullptr;
  }
  DreamKayTee* pair = new DreamKayTee();
  pair->SetSEmTDist(0, fSEmT[iPart1][iPart2]);
  pair->SetMEmTDist(0, fMEmT[iPart1][iPart2]);

  pair->SetSEmTDist(1, fSEmT[iPart1][iPart2]);
  pair->SetMEmTDist(1, fMEmT[iPart1][iPart2]);

  return pair;
}

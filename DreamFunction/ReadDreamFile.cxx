/*
 * ReadDreamFile.cxx
 *
 *  Created on: Aug 21, 2018
 *      Author: hohlweger
 */
#include "ReadDreamFile.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TLegend.h"
#include <iostream>
#include <iostream>
#include "stdlib.h"

ReadDreamFile::ReadDreamFile(int nPart1, int nPart2)
    : fQuiet(false),
      fNPart1(nPart1),
      fNPart2(nPart2),
      fSE(nullptr),
      fSEMult(nullptr),
      fSEkT(nullptr),
      fSEmT(nullptr),
      fSEmTProj(nullptr),
      fProjmT(nullptr),
      fMeanmT(nullptr),
      fSEdEtadPhimT(nullptr),
      fSEdEtadPhi(nullptr),
      fSEdEtadPhiAtRad(nullptr),
      fSEdEtadPhiAtRadSmallkStar(nullptr),
      fME(nullptr),
      fMEMult(nullptr),
      fMEkT(nullptr),
      fMEmT(nullptr),
      fMEdEtadPhimT(nullptr),
      fMEdEtadPhi(nullptr),
      fMEdEtadPhiAtRad(nullptr),
      fMEdEtadPhiAtRadSmallkStar(nullptr),
      fSECommon(nullptr),
      fSENonCommon(nullptr),
      fSEMultCommon(nullptr),
      fSEMultNonCommon(nullptr),
      fSEmTCommon(nullptr),
      fSEmTNonCommon(nullptr),
      fSEdEtadPhiCommon(nullptr),
      fSEdEtadPhiNonCommon(nullptr) {
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
      if (fSEmTProj && fSEmTProj[iPart1][iPart2])
        delete fSEmTProj[iPart1][iPart2];
      if (fProjmT && fProjmT[iPart1][iPart2])
        delete fProjmT[iPart1][iPart2];
      if (fMeanmT && fMeanmT[iPart1][iPart2])
        delete fMeanmT[iPart1][iPart2];
      if (fSEdEtadPhi && fSEdEtadPhi[iPart1][iPart2])
        delete fSEdEtadPhi[iPart1][iPart2];
      if (fME && fME[iPart1][iPart2])
        delete fME[iPart1][iPart2];
      if (fMEMult && fMEMult[iPart1][iPart2])
        delete fMEMult[iPart1][iPart2];
      if (fMEkT && fMEkT[iPart1][iPart2])
        delete fMEkT[iPart1][iPart2];
      if (fMEmT && fMEmT[iPart1][iPart2])
        delete fMEmT[iPart1][iPart2];
      if (fMEdEtadPhi && fMEdEtadPhi[iPart1][iPart2])
        delete fMEdEtadPhi[iPart1][iPart2];
      if (fSECommon && fSECommon[iPart1][iPart2])
        delete fSECommon[iPart1][iPart2];
      if (fSENonCommon && fSENonCommon[iPart1][iPart2])
        delete fSENonCommon[iPart1][iPart2];
      if (fSEMultCommon && fSEMultCommon[iPart1][iPart2])
        delete fSEMultCommon[iPart1][iPart2];
      if (fSEMultNonCommon && fSEMultNonCommon[iPart1][iPart2])
        delete fSEMultNonCommon[iPart1][iPart2];
      if (fSEdEtadPhiCommon && fSEdEtadPhiCommon[iPart1][iPart2])
        delete fSEdEtadPhiCommon[iPart1][iPart2];
      if (fSEdEtadPhiNonCommon && fSEdEtadPhiNonCommon[iPart1][iPart2])
        delete fSEdEtadPhiNonCommon[iPart1][iPart2];
    }
  }
}

void ReadDreamFile::SetAnalysisFile(const char* PathAnalysisFile,
                                    const char* Prefix, const char* Addon) {
  std::cout << "PathAnalysisFile: " << PathAnalysisFile << '\t' << " Prefix: "
            << Prefix << '\t' << " Addon: " << Addon << std::endl;
  TFile * _file0 = TFile::Open(PathAnalysisFile, "READ");
  TDirectoryFile *dirResults = (TDirectoryFile*) (_file0->FindObjectAny(
      TString::Format("%sResults%s", Prefix, Addon).Data()));
  if (!dirResults) {
    std::cout << "no dir results "
              << TString::Format("%sResults%s", Prefix, Addon).Data()
              << std::endl;
  }
  TList *Results = nullptr;
  dirResults->GetObject(TString::Format("%sResults%s", Prefix, Addon).Data(),
                        Results);
  if (!Results) {
    std::cout << "No results List for "
              << TString::Format("%sResults%s", Prefix, Addon).Data()
              << std::endl;
    dirResults->ls();
    return;
  }
  ExtractResults(Results);
  TIter next(Results);
  TObject *obj = nullptr;
  while (obj = next()) {
    TList *list = dynamic_cast<TList*>(obj);
    if (list)
      list->Delete();
  }
  Results->Delete();
  dirResults->Close();
  _file0->Close();
  delete _file0;
}

void ReadDreamFile::SetAnalysisFile(const char* PathAnalysisFile,
                                    const char* Path, const char* Prefix,
                                    const char* Addon) {
  TFile* _file0 = TFile::Open(PathAnalysisFile, "READ");
  TDirectoryFile *dirResults = (TDirectoryFile*) (_file0->FindObjectAny(
      Form("%sResults%s", Prefix, Addon)));
  TList *Results;
  dirResults->GetObject(Form("%sResults%s", Prefix, Addon), Results);
  auto listResults = (TList*) Results->FindObject(Path);
  ExtractResults(listResults);
  TIter next(listResults);
  TObject *obj = nullptr;
  while (obj = next()) {
    TList *list = dynamic_cast<TList*>(obj);
    if (list)
      list->Delete();
  }
  listResults->Delete();
  Results->Delete();
  dirResults->Close();
  _file0->Close();
  delete _file0;
}

void ReadDreamFile::SetAnalysisFileSample(const char* PathAnalysisFile,
                                    const char* Prefix, const char* Addon) {
  std::cout << "PathAnalysisFile: " << PathAnalysisFile << '\t' << " Prefix: "
            << Prefix << '\t' << " Addon: " << Addon << std::endl;
  TFile * _file0 = TFile::Open(PathAnalysisFile, "READ");
  TDirectoryFile *dirResults = (TDirectoryFile*) (_file0->FindObjectAny(
      TString::Format("%sResultsSample%s", Prefix, Addon).Data()));
  if (!dirResults) {
    std::cout << "no dir results "
              << TString::Format("%sResultsSample%s", Prefix, Addon).Data()
              << std::endl;
  }
  TList *Results = nullptr;
  dirResults->GetObject(TString::Format("%sResultsSample%s", Prefix, Addon).Data(),
                        Results);
  if (!Results) {
    std::cout << "No results List for "
              << TString::Format("%sResultsSample%s", Prefix, Addon).Data()
              << std::endl;
    dirResults->ls();
    return;
  }
  ExtractResults(Results);
  TIter next(Results);
  TObject *obj = nullptr;
  while (obj = next()) {
    TList *list = dynamic_cast<TList*>(obj);
    if (list)
      list->Delete();
  }
  Results->Delete();
  dirResults->Close();
  _file0->Close();
  delete _file0;
}

void ReadDreamFile::SetAnalysisFileAncestors(const char* PathAnalysisFile,
                                    const char* Prefix, const char* Addon) {
  std::cout << "PathAnalysisFile: " << PathAnalysisFile << '\t' << " Prefix: "
            << Prefix << '\t' << " Addon: " << Addon << std::endl;
  TFile * _file0 = TFile::Open(PathAnalysisFile, "READ");
  TDirectoryFile *dirResults = (TDirectoryFile*) (_file0->FindObjectAny(
      TString::Format("%sResults%s", Prefix, Addon).Data()));
  if (!dirResults) {
    std::cout << "no dir results "
              << TString::Format("%sResults%s", Prefix, Addon).Data()
              << std::endl;
  }
  TList *Results = nullptr;
  dirResults->GetObject(TString::Format("%sResults%s", Prefix, Addon).Data(),
                        Results);
  if (!Results) {
    std::cout << "No results List for "
              << TString::Format("%sResults%s", Prefix, Addon).Data()
              << std::endl;
    dirResults->ls();
    return;
  }
  ExtractResultsAncestors(Results);
  TIter next(Results);
  TObject *obj = nullptr;
  while (obj = next()) {
    TList *list = dynamic_cast<TList*>(obj);
    if (list)
      list->Delete();
  }
  Results->Delete();
  dirResults->Close();
  _file0->Close();
  delete _file0;
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
      if (!hist1D) {
        if (!fQuiet)
          std::cout << "SE Histogramm missing from " << FolderName.Data()
                    << std::endl;
      } else {
        fSE[iPart1][iPart2] = (TH1F*) hist1D->Clone(
            Form("%s_clone", hist1D->GetName()));
        fSE[iPart1][iPart2]->Sumw2();
      }
      fSEMult[iPart1][iPart2] = nullptr;
      auto hist2D = (TH2F*) PartList->FindObject(
          Form("SEMultDist_%s", FolderName.Data()));
      if (!hist2D) {
        if (!fQuiet)
          std::cout << "SEMult Histogramm missing from " << FolderName.Data()
                    << std::endl;
      } else {
        fSEMult[iPart1][iPart2] = (TH2F*) hist2D->Clone(
            Form("%s_clone", hist2D->GetName()));

        fSEMult[iPart1][iPart2]->Sumw2();
      }
      //instead start the fixed shifted binning at 8!
//      if (iPart1 == 1 && iPart2 == 5) {
//        fSE[iPart1][iPart2]->SetBinContent(1, 0);
//        fSEMult[iPart1][iPart2]->SetBinContent(
//            fSEMult[iPart1][iPart2]->GetXaxis()->FindBin(0.1),
//            fSEMult[iPart1][iPart2]->GetYaxis()->FindBin(12.1), 0);
//      }
      fME[iPart1][iPart2] = nullptr;
      hist1D = (TH1F*) PartList->FindObject(
          Form("MEDist_%s", FolderName.Data()));
      if (!hist1D) {
        if (!fQuiet)
          std::cout << "ME Histogramm missing from " << FolderName.Data()
                    << std::endl;
      } else {
        fME[iPart1][iPart2] = (TH1F*) hist1D->Clone(
            Form("%s_clone", hist1D->GetName()));

        fME[iPart1][iPart2]->Sumw2();
      }
      fMEMult[iPart1][iPart2] = nullptr;
      hist2D = (TH2F*) PartList->FindObject(
          Form("MEMultDist_%s", FolderName.Data()));
      if (!hist2D) {
        if (!fQuiet)
          std::cout << "ME Mult Histogramm missing from " << FolderName.Data()
                    << std::endl;
      } else {
        fMEMult[iPart1][iPart2] = (TH2F*) hist2D->Clone(
            Form("%s_clone", hist2D->GetName()));

        fMEMult[iPart1][iPart2]->Sumw2();
      }
    }
  }
}

void ReadDreamFile::ExtractResultsAncestors(const TList *Results) {
  TList *PartList;

  fSECommon = new TH1F**[fNPart1];
  fSEMultCommon = new TH2F**[fNPart1];

  fSENonCommon = new TH1F**[fNPart1];
  fSEMultNonCommon = new TH2F**[fNPart1];

  fME = new TH1F**[fNPart1];
  fMEMult = new TH2F**[fNPart1];

  for (int iPart1 = 0; iPart1 < fNPart1; ++iPart1) {
    fSECommon[iPart1] = new TH1F*[fNPart2];
    fSEMultCommon[iPart1] = new TH2F*[fNPart2];

    fSENonCommon[iPart1] = new TH1F*[fNPart2];
    fSEMultNonCommon[iPart1] = new TH2F*[fNPart2];

    fME[iPart1] = new TH1F*[fNPart2];
    fMEMult[iPart1] = new TH2F*[fNPart2];

    for (int iPart2 = iPart1; iPart2 < fNPart2; ++iPart2) {
      TString FolderName = Form("Particle%i_Particle%i", iPart1, iPart2);
      PartList = (TList*) Results->FindObject(FolderName.Data());
      fSECommon[iPart1][iPart2] = nullptr;
      fSENonCommon[iPart1][iPart2] = nullptr;
      auto hist1DCommon = (TH1F*) PartList->FindObject(
          Form("SEDistCommon_%s", FolderName.Data()));
      auto hist1DNonCommon = (TH1F*) PartList->FindObject(
          Form("SEDistNonCommon_%s", FolderName.Data()));
      if (!hist1DCommon | !hist1DNonCommon) {
        if (!fQuiet)
          std::cout << "SE Histogramm missing from " << FolderName.Data()
                    << std::endl;
      } else {
        fSECommon[iPart1][iPart2] = (TH1F*) hist1DCommon->Clone(
            Form("%s_clone", hist1DCommon->GetName()));
        fSECommon[iPart1][iPart2]->Sumw2();
        fSENonCommon[iPart1][iPart2] = (TH1F*) hist1DNonCommon->Clone(
            Form("%s_clone", hist1DNonCommon->GetName()));
        fSENonCommon[iPart1][iPart2]->Sumw2();
      }
      fSEMultCommon[iPart1][iPart2] = nullptr;
      fSEMultNonCommon[iPart1][iPart2] = nullptr;

      auto hist2DCommon = (TH2F*) PartList->FindObject(
          Form("SEMultDistCommon_%s", FolderName.Data()));
      auto hist2DNonCommon = (TH2F*) PartList->FindObject(
          Form("SEMultDistNonCommon_%s", FolderName.Data()));
      if (!hist2DCommon | !hist2DNonCommon) {
        if (!fQuiet)
          std::cout << "SEMult Histogramm missing from " << FolderName.Data()
                    << std::endl;
      } else {
        fSEMultCommon[iPart1][iPart2] = (TH2F*) hist2DCommon->Clone(
            Form("%s_clone", hist2DCommon->GetName()));
        fSEMultCommon[iPart1][iPart2]->Sumw2();
        fSEMultNonCommon[iPart1][iPart2] = (TH2F*) hist2DNonCommon->Clone(
            Form("%s_clone", hist2DNonCommon->GetName()));
        fSEMultNonCommon[iPart1][iPart2]->Sumw2();
      }
      //instead start the fixed shifted binning at 8!
//      if (iPart1 == 1 && iPart2 == 5) {
//        fSE[iPart1][iPart2]->SetBinContent(1, 0);
//        fSEMult[iPart1][iPart2]->SetBinContent(
//            fSEMult[iPart1][iPart2]->GetXaxis()->FindBin(0.1),
//            fSEMult[iPart1][iPart2]->GetYaxis()->FindBin(12.1), 0);
//      }
      fME[iPart1][iPart2] = nullptr;
      auto hist1D = (TH1F*) PartList->FindObject(
          Form("MEDist_%s", FolderName.Data()));
      if (!hist1D) {
        if (!fQuiet)
          std::cout << "ME Histogramm missing from " << FolderName.Data()
                    << std::endl;
      } else {
        fME[iPart1][iPart2] = (TH1F*) hist1D->Clone(
            Form("%s_clone", hist1D->GetName()));
        fME[iPart1][iPart2]->Sumw2();
      }

      fMEMult[iPart1][iPart2] = nullptr;
      auto hist2D = (TH2F*) PartList->FindObject(
          Form("MEMultDist_%s", FolderName.Data()));
      if (!hist2D) {
        if (!fQuiet)
          std::cout << "ME Mult Histogramm missing from " << FolderName.Data()
                    << std::endl;
      } else {
        fMEMult[iPart1][iPart2] = (TH2F*) hist2D->Clone(
            Form("%s_clone", hist2D->GetName()));
        fMEMult[iPart1][iPart2]->Sumw2();
    }
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
        if (!fQuiet)
          std::cout << "SEmT Histogramm missing from " << FolderName.Data()
                    << std::endl;
      }
      fMEmT[iPart1][iPart2] = nullptr;
      fMEmT[iPart1][iPart2] = (TH2F*) PartList->FindObject(
          Form("MEmTDist_%s", FolderName.Data()));
      if (!fMEmT[iPart1][iPart2]) {
        if (!fQuiet)
          std::cout << "MEmT Histogramm missing from " << FolderName.Data()
                    << std::endl;
      }
    }
  }
  return;
}

void ReadDreamFile::ReadAndProjectmTHistos(const char* AnalysisFile, const char* prefix,
                                 const char* addon, double kcut) {

  double kcutdummy = kcut/1000.;//transform in GeV
  double binwidth;

  fSEmT = new TH2F**[fNPart1];
  fMEmT = new TH2F**[fNPart1];
  fSEmTProj = new TH1F** [fNPart1];


  TFile* _file0 = TFile::Open(AnalysisFile, "READ");
  TDirectoryFile *dirResults = (TDirectoryFile*) (_file0->FindObjectAny(
      Form("%sResults%s", prefix, addon)));
  TList *Results;
  dirResults->GetObject(Form("%sResults%s", prefix, addon), Results);
  TList *PartList;

  for (int iPart1 = 0; iPart1 < fNPart1; ++iPart1) {
    fSEmT[iPart1] = new TH2F*[fNPart2];
    fMEmT[iPart1] = new TH2F*[fNPart2];
    fSEmTProj[iPart1] = new TH1F*[fNPart2];


    for (int iPart2 = iPart1; iPart2 < fNPart2; ++iPart2) {

      TString FolderName = Form("Particle%i_Particle%i", iPart1, iPart2);
      PartList = (TList*) Results->FindObject(FolderName.Data());

      fSEmT[iPart1][iPart2] = nullptr;
      fSEmT[iPart1][iPart2] = (TH2F*) PartList->FindObject(
          Form("SEmTDist_%s", FolderName.Data()));
	  fSEmTProj[iPart1][iPart2] = nullptr;

	  for(int ikstar = 1; ikstar < fSEmT[iPart1][iPart2]->GetNbinsX()+1; ikstar++)
      {
        binwidth = fSEmT[iPart1][iPart2]->GetXaxis()->GetBinUpEdge(ikstar)-
            fSEmT[iPart1][iPart2]->GetXaxis()->GetBinLowEdge(ikstar);

        fSEmTProj[iPart1][iPart2] = (TH1F*) fSEmT[iPart1][iPart2]->ProjectionY(
                TString::Format("fSEmTProj_part%i_part%i_%i",iPart1,iPart2,ikstar),ikstar,ikstar);

          if(fSEmT[iPart1][iPart2]->GetXaxis()->GetBinCenter(ikstar)>kcutdummy){
            break;
        }

        if (!fSEmTProj[iPart1][iPart2]) {
          if (!fQuiet)
            std::cout << "fSEmTProj Histogramm missing from " << FolderName.Data()
                      << std::endl;
        }
      }
      if (!fSEmT[iPart1][iPart2]) {
        if (!fQuiet)
          std::cout << "SEmT Histogramm missing from " << FolderName.Data()
                    << std::endl;
      }

    }
  }
  return;
}

void ReadDreamFile::ReadmTHistosAncestors(const char* AnalysisFile, const char* prefix,
                                 const char* addon) {
  fSEmTCommon = new TH2F**[fNPart1];
  fSEmTNonCommon = new TH2F**[fNPart1];

  fMEmT = new TH2F**[fNPart1];

  TFile* _file0 = TFile::Open(AnalysisFile, "READ");
  TDirectoryFile *dirResults = (TDirectoryFile*) (_file0->FindObjectAny(
      Form("%sResults%s", prefix, addon)));
  TList *Results;
  dirResults->GetObject(Form("%sResults%s", prefix, addon), Results);
  TList *PartList;
  for (int iPart1 = 0; iPart1 < fNPart1; ++iPart1) {
    fSEmTCommon[iPart1] = new TH2F*[fNPart2];
    fSEmTNonCommon[iPart1] = new TH2F*[fNPart2];

    fMEmT[iPart1] = new TH2F*[fNPart2];

    for (int iPart2 = iPart1; iPart2 < fNPart2; ++iPart2) {
      TString FolderName = Form("Particle%i_Particle%i", iPart1, iPart2);
      PartList = (TList*) Results->FindObject(FolderName.Data());

      fSEmTCommon[iPart1][iPart2] = nullptr;
      fSEmTNonCommon[iPart1][iPart2] = nullptr;

      fSEmTCommon[iPart1][iPart2] = (TH2F*) PartList->FindObject(
          Form("SEmTDistCommon_%s", FolderName.Data()));
      fSEmTNonCommon[iPart1][iPart2] = (TH2F*) PartList->FindObject(
          Form("SEmTDistNonCommon_%s", FolderName.Data()));
      if (!fSEmTCommon[iPart1][iPart2] | !fSEmTNonCommon[iPart1][iPart2]) {
        if (!fQuiet)
          std::cout << "SEmT Histogramm missing from " << FolderName.Data()
                    << std::endl;
      }
      fMEmT[iPart1][iPart2] = nullptr;
      fMEmT[iPart1][iPart2] = (TH2F*) PartList->FindObject(
          Form("MEmTDist_%s", FolderName.Data()));
      if (!fMEmT[iPart1][iPart2]) {
        if (!fQuiet)
          std::cout << "MEmT Histogramm missing from " << FolderName.Data()
                    << std::endl;
      }
    }
  }
  return;
}


void ReadDreamFile::ReadAndProjectkTHistos(const char* AnalysisFile, const char* prefix,
                                 const char* addon, double kcut) {

  double kcutdummy = kcut/1000.;//transform in GeV
  double binwidth;

  fSEmT = new TH2F**[fNPart1];
  fMEmT = new TH2F**[fNPart1];
  fSEmTProj = new TH1F** [fNPart1];


  TFile* _file0 = TFile::Open(AnalysisFile, "READ");
  TDirectoryFile *dirResults = (TDirectoryFile*) (_file0->FindObjectAny(
      Form("%sResults%s", prefix, addon)));
  TList *Results;
  dirResults->GetObject(Form("%sResults%s", prefix, addon), Results);
  TList *PartList;

  for (int iPart1 = 0; iPart1 < fNPart1; ++iPart1) {
    fSEmT[iPart1] = new TH2F*[fNPart2];
    fMEmT[iPart1] = new TH2F*[fNPart2];
    fSEmTProj[iPart1] = new TH1F*[fNPart2];


    for (int iPart2 = iPart1; iPart2 < fNPart2; ++iPart2) {

      TString FolderName = Form("Particle%i_Particle%i", iPart1, iPart2);
      PartList = (TList*) Results->FindObject(FolderName.Data());

      fSEmT[iPart1][iPart2] = nullptr;
      fSEmT[iPart1][iPart2] = (TH2F*) PartList->FindObject(
          Form("SEkTDist_%s", FolderName.Data()));
	  fSEmTProj[iPart1][iPart2] = nullptr;

	  for(int ikstar = 1; ikstar < fSEmT[iPart1][iPart2]->GetNbinsX()+1; ikstar++)
      {
      binwidth = fSEmT[iPart1][iPart2]->GetXaxis()->GetBinUpEdge(ikstar)-
    		  fSEmT[iPart1][iPart2]->GetXaxis()->GetBinLowEdge(ikstar);

      fSEmTProj[iPart1][iPart2] = (TH1F*) fSEmT[iPart1][iPart2]->ProjectionY(
    		  	  TString::Format("fSEkTProj_part%i_part%i_%i",iPart1,iPart2,ikstar),ikstar,ikstar);

      if(fSEmT[iPart1][iPart2]->GetXaxis()->GetBinCenter(ikstar)>kcutdummy){
    	  break;
      }

      if (!fSEmTProj[iPart1][iPart2]) {
        if (!fQuiet)
          std::cout << "fSEkTProj Histogramm missing from " << FolderName.Data()
                    << std::endl;
      }
      }

//====================================================
      if (!fSEmT[iPart1][iPart2]) {
        if (!fQuiet)
          std::cout << "SEkT Histogramm missing from " << FolderName.Data()
                    << std::endl;
      }

    }
  }
  return;
}


void ReadDreamFile::ReaddEtadPhiAtRadHists(const unsigned int nMaxMix,
                                           const char* AnalysisFile,
                                           const char* prefix,
                                           const char* Addon) {
  fSEdEtadPhiAtRad = new TH2F****[fNPart1];
  fSEdEtadPhiAtRadSmallkStar = new TH2F****[fNPart1];
  fMEdEtadPhiAtRad = new TH2F****[fNPart1];
  fMEdEtadPhiAtRadSmallkStar = new TH2F****[fNPart1];
  TFile* _file0 = TFile::Open(AnalysisFile, "READ");
  TDirectoryFile *dirResults = (TDirectoryFile*) (_file0->FindObjectAny(
      Form("%sResultQA%s", prefix, Addon)));
  TList *ResultsQA;
  dirResults->GetObject(Form("%sResultQA%s", prefix, Addon), ResultsQA);
  TList *PartList;

  for (int iPart1 = 0; iPart1 < fNPart1; ++iPart1) {
    fSEdEtadPhiAtRad[iPart1] = new TH2F***[fNPart2];
    fSEdEtadPhiAtRadSmallkStar[iPart1] = new TH2F***[fNPart2];
    fMEdEtadPhiAtRad[iPart1] = new TH2F***[fNPart2];
    fMEdEtadPhiAtRadSmallkStar[iPart1] = new TH2F***[fNPart2];

    for (int iPart2 = iPart1; iPart2 < fNPart2; ++iPart2) {
      TString FolderName = Form("QA_Particle%i_Particle%i", iPart1, iPart2);
      PartList = (TList*) ResultsQA->FindObject(FolderName.Data());
      fSEdEtadPhiAtRad[iPart1][iPart2] = new TH2F**[9];
      fSEdEtadPhiAtRadSmallkStar[iPart1][iPart2] = new TH2F**[9];
      fMEdEtadPhiAtRad[iPart1][iPart2] = new TH2F**[9];
      fMEdEtadPhiAtRadSmallkStar[iPart1][iPart2] = new TH2F**[9];
      for (int iRad = 0; iRad < 9; ++iRad) {
        fSEdEtadPhiAtRad[iPart1][iPart2][iRad] = new TH2F*[nMaxMix];
        fSEdEtadPhiAtRadSmallkStar[iPart1][iPart2][iRad] = new TH2F*[nMaxMix];
        fMEdEtadPhiAtRad[iPart1][iPart2][iRad] = new TH2F*[nMaxMix];
        fMEdEtadPhiAtRadSmallkStar[iPart1][iPart2][iRad] = new TH2F*[nMaxMix];
      }
      TIter next(PartList);
      TObject *obj = nullptr;
      while (obj = next()) {
        TString objName = obj->GetName();
        if (objName.Contains("Rad")) {
          TString beAChar = objName[objName.First('_') + 1];
          int iRad = atoi(beAChar.Data());
          TString beAnotherChar = objName[objName.First("x") + 1];
          int iMix = atoi(beAnotherChar.Data());
          if (objName.Contains("SE")) {
            if (objName.Contains("smallK")) {
              fSEdEtadPhiAtRadSmallkStar[iPart1][iPart2][iRad][iMix] =
                  (TH2F*) obj;
            } else {
              fSEdEtadPhiAtRad[iPart1][iPart2][iRad][iMix] = (TH2F*) obj;
            }
          } else if (objName.Contains("ME")) {
            if (objName.Contains("smallK")) {
              fMEdEtadPhiAtRadSmallkStar[iPart1][iPart2][iRad][iMix] =
                  (TH2F*) obj;
            } else {
              fMEdEtadPhiAtRad[iPart1][iPart2][iRad][iMix] = (TH2F*) obj;
            }
          }
        }
      }
    }
  }
  return;
}

void ReadDreamFile::ReaddEtadPhiHists(const unsigned int NBinsmT,
                                      const char* AnalysisFile,
                                      const char* prefix, const char* Addon) {
  fSEdEtadPhi = new TH2F**[fNPart1];
  fMEdEtadPhi = new TH2F**[fNPart1];
  if (NBinsmT > 0) {
    fSEdEtadPhimT = new TH2F***[fNPart1];
    fMEdEtadPhimT = new TH2F***[fNPart1];
  }
  TFile* _file0 = TFile::Open(AnalysisFile, "READ");
  TDirectoryFile *dirResults = (TDirectoryFile*) (_file0->FindObjectAny(
      Form("%sResults%s", prefix, Addon)));
  TList *Results;
  dirResults->GetObject(Form("%sResults%s", prefix, Addon), Results);
  TList *PartList;

  for (int iPart1 = 0; iPart1 < fNPart1; ++iPart1) {
    fSEdEtadPhi[iPart1] = new TH2F*[fNPart2];
    fMEdEtadPhi[iPart1] = new TH2F*[fNPart2];

    if (NBinsmT > 0) {
      fSEdEtadPhimT[iPart1] = new TH2F**[fNPart2];
      fMEdEtadPhimT[iPart1] = new TH2F**[fNPart2];
    }
    for (int iPart2 = iPart1; iPart2 < fNPart2; ++iPart2) {
      TString FolderName = Form("Particle%i_Particle%i", iPart1, iPart2);

      PartList = (TList*) Results->FindObject(FolderName.Data());
      if (NBinsmT > 0) {
        fSEdEtadPhimT[iPart1][iPart2] = new TH2F*[NBinsmT];
        fMEdEtadPhimT[iPart1][iPart2] = new TH2F*[NBinsmT];
      }
      int iSEmTCounter = 0;
      int iMEmTCounter = 0;
      TIter next(PartList);
      TObject *obj = nullptr;
      while (obj = next()) {
        TString objName = obj->GetName();
        if (NBinsmT > 0 && objName.Contains("imT")) {
          if (objName.Contains("SE")) {
            fSEdEtadPhimT[iPart1][iPart2][iSEmTCounter] = (TH2F*) obj;
//            std::cout << "SE: " << objName.Data() << std::endl;
            if (!fSEdEtadPhimT[iPart1][iPart2][iSEmTCounter]) {
              std::cout << objName.Data() << " failed to deliver an object \n";
            }
            iSEmTCounter++;
          } else if (objName.Contains("ME")) {
            fMEdEtadPhimT[iPart1][iPart2][iMEmTCounter] = (TH2F*) obj;
//            std::cout << "ME: " << objName.Data() << std::endl;
            if (!fMEdEtadPhimT[iPart1][iPart2][iMEmTCounter]) {
              std::cout << objName.Data() << " failed to deliver an object \n";
            }
            iMEmTCounter++;
          } else {
            std::cout << objName.Data()
                      << " contains imT but neither SE nor ME \n";
          }
        } else if (NBinsmT == 0 && objName.Contains("dPhidEtaDist")) {
          if (objName.Contains("SE")) {
            fSEdEtadPhi[iPart1][iPart2] = (TH2F*) obj;
          } else if (objName.Contains("ME")) {
            fMEdEtadPhi[iPart1][iPart2] = (TH2F*) obj;
          }
        }
      }

      if (NBinsmT > 0) {
        if (iSEmTCounter > 0 && iMEmTCounter) {
          std::cout << "Pair " << iPart1 << " & " << iPart2
                    << " Manually creating the summed SE and ME Hist ... \n";
          fSEdEtadPhi[iPart1][iPart2] = (TH2F*) fSEdEtadPhimT[iPart1][iPart2][0]
              ->Clone(
              Form("SEdPhidEtaDist_Particle%d_Particle%d", iPart1, iPart2));
          fMEdEtadPhi[iPart1][iPart2] = (TH2F*) fMEdEtadPhimT[iPart1][iPart2][0]
              ->Clone(
              Form("MEdPhidEtaDist_Particle%d_Particle%d", iPart1, iPart2));
          for (int imT = 1; imT < iSEmTCounter; ++imT) {
            fSEdEtadPhi[iPart1][iPart2]->Add(
                fSEdEtadPhimT[iPart1][iPart2][imT]);
            fMEdEtadPhi[iPart1][iPart2]->Add(
                fMEdEtadPhimT[iPart1][iPart2][imT]);
          }
        } else {
          std::cout << "Pair " << iPart1 << " & " << iPart2 << " not summed \n";
        }
      }
    }
  }
}

void ReadDreamFile::ReaddEtadPhiHistsAncestors(const unsigned int NBinsmT,
                                      const char* AnalysisFile,
                                      const char* prefix, const char* Addon) {
  fSEdEtadPhiCommon = new TH2F**[fNPart1];
  fSEdEtadPhiNonCommon = new TH2F**[fNPart1];
  fMEdEtadPhi = new TH2F**[fNPart1];
  TFile* _file0 = TFile::Open(AnalysisFile, "READ");
  TDirectoryFile *dirResults = (TDirectoryFile*) (_file0->FindObjectAny(
      Form("%sResults%s", prefix, Addon)));
  TList *Results;
  dirResults->GetObject(Form("%sResults%s", prefix, Addon), Results);
  TList *PartList;
  for (int iPart1 = 0; iPart1 < fNPart1; ++iPart1) {
    fSEdEtadPhiCommon[iPart1] = new TH2F*[fNPart2];
    fSEdEtadPhiNonCommon[iPart1] = new TH2F*[fNPart2];
    fMEdEtadPhi[iPart1] = new TH2F*[fNPart2];

    for (int iPart2 = iPart1; iPart2 < fNPart2; ++iPart2) {
      TString FolderName = Form("Particle%i_Particle%i", iPart1, iPart2);

      PartList = (TList*) Results->FindObject(FolderName.Data());

      TIter next(PartList);
      TObject *obj = nullptr;
      while (obj = next()) {
        TString objName = obj->GetName();
      if (NBinsmT == 0 && objName.Contains("dPhidEtaDist")) {
          if (objName.Contains("SE")) {
            if (objName.Contains("Common")) {
              if (objName.Contains("Non")) fSEdEtadPhiNonCommon[iPart1][iPart2] = (TH2F*) obj;
              else fSEdEtadPhiCommon[iPart1][iPart2] = (TH2F*) obj;
          }
          } else if (objName.Contains("ME")) {
            fMEdEtadPhi[iPart1][iPart2] = (TH2F*) obj;
          }
        }
      }
    }
  }
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

DreamDist* ReadDreamFile::GetPairDistributionsCommon(int iPart1, int iPart2,
                                               const char* name) {
//user needs to ensure deletion
  if (iPart2 < iPart1) {
    std::cout << "Particle Combination does not exist \n";
    return nullptr;
  }
  DreamDist* pair = new DreamDist();
  pair->SetSEDist(fSECommon[iPart1][iPart2], name);
  pair->SetSEMultDist(fSEMultCommon[iPart1][iPart2], name);
  pair->SetMEDist(fME[iPart1][iPart2], name);
  pair->SetMEMultDist(fMEMult[iPart1][iPart2], name);

  return pair;
}

DreamDist* ReadDreamFile::GetPairDistributionsNonCommon(int iPart1, int iPart2,
                                               const char* name) {
//user needs to ensure deletion
  if (iPart2 < iPart1) {
    std::cout << "Particle Combination does not exist \n";
    return nullptr;
  }
  DreamDist* pair = new DreamDist();
  pair->SetSEDist(fSENonCommon[iPart1][iPart2], name);
  pair->SetSEMultDist(fSEMultNonCommon[iPart1][iPart2], name);
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

  pair->SetSEkTDist(1, fSEkT[iAPart1][iAPart2]);
  pair->SetMEkTDist(1, fMEkT[iAPart1][iAPart2]);

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

  pair->SetSEmTDist(1, fSEmT[iAPart1][iAPart2]);
  pair->SetMEmTDist(1, fMEmT[iAPart1][iAPart2]);

  return pair;
}

DreamKayTee* ReadDreamFile::GetmTPairDistributionsCommon(int iPart1, int iPart2,
                                                   int iAPart1, int iAPart2) {
//user needs to ensure deletion
  if (iPart2 < iPart1) {
    std::cout << "Particle Combination does not exist \n";
    return nullptr;
  }
  DreamKayTee* pair = new DreamKayTee();
  pair->SetSEmTDist(0, fSEmTCommon[iPart1][iPart2]);
  pair->SetMEmTDist(0, fMEmT[iPart1][iPart2]);

  pair->SetSEmTDist(1, fSEmTCommon[iAPart1][iAPart2]);
  pair->SetMEmTDist(1, fMEmT[iAPart1][iAPart2]);

  return pair;
}

DreamKayTee* ReadDreamFile::GetmTPairDistributionsNonCommon(int iPart1, int iPart2,
                                                   int iAPart1, int iAPart2) {
//user needs to ensure deletion
  if (iPart2 < iPart1) {
    std::cout << "Particle Combination does not exist \n";
    return nullptr;
  }
  DreamKayTee* pair = new DreamKayTee();
  pair->SetSEmTDist(0, fSEmTNonCommon[iPart1][iPart2]);
  pair->SetMEmTDist(0, fMEmT[iPart1][iPart2]);

  pair->SetSEmTDist(1, fSEmTNonCommon[iAPart1][iAPart2]);
  pair->SetMEmTDist(1, fMEmT[iAPart1][iAPart2]);

  return pair;
}

DreamKayTee* ReadDreamFile::GetmTPairDistributionsCommon(int iPart1, int iPart2) {
//user needs to ensure deletion
  if (iPart2 < iPart1) {
    std::cout << "Particle Combination does not exist \n";
    return nullptr;
  }
  DreamKayTee* pair = new DreamKayTee();
  pair->SetSEmTDist(0, fSEmTCommon[iPart1][iPart2]);
  pair->SetMEmTDist(0, fMEmT[iPart1][iPart2]);

  return pair;
}

DreamKayTee* ReadDreamFile::GetmTPairDistributionsNonCommon(int iPart1, int iPart2) {
//user needs to ensure deletion
  if (iPart2 < iPart1) {
    std::cout << "Particle Combination does not exist \n";
    return nullptr;
  }
  DreamKayTee* pair = new DreamKayTee();
  pair->SetSEmTDist(0, fSEmTNonCommon[iPart1][iPart2]);
  pair->SetMEmTDist(0, fMEmT[iPart1][iPart2]);

  return pair;
}

DreamKayTee* ReadDreamFile::GetmTPairDistributionsBBar(int iPart1, int iPart2) {
//user needs to ensure deletion
  if (iPart2 < iPart1) {
    std::cout << "Particle Combination does not exist \n";
    return nullptr;
  }
  DreamKayTee* pair = new DreamKayTee();
  pair->SetSEmTDist(0, fSEmT[iPart1][iPart2]);
  pair->SetMEmTDist(0, fMEmT[iPart1][iPart2]);

  return pair;
}



DreamdEtadPhi* ReadDreamFile::GetdEtadPhiDistribution(int iPart1, int iPart2,
                                                      int iAPart1, int iAPart2,
                                                      int imT) {
  if (iPart2 < iPart1) {
    std::cout << "Particle Combination does not exist \n";
    return nullptr;
  }
  DreamdEtadPhi* outDist = new DreamdEtadPhi();
  if (imT < 0) {
    outDist->SetSEDistribution(fSEdEtadPhi[iPart1][iPart2], "");
    outDist->AddSEDistribution(fSEdEtadPhi[iAPart1][iAPart2]);
    outDist->SetMEDistribution(fMEdEtadPhi[iPart1][iPart2], "");
    outDist->AddMEDistribution(fMEdEtadPhi[iAPart1][iAPart2]);
  } else {
    outDist->SetSEDistribution(fSEdEtadPhimT[iPart1][iPart2][imT], "");
    outDist->AddSEDistribution(fSEdEtadPhimT[iAPart1][iAPart2][imT]);
    outDist->SetMEDistribution(fMEdEtadPhimT[iPart1][iPart2][imT], "");
    outDist->AddMEDistribution(fMEdEtadPhimT[iAPart1][iAPart2][imT]);
  }
  return outDist;
}

DreamdEtadPhi* ReadDreamFile::GetdEtadPhiDistributionCommon(int iPart1, int iPart2,
                                                      int iAPart1, int iAPart2,
                                                      int imT) {
  if (iPart2 < iPart1) {
    std::cout << "Particle Combination does not exist \n";
    return nullptr;
  }
  DreamdEtadPhi* outDistCommon = new DreamdEtadPhi();
  if (imT < 0) {
    outDistCommon->SetSEDistribution(fSEdEtadPhiCommon[iPart1][iPart2], "");
    outDistCommon->AddSEDistribution(fSEdEtadPhiCommon[iAPart1][iAPart2]);
    outDistCommon->SetMEDistribution(fMEdEtadPhi[iPart1][iPart2], "");
    outDistCommon->AddMEDistribution(fMEdEtadPhi[iAPart1][iAPart2]);
  }
  return outDistCommon;
}

DreamdEtadPhi* ReadDreamFile::GetdEtadPhiDistributionNonCommon(int iPart1, int iPart2,
                                                      int iAPart1, int iAPart2,
                                                      int imT) {
  if (iPart2 < iPart1) {
    std::cout << "Particle Combination does not exist \n";
    return nullptr;
  }
  DreamdEtadPhi* outDistNonCommon = new DreamdEtadPhi();

  if (imT < 0) {
    outDistNonCommon->SetSEDistribution(fSEdEtadPhiNonCommon[iPart1][iPart2], "A");
    outDistNonCommon->AddSEDistribution(fSEdEtadPhiNonCommon[iAPart1][iAPart2]);
    outDistNonCommon->SetMEDistribution(fMEdEtadPhi[iPart1][iPart2], "");
    outDistNonCommon->AddMEDistribution(fMEdEtadPhi[iAPart1][iAPart2]);
  }
  return outDistNonCommon;
}

DreamdEtadPhi* ReadDreamFile::GetdEtadPhiDistributionSingle(int iPart1,
                                                            int iPart2,
                                                            int imT) {
  if (iPart2 < iPart1) {
    std::cout << "Particle Combination does not exist \n";
    return nullptr;
  }
  DreamdEtadPhi* outDist = new DreamdEtadPhi();
  if (imT < 0) {
    outDist->SetSEDistribution(fSEdEtadPhi[iPart1][iPart2], "");
    outDist->SetMEDistribution(fMEdEtadPhi[iPart1][iPart2], "");
  } else {
    outDist->SetSEDistribution(fSEdEtadPhimT[iPart1][iPart2][imT], "");
    outDist->SetMEDistribution(fMEdEtadPhimT[iPart1][iPart2][imT], "");
  }
  return outDist;
}

DreamdEtadPhi* ReadDreamFile::GetdEtadPhiDistributionSingleCommon(int iPart1,
                                                            int iPart2,
                                                            int imT) {
  if (iPart2 < iPart1) {
    std::cout << "Particle Combination does not exist \n";
    return nullptr;
  }
  DreamdEtadPhi* outDist = new DreamdEtadPhi();
  if (imT < 0) {
    outDist->SetSEDistribution(fSEdEtadPhiCommon[iPart1][iPart2], "");
    outDist->SetMEDistribution(fMEdEtadPhi[iPart1][iPart2], "");
  }
  return outDist;
}

DreamdEtadPhi* ReadDreamFile::GetdEtadPhiDistributionSingleNonCommon(int iPart1,
                                                            int iPart2,
                                                            int imT) {
  if (iPart2 < iPart1) {
    std::cout << "Particle Combination does not exist \n";
    return nullptr;
  }
  DreamdEtadPhi* outDist = new DreamdEtadPhi();
  if (imT < 0) {
    outDist->SetSEDistribution(fSEdEtadPhiNonCommon[iPart1][iPart2], "");
    outDist->SetMEDistribution(fMEdEtadPhi[iPart1][iPart2], "");
  }
  return outDist;
}

DreamdEtadPhi* ReadDreamFile::GetdEtadPhiAtRadDistribution(int iPart1,
                                                           int iPart2,
                                                           int iMix1,
                                                           int iAPart1,
                                                           int iAPart2,
                                                           int iMix2, int iRad,
                                                           bool smallkStar) {
  DreamdEtadPhi* outDist = new DreamdEtadPhi();
  outDist->SetSEDistribution(
      smallkStar ?
          fSEdEtadPhiAtRadSmallkStar[iPart1][iPart2][iRad][iMix1] :
          fSEdEtadPhiAtRad[iPart1][iPart2][iRad][iMix1],
      "");
  outDist->AddSEDistribution(
      smallkStar ?
          fSEdEtadPhiAtRadSmallkStar[iAPart1][iAPart2][iRad][iMix2] :
          fSEdEtadPhiAtRad[iAPart1][iAPart2][iRad][iMix2]);

  outDist->SetMEDistribution(
      smallkStar ?
          fMEdEtadPhiAtRadSmallkStar[iPart1][iPart2][iRad][iMix1] :
          fMEdEtadPhiAtRad[iPart1][iPart2][iRad][iMix2],
      "");
  outDist->AddMEDistribution(
      smallkStar ?
          fMEdEtadPhiAtRadSmallkStar[iAPart1][iAPart2][iRad][iMix1] :
          fMEdEtadPhiAtRad[iAPart1][iAPart2][iRad][iMix2]);
  return outDist;
}

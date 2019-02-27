/*
 * CATSInput.cxx
 *
 *  Created on: Oct 30, 2018
 *      Author: hohlweger
 */

#include <CATSInput.h>
#include "TFile.h"
#include <iostream>
#include "TSystem.h"
#include "TROOT.h"
#include "TF1.h"

CATSInput::CATSInput()
    : fFixBinningExternal(false),
      fFixkMin(0.),
      fNameBasedir(),
      fNameMomResFile(),
      fNameSigmaFile(),
      fFraction_Res(0),
      fFraction_Sig(0),
      fUnitConv_Res(0),
      fUnitConv_Sig(0),
      fRes(),
      fSigma(),
      fDreamFile(nullptr),
      fCF_pp(nullptr),
      fCF_pL(nullptr),
      fCF_LL(nullptr),
      fCF_pXi(nullptr),
      fnormalizationLeft(-100),
      fnormalizationRight(-100) {
  // TODO Auto-generated constructor stub
  TH1::AddDirectory(kFALSE);
  TH2::AddDirectory(kFALSE);
}

CATSInput::~CATSInput() {
  if (fDreamFile)
    delete fDreamFile;
  if (fCF_pp)
    delete fCF_pp;
  if (fCF_pL)
    delete fCF_pL;
  if (fCF_LL)
    delete fCF_LL;
  if (fCF_pXi)
    delete fCF_pXi;

  for (auto it : fRes) {
    if (it)
      delete it;
  }
  for (auto it : fSigma) {
    if (it)
      delete it;
  }

}

void CATSInput::ReadResFile() {
  TString ResMatrixFileName = "";
  ResMatrixFileName += fNameBasedir;
  ResMatrixFileName += fNameMomResFile;
  auto FileRes = new TFile(ResMatrixFileName, "read");
  if (!FileRes) {
    std::cout << "No Resolution file found in " << ResMatrixFileName
              << std::endl;
  } else {
    FileRes->cd();
    std::vector<TH2F*> inputHist;
    std::vector<const char*> FileNames = { "hRes_pp_pL", "hRes_pL_pSigma0",
        "hRes_pL_pXim", "hRes_pXim_pXim1530" };
    for (auto it : FileNames) {
      auto hist = (TH2F*) FileRes->Get(it);
      if (hist) {
        inputHist.push_back(
            (TH2F*) hist->Clone(Form("%s_clone", hist->GetName())));
        delete hist;
      } else {
        std::cout << it << " histogram not found \n";
      }
    }
    if (inputHist.size() != FileNames.size()) {
      std::cout
          << "Something went horrible wrong while reading the input Resolution hists \n";
      return;
    }
    for (auto it : inputHist) {
      TString histName = Form("%s_MeV", it->GetName());
      TH2F* tmpResMeV = new TH2F(
          histName.Data(),
          histName.Data(),
          it->GetNbinsX() / fFraction_Res,
          it->GetXaxis()->GetBinLowEdge(1) * fUnitConv_Res,
          it->GetXaxis()->GetBinUpEdge(it->GetNbinsX() / fFraction_Res)
              * fUnitConv_Res,
          it->GetNbinsY() / fFraction_Res,
          it->GetYaxis()->GetBinLowEdge(1) * fUnitConv_Res,
          it->GetXaxis()->GetBinUpEdge(it->GetNbinsY() / fFraction_Res)
              * fUnitConv_Res);
      for (int iBinX = 1; iBinX <= it->GetNbinsX() / fFraction_Res; iBinX++) {
        for (int iBinY = 1; iBinY <= it->GetNbinsY() / fFraction_Res; iBinY++) {
          tmpResMeV->SetBinContent(iBinX, iBinY,
                                   it->GetBinContent(iBinX, iBinY));
        }
      }
      fRes.push_back(tmpResMeV);
      delete it;
    }
  }
  FileRes->Close();
  delete FileRes;
  return;
}

void CATSInput::ReadSigmaFile() {
  TString SigmaMatrixFileName = "";
  SigmaMatrixFileName += fNameBasedir;
  SigmaMatrixFileName += fNameSigmaFile;
  auto FileSigma = new TFile(SigmaMatrixFileName, "read");
  if (!FileSigma) {
    std::cout << "No Sigma file found in " << SigmaMatrixFileName << std::endl;
  } else {
    FileSigma->cd();
    std::vector<TH2F*> inputHist;
    std::vector<const char*> FileNames = { "hSigmaMeV_Proton_Proton",
        "hSigmaMeV_Proton_Lambda", "hSigmaMeV_Lambda_Lambda",
        "hSigmaMeV_Proton_Xim" };
    for (auto it : FileNames) {
      auto hist = (TH2F*) FileSigma->Get(it);
      if (hist) {
        inputHist.push_back(
            (TH2F*) hist->Clone(Form("%s_clone", hist->GetName())));
        delete hist;
      } else {
        std::cout << it << " histogram not found \n";
      }
    }
    if (inputHist.size() != FileNames.size()) {
      std::cout
          << "Something went horrible wrong while reading the input Resolution hists \n";
      return;
    }
    for (auto it : inputHist) {
      TString histName = Form("%s_MeV", it->GetName());
      TH2F* tmpSigmaMeV = new TH2F(
          histName.Data(),
          histName.Data(),
          it->GetNbinsX() / fFraction_Sig,
          it->GetXaxis()->GetBinLowEdge(1) * fUnitConv_Sig,
          it->GetXaxis()->GetBinUpEdge(it->GetNbinsX() / fFraction_Sig)
              * fUnitConv_Sig,
          it->GetNbinsY() / fFraction_Sig,
          it->GetYaxis()->GetBinLowEdge(1) * fUnitConv_Sig,
          it->GetXaxis()->GetBinUpEdge(it->GetNbinsY() / fFraction_Sig)
              * fUnitConv_Sig);
      for (int iBinX = 1; iBinX <= it->GetNbinsX() / fFraction_Sig; iBinX++) {
        for (int iBinY = 1; iBinY <= it->GetNbinsY() / fFraction_Sig; iBinY++) {
          tmpSigmaMeV->SetBinContent(iBinX, iBinY,
                                     it->GetBinContent(iBinX, iBinY));
        }
      }
      fSigma.push_back(tmpSigmaMeV);
      delete it;
    }
  }
  FileSigma->Close();
  delete FileSigma;
  return;
}

void CATSInput::ReadCorrelationFile(const char* path, const char* prefix,
                                    const char* suffix) {
  TString filename = Form("%s/AnalysisResults.root", path);
  fDreamFile = new ReadDreamFile(6, 6);
  fDreamFile->SetAnalysisFile(filename.Data(), prefix, suffix);
  return;
}

void CATSInput::ObtainCFs(int rebin, float normleft, float normright) {
  //normleft & right in MeV!
  normleft /= 1000.;
  normright /= 1000.;

  if (!fDreamFile) {
    std::cout << "No File was set via ReadCorrelationFile\n";
  } else {
    if (fnormalizationLeft != normleft || fnormalizationRight != normright) {
      if (fCF_pp) {
        delete fCF_pp;
      }
      if (fCF_pL) {
        delete fCF_pL;
      }
      if (fCF_LL) {
        delete fCF_LL;
      }
      if (fCF_pXi) {
        delete fCF_pXi;
      }
      fCF_pp = new DreamCF();
      DreamPair* pp = new DreamPair("Part", normleft, normright);
      DreamPair* ApAp = new DreamPair("AntiPart", normleft, normright);

      fCF_pL = new DreamCF();
      DreamPair* pL = new DreamPair("Part", normleft, normright);
      DreamPair* ApAL = new DreamPair("AntiPart", normleft, normright);

      fCF_LL = new DreamCF();
      DreamPair* LL = new DreamPair("Part", normleft, normright);
      DreamPair* ALAL = new DreamPair("AntiPart", normleft, normright);

      fCF_pXi = new DreamCF();
      DreamPair* pXi = new DreamPair("Part", normleft, normright);
      DreamPair* ApAXi = new DreamPair("AntiPart", normleft, normright);
      std::cout << "Set pair\n";
      pp->SetPair(fDreamFile->GetPairDistributions(0, 0, ""));
      ApAp->SetPair(fDreamFile->GetPairDistributions(1, 1, ""));

      pL->SetPair(fDreamFile->GetPairDistributions(0, 2, ""));
      ApAL->SetPair(fDreamFile->GetPairDistributions(1, 3, ""));

      LL->SetPair(fDreamFile->GetPairDistributions(2, 2, ""));
      ALAL->SetPair(fDreamFile->GetPairDistributions(3, 3, ""));

      pXi->SetPair(fDreamFile->GetPairDistributions(0, 4, ""));
      ApAXi->SetPair(fDreamFile->GetPairDistributions(1, 5, ""));
      pp->ShiftForEmpty(pp->GetPair());
      ApAp->ShiftForEmpty(ApAp->GetPair());

      pL->ShiftForEmpty(pL->GetPair());
      ApAL->ShiftForEmpty(ApAL->GetPair());

      LL->ShiftForEmpty(LL->GetPair());
      ALAL->ShiftForEmpty(ALAL->GetPair());

      pXi->ShiftForEmpty(pXi->GetPair());
      ApAXi->ShiftForEmpty(ApAXi->GetPair());
      pp->FixShift(pp->GetPairShiftedEmpty(0), ApAp->GetPairShiftedEmpty(0),
                   ApAp->GetFirstBin());
      ApAp->FixShift(ApAp->GetPairShiftedEmpty(0), pp->GetPairShiftedEmpty(0),
                     pp->GetFirstBin());

      pL->FixShift(pL->GetPairShiftedEmpty(0), ApAL->GetPairShiftedEmpty(0),
                   ApAL->GetFirstBin());
      ApAL->FixShift(ApAL->GetPairShiftedEmpty(0), pL->GetPairShiftedEmpty(0),
                     pL->GetFirstBin());

      LL->FixShift(LL->GetPairShiftedEmpty(0), ALAL->GetPairShiftedEmpty(0),
                   ALAL->GetFirstBin());
      ALAL->FixShift(ALAL->GetPairShiftedEmpty(0), LL->GetPairShiftedEmpty(0),
                     LL->GetFirstBin());

      pXi->FixShift(pXi->GetPairShiftedEmpty(0), ApAXi->GetPairShiftedEmpty(0),
                    ApAXi->GetFirstBin());
      ApAXi->FixShift(ApAXi->GetPairShiftedEmpty(0),
                      pXi->GetPairShiftedEmpty(0), pXi->GetFirstBin());
      pL->Rebin(pL->GetPairFixShifted(0), rebin);
      ApAL->Rebin(ApAL->GetPairFixShifted(0), rebin);

      LL->Rebin(LL->GetPairFixShifted(0), rebin);
      ALAL->Rebin(ALAL->GetPairFixShifted(0), rebin);

      pXi->Rebin(pXi->GetPairFixShifted(0), rebin);
      ApAXi->Rebin(ApAXi->GetPairFixShifted(0), rebin);

      pp->ReweightMixedEvent(pp->GetPairFixShifted(0), 0.2, 0.9);
      ApAp->ReweightMixedEvent(ApAp->GetPairFixShifted(0), 0.2, 0.9);

      pL->ReweightMixedEvent(pL->GetPairRebinned(0), 0.2, 0.9);
      ApAL->ReweightMixedEvent(ApAL->GetPairRebinned(0), 0.2, 0.9);

      LL->ReweightMixedEvent(LL->GetPairRebinned(0), 0.2, 0.9);
      ALAL->ReweightMixedEvent(ALAL->GetPairRebinned(0), 0.2, 0.9);

      pXi->ReweightMixedEvent(pXi->GetPairRebinned(0), 0.2, 0.9);
      ApAXi->ReweightMixedEvent(ApAXi->GetPairRebinned(0), 0.2, 0.9);

      fCF_pp->SetPairs(pp, ApAp);
      fCF_pp->GetCorrelations("pp");

      fCF_pL->SetPairs(pL, ApAL);
      fCF_pL->GetCorrelations("pL");

      fCF_LL->SetPairs(LL, ALAL);
      fCF_LL->GetCorrelations("LL");

      fCF_pXi->SetPairs(pXi, ApAXi);
      fCF_pXi->GetCorrelations("pXi");
      fnormalizationLeft = normleft;
      fnormalizationRight = normright;
    } else {
      std::cout << "Already existing normalization \n";
    }
  }
}

DreamCF* CATSInput::ObtainCFSyst(int rebin, const char* name, DreamDist* ppDist,
                                 DreamDist* ApApDist, DreamDist* ppFake,
                                 DreamDist* ApApFake) {
  //normleft & right in MeV!
  DreamCF* outCF = new DreamCF();
  DreamPair* pp = new DreamPair("Part", fnormalizationLeft,
                                fnormalizationRight);
  DreamPair* ApAp = new DreamPair("AntiPart", fnormalizationLeft,
                                  fnormalizationRight);
  if (fnormalizationLeft == 0 || fnormalizationRight == 0) {
    std::cout << "Normalization is 0! Bad results incoming! \n";
  }
  std::cout << "Set pair\n";
  pp->SetPair(ppDist);
  ApAp->SetPair(ApApDist);

  if (ppFake) {
    std::cout << "Faking SE Mult for " << name << std::endl;
    pp->GetPair()->SetSEMultDist(ppFake->GetSEMultDist(), "");
    pp->GetPair()->SetMEMultDist(ppFake->GetMEMultDist(), "");
  }

  if (ApApFake) {
    std::cout << "Faking ME Mult for " << name << std::endl;
    ApAp->GetPair()->SetSEMultDist(ApApFake->GetSEMultDist(), "");
    ApAp->GetPair()->SetMEMultDist(ApApFake->GetMEMultDist(), "");
  }

  pp->ShiftForEmpty(pp->GetPair());
  ApAp->ShiftForEmpty(ApAp->GetPair());

  if (fFixBinningExternal) {
    pp->FixShift(pp->GetPair(), ApAp->GetPair(),
                 fFixkMin, true);
    ApAp->FixShift(ApAp->GetPair(), pp->GetPair(),
                   fFixkMin, true);
  } else {
    pp->FixShift(pp->GetPairShiftedEmpty(0), ApAp->GetPairShiftedEmpty(0),
                 ApAp->GetFirstBin());
    ApAp->FixShift(ApAp->GetPairShiftedEmpty(0), pp->GetPairShiftedEmpty(0),
                   pp->GetFirstBin());
  }
  pp->ReweightMixedEvent(pp->GetPairFixShifted(0), 0.2, 0.9);
  ApAp->ReweightMixedEvent(ApAp->GetPairFixShifted(0), 0.2, 0.9);

  if (rebin != 1) {
    pp->Rebin(pp->GetPairFixShifted(0), rebin);
    ApAp->Rebin(ApAp->GetPairFixShifted(0), rebin);
    pp->ReweightMixedEvent(pp->GetPairRebinned(0), 0.2, 0.9);
    ApAp->ReweightMixedEvent(ApAp->GetPairRebinned(0), 0.2, 0.9);
  }

  outCF->SetPairs(pp, ApAp);
  outCF->GetCorrelations(name);
  return outCF;
}

TH1F* CATSInput::GetCF(TString pair, TString hist) {
  TH1F* output = nullptr;
  if (pair == TString("pp")) {
    for (auto it : fCF_pp->GetCorrelationFunctions()) {
      TString itName = it->GetName();
      std::cout << it->GetName() << std::endl;
      if (hist == itName) {
        output = it;
      }
    }
  } else if (pair == TString("pL")) {
    for (auto it : fCF_pL->GetCorrelationFunctions()) {
      TString itName = it->GetName();
      if (hist == itName) {
        std::cout << it->GetName() << std::endl;
        output = it;
      }
    }
  } else if (pair == TString("LL")) {
    for (auto it : fCF_LL->GetCorrelationFunctions()) {
      TString itName = it->GetName();
      if (hist == itName) {
        std::cout << it->GetName() << std::endl;
        output = it;
      }
    }
  } else if (pair == TString("pXi")) {
    for (auto it : fCF_pXi->GetCorrelationFunctions()) {
      TString itName = it->GetName();
      if (hist == itName) {
        std::cout << it->GetName() << std::endl;
        output = it;
      }
    }
  } else {
    std::cout << pair << " does not exist\n";
  }
  if (!output) {
    std::cout << "Danger! Histogram not set, maybe histname " << hist
              << " does not exist? \n";
  }
  return (TH1F*) output->Clone(Form("%sCloned", output->GetName()));
}

void CATSInput::AddSystematics(TString SysFile, TH1F* Hist) {
//!CHANGE PATH HERE
  TString SystErrFileName = TString::Format("%s%s", fNameBasedir.Data(),
                                            SysFile.Data());
  TFile* SystErrFile =
      SystErrFileName != "" ? new TFile(SystErrFileName, "read") : nullptr;
  if (SystErrFile) {
    TString sysName = SysFile;
    sysName.Replace(0, 10, "");
    sysName.Replace(2, 6, "");
    TH1F* outputParam = (TH1F*) SystErrFile->Get(
        Form("SysParam%s", sysName.Data()));
    if (outputParam) {
      TF1 *RelSyst = new TF1("sys", "pol2", 0, 3);
      RelSyst->SetParameter(0, outputParam->GetBinContent(1));
      RelSyst->SetParameter(1, outputParam->GetBinContent(2));
      RelSyst->SetParameter(2, outputParam->GetBinContent(3));

      int NumSEB = RelSyst == NULL ? 0 : Hist->FindBin(500);

      for (int iBin = 0; iBin < NumSEB; iBin++) {
        const float x = Hist->GetBinCenter(iBin + 1);
        const float y = Hist->GetBinContent(iBin + 1);
        Hist->SetBinError(
            iBin + 1,
            sqrt(
                pow(Hist->GetBinError(iBin + 1), 2.)
                    + pow(y * RelSyst->Eval(x / 1000.), 2.)));
      }
      delete RelSyst;
    }
  } else {
    std::cout << "Sys File not found " << SystErrFileName.Data() << std::endl;
  }
}

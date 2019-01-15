#include "CATSInputSigma0.h"
#include <iostream>

CATSInputSigma0::CATSInputSigma0()
    : fCF_pSigma(nullptr),
      fCF_SidebandUp(nullptr),
      fCF_SidebandLow(nullptr) {
}

CATSInputSigma0::~CATSInputSigma0() {
  delete fCF_pSigma;
  delete fCF_SidebandUp;
  delete fCF_SidebandLow;
}

void CATSInputSigma0::ReadSigma0CorrelationFile(const char* path,
                                                const char* appendix) {
  TString filename = Form("%s/AnalysisResults.root", path);
  fDreamFile = new ReadDreamFile(12, 12);
  fDreamFile->SetSigmaAnalysisFile(filename.Data(), appendix);
  return;
}

void CATSInputSigma0::ObtainCFs(int rebin, float normleft, float normright,
                                int rebinSyst) {
//normleft & right in MeV!
  normleft /= 1000.;
  normright /= 1000.;

  if (!fDreamFile) {
    std::cerr
        << "ERROR CATSInputSigma0: No File was set via ReadSigma0CorrelationFile\n";
    return;
  }
  if (fnormalizationLeft == normleft & fnormalizationRight == normright) {
    std::cerr << "ERROR CATSInputSigma0: Already existing normalization \n";
    return;
  }

  fCF_pp = new DreamCF();
  DreamPair* pp = new DreamPair("Part", normleft, normright);
  DreamPair* ApAp = new DreamPair("AntiPart", normleft, normright);

  fCF_pSigma = new DreamCF();
  DreamPair* pSigma = new DreamPair("Part", normleft, normright);
  DreamPair* ApASigma = new DreamPair("AntiPart", normleft, normright);

  fCF_SidebandUp = new DreamCF();
  DreamPair* pSiSBup = new DreamPair("Part_SB_up", normleft, normright);
  DreamPair* ApaSiSBup = new DreamPair("AntiPart_SB_up", normleft, normright);

  fCF_SidebandLow = new DreamCF();
  DreamPair* pSiSBlow = new DreamPair("Part_SB_low", normleft, normright);
  DreamPair* ApaSiSBlow = new DreamPair("AntiPart_SB_low", normleft, normright);

  pp->SetPair(fDreamFile->GetPairDistributions(0, 0, ""));
  ApAp->SetPair(fDreamFile->GetPairDistributions(1, 1, ""));

  pSigma->SetPair(fDreamFile->GetPairDistributions(0, 2, ""));
  ApASigma->SetPair(fDreamFile->GetPairDistributions(1, 3, ""));

  pSiSBup->SetPair(fDreamFile->GetPairDistributions(0, 4, ""));
  ApaSiSBup->SetPair(fDreamFile->GetPairDistributions(1, 5, ""));

  pSiSBlow->SetPair(fDreamFile->GetPairDistributions(0, 6, ""));
  ApaSiSBlow->SetPair(fDreamFile->GetPairDistributions(1, 7, ""));

  pp->ShiftForEmpty(pp->GetPair());
  ApAp->ShiftForEmpty(ApAp->GetPair());

  pSigma->ShiftForEmpty(pSigma->GetPair());
  ApASigma->ShiftForEmpty(ApASigma->GetPair());

  pSiSBup->ShiftForEmpty(pSiSBup->GetPair());
  ApaSiSBup->ShiftForEmpty(ApaSiSBup->GetPair());

  pSiSBlow->ShiftForEmpty(pSiSBlow->GetPair());
  ApaSiSBlow->ShiftForEmpty(ApaSiSBlow->GetPair());

  pp->FixShift(pp->GetPairShiftedEmpty(0), ApAp->GetPairShiftedEmpty(0),
               ApAp->GetFirstBin());
  ApAp->FixShift(ApAp->GetPairShiftedEmpty(0), pp->GetPairShiftedEmpty(0),
                 pp->GetFirstBin());

  pSigma->FixShift(pSigma->GetPairShiftedEmpty(0),
                   ApASigma->GetPairShiftedEmpty(0), ApASigma->GetFirstBin());
  ApASigma->FixShift(ApASigma->GetPairShiftedEmpty(0),
                     pSigma->GetPairShiftedEmpty(0), pSigma->GetFirstBin());

  pSiSBup->FixShift(pSiSBup->GetPairShiftedEmpty(0),
                    pSigma->GetPairShiftedEmpty(0),
                    ApASigma->GetPairShiftedEmpty(0), pSigma->GetFirstBin(),
                    ApASigma->GetFirstBin());
  ApaSiSBup->FixShift(ApaSiSBup->GetPairShiftedEmpty(0),
                      pSigma->GetPairShiftedEmpty(0),
                      ApASigma->GetPairShiftedEmpty(0), pSigma->GetFirstBin(),
                      ApASigma->GetFirstBin());

  pSiSBlow->FixShift(pSiSBlow->GetPairShiftedEmpty(0),
                     pSigma->GetPairShiftedEmpty(0),
                     ApASigma->GetPairShiftedEmpty(0), pSigma->GetFirstBin(),
                     ApASigma->GetFirstBin());
  ApaSiSBlow->FixShift(ApaSiSBlow->GetPairShiftedEmpty(0),
                       pSigma->GetPairShiftedEmpty(0),
                       ApASigma->GetPairShiftedEmpty(0), pSigma->GetFirstBin(),
                       ApASigma->GetFirstBin());

  if (rebinSyst != 1) {
    pp->Rebin(pp->GetPairFixShifted(0), rebinSyst);
    ApAp->Rebin(ApAp->GetPairFixShifted(0), rebinSyst);
  }

  pSigma->Rebin(pSigma->GetPairFixShifted(0), rebin * rebinSyst);
  ApASigma->Rebin(ApASigma->GetPairFixShifted(0), rebin * rebinSyst);

  pSiSBup->Rebin(pSiSBup->GetPairFixShifted(0), rebin * rebinSyst);
  ApaSiSBup->Rebin(ApaSiSBup->GetPairFixShifted(0), rebin * rebinSyst);

  pSiSBlow->Rebin(pSiSBlow->GetPairFixShifted(0), rebin * rebinSyst);
  ApaSiSBlow->Rebin(ApaSiSBlow->GetPairFixShifted(0), rebin * rebinSyst);

  if (rebinSyst != 1) {
    pp->ReweightMixedEvent(pp->GetPairRebinned(0), 0.2, 0.9);
    ApAp->ReweightMixedEvent(ApAp->GetPairRebinned(0), 0.2, 0.9);
  } else {
    pp->ReweightMixedEvent(pp->GetPairFixShifted(0), 0.2, 0.9);
    ApAp->ReweightMixedEvent(ApAp->GetPairFixShifted(0), 0.2, 0.9);
  }

  pSigma->ReweightMixedEvent(pSigma->GetPairRebinned(0), 0.2, 0.9);
  ApASigma->ReweightMixedEvent(ApASigma->GetPairRebinned(0), 0.2, 0.9);

  pSiSBup->ReweightMixedEvent(pSiSBup->GetPairRebinned(0), 0.2, 0.9);
  ApaSiSBup->ReweightMixedEvent(ApaSiSBup->GetPairRebinned(0), 0.2, 0.9);

  pSiSBlow->ReweightMixedEvent(pSiSBlow->GetPairRebinned(0), 0.2, 0.9);
  ApaSiSBlow->ReweightMixedEvent(ApaSiSBlow->GetPairRebinned(0), 0.2, 0.9);

  fCF_pp->SetPairs(pp, ApAp);
  fCF_pp->GetCorrelations("pp");

  fCF_pSigma->SetPairs(pSigma, ApASigma);
  fCF_pSigma->GetCorrelations("pSigma0");

  fCF_SidebandUp->SetPairs(pSiSBup, ApaSiSBup);
  fCF_SidebandUp->GetCorrelations("pSigmaSBUp");

  fCF_SidebandLow->SetPairs(pSiSBlow, ApaSiSBlow);
  fCF_SidebandLow->GetCorrelations("pSigmaSBLow");
  fnormalizationLeft = normleft;
  fnormalizationRight = normright;
}

TH1F* CATSInputSigma0::GetCF(TString pair, TString hist) {
  TH1F* output = nullptr;
  if (pair == TString("pp")) {
    for (const auto &it : fCF_pp->GetCorrelationFunctions()) {
      TString itName = it->GetName();
      if (hist == itName) {
        output = it;
      }
    }
  } else if (pair == TString("pSigma0")) {
    for (const auto &it : fCF_pSigma->GetCorrelationFunctions()) {
      TString itName = it->GetName();
      if (hist == itName) {
        output = it;
      }
    }
  } else if (pair == TString("pSigmaSBUp")) {
    for (const auto &it : fCF_SidebandUp->GetCorrelationFunctions()) {
      TString itName = it->GetName();
      if (hist == itName) {
        output = it;
      }
    }
  } else if (pair == TString("pSigmaSBLow")) {
    for (const auto &it : fCF_SidebandLow->GetCorrelationFunctions()) {
      TString itName = it->GetName();
      if (hist == itName) {
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

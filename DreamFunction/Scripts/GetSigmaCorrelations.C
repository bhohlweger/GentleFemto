#include <iostream>
#include "ReadDreamFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TROOT.h"

void GetSigmaCorrelations(const char* filename, const char* suffix) {
  // gStyle->SetOptStat(0);
  ReadDreamFile* DreamFile = new ReadDreamFile(12, 12);
  DreamFile->SetSigmaAnalysisFile(filename, suffix);

  DreamCF* CF_pp = new DreamCF();
  DreamPair* pp = new DreamPair("Part", 0.24, 0.34);
  DreamPair* ApAp = new DreamPair("AntiPart", 0.24, 0.34);

  DreamCF* CF_pSigma = new DreamCF();
  DreamPair* pSigma = new DreamPair("Part", 0.24, 0.34);
  DreamPair* ApASigma = new DreamPair("AntiPart", 0.24, 0.34);

  DreamCF* CF_SidebandUp = new DreamCF();
  DreamPair* pSiSBup = new DreamPair("Part_SB_up", 0.24, 0.34);
  DreamPair* ApaSiSBup = new DreamPair("AntiPart_SB_up", 0.24, 0.34);

  DreamCF* CF_SidebandLow = new DreamCF();
  DreamPair* pSiSBlow = new DreamPair("Part_SB_low", 0.24, 0.34);
  DreamPair* ApaSiSBlow = new DreamPair("AntiPart_SB_low", 0.24, 0.34);

  DreamCF* CF_SigmaLambda = new DreamCF();
  DreamPair* pLambdaSigma = new DreamPair("PartDaughter", 0.24, 0.34);
  DreamPair* ApALambdaSigma = new DreamPair("AntiPartDaughter", 0.24, 0.34);

  DreamCF* CF_SigmaPhoton = new DreamCF();
  DreamPair* pPhotonSigma = new DreamPair("PartDaughter", 0.24, 0.34);
  DreamPair* ApAPhotonSigma = new DreamPair("AntiPartDaughter", 0.24, 0.34);

  std::cout << "=========================" << std::endl;
  std::cout << "========Pair Set=========" << std::endl;
  std::cout << "=========================" << std::endl;
  pp->SetPair(DreamFile->GetPairDistributions(0, 0, ""));
  ApAp->SetPair(DreamFile->GetPairDistributions(1, 1, ""));

  pSigma->SetPair(DreamFile->GetPairDistributions(0, 2, ""));
  ApASigma->SetPair(DreamFile->GetPairDistributions(1, 3, ""));

  pSiSBup->SetPair(DreamFile->GetPairDistributions(0, 4, ""));
  ApaSiSBup->SetPair(DreamFile->GetPairDistributions(1, 5, ""));

  pSiSBlow->SetPair(DreamFile->GetPairDistributions(0, 6, ""));
  ApaSiSBlow->SetPair(DreamFile->GetPairDistributions(1, 7, ""));

  pLambdaSigma->SetPair(DreamFile->GetPairDistributions(0, 8, ""));
  ApALambdaSigma->SetPair(DreamFile->GetPairDistributions(1, 10, ""));

  pPhotonSigma->SetPair(DreamFile->GetPairDistributions(0, 9, ""));
  ApAPhotonSigma->SetPair(DreamFile->GetPairDistributions(1, 11, ""));

  std::cout << "=========================" << std::endl;
  std::cout << "======Pair Shifted=======" << std::endl;
  std::cout << "=========================" << std::endl;
  pp->ShiftForEmpty(pp->GetPair());
  ApAp->ShiftForEmpty(ApAp->GetPair());

  pSigma->ShiftForEmpty(pSigma->GetPair());
  ApASigma->ShiftForEmpty(ApASigma->GetPair());

  pSiSBup->ShiftForEmpty(pSiSBup->GetPair());
  ApaSiSBup->ShiftForEmpty(ApaSiSBup->GetPair());

  pSiSBlow->ShiftForEmpty(pSiSBlow->GetPair());
  ApaSiSBlow->ShiftForEmpty(ApaSiSBlow->GetPair());

  pLambdaSigma->ShiftForEmpty(pLambdaSigma->GetPair());
  ApALambdaSigma->ShiftForEmpty(ApALambdaSigma->GetPair());

  pPhotonSigma->ShiftForEmpty(pPhotonSigma->GetPair());
  ApAPhotonSigma->ShiftForEmpty(ApAPhotonSigma->GetPair());

  std::cout << "=========================" << std::endl;
  std::cout << "====Pair Fix Shifted=====" << std::endl;
  std::cout << "=========================" << std::endl;
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

  pLambdaSigma->FixShift(pLambdaSigma->GetPairShiftedEmpty(0),
                         pSigma->GetPairShiftedEmpty(0),
                         ApASigma->GetPairShiftedEmpty(0),
                         pSigma->GetFirstBin(), ApASigma->GetFirstBin());
  ApALambdaSigma->FixShift(ApALambdaSigma->GetPairShiftedEmpty(0),
                           pSigma->GetPairShiftedEmpty(0),
                           ApASigma->GetPairShiftedEmpty(0),
                           pSigma->GetFirstBin(), ApASigma->GetFirstBin());

  pPhotonSigma->FixShift(pPhotonSigma->GetPairShiftedEmpty(0),
                         pSigma->GetPairShiftedEmpty(0),
                         ApASigma->GetPairShiftedEmpty(0),
                         pSigma->GetFirstBin(), ApASigma->GetFirstBin());
  ApAPhotonSigma->FixShift(ApAPhotonSigma->GetPairShiftedEmpty(0),
                           pSigma->GetPairShiftedEmpty(0),
                           ApASigma->GetPairShiftedEmpty(0),
                           pSigma->GetFirstBin(), ApASigma->GetFirstBin());

  std::cout << "=========================" << std::endl;
  std::cout << "==Rebinning & Weighting==" << std::endl;
  std::cout << "=========================" << std::endl;

  for (int iReb = 4; iReb < 6; ++iReb) {
    std::cout << "==Rebinning==" << std::endl;
    pSigma->Rebin(pSigma->GetPairFixShifted(0), iReb);
    ApASigma->Rebin(ApASigma->GetPairFixShifted(0), iReb);
    pSiSBup->Rebin(pSiSBup->GetPairFixShifted(0), iReb);
    ApaSiSBup->Rebin(ApaSiSBup->GetPairFixShifted(0), iReb);
    pSiSBlow->Rebin(pSiSBlow->GetPairFixShifted(0), iReb);
    ApaSiSBlow->Rebin(ApaSiSBlow->GetPairFixShifted(0), iReb);
    pLambdaSigma->Rebin(pLambdaSigma->GetPairFixShifted(0), iReb);
    ApALambdaSigma->Rebin(ApALambdaSigma->GetPairFixShifted(0), iReb);
    pPhotonSigma->Rebin(pPhotonSigma->GetPairFixShifted(0), iReb);
    ApAPhotonSigma->Rebin(ApAPhotonSigma->GetPairFixShifted(0), iReb);
    std::cout << "==Weighting==" << std::endl;
    pSigma->ReweightMixedEvent(pSigma->GetPairRebinned(iReb - 4), 0.2, 0.9);
    ApASigma->ReweightMixedEvent(ApASigma->GetPairRebinned(iReb - 4), 0.2, 0.9);
    pSiSBup->ReweightMixedEvent(pSiSBup->GetPairRebinned(iReb - 4), 0.2, 0.9);
    ApaSiSBup->ReweightMixedEvent(ApaSiSBup->GetPairRebinned(iReb - 4), 0.2, 0.9);
    pSiSBlow->ReweightMixedEvent(pSiSBlow->GetPairRebinned(iReb - 4), 0.2, 0.9);
    ApaSiSBlow->ReweightMixedEvent(ApaSiSBlow->GetPairRebinned(iReb - 4), 0.2, 0.9);
    pLambdaSigma->ReweightMixedEvent(pLambdaSigma->GetPairRebinned(iReb - 4),
                                     0.2, 0.9);
    ApALambdaSigma->ReweightMixedEvent(
        ApALambdaSigma->GetPairRebinned(iReb - 4), 0.2, 0.9);
    pPhotonSigma->ReweightMixedEvent(pPhotonSigma->GetPairRebinned(iReb - 4),
                                     0.2, 0.9);
    ApAPhotonSigma->ReweightMixedEvent(
        ApAPhotonSigma->GetPairRebinned(iReb - 4), 0.2, 0.9);
  }
  pp->ReweightMixedEvent(pp->GetPairFixShifted(0), 0.2, 0.9);
  ApAp->ReweightMixedEvent(ApAp->GetPairFixShifted(0), 0.2, 0.9);

  pSigma->Rebin(pSigma->GetPair(), 4);
  pSigma->Rebin(pSigma->GetPair(), 5);
  ApASigma->Rebin(ApASigma->GetPair(), 4);
  ApASigma->Rebin(ApASigma->GetPair(), 5);
  pSiSBup->Rebin(pSiSBup->GetPair(), 4);
  pSiSBup->Rebin(pSiSBup->GetPair(), 5);
  ApaSiSBup->Rebin(ApaSiSBup->GetPair(), 4);
  ApaSiSBup->Rebin(ApaSiSBup->GetPair(), 5);
  pSiSBlow->Rebin(pSiSBlow->GetPair(), 4);
  pSiSBlow->Rebin(pSiSBlow->GetPair(), 5);
  ApaSiSBlow->Rebin(ApaSiSBlow->GetPair(), 4);
  ApaSiSBlow->Rebin(ApaSiSBlow->GetPair(), 5);
  pLambdaSigma->Rebin(pLambdaSigma->GetPair(), 4);
  pLambdaSigma->Rebin(pLambdaSigma->GetPair(), 5);
  ApALambdaSigma->Rebin(ApALambdaSigma->GetPair(), 4);
  ApALambdaSigma->Rebin(ApALambdaSigma->GetPair(), 5);
  pPhotonSigma->Rebin(pPhotonSigma->GetPair(), 4);
  pPhotonSigma->Rebin(pPhotonSigma->GetPair(), 5);
  ApAPhotonSigma->Rebin(ApAPhotonSigma->GetPair(), 4);
  ApAPhotonSigma->Rebin(ApAPhotonSigma->GetPair(), 5);

  std::cout << "=========================" << std::endl;
  std::cout << "=========CFs=============" << std::endl;
  std::cout << "=========================" << std::endl;

  TString foldername = filename;
  foldername.ReplaceAll("AnalysisResults.root", "");

  CF_pp->SetPairs(pp, ApAp);
  CF_pp->GetCorrelations();
  CF_pp->WriteOutput(Form("%s/CFOutput_pp.root", foldername.Data()));

  CF_pSigma->SetPairs(pSigma, ApASigma);
  CF_pSigma->GetCorrelations();
  CF_pSigma->WriteOutput(Form("%s/CFOutput_pSigma.root", foldername.Data()));

  CF_SidebandUp->SetPairs(pSiSBup, ApaSiSBup);
  CF_SidebandUp->GetCorrelations();
  CF_SidebandUp->WriteOutput(Form("%s/CFOutput_SB_up.root", foldername.Data()));

  CF_SidebandLow->SetPairs(pSiSBlow, ApaSiSBlow);
  CF_SidebandLow->GetCorrelations();
  CF_SidebandLow->WriteOutput(
      Form("%s/CFOutput_SB_low.root", foldername.Data()));

  CF_SigmaLambda->SetPairs(pLambdaSigma, ApALambdaSigma);
  CF_SigmaLambda->GetCorrelations();
  CF_SigmaLambda->WriteOutput(
      Form("%s/CFOutput_pLambda.root", foldername.Data()));

  CF_SigmaPhoton->SetPairs(pPhotonSigma, ApAPhotonSigma);
  CF_SigmaPhoton->GetCorrelations();
  CF_SigmaPhoton->WriteOutput(
      Form("%s/CFOutput_pPhoton.root", foldername.Data()));
}

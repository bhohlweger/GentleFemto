#include <iostream>
#include "ReadDreamFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TROOT.h"

void GetSigmaCorrelations(const char* filename, const char* trigger,
                          const char* suffixChar) {
  TString appendix = TString::Format("%s", trigger);
  TString suffix = TString::Format("%s", suffixChar);
  const int nPairs = (suffix == "0") ? 15 : 8;
  ReadDreamFile* DreamFile = new ReadDreamFile(nPairs, nPairs);
  DreamFile->SetAnalysisFile(filename, appendix.Data(), suffix.Data());

  DreamCF* CF_pp = new DreamCF();
  DreamPair* pp = new DreamPair("Part", 0.24, 0.34);
  DreamPair* ApAp = new DreamPair("AntiPart", 0.24, 0.34);

  DreamCF* CF_pSigma = new DreamCF();
  DreamPair* pSigma = new DreamPair("Part", 0.25, 0.4);
  DreamPair* ApASigma = new DreamPair("AntiPart", 0.25, 0.4);

  DreamCF* CF_SidebandUp = new DreamCF();
  DreamPair* pSiSBup = new DreamPair("Part_SB_up", 0.25, 0.4);
  DreamPair* ApaSiSBup = new DreamPair("AntiPart_SB_up", 0.25, 0.4);

  DreamCF* CF_SidebandLow = new DreamCF();
  DreamPair* pSiSBlow = new DreamPair("Part_SB_low", 0.25, 0.4);
  DreamPair* ApaSiSBlow = new DreamPair("AntiPart_SB_low", 0.25, 0.4);

  DreamCF* CF_pSigmaLambda = new DreamCF();
  DreamPair* pLambdaSigma = new DreamPair("PartDaughter", 0.25, 0.4);
  DreamPair* ApALambdaSigma = new DreamPair("AntiPartDaughter", 0.25, 0.4);

  DreamCF* CF_pLambda = new DreamCF();
  DreamPair* pLambda = new DreamPair("PartDaughter", 0.25, 0.4);
  DreamPair* ApALambda = new DreamPair("AntiPartDaughter", 0.25, 0.4);

  DreamCF* CF_pSigmaPhoton = new DreamCF();
  DreamPair* pPhotonSigma = new DreamPair("PartDaughter", 0.25, 0.4);
  DreamPair* ApPhotonSigma = new DreamPair("AntiPartDaughter", 0.25, 0.4);

  DreamCF* CF_pPhoton= new DreamCF();
  DreamPair* pPhoton = new DreamPair("PartDaughter", 0.25, 0.4);
  DreamPair* ApPhoton= new DreamPair("AntiPartDaughter", 0.25, 0.4);

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

  if (suffix == "0") {
    pLambdaSigma->SetPair(DreamFile->GetPairDistributions(0, 8, ""));
    ApALambdaSigma->SetPair(DreamFile->GetPairDistributions(1, 9, ""));

    pLambda->SetPair(DreamFile->GetPairDistributions(0, 10, ""));
    ApALambda->SetPair(DreamFile->GetPairDistributions(1, 11, ""));

    pPhotonSigma->SetPair(DreamFile->GetPairDistributions(0, 12, ""));
    ApPhotonSigma->SetPair(DreamFile->GetPairDistributions(1, 13, ""));

    pPhoton->SetPair(DreamFile->GetPairDistributions(0, 14, ""));
    ApPhoton->SetPair(DreamFile->GetPairDistributions(1, 14, ""));
  }

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

  if (suffix == "0") {
    pLambdaSigma->ShiftForEmpty(pLambdaSigma->GetPair());
    ApALambdaSigma->ShiftForEmpty(ApALambdaSigma->GetPair());

    pLambda->ShiftForEmpty(pLambda->GetPair());
    ApALambda->ShiftForEmpty(ApALambda->GetPair());

    pPhotonSigma->ShiftForEmpty(pPhotonSigma->GetPair());
    ApPhotonSigma->ShiftForEmpty(ApPhotonSigma->GetPair());

    pPhoton->ShiftForEmpty(pPhoton->GetPair());
    ApPhoton->ShiftForEmpty(ApPhoton->GetPair());
  }

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

  if (suffix == "0") {
    pLambdaSigma->FixShift(pLambdaSigma->GetPairShiftedEmpty(0),
                           ApALambdaSigma->GetPairShiftedEmpty(0),
                           ApALambdaSigma->GetFirstBin());

    ApALambdaSigma->FixShift(ApALambdaSigma->GetPairShiftedEmpty(0),
                             pLambdaSigma->GetPairShiftedEmpty(0),
                             pLambdaSigma->GetFirstBin());

    pLambda->FixShift(pLambda->GetPairShiftedEmpty(0),
                      ApALambda->GetPairShiftedEmpty(0),
                      ApALambda->GetFirstBin());

    ApALambda->FixShift(ApALambda->GetPairShiftedEmpty(0),
                        pLambda->GetPairShiftedEmpty(0),
                        pLambda->GetFirstBin());

    pPhotonSigma->FixShift(pPhotonSigma->GetPairShiftedEmpty(0),
                           ApPhotonSigma->GetPairShiftedEmpty(0),
                           ApPhotonSigma->GetFirstBin());

    ApPhotonSigma->FixShift(ApPhotonSigma->GetPairShiftedEmpty(0),
                            pPhotonSigma->GetPairShiftedEmpty(0),
                            pPhotonSigma->GetFirstBin());

    pPhoton->FixShift(pPhoton->GetPairShiftedEmpty(0),
                      ApPhoton->GetPairShiftedEmpty(0),
                      ApPhoton->GetFirstBin());

    ApPhoton->FixShift(ApPhoton->GetPairShiftedEmpty(0),
                       pPhoton->GetPairShiftedEmpty(0), pPhoton->GetFirstBin());
  }

  std::cout << "=========================" << std::endl;
  std::cout << "==Rebinning & Weighting==" << std::endl;
  std::cout << "=========================" << std::endl;

  std::vector<int> rebin = { { 1, 4, 5, 10 } };

  for (size_t iReb = 0; iReb < rebin.size(); ++iReb) {
    std::cout << "==Rebinning==" << std::endl;
    pSigma->Rebin(pSigma->GetPairFixShifted(0), rebin[iReb], true);
    ApASigma->Rebin(ApASigma->GetPairFixShifted(0), rebin[iReb], true);
    pSiSBup->Rebin(pSiSBup->GetPairFixShifted(0), rebin[iReb], true);
    ApaSiSBup->Rebin(ApaSiSBup->GetPairFixShifted(0), rebin[iReb], true);
    pSiSBlow->Rebin(pSiSBlow->GetPairFixShifted(0), rebin[iReb], true);
    ApaSiSBlow->Rebin(ApaSiSBlow->GetPairFixShifted(0), rebin[iReb], true);
    if (suffix == "0") {
      pLambdaSigma->Rebin(pLambdaSigma->GetPairFixShifted(0), rebin[iReb], true);
      ApALambdaSigma->Rebin(ApALambdaSigma->GetPairFixShifted(0), rebin[iReb], true);
      pLambda->Rebin(pLambda->GetPairFixShifted(0), rebin[iReb], true);
      ApALambda->Rebin(ApALambda->GetPairFixShifted(0), rebin[iReb], true);
      pPhotonSigma->Rebin(pPhotonSigma->GetPairFixShifted(0), rebin[iReb], true);
      ApPhotonSigma->Rebin(ApPhotonSigma->GetPairFixShifted(0), rebin[iReb], true);
      pPhoton->Rebin(pPhoton->GetPairFixShifted(0), rebin[iReb], true);
      ApPhoton->Rebin(ApPhoton->GetPairFixShifted(0), rebin[iReb], true);
    }
    std::cout << "==Weighting==" << std::endl;
    pSigma->ReweightMixedEvent(pSigma->GetPairRebinned(iReb), 0.2, 0.9, pSigma->GetPair());
    ApASigma->ReweightMixedEvent(ApASigma->GetPairRebinned(iReb), 0.2, 0.9, ApASigma->GetPair());
    pSiSBup->ReweightMixedEvent(pSiSBup->GetPairRebinned(iReb), 0.2, 0.9, pSiSBup->GetPair());
    ApaSiSBup->ReweightMixedEvent(ApaSiSBup->GetPairRebinned(iReb), 0.2, 0.9, ApaSiSBup->GetPair());
    pSiSBlow->ReweightMixedEvent(pSiSBlow->GetPairRebinned(iReb), 0.2, 0.9, pSiSBlow->GetPair());
    ApaSiSBlow->ReweightMixedEvent(ApaSiSBlow->GetPairRebinned(iReb), 0.2, 0.9, ApaSiSBlow->GetPair());
    if (suffix == "0") {
      pLambdaSigma->ReweightMixedEvent(pLambdaSigma->GetPairRebinned(iReb), 0.2,
                                       0.9, pLambdaSigma->GetPair());
      ApALambdaSigma->ReweightMixedEvent(ApALambdaSigma->GetPairRebinned(iReb),
                                         0.2, 0.9, ApALambdaSigma->GetPair());
      pLambda->ReweightMixedEvent(pLambda->GetPairRebinned(iReb), 0.2, 0.9, pLambda->GetPair());
      ApALambda->ReweightMixedEvent(ApALambda->GetPairRebinned(iReb), 0.2, 0.9, ApALambda->GetPair());
      pPhotonSigma->ReweightMixedEvent(pPhotonSigma->GetPairRebinned(iReb), 0.2,
                                       0.9, pPhotonSigma->GetPair());
      ApPhotonSigma->ReweightMixedEvent(ApPhotonSigma->GetPairRebinned(iReb),
                                        0.2, 0.9, ApPhotonSigma->GetPair());
      pPhoton->ReweightMixedEvent(pPhoton->GetPairRebinned(iReb), 0.2, 0.9, pPhoton->GetPair());
      ApPhoton->ReweightMixedEvent(ApPhoton->GetPairRebinned(iReb), 0.2, 0.9, ApPhoton->GetPair());
    }
  }
  pp->ReweightMixedEvent(pp->GetPairFixShifted(0), 0.2, 0.9, pp->GetPair());
  ApAp->ReweightMixedEvent(ApAp->GetPairFixShifted(0), 0.2, 0.9, ApAp->GetPair());

  pSigma->Rebin(pSigma->GetPair(), 4, true);
  pSigma->Rebin(pSigma->GetPair(), 5, true);
  ApASigma->Rebin(ApASigma->GetPair(), 4, true);
  ApASigma->Rebin(ApASigma->GetPair(), 5, true);
  pSiSBup->Rebin(pSiSBup->GetPair(), 4, true);
  pSiSBup->Rebin(pSiSBup->GetPair(), 5, true);
  ApaSiSBup->Rebin(ApaSiSBup->GetPair(), 4, true);
  ApaSiSBup->Rebin(ApaSiSBup->GetPair(), 5, true);
  pSiSBlow->Rebin(pSiSBlow->GetPair(), 4, true);
  pSiSBlow->Rebin(pSiSBlow->GetPair(), 5, true);
  ApaSiSBlow->Rebin(ApaSiSBlow->GetPair(), 4, true);
  ApaSiSBlow->Rebin(ApaSiSBlow->GetPair(), 5, true);
  if (suffix == "0") {
    pLambdaSigma->Rebin(pLambdaSigma->GetPair(), 4, true);
    pLambdaSigma->Rebin(pLambdaSigma->GetPair(), 5, true);
    ApALambdaSigma->Rebin(ApALambdaSigma->GetPair(), 4, true);
    ApALambdaSigma->Rebin(ApALambdaSigma->GetPair(), 5, true);
    pLambda->Rebin(pLambda->GetPair(), 4, true);
    pLambda->Rebin(pLambda->GetPair(), 5, true);
    ApALambda->Rebin(ApALambda->GetPair(), 4, true);
    ApALambda->Rebin(ApALambda->GetPair(), 5, true);

    pPhotonSigma->Rebin(pPhotonSigma->GetPair(), 4, true);
    pPhotonSigma->Rebin(pPhotonSigma->GetPair(), 5, true);
    ApPhotonSigma->Rebin(ApPhotonSigma->GetPair(), 4, true);
    ApPhotonSigma->Rebin(ApPhotonSigma->GetPair(), 5, true);
    pPhoton->Rebin(pPhoton->GetPair(), 4, true);
    pPhoton->Rebin(pPhoton->GetPair(), 5, true);
    ApPhoton->Rebin(ApPhoton->GetPair(), 4, true);
    ApPhoton->Rebin(ApPhoton->GetPair(), 5, true);
  }

  std::cout << "=========================" << std::endl;
  std::cout << "=========CFs=============" << std::endl;
  std::cout << "=========================" << std::endl;

  TString foldername = filename;
  foldername.ReplaceAll("AnalysisResults.root", "");

  TString fileAppendix =
      (suffix == "0") ? "" : TString::Format("_%s", suffix.Data());

  CF_pp->SetPairs(pp, ApAp);
  CF_pp->GetCorrelations();
  CF_pp->WriteOutput(
      Form("%s/CFOutput_pp%s.root", foldername.Data(), fileAppendix.Data()));

  CF_pSigma->SetPairs(pSigma, ApASigma);
  CF_pSigma->GetCorrelations();
  CF_pSigma->WriteOutput(
      Form("%s/CFOutput_pSigma%s.root", foldername.Data(),
           fileAppendix.Data()));

  CF_SidebandUp->SetPairs(pSiSBup, ApaSiSBup);
  CF_SidebandUp->GetCorrelations();
  CF_SidebandUp->WriteOutput(
      Form("%s/CFOutput_SB_up%s.root", foldername.Data(), fileAppendix.Data()));

  CF_SidebandLow->SetPairs(pSiSBlow, ApaSiSBlow);
  CF_SidebandLow->GetCorrelations();
  CF_SidebandLow->WriteOutput(
      Form("%s/CFOutput_SB_low%s.root", foldername.Data(),
           fileAppendix.Data()));

  if (suffix == "0") {
    CF_pSigmaLambda->SetPairs(pLambdaSigma, ApALambdaSigma);
    CF_pSigmaLambda->GetCorrelations();
    CF_pSigmaLambda->WriteOutput(
        Form("%s/CFOutput_pLambdaSigma.root", foldername.Data()));

    CF_pLambda->SetPairs(pLambda, ApALambda);
    CF_pLambda->GetCorrelations();
    CF_pLambda->WriteOutput(
        Form("%s/CFOutput_pLambda.root", foldername.Data()));

    CF_pSigmaPhoton->SetPairs(pPhotonSigma, ApPhotonSigma);
    CF_pSigmaPhoton->GetCorrelations();
    CF_pSigmaPhoton->WriteOutput(
        Form("%s/CFOutput_pPhotonSigma.root", foldername.Data()));

    CF_pPhoton->SetPairs(pPhoton, ApPhoton);
    CF_pPhoton->GetCorrelations();
    CF_pPhoton->WriteOutput(
        Form("%s/CFOutput_pPhoton.root", foldername.Data()));
  }
}

#include <iostream>
#include "ReadDreamFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TROOT.h"

void GetSigmaCorrelations(const char* filename, const char* trigger,
                          const char* suffixChar) {
  TString appendix = TString::Format("%s", trigger);
  TString suffix = TString::Format("%s", suffixChar);
  ReadDreamFile* DreamFile = new ReadDreamFile(8, 8);
  DreamFile->SetAnalysisFile(filename, appendix.Data(), suffix.Data());

  DreamCF* CF_pp = new DreamCF();
  DreamPair* pp = new DreamPair("Part", 0.24, 0.34);
  DreamPair* ApAp = new DreamPair("AntiPart", 0.24, 0.34);

  DreamCF* CF_pSigma = new DreamCF();
  DreamPair* pSigma = new DreamPair("Part", 0.25, 0.4);
  DreamPair* ApASigma = new DreamPair("AntiPart", 0.25, 0.4);

  DreamCF* CF_SidebandUp = new DreamCF();
  DreamPair* pSiSBup = new DreamPair("Part_SB_up", 0.3, 0.5);
  DreamPair* ApaSiSBup = new DreamPair("AntiPart_SB_up", 0.3, 0.5);

  DreamCF* CF_SidebandLow = new DreamCF();
  DreamPair* pSiSBlow = new DreamPair("Part_SB_low", 0.3, 0.5);
  DreamPair* ApaSiSBlow = new DreamPair("AntiPart_SB_low", 0.3, 0.5);

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
  std::cout << "=========================" << std::endl;
  std::cout << "==Rebinning & Weighting==" << std::endl;
  std::cout << "=========================" << std::endl;

  std::vector<int> rebin = { { 1, 4, 5, 10 } };

  for (size_t iReb = 0; iReb < rebin.size(); ++iReb) {
    std::cout << "==Rebinning==" << std::endl;
    pSigma->Rebin(pSigma->GetPairFixShifted(0), rebin[iReb]);
    ApASigma->Rebin(ApASigma->GetPairFixShifted(0), rebin[iReb]);
    pSiSBup->Rebin(pSiSBup->GetPairFixShifted(0), rebin[iReb]);
    ApaSiSBup->Rebin(ApaSiSBup->GetPairFixShifted(0), rebin[iReb]);
    pSiSBlow->Rebin(pSiSBlow->GetPairFixShifted(0), rebin[iReb]);
    ApaSiSBlow->Rebin(ApaSiSBlow->GetPairFixShifted(0), rebin[iReb]);
    std::cout << "==Weighting==" << std::endl;
    pSigma->ReweightMixedEvent(pSigma->GetPairRebinned(iReb), 0.2, 0.9);
    ApASigma->ReweightMixedEvent(ApASigma->GetPairRebinned(iReb), 0.2, 0.9);
    pSiSBup->ReweightMixedEvent(pSiSBup->GetPairRebinned(iReb), 0.2, 0.9);
    ApaSiSBup->ReweightMixedEvent(ApaSiSBup->GetPairRebinned(iReb), 0.2, 0.9);
    pSiSBlow->ReweightMixedEvent(pSiSBlow->GetPairRebinned(iReb), 0.2, 0.9);
    ApaSiSBlow->ReweightMixedEvent(ApaSiSBlow->GetPairRebinned(iReb), 0.2, 0.9);
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
}

#include "TROOT.h"
#include "ReadDreamFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <iostream>

void GetCorrelationsPd(const char* filename, const char* Path, const char* prefix,
                        const char* addon = "") {
  //gStyle->SetOptStat(0);
  ReadDreamFile* DreamFile = new ReadDreamFile(4,4);
  DreamFile->SetAnalysisFile(filename, prefix, addon);

  DreamCF* CF_pp = new DreamCF();
  DreamPair* pp = new DreamPair("Part", 0.24, 0.34);
  DreamPair* ApAp = new DreamPair("AntiPart", 0.24, 0.34);

  DreamCF* CF_pd = new DreamCF();
  DreamPair* pd = new DreamPair("Part", 0.24, 0.34);
  DreamPair* ApAd = new DreamPair("AntiPart", 0.24, 0.34);

  std::cout << "=========================" << std::endl;
  std::cout << "========Pair Set=========" << std::endl;
  std::cout << "=========================" << std::endl;
  pp->SetPair(DreamFile->GetPairDistributions(0, 0, ""));
  ApAp->SetPair(DreamFile->GetPairDistributions(1, 1, ""));

  pd->SetPair(DreamFile->GetPairDistributions(0, 2, ""));
  ApAd->SetPair(DreamFile->GetPairDistributions(1, 3, ""));

  std::cout << "=========================" << std::endl;
  std::cout << "======Pair Shifted=======" << std::endl;
  std::cout << "=========================" << std::endl;
  pp->ShiftForEmpty(pp->GetPair());
  ApAp->ShiftForEmpty(ApAp->GetPair());

  pd->ShiftForEmpty(pd->GetPair());
  ApAd->ShiftForEmpty(ApAd->GetPair());

  std::cout << "=========================" << std::endl;
  std::cout << "====Pair Fix Shifted=====" << std::endl;
  std::cout << "=========================" << std::endl;
  pp->FixShift(pp->GetPairShiftedEmpty(0), ApAp->GetPairShiftedEmpty(0),
               ApAp->GetFirstBin());
  ApAp->FixShift(ApAp->GetPairShiftedEmpty(0), pp->GetPairShiftedEmpty(0),
                 pp->GetFirstBin());
  pp->FixShift(pp->GetPair(), ApAp->GetPair(), 0.004, true);
  ApAp->FixShift(ApAp->GetPair(), pp->GetPair(), 0.004, true);
  pd->FixShift(pd->GetPairShiftedEmpty(0), ApAd->GetPairShiftedEmpty(0),
                 ApAd->GetFirstBin());
  ApAd->FixShift(ApAd->GetPairShiftedEmpty(0), pd->GetPairShiftedEmpty(0),
                  pd->GetFirstBin());
  std::cout << "=========================" << std::endl;
  std::cout << "==Rebinning & Weighting==" << std::endl;
  std::cout << "=========================" << std::endl;

  std::vector<int> rebinVec = { { 4, 5 } };
  for (size_t iReb = 0; iReb < rebinVec.size(); ++iReb) {
    pd->Rebin(pd->GetPairFixShifted(0), rebinVec[iReb]);
    pd->ReweightMixedEvent(pd->GetPairRebinned(iReb), 0.2, 0.9);
    ApAd->Rebin(ApAd->GetPairFixShifted(0), rebinVec[iReb]);
    ApAd->ReweightMixedEvent(ApAd->GetPairRebinned(iReb), 0.2, 0.9);
  }
  pp->ReweightMixedEvent(pp->GetPairFixShifted(0), 0.2, 0.9);
  ApAp->ReweightMixedEvent(ApAp->GetPairFixShifted(0), 0.2, 0.9);

  pp->ReweightMixedEvent(pp->GetPairFixShifted(1), 0.2, 0.9);
  ApAp->ReweightMixedEvent(ApAp->GetPairFixShifted(1), 0.2, 0.9);

  std::cout << "=========================" << std::endl;
  std::cout << "=========CFs=============" << std::endl;
  std::cout << "=========================" << std::endl;

  TString foldername = filename;
  foldername.ReplaceAll("AnalysisResults.root", "");

  CF_pp->SetPairs(pp, ApAp);
  CF_pp->GetCorrelations();
  CF_pp->WriteOutput(Form("%sCFOutput_pp.root", foldername.Data()));

  CF_pd->SetPairs(pd, ApAd);
  CF_pd->GetCorrelations();
  CF_pd->WriteOutput(Form("%sCFOutput_pd.root", foldername.Data()));
}

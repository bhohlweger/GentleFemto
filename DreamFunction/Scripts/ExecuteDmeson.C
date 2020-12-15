#include "ReadDreamFile.h"

int main(int argc, char* argv[]) {
  const char* filename = argv[1];
  TString appendix = TString::Format("%s", argv[2]);
  TString suffix = TString::Format("%s", argv[3]);
  ReadDreamFile* DreamFile = new ReadDreamFile(4, 4);
  DreamFile->SetAnalysisFile(filename, appendix.Data(), suffix.Data());

  const double normLower = 0.8;
  const double normUpper = 1.;

  DreamCF* CF_pDplus = new DreamCF();
  DreamPair* pDplus = new DreamPair("Part", normLower, normUpper);
  DreamPair* apDminus = new DreamPair("AntiPart", normLower, normUpper);

  DreamCF* CF_pDminus = new DreamCF();
  DreamPair* pDminus = new DreamPair("Part", normLower, normUpper);
  DreamPair* apDplus = new DreamPair("AntiPart", normLower, normUpper);

  pDplus->SetPair(DreamFile->GetPairDistributions(0, 2, ""));
  apDminus->SetPair(DreamFile->GetPairDistributions(1, 3, ""));
  pDminus->SetPair(DreamFile->GetPairDistributions(0, 3, ""));
  apDplus->SetPair(DreamFile->GetPairDistributions(1, 2, ""));

  pDplus->ShiftForEmpty(pDplus->GetPair());
  apDminus->ShiftForEmpty(apDminus->GetPair());
  pDminus->ShiftForEmpty(pDminus->GetPair());
  apDplus->ShiftForEmpty(apDplus->GetPair());

  pDplus->FixShift(pDplus->GetPairShiftedEmpty(0),
                   apDminus->GetPairShiftedEmpty(0), apDminus->GetFirstBin());
  apDminus->FixShift(apDminus->GetPairShiftedEmpty(0),
                     pDplus->GetPairShiftedEmpty(0), pDplus->GetFirstBin());
  pDminus->FixShift(pDminus->GetPairShiftedEmpty(0),
                    apDplus->GetPairShiftedEmpty(0), apDplus->GetFirstBin());
  apDplus->FixShift(apDplus->GetPairShiftedEmpty(0),
                    pDminus->GetPairShiftedEmpty(0), pDminus->GetFirstBin());

  std::vector<int> rebin = { { 1, 4, 5, 10, 20 } };

  for (size_t iReb = 0; iReb < rebin.size(); ++iReb) {
    pDplus->Rebin(pDplus->GetPairFixShifted(0), rebin[iReb], true);
    apDminus->Rebin(apDminus->GetPairFixShifted(0), rebin[iReb], true);
    pDminus->Rebin(pDminus->GetPairFixShifted(0), rebin[iReb], true);
    apDplus->Rebin(apDplus->GetPairFixShifted(0), rebin[iReb], true);

    pDplus->ReweightMixedEvent(pDplus->GetPairRebinned(iReb), 0.2, 0.9,
                               pDplus->GetPair());
    apDminus->ReweightMixedEvent(apDminus->GetPairRebinned(iReb), 0.2, 0.9,
                                 apDminus->GetPair());
    pDminus->ReweightMixedEvent(pDminus->GetPairRebinned(iReb), 0.2, 0.9,
                                pDminus->GetPair());
    apDplus->ReweightMixedEvent(apDplus->GetPairRebinned(iReb), 0.2, 0.9,
                                apDplus->GetPair());
  }


  TString foldername = filename;
  foldername.ReplaceAll("AnalysisResults.root", "");

  TString fileAppendix =
      (suffix == "0") ? "" : TString::Format("_%s", suffix.Data());

  CF_pDplus->SetPairs(pDplus, apDminus);
  CF_pDplus->GetCorrelations();
  CF_pDplus->WriteOutput(
      Form("%s/CFOutput_pDplus%s.root", foldername.Data(), fileAppendix.Data()));

  CF_pDminus->SetPairs(pDminus, apDplus);
  CF_pDminus->GetCorrelations();
  CF_pDminus->WriteOutput(
      Form("%s/CFOutput_pDminus%s.root", foldername.Data(),
           fileAppendix.Data()));

  return 1;
}

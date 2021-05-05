#include "ReadDreamFile.h"
#include "TLegend.h"
#include "METoSEReweighting.C"

int main(int argc, char* argv[]) {
  const char* filename = argv[1];
  TString appendix = TString::Format("%s", argv[2]);
  TString suffix = TString::Format("%s", argv[3]);
  ReadDreamFile* DreamFile = new ReadDreamFile(4, 4);
  DreamFile->SetAnalysisFile(filename, appendix.Data(), suffix.Data());

  const double normLower = 1.5;
  const double normUpper = 2.;

  DreamCF* CF_pDminus = new DreamCF();
  DreamPair* pDminus = new DreamPair("Part", normLower, normUpper);
  DreamPair* apDplus = new DreamPair("AntiPart", normLower, normUpper);

  pDminus->SetPair(DreamFile->GetPairDistributions(0, 3, ""));
  apDplus->SetPair(DreamFile->GetPairDistributions(1, 2, ""));

  pDminus->ShiftForEmpty(pDminus->GetPair());
  apDplus->ShiftForEmpty(apDplus->GetPair());

  pDminus->FixShift(pDminus->GetPairShiftedEmpty(0),
                    apDplus->GetPairShiftedEmpty(0), apDplus->GetFirstBin());
  apDplus->FixShift(apDplus->GetPairShiftedEmpty(0),
                    pDminus->GetPairShiftedEmpty(0), pDminus->GetFirstBin());

  std::vector<int> rebin = { { 1, 4, 5, 10, 20 } };

  for (size_t iReb = 0; iReb < rebin.size(); ++iReb) {
    pDminus->Rebin(pDminus->GetPairFixShifted(0), rebin[iReb], true);
    apDplus->Rebin(apDplus->GetPairFixShifted(0), rebin[iReb], true);

    pDminus->ReweightMixedEvent(pDminus->GetPairRebinned(iReb), 0.2, 0.9,
                                pDminus->GetPair());
    apDplus->ReweightMixedEvent(apDplus->GetPairRebinned(iReb), 0.2, 0.9,
                                apDplus->GetPair());
  }

  TString foldername = filename;
  foldername.ReplaceAll("AnalysisResults.root", "");

  TString fileAppendix =
      (suffix == "0") ? "" : TString::Format("_%s", suffix.Data());

  CF_pDminus->SetPairs(pDminus, apDplus);
  CF_pDminus->GetCorrelations();
  CF_pDminus->WriteOutput(
      Form("%s/CFOutput_pDminus%s.root", foldername.Data(),
           fileAppendix.Data()));

  TString FileName = Form("%s/CFOutput_pDminus%s.root", foldername.Data(),
                          fileAppendix.Data());
  TFile* file = TFile::Open(FileName, "update");
  TList* PairDist = (TList*) file->Get("PairDist");
  if (PairDist) {
    ReweightingQA(PairDist);
  } else {
    file->ls();
  }
  TList* AntiPairDist = (TList*) file->Get("AntiPairDist");
  if (AntiPairDist) {
    ReweightingQA(AntiPairDist);
  }
  return 1;
}

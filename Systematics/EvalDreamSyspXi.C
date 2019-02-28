#include "CATSInput.h"
#include "DreamPlot.h"
#include "DreamSystematics.h"
#include "ForgivingReader.h"
#include "CandidateCounter.h"
#include "TCanvas.h"
#include <iostream>

void EvalDreamSystematics(TString InputDir, TString prefix) {
  TString filename = Form("%s/AnalysisResults.root", InputDir.Data());
  DreamPlot::SetStyle();
  auto CATSinput = new CATSInput();
  CATSinput->SetNormalization(0.240, 0.340);
  CATSinput->SetFixedkStarMinBin(true, 0.008);

  auto counter = new CandidateCounter();

  ReadDreamFile* DreamFile = new ReadDreamFile(6, 6);
  DreamFile->SetAnalysisFile(filename.Data(), prefix);

  //Proton - Proton
  DreamSystematics protonproton(DreamSystematics::pp);
  const int protonVarStart = 1;

  //Proton - Xi

  DreamDist* pXi = DreamFile->GetPairDistributions(0, 4, "");
  DreamDist* ApAXi = DreamFile->GetPairDistributions(1, 5, "");
  DreamCF* CFpXiDef = CATSinput->ObtainCFSyst(20, "ppDef", pXi, ApAXi);
  DreamSystematics protonXi(DreamSystematics::pXi);
  protonXi.SetDefaultPair(CFpXiDef, "hCk_RebinnedppDefMeV_0");
  unsigned int nProtons = 1.799336e8+1.545516e8;
  unsigned int nCascade = 838588 + 792439;
  protonXi.SetNDefaultParticles(nProtons,nCascade);
  int iPXICounter = 0;
  const int XiVarStart = protonproton.GetNumberOfVars();
  for (int i = protonVarStart;
      i < protonVarStart + protonproton.GetNumberOfVars(); ++i) {
    ReadDreamFile* DreamVarFile = new ReadDreamFile(6, 6);
    DreamVarFile->SetAnalysisFile(filename.Data(), prefix, Form("%u", i));
    TString VarName = TString::Format("pXiVar%u", i);
    DreamCF* CFpXiVar = CATSinput->ObtainCFSyst(
        20, VarName.Data(), DreamVarFile->GetPairDistributions(0, 4, ""),
        DreamVarFile->GetPairDistributions(1, 5, ""), pXi, ApAXi);
    protonXi.SetVarPair(CFpXiVar,
                        TString::Format("Rebinned%sMeV", VarName.Data()));
    TString VarString = TString::Format("%u", i);
    ForgivingReader* ForgivingFile = new ForgivingReader(filename.Data(),
                                                         prefix,
                                                         VarString.Data());
    counter->SetNumberOfCandidates(ForgivingFile);
    protonXi.SetNVariationParticles(counter->GetNumberOfTracks(),
                                    counter->GetNumberOfCascades());
    iPXICounter++;
  }
  for (int i = XiVarStart;
      i
          < XiVarStart + protonXi.GetNumberOfVars()
              - protonproton.GetNumberOfVars(); ++i) {
    ReadDreamFile* DreamVarFile = new ReadDreamFile(6, 6);
    DreamVarFile->SetAnalysisFile(filename.Data(), prefix, Form("%u", i));
    TString VarName = TString::Format("pXiVar%u", i);
    DreamCF* CFpXiVar = CATSinput->ObtainCFSyst(
        20, VarName.Data(), DreamVarFile->GetPairDistributions(0, 4, ""),
        DreamVarFile->GetPairDistributions(1, 5, ""), pXi, ApAXi);
    protonXi.SetVarPair(CFpXiVar,
                        TString::Format("Rebinned%sMeV", VarName.Data()));
    TString VarString = TString::Format("%u", i);
    ForgivingReader* ForgivingFile = new ForgivingReader(filename.Data(),
                                                         prefix,
                                                         VarString.Data());
    counter->SetNumberOfCandidates(ForgivingFile);
    protonXi.SetNVariationParticles(counter->GetNumberOfTracks(),
                                    counter->GetNumberOfCascades());
    iPXICounter++;
  }
  protonXi.EvalSystematics();
  protonXi.EvalDifferenceInPairs();
  protonXi.EvalDifferenceInParticles();
  protonXi.WriteOutput();
  protonXi.DrawAllCF();
  std::cout << "Worked through " << iPXICounter << " variations" << std::endl;
}

int main(int argc, char* argv[]) {
  EvalDreamSystematics(argv[1], argv[2]);

  return 1;
}

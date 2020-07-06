#include "DreamKayTee.h"
#include "ReadDreamFile.h"
#include "DreamSystematics.h"
#include "TROOT.h"
#include "TSystem.h"
#include "CATSInput.h"
#include "CandidateCounter.h"
#include "ForgivingReader.h"

int main(int argc, char* argv[]) {

//to get systematics in each mT bin
TString InputFolderVar = "/Users/sartozza/cernbox/Analysis/BBbar/GentleFemto_Output/NanoOutput/Systematics_Nano/Syst_mT/LAL_pair/3_bins/";
TString InputFolderDef = "/Users/sartozza/cernbox/Analysis/BBbar/GentleFemto_Output/NanoOutput/mT_Analysis/Data/LAL_pair/mT_3bins/";
TFile* FileVar;
TFile* FileDef;
TH1F* histoDef;
TH1F* histoVar;

double upperFitRange = 400;
std::vector<int> mt_bins = {0,1,2};
for (int mbin=0;mbin<mt_bins.size(); mbin++){
  FileDef = new TFile(InputFolderDef+TString::Format("CFOutput_mT_LALDef_HMBBar_%i.root",mt_bins[mbin]),"read");
  TH1F* ClonehistoDef = (TH1F*)FileDef->Get("hCk_RebinnedLALDefMeV_0")->Clone(TString::Format("LALDef_mT%i", mt_bins[mbin]));

  DreamSystematics protonAL(DreamSystematics::LAL);
	if(ClonehistoDef) protonAL.SetDefaultHist(ClonehistoDef);
  else if (!ClonehistoDef) printf("Warning no existing Def Histo\n");
  protonAL.SetUpperFitRange(upperFitRange);

  for(int iVar=1;iVar<=44;iVar++){
    FileVar = new TFile(InputFolderVar+TString::Format("CFOutput_mT_LALVar%i_HMBBar_%i.root",iVar,mt_bins[mbin]),"read");
    TH1F* ClonehistoVar = (TH1F*)FileVar->Get(TString::Format("hCk_RebinnedLALVar%iMeV_0",iVar))->Clone(TString::Format("LALVar%i_mT%i", iVar,mt_bins[mbin]));
	  if(ClonehistoVar) protonAL.SetVarHist(ClonehistoVar);
    else if (!ClonehistoVar) printf("Warning no existing Def Histo\n");

  }
  protonAL.EvalSystematicsBBar(0);
  protonAL.WriteOutput(TString::Format("mt%i",mt_bins[mbin]));
  FileDef->Close();
  FileVar->Close();

}

}


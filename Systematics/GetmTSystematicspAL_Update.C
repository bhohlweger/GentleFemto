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
TString InputFolderVar = "/Users/sartozza/cernbox/Analysis/BBbar/GentleFemto_Output/NanoOutput/Systematics_Nano/Syst_mT/pAL_pair/";
TString InputFolderDef = "/Users/sartozza/cernbox/Analysis/BBbar/GentleFemto_Output/NanoOutput/mT_Analysis/Data/pAL_pair/";
TFile* FileVar;
TFile* FileDef;
TH1F* histoDef;
TH1F* histoVar;

double upperFitRange = 400;
std::vector<int> mt_bins = {0,1,2,3,4,5};
for (int mbin=0;mbin<mt_bins.size(); mbin++){
  printf("---mt Bin = %i\n",mt_bins[mbin]);
  FileDef = new TFile(InputFolderDef+TString::Format("CFOutput_mT_pALDef_HMBBar_%i.root",mt_bins[mbin]),"read");
  TH1F* ClonehistoDef = (TH1F*)FileDef->Get("hCk_RebinnedpALDefMeV_0")->Clone(TString::Format("pALDef_mT%i", mt_bins[mbin]));

  DreamSystematics protonAL(DreamSystematics::pAL);
	if(ClonehistoDef) protonAL.SetDefaultHist(ClonehistoDef);
  else if (!ClonehistoDef) printf("Warning no existing Def Histo\n");
  protonAL.SetUpperFitRange(upperFitRange);
  // FileDef->Close();

  for(int iVar=1;iVar<=42;iVar++){
    FileVar = new TFile(InputFolderVar+TString::Format("CFOutput_mT_pALVar%i_HMBBar_%i.root",iVar,mt_bins[mbin]),"read");
    TH1F* ClonehistoVar = (TH1F*)FileVar->Get(TString::Format("hCk_RebinnedpALVar%iMeV_0",iVar))->Clone(TString::Format("pALVar%i_mT%i", iVar,mt_bins[mbin]));
	  if(ClonehistoVar) protonAL.SetVarHist(ClonehistoVar);
    else if (!ClonehistoVar) printf("Warning no existing Def Histo\n");
    // FileVar->Close();
  }

  protonAL.EvalSystematicsBBar(0);
  protonAL.WriteOutput(TString::Format("mt%i",mt_bins[mbin]));
  FileDef->Close();
  FileVar->Close();

}

}


#include "TROOT.h"
#include "ReadDreamFile.h"
#include "DreamPlot.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <iostream>
#include "TPad.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TSystem.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TApplication.h"
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// void BbarB_MiniJet(const char* filename, const char* filenameMC) {
void BbarB_MiniJet(const char* STflag = "", const bool writefile = false) {
DreamPlot::SetStyle();

TString add1="1";
TString add2="2";
TString add3="3";
TString add4="4";
TString add5="5";

double maxfit = 0.4;

//Initializing the Lists

TList *listTrackCuts=0;
TList *listv0Cuts=0;
TList *listAntiTrackCuts=0;
TList *listAntiv0Cuts=0;

//Accessing all the directory in the root file
  TFile* file0data[4];
  TFile* file0MC[4];

  bool whichfolder = true;// true ->CorrectBin folder

  if(!whichfolder){
  if(strcmp(STflag, add1)==0)
  {
  std::cout << "Sphericity [0.,0.3] \n";
  file0data[0]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CFOutput_pAp_App_st1.root");
  file0MC[0]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/MCCFOutput_pAp_App_st1.root");
  //proton-AntiLambda index = 1
  file0data[1]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CFOutput_pAL_ApL_st1.root");
  file0MC[1]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/MCCFOutput_pAL_ApL_st1.root");
  //Antiproton-Lambda index = 3
  file0data[3]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CFOutput_pAL_ApL_st1.root");
  file0MC[3]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/MCCFOutput_pAL_ApL_st1.root");

  file0data[2]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CFOutput_LAL_ALL_st1.root");
  file0MC[2]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/MCCFOutput_LAL_ALL_st1.root");
  }
  else if(strcmp(STflag, add2)==0)
  {
  std::cout << "Sphericity [0.3,0.7] \n";
  file0data[0]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CFOutput_pAp_App_st2.root");
  file0MC[0]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/MCCFOutput_pAp_App_st2.root");
  //proton-AntiLambda index = 1
  file0data[1]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CFOutput_pAL_ApL_st2.root");
  file0MC[1]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/MCCFOutput_pAL_ApL_st2.root");
  //Antiproton-Lambda index = 3
  file0data[3]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CFOutput_pAL_ApL_st2.root");
  file0MC[3]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/MCCFOutput_pAL_ApL_st2.root");

  file0data[2]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CFOutput_LAL_ALL_st2.root");
  file0MC[2]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/MCCFOutput_LAL_ALL_st2.root");
  }
  else if(strcmp(STflag, add3)==0)
  {
  std::cout << "Sphericity [0.7,1.] \n";
  file0data[0]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CFOutput_pAp_App_st3.root");
  file0MC[0]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/MCCFOutput_pAp_App_st3.root");
  //proton-AntiLambda index = 1
  file0data[1]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CFOutput_pAL_ApL_st3.root");
  file0MC[1]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/MCCFOutput_pAL_ApL_st3.root");
  //Antiproton-Lambda index = 3
  file0data[3]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CFOutput_pAL_ApL_st3.root");
  file0MC[3]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/MCCFOutput_pAL_ApL_st3.root");

  file0data[2]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CFOutput_LAL_ALL_st3.root");
  file0MC[2]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/MCCFOutput_LAL_ALL_st3.root");
  }
  else if(strcmp(STflag, add4)==0)
  {
  file0data[0]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CFOutput_pAp_App_full.root");
  file0MC[0]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/MCCFOutput_pAp_App_full.root");
  //proton-AntiLambda index = 1
  file0data[1]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CFOutput_pAL_ApL_full.root");
  file0MC[1]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/MCCFOutput_pAL_ApL_full.root");
  //Antiproton-Lambda index = 3
  file0data[3]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CFOutput_pAL_ApL_full.root");
  file0MC[3]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/MCCFOutput_pAL_ApL_full.root");

  file0data[2]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CFOutput_LAL_ALL_full.root");
  file0MC[2]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/MCCFOutput_LAL_ALL_full.root");
  }
}
else if(whichfolder){
  if(strcmp(STflag, add1)==0)
  {
  std::cout << "Sphericity [0.,0.3] \n";
  file0data[0]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CorrectBin/CFOutput_pAp_App_st1.root");
  file0MC[0]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/CorrectBin/Pythia/MCCFOutput_pAp_App_st1.root");
  //proton-AntiLambda index = 1
  file0data[1]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CorrectBin/CFOutput_pAL_ApL_st1.root");
  file0MC[1]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/CorrectBin/Pythia/MCCFOutput_pAL_ApL_st1.root");
  //Antiproton-Lambda index = 3
  file0data[3]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CorrectBin/CFOutput_pAL_ApL_st1.root");
  file0MC[3]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/CorrectBin/Pythia/MCCFOutput_pAL_ApL_st1.root");

  file0data[2]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CorrectBin/CFOutput_LAL_ALL_st1.root");
  file0MC[2]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/CorrectBin/Pythia/MCCFOutput_LAL_ALL_st1.root");
  }
  else if(strcmp(STflag, add2)==0)
  {
  std::cout << "Sphericity [0.3,0.7] \n";
  file0data[0]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CorrectBin/CFOutput_pAp_App_st2.root");
  file0MC[0]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/CorrectBin/Pythia/MCCFOutput_pAp_App_st2.root");
  //proton-AntiLambda index = 1
  file0data[1]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CorrectBin/CFOutput_pAL_ApL_st2.root");
  file0MC[1]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/CorrectBin/Pythia/MCCFOutput_pAL_ApL_st2.root");
  //Antiproton-Lambda index = 3
  file0data[3]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CorrectBin/CFOutput_pAL_ApL_st2.root");
  file0MC[3]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/CorrectBin/Pythia/MCCFOutput_pAL_ApL_st2.root");

  file0data[2]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CorrectBin/CFOutput_LAL_ALL_st2.root");
  file0MC[2]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/CorrectBin/Pythia/MCCFOutput_LAL_ALL_st2.root");
  }
  else if(strcmp(STflag, add3)==0)
  {
  std::cout << "Sphericity [0.7,1.] \n";
  file0data[0]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CorrectBin/CFOutput_pAp_App_st3.root");
  file0MC[0]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/CorrectBin/Pythia/MCCFOutput_pAp_App_st3.root");
  //proton-AntiLambda index = 1
  file0data[1]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CorrectBin/CFOutput_pAL_ApL_st3.root");
  file0MC[1]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/CorrectBin/Pythia/MCCFOutput_pAL_ApL_st3.root");
  //Antiproton-Lambda index = 3
  file0data[3]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CorrectBin/CFOutput_pAL_ApL_st3.root");
  file0MC[3]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/CorrectBin/Pythia/MCCFOutput_pAL_ApL_st3.root");

  file0data[2]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CorrectBin/CFOutput_LAL_ALL_st3.root");
  file0MC[2]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/CorrectBin/Pythia/MCCFOutput_LAL_ALL_st3.root");
  }
  else if(strcmp(STflag, add4)==0)
  {
  file0data[0]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CorrectBin/CFOutput_pAp_App_full.root");
  file0MC[0]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/CorrectBin/Pythia/MCCFOutput_pAp_App_full.root");
  //proton-AntiLambda index = 1
  file0data[1]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CorrectBin/CFOutput_pAL_ApL_full.root");
  file0MC[1]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/CorrectBin/Pythia/MCCFOutput_pAL_ApL_full.root");
  //Antiproton-Lambda index = 3
  file0data[3]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CorrectBin/CFOutput_pAL_ApL_full.root");
  file0MC[3]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/CorrectBin/Pythia/MCCFOutput_pAL_ApL_full.root");

  file0data[2]=TFile::Open("../../../BBbar/GentleFemto_Output/data/CorrectBin/CFOutput_LAL_ALL_full.root");
  file0MC[2]=TFile::Open("../../../BBbar/GentleFemto_Output/mc/CorrectBin/Pythia/MCCFOutput_LAL_ALL_full.root");
  }
}
//REMEMBER THAT WHEN ACCESSING DIRECTLY FROM FILE USE FindObjectAny!!!!
   TH1F* hCFdata[3];
   TH1F* hCFmc[3];
   TList* listPairDistdata[4];
   TList* listPairDistMC[4];
   TList* listMEDistdata[4];
   TList* listMEDistMC[4];
   TH1F* hMEdata[4];
   TH1F* hMEmc[4];

   listPairDistdata[0]=(TList*)(file0data[0]->FindObjectAny("PairDist"));
   listPairDistMC[0]=(TList*)(file0MC[0]->FindObjectAny("PairDist"));
   //proton - AntiLambda index = 1
   listPairDistdata[1]=(TList*)(file0data[1]->FindObjectAny("PairDist"));
   listPairDistMC[1]=(TList*)(file0MC[1]->FindObjectAny("PairDist"));
   //Antiproton - Lambda index = 3
   listPairDistdata[3]=(TList*)(file0data[3]->FindObjectAny("AntiPairDist"));
   listPairDistMC[3]=(TList*)(file0MC[3]->FindObjectAny("AntiPairDist"));

   listPairDistdata[2]=(TList*)(file0data[2]->FindObjectAny("PairDist"));
   listPairDistMC[2]=(TList*)(file0MC[2]->FindObjectAny("PairDist"));


   listMEDistdata[0]=(TList*)listPairDistdata[0]->FindObject("PairReweighted");
   listMEDistMC[0]=(TList*)listPairDistMC[0]->FindObject("PairReweighted");
   //proton - AntiLambda index = 1
   listMEDistdata[1]=(TList*)listPairDistdata[1]->FindObject("PairReweighted");
   listMEDistMC[1]=(TList*)listPairDistMC[1]->FindObject("PairReweighted");
   //Antiproton - Lambda index = 3
   listMEDistdata[3]=(TList*)listPairDistdata[3]->FindObject("PairReweighted");
   listMEDistMC[3]=(TList*)listPairDistMC[3]->FindObject("PairReweighted");

   listMEDistdata[2]=(TList*)listPairDistdata[2]->FindObject("PairReweighted");
   listMEDistMC[2]=(TList*)listPairDistMC[2]->FindObject("PairReweighted");


    hMEdata[0]=(TH1F*)listMEDistdata[0]->FindObject("MEDist_Particle0_Particle1_Shifted_Reweighted");
    hMEdata[0]->Sumw2();
    hMEmc[0]=(TH1F*)listMEDistMC[0]->FindObject("MEDist_Particle0_Particle1_Shifted_Reweighted");
    hMEmc[0]->Sumw2();
    //proton-AntiLambda index = 1
    hMEdata[1]=(TH1F*)listMEDistdata[1]->FindObject("MEDist_Particle0_Particle3_Shifted_FixShifted_Rebinned_5_Reweighted");
    hMEdata[1]->Sumw2();
    hMEmc[1]=(TH1F*)listMEDistMC[1]->FindObject("MEDist_Particle0_Particle3_Shifted_FixShifted_Rebinned_5_Reweighted");
    hMEmc[1]->Sumw2();
    //Antiproton-Lambda index = 3
    hMEdata[3]=(TH1F*)listMEDistdata[3]->FindObject("MEDist_Particle1_Particle2_Shifted_FixShifted_Rebinned_5_Reweighted");
    hMEdata[3]->Sumw2();
    hMEmc[3]=(TH1F*)listMEDistMC[3]->FindObject("MEDist_Particle1_Particle2_Shifted_FixShifted_Rebinned_5_Reweighted");
    hMEmc[3]->Sumw2();
    //Lambda-AntiLambda index = 2
    hMEdata[2]=(TH1F*)listMEDistdata[2]->FindObject("MEDist_Particle2_Particle3_Shifted_FixShifted_Rebinned_5_Reweighted");
    hMEdata[2]->Sumw2();
    hMEmc[2]=(TH1F*)listMEDistMC[2]->FindObject("MEDist_Particle2_Particle3_Shifted_FixShifted_Rebinned_5_Reweighted");
    hMEmc[2]->Sumw2();


// Remember to call the Sumw2() for the error propagation, once done here is good for the whole
   hCFdata[0] = (TH1F*)(file0data[0]->FindObjectAny("hCk_ReweightedpApData_0"));
   hCFdata[0]->Sumw2();

   hCFmc[0] = (TH1F*)(file0MC[0]->FindObjectAny("hCk_ReweightedpApMC_0"));
   hCFmc[0]->Sumw2();

   hCFdata[1] = (TH1F*)(file0data[1]->FindObjectAny("hCk_ReweightedpALData_0"));
   hCFdata[1]->Sumw2();
   hCFmc[1] = (TH1F*)(file0MC[1]->FindObjectAny("hCk_ReweightedpALMC_0"));
   hCFmc[1]->Sumw2();
   hCFdata[2] = (TH1F*)(file0data[2]->FindObjectAny("hCk_ReweightedLALData_0"));
   hCFdata[2]->Sumw2();
   hCFmc[2] = (TH1F*)(file0MC[2]->FindObjectAny("hCk_ReweightedLALMC_0"));
   hCFmc[2]->Sumw2();

    TH1F *ratioCFDataMC[3];
    TH1F *ratioMEDataMC[4];
    TH1F *diffCF[3];
// Evaluating Ratio CF Data/MC
    std::cout<<"CF data/MC for Proton-Antiproton"<<std::endl;
    std::cout<<"number of bins data: "<<hCFdata[0]->GetNbinsX()<<std::endl;
    std::cout<<"number of bins MC: "<<hCFmc[0]->GetNbinsX()<<std::endl;
    ratioCFDataMC[0] = (TH1F*)hCFdata[0]->Clone(Form("%s_clone_ratioCFDataMC0",hCFdata[0]->GetName()));
    ratioCFDataMC[0]->Divide(hCFmc[0]);
    std::cout<<"CF data/MC for p-AL+Ap-L"<<std::endl;
    std::cout<<"number of bins data: "<<hCFdata[1]->GetNbinsX()<<std::endl;
    std::cout<<"number of bins MC: "<<hCFmc[1]->GetNbinsX()<<std::endl;
    ratioCFDataMC[1] = (TH1F*)hCFdata[1]->Clone(Form("%s_clone_ratioCFDataMC1",hCFdata[1]->GetName()));
    ratioCFDataMC[1]->Divide(hCFmc[1]);
    std::cout<<"CF data/MC for Lambda-AntiLambda"<<std::endl;
    std::cout<<"number of bins data: "<<hCFdata[2]->GetNbinsX()<<std::endl;
    std::cout<<"number of bins MC: "<<hCFmc[2]->GetNbinsX()<<std::endl;
    ratioCFDataMC[2] = (TH1F*)hCFdata[2]->Clone(Form("%s_clone_ratioCFDataMC2",hCFdata[2]->GetName()));
    ratioCFDataMC[2]->Divide(hCFmc[2]);

    std::cout<<"------------------------------"<<std::endl;

    // Evaluating Ratio ME Distributions Data/MC

    std::cout<<"ME Dist data/MC for Proton-Antiproton"<<std::endl;
     ratioMEDataMC[0] = (TH1F*)hMEdata[0]->Clone(Form("%s_clone_ratioMEDataMC0",hMEdata[0]->GetName()));
     ratioMEDataMC[0]->Divide(hMEmc[0]);
     //proton-AntiLambda index = 1
     std::cout<<"ME Dist data/MC for Proton-AntiLambda"<<std::endl;
     std::cout<<"number of bins data: "<<hMEdata[1]->GetNbinsX()<<std::endl;
     std::cout<<"number of bins MC: "<<hMEmc[1]->GetNbinsX()<<std::endl;

     ratioMEDataMC[1] = (TH1F*)hMEdata[1]->Clone(Form("%s_clone_ratioMEDataMC1",hMEdata[1]->GetName()));
     ratioMEDataMC[1]->Divide(hMEmc[1]);

     //Antiproton-Lambda index = 3
     std::cout<<"ME Dist data/MC for AntiProton-Lambda"<<std::endl;
     ratioMEDataMC[3] = (TH1F*)hMEdata[3]->Clone(Form("%s_clone_ratioMEDataMC3",hMEdata[3]->GetName()));
     ratioMEDataMC[3]->Divide(hMEmc[3]);

     std::cout<<"ME Dist data/MC for Lambda-AntiLambda"<<std::endl;
     std::cout<<"number of bins data: "<<hMEdata[2]->GetNbinsX()<<std::endl;
     std::cout<<"number of bins MC: "<<hMEmc[2]->GetNbinsX()<<std::endl;
     ratioMEDataMC[2] = (TH1F*)hMEdata[2]->Clone(Form("%s_clone_ratioMEDataMC2",hMEdata[2]->GetName()));
     ratioMEDataMC[2]->Divide(hMEmc[2]);

     std::cout<<"------------------------------"<<std::endl;

     // Evaluating Fit of MC with Gaus, Exp and Pol2
     TF1 *fBackgroundGaus[3];
     TF1 *fBackgroundExp[3];
     TF1 *fBackgroundPol[3];
     TF1 *fBackgroundPol3[3];

     TH1F *fitGaus[3];
     TH1F *fitExp[3];
     TH1F *fitPol[3];
     TH1F *fitPol3[3];

     Double_t chisqGaus[3];
     Double_t chisqExp[3];
     Double_t chisqPol[3];
     Double_t chisqPol3[3];

     Double_t pargaus[3];
     Double_t parexp[3];
     Double_t parpol[3];
     Double_t parpol3[3];

     TH1F *diffCFGaus[3];
     TH1F *diffCFExp[3];
     TH1F *diffCFPol[3];
     TH1F *diffCFPol3[3];

     TH1F *ratioUnity[3];// samefit/samefit
     TH1F *ratio1[3];
     TH1F *ratio2[3];
     TH1F *ratioHisto[3];  //divided by histo


     fBackgroundGaus[0] = new TF1("fBackgroundGaus0","gaus(0)",0.,maxfit);
     fBackgroundExp[0] = new TF1("fBackgroundExp0","expo(0)",0.,maxfit);
     fBackgroundPol[0] = new TF1("fBackgroundPol0","pol2(0)",0.,maxfit);
     fBackgroundPol3[0] = new TF1("fBackgroundPol30","pol3(0)",0.,maxfit);

     fBackgroundGaus[1] = new TF1("fBackgroundGaus1","gaus(0)",0.,maxfit);
     fBackgroundExp[1] = new TF1("fBackgroundExp1","expo(0)",0.,maxfit);
     fBackgroundPol[1] = new TF1("fBackgroundPol1","pol2(0)",0.,maxfit);
     fBackgroundPol3[1] = new TF1("fBackgroundPol31","pol3(0)",0.,maxfit);

     fBackgroundGaus[2] = new TF1("fBackgroundGaus2","gaus(0)",0.,maxfit);
     fBackgroundExp[2] = new TF1("fBackgroundExp2","expo(0)",0.,maxfit);
     fBackgroundPol[2] = new TF1("fBackgroundPol2","pol2(0)",0.,maxfit);
     fBackgroundPol3[2] = new TF1("fBackgroundPol32","pol3(0)",0.,maxfit);

     fitGaus[0] = (TH1F*)hCFmc[0]->Clone(Form("%s_clone_fitGaus0",hCFmc[0]->GetName()));
     fitExp[0] = (TH1F*)hCFmc[0]->Clone(Form("%s_clone_fitExp0",hCFmc[0]->GetName()));
     fitPol[0] = (TH1F*)hCFmc[0]->Clone(Form("%s_clone_fitPol0",hCFmc[0]->GetName()));
     fitPol3[0] = (TH1F*)hCFmc[0]->Clone(Form("%s_clone_fitPol30",hCFmc[0]->GetName()));

     fitGaus[1] = (TH1F*)hCFmc[1]->Clone(Form("%s_clone_fitGaus1",hCFmc[1]->GetName()));
     fitExp[1] = (TH1F*)hCFmc[1]->Clone(Form("%s_clone_fitExp1",hCFmc[1]->GetName()));
     fitPol[1] = (TH1F*)hCFmc[1]->Clone(Form("%s_clone_fitPol1",hCFmc[1]->GetName()));
     fitPol3[1] = (TH1F*)hCFmc[1]->Clone(Form("%s_clone_fitPol31",hCFmc[1]->GetName()));

     fitGaus[2] = (TH1F*)hCFmc[2]->Clone(Form("%s_clone_fitGaus2",hCFmc[2]->GetName()));
     fitExp[2] = (TH1F*)hCFmc[2]->Clone(Form("%s_clone_fitExp2",hCFmc[2]->GetName()));
     fitPol[2] = (TH1F*)hCFmc[2]->Clone(Form("%s_clone_fitPol2",hCFmc[2]->GetName()));
     fitPol3[2] = (TH1F*)hCFmc[2]->Clone(Form("%s_clone_fitPol32",hCFmc[2]->GetName()));

   for(unsigned i=0; i<3; i++){

     fBackgroundGaus[i]->SetParameter(0,1.);
     fBackgroundGaus[i]->SetParameter(1,1.);
     fBackgroundGaus[i]->SetParameter(2,1.);

     fBackgroundExp[i]->SetParameter(0,1.);
     fBackgroundExp[i]->SetParameter(1,1.);

     fBackgroundPol[i]->SetParameter(0,1.);
     fBackgroundPol[i]->SetParameter(1,1.);
     fBackgroundPol[i]->SetParameter(2,1.);

     fBackgroundPol3[i]->SetParameter(0,1.);
     fBackgroundPol3[i]->SetParameter(1,1.);
     fBackgroundPol3[i]->SetParameter(2,1.);
     fBackgroundPol3[i]->SetParameter(3,1.);

// std::cout<<"---------------------------------------"<<std::endl;
// std::cout<<"---------------GAUSS FIT-----------------"<<std::endl;
     fitGaus[i]->Fit(fBackgroundGaus[i],"SQR0","",0.,maxfit);
     fitExp[i]->Fit(fBackgroundExp[i],"SQR0","",0.,maxfit);
     fitPol[i]->Fit(fBackgroundPol[i],"SQR0","",0.,maxfit);
     fitPol3[i]->Fit(fBackgroundPol3[i],"SQR0","",0.,maxfit);

     chisqGaus[i] = fBackgroundGaus[i]->GetChisquare();
     pargaus[i] = fBackgroundGaus[i]->GetNDF();
     chisqExp[i] = fBackgroundExp[i]->GetChisquare();
     parexp[i] = fBackgroundExp[i]->GetNDF();
     chisqPol[i] = fBackgroundPol[i]->GetChisquare();
     parpol[i] = fBackgroundPol[i]->GetNDF();
     chisqPol3[i] = fBackgroundPol3[i]->GetChisquare();
     parpol3[i] = fBackgroundPol3[i]->GetNDF();
}
//     TFitResultPtr Fit(TF1 *function, Option_t *option, Option_t *goption,Axis_t xxmin, Axis_t  xxmax)
// option S = The result of the fit is returned in the TFitResultPtr
//        R = Use the range specified in the function range
//        Q = Quiet mode (minimum printing)
//        0 = Do not plot the result of the fit. By default the fitted function is drawn unless the option “N” above is specified
std::cout<<"-------------------------------"<<std::endl;
std::cout<<"Subtracting fit background to data"<<std::endl;
std::cout<<std::endl;
std::cout<<"CF data - MC for Proton-Antiproton"<<std::endl;
diffCFGaus[0] = (TH1F*)hCFdata[0]->Clone(Form("%s_clone_diffCFGaus0",hCFdata[0]->GetName()));//Append the name gaus or exp
diffCFGaus[0]->Divide(fBackgroundGaus[0]);
diffCFExp[0] = (TH1F*)hCFdata[0]->Clone(Form("%s_clone_diffCFExp0",hCFdata[0]->GetName()));//Append the name gaus or exp
diffCFExp[0]->Divide(fBackgroundExp[0]);
diffCFPol[0] = (TH1F*)hCFdata[0]->Clone(Form("%s_clone_diffCFPol0",hCFdata[0]->GetName()));//Append the name gaus or exp
diffCFPol[0]->Divide(fBackgroundPol[0]);
diffCFPol3[0] = (TH1F*)hCFdata[0]->Clone(Form("%s_clone_diffCFPol30",hCFdata[0]->GetName()));//Append the name gaus or exp
diffCFPol3[0]->Divide(fBackgroundPol3[0]);
std::cout<<"CF data - MC for Proton-AntiLambda + AntiProton-Lambda"<<std::endl;
diffCFGaus[1] = (TH1F*)hCFdata[1]->Clone(Form("%s_clone_diffCFGaus1",hCFdata[1]->GetName()));//Append the name gaus or exp
diffCFGaus[1]->Divide(fBackgroundGaus[1]);
diffCFExp[1] = (TH1F*)hCFdata[1]->Clone(Form("%s_clone_diffCFExp1",hCFdata[1]->GetName()));//Append the name gaus or exp
diffCFExp[1]->Divide(fBackgroundExp[1]);
diffCFPol[1] = (TH1F*)hCFdata[1]->Clone(Form("%s_clone_diffCFPol1",hCFdata[1]->GetName()));//Append the name gaus or exp
diffCFPol[1]->Divide(fBackgroundPol[1]);
diffCFPol3[1] = (TH1F*)hCFdata[1]->Clone(Form("%s_clone_diffCFPol31",hCFdata[1]->GetName()));//Append the name gaus or exp
diffCFPol3[1]->Divide(fBackgroundPol3[1]);
std::cout<<"CF data - MC for Lambda-AntiLambda"<<std::endl;
diffCFGaus[2] = (TH1F*)hCFdata[2]->Clone(Form("%s_clone_diffCFGaus2",hCFdata[2]->GetName()));//Append the name gaus or exp
diffCFGaus[2]->Divide(fBackgroundGaus[2]);
diffCFExp[2] = (TH1F*)hCFdata[2]->Clone(Form("%s_clone_diffCFExp2",hCFdata[2]->GetName()));//Append the name gaus or exp
diffCFExp[2]->Divide(fBackgroundExp[2]);
diffCFPol[2] = (TH1F*)hCFdata[2]->Clone(Form("%s_clone_diffCFPol2",hCFdata[2]->GetName()));//Append the name gaus or exp
diffCFPol[2]->Divide(fBackgroundPol[2]);
diffCFPol3[2] = (TH1F*)hCFdata[2]->Clone(Form("%s_clone_diffCFPol32",hCFdata[2]->GetName()));//Append the name gaus or exp
diffCFPol3[2]->Divide(fBackgroundPol3[2]);

std::cout<<"----------------------------------------------"<<std::endl;
std::cout<<"Proton-AntiProton"<<std::endl;
std::cout<<"ChiSquareGAUS/NDF= "<<chisqGaus[0]/pargaus[0]<<std::endl;
std::cout<<"ChiSquareEXP/NDF= "<<chisqExp[0]/parexp[0]<<std::endl;
std::cout<<"ChiSquarePOL2/NDF= "<<chisqPol[0]/parpol[0]<<std::endl;
std::cout<<"ChiSquarePOL3/NDF= "<<chisqPol3[0]/parpol3[0]<<std::endl;

std::cout<<"Proton-AntiLambda + Antiproton-Lambda"<<std::endl;
std::cout<<"ChiSquareGAUS/NDF= "<<chisqGaus[1]/pargaus[1]<<std::endl;
std::cout<<"ChiSquareEXP/NDF= "<<chisqExp[1]/parexp[1]<<std::endl;
std::cout<<"ChiSquarePOL2/NDF= "<<chisqPol[1]/parpol[1]<<std::endl;
std::cout<<"ChiSquarePOL3/NDF= "<<chisqPol3[1]/parpol3[1]<<std::endl;

std::cout<<"Lambda-AntiLambda"<<std::endl;
std::cout<<"ChiSquareGAUS/NDF= "<<chisqGaus[2]/pargaus[2]<<std::endl;
std::cout<<"ChiSquareEXP/NDF= "<<chisqExp[2]/parexp[2]<<std::endl;
std::cout<<"ChiSquarePOL2/NDF= "<<chisqPol[2]/parpol[2]<<std::endl;
std::cout<<"ChiSquarePOL3/NDF= "<<chisqPol3[2]/parpol3[2]<<std::endl;

std::cout<<"----------------------------------------------"<<std::endl;

std::cout<<"-------------------------------"<<std::endl;
std::cout<<"Making ratios of CF fitted"<<std::endl;
std::cout<<std::endl;
std::cout<<"Proton-Antiproton"<<std::endl;

ratioUnity[0] = (TH1F*)diffCFPol3[0]->Clone(Form("%s_clone_ratioUnity0",diffCFPol3[0]->GetName()));//Append the name gaus or exp
ratioUnity[0]->Divide(diffCFPol3[0]);
ratio1[0] = (TH1F*)diffCFPol3[0]->Clone(Form("%s_clone_ratio10",diffCFPol3[0]->GetName()));//Append the name gaus or exp
ratio1[0]->Divide(diffCFGaus[0]);
ratio2[0] = (TH1F*)diffCFPol3[0]->Clone(Form("%s_clone_ratio20",diffCFPol3[0]->GetName()));//Append the name gaus or exp
ratio2[0]->Divide(diffCFPol[0]);
ratioHisto[0] = (TH1F*)diffCFPol3[0]->Clone(Form("%s_clone_ratioHisto0",diffCFPol3[0]->GetName()));//Append the name gaus or exp
ratioHisto[0]->Divide(hCFmc[0]);

std::cout<<"Proton-AntiLambda + AntiProton-Lambda"<<std::endl;

ratioUnity[1] = (TH1F*)diffCFGaus[1]->Clone(Form("%s_clone_ratioUnity1",diffCFGaus[1]->GetName()));//Append the name gaus or exp
ratioUnity[1]->Divide(diffCFGaus[1]);
ratio1[1] = (TH1F*)diffCFGaus[1]->Clone(Form("%s_clone_ratio11",diffCFGaus[1]->GetName()));//Append the name gaus or exp
ratio1[1]->Divide(diffCFPol3[1]);
ratio2[1] = (TH1F*)diffCFGaus[1]->Clone(Form("%s_clone_ratio21",diffCFGaus[1]->GetName()));//Append the name gaus or exp
ratio2[1]->Divide(diffCFPol[1]);
ratioHisto[1] = (TH1F*)diffCFGaus[1]->Clone(Form("%s_clone_ratioHisto1",diffCFGaus[1]->GetName()));//Append the name gaus or exp
ratioHisto[1]->Divide(hCFmc[1]);

std::cout<<"Lambda-AntiLambda"<<std::endl;

ratioUnity[2] = (TH1F*)diffCFPol3[2]->Clone(Form("%s_clone_ratioUnity2",diffCFPol3[2]->GetName()));//Append the name gaus or exp
ratioUnity[2]->Divide(diffCFPol3[2]);
ratio1[2] = (TH1F*)diffCFPol3[2]->Clone(Form("%s_clone_ratio12",diffCFPol3[2]->GetName()));//Append the name gaus or exp
ratio1[2]->Divide(diffCFGaus[2]);
ratio2[2] = (TH1F*)diffCFPol3[2]->Clone(Form("%s_clone_ratio22",diffCFPol3[2]->GetName()));//Append the name gaus or exp
ratio2[2]->Divide(diffCFPol[2]);
ratioHisto[2] = (TH1F*)diffCFPol3[2]->Clone(Form("%s_clone_ratioHisto2",diffCFPol3[2]->GetName()));//Append the name gaus or exp
ratioHisto[2]->Divide(hCFmc[2]);

std::cout<<"----------------------------------------------"<<std::endl;
TString foldername;
TString foldernameplot;

if(maxfit==0.4){
foldername = "/Users/Valentina/cernbox/Analysis/BBbar/GentleFemto_Output/MiniJetCheck/FitRangeDefault04/";
foldernameplot = "/Users/Valentina/cernbox/Analysis/BBbar/GentleFemto_Output/MiniJetCheck/FitRangeDefault04/Plot/";
}
else if(maxfit==0.3){
   foldername = "/Users/Valentina/cernbox/Analysis/BBbar/GentleFemto_Output/MiniJetCheck/FitRange03/";
   foldernameplot = "/Users/Valentina/cernbox/Analysis/BBbar/GentleFemto_Output/MiniJetCheck/FitRange03/Plot/";
}
else if(maxfit == 0.35){
   foldername = "/Users/Valentina/cernbox/Analysis/BBbar/GentleFemto_Output/MiniJetCheck/FitRange035/";
   foldernameplot = "/Users/Valentina/cernbox/Analysis/BBbar/GentleFemto_Output/MiniJetCheck/FitRange035/Plot/";
}
else if(maxfit == 0.5){
   foldername = "/Users/Valentina/cernbox/Analysis/BBbar/GentleFemto_Output/MiniJetCheck/FitRange05/";
   foldernameplot = "/Users/Valentina/cernbox/Analysis/BBbar/GentleFemto_Output/MiniJetCheck/FitRange05/Plot/";
}
if(writefile){
FILE *fp = fopen(foldername+"Fitparameters_MC.dat","a");
if (fp!=NULL) {
std::cout<<"Writing to file"<<std::endl;
std::cout<<"Proton-AntiProton"<<std::endl;
if(strcmp(STflag, add1)==0){
  fprintf(fp,"Proton-AntiProton \n");
  fprintf(fp,"Sphericity [0-0.3]\n");
  fprintf(fp,"1) Gauss fit \n");
   for (unsigned igaus=0;igaus<fBackgroundGaus[0]->GetNpar();igaus++)
   {
      Float_t value = fBackgroundGaus[0]->GetParameter(igaus);
      Float_t valuerr = fBackgroundGaus[0]->GetParError(igaus);

      fprintf(fp,"p(%u) =  %f",igaus,value);
      fprintf(fp," +/-  %f",valuerr);
      fprintf(fp,"\n");

   }
   fprintf(fp,"\n");

   fprintf(fp,"2) Pol2 fit \n");
   for (unsigned ipol2=0;ipol2<fBackgroundPol[0]->GetNpar();ipol2++)
   {
      Float_t value = fBackgroundPol[0]->GetParameter(ipol2);
      Float_t valuerr = fBackgroundPol[0]->GetParError(ipol2);

      fprintf(fp,"p(%u) =  %f",ipol2,value);
      fprintf(fp," +/-  %f",valuerr);
      fprintf(fp,"\n");

   }
   fprintf(fp,"\n");

   fprintf(fp,"3) Pol3 fit \n");
   for (unsigned ipol3=0;ipol3<fBackgroundPol[0]->GetNpar();ipol3++)
   {
      Float_t value = fBackgroundPol3[0]->GetParameter(ipol3);
      Float_t valuerr = fBackgroundPol3[0]->GetParError(ipol3);

      fprintf(fp,"p(%u) =  %f",ipol3,value);
      fprintf(fp," +/-  %f",valuerr);

      fprintf(fp,"\n");
   }
   fprintf(fp, "\n");

   fprintf(fp,"GAUSS: chi2/ndf=  %f \n",chisqGaus[0]/pargaus[0]);
   fprintf(fp,"POL2: chi2/ndf=  %f \n",chisqPol[0]/parpol[0]);
   fprintf(fp,"POL3: chi2/ndf=  %f \n",chisqPol3[0]/parpol3[0]);
   fprintf(fp, "\n \n \n");
   fprintf(fp,"Proton-AntiLambda + Antiproton-Lambda \n");
   fprintf(fp,"Sphericity [0-0.3]\n");
   fprintf(fp,"1) Gauss fit \n");
    for (unsigned igaus=0;igaus<fBackgroundGaus[1]->GetNpar();igaus++)
    {
       Float_t value = fBackgroundGaus[1]->GetParameter(igaus);
       Float_t valuerr = fBackgroundGaus[1]->GetParError(igaus);

       fprintf(fp,"p(%u) =  %f",igaus,value);
       fprintf(fp," +/-  %f",valuerr);

       fprintf(fp,"\n");

    }
    fprintf(fp,"\n");

    fprintf(fp,"2) Pol2 fit \n");
    for (unsigned ipol2=0;ipol2<fBackgroundPol[1]->GetNpar();ipol2++)
    {
       Float_t value = fBackgroundPol[1]->GetParameter(ipol2);
       Float_t valuerr = fBackgroundPol[1]->GetParError(ipol2);

       fprintf(fp,"p(%u) =  %f",ipol2,value);
       fprintf(fp," +/-  %f",valuerr);

       fprintf(fp,"\n");

    }
    fprintf(fp,"\n");

    fprintf(fp,"3) Pol3 fit \n");
    for (unsigned ipol3=0;ipol3<fBackgroundPol[1]->GetNpar();ipol3++)
    {
       Float_t value = fBackgroundPol3[1]->GetParameter(ipol3);
       Float_t valuerr = fBackgroundPol3[1]->GetParError(ipol3);

       fprintf(fp,"p(%u) =  %f",ipol3,value);
       fprintf(fp," +/-  %f",valuerr);

       fprintf(fp,"\n");

    }
    fprintf(fp, "\n");

    fprintf(fp,"GAUSS: chi2/ndf=  %f \n",chisqGaus[1]/pargaus[1]);
    fprintf(fp,"POL2: chi2/ndf=  %f \n",chisqPol[1]/parpol[1]);
    fprintf(fp,"POL3: chi2/ndf=  %f \n",chisqPol3[1]/parpol3[1]);

    fprintf(fp, "\n \n \n");
    fprintf(fp,"Lambda-AntiLambda\n");
    fprintf(fp,"Sphericity [0-0.3]\n");
    fprintf(fp,"1) Gauss fit \n");
     for (unsigned igaus=0;igaus<fBackgroundGaus[2]->GetNpar();igaus++)
     {
        Float_t value = fBackgroundGaus[2]->GetParameter(igaus);
        Float_t valuerr = fBackgroundGaus[2]->GetParError(igaus);

        fprintf(fp,"p(%u) =  %f",igaus,value);
        fprintf(fp," +/-  %f",valuerr);

        fprintf(fp,"\n");

     }
     fprintf(fp,"\n");

     fprintf(fp,"2) Pol2 fit \n");
     for (unsigned ipol2=0;ipol2<fBackgroundPol[2]->GetNpar();ipol2++)
     {
        Float_t value = fBackgroundPol[2]->GetParameter(ipol2);
        Float_t valuerr = fBackgroundPol[2]->GetParError(ipol2);

        fprintf(fp,"p(%u) =  %f",ipol2,value);
        fprintf(fp," +/-  %f",valuerr);

        fprintf(fp,"\n");

     }
     fprintf(fp,"\n");

     fprintf(fp,"3) Pol3 fit \n");
     for (unsigned ipol3=0;ipol3<fBackgroundPol[2]->GetNpar();ipol3++)
     {
        Float_t value = fBackgroundPol3[2]->GetParameter(ipol3);
        Float_t valuerr = fBackgroundPol3[2]->GetParError(ipol3);

        fprintf(fp,"p(%u) =  %f",ipol3,value);
        fprintf(fp," +/-  %f",valuerr);

        fprintf(fp,"\n");

     }
     fprintf(fp, "\n");

     fprintf(fp,"GAUSS: chi2/ndf=  %f \n",chisqGaus[2]/pargaus[2]);
     fprintf(fp,"POL2: chi2/ndf=  %f \n",chisqPol[2]/parpol[2]);
     fprintf(fp,"POL3: chi2/ndf=  %f \n",chisqPol3[2]/parpol3[2]);
}
fprintf(fp, "\n \n \n");
fprintf(fp, "--------------------------------------------------\n");

if(strcmp(STflag, add2)==0){
  fprintf(fp,"Proton-AntiProton \n");
  fprintf(fp,"Sphericity [0.3-0.7]\n");
  fprintf(fp,"1) Gauss fit \n");
   for (unsigned igaus=0;igaus<fBackgroundGaus[0]->GetNpar();igaus++)
   {
      Float_t value = fBackgroundGaus[0]->GetParameter(igaus);
      Float_t valuerr = fBackgroundGaus[0]->GetParError(igaus);

      fprintf(fp,"p(%u) =  %f",igaus,value);
      fprintf(fp," +/-  %f",valuerr);

      fprintf(fp,"\n");

   }
   fprintf(fp,"\n");

   fprintf(fp,"2) Pol2 fit \n");
   for (unsigned ipol2=0;ipol2<fBackgroundPol[0]->GetNpar();ipol2++)
   {
      Float_t value = fBackgroundPol[0]->GetParameter(ipol2);
      Float_t valuerr = fBackgroundPol[0]->GetParError(ipol2);

      fprintf(fp,"p(%u) =  %f",ipol2,value);
      fprintf(fp," +/-  %f",valuerr);

      fprintf(fp,"\n");

   }
   fprintf(fp,"\n");

   fprintf(fp,"3) Pol3 fit \n");
   for (unsigned ipol3=0;ipol3<fBackgroundPol[0]->GetNpar();ipol3++)
   {
      Float_t value = fBackgroundPol3[0]->GetParameter(ipol3);
      Float_t valuerr = fBackgroundPol3[0]->GetParError(ipol3);

      fprintf(fp,"p(%u) =  %f",ipol3,value);
      fprintf(fp," +/-  %f",valuerr);

      fprintf(fp,"\n");

   }

   fprintf(fp, "\n");

   fprintf(fp,"GAUSS: chi2/ndf=  %f \n",chisqGaus[0]/pargaus[0]);
   fprintf(fp,"POL2: chi2/ndf=  %f \n",chisqPol[0]/parpol[0]);
   fprintf(fp,"POL3: chi2/ndf=  %f \n",chisqPol3[0]/parpol3[0]);
   fprintf(fp, "\n \n \n");
   fprintf(fp,"Proton-AntiLambda + Antiproton-Lambda \n");
   fprintf(fp,"Sphericity [0.3-0.7]\n");
   fprintf(fp,"1) Gauss fit \n");
    for (unsigned igaus=0;igaus<fBackgroundGaus[1]->GetNpar();igaus++)
    {
       Float_t value = fBackgroundGaus[1]->GetParameter(igaus);
       Float_t valuerr = fBackgroundGaus[1]->GetParError(igaus);

       fprintf(fp,"p(%u) =  %f",igaus,value);
       fprintf(fp," +/-  %f",valuerr);

       fprintf(fp,"\n");

    }
    fprintf(fp,"\n");

    fprintf(fp,"2) Pol2 fit \n");
    for (unsigned ipol2=0;ipol2<fBackgroundPol[1]->GetNpar();ipol2++)
    {
       Float_t value = fBackgroundPol[1]->GetParameter(ipol2);
       Float_t valuerr = fBackgroundPol[1]->GetParError(ipol2);

       fprintf(fp,"p(%u) =  %f",ipol2,value);
       fprintf(fp," +/-  %f",valuerr);

       fprintf(fp,"\n");

    }
    fprintf(fp,"\n");

    fprintf(fp,"3) Pol3 fit \n");
    for (unsigned ipol3=0;ipol3<fBackgroundPol[1]->GetNpar();ipol3++)
    {
       Float_t value = fBackgroundPol3[1]->GetParameter(ipol3);
       Float_t valuerr = fBackgroundPol3[1]->GetParError(ipol3);

       fprintf(fp,"p(%u) =  %f",ipol3,value);
       fprintf(fp," +/-  %f",valuerr);

       fprintf(fp,"\n");

    }
    fprintf(fp, "\n");

    fprintf(fp,"GAUSS: chi2/ndf=  %f \n",chisqGaus[1]/pargaus[1]);
    fprintf(fp,"POL2: chi2/ndf=  %f \n",chisqPol[1]/parpol[1]);
    fprintf(fp,"POL3: chi2/ndf=  %f \n",chisqPol3[1]/parpol3[1]);

    fprintf(fp, "\n \n \n");
    fprintf(fp,"Lambda-AntiLambda\n");
    fprintf(fp,"Sphericity [0.3-0.7]\n");
    fprintf(fp,"1) Gauss fit \n");
     for (unsigned igaus=0;igaus<fBackgroundGaus[2]->GetNpar();igaus++)
     {
        Float_t value = fBackgroundGaus[2]->GetParameter(igaus);
        Float_t valuerr = fBackgroundGaus[2]->GetParError(igaus);

        fprintf(fp,"p(%u) =  %f",igaus,value);
        fprintf(fp," +/-  %f",valuerr);

        fprintf(fp,"\n");

     }
     fprintf(fp,"\n");

     fprintf(fp,"2) Pol2 fit \n");
     for (unsigned ipol2=0;ipol2<fBackgroundPol[2]->GetNpar();ipol2++)
     {
        Float_t value = fBackgroundPol[2]->GetParameter(ipol2);
        Float_t valuerr = fBackgroundPol[2]->GetParError(ipol2);

        fprintf(fp,"p(%u) =  %f",ipol2,value);
        fprintf(fp," +/-  %f",valuerr);

        fprintf(fp,"\n");

     }
     fprintf(fp,"\n");

     fprintf(fp,"3) Pol3 fit \n");
     for (unsigned ipol3=0;ipol3<fBackgroundPol[2]->GetNpar();ipol3++)
     {
        Float_t value = fBackgroundPol3[2]->GetParameter(ipol3);
        Float_t valuerr = fBackgroundPol3[2]->GetParError(ipol3);

        fprintf(fp,"p(%u) =  %f",ipol3,value);
        fprintf(fp," +/-  %f",valuerr);

        fprintf(fp,"\n");

     }
     fprintf(fp, "\n");

     fprintf(fp,"GAUSS: chi2/ndf=  %f \n",chisqGaus[2]/pargaus[2]);
     fprintf(fp,"POL2: chi2/ndf=  %f \n",chisqPol[2]/parpol[2]);
     fprintf(fp,"POL3: chi2/ndf=  %f \n",chisqPol3[2]/parpol3[2]);
}
fprintf(fp, "\n \n \n");
fprintf(fp, "--------------------------------------------------\n");

if(strcmp(STflag, add3)==0){
  fprintf(fp,"Proton-AntiProton \n");
  fprintf(fp,"Sphericity [0.7-1]\n");
  fprintf(fp,"1) Gauss fit \n");
   for (unsigned igaus=0;igaus<fBackgroundGaus[0]->GetNpar();igaus++)
   {
      Float_t value = fBackgroundGaus[0]->GetParameter(igaus);
      Float_t valuerr = fBackgroundGaus[0]->GetParError(igaus);

      fprintf(fp,"p(%u) =  %f",igaus,value);
      fprintf(fp," +/-  %f",valuerr);

      fprintf(fp,"\n");

   }
   fprintf(fp,"\n");

   fprintf(fp,"2) Pol2 fit \n");
   for (unsigned ipol2=0;ipol2<fBackgroundPol[0]->GetNpar();ipol2++)
   {
      Float_t value = fBackgroundPol[0]->GetParameter(ipol2);
      Float_t valuerr = fBackgroundPol[0]->GetParError(ipol2);

      fprintf(fp,"p(%u) =  %f",ipol2,value);
      fprintf(fp," +/-  %f",valuerr);

      fprintf(fp,"\n");

   }
   fprintf(fp,"\n");

   fprintf(fp,"3) Pol3 fit \n");
   for (unsigned ipol3=0;ipol3<fBackgroundPol[0]->GetNpar();ipol3++)
   {
      Float_t value = fBackgroundPol3[0]->GetParameter(ipol3);
      Float_t valuerr = fBackgroundPol3[0]->GetParError(ipol3);

      fprintf(fp,"p(%u) =  %f",ipol3,value);
      fprintf(fp," +/-  %f",valuerr);

      fprintf(fp,"\n");

   }
   fprintf(fp, "\n");

   fprintf(fp,"GAUSS: chi2/ndf=  %f \n",chisqGaus[0]/pargaus[0]);
   fprintf(fp,"POL2: chi2/ndf=  %f \n",chisqPol[0]/parpol[0]);
   fprintf(fp,"POL3: chi2/ndf=  %f \n",chisqPol3[0]/parpol3[0]);
   fprintf(fp, "\n \n \n");
   fprintf(fp,"Proton-AntiLambda + Antiproton-Lambda \n");
   fprintf(fp,"Sphericity [0.7-1]\n");
   fprintf(fp,"1) Gauss fit \n");
    for (unsigned igaus=0;igaus<fBackgroundGaus[1]->GetNpar();igaus++)
    {
       Float_t value = fBackgroundGaus[1]->GetParameter(igaus);
       Float_t valuerr = fBackgroundGaus[1]->GetParError(igaus);

       fprintf(fp,"p(%u) =  %f",igaus,value);
       fprintf(fp," +/-  %f",valuerr);

       fprintf(fp,"\n");

    }
    fprintf(fp,"\n");

    fprintf(fp,"2) Pol2 fit \n");
    for (unsigned ipol2=0;ipol2<fBackgroundPol[1]->GetNpar();ipol2++)
    {
       Float_t value = fBackgroundPol[1]->GetParameter(ipol2);
       Float_t valuerr = fBackgroundPol[1]->GetParError(ipol2);

       fprintf(fp,"p(%u) =  %f",ipol2,value);
       fprintf(fp," +/-  %f",valuerr);

       fprintf(fp,"\n");

    }
    fprintf(fp,"\n");

    fprintf(fp,"3) Pol3 fit \n");
    for (unsigned ipol3=0;ipol3<fBackgroundPol[1]->GetNpar();ipol3++)
    {
       Float_t value = fBackgroundPol3[1]->GetParameter(ipol3);
       Float_t valuerr = fBackgroundPol3[1]->GetParError(ipol3);

       fprintf(fp,"p(%u) =  %f",ipol3,value);
       fprintf(fp," +/-  %f",valuerr);

       fprintf(fp,"\n");

    }
    fprintf(fp, "\n");

    fprintf(fp,"GAUSS: chi2/ndf=  %f \n",chisqGaus[1]/pargaus[1]);
    fprintf(fp,"POL2: chi2/ndf=  %f \n",chisqPol[1]/parpol[1]);
    fprintf(fp,"POL3: chi2/ndf=  %f \n",chisqPol3[1]/parpol3[1]);

    fprintf(fp, "\n \n \n");
    fprintf(fp,"Lambda-AntiLambda\n");
    fprintf(fp,"Sphericity [0.7-1]\n");
    fprintf(fp,"1) Gauss fit \n");
     for (unsigned igaus=0;igaus<fBackgroundGaus[2]->GetNpar();igaus++)
     {
        Float_t value = fBackgroundGaus[2]->GetParameter(igaus);
        Float_t valuerr = fBackgroundGaus[2]->GetParError(igaus);

        fprintf(fp,"p(%u) =  %f",igaus,value);
        fprintf(fp," +/-  %f",valuerr);

        fprintf(fp,"\n");

     }
     fprintf(fp,"\n");

     fprintf(fp,"2) Pol2 fit \n");
     for (unsigned ipol2=0;ipol2<fBackgroundPol[2]->GetNpar();ipol2++)
     {
        Float_t value = fBackgroundPol[2]->GetParameter(ipol2);
        Float_t valuerr = fBackgroundPol[2]->GetParError(ipol2);

        fprintf(fp,"p(%u) =  %f",ipol2,value);
        fprintf(fp," +/-  %f",valuerr);

        fprintf(fp,"\n");

     }
     fprintf(fp,"\n");

     fprintf(fp,"3) Pol3 fit \n");
     for (unsigned ipol3=0;ipol3<fBackgroundPol[2]->GetNpar();ipol3++)
     {
        Float_t value = fBackgroundPol3[2]->GetParameter(ipol3);
        Float_t valuerr = fBackgroundPol3[2]->GetParError(ipol3);

        fprintf(fp,"p(%u) =  %f",ipol3,value);
        fprintf(fp," +/-  %f",valuerr);

        fprintf(fp,"\n");

     }
     fprintf(fp, "\n");

     fprintf(fp,"GAUSS: chi2/ndf=  %f \n",chisqGaus[2]/pargaus[2]);
     fprintf(fp,"POL2: chi2/ndf=  %f \n",chisqPol[2]/parpol[2]);
     fprintf(fp,"POL3: chi2/ndf=  %f \n",chisqPol3[2]/parpol3[2]);
}
fprintf(fp, "\n \n \n");
fprintf(fp, "--------------------------------------------------\n");

if(strcmp(STflag, add4)==0){
  fprintf(fp,"Proton-AntiProton \n");
  fprintf(fp,"Sphericity [0-1]\n");
  fprintf(fp,"1) Gauss fit \n");
   for (unsigned igaus=0;igaus<fBackgroundGaus[0]->GetNpar();igaus++)
   {
      Float_t value = fBackgroundGaus[0]->GetParameter(igaus);
      Float_t valuerr = fBackgroundGaus[0]->GetParError(igaus);

      fprintf(fp,"p(%u) =  %f",igaus,value);
      fprintf(fp," +/-  %f",valuerr);

      fprintf(fp,"\n");

   }
   fprintf(fp,"\n");

   fprintf(fp,"2) Pol2 fit \n");
   for (unsigned ipol2=0;ipol2<fBackgroundPol[0]->GetNpar();ipol2++)
   {
      Float_t value = fBackgroundPol[0]->GetParameter(ipol2);
      Float_t valuerr = fBackgroundPol[0]->GetParError(ipol2);

      fprintf(fp,"p(%u) =  %f",ipol2,value);
      fprintf(fp," +/-  %f",valuerr);

      fprintf(fp,"\n");

   }
   fprintf(fp,"\n");

   fprintf(fp,"3) Pol3 fit \n");
   for (unsigned ipol3=0;ipol3<fBackgroundPol[0]->GetNpar();ipol3++)
   {
      Float_t value = fBackgroundPol3[0]->GetParameter(ipol3);
      Float_t valuerr = fBackgroundPol3[0]->GetParError(ipol3);

      fprintf(fp,"p(%u) =  %f",ipol3,value);
      fprintf(fp," +/-  %f",valuerr);

      fprintf(fp,"\n");

   }
   fprintf(fp, "\n");

   fprintf(fp,"GAUSS: chi2/ndf=  %f \n",chisqGaus[0]/pargaus[0]);
   fprintf(fp,"POL2: chi2/ndf=  %f \n",chisqPol[0]/parpol[0]);
   fprintf(fp,"POL3: chi2/ndf=  %f \n",chisqPol3[0]/parpol3[0]);
   fprintf(fp, "\n \n \n");
   fprintf(fp,"Proton-AntiLambda + Antiproton-Lambda \n");
   fprintf(fp,"Sphericity [0-1]\n");
   fprintf(fp,"1) Gauss fit \n");
    for (unsigned igaus=0;igaus<fBackgroundGaus[1]->GetNpar();igaus++)
    {
       Float_t value = fBackgroundGaus[1]->GetParameter(igaus);
       Float_t valuerr = fBackgroundGaus[1]->GetParError(igaus);

       fprintf(fp,"p(%u) =  %f",igaus,value);
       fprintf(fp," +/-  %f",valuerr);

       fprintf(fp,"\n");

    }
    fprintf(fp,"\n");

    fprintf(fp,"2) Pol2 fit \n");
    for (unsigned ipol2=0;ipol2<fBackgroundPol[1]->GetNpar();ipol2++)
    {
       Float_t value = fBackgroundPol[1]->GetParameter(ipol2);
       Float_t valuerr = fBackgroundPol[1]->GetParError(ipol2);

       fprintf(fp,"p(%u) =  %f",ipol2,value);
       fprintf(fp," +/-  %f",valuerr);

       fprintf(fp,"\n");

    }
    fprintf(fp,"\n");

    fprintf(fp,"3) Pol3 fit \n");
    for (unsigned ipol3=0;ipol3<fBackgroundPol[1]->GetNpar();ipol3++)
    {
       Float_t value = fBackgroundPol3[1]->GetParameter(ipol3);
       Float_t valuerr = fBackgroundPol3[1]->GetParError(ipol3);

       fprintf(fp,"p(%u) =  %f",ipol3,value);
       fprintf(fp," +/-  %f",valuerr);

       fprintf(fp,"\n");

    }
    fprintf(fp, "\n");

    fprintf(fp,"GAUSS: chi2/ndf=  %f \n",chisqGaus[1]/pargaus[1]);
    fprintf(fp,"POL2: chi2/ndf=  %f \n",chisqPol[1]/parpol[1]);
    fprintf(fp,"POL3: chi2/ndf=  %f \n",chisqPol3[1]/parpol3[1]);

    fprintf(fp, "\n \n \n");
    fprintf(fp,"Lambda-AntiLambda\n");
    fprintf(fp,"Sphericity [0-1]\n");
    fprintf(fp,"1) Gauss fit \n");
     for (unsigned igaus=0;igaus<fBackgroundGaus[2]->GetNpar();igaus++)
     {
        Float_t value = fBackgroundGaus[2]->GetParameter(igaus);
        Float_t valuerr = fBackgroundGaus[2]->GetParError(igaus);

        fprintf(fp,"p(%u) =  %f",igaus,value);
        fprintf(fp," +/-  %f",valuerr);

        fprintf(fp,"\n");

     }
     fprintf(fp,"\n");

     fprintf(fp,"2) Pol2 fit \n");
     for (unsigned ipol2=0;ipol2<fBackgroundPol[2]->GetNpar();ipol2++)
     {
        Float_t value = fBackgroundPol[2]->GetParameter(ipol2);
        Float_t valuerr = fBackgroundPol[2]->GetParError(ipol2);

        fprintf(fp,"p(%u) =  %f",ipol2,value);
        fprintf(fp," +/-  %f",valuerr);

        fprintf(fp,"\n");

     }
     fprintf(fp,"\n");

     fprintf(fp,"3) Pol3 fit \n");
     for (unsigned ipol3=0;ipol3<fBackgroundPol[2]->GetNpar();ipol3++)
     {
        Float_t value = fBackgroundPol3[2]->GetParameter(ipol3);
        Float_t valuerr = fBackgroundPol3[2]->GetParError(ipol3);

        fprintf(fp,"p(%u) =  %f",ipol3,value);
        fprintf(fp," +/-  %f",valuerr);

        fprintf(fp,"\n");

     }
     fprintf(fp, "\n");

     fprintf(fp,"GAUSS: chi2/ndf=  %f \n",chisqGaus[2]/pargaus[2]);
     fprintf(fp,"POL2: chi2/ndf=  %f \n",chisqPol[2]/parpol[2]);
     fprintf(fp,"POL3: chi2/ndf=  %f \n",chisqPol3[2]/parpol3[2]);
}

}

fclose(fp);
}
TCanvas* cGaus[3];
TCanvas* cRat[3];

for(unsigned i=0;i<3;i++){
cGaus[i] = new TCanvas();
cGaus[i]->Divide(1,2);
cGaus[i]->cd(1);
hCFmc[i]->GetXaxis()->SetRangeUser(0.,0.5);
hCFmc[i]->SetTitle("; k* [GeV/c]; C(k*)_{MC}");
hCFmc[i]->Draw();

fBackgroundGaus[i]->SetLineStyle(1);//Solid
fBackgroundGaus[i]->SetLineWidth(2);
fBackgroundGaus[i]->SetLineColor(1);// Black
fBackgroundExp[i]->SetLineStyle(7);//Dashed
fBackgroundExp[i]->SetLineWidth(2);
fBackgroundExp[i]->SetLineColor(2);// Red
fBackgroundPol[i]->SetLineStyle(2);//Dotted
fBackgroundPol[i]->SetLineWidth(2);
fBackgroundPol[i]->SetLineColor(4);// Blue
fBackgroundPol3[i]->SetLineStyle(10);//Dotted
fBackgroundPol3[i]->SetLineWidth(2);
fBackgroundPol3[i]->SetLineColor(6);// Magenta

TLatex letsWriteSomeTextToGranniesCanvas;
letsWriteSomeTextToGranniesCanvas.SetNDC(kTRUE);
letsWriteSomeTextToGranniesCanvas.DrawLatex(0.45, 0.4, TString::Format("#chi^{2}_{NDF-Gauss} = %.3f", chisqGaus[i]/pargaus[i]));
letsWriteSomeTextToGranniesCanvas.DrawLatex(0.45, 0.3, TString::Format("#chi^{2}_{NDF-Exp} = %.3f", chisqExp[i]/parexp[i]));
letsWriteSomeTextToGranniesCanvas.DrawLatex(0.45, 0.2, TString::Format("#chi^{2}_{NDF-Pol} = %.3f", chisqPol[i]/parpol[i]));
letsWriteSomeTextToGranniesCanvas.DrawLatex(0.45, 0.1, TString::Format("#chi^{2}_{NDF-Pol} = %.3f", chisqPol3[i]/parpol3[i]));

auto* leg0= new TLegend(0.65,0.65,0.8,0.8);
leg0->AddEntry(fBackgroundGaus[0], "Gaus", "l");
leg0->AddEntry(fBackgroundExp[0], "Exp", "l");
leg0->AddEntry(fBackgroundPol[0], "Pol2", "l");
leg0->AddEntry(fBackgroundPol3[0], "Pol3", "l");
leg0->Draw("same");
fBackgroundGaus[i]->Draw("same");
fBackgroundExp[i]->Draw("same");
fBackgroundPol[i]->Draw("same");
fBackgroundPol3[i]->Draw("same");

cGaus[i]->cd(2);
diffCFGaus[i]->SetTitle("; k* [GeV/c]; C(k*)");
diffCFGaus[i]->GetXaxis()->SetRangeUser(0.,0.25);
diffCFGaus[i]->SetFillColor(1);
diffCFExp[i]->GetXaxis()->SetRangeUser(0.,0.25);
diffCFExp[i]->SetFillColor(2);
diffCFPol[i]->GetXaxis()->SetRangeUser(0.,0.25);
diffCFPol[i]->SetFillColor(4);
diffCFPol3[i]->GetXaxis()->SetRangeUser(0.,0.25);
diffCFPol3[i]->SetFillColor(6);

diffCFGaus[i]->Draw();
diffCFExp[i]->Draw("same");
diffCFPol[i]->Draw("same");
diffCFPol3[i]->Draw("same");

cRat[i] = new TCanvas();

ratioUnity[i]->GetXaxis()->SetRangeUser(0.,0.25);
ratioUnity[i]->GetYaxis()->SetRangeUser(0.2,1.5);

ratioUnity[i]->SetTitle("; k* [GeV/c]; Ratio");
DreamPlot::SetStyleHisto(ratioUnity[i], 8, 4);//(histo,marker,color)
ratioUnity[i]->Draw();

ratio1[i]->GetXaxis()->SetRangeUser(0.,0.25);
DreamPlot::SetStyleHisto(ratio1[i], 8, 8);
ratio1[i]->Draw("same");

ratio2[i]->GetXaxis()->SetRangeUser(0.,0.25);
DreamPlot::SetStyleHisto(ratio2[i], 8, 7);
ratio2[i]->Draw("same");

ratioHisto[i]->GetXaxis()->SetRangeUser(0.,0.25);

DreamPlot::SetStyleHisto(ratioHisto[i], 8, 1);
ratioHisto[i]->Draw("same");

auto* legr= new TLegend(0.65,0.25,0.8,0.4);

if(i==0){
  legr->AddEntry(ratioUnity[0], "Pol3/Pol3", "l");
  legr->AddEntry(ratio1[0], "Pol3/Gaus", "l");
  legr->AddEntry(ratio2[0], "Pol3/Pol2", "l");
  legr->AddEntry(ratioHisto[0], "Pol3/Histo", "l");
  legr->Draw("same");
}

if(i==1){
  legr->AddEntry(ratioUnity[1], "Gaus/Gaus", "l");
  legr->AddEntry(ratio1[1], "Gaus/Pol3", "l");
  legr->AddEntry(ratio2[1], "Gaus/Pol2", "l");
  legr->AddEntry(ratioHisto[1], "Gaus/Histo", "l");
  legr->Draw("same");
}

if(i==2){
  legr->AddEntry(ratioUnity[2], "Pol3/Pol3", "l");
  legr->AddEntry(ratio1[2], "Pol3/Gaus", "l");
  legr->AddEntry(ratio2[2], "Pol3/Pol2", "l");
  legr->AddEntry(ratioHisto[2], "Pol3/Histo", "l");
  legr->Draw("same");
}

}
// std::cout<<"------------------------------"<<std::endl;

//------------------------------------------------

    //   // Making PDF file

   TFile *file = new TFile(foldername+"MiniJetCheck_test1.root","UPDATE");
   // cGaus[0]->SaveAs(foldername+"test1.pdf");
   // cGaus[0]->SaveAs(foldername+"test1.pdf");
   if(!file->FindObjectAny("PAP"))
   {
   file->mkdir("PAP");
   }
  if(!file->FindObjectAny("PAL+APL"))
  {
   file->mkdir("PAL+APL");
  }
  if(!file->FindObjectAny("LAL"))
  {
   file->mkdir("LAL");
  }
   for(unsigned i=0;i<3;i++){
     TCanvas* cCheck[i];
     cCheck[i] = new TCanvas();
   ratioCFDataMC[i]->SetTitle(" ; k* [GeV/c]; C(k*)_{data}/C(k*)_{MC}");
   ratioCFDataMC[i]->GetYaxis()->SetRangeUser(0.,6);
   ratioCFDataMC[i]->Draw();

   if(i==0){
     file->cd("PAP");
     if(strcmp(STflag, add1)==0)
     {
       auto* leg= new TLegend(0.25,0.65,0.45,0.8);
       leg->AddEntry(ratioCFDataMC[0], "0 < s_{T} < 0.3","");//"l" sets the legend as lines
       leg->Draw("same");
       cCheck[0]->SaveAs(foldernameplot+"CFDataMC_pAp_st1.pdf");
       cGaus[0]->SaveAs(foldernameplot+"FitCFDataMC_pAp_st1.pdf");


       ratioCFDataMC[0]->SetName("CFDataMC_pAp_st1");
       ratioCFDataMC[0]->Write("CFDataMC_pAp_st1", TObject::kOverwrite);

       diffCFGaus[0]->SetName("CFDataMCGaus_pAp_st1");
       diffCFGaus[0]->Write("CFDataMCGaus_pAp_st1", TObject::kOverwrite);

       hCFmc[0]->SetName("CFMCRaw_pAp_st1");
       hCFmc[0]->Write("CFMCRaw_pAp_st1", TObject::kOverwrite);

       fBackgroundGaus[0]->SetName("FitFuncGausMC_pAp_st1");
       fBackgroundGaus[0]->Write("FitFuncGausMC_pAp_st1", TObject::kOverwrite);
       // file->Write();
       diffCFExp[0]->SetName("CFDataMCExp_pAp_st1");
       diffCFExp[0]->Write("CFDataMCExp_pAp_st1", TObject::kOverwrite);
       // file->Write();
       diffCFPol[0]->SetName("CFDataMCPol_pAp_st1");
       diffCFPol[0]->Write("CFDataMCPol_pAp_st1", TObject::kOverwrite);
       fBackgroundPol[0]->SetName("FitFuncPol2MC_pAp_st1");
       fBackgroundPol[0]->Write("FitFuncPol2MC_pAp_st1", TObject::kOverwrite);
       // file->Write();
       diffCFPol3[0]->SetName("CFDataMCPol3_pAp_st1");
       diffCFPol3[0]->Write("CFDataMCPol3_pAp_st1", TObject::kOverwrite);
       fBackgroundPol3[0]->SetName("FitFuncPol3MC_pAp_st1");
       fBackgroundPol3[0]->Write("FitFuncPol3MC_pAp_st1", TObject::kOverwrite);
       file->Write();
     }
     else if(strcmp(STflag, add2)==0)
     {
       auto* leg= new TLegend(0.25,0.65,0.45,0.8);
       leg->AddEntry(ratioCFDataMC[0], "0.3 < s_{T} < 0.7","");//"l" sets the legend as lines
       leg->Draw("same");
       cCheck[0]->SaveAs(foldernameplot+"CFDataMC_pAp_st2.pdf");
       cGaus[0]->SaveAs(foldernameplot+"FitCFDataMC_pAp_st2.pdf");

       ratioCFDataMC[0]->SetName("CFDataMC_pAp_st2");
       ratioCFDataMC[0]->Write("CFDataMC_pAp_st2", TObject::kOverwrite);

       hCFmc[0]->SetName("CFMCRaw_pAp_st2");
       hCFmc[0]->Write("CFMCRaw_pAp_st2", TObject::kOverwrite);

       diffCFGaus[0]->SetName("CFDataMCGaus_pAp_st2");
       diffCFGaus[0]->Write("CFDataMCGaus_pAp_st2", TObject::kOverwrite);
       fBackgroundGaus[0]->SetName("FitFuncGausMC_pAp_st2");
       fBackgroundGaus[0]->Write("FitFuncGausMC_pAp_st2", TObject::kOverwrite);
       // file->Write();
       diffCFExp[0]->SetName("CFDataMCExp_pAp_st2");
       diffCFExp[0]->Write("CFDataMCExp_pAp_st2", TObject::kOverwrite);
       // file->Write();
       diffCFPol[0]->SetName("CFDataMCPol_pAp_st2");
       diffCFPol[0]->Write("CFDataMCPol_pAp_st2", TObject::kOverwrite);
       fBackgroundPol[0]->SetName("FitFuncPol2MC_pAp_st2");
       fBackgroundPol[0]->Write("FitFuncPol2MC_pAp_st2", TObject::kOverwrite);
       // file->Write();
       diffCFPol3[0]->SetName("CFDataMCPol3_pAp_st2");
       diffCFPol3[0]->Write("CFDataMCPol3_pAp_st2", TObject::kOverwrite);
       fBackgroundPol3[0]->SetName("FitFuncPol3MC_pAp_st2");
       fBackgroundPol3[0]->Write("FitFuncPol3MC_pAp_st2", TObject::kOverwrite);
       file->Write();
     }
     else if(strcmp(STflag, add3)==0)
     {
       auto* leg= new TLegend(0.25,0.65,0.45,0.8);
       leg->AddEntry(ratioCFDataMC[0], "0.7 < s_{T} < 1","");//"l" sets the legend as lines
       leg->Draw("same");
       cCheck[0]->SaveAs(foldernameplot+"CFDataMC_pAp_st3.pdf");
       cGaus[0]->SaveAs(foldernameplot+"FitCFDataMC_pAp_st3.pdf");
       cRat[0]->SaveAs(foldernameplot+"RatioFits_pAp_st3.pdf");


       ratioCFDataMC[0]->SetName("CFDataMC_pAp_st3");
       ratioCFDataMC[0]->Write("CFDataMC_pAp_st3", TObject::kOverwrite);

       diffCFGaus[0]->SetName("CFDataMCGaus_pAp_st3");
       diffCFGaus[0]->Write("CFDataMCGaus_pAp_st3", TObject::kOverwrite);

       hCFmc[0]->SetName("CFMCRaw_pAp_st3");
       hCFmc[0]->Write("CFMCRaw_pAp_st3", TObject::kOverwrite);

       fBackgroundGaus[0]->SetName("FitFuncGausMC_pAp_st3");
       fBackgroundGaus[0]->Write("FitFuncGausMC_pAp_st3", TObject::kOverwrite);
       // file->Write();
       diffCFExp[0]->SetName("CFDataMCExp_pAp_st3");
       diffCFExp[0]->Write("CFDataMCExp_pAp_st3", TObject::kOverwrite);
       // file->Write();
       diffCFPol[0]->SetName("CFDataMCPol_pAp_st3");
       diffCFPol[0]->Write("CFDataMCPol_pAp_st3", TObject::kOverwrite);
       fBackgroundPol[0]->SetName("FitFuncPol2MC_pAp_st3");
       fBackgroundPol[0]->Write("FitFuncPol2MC_pAp_st3", TObject::kOverwrite);
       // file->Write();
       diffCFPol3[0]->SetName("CFDataMCPol3_pAp_st3");
       diffCFPol3[0]->Write("CFDataMCPol3_pAp_st3", TObject::kOverwrite);
       fBackgroundPol3[0]->SetName("FitFuncPol3MC_pAp_st3");
       fBackgroundPol3[0]->Write("FitFuncPol3MC_pAp_st3", TObject::kOverwrite);
       file->Write();
     }
     else if(strcmp(STflag, add4)==0)
     {
       auto* leg= new TLegend(0.25,0.65,0.45,0.8);
       leg->AddEntry(ratioCFDataMC[0], "0 < s_{T} < 1","");//"l" sets the legend as lines
       leg->Draw("same");
       cCheck[0]->SaveAs(foldernameplot+"CFDataMC_pAp_full.pdf");
       cGaus[0]->SaveAs(foldernameplot+"FitCFDataMC_pAp_full.pdf");
       cRat[0]->SaveAs(foldernameplot+"RatioFits_pAp_full.pdf");


       ratioCFDataMC[0]->SetName("CFDataMC_pAp_full");
       ratioCFDataMC[0]->Write("CFDataMC_pAp_full", TObject::kOverwrite);

       diffCFGaus[0]->SetName("CFDataMCGaus_pAp_full");
       diffCFGaus[0]->Write("CFDataMCGaus_pAp_full", TObject::kOverwrite);

       hCFmc[0]->SetName("CFMCRaw_pAp_full");
       hCFmc[0]->Write("CFMCRaw_pAp_full", TObject::kOverwrite);

       fBackgroundGaus[0]->SetName("FitFuncGausMC_pAp_full");
       fBackgroundGaus[0]->Write("FitFuncGausMC_pAp_full", TObject::kOverwrite);

       diffCFExp[0]->SetName("CFDataMCExp_pAp_full");
       diffCFExp[0]->Write("CFDataMCExp_pAp_full", TObject::kOverwrite);
       // file->Write();
       diffCFPol[0]->SetName("CFDataMCPol_pAp_full");
       diffCFPol[0]->Write("CFDataMCPol_pAp_full", TObject::kOverwrite);
       fBackgroundPol[0]->SetName("FitFuncPol2MC_pAp_full");
       fBackgroundPol[0]->Write("FitFuncPol2MC_pAp_full", TObject::kOverwrite);
       // file->Write();
       diffCFPol3[0]->SetName("CFDataMCPol3_pAp_full");
       diffCFPol3[0]->Write("CFDataMCPol3_pAp_full", TObject::kOverwrite);
       fBackgroundPol3[0]->SetName("FitFuncPol3MC_pAp_full");
       fBackgroundPol3[0]->Write("FitFuncPol3MC_pAp_full", TObject::kOverwrite);
       file->Write();
     }
 }

 if (i==1) {
   file->cd("PAL+APL");

   if(strcmp(STflag, add1)==0)
   {
     auto* leg= new TLegend(0.25,0.65,0.45,0.8);
     leg->AddEntry(ratioCFDataMC[1], "0 < s_{T} < 0.3","");//"l" sets the legend as lines
     leg->Draw("same");
     cCheck[1]->SaveAs(foldernameplot+"CFDataMC_pAL_ApL_st1.pdf");
     cGaus[1]->SaveAs(foldernameplot+"FitCFDataMC_pAL_ApL_st1.pdf");


     ratioCFDataMC[1]->SetName("CFDataMC_pAL_ApL_st1");
     ratioCFDataMC[1]->Write("CFDataMC_pAL_ApL_st1", TObject::kOverwrite);
     // file->Write();
      diffCFGaus[1]->SetName("CFDataMCGaus_pAL_ApL_st1");
      diffCFGaus[1]->Write("CFDataMCGaus_pAL_ApL_st1", TObject::kOverwrite);

      hCFmc[1]->SetName("CFMCRaw_pAL_ApL_st1");
      hCFmc[1]->Write("CFMCRaw_pAL_ApL_st1", TObject::kOverwrite);

      fBackgroundGaus[1]->SetName("FitFuncGausMC_pAL_ApL_st1");
      fBackgroundGaus[1]->Write("FitFuncGausMC_pAL_ApL_st1", TObject::kOverwrite);
     // file->Write();
     diffCFExp[1]->SetName("CFDataMCExp_pAL_ApL_st1");
     diffCFExp[1]->Write("CFDataMCExp_pAL_ApL_st1", TObject::kOverwrite);
     // file->Write();
     diffCFPol[1]->SetName("CFDataMCPol_pAL_ApL_st1");
     diffCFPol[1]->Write("CFDataMCPol_pAL_ApL_st1", TObject::kOverwrite);
     fBackgroundPol[1]->SetName("FitFuncPolMC_pAL_ApL_st1");
     fBackgroundPol[1]->Write("FitFuncPolMC_pAL_ApL_st1", TObject::kOverwrite);
     // file->Write();
     diffCFPol3[1]->SetName("CFDataMCPol3_pAL_ApL_st1");
     diffCFPol3[1]->Write("CFDataMCPol3_pAL_ApL_st1", TObject::kOverwrite);
     fBackgroundPol3[1]->SetName("FitFuncPol3MC_pAL_ApL_st1");
     fBackgroundPol3[1]->Write("FitFuncPol3MC_pAL_ApL_st1", TObject::kOverwrite);
     file->Write();

   }
   else if(strcmp(STflag, add2)==0)
   {
     auto* leg= new TLegend(0.25,0.65,0.45,0.8);
     leg->AddEntry(ratioCFDataMC[1], "0.3 < s_{T} < 0.7","");//"l" sets the legend as lines
     leg->Draw("same");
     cCheck[1]->SaveAs(foldernameplot+"CFDataMC_pAL_ApL_st2.pdf");
     cGaus[1]->SaveAs(foldernameplot+"FitCFDataMC_pAL_ApL_st2.pdf");

     ratioCFDataMC[1]->SetName("CFDataMC_pAL_ApL_st2");
     ratioCFDataMC[1]->Write("CFDataMC_pAL_ApL_st2", TObject::kOverwrite);
     // file->Write();
     diffCFGaus[1]->SetName("CFDataMCGaus_pAL_ApL_st2");
     diffCFGaus[1]->Write("CFDataMCGaus_pAL_ApL_st2", TObject::kOverwrite);

     hCFmc[1]->SetName("CFMCRaw_pAL_ApL_st2");
     hCFmc[1]->Write("CFMCRaw_pAL_ApL_st2", TObject::kOverwrite);

     fBackgroundGaus[1]->SetName("FitFuncGausMC_pAL_ApL_st2");
     fBackgroundGaus[1]->Write("FitFuncGausMC_pAL_ApL_st2", TObject::kOverwrite);
     // file->Write();
     diffCFExp[1]->SetName("CFDataMCExp_pAL_ApL_st2");
     diffCFExp[1]->Write("CFDataMCExp_pAL_ApL_st2", TObject::kOverwrite);
     // file->Write();
     diffCFPol[1]->SetName("CFDataMCPol_pAL_ApL_st2");
     diffCFPol[1]->Write("CFDataMCPol_pAL_ApL_st2", TObject::kOverwrite);
     fBackgroundPol[1]->SetName("FitFuncPolMC_pAL_ApL_st2");
     fBackgroundPol[1]->Write("FitFuncPolMC_pAL_ApL_st2", TObject::kOverwrite);
     // file->Write();
     diffCFPol3[1]->SetName("CFDataMCPol3_pAL_ApL_st2");
     diffCFPol3[1]->Write("CFDataMCPol3_pAL_ApL_st2", TObject::kOverwrite);
     fBackgroundPol3[1]->SetName("FitFuncPol3MC_pAL_ApL_st2");
     fBackgroundPol3[1]->Write("FitFuncPol3MC_pAL_ApL_st2", TObject::kOverwrite);
     file->Write();
   }
   else if(strcmp(STflag, add3)==0)
   {
     auto* leg= new TLegend(0.25,0.65,0.45,0.8);
     leg->AddEntry(ratioCFDataMC[1], "0.7 < s_{T} < 1","");//"l" sets the legend as lines
     leg->Draw("same");
     cCheck[1]->SaveAs(foldernameplot+"CFDataMC_pAL_ApL_st3.pdf");
     cGaus[1]->SaveAs(foldernameplot+"FitCFDataMC_pAL_ApL_st3.pdf");
     cRat[1]->SaveAs(foldernameplot+"RatioFits_pAL_ApL_st3.pdf");

     ratioCFDataMC[1]->SetName("CFDataMC_pAL_ApL_st3");
     ratioCFDataMC[1]->Write("CFDataMC_pAL_ApL_st3", TObject::kOverwrite);
     // file->Write();
     diffCFGaus[1]->SetName("CFDataMCGaus_pAL_ApL_st3");
     diffCFGaus[1]->Write("CFDataMCGaus_pAL_ApL_st3", TObject::kOverwrite);

     hCFmc[1]->SetName("CFMCRaw_pAL_ApL_st3");
     hCFmc[1]->Write("CFMCRaw_pAL_ApL_st3", TObject::kOverwrite);

     fBackgroundGaus[1]->SetName("FitFuncGausMC_pAL_ApL_st3");
     fBackgroundGaus[1]->Write("FitFuncGausMC_pAL_ApL_st3", TObject::kOverwrite);
     // file->Write();
     diffCFExp[1]->SetName("CFDataMCExp_pAL_ApL_st3");
     diffCFExp[1]->Write("CFDataMCExp_pAL_ApL_st3", TObject::kOverwrite);
     // file->Write();
     diffCFPol[1]->SetName("CFDataMCPol_pAL_ApL_st3");
     diffCFPol[1]->Write("CFDataMCPol_pAL_ApL_st3", TObject::kOverwrite);
     fBackgroundPol[1]->SetName("FitFuncPolMC_pAL_ApL_st3");
     fBackgroundPol[1]->Write("FitFuncPolMC_pAL_ApL_st3", TObject::kOverwrite);
     // file->Write();
     diffCFPol3[1]->SetName("CFDataMCPol3_pAL_ApL_st3");
     diffCFPol3[1]->Write("CFDataMCPol3_pAL_ApL_st3", TObject::kOverwrite);
     fBackgroundPol3[1]->SetName("FitFuncPol3MC_pAL_ApL_st3");
     fBackgroundPol3[1]->Write("FitFuncPol3MC_pAL_ApL_st3", TObject::kOverwrite);
     file->Write();
   }
   else if(strcmp(STflag, add4)==0)
   {
     auto* leg= new TLegend(0.25,0.65,0.45,0.8);
     leg->AddEntry(ratioCFDataMC[1], "0 < s_{T} < 1","");//"l" sets the legend as lines
     leg->Draw("same");
     cCheck[1]->SaveAs(foldernameplot+"CFDataMC_pAL_ApL_full.pdf");
     cGaus[1]->SaveAs(foldernameplot+"FitCFDataMC_pAL_ApL_full.pdf");
     cRat[1]->SaveAs(foldernameplot+"RatioFits_pAL_ApL_full.pdf");


     ratioCFDataMC[1]->SetName("CFDataMC_pAL_ApL_full");
     ratioCFDataMC[1]->Write("CFDataMC_pAL_ApL_full", TObject::kOverwrite);
     // file->Write();
     diffCFGaus[1]->SetName("CFDataMCGaus_pAL_ApL_full");
     diffCFGaus[1]->Write("CFDataMCGaus_pAL_ApL_full", TObject::kOverwrite);

     hCFmc[1]->SetName("CFMCRaw_pAL_ApL_full");
     hCFmc[1]->Write("CFMCRaw_pAL_ApL_full", TObject::kOverwrite);

     fBackgroundGaus[1]->SetName("FitFuncGausMC_pAL_ApL_full");
     fBackgroundGaus[1]->Write("FitFuncGausMC_pAL_ApL_full", TObject::kOverwrite);
     // file->Write();
     diffCFExp[1]->SetName("CFDataMCExp_pAL_ApL_full");
     diffCFExp[1]->Write("CFDataMCExp_pAL_ApL_full", TObject::kOverwrite);
     // file->Write();
     diffCFPol[1]->SetName("CFDataMCPol_pAL_ApL_full");
     diffCFPol[1]->Write("CFDataMCPol_pAL_ApL_full", TObject::kOverwrite);
     fBackgroundPol[1]->SetName("FitFuncPolMC_pAL_ApL_full");
     fBackgroundPol[1]->Write("FitFuncPolMC_pAL_ApL_full", TObject::kOverwrite);
     // file->Write();
     diffCFPol3[1]->SetName("CFDataMCPol3_pAL_ApL_full");
     diffCFPol3[1]->Write("CFDataMCPol3_pAL_ApL_full", TObject::kOverwrite);
     fBackgroundPol3[1]->SetName("FitFuncPol3MC_pAL_ApL_full");
     fBackgroundPol3[1]->Write("FitFuncPol3MC_pAL_ApL_full", TObject::kOverwrite);
     file->Write();
   }
}
  if(i==2){
    file->cd("LAL");

    if(strcmp(STflag, add1)==0)
    {
      auto* leg= new TLegend(0.25,0.65,0.45,0.8);
      leg->AddEntry(ratioCFDataMC[2], "0 < s_{T} < 0.3","");//"l" sets the legend as lines
      leg->Draw("same");
      cCheck[2]->SaveAs(foldernameplot+"CFDataMC_LAL_st1.pdf");
      cGaus[2]->SaveAs(foldernameplot+"FitCFDataMC_LAL_st1.pdf");

      ratioCFDataMC[2]->SetName("CFDataMC_LAL_st1");
      ratioCFDataMC[2]->Write("CFDataMC_LAL_st1", TObject::kOverwrite);
      // file->Write();
      diffCFGaus[2]->SetName("CFDataMCGaus_LAL_st1");
      diffCFGaus[2]->Write("CFDataMCGaus_LAL_st1", TObject::kOverwrite);

      hCFmc[2]->SetName("CFMCRaw_LAL_st1");
      hCFmc[2]->Write("CFMCRaw_LAL_st1", TObject::kOverwrite);

      fBackgroundGaus[2]->SetName("FitFuncGausMC_LAL_st1");
      fBackgroundGaus[2]->Write("FitFuncGausMC_LAL_st1", TObject::kOverwrite);
      // file->Write();
      diffCFExp[2]->SetName("CFDataMCExp_LAL_st1");
      diffCFExp[2]->Write("CFDataMCExp_LAL_st1", TObject::kOverwrite);
      // file->Write();
      diffCFPol[2]->SetName("CFDataMCPol_LAL_st1");
      diffCFPol[2]->Write("CFDataMCPol_LAL_st1", TObject::kOverwrite);
      fBackgroundPol[2]->SetName("FitFuncPolMC_LAL_st1");
      fBackgroundPol[2]->Write("FitFuncPolMC_LAL_st1", TObject::kOverwrite);
      // file->Write();
      diffCFPol3[2]->SetName("CFDataMCPol3_LAL_st1");
      diffCFPol3[2]->Write("CFDataMCPol3_LAL_st1", TObject::kOverwrite);
      fBackgroundPol3[2]->SetName("FitFuncPol3MC_LAL_st1");
      fBackgroundPol3[2]->Write("FitFuncPol3MC_LAL_st1", TObject::kOverwrite);
      file->Write();

    }
    else if(strcmp(STflag, add2)==0)
    {
      auto* leg= new TLegend(0.25,0.65,0.45,0.8);
      leg->AddEntry(ratioCFDataMC[2], "0.3 < s_{T} < 0.7","");//"l" sets the legend as lines
      leg->Draw("same");
      cCheck[2]->SaveAs(foldernameplot+"CFDataMC_LAL_st2.pdf");
      cGaus[2]->SaveAs(foldernameplot+"FitCFDataMC_LAL_st2.pdf");

      ratioCFDataMC[2]->SetName("CFDataMC_LAL_st2");
      ratioCFDataMC[2]->Write("CFDataMC_LAL_st2", TObject::kOverwrite);
      // file->Write();
      diffCFGaus[2]->SetName("CFDataMCGaus_LAL_st2");
      diffCFGaus[2]->Write("CFDataMCGaus_LAL_st2", TObject::kOverwrite);

      hCFmc[2]->SetName("CFMCRaw_LAL_st2");
      hCFmc[2]->Write("CFMCRaw_LAL_st2", TObject::kOverwrite);

      fBackgroundGaus[2]->SetName("FitFuncGausMC_LAL_st2");
      fBackgroundGaus[2]->Write("FitFuncGausMC_LAL_st2", TObject::kOverwrite);
      // file->Write();
      diffCFExp[2]->SetName("CFDataMCExp_LAL_st2");
      diffCFExp[2]->Write("CFDataMCExp_LAL_st2", TObject::kOverwrite);
      // file->Write();
      diffCFPol[2]->SetName("CFDataMCPol_LAL_st2");
      diffCFPol[2]->Write("CFDataMCPol_LAL_st2", TObject::kOverwrite);
      fBackgroundPol[2]->SetName("FitFuncPolMC_LAL_st2");
      fBackgroundPol[2]->Write("FitFuncPolMC_LAL_st2", TObject::kOverwrite);
      // file->Write();
      diffCFPol3[2]->SetName("CFDataMCPol3_LAL_st2");
      diffCFPol3[2]->Write("CFDataMCPol3_LAL_st2", TObject::kOverwrite);
      fBackgroundPol3[2]->SetName("FitFuncPol3MC_LAL_st2");
      fBackgroundPol3[2]->Write("FitFuncPol3MC_LAL_st2", TObject::kOverwrite);
      file->Write();
    }
    else if(strcmp(STflag, add3)==0)
    {
      auto* leg= new TLegend(0.25,0.65,0.45,0.8);
      leg->AddEntry(ratioCFDataMC[2], "0.7 < s_{T} < 1","");//"l" sets the legend as lines
      leg->Draw("same");
      cCheck[2]->SaveAs(foldernameplot+"CFDataMC_LAL_st3.pdf");
      cGaus[2]->SaveAs(foldernameplot+"FitCFDataMC_LAL_st3.pdf");
      cRat[2]->SaveAs(foldernameplot+"RatioFits_LAL_st3.pdf");

      ratioCFDataMC[2]->SetName("CFDataMC_LAL_st3");
      ratioCFDataMC[2]->Write("CFDataMC_LAL_st3", TObject::kOverwrite);
      // file->Write();
      diffCFGaus[2]->SetName("CFDataMCGaus_LAL_st3");
      diffCFGaus[2]->Write("CFDataMCGaus_LAL_st3", TObject::kOverwrite);

      hCFmc[2]->SetName("CFMCRaw_LAL_st3");
      hCFmc[2]->Write("CFMCRaw_LAL_st3", TObject::kOverwrite);

      fBackgroundGaus[2]->SetName("FitFuncGausMC_LAL_st3");
      fBackgroundGaus[2]->Write("FitFuncGausMC_LAL_st3", TObject::kOverwrite);
      // file->Write();
      diffCFExp[2]->SetName("CFDataMCExp_LAL_st3");
      diffCFExp[2]->Write("CFDataMCExp_LAL_st3", TObject::kOverwrite);
      // file->Write();
      diffCFPol[2]->SetName("CFDataMCPol_LAL_st3");
      diffCFPol[2]->Write("CFDataMCPol_LAL_st3", TObject::kOverwrite);
      fBackgroundPol[2]->SetName("FitFuncPolMC_LAL_st3");
      fBackgroundPol[2]->Write("FitFuncPolMC_LAL_st3", TObject::kOverwrite);
      // file->Write();
      diffCFPol3[2]->SetName("CFDataMCPol3_LAL_st3");
      diffCFPol3[2]->Write("CFDataMCPol3_LAL_st3", TObject::kOverwrite);
      fBackgroundPol3[2]->SetName("FitFuncPol3MC_LAL_st3");
      fBackgroundPol3[2]->Write("FitFuncPol3MC_LAL_st3", TObject::kOverwrite);
      file->Write();
    }
    else if(strcmp(STflag, add4)==0)
    {
      auto* leg= new TLegend(0.25,0.65,0.45,0.8);
      leg->AddEntry(ratioCFDataMC[2], "0 < s_{T} < 1","");//"l" sets the legend as lines
      leg->Draw("same");
      cCheck[2]->SaveAs(foldernameplot+"CFDataMC_LAL_full.pdf");
      cGaus[2]->SaveAs(foldernameplot+"FitCFDataMC_LAL_full.pdf");
      cRat[2]->SaveAs(foldernameplot+"RatioFits_LAL_full.pdf");

      ratioCFDataMC[2]->SetName("CFDataMC_LAL_full");
      ratioCFDataMC[2]->Write("CFDataMC_LAL_full", TObject::kOverwrite);
      // file->Write();
      diffCFGaus[2]->SetName("CFDataMCGaus_LAL_full");
      diffCFGaus[2]->Write("CFDataMCGaus_LAL_full", TObject::kOverwrite);

      hCFmc[2]->SetName("CFMCRaw_LAL_full");
      hCFmc[2]->Write("CFMCRaw_LAL_full", TObject::kOverwrite);

      fBackgroundGaus[2]->SetName("FitFuncGausMC_LAL_full");
      fBackgroundGaus[2]->Write("FitFuncGausMC_LAL_full", TObject::kOverwrite);
      // file->Write();
      diffCFExp[2]->SetName("CFDataMCExp_LAL_full");
      diffCFExp[2]->Write("CFDataMCExp_LAL_full", TObject::kOverwrite);
      // file->Write();
      diffCFPol[2]->SetName("CFDataMCPol_LAL_full");
      diffCFPol[2]->Write("CFDataMCPol_LAL_full", TObject::kOverwrite);
      fBackgroundPol[2]->SetName("FitFuncPolMC_LAL_full");
      fBackgroundPol[2]->Write("FitFuncPolMC_LAL_full", TObject::kOverwrite);
      // file->Write();
      diffCFPol3[2]->SetName("CFDataMCPol3_LAL_full");
      diffCFPol3[2]->Write("CFDataMCPol3_LAL_full", TObject::kOverwrite);
      fBackgroundPol3[2]->SetName("FitFuncPol3MC_LAL_full");
      fBackgroundPol3[2]->Write("FitFuncPol3MC_LAL_full", TObject::kOverwrite);
      file->Write();
    }
  }
   }


   for(unsigned j=0;j<4;j++){
   TCanvas* cMECheck[j];
   cMECheck[j] = new TCanvas();
   ratioMEDataMC[j]->SetTitle(" ; k* [GeV/c]; ME(k*)_{data}/ ME(k*)_{MC}");
   ratioMEDataMC[j]->Draw();
   if(j==0)
  {
     if(strcmp(STflag, add1)==0)
     {
       auto* leg= new TLegend(0.25,0.65,0.45,0.8);
       leg->AddEntry(ratioMEDataMC[0], "0 < s_{T} < 0.3","");//"l" sets the legend as lines
       leg->Draw("same");
       cMECheck[0]->SaveAs(foldernameplot+"MEDataMC_pAp_st1.pdf");
     }
     else if(strcmp(STflag, add2)==0)
     {
       auto* leg= new TLegend(0.25,0.65,0.45,0.8);
       leg->AddEntry(ratioMEDataMC[0], "0.3 < s_{T} < 0.7","");//"l" sets the legend as lines
       leg->Draw("same");
       cMECheck[0]->SaveAs(foldernameplot+"MEDataMC_pAp_st2.pdf");
     }
     else if(strcmp(STflag, add3)==0)
     {
       auto* leg= new TLegend(0.25,0.65,0.45,0.8);
       leg->AddEntry(ratioMEDataMC[0], "0.7 < s_{T} < 1","");//"l" sets the legend as lines
       leg->Draw("same");
       cMECheck[0]->SaveAs(foldernameplot+"MEDataMC_pAp_st3.pdf");
     }
     else if(strcmp(STflag, add4)==0)
     {
       auto* leg= new TLegend(0.25,0.2,0.45,0.4);
       leg->AddEntry(ratioMEDataMC[0], "0 < s_{T} < 1","");//"l" sets the legend as lines
       leg->Draw("same");
       cMECheck[0]->SaveAs(foldernameplot+"MEDataMC_pAp_full.pdf");
     }
  }
 if (j==1)
 {

   if(strcmp(STflag, add1)==0)
   {
     auto* leg= new TLegend(0.25,0.65,0.45,0.8);
     leg->AddEntry(ratioMEDataMC[1], "0 < s_{T} < 0.3","");//"l" sets the legend as lines
     leg->Draw("same");
     cMECheck[1]->SaveAs(foldernameplot+"MEDataMC_pAL_st1.pdf");
   }
   else if(strcmp(STflag, add2)==0)
   {
     auto* leg= new TLegend(0.25,0.65,0.45,0.8);
     leg->AddEntry(ratioMEDataMC[1], "0.3 < s_{T} < 0.7","");//"l" sets the legend as lines
     leg->Draw("same");
     cMECheck[1]->SaveAs(foldernameplot+"MEDataMC_pAL_st2.pdf");
   }
   else if(strcmp(STflag, add3)==0)
   {
     auto* leg= new TLegend(0.25,0.65,0.45,0.8);
     leg->AddEntry(ratioMEDataMC[1], "0.7 < s_{T} < 1","");//"l" sets the legend as lines
     leg->Draw("same");
     cMECheck[1]->SaveAs(foldernameplot+"MEDataMC_pAL_st3.pdf");
   }
   else if(strcmp(STflag, add4)==0)
   {
     auto* leg= new TLegend(0.25,0.65,0.45,0.8);
     leg->AddEntry(ratioMEDataMC[1], "0 < s_{T} < 1","");//"l" sets the legend as lines
     leg->Draw("same");
     cMECheck[1]->SaveAs(foldernameplot+"MEDataMC_pAL_full.pdf");
   }
 }
   if(j==2)
 {

     if(strcmp(STflag, add1)==0)
     {
       auto* leg= new TLegend(0.25,0.65,0.45,0.8);
       leg->AddEntry(ratioMEDataMC[2], "0 < s_{T} < 0.3","");//"l" sets the legend as lines
       leg->Draw("same");
       cMECheck[2]->SaveAs(foldernameplot+"MEDataMC_LAL_st1.pdf");
     }
     else if(strcmp(STflag, add2)==0)
     {
       auto* leg= new TLegend(0.25,0.65,0.45,0.8);
       leg->AddEntry(ratioMEDataMC[2], "0.3 < s_{T} < 0.7","");//"l" sets the legend as lines
       leg->Draw("same");
       cMECheck[2]->SaveAs(foldernameplot+"MEDataMC_LAL_st2.pdf");
     }
     else if(strcmp(STflag, add3)==0)
     {
       auto* leg= new TLegend(0.25,0.65,0.45,0.8);
       leg->AddEntry(ratioMEDataMC[2], "0.7 < s_{T} < 1","");//"l" sets the legend as lines
       leg->Draw("same");
       cMECheck[2]->SaveAs(foldernameplot+"MEDataMC_LAL_st3.pdf");
     }
     else if(strcmp(STflag, add4)==0)
     {
       auto* leg= new TLegend(0.25,0.65,0.45,0.8);
       leg->AddEntry(ratioMEDataMC[2], "0 < s_{T} < 1","");//"l" sets the legend as lines
       leg->Draw("same");
       cMECheck[2]->SaveAs(foldernameplot+"MEDataMC_LAL_full.pdf");
     }
 }
   if(j==3)
 {

     if(strcmp(STflag, add1)==0)
     {
       auto* leg= new TLegend(0.25,0.65,0.45,0.8);
       leg->AddEntry(ratioMEDataMC[3], "0 < s_{T} < 0.3","");//"l" sets the legend as lines
       leg->Draw("same");
       cMECheck[3]->SaveAs(foldernameplot+"MEDataMC_apL_st1.pdf");
     }
     else if(strcmp(STflag, add2)==0)
     {
       auto* leg= new TLegend(0.25,0.65,0.45,0.8);
       leg->AddEntry(ratioMEDataMC[3], "0.3 < s_{T} < 0.7","");//"l" sets the legend as lines
       leg->Draw("same");
       cMECheck[3]->SaveAs(foldernameplot+"MEDataMC_apL_st2.pdf");
     }
     else if(strcmp(STflag, add3)==0)
     {
       auto* leg= new TLegend(0.25,0.65,0.45,0.8);
       leg->AddEntry(ratioMEDataMC[3], "0.7 < s_{T} < 1","");//"l" sets the legend as lines
       leg->Draw("same");
       cMECheck[3]->SaveAs(foldernameplot+"MEDataMC_apL_st3.pdf");
     }
     else if(strcmp(STflag, add4)==0)
     {
       auto* leg= new TLegend(0.25,0.65,0.45,0.8);
       leg->AddEntry(ratioMEDataMC[3], "0 < s_{T} < 1","");//"l" sets the legend as lines
       leg->Draw("same");
       cMECheck[3]->SaveAs(foldernameplot+"MEDataMC_apL_full.pdf");
     }
 }
}

  std::cout << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << std::endl;
  std::cout<<"The system only dreams in total darkness"<<std::endl;

}

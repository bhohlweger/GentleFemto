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

int main(int argc, char* argv[]) {
  const char* filename = argv[1];
  const char* prefix = argv[2];
  const char* addon = (argv[3]) ? argv[3] : "";

gStyle->SetOptStat(0);
gStyle->SetOptTitle(kFALSE);

DreamPlot::SetStyle();

//Accessing all the directory in the root file
  TFile* _file0=TFile::Open(filename);

  TList *list1 = 0;

  TDirectoryFile *dir=(TDirectoryFile*)(_file0->FindObjectAny(Form("%sResults%s", prefix, addon)));
  dir->GetObject(Form("%sResults%s",prefix,addon),list1);


  TList* list_p0_p1=(TList*)list1->FindObject("Particle0_Particle1");
  TList* list_p0_p3=(TList*)list1->FindObject("Particle0_Particle3");
  TList* list_p1_p2=(TList*)list1->FindObject("Particle1_Particle2");
  TList* list_p2_p3=(TList*)list1->FindObject("Particle2_Particle3");

  auto* SEmT_p0_p1_2D = (TH2F*)list_p0_p1->FindObject("SEmTDist_Particle0_Particle1");
  auto* SEmT_p0_p3_2D = (TH2F*)list_p0_p3->FindObject("SEmTDist_Particle0_Particle3");
  auto* SEmT_p1_p2_2D = (TH2F*)list_p1_p2->FindObject("SEmTDist_Particle1_Particle2");
  auto* SEmT_p2_p3_2D = (TH2F*)list_p2_p3->FindObject("SEmTDist_Particle2_Particle3");



//run up to k* bin 50, so up to k*=200 MeV
int Num_kstar_bins = 50;
TH1F* SEmT_p0_p1;float mean_p0_p1[Num_kstar_bins];
TH1F* SEmT_p0_p3;float mean_p0_p3[Num_kstar_bins];
TH1F* SEmT_p1_p2;float mean_p1_p2[Num_kstar_bins];
TH1F* SEmT_p2_p3;float mean_p2_p3[Num_kstar_bins];


for(int i_k = 1; i_k < Num_kstar_bins+1; i_k++)
{
  SEmT_p0_p1 = (TH1F*) SEmT_p0_p1_2D->ProjectionY(TString::Format("SEmT_p0_p1_%i",i_k),i_k,i_k,"e")->Clone(Form("SEmT_p0_p1_%i",i_k));
  mean_p0_p1[i_k] = SEmT_p0_p1->GetMean(1);
  printf("# bins in k* = %i\n", SEmT_p0_p1->GetNbinsX());
  Double_t binCenter = SEmT_p0_p1->GetXaxis()->GetBinCenter(i_k);
  printf("Mean for p-antip\n");
  std::cout<<"bin in k* = "<<i_k<< "value of k* = "<<binCenter<<"<m_T> (p-antip) = "<<mean_p0_p1[i_k]<<std::endl;

  SEmT_p0_p3 = (TH1F*) SEmT_p0_p3_2D->ProjectionY(TString::Format("SEmT_p0_p3_%i",i_k),i_k,i_k,"e")->Clone(Form("SEmT_p0_p3_%i",i_k));
  mean_p0_p3[i_k] = SEmT_p0_p3->GetMean(1);
  SEmT_p1_p2 = (TH1F*) SEmT_p1_p2_2D->ProjectionY(TString::Format("SEmT_p1_p2_%i",i_k),i_k,i_k,"e")->Clone(Form("SEmT_p1_p2_%i",i_k));
  mean_p1_p2[i_k] = SEmT_p1_p2->GetMean(1);
  SEmT_p2_p3 = (TH1F*) SEmT_p2_p3_2D->ProjectionY(TString::Format("SEmT_p2_p3_%i",i_k),i_k,i_k,"e")->Clone(Form("SEmT_p2_p3_%i",i_k));
  mean_p2_p3[i_k] = SEmT_p0_p1->GetMean(1);

  //h2->Fill(x,y);
}

return 0;
}

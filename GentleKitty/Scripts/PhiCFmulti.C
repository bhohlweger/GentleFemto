#include <iostream>
#include "CATS.h"
#include "CATSInputSigma0.h"
#include "CATSLambdaParam.h"
#include "DLM_CkDecomposition.h"
#include "DLM_CkModels.h"
#include "DLM_Potentials.h"
#include "DLM_Source.h"
#include "DreamPlot.h"
#include "SidebandSigma.h"
#include "TCanvas.h"
#include "TDatabasePDG.h"
#include "TGenPhaseSpace.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TNtuple.h"
#include "TROOT.h"
#include "TRandom3.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TidyCats.h"

/// =====================================================================================
void FillWaveGraph(CATS& kitty, TGraph* gr) {
  for (unsigned int i = 0; i < kitty.GetNumMomBins(); ++i) {
    const float mom = kitty.GetMomentum(i);
    gr->SetPoint(i, mom, kitty.GetCorrFun(i));
  }
}

double integrate(TGraph* CF, double min, double max) {
  double sum = 0;
  int diff = max - min;
  int a = 0;
  for (int i = min; i <= max; i++) {
    double ycf = CF->Eval(i);
  //  cout << "value: " << i << " ----> " << ycf << endl;
    sum = sum + ycf;
    a++;
  }
 // cout << sum << endl;
  return sum / a;
 // cout << sum / a << endl;
}

TGraph* GetSmearedCFscaled(TGraph* CF, TH2* matrix, double scale) {
  // Define new Histogram which have dimension according to the yaxis (new
  // momentum axis):
  const int nbins_original = matrix->GetXaxis()->GetNbins();
  const Int_t nbins_transformed = matrix->GetYaxis()->GetNbins();

  TGraph* smearedCF = new TGraph();

  int countPoint = 0;
  for (int momTrans = 0; momTrans < nbins_transformed; momTrans++) {
    Double_t matrixvalues_sum = 0.;
    Double_t weighted_matrixvalues_sum = 0.;

    for (int momOri = 0; momOri < nbins_original; momOri++) {
      Double_t momentum_original = matrix->GetXaxis()->GetBinCenter(momOri + 1);
      matrixvalues_sum += matrix->GetBinContent(momOri + 1, momTrans + 1);
      weighted_matrixvalues_sum +=
          CF->Eval(momentum_original) *
          matrix->GetBinContent(momOri + 1, momTrans + 1);
    }
    Double_t transformed_CF = 0.;
    if (matrixvalues_sum != 0.)
      transformed_CF = weighted_matrixvalues_sum / matrixvalues_sum;

    smearedCF->SetPoint(countPoint++,
                        matrix->GetYaxis()->GetBinCenter(momTrans + 1),
                        transformed_CF * scale);
  }
  return smearedCF;
}

/// =====================================================================================
TGraph* GetSmearedCF(TGraph* CF, TH2* matrix) {
  // Define new Histogram which have dimension according to the yaxis (new
  // momentum axis):
  const int nbins_original = matrix->GetXaxis()->GetNbins();
  const Int_t nbins_transformed = matrix->GetYaxis()->GetNbins();

  TGraph* smearedCF = new TGraph();

  int countPoint = 0;
  for (int momTrans = 0; momTrans < nbins_transformed; momTrans++) {
    Double_t matrixvalues_sum = 0.;
    Double_t weighted_matrixvalues_sum = 0.;

    for (int momOri = 0; momOri < nbins_original; momOri++) {
      Double_t momentum_original = matrix->GetXaxis()->GetBinCenter(momOri + 1);
      matrixvalues_sum += matrix->GetBinContent(momOri + 1, momTrans + 1);
      weighted_matrixvalues_sum +=
          CF->Eval(momentum_original) *
          matrix->GetBinContent(momOri + 1, momTrans + 1);
    }
    Double_t transformed_CF = 0.;
    if (matrixvalues_sum != 0.)
      transformed_CF = weighted_matrixvalues_sum / matrixvalues_sum;

    smearedCF->SetPoint(countPoint++,
                        matrix->GetYaxis()->GetBinCenter(momTrans + 1),
                        transformed_CF);
  }
  return smearedCF;
}

TGraph* GetSmearedCF2(TGraph* CF, TH2* matrix) {
  // Define new Histogram which have dimension according to the yaxis (new
  // momentum axis):
  const int nbins_original = matrix->GetXaxis()->GetNbins();
  const Int_t nbins_transformed = matrix->GetYaxis()->GetNbins();

  TGraph* smearedCF = new TGraph();

  int countPoint = 0;
  for (int momTrans = 0; momTrans < nbins_transformed; momTrans++) {
    Double_t matrixvalues_sum = 0.;
    Double_t weighted_matrixvalues_sum = 0.;

    for (int momOri = 0; momOri < nbins_original; momOri++) {
      Double_t momentum_original = matrix->GetXaxis()->GetBinCenter(momOri + 1);
      matrixvalues_sum += matrix->GetBinContent(momOri + 1, momTrans + 1);
      weighted_matrixvalues_sum +=
          CF->Eval(momentum_original) *
          matrix->GetBinContent(momOri + 1, momTrans + 1);
    }
    Double_t transformed_CF = 0.;
    if (matrixvalues_sum != 0.)
      transformed_CF = weighted_matrixvalues_sum / matrixvalues_sum;

    smearedCF->SetPoint(countPoint++,
                        (matrix->GetYaxis()->GetBinCenter(momTrans + 1)),
                        (transformed_CF));
  }
  return smearedCF;
}

TGraph* GetSmearedCFmoved(TGraph* CF, TH2* matrix, double max_x, double max_y) {
  // Define new Histogram which have dimension according to the yaxis (new
  // momentum axis):

  //  double max_cf_x=CF->GetMaximum();

  const int nbins_original = matrix->GetXaxis()->GetNbins();
  const Int_t nbins_transformed = matrix->GetYaxis()->GetNbins();

  TGraph* smearedCF = new TGraph();

  int countPoint = 0;
  for (int momTrans = 0; momTrans < nbins_transformed; momTrans++) {
    Double_t matrixvalues_sum = 0.;
    Double_t weighted_matrixvalues_sum = 0.;

    for (int momOri = 0; momOri < nbins_original; momOri++) {
      Double_t momentum_original = matrix->GetXaxis()->GetBinCenter(momOri + 1);
      matrixvalues_sum += matrix->GetBinContent(momOri + 1, momTrans + 1);
      weighted_matrixvalues_sum +=
          CF->Eval(momentum_original) *
          matrix->GetBinContent(momOri + 1, momTrans + 1);
    }
    Double_t transformed_CF = 0.;
    if (matrixvalues_sum != 0.)
      transformed_CF = weighted_matrixvalues_sum / matrixvalues_sum;

    smearedCF->SetPoint(
        countPoint++, (matrix->GetYaxis()->GetBinCenter(momTrans + 1) - max_x),
        (transformed_CF - max_y));
  }
  return smearedCF;
}

double getxmax(TGraph* CF) {
  double ymaxcf = 0;
  double xmaxcf = 0;
  for (int i = 0; i < 1000; i++) {
    double ycf = CF->Eval(i);
    // cout<<i<< " yval: "<<ycf<<endl;
    if (ycf >= ymaxcf) {
      ymaxcf = ycf;
      xmaxcf = i;
    }
  }
 // cout << xmaxcf << endl;
  return xmaxcf;
}

double getymax(TGraph* CF) {
  double ymaxcf = 0;
  double xmaxcf = 0;
  for (int i = 0; i < 1000; i++) {
    double ycf = CF->Eval(i);
    // cout<<i<< " yval: "<<ycf<<endl;
    if (ycf >= ymaxcf) {
      ymaxcf = ycf;
      xmaxcf = i;
    }
  }
// cout << ymaxcf << endl;

  return ymaxcf;
}



//double integrate(TGraph* CF, double min, double max){
//double sum=0;
//int diff=max-min;
//int a=0;
//for (int i=min;i<=max;i++){
//double ycf= CF->Eval(i);
//cout<<"value: "<< i<< " ----> "<<ycf<<endl;
//sum = sum+ycf;
//a++;
//}
//cout<<sum<<endl;
//return sum/a;
//cout<<sum/a<<endl;

//}

// void setmaxima(TGraph* CF,double ymaxcf, double xmaxcf ){
// ymaxcf =0;
// xmaxcf=0;
// for (int i=0;i<CF->GetN();i++){
// double ycf= CF->Eval(i);
// cout<<i<< " yval: "<<ycf<<endl;
// while (ycf>=ymaxcf){
//    ymaxcf=ycf;
//    xmaxcf=i;
//}
//}
//}

TH1F* CalcCF(const char* prefix, int a) {
  const char* addon = Form("%d", a);
//  const char* name =
//      Form("/home/emma/FemtoPhiHM_checks2/800900/CFOutput_pPhi_%s_%s.root",
//           prefix, addon);
  const char* name =Form("/home/emma/FemtoPhiHM_allSB/1000-1200/CFOutput_pPhi_%s_%s.root",prefix, addon);
  auto file = TFile::Open(name);
  // file->ls();
  TH1F* corr = (TH1F*)file->FindObjectAny("hCk_ReweightedMeV_2");
  return corr;
}

TH1F* CalcCF2(const char* prefix, int a) {
  const char* addon = Form("%d", a);
//  const char* name =
//      Form("/home/emma/FemtoPhiHM_allSB/800900/CFOutput_pPhi_%s_%s.root",
//           prefix, addon);
  const char* name =Form("/home/emma/FemtoPhiHM_allSB/1000-1200/CFOutput_pPhi_%s_%s.root",prefix, addon);
  auto file = TFile::Open(name);
  // file->ls();
  TH1F* corr = (TH1F*)file->FindObjectAny("hCk_ReweightedMeV_2");
  return corr;
}

// TH2F* GetMatrix( const char* prefix,  int a) {
//  const char* addon=Form("%d",a);
//  const char*
//  name=Form("/home/emma/PhiTRUTH_2/phitruth/CFOutput_pPhi_%s_%s.root",prefix,
//  addon);
//  auto file = TFile::Open(name);
//  //file->ls();
//  TH1F* corr = (TH1F*)file->FindObjectAny("hCk_ReweightedMeV_2");
//  return corr;

//}

double evalMomOri(double kstar, TH2* matrix) {
  // Define new Histogram which have dimension according to the yaxis (new
  // momentum axis):
  const int nbins_original = matrix->GetYaxis()->GetNbins();

  int binori = 0;
  for (int momOri = 0; momOri < nbins_original; momOri++) {
    Double_t momentum_original = matrix->GetYaxis()->GetBinCenter(momOri + 1);
    if (momentum_original >= kstar) {
      binori = (momOri + 1);
      break;
    }
  }
  return binori;
}

double evalMomTrafo(double bin, TH2* matrix) {
  // Define new Histogram which have dimension according to the yaxis (new
  // momentum axis):
  const Int_t nbins_transformed = matrix->GetXaxis()->GetNbins();
  double max = 0;
  double maxbin = 0;
  for (int momTrans = 0; momTrans < nbins_transformed; momTrans++) {
    double matrixvalues = 0.;
    matrixvalues = matrix->GetBinContent(momTrans + 1, bin);
    // cout<<"matrixval"<<momTrans + 1<< ": "<<matrixvalues<<endl;
    if (matrixvalues >= max) {
      max = matrixvalues;
      maxbin = momTrans + 1;
    }
    // cout<<max<<endl;
  }
  Double_t momentum_transformed = matrix->GetXaxis()->GetBinCenter(maxbin + 1);

  return momentum_transformed;
}

double GetWidth(double bin, TH2* matrix) {
  const Int_t nbins_transformed = matrix->GetXaxis()->GetNbins();
  double bin_width = matrix->GetXaxis()->GetBinWidth(1);
  int count = 0;
  for (int momTrans = 0; momTrans < nbins_transformed; momTrans++) {
    Double_t matrixvalues = 0.;
    matrixvalues = matrix->GetBinContent(momTrans + 1, bin);
    if (matrixvalues != 0) {
      count++;
    }
  }
  double width = count * bin_width;
  return width;
}



double GetProjWidth(TH1* proj){
  const Int_t nbins_transformed = proj->GetNbinsX();
  double bin_width=proj->GetXaxis()->GetBinWidth(1);
  int count = 0;
  for (int momTrans = 0; momTrans < nbins_transformed; momTrans++) {
    Double_t matrixvalues = 0.;
    matrixvalues = proj->GetBinContent(momTrans + 1);
    if (matrixvalues!=0){
        count++;
    }
    }
  double width=count*bin_width;
  return width;
}


double GetExright(TH1* proj, double mean){
  const Int_t nbins_transformed = proj->GetNbinsX();
  double bin_width=proj->GetXaxis()->GetBinWidth(1);
  int count = 0;
  int maxbin=0;
  for (int momTrans = 0; momTrans < nbins_transformed; momTrans++) {
    Double_t matrixvalues = 0.;
    matrixvalues = proj->GetBinContent(momTrans + 1);
    if (matrixvalues!=0){
        count++;
        maxbin=momTrans;

    }
    }
  double width=count*bin_width;
  double max= proj->GetBinCenter(maxbin)+1;
 // cout<<"wifht: "<< width<<" maxbin: "<<maxbin<<"maxvaue: "<<max<<"mean: "<< mean<<endl;

  double er= max-mean;



  return er;

}



double GetExleft(TH1* proj, double mean){
  const Int_t nbins_transformed = proj->GetNbinsX();
  double bin_width=proj->GetXaxis()->GetBinWidth(1);
  int count = 0;
  int maxbin=0;
  for (int momTrans = 0; momTrans < nbins_transformed; momTrans++) {
    Double_t matrixvalues = 0.;
    matrixvalues = proj->GetBinContent(momTrans + 1);
    if (matrixvalues!=0){
        count++;
        maxbin=momTrans;
    }
    }
  double width=count*bin_width;

  double max= proj->GetBinCenter(maxbin)+1;
  //cout<<"wifht: "<< width<<" maxbin: "<<maxbin<<"maxvaue: "<<max<<"mean: "<< mean<<endl;

  double el=mean-(max-width);



  return el;

}


float getchisquared(TH1F *hist, TGraph *func)
{
    int nbins=hist->GetNbinsX();
    double chi2=0;
    double tmp=0;
    cout<< "nbins="<< nbins<<endl;
    double momentum= 0;
    for (int i = 0; i <nbins; ++i ) {
    momentum= hist->GetBinCenter(i + 1);
    if (momentum<=1000){


    double funval= func->Eval(momentum);
    double histval= hist->GetBinContent(i + 1);
    double error=hist->GetBinError(i + 1);
    tmp=(histval-funval)/error;
    chi2+=tmp*tmp;

    //cout<< "momentum: "<< momentum<< "   funval: "<<funval<<"   histval: "<< histval<< "   error: "<< error<<endl;

    }
    }
    cout<< "chi2: ..."<< chi2<< endl;
    return chi2;
}


// void SetCFPoint(Int_t i, TH1F* CF, TH2F* matrix, int binnr, TGraphErrors
// Graph ){
//    double kstar= CF->GetBinCenter(binnr);
//    double ey= CF->GetBinError(binnr);
//    double val= CF->GetBinContent(binnr);

//    int biny=evalMomOri(kstar,matrix);
//    double max=evalMomTrafo(biny,matrix);
//    double ex=GetWidth(biny,matrix);

//    Graph.SetPoint(i,max,val);
//    Graph.SetPointError(i, ex, ey);

//}

/// =====================================================================================
int main(int argc, char* argv[]) {
  gROOT->ProcessLine("gErrorIgnoreLevel = 3001");
  DreamPlot::SetStyle();

  TString InputDir = argv[1];
  const char* trigger = argv[2];
  int suffix = atoi(argv[3]);

  auto filename = TString::Format("%s/SidebandMix.root", InputDir.Data());
  auto outfile = new TFile(filename, "RECREATE");

  auto filenameSmear =
      TString::Format("%s/SmearSideband.root", InputDir.Data());
  auto infile = TFile::Open(filenameSmear);
  if (!infile) {
    std::cout
        << "No smearing matrix found - start the sideband computation task!\n";
    return 0;
  }

  // Get Matrices

  std::cout << "get matrices" << std::endl;
  auto histSmearTOTSBleft = (TH2F*)infile->Get("histSmearTOTSBleft");
  auto histSmearTOTSBright = (TH2F*)infile->Get("histSmearTOTSBright");
  auto histSmearTOTSB12 = (TH2F*)infile->Get("histSmearTOTSB12");
  auto histSmearTOTSB23 = (TH2F*)infile->Get("histSmearTOTSB23");
  auto histSmearTOTSB34 = (TH2F*)infile->Get("histSmearTOTSB34");
  auto histSmearTOTSB45 = (TH2F*)infile->Get("histSmearTOTSB45");
  auto histSmearTOTSB56 = (TH2F*)infile->Get("histSmearTOTSB56");
  auto histSmearTOTSB67 = (TH2F*)infile->Get("histSmearTOTSB67");
  auto histSmearTOTSB78 = (TH2F*)infile->Get("histSmearTOTSB78");
  auto histSmearTOTSB89 = (TH2F*)infile->Get("histSmearTOTSB89");
  auto histSmearTOTSB90 = (TH2F*)infile->Get("histSmearTOTSB90");
  auto histSmearTOTPeak = (TH2F*)infile->Get("histSmearTOTPeak");

  auto histSmearTOTSB2762 = (TH2F*)infile->Get("histSmearTOT2762");
  auto histSmearTOTSB6210 = (TH2F*)infile->Get("histSmearTOT6210");
  auto histSmearTOTSB1015 = (TH2F*)infile->Get("histSmearTOT1015");
  auto histSmearTOTSB1520 = (TH2F*)infile->Get("histSmearTOT1520");
  auto histSmearTOTSB2025 = (TH2F*)infile->Get("histSmearTOT2025");
  auto histSmearTOTSB2530 = (TH2F*)infile->Get("histSmearTOT2530");
  auto histSmearTOTSB3035 = (TH2F*)infile->Get("histSmearTOT3035");
  auto histSmearTOTSB3540 = (TH2F*)infile->Get("histSmearTOT3540");
  auto histSmearTOTSB4045 = (TH2F*)infile->Get("histSmearTOT4045");
  auto histSmearTOTSB4550 = (TH2F*)infile->Get("histSmearTOT4550");
  auto histSmearTOTSB5055 = (TH2F*)infile->Get("histSmearTOT5055");
  auto histSmearTOTSB5560 = (TH2F*)infile->Get("histSmearTOT5560");

  auto histSmearTOTSB2223 = (TH2F*)infile->Get("histSmearTOT2223");
  auto histSmearTOTSB2425 = (TH2F*)infile->Get("histSmearTOT2425");
  auto histSmearTOTSB2627 = (TH2F*)infile->Get("histSmearTOT2627");
  auto histSmearTOTSB2829 = (TH2F*)infile->Get("histSmearTOT2829");
  auto histSmearTOTSB3031 = (TH2F*)infile->Get("histSmearTOT3031");
  auto histSmearTOTSB3334 = (TH2F*)infile->Get("histSmearTOT3334");
  auto histSmearTOTSB3738 = (TH2F*)infile->Get("histSmearTOT3738");
  auto histSmearTOTSB4041 = (TH2F*)infile->Get("histSmearTOT4041");

  const float right = 0.04;
  const float top = 0.025;

  auto dummyHist = new TH1F(
      "dummyHist", ";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)", 100, 4, 1000);
  DreamPlot::SetStyleHisto(dummyHist, 20, kGreen + 2);

  std::cout << "setstyle" << std::endl;
  DreamPlot::SetStyleHisto(histSmearTOTSBleft);
  DreamPlot::SetStyleHisto(histSmearTOTSBright);
  DreamPlot::SetStyleHisto(histSmearTOTSB12);
  DreamPlot::SetStyleHisto(histSmearTOTSB23);
  DreamPlot::SetStyleHisto(histSmearTOTSB34);
  DreamPlot::SetStyleHisto(histSmearTOTSB45);
  DreamPlot::SetStyleHisto(histSmearTOTSB56);
  DreamPlot::SetStyleHisto(histSmearTOTSB67);
  DreamPlot::SetStyleHisto(histSmearTOTSB78);
  DreamPlot::SetStyleHisto(histSmearTOTSB89);
  DreamPlot::SetStyleHisto(histSmearTOTSB90);
  DreamPlot::SetStyleHisto(histSmearTOTPeak);

  //  DreamPlot::SetStyleHisto(histSmearTOTSB2762);
  //  DreamPlot::SetStyleHisto(histSmearTOTSB6210);
  //  DreamPlot::SetStyleHisto(histSmearTOTSB1015);
  //  DreamPlot::SetStyleHisto(histSmearTOTSB1520);
  //  DreamPlot::SetStyleHisto(histSmearTOTSB2025);
  //  DreamPlot::SetStyleHisto(histSmearTOTSB2530);
  //  DreamPlot::SetStyleHisto(histSmearTOTSB3035);
  //  DreamPlot::SetStyleHisto(histSmearTOTSB3540);
  //  DreamPlot::SetStyleHisto(histSmearTOTSB4045);
  //  DreamPlot::SetStyleHisto(histSmearTOTSB4550);
  //  DreamPlot::SetStyleHisto(histSmearTOTSB5055);
  //  DreamPlot::SetStyleHisto(histSmearTOTSB5560);

  //  DreamPlot::SetStyleHisto(histSmearTOTSB2223);
  //  DreamPlot::SetStyleHisto(histSmearTOTSB2425);
  //  DreamPlot::SetStyleHisto(histSmearTOTSB2627);
  //  DreamPlot::SetStyleHisto(histSmearTOTSB2829);
  //  DreamPlot::SetStyleHisto(histSmearTOTSB3031);
  //  DreamPlot::SetStyleHisto(histSmearTOTSB3334);
  //  DreamPlot::SetStyleHisto(histSmearTOTSB3738);
  //  DreamPlot::SetStyleHisto(histSmearTOTSB4041);

  std::cout << "end style" << std::endl;

  auto m1 = new TCanvas("m1", "m1", 650, 550);
  histSmearTOTPeak->Draw("col");
  histSmearTOTPeak->GetXaxis()->SetRangeUser(0, 500);
  histSmearTOTPeak->GetYaxis()->SetRangeUser(0, 500);
  histSmearTOTPeak->GetXaxis()->SetNdivisions(505);
  histSmearTOTPeak->GetYaxis()->SetNdivisions(505);
  histSmearTOTPeak->SetTitle(
      "; #it{k}*_{p#minus K^{+}} (MeV/#it{c}); #it{k}*_{p#minus#phi} "
      "(MeV/#it{c})");
  m1->Print("histSmearTOTPeak.pdf");
  m1->Print("histSmearTOTPeak.png");
  delete m1;

  auto m2 = new TCanvas("m2", "m2", 650, 550);
  histSmearTOTSBleft->Draw("col");
  histSmearTOTSBleft->GetXaxis()->SetRangeUser(0, 500);
  histSmearTOTSBleft->GetYaxis()->SetRangeUser(0, 500);
  histSmearTOTSBleft->GetXaxis()->SetNdivisions(505);
  histSmearTOTSBleft->GetYaxis()->SetNdivisions(505);
  histSmearTOTSBleft->SetTitle(
      "; #it{k}*_{p#minus K^{+}} (MeV/#it{c}); #it{k}*_{p#minus (K^{+}K^{-})} "
      "(MeV/#it{c})");
  m2->Print("histSmearTOTSBleft.pdf");
  m2->Print("histSmearTOTSBleft.png");
  delete m2;

  TH1F* SBleftc = CalcCF("HMPhi", 1);
  TH1F* SBrightc = CalcCF("HMPhi", 2);
  TH1F* SB12c = CalcCF("HMPhi", 3);
  TH1F* SB23c = CalcCF("HMPhi", 4);
  TH1F* SB34c = CalcCF("HMPhi", 5);
  TH1F* SB45c = CalcCF("HMPhi", 6);
  TH1F* SB56c = CalcCF("HMPhi", 7);
  TH1F* SB67c = CalcCF("HMPhi", 8);
  TH1F* SB78c = CalcCF("HMPhi", 9);
  TH1F* SB89c = CalcCF("HMPhi", 10);
  TH1F* SB90c = CalcCF("HMPhi", 11);
  TH1F* pPhiCFc = CalcCF("HMPhi", 0);

  pPhiCFc->SetMarkerStyle(20);
  pPhiCFc->SetMarkerSize(1);
  pPhiCFc->SetMarkerColor(kBlack);
  pPhiCFc->SetLineColor(kBlack);
  pPhiCFc->SetLineWidth(2);

  SBleftc->SetMarkerStyle(20);
  SBleftc->SetMarkerSize(1);
  SBleftc->SetMarkerColor(kBlack);
  SBleftc->SetLineColor(kBlack);
  SBleftc->SetLineWidth(2);

  SBrightc->SetMarkerStyle(20);
  SBrightc->SetMarkerSize(1);
  SBrightc->SetMarkerColor(kBlack);
  SBrightc->SetLineColor(kBlack);
  SBrightc->SetLineWidth(2);

  SB12c->SetMarkerStyle(20);
  SB12c->SetMarkerSize(1);
  SB12c->SetMarkerColor(kBlack);
  SB12c->SetLineColor(kBlack);
  SB12c->SetLineWidth(2);

  SB23c->SetMarkerStyle(20);
  SB23c->SetMarkerSize(1);
  SB23c->SetMarkerColor(kBlack);
  SB23c->SetLineColor(kBlack);
  SB23c->SetLineWidth(2);

  SB34c->SetMarkerStyle(20);
  SB34c->SetMarkerSize(1);
  SB34c->SetMarkerColor(kBlack);
  SB34c->SetLineColor(kBlack);
  SB34c->SetLineWidth(2);

  SB45c->SetMarkerStyle(20);
  SB45c->SetMarkerSize(1);
  SB45c->SetMarkerColor(kBlack);
  SB45c->SetLineColor(kBlack);
  SB45c->SetLineWidth(2);

  SB56c->SetMarkerStyle(20);
  SB56c->SetMarkerSize(1);
  SB56c->SetMarkerColor(kBlack);
  SB56c->SetLineColor(kBlack);
  SB56c->SetLineWidth(2);

  SB67c->SetMarkerStyle(20);
  SB67c->SetMarkerSize(1);
  SB67c->SetMarkerColor(kBlack);
  SB67c->SetLineColor(kBlack);
  SB67c->SetLineWidth(2);

  SB78c->SetMarkerStyle(20);
  SB78c->SetMarkerSize(1);
  SB78c->SetMarkerColor(kBlack);
  SB78c->SetLineColor(kBlack);
  SB78c->SetLineWidth(2);

  SB89c->SetMarkerStyle(20);
  SB89c->SetMarkerSize(1);
  SB89c->SetMarkerColor(kBlack);
  SB89c->SetLineColor(kBlack);
  SB89c->SetLineWidth(2);

  SB90c->SetMarkerStyle(20);
  SB90c->SetMarkerSize(1);
  SB90c->SetMarkerColor(kBlack);
  SB90c->SetLineColor(kBlack);
  SB90c->SetLineWidth(2);



  TH1F* pPhiCF = CalcCF2("HMPhi", 0);
  TH1F* SBleft = CalcCF2("HMPhi", 1);
  TH1F* SB2762 = CalcCF2("HMPhi", 2);
  TH1F* SB6210 = CalcCF2("HMPhi", 3);
  TH1F* SB1015 = CalcCF2("HMPhi", 4);
  TH1F* SB1520 = CalcCF2("HMPhi", 5);
  TH1F* SB2025 = CalcCF2("HMPhi", 6);
  TH1F* SB2530 = CalcCF2("HMPhi", 7);
  TH1F* SB3035 = CalcCF2("HMPhi", 8);
  TH1F* SB3540 = CalcCF2("HMPhi", 9);
  TH1F* SB4045 = CalcCF2("HMPhi", 10);
  TH1F* SB4550 = CalcCF2("HMPhi", 11);
  TH1F* SB5055 = CalcCF2("HMPhi", 12);
  TH1F* SB5560 = CalcCF2("HMPhi", 13);
  TH1F* SB67 = CalcCF2("HMPhi", 14);
  TH1F* SB78 = CalcCF2("HMPhi", 15);
  TH1F* SB89 = CalcCF2("HMPhi", 16);
  TH1F* SB90 = CalcCF2("HMPhi", 17);
  TH1F* SB2223 = CalcCF2("HMPhi", 18);
  TH1F* SB2425 = CalcCF2("HMPhi", 19);
  TH1F* SB2627 = CalcCF2("HMPhi", 20);
  TH1F* SB2829 = CalcCF2("HMPhi", 21);
  TH1F* SB3031 = CalcCF2("HMPhi", 22);
  TH1F* SB3334 = CalcCF2("HMPhi", 23);
  TH1F* SB3738 = CalcCF2("HMPhi", 24);
  TH1F* SB4041 = CalcCF2("HMPhi", 25);

  pPhiCF->SetMarkerStyle(20);
  pPhiCF->SetMarkerSize(1);
  pPhiCF->SetMarkerColor(kBlack);
  pPhiCF->SetLineColor(kBlack);
  pPhiCF->SetLineWidth(2);

  SBleft->SetMarkerStyle(20);
  SBleft->SetMarkerSize(1);
  SBleft->SetMarkerColor(kBlack);
  SBleft->SetLineColor(kBlack);
  SBleft->SetLineWidth(2);

  SB2762->SetMarkerStyle(20);
  SB2762->SetMarkerSize(1);
  SB2762->SetMarkerColor(kBlack);
  SB2762->SetLineColor(kBlack);
  SB2762->SetLineWidth(2);

  SB6210->SetMarkerStyle(20);
  SB6210->SetMarkerSize(1);
  SB6210->SetMarkerColor(kBlack);
  SB6210->SetLineColor(kBlack);
  SB6210->SetLineWidth(2);

  SB1015->SetMarkerStyle(20);
  SB1015->SetMarkerSize(1);
  SB1015->SetMarkerColor(kBlack);
  SB1015->SetLineColor(kBlack);
  SB1015->SetLineWidth(2);

  SB1520->SetMarkerStyle(20);
  SB1520->SetMarkerSize(1);
  SB1520->SetMarkerColor(kBlack);
  SB1520->SetLineColor(kBlack);
  SB1520->SetLineWidth(2);

  SB2025->SetMarkerStyle(20);
  SB2025->SetMarkerSize(1);
  SB2025->SetMarkerColor(kBlack);
  SB2025->SetLineColor(kBlack);
  SB2025->SetLineWidth(2);

  SB2530->SetMarkerStyle(20);
  SB2530->SetMarkerSize(1);
  SB2530->SetMarkerColor(kBlack);
  SB2530->SetLineColor(kBlack);
  SB2530->SetLineWidth(2);

  SB3035->SetMarkerStyle(20);
  SB3035->SetMarkerSize(1);
  SB3035->SetMarkerColor(kBlack);
  SB3035->SetLineColor(kBlack);
  SB3035->SetLineWidth(2);

  SB3540->SetMarkerStyle(20);
  SB3540->SetMarkerSize(1);
  SB3540->SetMarkerColor(kBlack);
  SB3540->SetLineColor(kBlack);
  SB3540->SetLineWidth(2);

  SB4045->SetMarkerStyle(20);
  SB4045->SetMarkerSize(1);
  SB4045->SetMarkerColor(kBlack);
  SB4045->SetLineColor(kBlack);
  SB4045->SetLineWidth(2);

  SB4550->SetMarkerStyle(20);
  SB4550->SetMarkerSize(1);
  SB4550->SetMarkerColor(kBlack);
  SB4550->SetLineColor(kBlack);
  SB4550->SetLineWidth(2);

  SB5055->SetMarkerStyle(20);
  SB5055->SetMarkerSize(1);
  SB5055->SetMarkerColor(kBlack);
  SB5055->SetLineColor(kBlack);
  SB5055->SetLineWidth(2);

  SB5560->SetMarkerStyle(20);
  SB5560->SetMarkerSize(1);
  SB5560->SetMarkerColor(kBlack);
  SB5560->SetLineColor(kBlack);
  SB5560->SetLineWidth(2);

  SB67->SetMarkerStyle(20);
  SB67->SetMarkerSize(1);
  SB67->SetMarkerColor(kBlack);
  SB67->SetLineColor(kBlack);
  SB67->SetLineWidth(2);

  SB78->SetMarkerStyle(20);
  SB78->SetMarkerSize(1);
  SB78->SetMarkerColor(kBlack);
  SB78->SetLineColor(kBlack);
  SB78->SetLineWidth(2);

  SB89->SetMarkerStyle(20);
  SB89->SetMarkerSize(1);
  SB89->SetMarkerColor(kBlack);
  SB89->SetLineColor(kBlack);
  SB89->SetLineWidth(2);

  SB90->SetMarkerStyle(20);
  SB90->SetMarkerSize(1);
  SB90->SetMarkerColor(kBlack);
  SB90->SetLineColor(kBlack);
  SB90->SetLineWidth(2);

  SB2223->SetMarkerStyle(20);
  SB2223->SetMarkerSize(1);
  SB2223->SetMarkerColor(kBlack);
  SB2223->SetLineColor(kBlack);
  SB2223->SetLineWidth(2);

  SB2425->SetMarkerStyle(20);
  SB2425->SetMarkerSize(1);
  SB2425->SetMarkerColor(kBlack);
  SB2425->SetLineColor(kBlack);
  SB2425->SetLineWidth(2);

  SB2627->SetMarkerStyle(20);
  SB2627->SetMarkerSize(1);
  SB2627->SetMarkerColor(kBlack);
  SB2627->SetLineColor(kBlack);
  SB2627->SetLineWidth(2);

  SB2829->SetMarkerStyle(20);
  SB2829->SetMarkerSize(1);
  SB2829->SetMarkerColor(kBlack);
  SB2829->SetLineColor(kBlack);
  SB2829->SetLineWidth(2);

  SB3031->SetMarkerStyle(20);
  SB3031->SetMarkerSize(1);
  SB3031->SetMarkerColor(kBlack);
  SB3031->SetLineColor(kBlack);
  SB3031->SetLineWidth(2);

  SB3334->SetMarkerStyle(20);
  SB3334->SetMarkerSize(1);
  SB3334->SetMarkerColor(kBlack);
  SB3334->SetLineColor(kBlack);
  SB3334->SetLineWidth(2);

  SB3738->SetMarkerStyle(20);
  SB3738->SetMarkerSize(1);
  SB3738->SetMarkerColor(kBlack);
  SB3738->SetLineColor(kBlack);
  SB3738->SetLineWidth(2);

  SB4041->SetMarkerStyle(20);
  SB4041->SetMarkerSize(1);
  SB4041->SetMarkerColor(kBlack);
  SB4041->SetLineColor(kBlack);
  SB4041->SetLineWidth(2);

  //  void SetCFPoint(Int_t i, TH1F* CF, TH2F* matrix, int binnr, TGraphErrors
  //  Graph ){
  //      double kstar= CF->GetBinCenter(binnr);
  //      double ey= CF->GetBinError(binnr);
  //      double val= CF->GetBinContent(binnr);

  //      int biny=evalMomOri(kstar,matrix);
  //      double max=evalMomTrafo(biny,matrix);
  //      double ex=GetWidth(biny,matrix);

  //      Graph.SetPoint(i,max,val);
  //      Graph.SetPointError(i, ex, ey);

  //  }

  // auto pKCF = new TGraphErrors();

  // double k1= SBleft->GetBinCenter(1);
  // double ey1= SBleft->GetBinError(1);
  // double C1= SBleft->GetBinContent(1);

  // int biny1=evalMomOri(k1,histSmearTOTSBleft);
  // double max1=evalMomTrafo(biny1,histSmearTOTSBleft);
  // double ex1=GetWidth(biny1,histSmearTOTSBleft);
  // pKCF.SetPoint(1,max1,C1);
  // pKCF.SetPointError(1, ex1, ey1);

  const Int_t n = 25;
  Double_t x[n] = {0};
  Double_t y[n] = {0};
  Double_t ex[n] = {0};
  Double_t ey[n] = {0};
  Double_t exl[n] = {0};
  Double_t exh[n] = {0};


  const Int_t n1 = 9;
  Double_t xx[n1] = {0};
  Double_t yy[n1] = {0};
  Double_t exx[n1] = {0};
  Double_t eyy[n1] = {0};

//  int i = 2;  // cf bin nr

//  y[0] = SBleft->GetBinContent(i);
//  cout << y[0] << endl;
//  ey[0] = SBleft->GetBinError(i);
//  double k0 = SBleft->GetBinCenter(i);
//  int biny0 = evalMomOri(k0, histSmearTOTSBleft);
//  double max0 = evalMomTrafo(biny0, histSmearTOTSBleft);
//  double ex0 = GetWidth(biny0, histSmearTOTSBleft) / 2;
//  x[0] = max0;
//  ex[0] = ex0;

//  yy[0] = y[0];
//  xx[0] = x[0];
//  eyy[0] = ey[0];
//  exx[0] = ex[0];

//  y[1] = SB2762->GetBinContent(i);
//  ey[1] = SB2762->GetBinError(i);
//  double k1 = SB2762->GetBinCenter(i);
//  int biny1 = evalMomOri(k1, histSmearTOTSB2762);
//  double max1 = evalMomTrafo(biny1, histSmearTOTSB2762);
//  double ex1 = GetWidth(biny1, histSmearTOTSB2762) / 2;
//  x[1] = max1;
//  ex[1] = ex1;

//  y[2] = SB6210->GetBinContent(i);
//  ey[2] = SB6210->GetBinError(i);
//  double k2 = SB6210->GetBinCenter(i);
//  int biny2 = evalMomOri(k2, histSmearTOTSB6210);
//  double max2 = evalMomTrafo(biny2, histSmearTOTSB6210);
//  double ex2 = GetWidth(biny2, histSmearTOTSB6210) / 2;
//  x[2] = max2;
//  ex[2] = ex2;

//  y[3] = SB1015->GetBinContent(i);
//  ey[3] = SB1015->GetBinError(i);
//  double k3 = SB1015->GetBinCenter(i);
//  int biny3 = evalMomOri(k2, histSmearTOTSB1015);
//  double max3 = evalMomTrafo(biny2, histSmearTOTSB1015);
//  double ex3 = GetWidth(biny2, histSmearTOTSB1015) / 2;
//  x[3] = max3;
//  ex[3] = ex3;

//  y[4] = SB1520->GetBinContent(i);
//  ey[4] = SB1520->GetBinError(i);
//  double k4 = SB1520->GetBinCenter(i);
//  int biny4 = evalMomOri(k4, histSmearTOTSB1520);
//  double max4 = evalMomTrafo(biny4, histSmearTOTSB1520);
//  double ex4 = GetWidth(biny4, histSmearTOTSB1520) / 2;
//  x[4] = max4;
//  ex[4] = ex4;

//  y[5] = SB2025->GetBinContent(i);
//  ey[5] = SB2025->GetBinError(i);
//  double k5 = SB2025->GetBinCenter(i);
//  int biny5 = evalMomOri(k5, histSmearTOTSB2025);
//  double max5 = evalMomTrafo(biny5, histSmearTOTSB2025);
//  double ex5 = GetWidth(biny5, histSmearTOTSB2025) / 2;
//  x[5] = max5;
//  ex[5] = ex5;

//  y[6] = SB2530->GetBinContent(i);
//  ey[6] = SB2530->GetBinError(i);
//  double k6 = SB2530->GetBinCenter(i);
//  int biny6 = evalMomOri(k6, histSmearTOTSB2530);
//  double max6 = evalMomTrafo(biny6, histSmearTOTSB2530);
//  double ex6 = GetWidth(biny6, histSmearTOTSB2530) / 2;
//  x[6] = max6;
//  ex[6] = ex6;

//  y[7] = SB3035->GetBinContent(i);
//  ey[7] = SB3035->GetBinError(i);
//  double k7 = SB3035->GetBinCenter(i);
//  int biny7 = evalMomOri(k7, histSmearTOTSB3035);
//  double max7 = evalMomTrafo(biny7, histSmearTOTSB3035);
//  double ex7 = GetWidth(biny7, histSmearTOTSB3035) / 2;
//  x[7] = max7;
//  ex[7] = ex7;

//  y[8] = SB3540->GetBinContent(i);
//  ey[8] = SB3540->GetBinError(i);
//  double k8 = SB3540->GetBinCenter(i);
//  int biny8 = evalMomOri(k8, histSmearTOTSB3540);
//  double max8 = evalMomTrafo(biny8, histSmearTOTSB3540);
//  double ex8 = GetWidth(biny8, histSmearTOTSB3540) / 2;
//  x[8] = max8;
//  ex[8] = ex8;

//  y[9] = SB4045->GetBinContent(i);
//  // y[9]=1.01;
//  ey[9] = SB4045->GetBinError(i);
//  double k9 = SB4045->GetBinCenter(i);
//  int biny9 = evalMomOri(k9, histSmearTOTSB4045);
//  double max9 = evalMomTrafo(biny9, histSmearTOTSB4045);
//  double ex9 = GetWidth(biny9, histSmearTOTSB4045) / 2;
//  x[9] = max9;
//  ex[9] = ex9;

//  y[10] = SB4550->GetBinContent(i);
//  ey[10] = SB4550->GetBinError(i);
//  double k10 = SB4550->GetBinCenter(i);
//  int biny10 = evalMomOri(k10, histSmearTOTSB4550);
//  double max10 = evalMomTrafo(biny10, histSmearTOTSB4550);
//  double ex10 = GetWidth(biny10, histSmearTOTSB4550) / 2;
//  x[10] = max10;
//  ex[10] = ex10;

//  y[11] = SB5055->GetBinContent(i);
//  ey[11] = SB5055->GetBinError(i);
//  double k11 = SB5055->GetBinCenter(i);
//  int biny11 = evalMomOri(k11, histSmearTOTSB5055);
//  double max11 = evalMomTrafo(biny11, histSmearTOTSB5055);
//  double ex11 = GetWidth(biny11, histSmearTOTSB5055) / 2;
//  x[11] = max11;
//  ex[11] = ex11;

//  y[12] = SB5560->GetBinContent(i);
//  ey[12] = SB5560->GetBinError(i);
//  double k12 = SB5560->GetBinCenter(i);
//  int biny12 = evalMomOri(k12, histSmearTOTSB5560);
//  double max12 = evalMomTrafo(biny12, histSmearTOTSB5560);
//  double ex12 = GetWidth(biny12, histSmearTOTSB5560) / 2;
//  x[12] = max12;
//  ex[12] = ex12;

//  y[13] = SB67->GetBinContent(i);
//  ey[13] = SB67->GetBinError(i);
//  double k13 = SB67->GetBinCenter(i);
//  int biny13 = evalMomOri(k13, histSmearTOTSB67);
//  double max13 = evalMomTrafo(biny13, histSmearTOTSB67);
//  double ex13 = GetWidth(biny13, histSmearTOTSB67) / 2;
//  x[13] = max13;
//  ex[13] = ex13;

//  y[14] = SB78->GetBinContent(i);
//  ey[14] = SB78->GetBinError(i);
//  double k14 = SB78->GetBinCenter(i);
//  int biny14 = evalMomOri(k14, histSmearTOTSB78);
//  double max14 = evalMomTrafo(biny14, histSmearTOTSB78);
//  double ex14 = GetWidth(biny14, histSmearTOTSB78) / 2;
//  x[14] = max14;
//  ex[14] = ex14;

//  y[15] = SB89->GetBinContent(i);
//  ey[15] = SB89->GetBinError(i);
//  double k15 = SB89->GetBinCenter(i);
//  int biny15 = evalMomOri(k15, histSmearTOTSB89);
//  double max15 = evalMomTrafo(biny15, histSmearTOTSB89);
//  double ex15 = GetWidth(biny15, histSmearTOTSB89) / 2;
//  x[15] = max15;
//  ex[15] = ex15;

//  y[16] = SB90->GetBinContent(i);
//  ey[16] = SB90->GetBinError(i);
//  double k16 = SB90->GetBinCenter(i);
//  int biny16 = evalMomOri(k16, histSmearTOTSB90);
//  double max16 = evalMomTrafo(biny16, histSmearTOTSB90);
//  double ex16 = GetWidth(biny16, histSmearTOTSB90) / 2;
//  x[16] = max16;
//  ex[16] = ex16;

//  y[17] = SB2223->GetBinContent(i);
//  ey[17] = SB2223->GetBinError(i);
//  double k17 = SB2223->GetBinCenter(i);
//  int biny17 = evalMomOri(k17, histSmearTOTSB2223);
//  double max17 = evalMomTrafo(biny17, histSmearTOTSB2223);
//  double ex17 = GetWidth(biny17, histSmearTOTSB2223) / 2;
//  x[17] = max17;
//  ex[17] = ex17;

//  yy[1] = y[17];
//  xx[1] = x[17];
//  eyy[1] = ey[17];
//  exx[1] = ex[17];

//  y[18] = SB2425->GetBinContent(i);
//  ey[18] = SB2425->GetBinError(i);
//  double k18 = SB2425->GetBinCenter(i);
//  int biny18 = evalMomOri(k18, histSmearTOTSB2425);
//  double max18 = evalMomTrafo(biny18, histSmearTOTSB2425);
//  double ex18 = GetWidth(biny18, histSmearTOTSB2425) / 2;
//  x[18] = max18;
//  ex[18] = ex18;

//  yy[2] = y[18];
//  xx[2] = x[18];
//  eyy[2] = ey[18];
//  exx[2] = ex[18];

//  y[19] = SB2627->GetBinContent(i);
//  ey[19] = SB2627->GetBinError(i);
//  double k19 = SB2627->GetBinCenter(i);
//  int biny19 = evalMomOri(k19, histSmearTOTSB2627);
//  double max19 = evalMomTrafo(biny19, histSmearTOTSB2627);
//  double ex19 = GetWidth(biny19, histSmearTOTSB2627) / 2;
//  x[19] = max19;
//  ex[19] = ex19;

//  yy[3] = y[19];
//  xx[3] = x[19];
//  eyy[3] = ey[19];
//  exx[3] = ex[19];

//  y[20] = SB2829->GetBinContent(i);
//  ey[20] = SB2829->GetBinError(i);
//  double k20 = SB2829->GetBinCenter(i);
//  int biny20 = evalMomOri(k20, histSmearTOTSB2829);
//  double max20 = evalMomTrafo(biny20, histSmearTOTSB2829);
//  double ex20 = GetWidth(biny20, histSmearTOTSB2829) / 2;
//  x[20] = max20;
//  ex[20] = ex20;

//  yy[4] = y[20];
//  xx[4] = x[20];
//  eyy[4] = ey[20];
//  exx[4] = ex[20];

//  y[21] = SB3031->GetBinContent(i);
//  ey[21] = SB3031->GetBinError(i);
//  double k21 = SB3031->GetBinCenter(i);
//  int biny21 = evalMomOri(k21, histSmearTOTSB3031);
//  double max21 = evalMomTrafo(biny21, histSmearTOTSB3031);
//  double ex21 = GetWidth(biny21, histSmearTOTSB3031) / 2;
//  x[21] = max21;
//  ex[21] = ex21;

//  yy[5] = y[21];
//  xx[5] = x[21];
//  eyy[5] = ey[21];
//  exx[5] = ex[21];

//  y[22] = SB3334->GetBinContent(i);
//  ey[22] = SB3334->GetBinError(i);
//  double k22 = SB3334->GetBinCenter(i);
//  int biny22 = evalMomOri(k22, histSmearTOTSB3334);
//  double max22 = evalMomTrafo(biny22, histSmearTOTSB3334);
//  double ex22 = GetWidth(biny22, histSmearTOTSB3334) / 2;
//  x[22] = max22;
//  ex[22] = ex22;

//  yy[6] = y[22];
//  xx[6] = x[22];
//  eyy[6] = ey[22];
//  exx[6] = ex[22];

//  y[23] = SB3738->GetBinContent(i);
//  ey[23] = SB3738->GetBinError(i);
//  double k23 = SB3738->GetBinCenter(i);
//  int biny23 = evalMomOri(k23, histSmearTOTSB3738);
//  double max23 = evalMomTrafo(biny23, histSmearTOTSB3738);
//  double ex23 = GetWidth(biny23, histSmearTOTSB3738) / 2;
//  x[23] = max23;
//  ex[23] = ex23;

//  yy[7] = y[23];
//  xx[7] = x[23];
//  eyy[7] = ey[23];
//  exx[7] = ex[23];

//  y[24] = SB4041->GetBinContent(i);
//  ey[24] = SB4041->GetBinError(i);
//  double k24 = SB4041->GetBinCenter(i);
//  int biny24 = evalMomOri(k24, histSmearTOTSB4041);
//  double max24 = evalMomTrafo(biny24, histSmearTOTSB4041);
//  double ex24 = GetWidth(biny24, histSmearTOTSB4041) / 2;
//  x[24] = max24;
//  ex[24] = ex24;

//  yy[8] = y[24];
//  xx[8] = x[24];
//  eyy[8] = ey[24];
//  exx[8] = ex[24];

//  // y[666]=SB12->GetBinContent(i);
//  // ey[666]=SB12->GetBinError(i);
//  // double k666= SB12->GetBinCenter(i);
//  // int biny666=evalMomOri(k666,histSmearTOTSB12);
//  // double max666=evalMomTrafo(biny666,histSmearTOTSB12);
//  // double ex666=GetWidth(biny666,histSmearTOTSB12)/2;
//  // x[666]=max666;
//  // ex[666]=ex666;

//  auto pKCF = new TGraphErrors(n, x, y, ex, ey);
//  int lineWidth2 = 3;
//  DreamPlot::SetStyleGraph(pKCF, 20, kBlack, 0.7);
//  pKCF->SetLineWidth(lineWidth2);

//  auto pKCF2 = new TGraphErrors(n, x, y, ex, ey);
//  DreamPlot::SetStyleGraph(pKCF2, 20, kBlack, 0.7);
//  pKCF2->SetLineWidth(lineWidth2);

//  auto pKBG = new TGraphErrors(n1, xx, yy, exx, eyy);
//  DreamPlot::SetStyleGraph(pKBG, 20, kBlack, 0.7);
//  pKBG->SetLineWidth(lineWidth2);


  int i=2; //cf bin nr

  y[0]=SBleft->GetBinContent(i);
  ey[0]=SBleft->GetBinError(i);
  double k0= SBleft->GetBinCenter(i);
  int bw0=(SBleft->GetBinWidth(i))/2;
  int binminy0=evalMomOri(k0-bw0,histSmearTOTSBleft);
  int binmaxy0=evalMomOri(k0+bw0,histSmearTOTSBleft);
  auto projection0 = (TH1F*)((TH2F*)histSmearTOTSBleft->ProjectionX("ll", binminy0, binmaxy0-1));
  projection0->SetFillColor(kBlue+1);
  DreamPlot::SetStyleHisto(projection0);
  outfile->cd();
  projection0->Write("proj0");
  //int maxb0=projection0->GetMaximumBin();
  //double max0=projection0->GetBinCenter( maxb0);
  //int biny0=evalMomOri(k0,histSmearTOTSBleft);
  //double max0=evalMomTrafo(biny0,histSmearTOTSBleft);

  int max0=projection0->GetMean();
  double ex0=GetProjWidth(projection0)/2;
  double exvol0=GetProjWidth(projection0);
  double ex0_l=GetExleft(projection0,max0);
  double ex0_r=GetExright(projection0,max0);

  cout<<"halferror: "<<ex0<<"  errorleft: "<<ex0_l<<"  errorright: "<<ex0_r<<endl;


  x[0]=max0;
  ex[0]=ex0;
  exl[0] =ex0_l;
  exh[0] =ex0_r;



  yy[0]=y[0];
  xx[0]=x[0];
  eyy[0]=ey[0];
  exx[0]=ex[0];






  //y[0]=SBleft->GetBinContent(i);
  //cout<<y[0]<<endl;
  //ey[0]=SBleft->GetBinError(i);
  //double k0= SBleft->GetBinCenter(i);
  //int biny0=evalMomOri(k0,histSmearTOTSBleft);
  //double max0=evalMomTrafo(biny0,histSmearTOTSBleft);
  //double ex0=GetWidth(biny0,histSmearTOTSBleft)/2;
  //x[0]=max0;
  //ex[0]=ex0;

  //yy[0]=y[0];
  //xx[0]=x[0];
  //eyy[0]=ey[0];
  //exx[0]=ex[0];

  y[1]=SB2762->GetBinContent(i);
  ey[1]=SB2762->GetBinError(i);
  double k1= SB2762->GetBinCenter(i);
  int bw1=(SBleft->GetBinWidth(i))/2;
  int binminy1=evalMomOri(k1-bw1,histSmearTOTSB2762);
  int binmaxy1=evalMomOri(k1+bw1,histSmearTOTSB2762);
  auto projection1 = (TH1F*)((TH2F*)histSmearTOTSB2762->ProjectionX("ll", binminy1, binmaxy1-1));
  projection1->SetFillColor(kBlue+1);
  DreamPlot::SetStyleHisto(projection1);
  outfile->cd();
  projection1->Write("proj1");

  //int maxb0=projection0->GetMaximumBin();
  //double max0=projection0->GetBinCenter( maxb0);
  //int biny0=evalMomOri(k0,histSmearTOTSBleft);
  //double max0=evalMomTrafo(biny0,histSmearTOTSBleft);

  int max1=projection1->GetMean();
  double ex1=GetProjWidth(projection1)/2;
  double exvol1=GetProjWidth(projection1);
  double ex1_l=GetExleft(projection1,max1);
  double ex1_r=GetExright(projection1,max1);

  cout<<"halferror: "<<ex1<<"  errorleft: "<<ex1_l<<"  errorright: "<<ex1_r<<endl;

  exl[1] =ex1_l;
  exh[1] =ex1_r;

  x[1]=max1;
  ex[1]=ex1;




  y[2]=SB6210->GetBinContent(i);
  ey[2]=SB6210->GetBinError(i);
  double k2= SB6210->GetBinCenter(i);
  int bw2=(SBleft->GetBinWidth(i))/2;
  int binminy2=evalMomOri(k2-bw2,histSmearTOTSB6210);
  int binmaxy2=evalMomOri(k2+bw2,histSmearTOTSB6210);
  auto projection2 = (TH2F*)((TH2F*)histSmearTOTSB6210->ProjectionX("ll", binminy2, binmaxy2-1));
  projection2->SetFillColor(kBlue+2);
  DreamPlot::SetStyleHisto(projection2);
  outfile->cd();
  projection2->Write("proj2");

  //int maxb2=projection2->GetMaximumBin();
  //double max2=projection2->GetBinCenter(maxb2);
  //int biny2=evalMomOri(k2,histSmearTOTSB6210);
  //double max2=evalMomTrafo(biny2,histSmearTOTSB6210);

  int max2=projection2->GetMean();
  double ex2=GetProjWidth(projection2)/2;
  double exvol2=GetProjWidth(projection2);
  double ex2_l=GetExleft(projection2,max2);
  double ex2_r=GetExright(projection2,max2);

  cout<<"halferror: "<<ex2<<"  errorleft: "<<ex2_l<<"  errorright: "<<ex2_r<<endl;

  exl[2] =ex2_l;
  exh[2] =ex2_r;
  //double ex2=GetWidth(biny2,histSmearTOTSB6210)/2;
  x[2]=max2;
  ex[2]=ex2;





  y[3]=SB1015->GetBinContent(i);
  ey[3]=SB1015->GetBinError(i);
  double k3= SB1015->GetBinCenter(i);
  int bw3=(SBleft->GetBinWidth(i))/2;
  int binminy3=evalMomOri(k3-bw3,histSmearTOTSB1015);
  int binmaxy3=evalMomOri(k3+bw3,histSmearTOTSB1015);
  auto projection3 = (TH2F*)((TH2F*)histSmearTOTSB1015->ProjectionX("ll", binminy3, binmaxy3-1));
  projection3->SetFillColor(kBlue+2);
  DreamPlot::SetStyleHisto(projection3);
  outfile->cd();
  projection3->Write("proj3");

  //int maxb3=projection3->GetMaximumBin();
  //double max3=projection3->GetBinCenter(maxb3);
  //int biny3=evalMomOri(k3,histSmearTOTSB1015);
  //double max3=evalMomTrafo(biny2,histSmearTOTSB1015);

  int max3=projection3->GetMean();
  double ex3=GetProjWidth(projection3)/2;
  double exvol3=GetProjWidth(projection3);
  double ex3_l=GetExleft(projection3,max3);
  double ex3_r=GetExright(projection3,max3);

  cout<<"halferror: "<<ex3<<"  errorleft: "<<ex3_l<<"  errorright: "<<ex3_r<<endl;

  exl[3] =ex3_l;
  exh[3] =ex3_r;
  //double ex3=GetWidth(biny2,histSmearTOTSB1015)/2;
  x[3]=max3;
  ex[3]=ex3;



  y[4]=SB1520->GetBinContent(i);
  ey[4]=SB1520->GetBinError(i);
  double k4= SB1520->GetBinCenter(i);
  int bw4=(SBleft->GetBinWidth(i))/2;
  int binminy4=evalMomOri(k4-bw4,histSmearTOTSB1520);
  int binmaxy4=evalMomOri(k4+bw4,histSmearTOTSB1520);
  auto projection4 = (TH2F*)((TH2F*)histSmearTOTSB1520->ProjectionX("ll", binminy4, binmaxy4-1));
  projection4->SetFillColor(kBlue+2);
  DreamPlot::SetStyleHisto(projection4);
  outfile->cd();
  projection4->Write("proj4");

  //int maxb4=projection4->GetMaximumBin();
  //double max4=projection4->GetBinCenter(maxb4);
  //int biny4=evalMomOri(k4,histSmearTOTSB1520);
  //double max4=evalMomTrafo(biny4,histSmearTOTSB1520);

  int max4=projection4->GetMean();
  double ex4=GetProjWidth(projection4)/2;
  double exvol4=GetProjWidth(projection4);
  double ex4_l=GetExleft(projection4,max4);
  double ex4_r=GetExright(projection4,max4);

  cout<<"halferror: "<<ex4<<"  errorleft: "<<ex4_l<<"  errorright: "<<ex4_r<<endl;

  exl[4] =ex4_l;
  exh[4] =ex4_r;
  //double ex4=GetWidth(biny4,histSmearTOTSB1520)/2;
  x[4]=max4;
  ex[4]=ex4;


  y[5]=SB2025->GetBinContent(i);
  ey[5]=SB2025->GetBinError(i);
  double k5= SB2025->GetBinCenter(i);
  int bw5=(SBleft->GetBinWidth(i))/2;
  int binminy5=evalMomOri(k5-bw5,histSmearTOTSB2025);
  int binmaxy5=evalMomOri(k5+bw5,histSmearTOTSB2025);
  auto projection5 = (TH2F*)((TH2F*)histSmearTOTSB2025->ProjectionX("ll", binminy5, binmaxy5-1));
  projection5->SetFillColor(kBlue+2);
  DreamPlot::SetStyleHisto(projection5);
  outfile->cd();
  projection5->Write("proj5");

  //int maxb5=projection5->GetMaximumBin();
  //double max5=projection5->GetBinCenter(maxb5);
  //int biny5=evalMomOri(k5,histSmearTOTSB2025);
  //double max5=evalMomTrafo(biny5,histSmearTOTSB2025);

  int max5=projection5->GetMean();
  double ex5=GetProjWidth(projection5)/2;
  double exvol5=GetProjWidth(projection5);
  double ex5_l=GetExleft(projection5,max5);
  double ex5_r=GetExright(projection5,max5);

  cout<<"halferror: "<<ex5<<"  errorleft: "<<ex5_l<<"  errorright: "<<ex5_r<<endl;

  exl[5] =ex5_l;
  exh[5] =ex5_r;
  //double ex5=GetWidth(biny5,histSmearTOTSB2025)/2;
  x[5]=max5;
  ex[5]=ex5;

  y[6]=SB2530->GetBinContent(i);
  ey[6]=SB2530->GetBinError(i);
  double k6= SB2530->GetBinCenter(i);
  int bw6=(SBleft->GetBinWidth(i))/2;
  int binminy6=evalMomOri(k6-bw6,histSmearTOTSB2530);
  int binmaxy6=evalMomOri(k6+bw6,histSmearTOTSB2530);
  auto projection6 = (TH2F*)((TH2F*)histSmearTOTSB2530->ProjectionX("ll", binminy6, binmaxy6-1));
  projection6->SetFillColor(kBlue+2);
  DreamPlot::SetStyleHisto(projection6);
  outfile->cd();
  projection6->Write("projjj6");

  //int maxb6=projection6->GetMaximumBin();
  //double max6=projection6->GetBinCenter(maxb6);
  //int biny6=evalMomOri(k6,histSmearTOTSB2530);
  //double max6=evalMomTrafo(biny6,histSmearTOTSB2530);

  int max6=projection6->GetMean();
  double ex6=GetProjWidth(projection6)/2;
  double exvol6=GetProjWidth(projection6);
  double ex6_l=GetExleft(projection6,max6);
  double ex6_r=GetExright(projection6,max6);

  cout<<"halferror: "<<ex6<<"  errorleft: "<<ex6_l<<"  errorright: "<<ex6_r<<endl;

  exl[6] =ex6_l;
  exh[6] =ex6_r;
  //double ex6=GetWidth(biny6,histSmearTOTSB2530)/2;
  x[6]=max6;
  ex[6]=ex6;

  y[7]=SB3035->GetBinContent(i);
  ey[7]=SB3035->GetBinError(i);
  double k7= SB3035->GetBinCenter(i);
  int bw7=(SBleft->GetBinWidth(i))/2;
  int binminy7=evalMomOri(k7-bw7,histSmearTOTSB3035);
  int binmaxy7=evalMomOri(k7+bw7,histSmearTOTSB3035);
  auto projection7 = (TH2F*)((TH2F*)histSmearTOTSB3035->ProjectionX("ll", binminy7, binmaxy7-1));
  projection2->SetFillColor(kBlue+2);
  DreamPlot::SetStyleHisto(projection7);
  outfile->cd();
  projection7->Write("proj7");

  //int maxb7=projection7->GetMaximumBin();
  //double max7=projection7->GetBinCenter(maxb7);
  //int biny7=evalMomOri(k7,histSmearTOTSB3035);
  //double max7=evalMomTrafo(biny7,histSmearTOTSB3035);

  int max7=projection7->GetMean();
  double ex7=GetProjWidth(projection7)/2;
  double exvol7=GetProjWidth(projection7);
  double ex7_l=GetExleft(projection7,max7);
  double ex7_r=GetExright(projection7,max7);

  cout<<"halferror: "<<ex7<<"  errorleft: "<<ex7_l<<"  errorright: "<<ex7_r<<endl;

  exl[7] =ex7_l;
  exh[7] =ex7_r;
  //double ex7=GetWidth(biny7,histSmearTOTSB3035)/2;
  x[7]=max7;
  ex[7]=ex7;

  y[8]=SB3540->GetBinContent(i);
  ey[8]=SB3540->GetBinError(i);
  double k8= SB3540->GetBinCenter(i);
  int bw8=(SBleft->GetBinWidth(i))/2;
  int binminy8=evalMomOri(k8-bw8,histSmearTOTSB3540);
  int binmaxy8=evalMomOri(k8+bw8,histSmearTOTSB3540);
  auto projection8 = (TH2F*)((TH2F*)histSmearTOTSB3540->ProjectionX("ll", binminy8, binmaxy8-1));
  projection8->SetFillColor(kBlue+2);
  DreamPlot::SetStyleHisto(projection8);
  outfile->cd();
  projection8->Write("proj8");

  //int maxb8=projection8->GetMaximumBin();
  //double max8=projection8->GetBinCenter(maxb8);
  //int biny8=evalMomOri(k8,histSmearTOTSB3540);
  //double max8=evalMomTrafo(biny8,histSmearTOTSB3540);

  int max8=projection8->GetMean();
  double ex8=GetProjWidth(projection8)/2;
  double exvol8=GetProjWidth(projection8);
  double ex8_l=GetExleft(projection8,max8);
  double ex8_r=GetExright(projection8,max8);

  cout<<"halferror: "<<ex8<<"  errorleft: "<<ex8_l<<"  errorright: "<<ex8_r<<endl;

  exl[8] =ex8_l;
  exh[8] =ex8_r;
  //double ex8=GetWidth(biny8,histSmearTOTSB3540)/2;
  x[8]=max8;
  ex[8]=ex8;



  y[9]=SB4045->GetBinContent(i);
  ey[9]=SB4045->GetBinError(i);
  double k9= SB4045->GetBinCenter(i);
  int bw9=(SBleft->GetBinWidth(i))/2;
  int binminy9=evalMomOri(k9-bw9,histSmearTOTSB4045);
  int binmaxy9=evalMomOri(k9+bw9,histSmearTOTSB4045);
  auto projection9 = (TH2F*)((TH2F*)histSmearTOTSB4045->ProjectionX("ll", binminy9, binmaxy9-1));
  projection9->SetFillColor(kBlue+2);
  DreamPlot::SetStyleHisto(projection9);
  outfile->cd();
  projection9->Write("proj9");

  //int maxb9=projection9->GetMaximumBin();
  //double max9=projection9->GetBinCenter(maxb9);
  //int biny9=evalMomOri(k9,histSmearTOTSB4045);
  //double max9=evalMomTrafo(biny9,histSmearTOTSB4045);

  int max9=projection9->GetMean();
  double ex9=GetProjWidth(projection9)/2;
  double exvol9=GetProjWidth(projection9);
  double ex9_l=GetExleft(projection9,max9);
  double ex9_r=GetExright(projection9,max9);

  cout<<"halferror: "<<ex9<<"  errorleft: "<<ex9_l<<"  errorright: "<<ex9_r<<endl;

  exl[9] =ex9_l;
  exh[9] =ex9_r;
  //double ex9=GetWidth(biny9,histSmearTOTSB4045)/2;
  x[9]=max9;
  ex[9]=ex9;




  y[10]=SB4550->GetBinContent(i);
  ey[10]=SB4550->GetBinError(i);
  double k10= SB4550->GetBinCenter(i);
  int bw10=(SBleft->GetBinWidth(i))/2;
  int binminy10=evalMomOri(k10-bw10,histSmearTOTSB4550);
  int binmaxy10=evalMomOri(k10+bw10,histSmearTOTSB4550);
  auto projection10 = (TH2F*)((TH2F*)histSmearTOTSB4550->ProjectionX("ll", binminy10, binmaxy10-1));
  projection10->SetFillColor(kBlue+2);
  DreamPlot::SetStyleHisto(projection10);
  outfile->cd();
  projection10->Write("proj10");

  //int maxb10=projection10->GetMaximumBin();
  //double max10=projection10->GetBinCenter(maxb10);
  //int biny10=evalMomOri(k10,histSmearTOTSB4550);
  //double max10=evalMomTrafo(biny10,histSmearTOTSB4550);

  int max10=projection10->GetMean();
  double ex10=GetProjWidth(projection10)/2;
  double exvol10=GetProjWidth(projection10);
  double ex10_l=GetExleft(projection10,max10);
  double ex10_r=GetExright(projection10,max10);

  cout<<"halferror: "<<ex10<<"  errorleft: "<<ex10_l<<"  errorright: "<<ex10_r<<endl;

  exl[10] =ex10_l;
  exh[10] =ex10_r;
  //double ex10=GetWidth(biny10,histSmearTOTSB4550)/2;
  x[10]=max10;
  ex[10]=ex10;



  y[11]=SB5055->GetBinContent(i);
  ey[11]=SB5055->GetBinError(i);
  double k11= SB5055->GetBinCenter(i);
  int bw11=(SBleft->GetBinWidth(i))/2;
  int binminy11=evalMomOri(k11-bw11,histSmearTOTSB5055);
  int binmaxy11=evalMomOri(k11+bw11,histSmearTOTSB5055);
  auto projection11 = (TH2F*)((TH2F*)histSmearTOTSB5055->ProjectionX("ll", binminy11, binmaxy11-1));
  projection11->SetFillColor(kBlue+2);
  DreamPlot::SetStyleHisto(projection11);
  outfile->cd();
  projection11->Write("proj11");

  //int maxb11=projection11->GetMaximumBin();
  //double max11=projection11->GetBinCenter(maxb11);
  //int biny11=evalMomOri(k11,histSmearTOTSB5055);
  //double max11=evalMomTrafo(biny11,histSmearTOTSB5055);

  int max11=projection11->GetMean();
  double ex11=GetProjWidth(projection11)/2;
  double exvol11=GetProjWidth(projection11);
  double ex11_l=GetExleft(projection11,max11);
  double ex11_r=GetExright(projection11,max11);

  cout<<"halferror: "<<ex11<<"  errorleft: "<<ex11_l<<"  errorright: "<<ex11_r<<endl;

  exl[11] =ex11_l;
  exh[11] =ex11_r;
  //double ex11=GetWidth(biny11,histSmearTOTSB5055)/2;
  x[11]=max11;
  ex[11]=ex11;



  y[12]=SB5560->GetBinContent(i);
  ey[12]=SB5560->GetBinError(i);
  double k12= SB5560->GetBinCenter(i);
  int bw12=(SBleft->GetBinWidth(i))/2;
  int binminy12=evalMomOri(k12-bw12,histSmearTOTSB5560);
  int binmaxy12=evalMomOri(k12+bw12,histSmearTOTSB5560);
  auto projection12 = (TH2F*)((TH2F*)histSmearTOTSB5560->ProjectionX("ll", binminy12, binmaxy12-1));
  projection12->SetFillColor(kBlue+2);
  DreamPlot::SetStyleHisto(projection12);
  outfile->cd();
  projection12->Write("proj12");

  //int maxb12=projection12->GetMaximumBin();
  //double max12=projection12->GetBinCenter(maxb12);
  //int biny12=evalMomOri(k12,histSmearTOTSB5560);
  //double max12=evalMomTrafo(biny12,histSmearTOTSB5560);

  int max12=projection12->GetMean();
  double ex12=GetProjWidth(projection12)/2;
  double exvol12=GetProjWidth(projection12);
  double ex12_l=GetExleft(projection12,max12);
  double ex12_r=GetExright(projection12,max12);

  cout<<"halferror: "<<ex12<<"  errorleft: "<<ex12_l<<"  errorright: "<<ex12_r<<endl;

  exl[12] =ex12_l;
  exh[12] =ex12_r;
  //double ex12=GetWidth(biny12,histSmearTOTSB5560)/2;
  x[12]=max12;
  ex[12]=ex12;



  y[13]=SB67->GetBinContent(i);
  ey[13]=SB67->GetBinError(i);
  double k13= SB67->GetBinCenter(i);
  int bw13=(SBleft->GetBinWidth(i))/2;
  int binminy13=evalMomOri(k13-bw13,histSmearTOTSB67);
  int binmaxy13=evalMomOri(k13+bw13,histSmearTOTSB67);
  auto projection13 = (TH2F*)((TH2F*)histSmearTOTSB67->ProjectionX("ll", binminy13, binmaxy13-1));
  projection13->SetFillColor(kBlue+2);
  DreamPlot::SetStyleHisto(projection13 );
  outfile->cd();
  projection13 ->Write("proj13 ");

  //int maxb13 =projection13 ->GetMaximumBin();
  //double max13 =projection13 ->GetBinCenter(maxb13 );
  //int biny13=evalMomOri(k13,histSmearTOTSB67);
  //double max13=evalMomTrafo(biny13,histSmearTOTSB67);
  int max13=projection13->GetMean();
  double ex13=GetProjWidth(projection13)/2;
  double exvol13=GetProjWidth(projection13);
  double ex13_l=GetExleft(projection13,max13);
  double ex13_r=GetExright(projection13,max13);

  cout<<"halferror: "<<ex13<<"  errorleft: "<<ex13_l<<"  errorright: "<<ex13_r<<endl;

  exl[13] =ex13_l;
  exh[13] =ex13_r;
  //double ex13=GetWidth(biny13,histSmearTOTSB67)/2;
  x[13]=max13;
  ex[13]=ex13;



  y[14]=SB78->GetBinContent(i);
  ey[14]=SB78->GetBinError(i);
  double k14= SB78->GetBinCenter(i);
  int bw14=(SBleft->GetBinWidth(i))/2;
  int binminy14=evalMomOri(k14-bw14,histSmearTOTSB78);
  int binmaxy14=evalMomOri(k14+bw14,histSmearTOTSB78);
  auto projection14 = (TH2F*)((TH2F*)histSmearTOTSB78->ProjectionX("ll", binminy14, binmaxy14-1));
  projection14->SetFillColor(kBlue+2);
  DreamPlot::SetStyleHisto(projection14);
  outfile->cd();
  projection14->Write("proj14");

  //int maxb14=projection14->GetMaximumBin();
  //double max14=projection14->GetBinCenter(maxb14);
  //int biny14=evalMomOri(k14,histSmearTOTSB78);
  //double max14=evalMomTrafo(biny14,histSmearTOTSB78);

  int max14=projection14->GetMean();
  double ex14=GetProjWidth(projection14)/2;
  double exvol14=GetProjWidth(projection14);
  double ex14_l=GetExleft(projection14,max14);
  double ex14_r=GetExright(projection14,max14);

  cout<<"halferror: "<<ex14<<"  errorleft: "<<ex14_l<<"  errorright: "<<ex14_r<<endl;

  exl[14] =ex14_l;
  exh[14] =ex14_r;
  //double ex14=GetWidth(biny14,histSmearTOTSB78)/2;
  x[14]=max14;
  ex[14]=ex14;






  y[15]=SB89->GetBinContent(i);
  ey[15]=SB89->GetBinError(i);
  double k15= SB89->GetBinCenter(i);
  int bw15=(SBleft->GetBinWidth(i))/2;
  int binminy15=evalMomOri(k15-bw15,histSmearTOTSB89);
  int binmaxy15=evalMomOri(k15+bw15,histSmearTOTSB89);
  auto projection15 = (TH2F*)((TH2F*)histSmearTOTSB89->ProjectionX("ll", binminy15, binmaxy15-1));
  projection15->SetFillColor(kBlue+2);
  DreamPlot::SetStyleHisto(projection15);
  outfile->cd();
  projection15->Write("proj15");

  //int maxb15=projection15->GetMaximumBin();
  //double max15=projection15->GetBinCenter(maxb15);
  //int biny15=evalMomOri(k15,histSmearTOTSB89);
  //double max15=evalMomTrafo(biny15,histSmearTOTSB89);

  int max15=projection15->GetMean();
  double ex15=GetProjWidth(projection15)/2;
  double exvol15=GetProjWidth(projection15);
  double ex15_l=GetExleft(projection15,max15);
  double ex15_r=GetExright(projection15,max15);

  cout<<"halferror: "<<ex15<<"  errorleft: "<<ex15_l<<"  errorright: "<<ex15_r<<endl;

  exl[15] =ex15_l;
  exh[15] =ex15_r;
  //double ex15=GetWidth(biny15,histSmearTOTSB89)/2;
  x[15]=max15;
  ex[15]=ex15;


  y[16]=SB90->GetBinContent(i);
  ey[16]=SB90->GetBinError(i);
  double k16= SB90->GetBinCenter(i);
  int bw16=(SBleft->GetBinWidth(i))/2;
  int binminy16=evalMomOri(k16-bw16,histSmearTOTSB90);
  int binmaxy16=evalMomOri(k16+bw16,histSmearTOTSB90);
  auto projection16 = (TH2F*)((TH2F*)histSmearTOTSB90->ProjectionX("ll", binminy16, binmaxy16-1));
  projection16->SetFillColor(kBlue+2);
  DreamPlot::SetStyleHisto(projection16);
  outfile->cd();
  projection16->Write("proj16");

  //int maxb16=projection16->GetMaximumBin();
  //double max16=projection16->GetBinCenter(maxb16);
  //int biny16=evalMomOri(k16,histSmearTOTSB90);
  //double max16=evalMomTrafo(biny16,histSmearTOTSB90);

  int max16=projection16->GetMean();
  double ex16=GetProjWidth(projection16)/2;
  double exvol16=GetProjWidth(projection16);
  double ex16_l=GetExleft(projection16,max16);
  double ex16_r=GetExright(projection16,max16);

  cout<<"halferror: "<<ex16<<"  errorleft: "<<ex16_l<<"  errorright: "<<ex16_r<<endl;

  exl[16] =ex16_l;
  exh[16] =ex16_r;
  //double ex16=GetWidth(biny16,histSmearTOTSB90)/2;
  x[16]=max16;
  ex[16]=ex16;


  y[17]=SB2223->GetBinContent(i);
  ey[17]=SB2223->GetBinError(i);
  double k17= SB2223->GetBinCenter(i);
  int bw17=(SBleft->GetBinWidth(i))/2;
  int binminy17=evalMomOri(k17-bw17,histSmearTOTSB2223);
  int binmaxy17=evalMomOri(k17+bw17,histSmearTOTSB2223);
  auto projection17 = (TH2F*)((TH2F*)histSmearTOTSB2223->ProjectionX("ll", binminy17, binmaxy17-1));
  projection17->SetFillColor(kBlue+2);
  DreamPlot::SetStyleHisto(projection17);
  outfile->cd();
  projection17->Write("projjj17");

  //int maxb17=projection17->GetMaximumBin();
  //double max17=projection17->GetBinCenter(maxb17);
  //int biny17=evalMomOri(k17,histSmearTOTSB2223);
  //double max17=evalMomTrafo(biny17,histSmearTOTSB2223);

  int max17=projection17->GetMean();
  double ex17=GetProjWidth(projection17)/2;
  double exvol17=GetProjWidth(projection17);
  double ex17_l=GetExleft(projection17,max17);
  double ex17_r=GetExright(projection17,max17);

  cout<<"halferror: "<<ex17<<"  errorleft: "<<ex17_l<<"  errorright: "<<ex17_r<<endl;

  exl[17] =ex17_l;
  exh[17] =ex17_r;
  //double ex17=GetWidth(biny17,histSmearTOTSB2223)/2;
  x[17]=max17;
  ex[17]=ex17;

  yy[1]=y[17];
  xx[1]=x[17];
  eyy[1]=ey[17];
  exx[1]=ex[17];

  y[18]=SB2425->GetBinContent(i);
  ey[18]=SB2425->GetBinError(i);
  double k18= SB2425->GetBinCenter(i);
  int bw18=(SBleft->GetBinWidth(i))/2;
  int binminy18=evalMomOri(k18-bw18,histSmearTOTSB2425);
  int binmaxy18=evalMomOri(k18+bw18,histSmearTOTSB2425);
  auto projection18 = (TH2F*)((TH2F*)histSmearTOTSB2425->ProjectionX("ll", binminy18, binmaxy18-1));
  projection18->SetFillColor(kBlue+2);
  DreamPlot::SetStyleHisto(projection18);
  outfile->cd();
  projection18->Write("proj18");

  //int maxb18=projection18->GetMaximumBin();
  //double max18=projection18->GetBinCenter(maxb18);
  //int biny18=evalMomOri(k18,histSmearTOTSB2425);
  //double max18=evalMomTrafo(biny18,histSmearTOTSB2425);

  int max18=projection18->GetMean();
  double ex18=GetProjWidth(projection18)/2;
  double exvol18=GetProjWidth(projection18);
  double ex18_l=GetExleft(projection18,max18);
  double ex18_r=GetExright(projection18,max18);

  cout<<"halferror: "<<ex18<<"  errorleft: "<<ex18_l<<"  errorright: "<<ex18_r<<endl;

  exl[18] =ex18_l;
  exh[18] =ex18_r;
  //double ex18=GetWidth(biny18,histSmearTOTSB2425)/2;
  x[18]=max18;
  ex[18]=ex18;

  yy[2]=y[18];
  xx[2]=x[18];
  eyy[2]=ey[18];
  exx[2]=ex[18];

  y[19]=SB2627->GetBinContent(i);
  ey[19]=SB2627->GetBinError(i);
  double k19= SB2627->GetBinCenter(i);
  int bw19=(SBleft->GetBinWidth(i))/2;
  int binminy19=evalMomOri(k19-bw19,histSmearTOTSB2627);
  int binmaxy19=evalMomOri(k19+bw19,histSmearTOTSB2627);
  auto projection19 = (TH2F*)((TH2F*)histSmearTOTSB2627->ProjectionX("ll", binminy19, binmaxy19-1));
  projection19->SetFillColor(kBlue+2);
  DreamPlot::SetStyleHisto(projection19);
  outfile->cd();
  projection19->Write("proj19");

  //int maxb19=projection19->GetMaximumBin();
  //double max19=projection19->GetBinCenter(maxb19);
  //int biny19=evalMomOri(k19,histSmearTOTSB2627);
  //double max19=evalMomTrafo(biny19,histSmearTOTSB2627);

  int max19=projection19->GetMean();
  double ex19=GetProjWidth(projection19)/2;
  double exvol19=GetProjWidth(projection19);
  double ex19_l=GetExleft(projection19,max19);
  double ex19_r=GetExright(projection19,max19);

  cout<<"halferror: "<<ex19<<"  errorleft: "<<ex19_l<<"  errorright: "<<ex19_r<<endl;

  exl[19] =ex19_l;
  exh[19] =ex19_r;
  //double ex19=GetWidth(biny19,histSmearTOTSB2627)/2;
  x[19]=max19;
  ex[19]=ex19;

  yy[3]=y[19];
  xx[3]=x[19];
  eyy[3]=ey[19];
  exx[3]=ex[19];

  y[20]=SB2829->GetBinContent(i);
  ey[20]=SB2829->GetBinError(i);
  double k20= SB2829->GetBinCenter(i);
  int bw20=(SBleft->GetBinWidth(i))/2;
  int binminy20=evalMomOri(k20-bw20,histSmearTOTSB2829);
  int binmaxy20=evalMomOri(k20+bw20,histSmearTOTSB2829);
  auto projection20 = (TH2F*)((TH2F*)histSmearTOTSB2829->ProjectionX("ll", binminy20, binmaxy20-1));
  projection20->SetFillColor(kBlue+2);
  DreamPlot::SetStyleHisto(projection20);
  outfile->cd();
  projection20->Write("proj20");

  //int maxb20=projection20->GetMaximumBin();
  //double max20=projection20->GetBinCenter(maxb20);
  //int biny20=evalMomOri(k20,histSmearTOTSB2829);
  //double max20=evalMomTrafo(biny20,histSmearTOTSB2829);

  int max20=projection20->GetMean();
  double ex20=GetProjWidth(projection20)/2;
  double exvol20=GetProjWidth(projection20);
  double ex20_l=GetExleft(projection20,max20);
  double ex20_r=GetExright(projection20,max20);

  cout<<"halferror: "<<ex20<<"  errorleft: "<<ex20_l<<"  errorright: "<<ex20_r<<endl;

  exl[20] =ex20_l;
  exh[20] =ex20_r;
  //double ex20=GetWidth(biny20,histSmearTOTSB2829)/2;
  x[20]=max20;
  ex[20]=ex20;

  yy[4]=y[20];
  xx[4]=x[20];
  eyy[4]=ey[20];
  exx[4]=ex[20];


  y[21]=SB3031->GetBinContent(i);
  ey[21]=SB3031->GetBinError(i);
  double k21= SB3031->GetBinCenter(i);
  int bw21=(SBleft->GetBinWidth(i))/2;
  int binminy21=evalMomOri(k21-bw21,histSmearTOTSB3031);
  int binmaxy21=evalMomOri(k21+bw21,histSmearTOTSB3031);
  auto projection21 = (TH2F*)((TH2F*)histSmearTOTSB3031->ProjectionX("ll", binminy21, binmaxy21-1));
  projection21->SetFillColor(kBlue+2);
  DreamPlot::SetStyleHisto(projection21);
  outfile->cd();
  projection21->Write("proj21");

  //int maxb21=projection21->GetMaximumBin();
  //double max21=projection21->GetBinCenter(maxb21);
  //int biny21=evalMomOri(k21,histSmearTOTSB3031);
  //double max21=evalMomTrafo(biny21,histSmearTOTSB3031);

  int max21=projection21->GetMean();
  double ex21=GetProjWidth(projection21)/2;
  double exvol21=GetProjWidth(projection21);
  double ex21_l=GetExleft(projection21,max21);
  double ex21_r=GetExright(projection21,max21);

  cout<<"halferror: "<<ex21<<"  errorleft: "<<ex21_l<<"  errorright: "<<ex21_r<<endl;

  exl[21] =ex21_l;
  exh[21] =ex21_r;
  //double ex21=GetWidth(biny21,histSmearTOTSB3031)/2;
  x[21]=max21;
  ex[21]=ex21;

  yy[5]=y[21];
  xx[5]=x[21];
  eyy[5]=ey[21];
  exx[5]=ex[21];

  y[22]=SB3334->GetBinContent(i);
  ey[22]=SB3334->GetBinError(i);
  double k22= SB3334->GetBinCenter(i);
  int bw22=(SBleft->GetBinWidth(i))/2;
  int binminy22=evalMomOri(k22-bw22,histSmearTOTSB3334);
  int binmaxy22=evalMomOri(k22+bw22,histSmearTOTSB3334);
  auto projection22 = (TH2F*)((TH2F*)histSmearTOTSB3334->ProjectionX("ll", binminy22, binmaxy22-1));
  projection22->SetFillColor(kBlue+2);
  DreamPlot::SetStyleHisto(projection22);
  outfile->cd();
  projection22->Write("proj22");

  //int maxb22=projection22->GetMaximumBin();
  //double max22=projection22->GetBinCenter(maxb22);
  //int biny22=evalMomOri(k22,histSmearTOTSB3334);
  //double max22=evalMomTrafo(biny22,histSmearTOTSB3334);

  int max22=projection22->GetMean();
  double ex22=GetProjWidth(projection22)/2;
  double exvol22=GetProjWidth(projection22);
  double ex22_l=GetExleft(projection22,max22);
  double ex22_r=GetExright(projection22,max22);

  cout<<"halferror: "<<ex22<<"  errorleft: "<<ex22_l<<"  errorright: "<<ex22_r<<endl;

  exl[22] =ex22_l;
  exh[22] =ex22_r;
  //double ex22=GetWidth(biny22,histSmearTOTSB3334)/2;
  x[22]=max22;
  ex[22]=ex22;

  yy[6]=y[22];
  xx[6]=x[22];
  eyy[6]=ey[22];
  exx[6]=ex[22];

  y[23]=SB3738->GetBinContent(i);
  ey[23]=SB3738->GetBinError(i);
  double k23= SB3738->GetBinCenter(i);
  int bw23=(SBleft->GetBinWidth(i))/2;
  int binminy23=evalMomOri(k23-bw23,histSmearTOTSB3738);
  int binmaxy23=evalMomOri(k23+bw23,histSmearTOTSB3738);
  auto projection23 = (TH2F*)((TH2F*)histSmearTOTSB3738->ProjectionX("ll", binminy23, binmaxy23-1));
  projection23->SetFillColor(kBlue+2);
  DreamPlot::SetStyleHisto(projection23);
  outfile->cd();
  projection23->Write("proj23");

  //int maxb23=projection23->GetMaximumBin();
  //double max23=projection23->GetBinCenter(maxb23);
  //int biny23=evalMomOri(k23,histSmearTOTSB3738);
  //double max23=evalMomTrafo(biny23,histSmearTOTSB3738);

  int max23=projection23->GetMean();
  double ex23=GetProjWidth(projection23)/2;
  double exvol23=GetProjWidth(projection23);
  double ex23_l=GetExleft(projection23,max23);
  double ex23_r=GetExright(projection23,max23);

  cout<<"halferror: "<<ex23<<"  errorleft: "<<ex23_l<<"  errorright: "<<ex23_r<<endl;

  exl[23] =ex23_l;
  exh[23] =ex23_r;
  //double ex23=GetWidth(biny23,histSmearTOTSB3738)/2;
  x[23]=max23;
  ex[23]=ex23;


  yy[7]=y[23];
  xx[7]=x[23];
  eyy[7]=ey[23];
  exx[7]=ex[23];

  y[24]=SB4041->GetBinContent(i);
  ey[24]=SB4041->GetBinError(i);
  double k24= SB4041->GetBinCenter(i);
  int bw24=(SBleft->GetBinWidth(i))/2;
  int binminy24=evalMomOri(k24-bw24,histSmearTOTSB4041);
  int binmaxy24=evalMomOri(k24+bw24,histSmearTOTSB4041);
  auto projection24 = (TH2F*)((TH2F*)histSmearTOTSB4041->ProjectionX("ll", binminy24, binmaxy24-1));
  projection24->SetFillColor(kBlue+2);
  DreamPlot::SetStyleHisto(projection24);
  outfile->cd();
  projection24->Write("proj24");

  //int maxb24=projection24->GetMaximumBin();
  //double max24=projection24->GetBinCenter(maxb24);
  //int biny24=evalMomOri(k24,histSmearTOTSB4041);
  //double max24=evalMomTrafo(biny24,histSmearTOTSB4041);

  int max24=projection24->GetMean();
  double ex24=GetProjWidth(projection24)/2;
  double exvol24=GetProjWidth(projection24);
  double ex24_l=GetExleft(projection24,max24);
  double ex24_r=GetExright(projection24,max24);

  cout<<"halferror: "<<ex24<<"  errorleft: "<<ex24_l<<"  errorright: "<<ex24_r<<endl;

  exl[24] =ex24_l;
  exh[24] =ex24_r;
  //double ex24=GetWidth(biny24,histSmearTOTSB4041)/2;
  x[24]=max24;
  ex[24]=ex24;

  yy[8]=y[24];
  xx[8]=x[24];
  eyy[8]=ey[24];
  exx[8]=ex[24];

  //y[666]=SB12->GetBinContent(i);
  //ey[666]=SB12->GetBinError(i);
  //double k666= SB12->GetBinCenter(i);
  //int biny666=evalMomOri(k666,histSmearTOTSB12);
  //double max666=evalMomTrafo(biny666,histSmearTOTSB12);
  //double ex666=GetWidth(biny666,histSmearTOTSB12)/2;
  //x[666]=max666;
  //ex[666]=ex666;








  auto pKCF = new TGraphAsymmErrors(n,x,y,exl, exh,ey, ey);
  int lineWidth2 = 3;
  DreamPlot::SetStyleGraph(pKCF, 20, kBlack, 0.7);
  //pKCF->SetLineWidth(lineWidth2);


  //auto pKCF = new TGraphErrors(n,x,y,ex,ey);
  //int lineWidth2 = 3;
  //DreamPlot::SetStyleGraph(pKCF, 20, kBlack, 0.7);
  //pKCF->SetLineWidth(lineWidth2);


  auto pKCF2 = new TGraphErrors(n,x,y,ex,ey);
  DreamPlot::SetStyleGraph(pKCF2, 20, kBlack, 0.7);
  pKCF2->SetLineWidth(lineWidth2);



  auto pKBG = new TGraphErrors(n1,xx,yy,exx,eyy);
  DreamPlot::SetStyleGraph(pKBG, 20, kBlack, 0.7);
  pKBG->SetLineWidth(lineWidth2);

  auto pKCF666 = new TGraphAsymmErrors(n,x,y,exl, exh,ey, ey);
  DreamPlot::SetStyleGraph(pKCF666, 20, kBlack, 0.7);
  pKCF666->SetLineWidth(lineWidth2);

  pKCF666->SetLineColor(4);
  // pKCF->SetLineWidth(4);
  pKCF666->SetMarkerColor(4);
  pKCF666->SetMarkerSize(1.5);
  pKCF666->SetMarkerStyle(21);




  // fitfun->Draw("L3same");

  for (int a = 0; a < 11; a++) {
    cout << x[a] << " " << ex[a] << " " << y[a] << " " << ey[a] << endl;
  }

  // SetCFPoint(int i, TH1F* CF, TH2F* matrix, int binnr, TGraphErrors Graph )

  pKCF->SetLineColor(4);
  // pKCF->SetLineWidth(4);
  pKCF->SetMarkerColor(4);
  pKCF->SetMarkerSize(1.5);
  pKCF->SetMarkerStyle(21);

  pKCF2->SetLineColor(4);
  // pKCF->SetLineWidth(4);
  pKCF2->SetMarkerColor(4);
  pKCF2->SetMarkerSize(1.5);
  pKCF2->SetMarkerStyle(21);

  pKBG->SetLineColor(4);
  // pKCF->SetLineWidth(4);
  pKBG->SetMarkerColor(4);
  pKBG->SetMarkerSize(1.5);
  pKBG->SetMarkerStyle(21);

  auto smearedCF = GetSmearedCF(pKCF, histSmearTOTPeak);
  smearedCF->SetLineColor(2);
  // smearedCF->SetLineWidth(4);
  smearedCF->SetMarkerColor(2);
  smearedCF->SetMarkerSize(1.5);
  smearedCF->SetMarkerStyle(21);

  auto smearedSBleft = GetSmearedCF(pKCF, histSmearTOTSBleft);
  smearedSBleft->SetLineColor(2);
  // smearedCF->SetLineWidth(4);
  smearedSBleft->SetMarkerColor(2);
  smearedSBleft->SetMarkerSize(1.5);
  smearedSBleft->SetMarkerStyle(21);

  auto smearedSBright = GetSmearedCF(pKCF, histSmearTOTSBright);
  smearedSBright->SetLineColor(2);
  // smearedCF->SetLineWidth(4);
  smearedSBright->SetMarkerColor(2);
  smearedSBright->SetMarkerSize(1.5);
  smearedSBright->SetMarkerStyle(21);

  auto smearedSB12 = GetSmearedCF(pKCF, histSmearTOTSB12);
  smearedSB12->SetLineColor(2);
  // smearedCF->SetLineWidth(4);
  smearedSB12->SetMarkerColor(2);
  smearedSB12->SetMarkerSize(1.5);
  smearedSB12->SetMarkerStyle(21);

  auto smearedSB23 = GetSmearedCF(pKCF, histSmearTOTSB23);
  smearedSB23->SetLineColor(2);
  // smearedCF->SetLineWidth(4);
  smearedSB23->SetMarkerColor(2);
  smearedSB23->SetMarkerSize(1.5);
  smearedSB23->SetMarkerStyle(21);

  auto smearedSB34 = GetSmearedCF(pKCF, histSmearTOTSB34);
  smearedSB34->SetLineColor(2);
  // smearedCF->SetLineWidth(4);
  smearedSB34->SetMarkerColor(2);
  smearedSB34->SetMarkerSize(1.5);
  smearedSB34->SetMarkerStyle(21);

  auto smearedSB45 = GetSmearedCF(pKCF, histSmearTOTSB45);
  smearedSB45->SetLineColor(2);
  // smearedCF->SetLineWidth(4);
  smearedSB45->SetMarkerColor(2);
  smearedSB45->SetMarkerSize(1.5);
  smearedSB45->SetMarkerStyle(21);

  auto smearedSB56 = GetSmearedCF(pKCF, histSmearTOTSB56);
  smearedSB56->SetLineColor(2);
  // smearedCF->SetLineWidth(4);
  smearedSB56->SetMarkerColor(2);
  smearedSB56->SetMarkerSize(1.5);
  smearedSB56->SetMarkerStyle(21);

  auto smearedSB67 = GetSmearedCF(pKCF, histSmearTOTSB67);
  smearedSB67->SetLineColor(2);
  // smearedCF->SetLineWidth(4);
  smearedSB67->SetMarkerColor(2);
  smearedSB67->SetMarkerSize(1.5);
  smearedSB67->SetMarkerStyle(21);

  auto smearedSB78 = GetSmearedCF(pKCF, histSmearTOTSB78);
  smearedSB78->SetLineColor(2);
  // smearedCF->SetLineWidth(4);
  smearedSB78->SetMarkerColor(2);
  smearedSB78->SetMarkerSize(1.5);
  smearedSB78->SetMarkerStyle(21);

  auto smearedSB89 = GetSmearedCF(pKCF, histSmearTOTSB89);
  smearedSB89->SetLineColor(2);
  // smearedCF->SetLineWidth(4);
  smearedSB89->SetMarkerColor(2);
  smearedSB89->SetMarkerSize(1.5);
  smearedSB89->SetMarkerStyle(21);

  auto smearedSB90 = GetSmearedCF(pKCF, histSmearTOTSB90);
  smearedSB90->SetLineColor(2);
  // smearedCF->SetLineWidth(4);
  smearedSB90->SetMarkerColor(2);
  smearedSB90->SetMarkerSize(1.5);
  smearedSB90->SetMarkerStyle(21);

//  Double_t par[6];
//  TF1* f1 = new TF1("f1", "gaus", 0, 450);
//  TF1* f2 = new TF1("f2", "pol2", 0, 850);

//  TF1* total = new TF1("total", "gaus(0)+pol2(3)", 0, 850);

//  f1->SetParLimits(3, 90, 98);  // 0.072,0.95

//  f2->SetParLimits(0, 0.85, 1.00);  // 0.072,0.95

//  total->SetLineColor(2);
//  pKCF->Fit(f1, "R");
//  pKBG->Fit(f2, "R");

//  f1->GetParameters(&par[0]);
//  f2->GetParameters(&par[3]);

//  total->SetParameters(par);
//  pKCF->Fit(total, "F", "", 0, 900);

//  float pg0 = total->GetParameter(0);
//  float pg1 = total->GetParameter(1);
//  float pg2 = total->GetParameter(2);
//  float pg3 = total->GetParameter(3);
//  float pg4 = total->GetParameter(4);
//  float pg5 = total->GetParameter(5);

//  float eg0 = total->GetParError(0);
//  float eg1 = total->GetParError(1);
//  float eg2 = total->GetParError(2);
//  float eg3 = total->GetParError(3);
//  float eg4 = total->GetParError(4);
//  float eg5 = total->GetParError(5);

  // std::vector<float> param0={{pg0-eg0,pg0,pg0+eg0}};
  // std::vector<float> param1={{pg1-eg1,pg1,pg1+eg1}};
  // std::vector<float> param2={{pg2-eg2,pg2,pg2+eg2}};
  // std::vector<float> param3={{pg3-eg3,pg3,pg3+eg3}};
  // std::vector<float> param4={{pg4-eg4,pg4,pg4+eg4}};
  // std::vector<float> param5={{pg5-eg5,pg5,pg5+eg5}};


//-----------------------------------------------------old fit---------------------
//  Double_t par[6];
//  TF1 *f1 = new TF1("f1","gaus",0,450);
//  TF1 *f11 = new TF1("f11","gaus(0)+pol0(3)",0,850);
//  TF1 *f12 = new TF1("f12","gaus(0)+pol1(3)",0,850);
//  TF1 *f13 = new TF1("f13","gaus(0)+pol2(3)",0,850);

//  TF1 *f2 = new TF1("f2","pol2",0,850);

//  TF1 *total = new TF1("total","gaus(0)+pol2(3)",0,850);


//  pKCF->Fit(f1,"R");
//  f11->SetParameter(0,f1->GetParameter(0));
//  f11->SetParameter(1,f1->GetParameter(1));
//  f11->SetParameter(2,f1->GetParameter(2));
//  f11->SetParLimits(3,0.9,1.2); //0.072,0.95

//  pKCF->Fit(f11,"R");
//  f11->SetParLimits(0,0.08,0.10); //0.072,0.95
//  f11->SetParLimits(1,220,300); //0.072,0.95
//  pKCF->Fit(f11,"0", 0,850);


//  f12->SetParameter(0,f11->GetParameter(0));
//  f12->SetParameter(1,f11->GetParameter(1));
//  f12->SetParameter(2,f11->GetParameter(2));
//  f12->SetParameter(3,f11->GetParameter(3));
//  f12->SetParLimits(3,0.5,0.80); //0.072,0.95

//  pKCF->Fit(f12,"R");

//  f13->SetParameter(0,f12->GetParameter(0));
//  f13->SetParameter(1,f12->GetParameter(1));
//  f13->SetParameter(2,f12->GetParameter(2));
//  f13->SetParameter(3,f12->GetParameter(3));
//  f13->SetParLimits(3,0.6,0.80); //0.072,0.95
//  f13->SetParameter(4,f12->GetParameter(4));


//  pKCF->Fit(f13,"R");


//  pKCF->Fit(f13,"R");
//  f13->SetParLimits(4,0.0003,0.0007); //0.072,0.95
//  f13->SetParLimits(3,0.7,0.75); //0.072,0.95


//  pKCF->Fit(f13,"R");


//  f2->SetParLimits(0,0.85,1.00); //0.072,0.95

//  total->SetLineColor(2);


//  pKBG->Fit(f2,"R");



//  f13->GetParameters(&par[0]);


  //-----------------------------------------------------new fit---------------------


  cout<<"GAUSS FITTTTT----------------"<<endl;
  TF1 *f13 = new TF1("fgp","gaus(0)+pol2(3)",0,650);

  f13->SetParameter(0,0.214);
  //fgp->FixParameter(1,188.117);
  //fgp->FixParameter(2,141.988);
  f13->SetParLimits(1,175,195);
  f13->SetParLimits(2,120,160);
  f13->SetParameter(3,0.075);
  f13->SetParLimits(3,0.6,0.80); //0.072,0.95
  f13->SetParameter(4,0.00069999);
  f13->SetParLimits(4,0.0007,0.00085);
  f13->SetParLimits(5,-0.0000006,-0.0000005);
  //fgp->SetParameter(5,0.00069999);
  pKCF->Fit(f13,"R");
  pKCF->Fit(f13," "," ", 0,750);
  //fgp->FixParameter(1,188.117);
  cout<<"GAUSS FINALL----------------"<<endl;

  pKCF->Fit(f13," "," ", 0,850);




  float pg0 = f13->GetParameter(0);
  float pg1 = f13->GetParameter(1);
  float pg2 = f13->GetParameter(2);
  float pg3 = f13->GetParameter(3);
  float pg4 = f13->GetParameter(4);
  float pg5 = f13->GetParameter(5);

  float eg0 = f13->GetParError(0);
  float eg1 = f13->GetParError(1);
  float eg2 = f13->GetParError(2);
  float eg3 = f13->GetParError(3);
  float eg4 = f13->GetParError(4);
  float eg5 = f13->GetParError(5);

  float perer = 1.0;

  std::vector<float> param0 = {{pg0 - (eg0 * perer), pg0, pg0 + (eg0 * perer)}};
  std::vector<float> param1 = {{pg1 - (eg1 * perer), pg1, pg1 + (eg1 * perer)}};
  std::vector<float> param2 = {{pg2 - (eg2 * perer), pg2, pg2 + (eg2 * perer)}};
  std::vector<float> param3 = {{pg3 - (eg3 * perer), pg3, pg3 + (eg3 * perer)}};
  std::vector<float> param4 = {{pg4 - (eg4 * perer), pg4, pg4 + (eg4 * perer)}};
  std::vector<float> param5 = {{pg5 - (eg5 * perer), pg5, pg5 + (eg5 * perer)}};



  int histnr = 730;
  TF1* fits[histnr];
  int counter = 0;
  float chindfmax=2.0;
  float minimumchindf=10;
  int mincounter=0;
  float minTOTSB=10000.0;
  int minTOTcounterSB=0;
  float minTOT=10000.0;
  int minTOTcounter=0;
  float minTOTSBleft=10000.0;
  int minTOTcounterSBleft=0;


  outfile->cd();

  int ndfSB=0;
  int SBCFbins=SBrightc->GetNbinsX();
  for(int i=0;i<SBCFbins;i++){
  double momentum= SBrightc->GetBinCenter(i + 1);
  if(momentum<=1000){
      ndfSB++;
  }
  }

  int ndfSBleft=0;
  int SBCFbinsleft=SBleftc->GetNbinsX();
  for(int i=0;i<SBCFbinsleft;i++){
  double momentum= SBleftc->GetBinCenter(i + 1);
  if(momentum<=1000){
      ndfSBleft++;
  }
  }

  int ndfFIT=0;

  cout<<"SBNDFright: "<< ndfSB<<endl;
  cout<<"SBNDFleft: "<< ndfSBleft<<endl;

  int a = 0;
  for (int p0 = 0; p0 < param0.size(); p0++) {
    for (int p1 = 0; p1 < param1.size(); p1++) {
      for (int p2 = 0; p2 < param2.size(); p2++) {
        for (int p3 = 0; p3 < param3.size(); p3++) {
          for (int p4 = 0; p4 < param4.size(); p4++) {
            for (int p5 = 0; p5 < param5.size(); p5++) {

//     int p1=0;
//     int p2=0;
//     int p3=0;
//     int p4=0;
//     int p5=0;
              int TotalChi2=0;


              // while(counter<histnr){
              cout << "counter nr: " << counter << endl;

              fits[counter] = nullptr;
              fits[counter] =
                  new TF1(Form("fit_%d", counter), "gaus(0)+pol2(3)", 0, 850);
              fits[counter]->FixParameter(0, param0[p0]);
              fits[counter]->FixParameter(1, param1[p1]);
              fits[counter]->FixParameter(2, param2[p2]);
              fits[counter]->FixParameter(3, param3[p3]);
              fits[counter]->FixParameter(4, param4[p4]);
              fits[counter]->FixParameter(5, param5[p5]);

              pKCF->Fit(fits[counter]);

//              fits[counter]->SetLineColor(2);
//              // smearedCF->SetLineWidth(4);
//              fits[counter]->SetMarkerColor(2);
//              fits[counter]->SetMarkerSize(1.5);
//              fits[counter]->SetMarkerStyle(21);

             // fits[counter]->Write(Form("fit_%d", counter));


              fits[counter]->SetLineColorAlpha(kGreen, 0.2);
              double chi = fits[counter]->GetChisquare();
              double bin1 = fits[counter]->Eval(0);
              double ndf= fits[counter]->GetNDF();
              ndfFIT=ndf;
              double cndf= chi/ndf;
             // cout<<"chi2/ndf FIT fun: "<<cndf<<endl;

              if (cndf < chindfmax) {

                  fits[counter]->SetLineColor(2);
                  // smearedCF->SetLineWidth(4);
                  fits[counter]->SetMarkerColor(2);
                  fits[counter]->SetMarkerSize(1.5);
                  fits[counter]->SetMarkerStyle(21);

                  fits[counter]->Write(Form("chifit_%d", counter));

              }

              double para0=fits[counter]->GetParameter(0);
              double para1=fits[counter]->GetParameter(1);
              double para2=fits[counter]->GetParameter(2);
              double para3=fits[counter]->GetParameter(3);
              double para4=fits[counter]->GetParameter(4);
              double para5=fits[counter]->GetParameter(5);

              //float chiSBsmear=0;

              if ((cndf < chindfmax)) {
                  cout<<cndf<<endl;
                  const int binwidthg = 1;
                  const int NumMomBinsg = 800;

                  Double_t x2g[NumMomBinsg], y2g[NumMomBinsg];

                for (int i = 0; i < NumMomBinsg; i++) {
                  x2g[i] = i * binwidthg;
                  y2g[i] = para0 * TMath::Gaus(x2g[i], para1,
                                                    para2, false) +
                           para3 + para4 * x2g[i] +
                           para5 * x2g[i] * x2g[i];
                }
                TGraph* gausfit = new TGraph(NumMomBinsg, x2g, y2g);
                //gausfit->Write(Form("gaussfit_%d", counter));

                auto smearedfSB_CF = GetSmearedCF(gausfit, histSmearTOTSBright);

               // smearedfSB_CF->Write(Form("smearedSB_%d", counter));
                double IntegralCFright =integrate(smearedfSB_CF,800,900);
                Double_t scaleright = 1/IntegralCFright;

                auto normSBright=GetSmearedCFscaled(gausfit, histSmearTOTSBright, scaleright);
                normSBright->SetLineColor(2);
                normSBright->SetMarkerColor(2);
                normSBright->SetMarkerSize(1.5);
                normSBright->SetMarkerStyle(21);
                normSBright->Write(Form("normsmearedSB_%d", counter));
                int ndfsbright=ndfSB;
                float chiSBsmear= getchisquared(SBrightc,normSBright);
                float chisqTOTSB=chiSBsmear/ndfsbright;
                cout<<"chi2 pK: "<<chi<<endl;
                cout<<"chi2/ndf pK: "<<cndf<<endl;
                cout<<"chi2 SB right: "<<chiSBsmear<<endl;
                cout<<"chi2/ndf SB right: "<<chisqTOTSB<<endl;







                //                float chiTOT=chisqTOTSB+cndf;
                                auto smearedfSBl_CF = GetSmearedCF(gausfit, histSmearTOTSBleft);

                               // smearedfSB_CF->Write(Form("smearedSB_%d", counter));
                                double IntegralCFleft =integrate(smearedfSBl_CF,800,900);
                                Double_t scaleleft = 1/IntegralCFleft;

                                auto normSBleft=GetSmearedCFscaled(gausfit, histSmearTOTSBleft, scaleleft);
                                normSBleft->SetLineColor(2);
                                normSBleft->SetMarkerColor(2);
                                normSBleft->SetMarkerSize(1.5);
                                normSBleft->SetMarkerStyle(21);
                                normSBleft->Write(Form("normsmearedSBleft_%d", counter));

                                int ndfsbleft=ndfSBleft;
                                float chiSBleftsmear= getchisquared(SBleftc,normSBleft);
                                float chisqTOTSBleft=chiSBleftsmear/ndfsbleft;

                                cout<<"chi2 SB left: "<<chiSBleftsmear<<endl;
                                cout<<"chi2/ndf SB left: "<<chisqTOTSBleft<<endl;



                                float chiTOT=(chiSBsmear+chi+chiSBleftsmear)/(ndf+ndfsbright+ndfsbleft);

                                cout<<"chiTOT/ndf: "<<chiTOT<< "   chi2 sum: "<< chiSBsmear+chi+chiSBleftsmear<<"   ndfsum: "<<ndf+ndfsbright+ndfsbleft<<endl;


                if(chiTOT<=minTOT)
                {
                    minTOT=chiTOT;
                    minTOTcounter=counter;


                }

                if(chisqTOTSB<=minTOTSB)
                {
                    minTOTSB=chisqTOTSB;
                    minTOTcounterSB=counter;
                }


                if(chisqTOTSBleft<=minTOTSBleft)
                {
                    minTOTSBleft=chisqTOTSBleft;
                    minTOTcounterSBleft=counter;
                }



             }



              if(cndf<=minimumchindf)
              {
                  minimumchindf=cndf;
                  mincounter=counter;
              }

//              if(chiSBsmear<=minTOT)
//              {
//                  minTOT=chiSBsmear;
//                  minTOTcounter=counter;
//              }

              counter++;
            }
          }
        }
      }
    }
  }

  cout<<"ndf FIT: "<< ndfFIT<<"  ndf SB: "<< ndfSB<< endl;

  cout<<"mini cndf FIT: "<< minimumchindf<<endl;
  cout<<"mini counter FIT: "<< mincounter<<endl;

  cout<<"mini cndf SB right: "<< minTOTSB<<endl;
  cout<<"mini counter SB right: "<< minTOTcounterSB<<endl;

  cout<<"mini cndf SB left: "<< minTOTSBleft<<endl;
  cout<<"mini counter SB left: "<< minTOTcounterSBleft<<endl;


  cout<<"mini cndf TOT: "<< minTOT<<endl;
  cout<<"mini counter TOT: "<< minTOTcounter<<endl;




  int toto=3*3*3*3*3*3;

  auto gcy = new TCanvas("pKCFall", "pKCFall", 0, 0, 650, 550);
  gcy->SetRightMargin(right);
  gcy->SetTopMargin(top);

  dummyHist->Draw();
  dummyHist->GetYaxis()->SetRangeUser(0.85, 1.2);
  pKCF666->GetYaxis()->SetTitle(" #it{C}(#it{k}*)");
  pKCF666->GetXaxis()->SetTitle(" #it{k}* (MeV/#it{c})");

  pKCF666->SetTitle("p-K function with gaus+pol2 fit");
  pKCF666->SetMarkerStyle(21);
  pKCF666->Draw("AP");
  pKCF666->GetYaxis()->SetRangeUser(0.9, 1.2);
  outfile->cd();
//  for (int i =0;i<toto;i++)
//  {
//      auto fitfun= (TF1*)outfile->Get(Form("chifit_%d", i));
//      if(fitfun){
//      fitfun->SetLineColorAlpha(kGreen, 0.2);
//      fitfun->Draw("same");
//      }
//  }
  auto orifit= (TF1*)outfile->Get("chifit_364");
  if(orifit){
  orifit->Draw("same");
  }
//  auto bestfit= (TF1*)outfile->Get("chifit_498");
  auto bestfit= (TF1*)outfile->Get("chifit_526");
//  auto bestfit= (TF1*)outfile->Get("chifit_290");

  if(bestfit){
  bestfit->SetLineColor(6);
  bestfit->Draw("same");
  }

  pKCF666->SetMarkerStyle(21);
  pKCF666->SetMarkerColor(4);

  pKCF666->Draw("Psame");


  gcy->Print("pKCFallhh.pdf");
  gcy->Print("pKCFallhh.png");



  auto gcy2 = new TCanvas("pKSBall", "pKSBall", 0, 0, 650, 550);
  gcy2->SetRightMargin(right);
  gcy2->SetTopMargin(top);

      dummyHist->Draw();
      dummyHist->GetYaxis()->SetRangeUser(0.85, 1.15);
      SBrightc->SetTitle(
          "SB right 1.027<M_{K^{+}K^{-}}<1.1");

  //    SBrightc->GetYaxis()->SetRangeUser(0.95, 1.15);
      SBrightc->GetYaxis()->SetRangeUser(0.95, 1.1);

      SBrightc->GetXaxis()->SetRangeUser(5, 1000);
      SBrightc->GetYaxis()->SetTitle(" #it{C}(#it{k}*)");
      SBrightc->GetXaxis()->SetTitle(" #it{k}* (MeV/#it{c})");
      SBrightc->Draw("");

  outfile->cd();
//  for (int i =0;i<toto;i++)
//  {
//      auto fitfun= (TGraph*)outfile->Get(Form("normsmearedSB_%d", i));
//      if(fitfun){
//      fitfun->SetLineColorAlpha(kGreen, 0.2);
//      fitfun->Draw("same");
//      }
//  }
  auto orifitsb= (TF1*)outfile->Get("normsmearedSB_364");
  if(orifitsb){
  orifitsb->Draw("same");
  orifitsb->SetMarkerStyle(21);
  orifitsb->SetMarkerColor(4);
  }

//  auto bestfitsb= (TF1*)outfile->Get("normsmearedSB_498");
  auto bestfitsb= (TF1*)outfile->Get("normsmearedSB_526");
//  auto bestfitsb= (TF1*)outfile->Get("normsmearedSB_290");

  if(bestfitsb){
  bestfitsb->SetLineColor(6);
  bestfitsb->Draw("same");

  }

  SBrightc->Draw("same");
  gcy2->Print("rightsballhh.pdf");
  gcy2->Print("rightsballhh.png");



  auto gcy3 = new TCanvas("pKSBleall", "pKSBleall", 0, 0, 650, 550);
  gcy3->SetRightMargin(right);
  gcy3->SetTopMargin(top);

      dummyHist->Draw();
      dummyHist->GetYaxis()->SetRangeUser(0.85, 1.15);
      SBleftc->SetTitle(
          "SB right 1.027<M_{K^{+}K^{-}}<1.1");

  //    SBrightc->GetYaxis()->SetRangeUser(0.95, 1.15);
      SBleftc->GetYaxis()->SetRangeUser(0.85, 1.15);

      SBleftc->GetXaxis()->SetRangeUser(5, 1000);
      SBleftc->GetYaxis()->SetTitle(" #it{C}(#it{k}*)");
      SBleftc->GetXaxis()->SetTitle(" #it{k}* (MeV/#it{c})");
      SBleftc->Draw("");

  outfile->cd();
//  for (int i =0;i<toto;i++)
//  {
//      auto fitfunle= (TGraph*)outfile->Get(Form("normsmearedSBleft_%d", i));
//      if(fitfunle){
//      fitfunle->SetLineColorAlpha(kGreen, 0.2);
//      fitfunle->Draw("same");
//      }
//  }
  auto orifitsble= (TF1*)outfile->Get("normsmearedSBleft_364");
  if(orifitsble){
  orifitsble->Draw("same");
  orifitsble->SetMarkerStyle(21);
  orifitsble->SetMarkerColor(4);
  }

//  auto bestfitsb= (TF1*)outfile->Get("normsmearedSB_498");
  auto bestfitsble= (TF1*)outfile->Get("normsmearedSBleft_526");
// auto bestfitsble= (TF1*)outfile->Get("normsmearedSBleft_290");

  if(bestfitsble){
  bestfitsble->SetLineColor(6);
  bestfitsble->Draw("same");

  }

  SBleftc->Draw("same");
  gcy3->Print("leftsballhh.pdf");
  gcy3->Print("leftsballhh.png");





  cout << "nana:  " << a << endl;

  TF1* f3 = new TF1("f3", "pol3", 0, 700);

  f3->SetParLimits(0, 0.0, 0.95);  // 0.002,0.005

  // f1->SetParLimits(1,0.001,0.005); //0.002,0.005
  // f1->SetParLimits(2,-0.00003,0.0);
  // f1->SetParLimits(3,-0.00000001,0.00000001);
  // f1->SetParLimits(4,0.00000000001,0.0000000001);
  // f1->SetParLimits(5,0.000000000000001,0.00000000001);
  // f1->SetParLimits(6,-1,1);

  pKCF2->Fit(f3, "F", "", 0, 700);
  f3->SetParLimits(0, 0.5, 0.95);
  pKCF2->Fit(f3, "F", "", 0, 700);

  f3->SetParLimits(0, 0.6, 0.95);

  f3->SetParLimits(1, 0.002, 0.003);
  f3->SetParLimits(2, -0.000007, -0.000003);
  f3->SetParLimits(3, 0.0, 0.000000007);

  pKCF2->Fit(f3, "F", "", 0, 910);  // 910

  TF1* f4 = new TF1("f4", "pol4", 0, 950);
  f4->SetParLimits(0, 0.6, 0.95);

  f4->SetParameter(0, f3->GetParameter(0));
  f4->SetParameter(1, f3->GetParameter(1));
  f4->SetParameter(2, f3->GetParameter(2));
  f4->SetParameter(3, f3->GetParameter(3));

  pKCF2->Fit(f4, "F", "", 0, 910);

  TF1* f5 = new TF1("f5", "pol5", 0, 780);
  f5->SetParLimits(0, 0.6, 0.95);

  f5->SetParameter(0, f4->GetParameter(0));
  f5->SetParameter(1, f4->GetParameter(1));
  f5->SetParameter(2, f4->GetParameter(2));
  f5->SetParameter(3, f4->GetParameter(3));
  f5->SetParameter(4, f4->GetParameter(4));

  pKCF2->Fit(f5, "F", "", 0, 780);

  //  float ppo0 = f5->GetParameter(0);
  //  float ppo1 = f5->GetParameter(1);
  //  float ppo2 = f5->GetParameter(2);
  //  float ppo3 = f5->GetParameter(3);
  //  float ppo4 = f5->GetParameter(4);
  //  float ppo5 = f5->GetParameter(5);

  //  float epo0 = f5->GetParError(0);
  //  float epo1 = f5->GetParError(1);
  //  float epo2 = f5->GetParError(2);
  //  float epo3 = f5->GetParError(3);
  //  float epo4 = f5->GetParError(4);
  //  float epo5 = f5->GetParError(5);

  //  // std::vector<float> param0={{pg0-eg0,pg0,pg0+eg0}};
  //  // std::vector<float> param1={{pg1-eg1,pg1,pg1+eg1}};
  //  // std::vector<float> param2={{pg2-eg2,pg2,pg2+eg2}};
  //  // std::vector<float> param3={{pg3-eg3,pg3,pg3+eg3}};
  //  // std::vector<float> param4={{pg4-eg4,pg4,pg4+eg4}};
  //  // std::vector<float> param5={{pg5-eg5,pg5,pg5+eg5}};

  //  float perpoer = 0.1;

  //  std::vector<float> param0po = {
  //      {ppo0 - (epo0 * perpoer), ppo0, ppo0 + (epo0 * perpoer)}};
  //  std::vector<float> param1po = {
  //      {ppo1 - (epo1 * perpoer), ppo1, ppo1 + (epo1 * perpoer)}};
  //  std::vector<float> param2po = {
  //      {ppo2 - (epo2 * perpoer), ppo2, ppo2 + (epo2 * perpoer)}};
  //  std::vector<float> param3po = {
  //      {ppo3 - (epo3 * perpoer), ppo3, ppo3 + (epo3 * perpoer)}};
  //  std::vector<float> param4po = {
  //      {ppo4 - (epo4 * perpoer), ppo4, ppo4 + (epo4 * perpoer)}};
  //  std::vector<float> param5po = {
  //      {ppo5 - (epo5 * perpoer), ppo5, ppo5 + (epo5 * perpoer)}};

  //  // for(int p0=0;p0<param0.size(); p0++){
  //  // std::cout<< "p1: "<< param0[p0]<<endl;
  //  //}
  //  // cout<<"para0 size: "<<param0.size()<<endl;

  //  auto gcpol = new TCanvas("pKCFpo", "pKCFpo", 0, 0, 650, 550);
  //  gcpol->SetRightMargin(right);
  //  gcpol->SetTopMargin(top);

  //  dummyHist->Draw();
  //  dummyHist->GetYaxis()->SetRangeUser(0.85, 1.2);
  //  pKCF2->SetTitle("p-K function with pol5 fit");
  //  pKCF2->SetMarkerStyle(21);
  //  pKCF2->Draw("ALP");

  //  int histnr2 = 730;
  //  TF1* fits2[histnr2];

  //  // for(int counter=0; counter<histnr; counter++){

  //  int counter2 = 0;
  //  outfile->cd();

  //  // histSmearTOTSBleft->Write("histSmearTOTSBleft");
  //  // for(int counter2=0; counter2<histnr; counter2++){
  //  // while(counter2<histnr){
  //  for (int p0 = 0; p0 < param0po.size(); p0++) {
  //    for (int p1 = 0; p1 < param1po.size(); p1++) {
  //      for (int p2 = 0; p2 < param2po.size(); p2++) {
  //        for (int p3 = 0; p3 < param3po.size(); p3++) {
  //          for (int p4 = 0; p4 < param4po.size(); p4++) {
  //            for (int p5 = 0; p5 < param5po.size(); p5++) {
  //              // while(counter2<histnr){
  //              cout << "counter2 nr: " << counter2 << endl;

  //              fits2[counter2] = nullptr;
  //              fits2[counter2] =
  //                  new TF1(Form("polfit_%d", counter2), "pol5", 0, 850);
  //              // std::cout<< "p1: "<< param0[p0]<<" p2: "<<
  //              param1[p1]<<endl;
  //              fits2[counter2]->SetParameter(0, param0po[p0]);
  //              fits2[counter2]->SetParameter(1, param1po[p1]);
  //              fits2[counter2]->SetParameter(2, param2po[p2]);
  //              fits2[counter2]->SetParameter(3, param3po[p3]);
  //              fits2[counter2]->SetParameter(4, param4po[p4]);
  //              fits2[counter2]->SetParameter(5, param5po[p5]);

  //              pKCF2->Fit(fits2[counter2]);
  //              fits2[counter2]->SetLineColorAlpha(kGreen, 0.2);
  //              fits2[counter2]->Draw("same");
  //              // pKCF->Draw("same");
  //              fits2[counter2]->Write(Form("polfit_%d", counter2));
  //              counter2++;
  //            }
  //          }
  //          //}
  //        }
  //      }
  //    }
  //  }

  //  gcpol->Print("pKCFall2.pdf");
  //  gcpol->Print("pKCFall2.png");

  // f3->SetParLimits(0,0.7,0.95);
  // pKCF2->Fit(f3,"F", "", 0, 910);
  // f3->SetParLimits(0,0.75,0.95);
  // pKCF2->Fit(f3,"F", "", 0, 910);

  // TF1 *f4 = new TF1("f4","pol5",0,700);
  // f3->SetParLimits(0,0.5,0.95);

//  const int binwidth = 1;
//  const int NumMomBins = 800;
//  float p1 = total->GetParameter(0);
//  float p2 = total->GetParameter(1);
//  float p3 = total->GetParameter(2);
//  float p4 = total->GetParameter(3);
//  float p5 = total->GetParameter(4);
//  float p6 = total->GetParameter(5);


  const int binwidth = 1;
  const int NumMomBins = 800;
  float p1 = f13->GetParameter(0);
  float p2 = f13->GetParameter(1);
  float p3 = f13->GetParameter(2);
  float p4 = f13->GetParameter(3);
  float p5 = f13->GetParameter(4);
  float p6 = f13->GetParameter(5);

  std::cout<<p1<< " "<< p2<< " "<< p3<< " "<< p4<< " "<< p5<< " "<< p6<<endl;
  // float p6=f1->GetParameter(5);

  float p1s = f5->GetParameter(0);
  float p2s = f5->GetParameter(1);
  float p3s = f5->GetParameter(2);
  float p4s = f5->GetParameter(3);
  float p5s = f5->GetParameter(4);
  float p6s = f5->GetParameter(5);

  // float dd=fit1->GetParameter(3);

  Int_t nn = NumMomBins;
  Double_t x2[nn], y2[nn];
  for (Int_t i = 0; i < nn; i++) {
    x2[i] = i * binwidth;
    y2[i] = p1s + p2s * x2[i] + p3s * x2[i] * x2[i] +
            p4s * x2[i] * x2[i] * x2[i] + p5s * x2[i] * x2[i] * x2[i] * x2[i] +
            p6s * x2[i] * x2[i] * x2[i] * x2[i] * x2[i];

    // cout<< x2[i]<< "    y= "<<y2[i]<<endl;
  }

  Double_t x3[nn], y3[nn];

  for (Int_t i = 0; i < nn; i++) {
    x3[i] = i * binwidth;
    y3[i] = p1 * TMath::Gaus(x3[i], p2, p3, false) + p4 + p5 * x3[i] +
            p6 * x3[i] * x3[i];

    // cout<< x2[i]<< "    y= "<<y2[i]<<endl;
  }

  TGraph* fitfunpol = new TGraph(nn, x2, y2);
  TGraph* fitfungaus = new TGraph(nn, x3, y3);

  auto smearedfitCF = GetSmearedCF(fitfungaus, histSmearTOTPeak);
  auto smearedfitSBleft = GetSmearedCF(fitfungaus, histSmearTOTSBleft);
  auto smearedfitSBright = GetSmearedCF(fitfungaus, histSmearTOTSBright);

  smearedfitCF->SetLineColor(2);
  // smearedCF->SetLineWidth(4);
  smearedfitCF->SetMarkerColor(2);
  smearedfitCF->SetMarkerSize(1.5);
  smearedfitCF->SetMarkerStyle(21);

  smearedfitSBleft->SetLineColor(2);
  // smearedCF->SetLineWidth(4);
  smearedfitSBleft->SetMarkerColor(2);
  smearedfitSBleft->SetMarkerSize(1.5);
  smearedfitSBleft->SetMarkerStyle(21);

  smearedfitSBright->SetLineColor(2);
  // smearedCF->SetLineWidth(4);
  smearedfitSBright->SetMarkerColor(2);
  smearedfitSBright->SetMarkerSize(1.5);
  smearedfitSBright->SetMarkerStyle(21);

  auto smearedfitCFpol = GetSmearedCF(fitfunpol, histSmearTOTPeak);
  auto smearedfitSBleftpol = GetSmearedCF(fitfunpol, histSmearTOTSBleft);
  auto smearedfitSBrightpol = GetSmearedCF(fitfunpol, histSmearTOTSBright);

  smearedfitCFpol->SetLineColor(2);
  // smearedCF->SetLineWidth(4);
  smearedfitCFpol->SetMarkerColor(2);
  smearedfitCFpol->SetMarkerSize(1.5);
  smearedfitCFpol->SetMarkerStyle(21);

  smearedfitSBleftpol->SetLineColor(2);
  // smearedCF->SetLineWidth(4);
  smearedfitSBleftpol->SetMarkerColor(2);
  smearedfitSBleftpol->SetMarkerSize(1.5);
  smearedfitSBleftpol->SetMarkerStyle(21);

  smearedfitSBrightpol->SetLineColor(2);
  // smearedCF->SetLineWidth(4);
  smearedfitSBrightpol->SetMarkerColor(2);
  smearedfitSBrightpol->SetMarkerSize(1.5);
  smearedfitSBrightpol->SetMarkerStyle(21);

  double max_x = getxmax(smearedfitCF);
  double max_y = getymax(smearedfitCF);

  double reX = std::abs(260 - max_x);
  double reY = std::abs(1.035 - max_y);
  cout << "rely : " << reY << endl;
  auto movedCF = GetSmearedCFmoved(fitfungaus, histSmearTOTPeak, reX, reY);
  movedCF->SetLineColor(2);
  // smearedCF->SetLineWidth(4);
  movedCF->SetMarkerColor(2);
  movedCF->SetMarkerSize(1.5);
  movedCF->SetMarkerStyle(21);

  double max_x2 = getxmax(smearedCF);
  double max_y2 = getymax(smearedCF);

  double reX2 = std::abs(260 - max_x2);
  double reY2 = std::abs(1.035 - max_y2);
  cout << "rely : " << reY << endl;
  auto movedCF2 = GetSmearedCFmoved(pKCF, histSmearTOTPeak, reX2, reY2);
  movedCF2->SetLineColor(2);
  // smearedCF->SetLineWidth(4);
  movedCF2->SetMarkerColor(2);
  movedCF2->SetMarkerSize(1.5);
  movedCF2->SetMarkerStyle(21);

  double max_xl = getxmax(smearedfitSBleft);
  double max_yl = getymax(smearedfitSBleft);

  double reXl = std::abs(340 - max_xl);
  double reYl = std::abs(1.071 - max_yl);
  auto movedleft =
      GetSmearedCFmoved(fitfungaus, histSmearTOTSBleft, reXl, reYl);
  movedleft->SetLineColor(2);
  // smearedCF->SetLineWidth(4);
  movedleft->SetMarkerColor(2);
  movedleft->SetMarkerSize(1.5);
  movedleft->SetMarkerStyle(21);

  double max_xl2 = getxmax(smearedSBleft);
  double max_yl2 = getymax(smearedSBleft);

  double reXl2 = std::abs(340 - max_xl2);
  double reYl2 = std::abs(1.071 - max_yl2);
  auto movedleft2 = GetSmearedCFmoved(pKCF, histSmearTOTSBleft, reXl2, reYl2);
  movedleft2->SetLineColor(2);
  // smearedCF->SetLineWidth(4);
  movedleft2->SetMarkerColor(2);
  movedleft2->SetMarkerSize(1.5);
  movedleft2->SetMarkerStyle(21);

  double max_xr = getxmax(smearedfitSBright);
  double max_yr = getymax(smearedfitSBright);

  double reXr = std::abs(300 - max_xr);
  double reYr = std::abs(1.056 - max_yr);
  auto movedright =
      GetSmearedCFmoved(fitfungaus, histSmearTOTSBright, reXr, reYr);
  movedright->SetLineColor(2);
  movedright->SetMarkerColor(2);
  movedright->SetMarkerSize(1.5);
  movedright->SetMarkerStyle(21);

  double max_xr2 = getxmax(smearedSBright);
  double max_yr2 = getymax(smearedSBright);

  double reXr2 = std::abs(300 - max_xr2);
  double reYr2 = std::abs(1.056 - max_yr2);
  auto movedright2 = GetSmearedCFmoved(pKCF, histSmearTOTSBright, reXr2, reYr2);
  movedright2->SetLineColor(2);
  movedright2->SetMarkerColor(2);
  movedright2->SetMarkerSize(1.5);
  movedright2->SetMarkerStyle(21);

  int lineWidth = 3;

  outfile->cd();

  histSmearTOTSBleft->Write("histSmearTOTSBleft");
  histSmearTOTPeak->Write("histSmearTOTPeak");
  histSmearTOTSB2762->Write("histSmearTOTSB2762");
  histSmearTOTSB6210->Write("histSmearTOTSB6210");
  histSmearTOTSB1015->Write("histSmearTOTSB1015");
  histSmearTOTSB1520->Write("histSmearTOTSB1520");
  histSmearTOTSB2025->Write("histSmearTOTSB2025");
  histSmearTOTSB2530->Write("histSmearTOTSB2530");
  histSmearTOTSB3035->Write("histSmearTOTSB3035");
  histSmearTOTSB3540->Write("histSmearTOTSB3540");
  histSmearTOTSB4045->Write("histSmearTOTSB4045");
  histSmearTOTSB4550->Write("histSmearTOTSB4550");
  histSmearTOTSB5055->Write("histSmearTOTSB5055");
  histSmearTOTSB5560->Write("histSmearTOTSB5560");
  histSmearTOTSB67->Write("histSmearTOTSB67");
  histSmearTOTSB78->Write("histSmearTOTSB78");
  histSmearTOTSB89->Write("histSmearTOTSB89");
  histSmearTOTSB90->Write("histSmearTOTSB90");
  histSmearTOTSB2223->Write("histSmearTOTSB2223");
  histSmearTOTSB2425->Write("histSmearTOTSB2425");
  histSmearTOTSB2627->Write("histSmearTOTSB2627");
  histSmearTOTSB2829->Write("histSmearTOTSB2829");
  histSmearTOTSB3031->Write("histSmearTOTSB3031");
  histSmearTOTSB3334->Write("histSmearTOTSB3334");
  histSmearTOTSB3738->Write("histSmearTOTSB3738");
  histSmearTOTSB4041->Write("histSmearTOTSB4041");

  SBleftc->Write("SBleftc");
  pPhiCFc->Write("Peakc");
  SBrightc->Write("SBrightc");
  SB12c->Write("SB12c");
  SB23c->Write("SB23c");
  SB34c->Write("SB34c");
  SB45c->Write("SB45c");
  SB56c->Write("SB56c");
  SB67c->Write("SB67c");
  SB78c->Write("SB78c");
  SB89c->Write("SB89c");
  SB90c->Write("SB90c");

  SBleft->Write("SBleft");
  pPhiCF->Write("Peak");
  SB2762->Write("SB2762");
  SB6210->Write("SB6210");
  SB1015->Write("SB1015");
  SB1520->Write("SB1520");
  SB2025->Write("SB2025");
  SB2530->Write("SB2530");
  SB3035->Write("SB3035");
  SB3540->Write("SB3540");
  SB4045->Write("SB4045");
  SB4550->Write("SB4550");
  SB5055->Write("SB5055");
  SB5560->Write("SB5560");
  SB67->Write("SB67");
  SB78->Write("SB78");
  SB89->Write("SB89");
  SB90->Write("SB90");
  SB2223->Write("SB2223");
  SB2425->Write("SB2425");
  SB2627->Write("SB2627");
  SB2829->Write("SB2829");
  SB3031->Write("SB3031");
  SB3334->Write("SB3334");
  SB3738->Write("SB3738");
  SB4041->Write("SB4041");

  pKCF->Write("pKCFgaus");
  pKBG->Write("pKBG");
  pKCF2->Write("pKCFpol");

  //  smearedCF->Write("smearedCF");
  //  smearedSBleft->Write("smearedSBleft");
  //  smearedSBright->Write("smearedSBright");

  //  fitfunpol->Write("fitfunpol");
  //  fitfungaus->Write("fitfungaus");

  //  smearedfitCF->Write("smearedfitCF");
  //  smearedfitSBleft->Write("smearedfitSBleft");
  //  smearedfitSBright->Write("smearedfitSBright");

  //  smearedfitCFpol->Write("smearedfitCFpol");
  //  smearedfitSBleftpol->Write("smearedfitSBleftpol");
  //  smearedfitSBrightpol->Write("smearedfitSBrightpol");

  //  movedCF->Write("movedCF");
  //  movedleft->Write("movedleft");
  //  movedright->Write("movedright");

  //  movedCF2->Write("movedCF2");
  //  movedleft2->Write("movedleft2");
  //  movedright2->Write("movedright2");

  outfile->Close();
}

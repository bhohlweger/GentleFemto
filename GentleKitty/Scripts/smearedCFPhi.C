#include "TidyCats.h"
#include "TRandom3.h"
#include "SidebandSigma.h"
#include "CATSInputSigma0.h"
#include "CATSLambdaParam.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "DLM_Source.h"
#include "DLM_Potentials.h"
#include "DLM_CkModels.h"
#include "CATS.h"
#include "DLM_CkDecomposition.h"
#include "TDatabasePDG.h"
#include <iostream>
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TNtuple.h"
#include "TSystem.h"
#include "TRandom3.h"
#include "TGenPhaseSpace.h"
#include "TStopwatch.h"
#include "DreamPlot.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TVirtualFitter.h"
#include "TGraph.h"
#include "TGraphErrors.h"

/// =====================================================================================
void FillWaveGraph(CATS &kitty, TGraph *gr) {
  for (unsigned int i = 0; i < kitty.GetNumMomBins(); ++i) {
    const float mom = kitty.GetMomentum(i);
    gr->SetPoint(i, mom, kitty.GetCorrFun(i));
  }
}

/// =====================================================================================
TGraph *GetSmearedCF(TGraph* CF, TH2* matrix) {
  //Define new Histogram which have dimension according to the yaxis (new momentum axis):
  const int nbins_original = matrix->GetXaxis()->GetNbins();
  const Int_t nbins_transformed = matrix->GetYaxis()->GetNbins();

  TGraph *smearedCF = new TGraph();

  int countPoint = 0;
  for (int momTrans = 0; momTrans < nbins_transformed; momTrans++) {
    Double_t matrixvalues_sum = 0.;
    Double_t weighted_matrixvalues_sum = 0.;

    for (int momOri = 0; momOri < nbins_original; momOri++) {
      Double_t momentum_original = matrix->GetXaxis()->GetBinCenter(momOri + 1);
      matrixvalues_sum += matrix->GetBinContent(momOri + 1, momTrans + 1);
      weighted_matrixvalues_sum += CF->Eval(momentum_original)
          * matrix->GetBinContent(momOri + 1, momTrans + 1);
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


TGraph *GetSmearedCFscaled(TGraph* CF, TH2* matrix, double scale) {
  //Define new Histogram which have dimension according to the yaxis (new momentum axis):
  const int nbins_original = matrix->GetXaxis()->GetNbins();
  const Int_t nbins_transformed = matrix->GetYaxis()->GetNbins();

  TGraph *smearedCF = new TGraph();

  int countPoint = 0;
  for (int momTrans = 0; momTrans < nbins_transformed; momTrans++) {
    Double_t matrixvalues_sum = 0.;
    Double_t weighted_matrixvalues_sum = 0.;

    for (int momOri = 0; momOri < nbins_original; momOri++) {
      Double_t momentum_original = matrix->GetXaxis()->GetBinCenter(momOri + 1);
      matrixvalues_sum += matrix->GetBinContent(momOri + 1, momTrans + 1);
      weighted_matrixvalues_sum += CF->Eval(momentum_original)
          * matrix->GetBinContent(momOri + 1, momTrans + 1);
    }
    Double_t transformed_CF = 0.;
    if (matrixvalues_sum != 0.)
      transformed_CF = weighted_matrixvalues_sum / matrixvalues_sum;

    smearedCF->SetPoint(countPoint++,
                        matrix->GetYaxis()->GetBinCenter(momTrans + 1),
                        transformed_CF*scale);
  }
  return smearedCF;
}


TGraph *GetSmearedCF2(TGraph* CF, TH2* matrix) {
  //Define new Histogram which have dimension according to the yaxis (new momentum axis):
  const int nbins_original = matrix->GetXaxis()->GetNbins();
  const Int_t nbins_transformed = matrix->GetYaxis()->GetNbins();

  TGraph *smearedCF = new TGraph();

  int countPoint = 0;
  for (int momTrans = 0; momTrans < nbins_transformed; momTrans++) {
    Double_t matrixvalues_sum = 0.;
    Double_t weighted_matrixvalues_sum = 0.;

    for (int momOri = 0; momOri < nbins_original; momOri++) {
      Double_t momentum_original = matrix->GetXaxis()->GetBinCenter(momOri + 1);
      matrixvalues_sum += matrix->GetBinContent(momOri + 1, momTrans + 1);
      weighted_matrixvalues_sum += CF->Eval(momentum_original)
          * matrix->GetBinContent(momOri + 1, momTrans + 1);
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


TGraph *GetSmearedCFmoved(TGraph* CF, TH2* matrix, double max_x, double max_y) {
  //Define new Histogram which have dimension according to the yaxis (new momentum axis):

//  double max_cf_x=CF->GetMaximum();

  const int nbins_original = matrix->GetXaxis()->GetNbins();
  const Int_t nbins_transformed = matrix->GetYaxis()->GetNbins();

  TGraph *smearedCF = new TGraph();

  int countPoint = 0;
  for (int momTrans = 0; momTrans < nbins_transformed; momTrans++) {
    Double_t matrixvalues_sum = 0.;
    Double_t weighted_matrixvalues_sum = 0.;

    for (int momOri = 0; momOri < nbins_original; momOri++) {
      Double_t momentum_original = matrix->GetXaxis()->GetBinCenter(momOri + 1);
      matrixvalues_sum += matrix->GetBinContent(momOri + 1, momTrans + 1);
      weighted_matrixvalues_sum += CF->Eval(momentum_original)
          * matrix->GetBinContent(momOri + 1, momTrans + 1);
    }
    Double_t transformed_CF = 0.;
    if (matrixvalues_sum != 0.)
      transformed_CF = weighted_matrixvalues_sum / matrixvalues_sum;

    smearedCF->SetPoint(countPoint++,
                        (matrix->GetYaxis()->GetBinCenter(momTrans + 1)-max_x),
                        (transformed_CF-max_y));

  }
  return smearedCF;
}


double getxmax(TGraph* CF){
double ymaxcf =0;
double xmaxcf=0;
for (int i=0;i<1000;i++){
double ycf= CF->Eval(i);
//cout<<i<< " yval: "<<ycf<<endl;
if (ycf>=ymaxcf){
    ymaxcf=ycf;
    xmaxcf=i;
}
}
//cout <<xmaxcf<<endl;
return xmaxcf;

}

double integrate(TGraph* CF, double min, double max){
double sum=0;
int diff=max-min;
int a=0;
for (int i=min;i<=max;i++){
double ycf= CF->Eval(i);
//cout<<"value: "<< i<< " ----> "<<ycf<<endl;
sum = sum+ycf;
a++;
}
//cout<<sum<<endl;
return sum/a;
//cout<<sum/a<<endl;

}

double getymax(TGraph* CF){
double ymaxcf =0;
double xmaxcf=0;
for (int i=0;i<1000;i++){
double ycf= CF->Eval(i);
//cout<<i<< " yval: "<<ycf<<endl;
if (ycf>=ymaxcf){
    ymaxcf=ycf;
    xmaxcf=i;
}
}
//cout <<ymaxcf<<endl;

return ymaxcf;
}



//double getpointtoval(TGraph* CF, double value){
//for (int i=0;i<CF->GetN();i++){
//double ycf= CF->Eval(i);
////cout<<i<< " yval: "<<ycf<<endl;
//if (ycf>=ymaxcf){
//    ymaxcf=ycf;
//    xmaxcf=i;
//}
//}
//cout <<ymaxcf<<endl;

//return ymaxcf;
//}

//void setmaxima(TGraph* CF,double ymaxcf, double xmaxcf ){
//ymaxcf =0;
//xmaxcf=0;
//for (int i=0;i<CF->GetN();i++){
//double ycf= CF->Eval(i);
//cout<<i<< " yval: "<<ycf<<endl;
//while (ycf>=ymaxcf){
//    ymaxcf=ycf;
//    xmaxcf=i;
//}
//}
//}


TH1F* CalcCF( const char* prefix,  int a) {
  const char* addon=Form("%d",a);
//  const char* name=Form("/home/emma/FemtoPhiHM_checks2/1000-1200/CFOutput_pPhi_%s_%s.root",prefix, addon);
//  const char* name=Form("/home/emma/FemtoPhiHM_checks2/800900/CFOutput_pPhi_%s_%s.root",prefix, addon);
  const char* name=Form("/home/emma/FemtoPhiHM_ROTMC/phinormal10001200/CFOutput_pPhi_%s_%s.root",prefix, addon);

  auto file = TFile::Open(name);
 // file->ls();
  TH1F* corr = (TH1F*)file->FindObjectAny("hCk_ReweightedMeV_2");
  return corr;

}

TH1F* CalcCF3( const char* prefix,  int a) {
  const char* addon=Form("%d",a);
  const char* name=Form("/home/emma/FemtoPhiHM_checks2/1000-1200/CFOutput_pPhi_%s_%s.root",prefix, addon);
//  const char* name=Form("/home/emma/FemtoPhiHM_checks2/800900/CFOutput_pPhi_%s_%s.root",prefix, addon);
//  const char* name=Form("/home/emma/FemtoPhiHM_checks2/1000-1200/CFOutput_pPhi_%s_%s.root",prefix, addon);

  auto file = TFile::Open(name);
 // file->ls();
  TH1F* corr = (TH1F*)file->FindObjectAny("hCk_ReweightedMeV_2");
  return corr;

}

TH1F* CalcCF2( const char* prefix,  int a) {
  const char* addon=Form("%d",a);
//  const char* name=Form("/home/emma/FemtoPhiHM_allSB/1000-1200/CFOutput_pPhi_%s_%s.root",prefix, addon);
//  const char* name=Form("/home/emma/FemtoPhiHM_allSB/800900/CFOutput_pPhi_%s_%s.root",prefix, addon);
  const char* name=Form("/home/emma/FemtoPhiHM_allSBMC/CFOutput_pPhi_%s_%s.root",prefix, addon);


  auto file = TFile::Open(name);
  //file->ls();
  TH1F* corr = (TH1F*)file->FindObjectAny("hCk_ReweightedMeV_2");
  return corr;

}

//TH2F* GetMatrix( const char* prefix,  int a) {
//  const char* addon=Form("%d",a);
//  const char* name=Form("/home/emma/PhiTRUTH_2/phitruth/CFOutput_pPhi_%s_%s.root",prefix, addon);
//  auto file = TFile::Open(name);
//  //file->ls();
//  TH1F* corr = (TH1F*)file->FindObjectAny("hCk_ReweightedMeV_2");
//  return corr;

//}


double evalMomOri(double kstar, TH2* matrix) {
  //Define new Histogram which have dimension according to the yaxis (new momentum axis):
  const int nbins_original = matrix->GetYaxis()->GetNbins();

  int binori=0;
  for (int momOri = 0; momOri < nbins_original; momOri++) {
    Double_t momentum_original = matrix->GetYaxis()->GetBinCenter(momOri + 1);
    if (momentum_original>=kstar)
    {
        binori=(momOri+1);
        break;
    }
  }
  return binori;

}

double evalMomTrafo(double bin, TH2* matrix) {
  //Define new Histogram which have dimension according to the yaxis (new momentum axis):
  const Int_t nbins_transformed = matrix->GetXaxis()->GetNbins();
  double max=0;
  double maxbin=0;
  for (int momTrans = 0; momTrans < nbins_transformed; momTrans++) {
    double matrixvalues = 0.;
    matrixvalues = matrix->GetBinContent( momTrans + 1, bin);
    //cout<<"matrixval"<<momTrans + 1<< ": "<<matrixvalues<<endl;
    if (matrixvalues>=max){
        max=matrixvalues;
        maxbin=momTrans+1;
    }
    //cout<<max<<endl;
}
  Double_t momentum_transformed = matrix->GetXaxis()->GetBinCenter(maxbin + 1);

   return momentum_transformed;
}

double GetWidth(double bin, TH2* matrix){
  const Int_t nbins_transformed = matrix->GetXaxis()->GetNbins();
  double bin_width=matrix->GetXaxis()->GetBinWidth(1);
  int count = 0;
  for (int momTrans = 0; momTrans < nbins_transformed; momTrans++) {
    Double_t matrixvalues = 0.;
    matrixvalues = matrix->GetBinContent(momTrans + 1,bin);
    if (matrixvalues!=0){
        count++;
    }
    }
  double width=count*bin_width;
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
  //cout<<"wifht: "<< width<<" maxbin: "<<maxbin<<"maxvaue: "<<max<<"mean: "<< mean<<endl;

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


double PolN_flat_end(double* x, double* par){
     double& lim = par[1];
     unsigned Order = par[0];
     double* pol = &par[2];
     if(lim!=-1e6){
         pol[1] = 0;
         for(unsigned uOrd=2; uOrd<=Order; uOrd++){
             pol[1] -= uOrd*pol[uOrd]*(pow(lim,uOrd-1));
         }
     }
     double Result=0;
     for(unsigned uOrd=0; uOrd<=Order; uOrd++){
         Result += pol[uOrd]*pow(*x,uOrd);
     }
     return Result;
}

//void myFcn(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */  )
//{
//  int n = coords.size();
//  double chi2 = 0;
//  double tmp,x[2];
//  for (int i = 0; i <n; ++i ) {
//    x[0] = coords[i].first;
//    x[1] = coords[i].second;
//    tmp = ( values[i] - my2Dfunc(x,p))/errors[i];
//    chi2 += tmp*tmp;
//  }
//  fval = chi2;
//}


float getchisquared(TH1F *hist, TGraph *func)
{
    int nbins=hist->GetNbinsX();
    double chi2=0;
    double tmp=0;
    cout<< "nbins="<< nbins<<endl;
    for (int i = 0; i <nbins; ++i ) {
    double momentum= hist->GetBinCenter(i + 1);
    double funval= func->Eval(momentum);
    double histval= hist->GetBinContent(i + 1);
    double error=hist->GetBinError(i + 1);
    tmp=(histval-funval)/error;
    chi2+=tmp*tmp;
    cout<< "funval: "<<funval<<"   histval: "<< histval<< "   error: "<< error<<endl;
    }
    return chi2;
    cout<< "chi2: ..."<< chi2<< endl;
}

//void SetCFPoint(Int_t i, TH1F* CF, TH2F* matrix, int binnr, TGraphErrors Graph ){
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
int main(int argc, char *argv[]) {
  gROOT->ProcessLine("gErrorIgnoreLevel = 3001");
  DreamPlot::SetStyle();

  TString InputDir = argv[1];
  const char* trigger = argv[2];
  int suffix= atoi(argv[3]);

  auto filename = TString::Format("%s/SidebandInspect.root", InputDir.Data());
  auto outfile = new TFile(filename, "RECREATE");


  auto filenameSmear = TString::Format("%s/SmearSideband.root", InputDir.Data());
  auto infile = TFile::Open(filenameSmear);
  if (!infile) {
    std::cout
        << "No smearing matrix found - start the sideband computation task!\n";
    return 0;
  }

// Get Matrices

  std::cout<<"get matrices"<<std::endl;
  auto histSmearTOTSBleft = (TH2F*) infile->Get("histSmearTOTSBleft");
  auto histSmearTOTSBright = (TH2F*) infile->Get("histSmearTOTSBright");
  auto histSmearTOTSB12 = (TH2F*) infile->Get("histSmearTOTSB12");
  auto histSmearTOTSB23 = (TH2F*) infile->Get("histSmearTOTSB23");
  auto histSmearTOTSB34 = (TH2F*) infile->Get("histSmearTOTSB34");
  auto histSmearTOTSB45 = (TH2F*) infile->Get("histSmearTOTSB45");
  auto histSmearTOTSB56 = (TH2F*) infile->Get("histSmearTOTSB56");
  auto histSmearTOTSB67 = (TH2F*) infile->Get("histSmearTOTSB67");
  auto histSmearTOTSB78 = (TH2F*) infile->Get("histSmearTOTSB78");
  auto histSmearTOTSB89 = (TH2F*) infile->Get("histSmearTOTSB89");
  auto histSmearTOTSB90 = (TH2F*) infile->Get("histSmearTOTSB90");
  auto histSmearTOTPeak = (TH2F*) infile->Get("histSmearTOTPeak");

  auto histSmearTOTSB2762 = (TH2F*) infile->Get("histSmearTOT2762");
  auto histSmearTOTSB6210 = (TH2F*) infile->Get("histSmearTOT6210");
  auto histSmearTOTSB1015 = (TH2F*) infile->Get("histSmearTOT1015");
  auto histSmearTOTSB1520 = (TH2F*) infile->Get("histSmearTOT1520");
  auto histSmearTOTSB2025 = (TH2F*) infile->Get("histSmearTOT2025");
  auto histSmearTOTSB2530 = (TH2F*) infile->Get("histSmearTOT2530");
  auto histSmearTOTSB3035 = (TH2F*) infile->Get("histSmearTOT3035");
  auto histSmearTOTSB3540 = (TH2F*) infile->Get("histSmearTOT3540");
  auto histSmearTOTSB4045 = (TH2F*) infile->Get("histSmearTOT4045");
  auto histSmearTOTSB4550 = (TH2F*) infile->Get("histSmearTOT4550");
  auto histSmearTOTSB5055 = (TH2F*) infile->Get("histSmearTOT5055");
  auto histSmearTOTSB5560 = (TH2F*) infile->Get("histSmearTOT5560");

  auto histSmearTOTSB2223 = (TH2F*) infile->Get("histSmearTOT2223");
  auto histSmearTOTSB2425 = (TH2F*) infile->Get("histSmearTOT2425");
  auto histSmearTOTSB2627 = (TH2F*) infile->Get("histSmearTOT2627");
  auto histSmearTOTSB2829 = (TH2F*) infile->Get("histSmearTOT2829");
  auto histSmearTOTSB3031 = (TH2F*) infile->Get("histSmearTOT3031");
  auto histSmearTOTSB3334 = (TH2F*) infile->Get("histSmearTOT3334");
  auto histSmearTOTSB3738 = (TH2F*) infile->Get("histSmearTOT3738");
  auto histSmearTOTSB4041 = (TH2F*) infile->Get("histSmearTOT4041");


  const float right = 0.04;
  const float top = 0.025;

  auto dummyHist = new TH1F("dummyHist",
                            ";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)", 100, 4,
                            1000);
  DreamPlot::SetStyleHisto(dummyHist, 20, kGreen + 2);


  std::cout<<"setstyle"<<std::endl;
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

  std::cout<<"end style"<<std::endl;


  auto m1 = new TCanvas("m1", "m1", 650, 550);
  histSmearTOTPeak->Draw("col");
  histSmearTOTPeak->GetXaxis()->SetRangeUser(0, 500);
  histSmearTOTPeak->GetYaxis()->SetRangeUser(0, 500);
  histSmearTOTPeak->GetXaxis()->SetNdivisions(505);
  histSmearTOTPeak->GetYaxis()->SetNdivisions(505);
  histSmearTOTPeak->SetTitle(
      "; #it{k}*_{p#minus K^{+}} (MeV/#it{c}); #it{k}*_{p#minus#phi} (MeV/#it{c})");
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
      "; #it{k}*_{p#minus K^{+}} (MeV/#it{c}); #it{k}*_{p#minus (K^{+}K^{-})} (MeV/#it{c})");
  m2->Print("histSmearTOTSBleft.pdf");
  m2->Print("histSmearTOTSBleft.png");
delete m2;


 TH1F* SBleftc=CalcCF("HMPhi",1);
 TH1F* SBrightc=CalcCF("HMPhi",2);
 TH1F* SB12c=CalcCF3("HMPhi",3);
 TH1F* SB23c=CalcCF3("HMPhi",4);
 TH1F* SB34c=CalcCF3("HMPhi",5);
 TH1F* SB45c=CalcCF3("HMPhi",6);
 TH1F* SB56c=CalcCF3("HMPhi",7);
 TH1F* SB67c=CalcCF3("HMPhi",8);
 TH1F* SB78c=CalcCF3("HMPhi",9);
 TH1F* SB89c=CalcCF3("HMPhi",10);
 TH1F* SB90c=CalcCF3("HMPhi",11);
 TH1F* pPhiCFc=CalcCF("HMPhi",0);




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


// TH1F* SBleft=CalcCF("HMPhi",1);
// TH1F* SBright=CalcCF("HMPhi",2);
// TH1F* SB12=CalcCF("HMPhi",3);
// TH1F* SB23=CalcCF("HMPhi",4);
// TH1F* SB34=CalcCF("HMPhi",5);
// TH1F* SB45=CalcCF("HMPhi",6);
// TH1F* SB56=CalcCF("HMPhi",7);
// TH1F* SB67=CalcCF("HMPhi",8);
// TH1F* SB78=CalcCF("HMPhi",9);
// TH1F* SB89=CalcCF("HMPhi",10);
// TH1F* SB90=CalcCF("HMPhi",11);
// TH1F* pPhiCF=CalcCF("HMPhi",0);

 TH1F* pPhiCF=CalcCF2("HMPhi",0);
 TH1F* SBleft=CalcCF2("HMPhi",1);
 TH1F* SB2762=CalcCF2("HMPhi",2);
 TH1F* SB6210=CalcCF2("HMPhi",3);
 TH1F* SB1015=CalcCF2("HMPhi",4);
 TH1F* SB1520=CalcCF2("HMPhi",5);
 TH1F* SB2025=CalcCF2("HMPhi",6);
 TH1F* SB2530=CalcCF2("HMPhi",7);
 TH1F* SB3035=CalcCF2("HMPhi",8);
 TH1F* SB3540=CalcCF2("HMPhi",9);
 TH1F* SB4045=CalcCF2("HMPhi",10);
 TH1F* SB4550=CalcCF2("HMPhi",11);
 TH1F* SB5055=CalcCF2("HMPhi",12);
 TH1F* SB5560=CalcCF2("HMPhi",13);
 TH1F* SB67=CalcCF2("HMPhi",14);
 TH1F* SB78=CalcCF2("HMPhi",15);
 TH1F* SB89=CalcCF2("HMPhi",16);
 TH1F* SB90=CalcCF2("HMPhi",17);
 TH1F* SB2223=CalcCF2("HMPhi",18);
 TH1F* SB2425=CalcCF2("HMPhi",19);
 TH1F* SB2627=CalcCF2("HMPhi",20);
 TH1F* SB2829=CalcCF2("HMPhi",21);
 TH1F* SB3031=CalcCF2("HMPhi",22);
 TH1F* SB3334=CalcCF2("HMPhi",23);
 TH1F* SB3738=CalcCF2("HMPhi",24);
 TH1F* SB4041=CalcCF2("HMPhi",25);



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


//  void SetCFPoint(Int_t i, TH1F* CF, TH2F* matrix, int binnr, TGraphErrors Graph ){
//      double kstar= CF->GetBinCenter(binnr);
//      double ey= CF->GetBinError(binnr);
//      double val= CF->GetBinContent(binnr);

//      int biny=evalMomOri(kstar,matrix);
//      double max=evalMomTrafo(biny,matrix);
//      double ex=GetWidth(biny,matrix);

//      Graph.SetPoint(i,max,val);
//      Graph.SetPointError(i, ex, ey);

//  }

//auto pKCF = new TGraphErrors();

//double k1= SBleft->GetBinCenter(1);
//double ey1= SBleft->GetBinError(1);
//double C1= SBleft->GetBinContent(1);

//int biny1=evalMomOri(k1,histSmearTOTSBleft);
//double max1=evalMomTrafo(biny1,histSmearTOTSBleft);
//double ex1=GetWidth(biny1,histSmearTOTSBleft);
//pKCF.SetPoint(1,max1,C1);
//pKCF.SetPointError(1, ex1, ey1);


const Int_t n = 25;
Double_t x[n]  = {0};
Double_t y[n]  = {0};
Double_t ex[n] = {0};
Double_t ey[n] = {0};
Double_t exl[n] = {0};
Double_t exh[n] = {0};


const Int_t n1= 9;
Double_t xx[n1]  = {0};
Double_t yy[n1]  = {0};
Double_t exx[n1] = {0};
Double_t exx2[n1] = {0};
Double_t eyy[n1] = {0};

//y[1]=SBleft->GetBinContent(1);

//std::cout<<"cfvalue "<<y[1]<<endl;

//ey[1]=SBleft->GetBinError(1);

//std::cout<<"error "<<ey[1]<<endl;

//double k1= SBleft->GetBinCenter(1);

//std::cout<<"kstarcf "<<k1<<endl;

//int biny1=evalMomOri(k1,histSmearTOTSBleft);

//std::cout<<"biny  "<<biny1<<endl;

//double max1=evalMomTrafo(biny1,histSmearTOTSBleft);

//std::cout<<"valuemax k* bin  "<<max1<<endl;

//double ex1=GetWidth(biny1,histSmearTOTSBleft);

//std::cout<<"errorx  "<<ex1<<endl;

//x[1]=max1;
//ex[1]=ex1;





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
exx[0]=exl[0];
exx2[0]=exh[0];






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
exx[1]=exl[17];
exx2[1]=exh[17];


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

//yy[2]=y[18];
//xx[2]=x[18];
//eyy[2]=ey[18];
//exx[2]=ex[18];

yy[2]=y[18];
xx[2]=x[18];
eyy[2]=ey[18];
exx[2]=exl[18];
exx2[2]=exh[18];


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

//yy[3]=y[19];
//xx[3]=x[19];
//eyy[3]=ey[19];
//exx[3]=ex[19];

yy[3]=y[19];
xx[3]=x[19];
eyy[3]=ey[19];
exx[3]=exl[19];
exx2[3]=exh[19];


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
exx[4]=exl[20];
exx2[4]=exh[20];


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
exx[5]=exl[21];
exx2[5]=exh[21];

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
exx[6]=exl[22];
exx2[6]=exh[22];

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
exx[7]=exl[23];
exx2[7]=exh[23];

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
exx[8]=exl[24];
exx2[8]=exh[24];

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
pKCF->SetLineWidth(lineWidth2);


//auto pKCF = new TGraphErrors(n,x,y,ex,ey);
//int lineWidth2 = 3;
//DreamPlot::SetStyleGraph(pKCF, 20, kBlack, 0.7);
//pKCF->SetLineWidth(lineWidth2);


auto pKCF2 = new TGraphAsymmErrors(n,x,y,exl, exh ,ey,ey);
DreamPlot::SetStyleGraph(pKCF2, 20, kBlack, 0.7);
pKCF2->SetLineWidth(lineWidth2);

auto pKCFpol5 = new TGraphAsymmErrors(n,x,y,exl, exh ,ey,ey);
DreamPlot::SetStyleGraph(pKCFpol5, 20, kBlack, 0.7);
pKCFpol5->SetLineWidth(lineWidth2);


auto pKCFDimi = new TGraphAsymmErrors(n,x,y,exl, exh ,ey,ey);
DreamPlot::SetStyleGraph(pKCFDimi, 20, kBlack, 0.7);
pKCFDimi->SetLineWidth(lineWidth2);

auto pKCFDimipol5 = new TGraphAsymmErrors(n,x,y,exl, exh ,ey,ey);
DreamPlot::SetStyleGraph(pKCFDimipol5, 20, kBlack, 0.7);
pKCFDimipol5->SetLineWidth(lineWidth2);



auto pKBG = new TGraphAsymmErrors(n1,xx,yy,exx, exx2, eyy,eyy);
DreamPlot::SetStyleGraph(pKBG, 20, kBlack, 0.7);
pKBG->SetLineWidth(lineWidth2);


auto pKCF66 = new TGraphAsymmErrors(n,x,y,exl, exh,ey, ey);
DreamPlot::SetStyleGraph(pKCF66, 20, kBlack, 0.7);
pKCF66->SetLineWidth(lineWidth2);



//fitfun->Draw("L3same");

for (int a=0;a<11;a++){
    cout<< x[a]<<" "<< ex[a]<<" "<< y[a]<< " "<< ey[a]<<endl;
}

//SetCFPoint(int i, TH1F* CF, TH2F* matrix, int binnr, TGraphErrors Graph )


pKCF->SetLineColor(4);
//pKCF->SetLineWidth(4);
pKCF->SetMarkerColor(4);
pKCF->SetMarkerSize(1.5);
pKCF->SetMarkerStyle(21);

pKCF66->SetLineColor(4);
//pKCF->SetLineWidth(4);
pKCF66->SetMarkerColor(4);
pKCF66->SetMarkerSize(1.5);
pKCF66->SetMarkerStyle(21);


pKCF2->SetLineColor(4);
//pKCF->SetLineWidth(4);
pKCF2->SetMarkerColor(4);
pKCF2->SetMarkerSize(1.5);
pKCF2->SetMarkerStyle(21);


pKCFDimi->SetLineColor(4);
//pKCF->SetLineWidth(4);
pKCFDimi->SetMarkerColor(4);
pKCFDimi->SetMarkerSize(1.5);
pKCFDimi->SetMarkerStyle(21);


pKCFDimipol5->SetLineColor(4);
//pKCF->SetLineWidth(4);
pKCFDimipol5->SetMarkerColor(4);
pKCFDimipol5->SetMarkerSize(1.5);
pKCFDimipol5->SetMarkerStyle(21);

pKBG->SetLineColor(4);
//pKCF->SetLineWidth(4);
pKBG->SetMarkerColor(4);
pKBG->SetMarkerSize(1.5);
pKBG->SetMarkerStyle(21);

auto smearedCF=GetSmearedCF(pKCF,histSmearTOTPeak);
smearedCF->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedCF->SetMarkerColor(2);
smearedCF->SetMarkerSize(1.5);
smearedCF->SetMarkerStyle(21);

auto smearedSBleft=GetSmearedCF(pKCF,histSmearTOTSBleft);
smearedSBleft->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedSBleft->SetMarkerColor(2);
smearedSBleft->SetMarkerSize(1.5);
smearedSBleft->SetMarkerStyle(21);

auto smearedSBright=GetSmearedCF(pKCF,histSmearTOTSBright);
smearedSBright->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedSBright->SetMarkerColor(2);
smearedSBright->SetMarkerSize(1.5);
smearedSBright->SetMarkerStyle(21);

auto smearedSB12=GetSmearedCF(pKCF,histSmearTOTSB12);
smearedSB12->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedSB12->SetMarkerColor(2);
smearedSB12->SetMarkerSize(1.5);
smearedSB12->SetMarkerStyle(21);

auto smearedSB23=GetSmearedCF(pKCF,histSmearTOTSB23);
smearedSB23->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedSB23->SetMarkerColor(2);
smearedSB23->SetMarkerSize(1.5);
smearedSB23->SetMarkerStyle(21);

auto smearedSB34=GetSmearedCF(pKCF,histSmearTOTSB34);
smearedSB34->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedSB34->SetMarkerColor(2);
smearedSB34->SetMarkerSize(1.5);
smearedSB34->SetMarkerStyle(21);

auto smearedSB45=GetSmearedCF(pKCF,histSmearTOTSB45);
smearedSB45->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedSB45->SetMarkerColor(2);
smearedSB45->SetMarkerSize(1.5);
smearedSB45->SetMarkerStyle(21);

auto smearedSB56=GetSmearedCF(pKCF,histSmearTOTSB56);
smearedSB56->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedSB56->SetMarkerColor(2);
smearedSB56->SetMarkerSize(1.5);
smearedSB56->SetMarkerStyle(21);

auto smearedSB67=GetSmearedCF(pKCF,histSmearTOTSB67);
smearedSB67->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedSB67->SetMarkerColor(2);
smearedSB67->SetMarkerSize(1.5);
smearedSB67->SetMarkerStyle(21);

auto smearedSB78=GetSmearedCF(pKCF,histSmearTOTSB78);
smearedSB78->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedSB78->SetMarkerColor(2);
smearedSB78->SetMarkerSize(1.5);
smearedSB78->SetMarkerStyle(21);

auto smearedSB89=GetSmearedCF(pKCF,histSmearTOTSB89);
smearedSB89->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedSB89->SetMarkerColor(2);
smearedSB89->SetMarkerSize(1.5);
smearedSB89->SetMarkerStyle(21);

auto smearedSB90=GetSmearedCF(pKCF,histSmearTOTSB90);
smearedSB90->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedSB90->SetMarkerColor(2);
smearedSB90->SetMarkerSize(1.5);
smearedSB90->SetMarkerStyle(21);


//int mPairs = 4;
//auto fit1 = new TF1("Fit", [&](double *x, double *p) {

//         return p[0]+p[1]*x[0]+p[2]*TMath::Power(x[0],2)+p[3]*TMath::Power(x[0],3);

//        }, 0, 500, mPairs);
//auto fit2 = new TF1("Fit", [&](double *x, double *p) {

//         return p[0]+p[1]*x[0]+p[2]*TMath::Power(x[0],2)+p[3]*TMath::Power(x[0],3);

//        }, 0, 500, mPairs);
//auto fit3 = new TF1("Fit", [&](double *x, double *p) {

//         return p[0]+p[1]*x[0]+p[2]*TMath::Power(x[0],2)+p[3]*TMath::Power(x[0],3);

//        }, 0, 500, mPairs);


// phitruth->Fit(fit1,"F", "", 5, 400);

// SBleft->Fit(fit2,"F", "", 5, 360);
// SBright->Fit(fit3,"F", "", 5, 360);


//TF1 *f1 = new TF1("f1","landau",0,500);






//--------------------------fit gauss--------------------------------------

//TF1 *f1 = new TF1("f1","gaus",0,450);
//TF1 *f11 = new TF1("f11","gaus(0)+pol0(3)",0,850);
//TF1 *f12 = new TF1("f12","gaus(0)+pol1(3)",0,850);
//TF1 *f13 = new TF1("f13","gaus(0)+pol2(3)",0,850);


//TF1 *total = new TF1("total","gaus(0)+pol2(3)",0,850);


//pKCF->Fit(f1,"R");
//f11->SetParameter(0,f1->GetParameter(0));
//f11->SetParameter(1,f1->GetParameter(1));
//f11->SetParameter(2,f1->GetParameter(2));
//f11->SetParLimits(3,0.9,1.2); //0.072,0.95


//pKCF->Fit(f11,"R");
//f11->SetParLimits(0,0.08,0.10); //0.072,0.95
//f11->SetParLimits(1,220,300); //0.072,0.95
//pKCF->Fit(f11,"R");


//f12->SetParameter(0,f11->GetParameter(0));
//f12->SetParameter(1,f11->GetParameter(1));
//f12->SetParameter(2,f11->GetParameter(2));
//f12->SetParameter(3,f11->GetParameter(3));
//f12->SetParLimits(3,0.5,0.80); //0.072,0.95

//pKCF->Fit(f12,"R");

//f13->SetParameter(0,f12->GetParameter(0));
//f13->SetParameter(1,f12->GetParameter(1));
//f13->SetParameter(2,f12->GetParameter(2));
//f13->SetParameter(3,f12->GetParameter(3));
//f13->SetParLimits(3,0.6,0.80); //0.072,0.95
//f13->SetParameter(4,f12->GetParameter(4));


//pKCF->Fit(f13,"R");


//pKCF->Fit(f13,"R");
//f13->SetParLimits(4,0.0003,0.0007); //0.072,0.95
//f13->SetParLimits(3,0.7,0.75); //0.072,0.95


//pKCF->Fit(f13,"R");


//cout<<"GAUSS FITTTTT----------------"<<endl;
//TF1 *fgp = new TF1("fgp","gaus(0)+pol2(3)",0,650);

//fgp->SetParameter(0,0.214);
//fgp->SetParLimits(1,175,195); //188.117
//fgp->SetParLimits(2,120,160); //141.988
//fgp->SetParameter(3,0.075);
//fgp->SetParLimits(3,0.6,0.80); //0.072,0.95
//fgp->SetParameter(4,0.00069999);
//fgp->SetParLimits(4,0.0007,0.00085);
//fgp->SetParLimits(5,-0.0000006,-0.0000005);
//pKCF66->Fit(fgp,"R");
//pKCF66->Fit(fgp," "," ", 0,750);
////fgp->FixParameter(1,188.117);
//cout<<"GAUSS FINALL----------------"<<endl;
//pKCF66->Fit(fgp," "," ", 0,850);


//----------------------newfit-----------------------------------------------------

TF1 *f13 = new TF1("fgp","gaus(0)+pol2(3)",0,650);

f13->SetParameter(0,0.214);
f13->SetParLimits(1,175,195); //188.117
f13->SetParLimits(2,120,160); //141.988
f13->SetParameter(3,0.075);
f13->SetParLimits(3,0.6,0.80); //0.072,0.95
f13->SetParameter(4,0.00069999);
f13->SetParLimits(4,0.0007,0.00085);
f13->SetParLimits(5,-0.0000006,-0.0000005);
pKCF->Fit(f13,"R");
pKCF->Fit(f13," "," ", 0,750);
//fgp->FixParameter(1,188.117);
cout<<"GAUSS FINALL----------------"<<endl;
pKCF->Fit(f13," "," ", 0,850);



//--------------------------confidencd intervals--------------------------------------

Int_t ngr = 25;
TGraphErrors *grint = new TGraphErrors(ngr);
grint->SetTitle("Fitted line with .95 conf. band");
for (i=0; i<ngr; i++)
   grint->SetPoint(i, pKCF->GetX()[i], 0);
/*Compute the confidence intervals at the x points of the created graph*/
(TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint, 0.68);
//Now the "grint" graph contains function values as its y-coordinates
//and confidence intervals as the errors on these coordinates
//Draw the graph, the function and the confidence intervals
//myc->cd(1);
//grint->SetLineColor(kRed);
//grint->Draw("ap");
//gr->SetMarkerStyle(5);
//gr->SetMarkerSize(0.7);
//gr->Draw("psame");


// smeard intrvals


Double_t cx[ngr]={0};
Double_t cy[ngr]={0};
//Double_t cexu[ngr]={0};
//Double_t cexl[ngr]={0};
Double_t ceyu[ngr]={0};
Double_t ceyl[ngr]={0};
Double_t cUP[ngr]={0};
Double_t cLOW[ngr]={0};


for(int i=0;i<ngr;i++){
    Double_t x,y;
//    cx[i]=grint->GetPointX(i);
//    cy[i]=grint->GetPointY(i);
      grint->GetPoint(i,x,y);
       cx[i]=x;
       cy[i]=y;
//    cexu[i]=grint->GetErrorXhigh(i);
//    cexl[i]=grint->GetErrorXlow(i);
    ceyu[i]=grint->GetErrorYhigh(i);
    ceyl[i]=grint->GetErrorYlow(i);
    cUP[i]=cy[i]+ceyu[i];
    cLOW[i]=cy[i]-ceyl[i];
    std::cout<< "x:"<<cx[i]<<"  y:"<<cy[i]<<endl;
}

auto ConfintUP = new TGraph(ngr,cx,cUP);
DreamPlot::SetStyleGraph(ConfintUP, 20, kBlack, 0.7);
ConfintUP->SetLineWidth(lineWidth2);


auto ConfintLOW = new TGraph(ngr,cx,cLOW);
DreamPlot::SetStyleGraph(ConfintLOW, 20, kBlack, 0.7);
ConfintLOW->SetLineWidth(lineWidth2);


//--------------------------fit BG--------------------------------------

TF1 *f2 = new TF1("f2","pol2",0,850);


f2->SetParLimits(0,0.85,1.00); //0.072,0.95

pKBG->Fit(f2,"R");



//Double_t par[6];
//f1->GetParameters(&par[0]);
//f2->GetParameters(&par[3]);

//f13->GetParameters(&par[0]);


//total->SetLineColor(2);

//total->SetParameters(par);
//total->FixParameter(1,180); //0.072,0.95
//total->SetParLimits(0,1.08,1.14); //0.072,0.95
//total->SetParLimits(1,220,320); //0.072,0.95


//pKCF->Fit(total,"F", "", 0,900);



//--------------------------fit pol5--------------------------------------

TF1 *f3 = new TF1("f3","pol3",0,700);

f3->SetParLimits(0,0.0,0.95); //0.002,0.005

//f1->SetParLimits(1,0.001,0.005); //0.002,0.005
//f1->SetParLimits(2,-0.00003,0.0);
//f1->SetParLimits(3,-0.00000001,0.00000001);
//f1->SetParLimits(4,0.00000000001,0.0000000001);
//f1->SetParLimits(5,0.000000000000001,0.00000000001);
//f1->SetParLimits(6,-1,1);



pKCF2->Fit(f3,"F", "", 0, 700);
f3->SetParLimits(0,0.5,0.95);
pKCF2->Fit(f3,"F", "", 0, 700);

f3->SetParLimits(0,0.6,0.95);

f3->SetParLimits(1,0.002,0.003);
f3->SetParLimits(2,-0.000007,-0.000003);
f3->SetParLimits(3,0.0,0.000000007);

pKCF2->Fit(f3,"F", "", 0, 910); //910

TF1 *f4 = new TF1("f4","pol4",0,950);
f4->SetParLimits(0,0.6,0.95);

f4->SetParameter(0,f3->GetParameter(0));
f4->SetParameter(1,f3->GetParameter(1));
f4->SetParameter(2,f3->GetParameter(2));
f4->SetParameter(3,f3->GetParameter(3));


pKCF2->Fit(f4,"F", "", 0, 910);

TF1 *f5 = new TF1("f5","pol5",0,780);
f5->SetParLimits(0,0.6,0.95);

f5->SetParameter(0,f4->GetParameter(0));
f5->SetParameter(1,f4->GetParameter(1));
f5->SetParameter(2,f4->GetParameter(2));
f5->SetParameter(3,f4->GetParameter(3));
f5->SetParameter(4,f4->GetParameter(4));


pKCF2->Fit(f5,"F", "", 0, 780);



//f3->SetParLimits(0,0.7,0.95);
//pKCF2->Fit(f3,"F", "", 0, 910);
//f3->SetParLimits(0,0.75,0.95);
//pKCF2->Fit(f3,"F", "", 0, 910);

//TF1 *f4 = new TF1("f4","pol5",0,700);
//f3->SetParLimits(0,0.5,0.95);


//--------------------------fit Dimi fitfunction pol3--------------------------------------


TF1* fun2 = new TF1("fun2",PolN_flat_end,0,780);

fun2->FixParameter(0,3);
fun2->FixParameter(1,600);
fun2->SetParameter(2,f3->GetParameter(0));//0
fun2->SetParameter(3,f3->GetParameter(1));//1
fun2->SetParameter(4,f3->GetParameter(2));//2
fun2->SetParameter(5,f3->GetParameter(3));//3
fun2->SetParLimits(2,0.0,0.95); //0.002,0.005


pKCFDimi->Fit(fun2,"F", "", 0, 780);




//--------------------------fit Dimi fitfunction pol5--------------------------------------


TF1* funDpol5 = new TF1("funDpol5",PolN_flat_end,0,780);

fun2->FixParameter(0,5);
fun2->FixParameter(1,600);
fun2->SetParameter(2,f5->GetParameter(0));//0
fun2->SetParameter(3,f5->GetParameter(1));//1
fun2->SetParameter(4,f5->GetParameter(2));//2
fun2->SetParameter(5,f5->GetParameter(3));//3
fun2->SetParameter(6,f5->GetParameter(4));//4
fun2->SetParameter(7,f5->GetParameter(5));//5

fun2->SetParLimits(2,0.0,0.95); //0.002,0.005


pKCFDimipol5->Fit(funDpol5,"F", "", 0, 780);



//-------------------------creat fit graphs--------------------------------------

const int binwidth = 1;
const int NumMomBins = 1000;
float p1=f13->GetParameter(0);
float p2=f13->GetParameter(1);
float p3=f13->GetParameter(2);
float p4=f13->GetParameter(3);
float p5=f13->GetParameter(4);
float p6=f13->GetParameter(5);

//float p6=f1->GetParameter(5);


float p1s=f5->GetParameter(0);
float p2s=f5->GetParameter(1);
float p3s=f5->GetParameter(2);
float p4s=f5->GetParameter(3);
float p5s=f5->GetParameter(4);
float p6s=f5->GetParameter(5);


//float dd=fit1->GetParameter(3);

Int_t nn =NumMomBins;
Double_t x2[nn], y2[nn];
for (Int_t i=0; i<nn; i++) {
   x2[i] = i*binwidth;
   y2[i]=p1s+p2s*x2[i]+p3s*x2[i]*x2[i]+p4s*x2[i]*x2[i]*x2[i]+p5s*x2[i]*x2[i]*x2[i]*x2[i]+p6s*x2[i]*x2[i]*x2[i]*x2[i]*x2[i];


  // cout<< x2[i]<< "    y= "<<y2[i]<<endl;
}


Double_t x3[nn], y3[nn];

for (Int_t i=0; i<nn; i++) {
   x3[i] = i*binwidth;
   y3[i] = p1*TMath::Gaus(x3[i], p2, p3, false)+p4+p5*x3[i]+p6*x3[i]*x3[i];


  // cout<< x2[i]<< "    y= "<<y2[i]<<endl;
}




TGraph *fitfunpol = new TGraph (nn, x2, y2);
TGraph *fitfungaus = new TGraph (nn, x3, y3);


auto smearedfitCF=GetSmearedCF(fitfungaus,histSmearTOTPeak);
auto smearedfitSBleft=GetSmearedCF(fitfungaus,histSmearTOTSBleft);
auto smearedfitSBright=GetSmearedCF(fitfungaus,histSmearTOTSBright);
auto smearedfitSB12=GetSmearedCF(fitfungaus,histSmearTOTSB12);
auto smearedfitSB23=GetSmearedCF(fitfungaus,histSmearTOTSB23);
auto smearedfitSB34=GetSmearedCF(fitfungaus,histSmearTOTSB34);
auto smearedfitSB45=GetSmearedCF(fitfungaus,histSmearTOTSB45);
auto smearedfitSB56=GetSmearedCF(fitfungaus,histSmearTOTSB56);
auto smearedfitSB67=GetSmearedCF(fitfungaus,histSmearTOTSB67);
auto smearedfitSB78=GetSmearedCF(fitfungaus,histSmearTOTSB78);
auto smearedfitSB89=GetSmearedCF(fitfungaus,histSmearTOTSB89);
auto smearedfitSB90=GetSmearedCF(fitfungaus,histSmearTOTSB90);


auto smearedconfSBleftup=GetSmearedCF(ConfintUP,histSmearTOTSBleft);
auto smearedconfSBleftlow=GetSmearedCF(ConfintLOW,histSmearTOTSBleft);

auto smearedconfSBrightup=GetSmearedCF(ConfintUP,histSmearTOTSBright);
auto smearedconfSBrightlow=GetSmearedCF(ConfintLOW,histSmearTOTSBright);

smearedconfSBleftup->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedconfSBleftup->SetMarkerColor(2);
smearedconfSBleftup->SetMarkerSize(1.5);
smearedconfSBleftup->SetMarkerStyle(21);

smearedconfSBleftlow->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedconfSBleftlow->SetMarkerColor(2);
smearedconfSBleftlow->SetMarkerSize(1.5);
smearedconfSBleftlow->SetMarkerStyle(21);

smearedconfSBrightup->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedconfSBrightup->SetMarkerColor(2);
smearedconfSBrightup->SetMarkerSize(1.5);
smearedconfSBrightup->SetMarkerStyle(21);

smearedconfSBrightlow->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedconfSBrightlow->SetMarkerColor(2);
smearedconfSBrightlow->SetMarkerSize(1.5);
smearedconfSBrightlow->SetMarkerStyle(21);


smearedfitCF->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedfitCF->SetMarkerColor(2);
smearedfitCF->SetMarkerSize(1.5);
smearedfitCF->SetMarkerStyle(21);


smearedfitCF->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedfitCF->SetMarkerColor(2);
smearedfitCF->SetMarkerSize(1.5);
smearedfitCF->SetMarkerStyle(21);

smearedfitSBleft->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedfitSBleft->SetMarkerColor(2);
smearedfitSBleft->SetMarkerSize(1.5);
smearedfitSBleft->SetMarkerStyle(21);

smearedfitSBright->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedfitSBright->SetMarkerColor(2);
smearedfitSBright->SetMarkerSize(1.5);
smearedfitSBright->SetMarkerStyle(21);

smearedfitSB12->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedfitSB12->SetMarkerColor(2);
smearedfitSB12->SetMarkerSize(1.5);
smearedfitSB12->SetMarkerStyle(21);

smearedfitSB23->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedfitSB23->SetMarkerColor(2);
smearedfitSB23->SetMarkerSize(1.5);
smearedfitSB23->SetMarkerStyle(21);

smearedfitSB34->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedfitSB34->SetMarkerColor(2);
smearedfitSB34->SetMarkerSize(1.5);
smearedfitSB34->SetMarkerStyle(21);


smearedfitSB45->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedfitSB45->SetMarkerColor(2);
smearedfitSB45->SetMarkerSize(1.5);
smearedfitSB45->SetMarkerStyle(21);


smearedfitSB56->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedfitSB56->SetMarkerColor(2);
smearedfitSB56->SetMarkerSize(1.5);
smearedfitSB56->SetMarkerStyle(21);


smearedfitSB67->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedfitSB67->SetMarkerColor(2);
smearedfitSB67->SetMarkerSize(1.5);
smearedfitSB67->SetMarkerStyle(21);


smearedfitSB78->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedfitSB78->SetMarkerColor(2);
smearedfitSB78->SetMarkerSize(1.5);
smearedfitSB78->SetMarkerStyle(21);


smearedfitSB89->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedfitSB89->SetMarkerColor(2);
smearedfitSB89->SetMarkerSize(1.5);
smearedfitSB89->SetMarkerStyle(21);


smearedfitSB90->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedfitSB90->SetMarkerColor(2);
smearedfitSB90->SetMarkerSize(1.5);
smearedfitSB90->SetMarkerStyle(21);



auto smearedfitCFpol=GetSmearedCF(fitfunpol,histSmearTOTPeak);
auto smearedfitSBleftpol=GetSmearedCF(fitfunpol,histSmearTOTSBleft);
auto smearedfitSBrightpol=GetSmearedCF(fitfunpol,histSmearTOTSBright);
auto smearedfitSB12pol=GetSmearedCF(fitfunpol,histSmearTOTSB12);
auto smearedfitSB23pol=GetSmearedCF(fitfunpol,histSmearTOTSB23);
auto smearedfitSB34pol=GetSmearedCF(fitfunpol,histSmearTOTSB34);
auto smearedfitSB45pol=GetSmearedCF(fitfunpol,histSmearTOTSB45);
auto smearedfitSB56pol=GetSmearedCF(fitfunpol,histSmearTOTSB56);
auto smearedfitSB67pol=GetSmearedCF(fitfunpol,histSmearTOTSB67);
auto smearedfitSB78pol=GetSmearedCF(fitfunpol,histSmearTOTSB78);
auto smearedfitSB89pol=GetSmearedCF(fitfunpol,histSmearTOTSB89);
auto smearedfitSB90pol=GetSmearedCF(fitfunpol,histSmearTOTSB90);

smearedfitCFpol->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedfitCFpol->SetMarkerColor(2);
smearedfitCFpol->SetMarkerSize(1.5);
smearedfitCFpol->SetMarkerStyle(21);

smearedfitSBleftpol->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedfitSBleftpol->SetMarkerColor(2);
smearedfitSBleftpol->SetMarkerSize(1.5);
smearedfitSBleftpol->SetMarkerStyle(21);

smearedfitSBrightpol->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedfitSBrightpol->SetMarkerColor(2);
smearedfitSBrightpol->SetMarkerSize(1.5);
smearedfitSBrightpol->SetMarkerStyle(21);

smearedfitSB12pol->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedfitSB12pol->SetMarkerColor(2);
smearedfitSB12pol->SetMarkerSize(1.5);
smearedfitSB12pol->SetMarkerStyle(21);

smearedfitSB23pol->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedfitSB23pol->SetMarkerColor(2);
smearedfitSB23pol->SetMarkerSize(1.5);
smearedfitSB23pol->SetMarkerStyle(21);

smearedfitSB34pol->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedfitSB34pol->SetMarkerColor(2);
smearedfitSB34pol->SetMarkerSize(1.5);
smearedfitSB34pol->SetMarkerStyle(21);


smearedfitSB45pol->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedfitSB45pol->SetMarkerColor(2);
smearedfitSB45pol->SetMarkerSize(1.5);
smearedfitSB45pol->SetMarkerStyle(21);


smearedfitSB56pol->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedfitSB56pol->SetMarkerColor(2);
smearedfitSB56pol->SetMarkerSize(1.5);
smearedfitSB56pol->SetMarkerStyle(21);


smearedfitSB67pol->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedfitSB67pol->SetMarkerColor(2);
smearedfitSB67pol->SetMarkerSize(1.5);
smearedfitSB67pol->SetMarkerStyle(21);


smearedfitSB78pol->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedfitSB78pol->SetMarkerColor(2);
smearedfitSB78pol->SetMarkerSize(1.5);
smearedfitSB78pol->SetMarkerStyle(21);


smearedfitSB89pol->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedfitSB89pol->SetMarkerColor(2);
smearedfitSB89pol->SetMarkerSize(1.5);
smearedfitSB89pol->SetMarkerStyle(21);


smearedfitSB90pol->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedfitSB90pol->SetMarkerColor(2);
smearedfitSB90pol->SetMarkerSize(1.5);
smearedfitSB90pol->SetMarkerStyle(21);



//auto smearedCF3=GetSmearedCF2(fitfun,histSmearTOTPeak);
//smearedCF3->SetLineColor(3);
////smearedCF->SetLineWidth(4);
//smearedCF3->SetMarkerColor(3);
//smearedCF3->SetMarkerSize(1.5);
//smearedCF3->SetMarkerStyle(21);



//auto smearedCF2=GetSmearedCF(total,histSmearTOTPeak);
//smearedCF2->SetLineColor(2);
////smearedCF->SetLineWidth(4);
//smearedCF2->SetMarkerColor(2);
//smearedCF2->SetMarkerSize(1.5);
//smearedCF2->SetMarkerStyle(21);

double max_x=getxmax(smearedfitCF);
double max_y=getymax(smearedfitCF);

double reX=std::abs(260-max_x);
//double reY=std::abs(1.035-max_y);
double reY=std::abs(1.041-max_y);
cout<<"rely : " << reY<<endl;
auto movedCF=GetSmearedCFmoved(fitfungaus,histSmearTOTPeak,reX,reY);
movedCF->SetLineColor(2);
//smearedCF->SetLineWidth(4);
movedCF->SetMarkerColor(2);
movedCF->SetMarkerSize(1.5);
movedCF->SetMarkerStyle(21);

double max_x2=getxmax(smearedCF);
double max_y2=getymax(smearedCF);

double reX2=std::abs(260-max_x2);
//double reY2=std::abs(1.035-max_y2);
double reY2=std::abs(1.041-max_y2);

cout<<"rely : " << reY<<endl;
auto movedCF2=GetSmearedCFmoved(pKCF,histSmearTOTPeak,reX2,reY2);
movedCF2->SetLineColor(2);
//smearedCF->SetLineWidth(4);
movedCF2->SetMarkerColor(2);
movedCF2->SetMarkerSize(1.5);
movedCF2->SetMarkerStyle(21);





double max_xl=getxmax(smearedfitSBleft);
double max_yl=getymax(smearedfitSBleft);

double reXl=std::abs(340-max_xl);
double reYl=std::abs(1.073-max_yl);
//double reYl=std::abs(1.071-max_yl);

auto movedleft=GetSmearedCFmoved(fitfungaus,histSmearTOTSBleft,reXl,reYl);
movedleft->SetLineColor(2);
//smearedCF->SetLineWidth(4);
movedleft->SetMarkerColor(2);
movedleft->SetMarkerSize(1.5);
movedleft->SetMarkerStyle(21);

double max_xl2=getxmax(smearedSBleft);
double max_yl2=getymax(smearedSBleft);

double reXl2=std::abs(340-max_xl2);
double reYl2=std::abs(1.073-max_yl2);
//double reYl2=std::abs(1.071-max_yl2);

auto movedleft2=GetSmearedCFmoved(pKCF,histSmearTOTSBleft,reXl2,reYl2);
movedleft2->SetLineColor(2);
//smearedCF->SetLineWidth(4);
movedleft2->SetMarkerColor(2);
movedleft2->SetMarkerSize(1.5);
movedleft2->SetMarkerStyle(21);






double max_xr=getxmax(smearedfitSBright);
double max_yr=getymax(smearedfitSBright);

double reXr=std::abs(300-max_xr);
//double reYr=std::abs(1.056-max_yr);
double reYr=std::abs(1.058-max_yr);

auto movedright=GetSmearedCFmoved(fitfungaus,histSmearTOTSBright,reXr,reYr);
movedright->SetLineColor(2);
movedright->SetMarkerColor(2);
movedright->SetMarkerSize(1.5);
movedright->SetMarkerStyle(21);


double max_xr2=getxmax(smearedSBright);
double max_yr2=getymax(smearedSBright);

double reXr2=std::abs(300-max_xr2);
//double reYr2=std::abs(1.056-max_yr2);
double reYr2=std::abs(1.058-max_yr2);

auto movedright2=GetSmearedCFmoved(pKCF,histSmearTOTSBright,reXr2,reYr2);
movedright2->SetLineColor(2);
movedright2->SetMarkerColor(2);
movedright2->SetMarkerSize(1.5);
movedright2->SetMarkerStyle(21);







double max_x12=getxmax(smearedfitSB12);
double max_y12=getymax(smearedfitSB12);

double reX12=std::abs(300-max_x12);
double reY12=std::abs(1.056-max_y12);
auto moved12=GetSmearedCFmoved(fitfungaus,histSmearTOTSB12,reX12,reY12);
moved12->SetLineColor(2);
moved12->SetMarkerColor(2);
moved12->SetMarkerSize(1.5);
moved12->SetMarkerStyle(21);


double max_x23=getxmax(smearedfitSB23);
double max_y23=getymax(smearedfitSB23);

double reX23=std::abs(300-max_x23);
double reY23=std::abs(1.056-max_y23);
auto moved23=GetSmearedCFmoved(fitfungaus,histSmearTOTSB23,reX23,reY23);
moved23->SetLineColor(2);
moved23->SetMarkerColor(2);
moved23->SetMarkerSize(1.5);
moved23->SetMarkerStyle(21);


double max_x34=getxmax(smearedfitSB34);
double max_y34=getymax(smearedfitSB34);

double reX34=std::abs(300-max_x34);
double reY34=std::abs(1.056-max_y34);
auto moved34=GetSmearedCFmoved(fitfungaus,histSmearTOTSB34,reX34,reY34);
moved34->SetLineColor(2);
moved34->SetMarkerColor(2);
moved34->SetMarkerSize(1.5);
moved34->SetMarkerStyle(21);


double max_x45=getxmax(smearedfitSB45);
double max_y45=getymax(smearedfitSB45);

double reX45=std::abs(300-max_x45);
double reY45=std::abs(1.056-max_y45);
auto moved45=GetSmearedCFmoved(fitfungaus,histSmearTOTSB45,reX45,reY45);
moved45->SetLineColor(2);
moved45->SetMarkerColor(2);
moved45->SetMarkerSize(1.5);
moved45->SetMarkerStyle(21);


double max_x56=getxmax(smearedfitSB56);
double max_y56=getymax(smearedfitSB56);

double reX56=std::abs(300-max_x56);
double reY56=std::abs(1.056-max_y56);
auto moved56=GetSmearedCFmoved(fitfungaus,histSmearTOTSB56,reX56,reY56);
moved56->SetLineColor(2);
moved56->SetMarkerColor(2);
moved56->SetMarkerSize(1.5);
moved56->SetMarkerStyle(21);


double max_x67=getxmax(smearedfitSB67);
double max_y67=getymax(smearedfitSB67);

double reX67=std::abs(300-max_x67);
double reY67=std::abs(1.056-max_y67);
auto moved67=GetSmearedCFmoved(fitfungaus,histSmearTOTSB67,reX67,reY67);
moved67->SetLineColor(2);
moved67->SetMarkerColor(2);
moved67->SetMarkerSize(1.5);
moved67->SetMarkerStyle(21);

double max_x78=getxmax(smearedfitSB78);
double max_y78=getymax(smearedfitSB78);

double reX78=std::abs(300-max_x78);
double reY78=std::abs(1.056-max_y78);
auto moved78=GetSmearedCFmoved(fitfungaus,histSmearTOTSB78,reX78,reY78);
moved78->SetLineColor(2);
moved78->SetMarkerColor(2);
moved78->SetMarkerSize(1.5);
moved78->SetMarkerStyle(21);


double max_x89=getxmax(smearedfitSB89);
double max_y89=getymax(smearedfitSB89);

double reX89=std::abs(300-max_x89);
double reY89=std::abs(1.056-max_y89);
auto moved89=GetSmearedCFmoved(fitfungaus,histSmearTOTSB89,reX89,reY89);
moved89->SetLineColor(2);
moved89->SetMarkerColor(2);
moved89->SetMarkerSize(1.5);
moved89->SetMarkerStyle(21);


double max_x90=getxmax(smearedfitSB90);
double max_y90=getymax(smearedfitSB90);

double reX90=std::abs(300-max_x90);
double reY90=std::abs(1.056-max_y90);
auto moved90=GetSmearedCFmoved(fitfungaus,histSmearTOTSB90,reX90,reY90);
moved90->SetLineColor(2);
moved90->SetMarkerColor(2);
moved90->SetMarkerSize(1.5);
moved90->SetMarkerStyle(21);



cout<<"npoints: "<< smearedfitCF->GetN()<< endl;

//double IntegralCF = smearedfitCF->Integral(400,
//                                  450);

double IntegralCF =integrate(smearedfitCF,800,900);
Double_t scale = 1/IntegralCF;

//cout<<"integral:"<< IntegralCF<< "scale: "<< scale<< endl;
//smearedfitCF->Scale(scale);

auto normCF=GetSmearedCFscaled(fitfungaus, histSmearTOTPeak, scale);
normCF->SetLineColor(2);
normCF->SetMarkerColor(2);
normCF->SetMarkerSize(1.5);
normCF->SetMarkerStyle(21);


double IntegralCF2 =integrate(smearedCF,800,900);
Double_t scale2 = 1/IntegralCF;

//cout<<"integral:"<< IntegralCF<< "scale: "<< scale<< endl;
//smearedfitCF->Scale(scale);

auto normCF2=GetSmearedCFscaled(pKCF, histSmearTOTPeak, scale2);
normCF2->SetLineColor(2);
normCF2->SetMarkerColor(2);
normCF2->SetMarkerSize(1.5);
normCF2->SetMarkerStyle(21);







double IntegralCFleft =integrate(smearedfitSBleft,800,900);
Double_t scaleleft = 1/IntegralCFleft;

auto normSBleft=GetSmearedCFscaled(fitfungaus, histSmearTOTSBleft, scaleleft);
normSBleft->SetLineColor(2);
normSBleft->SetMarkerColor(2);
normSBleft->SetMarkerSize(1.5);
normSBleft->SetMarkerStyle(21);

double IntegralCFleft2 =integrate(smearedSBleft,800,900);
Double_t scaleleft2 = 1/IntegralCFleft2;

auto normSBleft2=GetSmearedCFscaled(pKCF, histSmearTOTSBleft, scaleleft2);
normSBleft2->SetLineColor(2);
normSBleft2->SetMarkerColor(2);
normSBleft2->SetMarkerSize(1.5);
normSBleft2->SetMarkerStyle(21);



double IntegralCFright =integrate(smearedfitSBright,800,900);
Double_t scaleright = 1/IntegralCFright;

auto normSBright=GetSmearedCFscaled(fitfungaus, histSmearTOTSBright, scaleright);
normSBright->SetLineColor(2);
normSBright->SetMarkerColor(2);
normSBright->SetMarkerSize(1.5);
normSBright->SetMarkerStyle(21);

double IntegralCFright2 =integrate(smearedSBright,800,900);
Double_t scaleright2 = 1/IntegralCFright2;

auto normSBright2=GetSmearedCFscaled(pKCF, histSmearTOTSBright, scaleright2);
normSBright2->SetLineColor(2);
normSBright2->SetMarkerColor(2);
normSBright2->SetMarkerSize(1.5);
normSBright2->SetMarkerStyle(21);




Double_t xnor,ynor, xo,yo;

normSBright->GetPoint(300, xnor, ynor);
smearedfitSBleft->GetPoint(300, xo, yo);


float normghtsb=ynor;
float orirightsb=yo;

Double_t xnorle,ynorle, xole,yole;

normSBright->GetPoint(300, xnorle, ynorle);
smearedfitSBleft->GetPoint(300, xole, yole);


float normleftsb=ynorle;
float orileftsb=yole;

float yshiftright=abs(normghtsb-orirightsb);

float yshiftleft=abs(normleftsb-orileftsb);

cout<< "yshift right:"<<yshiftright<<endl;

double IntegralSBleftUP =integrate(smearedconfSBleftup,800,900);
Double_t scaleleftUP= 1/IntegralSBleftUP;

auto normSBleftUP=GetSmearedCFmoved(ConfintUP, histSmearTOTSBleft, 0,yshiftleft);
normSBleftUP->SetLineColor(2);
normSBleftUP->SetMarkerColor(2);
normSBleftUP->SetMarkerSize(1.5);
normSBleftUP->SetMarkerStyle(21);



double IntegralSBleftLOW =integrate(smearedconfSBleftlow,800,900);
Double_t scaleleftLOW= 1/IntegralSBleftLOW;

auto normSBleftLOW=GetSmearedCFmoved(ConfintLOW, histSmearTOTSBleft,0, yshiftleft);
normSBleftLOW->SetLineColor(2);
normSBleftLOW->SetMarkerColor(2);
normSBleftLOW->SetMarkerSize(1.5);
normSBleftLOW->SetMarkerStyle(21);



double IntegralSBrightUP =integrate(smearedconfSBrightup,800,900);
Double_t scalerightup= 1/IntegralSBrightUP;

auto normSBrightUP=GetSmearedCFmoved(ConfintUP, histSmearTOTSBright, 0, yshiftright);
normSBrightUP->SetLineColor(2);
normSBrightUP->SetMarkerColor(2);
normSBrightUP->SetMarkerSize(1.5);
normSBrightUP->SetMarkerStyle(21);



double IntegralSBrightLOW =integrate(smearedconfSBrightlow,800,900);
Double_t scalerightLOW= 1/IntegralSBrightLOW;

auto normSBrightLOW=GetSmearedCFmoved(ConfintLOW, histSmearTOTSBright, 0,yshiftright);
normSBrightLOW->SetLineColor(2);
normSBrightLOW->SetMarkerColor(2);
normSBrightLOW->SetMarkerSize(1.5);
normSBrightLOW->SetMarkerStyle(21);




int lineWidth = 3;


DreamPlot::SetStyleGraph(smearedCF, 20, kGray, 0.7);
smearedCF->SetLineWidth(lineWidth);
DreamPlot::SetStyleGraph(smearedSBleft, 20, kGray, 0.7);
smearedSBleft->SetLineWidth(4);
DreamPlot::SetStyleGraph(smearedSBright, 20, kGray, 0.7);
smearedSBright->SetLineWidth(4);
DreamPlot::SetStyleGraph(smearedSB12, 20, kGray, 0.7);
smearedSB12->SetLineWidth(4);
DreamPlot::SetStyleGraph(smearedSB23, 20, kGray, 0.7);
smearedSB23->SetLineWidth(4);
DreamPlot::SetStyleGraph(smearedSB34, 20, kGray, 0.7);
smearedSB34->SetLineWidth(4);
DreamPlot::SetStyleGraph(smearedSB45, 20, kGray, 0.7);
smearedSB45->SetLineWidth(4);
DreamPlot::SetStyleGraph(smearedSB56, 20, kGray, 0.7);
smearedSB56->SetLineWidth(4);
DreamPlot::SetStyleGraph(smearedSB67, 20, kGray, 0.7);
smearedSB67->SetLineWidth(4);
DreamPlot::SetStyleGraph(smearedSB78, 20, kGray, 0.7);
smearedSB78->SetLineWidth(4);
DreamPlot::SetStyleGraph(smearedSB89, 20, kGray, 0.7);
smearedSB89->SetLineWidth(4);
DreamPlot::SetStyleGraph(smearedSB90, 20, kGray, 0.7);
smearedSB90->SetLineWidth(4);


DreamPlot::SetStyleGraph(fitfungaus, 20, kGreen, 0.7);
fitfungaus->SetLineWidth(lineWidth);
DreamPlot::SetStyleGraph(smearedfitCF, 20, kRed, 0.7);
smearedfitCF->SetLineWidth(lineWidth);
DreamPlot::SetStyleGraph(smearedfitSBleft, 20, kRed, 0.7);
smearedfitSBleft->SetLineWidth(4);
DreamPlot::SetStyleGraph(smearedfitSBright, 20, kRed, 0.7);
smearedfitSBright->SetLineWidth(4);
DreamPlot::SetStyleGraph(smearedfitSB12, 20, kRed, 0.7);
smearedfitSB12->SetLineWidth(4);
DreamPlot::SetStyleGraph(smearedfitSB23, 20, kRed, 0.7);
smearedfitSB23->SetLineWidth(4);
DreamPlot::SetStyleGraph(smearedfitSB34, 20, kRed, 0.7);
smearedfitSB34->SetLineWidth(4);
DreamPlot::SetStyleGraph(smearedfitSB45, 20, kRed, 0.7);
smearedfitSB45->SetLineWidth(4);
DreamPlot::SetStyleGraph(smearedfitSB56, 20, kRed, 0.7);
smearedfitSB56->SetLineWidth(4);
DreamPlot::SetStyleGraph(smearedfitSB67, 20, kRed, 0.7);
smearedfitSB67->SetLineWidth(4);
DreamPlot::SetStyleGraph(smearedfitSB78, 20, kRed, 0.7);
smearedfitSB78->SetLineWidth(4);
DreamPlot::SetStyleGraph(smearedfitSB89, 20, kRed, 0.7);
smearedfitSB89->SetLineWidth(4);
DreamPlot::SetStyleGraph(smearedfitSB90, 20, kRed, 0.7);
smearedfitSB90->SetLineWidth(4);




DreamPlot::SetStyleGraph(fitfunpol, 20, kGreen, 0.7);
fitfunpol->SetLineWidth(lineWidth);
DreamPlot::SetStyleGraph(smearedfitCFpol, 20, kBlue, 0.7);
smearedfitCFpol->SetLineWidth(4);
DreamPlot::SetStyleGraph(smearedfitSBleftpol, 20, kBlue, 0.7);
smearedfitSBleftpol->SetLineWidth(4);
DreamPlot::SetStyleGraph(smearedfitSBrightpol, 20, kBlue, 0.7);
smearedfitSBrightpol->SetLineWidth(4);
DreamPlot::SetStyleGraph(smearedfitSB12pol, 20, kBlue, 0.7);
smearedfitSB12pol->SetLineWidth(4);
DreamPlot::SetStyleGraph(smearedfitSB23pol, 20, kBlue, 0.7);
smearedfitSB23pol->SetLineWidth(4);
DreamPlot::SetStyleGraph(smearedfitSB34pol, 20, kBlue, 0.7);
smearedfitSB34pol->SetLineWidth(4);
DreamPlot::SetStyleGraph(smearedfitSB45pol, 20, kBlue, 0.7);
smearedfitSB45pol->SetLineWidth(4);
DreamPlot::SetStyleGraph(smearedfitSB56pol, 20, kBlue, 0.7);
smearedfitSB56pol->SetLineWidth(4);
DreamPlot::SetStyleGraph(smearedfitSB67pol, 20, kBlue, 0.7);
smearedfitSB67pol->SetLineWidth(4);
DreamPlot::SetStyleGraph(smearedfitSB78pol, 20, kBlue, 0.7);
smearedfitSB78pol->SetLineWidth(4);
DreamPlot::SetStyleGraph(smearedfitSB89pol, 20, kBlue, 0.7);
smearedfitSB89pol->SetLineWidth(4);
DreamPlot::SetStyleGraph(smearedfitSB90pol, 20, kBlue, 0.7);
smearedfitSB90pol->SetLineWidth(4);







auto gc3 = new TCanvas("pKCF", "pKCF", 0, 0, 650, 550);
gc3->SetRightMargin(right);
gc3->SetTopMargin(top);
dummyHist->Draw();
dummyHist->GetYaxis()->SetRangeUser(0.85, 1.2);
//dummyHist->GetXaxis()->SetRangeUser(0.0, 200);

//dummyHist->GetXaxis()->SetNdivisions(504);

pKCF->SetTitle("p-K function with gaus+pol2 fit");
//pKCF->SetMarkerColor(4);
pKCF->SetMarkerStyle(21);
pKCF->Draw("ALPsame");
fitfungaus->Draw("same");
//pKCF->Draw("L3same");
gc3->Print("pKCF.pdf");
gc3->Print("pKCF.png");


auto gc33 = new TCanvas("pKCF2", "pKCF2", 0, 0, 650, 550);
gc33->SetRightMargin(right);
gc33->SetTopMargin(top);
dummyHist->Draw();
//dummyHist->GetYaxis()->SetRangeUser(0.8, 1.1);
//dummyHist->GetXaxis()->SetRangeUser(0.0, 200);

//dummyHist->GetXaxis()->SetNdivisions(504);

pKCF2->SetTitle("p-K function with pol5 fit");
pKCF2->SetMarkerColor(4);
pKCF2->SetMarkerStyle(21);
pKCF2->Draw("ALPsame");
fitfunpol->Draw("same");
//pKCF->Draw("L3same");
gc33->Print("pKCF2.pdf");
gc33->Print("pKCF2.png");



auto gc335 = new TCanvas("pKCFDimi", "pKCFDimi", 0, 0, 650, 550);
gc335->SetRightMargin(right);
gc335->SetTopMargin(top);
dummyHist->Draw();
//dummyHist->GetYaxis()->SetRangeUser(0.8, 1.1);
//dummyHist->GetXaxis()->SetRangeUser(0.0, 200);

//dummyHist->GetXaxis()->SetNdivisions(504);

pKCFDimi->SetTitle("p-K function with pol5 fit");
pKCFDimi->SetMarkerColor(4);
pKCFDimi->SetMarkerStyle(21);
pKCFDimi->Draw("ALPsame");
fitfunpol->Draw("same");
//pKCF->Draw("L3same");
gc335->Print("pKCFDimi.pdf");
gc335->Print("pKCFDimi.png");


auto gc33522 = new TCanvas("pKCFDimipol5", "pKCFDimipol5", 0, 0, 650, 550);
gc33522->SetRightMargin(right);
gc33522->SetTopMargin(top);
dummyHist->Draw();
//dummyHist->GetYaxis()->SetRangeUser(0.8, 1.1);
//dummyHist->GetXaxis()->SetRangeUser(0.0, 200);

//dummyHist->GetXaxis()->SetNdivisions(504);

pKCFDimipol5->SetTitle("p-K function with pol5 fit");
pKCFDimipol5->SetMarkerColor(4);
pKCFDimipol5->SetMarkerStyle(21);
pKCFDimipol5->Draw("ALPsame");
fitfunpol->Draw("same");
//pKCF->Draw("L3same");
gc33522->Print("pKCFDimipol5.pdf");
gc33522->Print("pKCFDimipol5.png");






auto gc333 = new TCanvas("pKBG", "pKBG", 0, 0, 650, 550);
gc333->SetRightMargin(right);
gc333->SetTopMargin(top);
dummyHist->Draw();
//dummyHist->GetYaxis()->SetRangeUser(0.8, 1.1);
//dummyHist->GetXaxis()->SetRangeUser(0.0, 200);

//dummyHist->GetXaxis()->SetNdivisions(504);

pKBG->SetTitle("pol2 BG fit");
pKBG->SetMarkerColor(4);
pKBG->SetMarkerStyle(21);
pKBG->Draw("ALPsame");
//pKCF->Draw("L3same");
gc333->Print("pKBG.pdf");
gc333->Print("pKBG.png");


//dummyHist->GetXaxis()->SetRangeUser(5.0, 1000);


auto gc3336 = new TCanvas("pKBG", "pKBG", 0, 0, 650, 550);
gc3336->SetRightMargin(right);
gc3336->SetTopMargin(top);
dummyHist->Draw();
//dummyHist->GetYaxis()->SetRangeUser(0.8, 1.1);
//dummyHist->GetXaxis()->SetRangeUser(0.0, 200);

//dummyHist->GetXaxis()->SetNdivisions(504);
grint->SetLineColor(kRed);
grint->Draw("ap");
pKCF->SetMarkerStyle(5);
pKCF->SetMarkerSize(0.7);
pKCF->Draw("ALPsame");
gc3336->Print("confint.pdf");
gc3336->Print("confint.png");

ConfintUP->SetLineColor(kRed);
ConfintLOW->SetLineColor(kRed);


auto totPeak = new TCanvas("TOTsmearingPeak", "TOTsmearingPeak", 0, 0, 650, 550);
totPeak->SetRightMargin(right);
totPeak->SetTopMargin(top);

dummyHist->Draw();
//dummyHist->GetYaxis()->SetRangeUser(0.85, 1.05);
//dummyHist->GetXaxis()->SetNdivisions(504);
pPhiCFc->SetTitle("p-K#rightarrow p-#phi");
pPhiCFc->Draw("same");
smearedfitCF->Draw("L3same");
//smearedfitCFpol->Draw("L3same");
//smearedCF->Draw("L3same");

//fitfun2->Draw("L3same");
auto totPeakleg = new TLegend(0.65, 0.2, 0.65 + 0.25, 0.45);
totPeakleg->SetBorderSize(0);
totPeakleg->SetTextFont(42);
totPeakleg->SetTextSize(gStyle->GetTextSize() * 0.7);
totPeakleg->SetFillStyle(0);
totPeakleg->AddEntry(pPhiCFc,"p-#phi", "pez");
totPeakleg->AddEntry(smearedfitCF,"gaus+pol2","l");
//gle1g55->AddEntry(fitfun2,"Fit function", "l");
totPeakleg->Draw("same");
totPeak->Print("TOTsmearingPeak.pdf");
totPeak->Print("TOTsmearingPeak.png");

auto totPeak2 = new TCanvas("TOTsmearingPeak2", "TOTsmearingPeak2", 0, 0, 650, 550);
totPeak2->SetRightMargin(right);
totPeak2->SetTopMargin(top);

dummyHist->Draw();
//dummyHist->GetYaxis()->SetRangeUser(0.85, 1.05);
//dummyHist->GetXaxis()->SetNdivisions(504);
pPhiCFc->SetTitle("p-K#rightarrow p-#phi");
pPhiCFc->Draw("same");
//smearedfitCF->Draw("L3same");
smearedfitCFpol->Draw("L3same");
//smearedCF->Draw("L3same");

//fitfun2->Draw("L3same");
auto totPeak2leg = new TLegend(0.65, 0.2, 0.65 + 0.25, 0.45);
totPeak2leg->SetBorderSize(0);
totPeak2leg->SetTextFont(42);
totPeak2leg->SetTextSize(gStyle->GetTextSize() * 0.7);
totPeak2leg->SetFillStyle(0);
totPeak2leg->AddEntry(pPhiCFc,"p-#phi", "pez");
totPeak2leg->AddEntry(smearedfitCFpol,"pol5","l");
//gle1g55->AddEntry(fitfun2,"Fit function", "l");
totPeak2leg->Draw("same");
totPeak2->Print("TOTsmearingPeak2.pdf");
totPeak2->Print("TOTsmearingPeak2.png");

auto totPeak3 = new TCanvas("TOTsmearingPeak3", "TOTsmearingPeak3", 0, 0, 650, 550);
totPeak3->SetRightMargin(right);
totPeak3->SetTopMargin(top);

dummyHist->Draw();
//dummyHist->GetYaxis()->SetRangeUser(0.85, 1.05);
//dummyHist->GetXaxis()->SetNdivisions(504);
pPhiCFc->SetTitle("p-K#rightarrow p-#phi");
pPhiCFc->Draw("same");
//smearedfitCF->Draw("L3same");
//smearedfitCFpol->Draw("L3same");
smearedCF->Draw("L3same");

//fitfun2->Draw("L3same");
auto totPeak3leg = new TLegend(0.65, 0.2, 0.65 + 0.25, 0.45);
totPeak3leg->SetBorderSize(0);
totPeak3leg->SetTextFont(42);
totPeak3leg->SetTextSize(gStyle->GetTextSize() * 0.7);
totPeak3leg->SetFillStyle(0);
totPeak3leg->AddEntry(pPhiCFc,"p-#phi", "pez");
totPeak3leg->AddEntry(smearedCF,"smeared data","l");
//gle1g55->AddEntry(fitfun2,"Fit function", "l");
totPeak3leg->Draw("same");
totPeak3->Print("TOTsmearingPeak3.pdf");
totPeak3->Print("TOTsmearingPeak3.png");




auto totPeak4 = new TCanvas("TOTsmearingPeak4", "TOTsmearingPeak4", 0, 0, 650, 550);
totPeak4->SetRightMargin(right);
totPeak4->SetTopMargin(top);

dummyHist->Draw();
//dummyHist->GetYaxis()->SetRangeUser(0.85, 1.05);
//dummyHist->GetXaxis()->SetNdivisions(504);
pPhiCFc->SetTitle("p-K#rightarrow#phi");
pPhiCFc->Draw("same");
smearedfitCF->Draw("L3same");
smearedfitCFpol->Draw("L3same");
smearedCF->Draw("L3same");
//normCF->Draw("same");
//fitfun2->Draw("L3same");
auto totPeak4leg = new TLegend(0.65, 0.2, 0.65 + 0.25, 0.45);
totPeak4leg->SetBorderSize(0);
totPeak4leg->SetTextFont(42);
totPeak4leg->SetTextSize(gStyle->GetTextSize() * 0.7);
totPeak4leg->SetFillStyle(0);
totPeak4leg->AddEntry(pPhiCFc,"p-#phi", "pez");
totPeak4leg->AddEntry(smearedfitCF,"gaus+pol2","l");
totPeak4leg->AddEntry(smearedfitCFpol,"pol5","l");
totPeak4leg->AddEntry(smearedCF,"smeared data","l");
//gle1g55->AddEntry(fitfun2,"Fit function", "l");
totPeak4leg->Draw("same");
totPeak4->Print("TOTsmearingPeak4.pdf");
totPeak4->Print("TOTsmearingPeak4.png");

DreamPlot::SetStyleGraph(smearedfitCF, 20, kGray, 0.7);
smearedfitCF->SetLineWidth(lineWidth);

auto totPeakmove = new TCanvas("TOTsmearingPeakmove", "TOTsmearingPeakmove", 0, 0, 650, 550);
totPeakmove->SetRightMargin(right);
totPeakmove->SetTopMargin(top);

dummyHist->Draw();
//dummyHist->GetYaxis()->SetRangeUser(0.85, 1.05);
//dummyHist->GetXaxis()->SetNdivisions(504);
pPhiCFc->SetTitle("p-K#rightarrow#phi");
pPhiCFc->Draw("same");
smearedfitCF->Draw("L3same");
movedCF->Draw("L3same");

//fitfun2->Draw("L3same");
auto totPeakmoveleg = new TLegend(0.65, 0.2, 0.65 + 0.25, 0.45);
totPeakmoveleg->SetBorderSize(0);
totPeakmoveleg->SetTextFont(42);
totPeakmoveleg->SetTextSize(gStyle->GetTextSize() * 0.7);
totPeakmoveleg->SetFillStyle(0);
totPeakmoveleg->AddEntry(pPhiCFc,"p-#phi", "pez");
totPeakmoveleg->AddEntry(smearedfitCF,"gaus+pol2","l");
totPeakmoveleg->AddEntry(movedCF,"moved CF","l");
//totPeakmoveleg->AddEntry(smearedCF,"smeared data","l");
//gle1g55->AddEntry(fitfun2,"Fit function", "l");
totPeakmoveleg->Draw("same");
totPeakmove->Print("TOTsmearingPeakmove.pdf");
totPeakmove->Print("TOTsmearingPeakmove.png");

DreamPlot::SetStyleGraph(smearedCF, 20, kGray, 0.7);
smearedCF->SetLineWidth(lineWidth);

auto totPeakmove2 = new TCanvas("TOTsmearingPeakmove2", "TOTsmearingPeakmove2", 0, 0, 650, 550);
totPeakmove2->SetRightMargin(right);
totPeakmove2->SetTopMargin(top);

dummyHist->Draw();
//dummyHist->GetYaxis()->SetRangeUser(0.85, 1.05);
//dummyHist->GetXaxis()->SetNdivisions(504);
pPhiCFc->SetTitle("p-K#rightarrow#phi");
pPhiCFc->Draw("same");
smearedCF->Draw("L3same");
movedCF2->Draw("L3same");

//fitfun2->Draw("L3same");
auto totPeakmoveleg2 = new TLegend(0.65, 0.2, 0.65 + 0.25, 0.45);
totPeakmoveleg2->SetBorderSize(0);
totPeakmoveleg2->SetTextFont(42);
totPeakmoveleg2->SetTextSize(gStyle->GetTextSize() * 0.7);
totPeakmoveleg2->SetFillStyle(0);
totPeakmoveleg2->AddEntry(pPhiCFc,"p-#phi", "pez");
totPeakmoveleg2->AddEntry(smearedCF,"smeared CF","l");
totPeakmoveleg2->AddEntry(movedCF2,"moved CF","l");
//totPeakmoveleg->AddEntry(smearedCF,"smeared data","l");
//gle1g55->AddEntry(fitfun2,"Fit function", "l");
totPeakmoveleg2->Draw("same");
totPeakmove2->Print("TOTsmearingPeakmove2.pdf");
totPeakmove2->Print("TOTsmearingPeakmove2.png");





auto totPeaknorm = new TCanvas("TOTsmearingPeaknorm", "TOTsmearingPeaknorm", 0, 0, 650, 550);
totPeaknorm->SetRightMargin(right);
totPeaknorm->SetTopMargin(top);

dummyHist->Draw();
//dummyHist->GetYaxis()->SetRangeUser(0.85, 1.05);
//dummyHist->GetXaxis()->SetNdivisions(504);
pPhiCFc->SetTitle("p-K#rightarrow#phi");
pPhiCFc->Draw("same");
smearedfitCF->Draw("L3same");
normCF->Draw("L3same");

//fitfun2->Draw("L3same");
auto totPeaknormleg = new TLegend(0.65, 0.2, 0.65 + 0.25, 0.45);
totPeaknormleg->SetBorderSize(0);
totPeaknormleg->SetTextFont(42);
totPeaknormleg->SetTextSize(gStyle->GetTextSize() * 0.7);
totPeaknormleg->SetFillStyle(0);
totPeaknormleg->AddEntry(pPhiCFc,"p-#phi", "pez");
totPeaknormleg->AddEntry(smearedfitCF,"gaus+pol2","l");
totPeaknormleg->AddEntry(normCF,"norm CF","l");
//totPeaknormleg->AddEntry(smearedCF,"smeared data","l");
//gle1g55->AddEntry(fitfun2,"Fit function", "l");
totPeaknormleg->Draw("same");
totPeaknorm->Print("TOTsmearingPeaknorm.pdf");
totPeaknorm->Print("TOTsmearingPeaknorm.png");

DreamPlot::SetStyleGraph(smearedCF, 20, kGray, 0.7);
smearedCF->SetLineWidth(lineWidth);

auto totPeaknorm2 = new TCanvas("TOTsmearingPeaknorm2", "TOTsmearingPeaknorm2", 0, 0, 650, 550);
totPeaknorm2->SetRightMargin(right);
totPeaknorm2->SetTopMargin(top);

dummyHist->Draw();
//dummyHist->GetYaxis()->SetRangeUser(0.85, 1.05);
//dummyHist->GetXaxis()->SetNdivisions(504);
pPhiCFc->SetTitle("p-K#rightarrow#phi");
pPhiCFc->Draw("same");
smearedCF->Draw("L3same");
normCF2->Draw("L3same");

//fitfun2->Draw("L3same");
auto totPeaknormleg2 = new TLegend(0.65, 0.2, 0.65 + 0.25, 0.45);
totPeaknormleg2->SetBorderSize(0);
totPeaknormleg2->SetTextFont(42);
totPeaknormleg2->SetTextSize(gStyle->GetTextSize() * 0.7);
totPeaknormleg2->SetFillStyle(0);
totPeaknormleg2->AddEntry(pPhiCFc,"p-#phi", "pez");
totPeaknormleg2->AddEntry(smearedCF,"smeared CF","l");
totPeaknormleg2->AddEntry(normCF2,"norm CF","l");
//totPeaknormleg->AddEntry(smearedCF,"smeared data","l");
//gle1g55->AddEntry(fitfun2,"Fit function", "l");
totPeaknormleg2->Draw("same");
totPeaknorm2->Print("TOTsmearingPeaknorm2.pdf");
totPeaknorm2->Print("TOTsmearingPeaknorm2.png");





auto totleft = new TCanvas("TOTsmearingSBleft", "TOTsmearingSBleft", 0, 0, 650, 550);
totleft->SetRightMargin(right);
totleft->SetTopMargin(top);

dummyHist->Draw();
//dummyHist->GetYaxis()->SetRangeUser(0.85, 1.05);
//dummyHist->GetXaxis()->SetNdivisions(504);
SBleftc->SetTitle("p-K#rightarrow p-(K^{+}K^{-}) SB left 0.987<M_{K^{+}K^{-}}<1.011");
SBleftc->Draw("same");
smearedfitSBleft->Draw("L3same");
smearedfitSBleftpol->Draw("L3same");
smearedSBleft->Draw("L3same");

//fitfun2->Draw("L3same");
auto totleftleg = new TLegend(0.65, 0.2, 0.65 + 0.25, 0.45);
totleftleg->SetBorderSize(0);
totleftleg->SetTextFont(42);
totleftleg->SetTextSize(gStyle->GetTextSize() * 0.7);
totleftleg->SetFillStyle(0);
totleftleg->AddEntry(SBleftc,"SB left", "pez");
totleftleg->AddEntry(smearedfitSBleft,"gaus+pol2","l");
totleftleg->AddEntry(smearedfitSBleftpol,"pol5","l");
totleftleg->AddEntry(smearedSBleft,"smeared data","l");
//gle1g55->AddEntry(fitfun2,"Fit function", "l");
totleftleg->Draw("same");
totleft->Print("TOTsmearingSBleft.pdf");
totleft->Print("TOTsmearingSBleft.png");

auto totleft2 = new TCanvas("TOTsmearingSBleft2", "TOTsmearingSBleft2", 0, 0, 650, 550);
totleft2->SetRightMargin(right);
totleft2->SetTopMargin(top);

dummyHist->Draw();
//dummyHist->GetYaxis()->SetRangeUser(0.85, 1.05);
//dummyHist->GetXaxis()->SetNdivisions(504);
SBleftc->SetTitle("p-K#rightarrow p-(K^{+}K^{-}) SB left 0.987<M_{K^{+}K^{-}}<1.011");
SBleftc->Draw("same");
smearedfitSBleft->Draw("L3same");
//smearedfitSBleftpol->Draw("L3same");
//smearedSBleft->Draw("L3same");

//fitfun2->Draw("L3same");
auto totleftleg2 = new TLegend(0.65, 0.2, 0.65 + 0.25, 0.45);
totleftleg2->SetBorderSize(0);
totleftleg2->SetTextFont(42);
totleftleg2->SetTextSize(gStyle->GetTextSize() * 0.7);
totleftleg2->SetFillStyle(0);
totleftleg2->AddEntry(SBleftc,"SB left", "pez");
totleftleg2->AddEntry(smearedfitSBleft,"gaus+pol2","l");
//gle1g55->AddEntry(fitfun2,"Fit function", "l");
totleftleg2->Draw("same");
totleft2->Print("TOTsmearingSBleft2.pdf");
totleft2->Print("TOTsmearingSBleft2.png");

auto totleft3 = new TCanvas("TOTsmearingSBleft", "TOTsmearingSBleft", 0, 0, 650, 550);
totleft3->SetRightMargin(right);
totleft3->SetTopMargin(top);

dummyHist->Draw();
//dummyHist->GetYaxis()->SetRangeUser(0.85, 1.05);
//dummyHist->GetXaxis()->SetNdivisions(504);
SBleftc->SetTitle("p-K#rightarrow p-(K^{+}K^{-}) SB left 0.987<M_{K^{+}K^{-}}<1.011");
SBleftc->Draw("same");
//smearedfitSBleft->Draw("L3same");
smearedfitSBleftpol->Draw("L3same");
//smearedSBleft->Draw("L3same");

//fitfun2->Draw("L3same");
auto totleftleg3 = new TLegend(0.65, 0.2, 0.65 + 0.25, 0.45);
totleftleg3->SetBorderSize(0);
totleftleg3->SetTextFont(42);
totleftleg3->SetTextSize(gStyle->GetTextSize() * 0.7);
totleftleg3->SetFillStyle(0);
totleftleg3->AddEntry(SBleftc,"SB left", "pez");
totleftleg3->AddEntry(smearedfitSBleftpol,"pol5","l");
//gle1g55->AddEntry(fitfun2,"Fit function", "l");
totleftleg3->Draw("same");
totleft3->Print("TOTsmearingSBleft3.pdf");
totleft3->Print("TOTsmearingSBleft3.png");

auto totleft4 = new TCanvas("TOTsmearingSBleft", "TOTsmearingSBleft", 0, 0, 650, 550);
totleft4->SetRightMargin(right);
totleft4->SetTopMargin(top);

dummyHist->Draw();
//dummyHist->GetYaxis()->SetRangeUser(0.85, 1.05);
//dummyHist->GetXaxis()->SetNdivisions(504);
SBleftc->SetTitle("p-K#rightarrow p-(K^{+}K^{-}) SB left 0.987<M_{K^{+}K^{-}}<1.011");
SBleftc->Draw("same");
//smearedfitSBleft->Draw("L3same");
//smearedfitSBleftpol->Draw("L3same");
smearedSBleft->Draw("L3same");

//fitfun2->Draw("L3same");
auto totleftleg4 = new TLegend(0.65, 0.2, 0.65 + 0.25, 0.45);
totleftleg4->SetBorderSize(0);
totleftleg4->SetTextFont(42);
totleftleg4->SetTextSize(gStyle->GetTextSize() * 0.7);
totleftleg4->SetFillStyle(0);
totleftleg4->AddEntry(SBleftc,"SB left", "pez");
totleftleg4->AddEntry(smearedSBleft,"smeared data","l");
//gle1g55->AddEntry(fitfun2,"Fit function", "l");
totleftleg4->Draw("same");
totleft4->Print("TOTsmearingSBleft4.pdf");
totleft4->Print("TOTsmearingSBleft4.png");

DreamPlot::SetStyleGraph(smearedfitSBleft, 20, kGray, 0.7);
smearedfitSBleft->SetLineWidth(lineWidth);

auto totleftmove = new TCanvas("TOTsmearingSBleftmove", "TOTsmearingSBleftmove", 0, 0, 650, 550);
totleftmove->SetRightMargin(right);
totleftmove->SetTopMargin(top);

dummyHist->Draw();
//dummyHist->GetYaxis()->SetRangeUser(0.85, 1.05);
//dummyHist->GetXaxis()->SetNdivisions(504);
SBleftc->SetTitle("p-K#rightarrow p-(K^{+}K^{-}) SB left 0.987<M_{K^{+}K^{-}}<1.011");
SBleftc->Draw("same");
smearedfitSBleft->Draw("L3same");
movedleft->Draw("L3same");
//smearedSBleft->Draw("L3same");

//fitfun2->Draw("L3same");
auto totleftlegmove = new TLegend(0.65, 0.2, 0.65 + 0.25, 0.45);
totleftlegmove->SetBorderSize(0);
totleftlegmove->SetTextFont(42);
totleftlegmove->SetTextSize(gStyle->GetTextSize() * 0.7);
totleftlegmove->SetFillStyle(0);
totleftlegmove->AddEntry(SBleftc,"SB left", "pez");
totleftlegmove->AddEntry(smearedfitSBleft,"gaus+pol2","l");
totleftlegmove->AddEntry(movedleft,"moved CF","l");

//gle1g55->AddEntry(fitfun2,"Fit function", "l");
totleftlegmove->Draw("same");
totleftmove->Print("TOTsmearingSBleftmove.pdf");
totleftmove->Print("TOTsmearingSBleftmove.png");

DreamPlot::SetStyleGraph(smearedSBleft, 20, kGray, 0.7);
smearedSBleft->SetLineWidth(lineWidth);

auto totleftmove2 = new TCanvas("TOTsmearingSBleftmove2", "TOTsmearingSBleftmove2", 0, 0, 650, 550);
totleftmove2->SetRightMargin(right);
totleftmove2->SetTopMargin(top);

dummyHist->Draw();
//dummyHist->GetYaxis()->SetRangeUser(0.85, 1.05);
//dummyHist->GetXaxis()->SetNdivisions(504);
SBleftc->SetTitle("p-K#rightarrow p-(K^{+}K^{-}) SB left 0.987<M_{K^{+}K^{-}}<1.011");
SBleftc->Draw("same");
smearedSBleft->Draw("L3same");
movedleft2->Draw("L3same");
//smearedSBleft->Draw("L3same");

//fitfun2->Draw("L3same");
auto totleftlegmove2 = new TLegend(0.65, 0.2, 0.65 + 0.25, 0.45);
totleftlegmove2->SetBorderSize(0);
totleftlegmove2->SetTextFont(42);
totleftlegmove2->SetTextSize(gStyle->GetTextSize() * 0.7);
totleftlegmove2->SetFillStyle(0);
totleftlegmove2->AddEntry(SBleftc,"SB left", "pez");
totleftlegmove2->AddEntry(smearedSBleft,"smeared CF","l");
totleftlegmove2->AddEntry(movedleft2,"moved CF","l");

//gle1g55->AddEntry(fitfun2,"Fit function", "l");
totleftlegmove2->Draw("same");
totleftmove2->Print("TOTsmearingSBleftmove2.pdf");
totleftmove2->Print("TOTsmearingSBleftmove2.png");




auto totleftnorm = new TCanvas("TOTsmearingSBleftnorm", "TOTsmearingSBleftnorm", 0, 0, 650, 550);
totleftnorm->SetRightMargin(right);
totleftnorm->SetTopMargin(top);

dummyHist->Draw();
//dummyHist->GetYaxis()->SetRangeUser(0.85, 1.05);
//dummyHist->GetXaxis()->SetNdivisions(504);
SBleftc->SetTitle("p-K#rightarrow p-(K^{+}K^{-}) SB left 0.987<M_{K^{+}K^{-}}<1.011");
SBleftc->Draw("same");
smearedfitSBleft->Draw("L3same");
normSBleft->Draw("L3same");
//smearedSBleft->Draw("L3same");

//fitfun2->Draw("L3same");
auto totleftlegnorm = new TLegend(0.65, 0.2, 0.65 + 0.25, 0.45);
totleftlegnorm->SetBorderSize(0);
totleftlegnorm->SetTextFont(42);
totleftlegnorm->SetTextSize(gStyle->GetTextSize() * 0.7);
totleftlegnorm->SetFillStyle(0);
totleftlegnorm->AddEntry(SBleftc,"SB left", "pez");
totleftlegnorm->AddEntry(smearedfitSBleft,"gaus+pol2","l");
totleftlegnorm->AddEntry(normSBleft,"normd CF","l");

//gle1g55->AddEntry(fitfun2,"Fit function", "l");
totleftlegnorm->Draw("same");
totleftnorm->Print("TOTsmearingSBleftnorm.pdf");
totleftnorm->Print("TOTsmearingSBleftnorm.png");

DreamPlot::SetStyleGraph(smearedSBleft, 20, kGray, 0.7);
smearedSBleft->SetLineWidth(lineWidth);

auto totleftnorm2 = new TCanvas("TOTsmearingSBleftnorm2", "TOTsmearingSBleftnorm2", 0, 0, 650, 550);
totleftnorm2->SetRightMargin(right);
totleftnorm2->SetTopMargin(top);

dummyHist->Draw();
//dummyHist->GetYaxis()->SetRangeUser(0.85, 1.05);
//dummyHist->GetXaxis()->SetNdivisions(504);
SBleftc->SetTitle("p-K#rightarrow p-(K^{+}K^{-}) SB left 0.987<M_{K^{+}K^{-}}<1.011");
SBleftc->Draw("same");
smearedSBleft->Draw("L3same");
normSBleft2->Draw("L3same");
//smearedSBleft->Draw("L3same");

//fitfun2->Draw("L3same");
auto totleftlegnorm2 = new TLegend(0.65, 0.2, 0.65 + 0.25, 0.45);
totleftlegnorm2->SetBorderSize(0);
totleftlegnorm2->SetTextFont(42);
totleftlegnorm2->SetTextSize(gStyle->GetTextSize() * 0.7);
totleftlegnorm2->SetFillStyle(0);
totleftlegnorm2->AddEntry(SBleftc,"SB left", "pez");
totleftlegnorm2->AddEntry(smearedSBleft,"smeared CF","l");
totleftlegnorm2->AddEntry(normSBleft2,"normd CF","l");

//gle1g55->AddEntry(fitfun2,"Fit function", "l");
totleftlegnorm2->Draw("same");
totleftnorm2->Print("TOTsmearingSBleftnorm2.pdf");
totleftnorm2->Print("TOTsmearingSBleftnorm2.png");














auto totright = new TCanvas("TOTsmearingSBright", "TOTsmearingSBright", 0, 0, 650, 550);
totright->SetRightMargin(right);
totright->SetTopMargin(top);

dummyHist->Draw();
//dummyHist->GetYaxis()->SetRangeUser(0.85, 1.05);
//dummyHist->GetXaxis()->SetNdivisions(504);
SBrightc->SetTitle("p-K#rightarrow p-(K^{+}K^{-}) SB right 1.027<M_{K^{+}K^{-}}<1.1");
SBrightc->Draw("same");
smearedfitSBright->Draw("L3same");
smearedfitSBrightpol->Draw("L3same");
smearedSBright->Draw("L3same");

//fitfun2->Draw("L3same");
auto totrightleg = new TLegend(0.65, 0.2, 0.65 + 0.25, 0.45);
totrightleg->SetBorderSize(0);
totrightleg->SetTextFont(42);
totrightleg->SetTextSize(gStyle->GetTextSize() * 0.7);
totrightleg->SetFillStyle(0);
totrightleg->AddEntry(SBrightc,"SB right", "pez");
totrightleg->AddEntry(smearedfitSBright,"gaus+pol2","l");
totrightleg->AddEntry(smearedfitSBrightpol,"pol5","l");
totrightleg->AddEntry(smearedSBright,"smeared data","l");
//gle1g55->AddEntry(fitfun2,"Fit function", "l");
totrightleg->Draw("same");
totright->Print("TOTsmearingSBright.pdf");
totright->Print("TOTsmearingSBright.png");

auto totright2 = new TCanvas("TOTsmearingSBright", "TOTsmearingSBright", 0, 0, 650, 550);
totright2->SetRightMargin(right);
totright2->SetTopMargin(top);

dummyHist->Draw();
//dummyHist->GetYaxis()->SetRangeUser(0.85, 1.05);
//dummyHist->GetXaxis()->SetNdivisions(504);
SBrightc->SetTitle("p-K#rightarrow p-(K^{+}K^{-}) SB right 1.027<M_{K^{+}K^{-}}<1.1");
SBrightc->Draw("same");
smearedfitSBright->Draw("L3same");
//smearedfitSBrightpol->Draw("L3same");
//smearedSBright->Draw("L3same");

//fitfun2->Draw("L3same");
auto totrightleg2 = new TLegend(0.65, 0.2, 0.65 + 0.25, 0.45);
totrightleg2->SetBorderSize(0);
totrightleg2->SetTextFont(42);
totrightleg2->SetTextSize(gStyle->GetTextSize() * 0.7);
totrightleg2->SetFillStyle(0);
totrightleg2->AddEntry(SBrightc,"SB right", "pez");
totrightleg2->AddEntry(smearedfitSBright,"gaus+pol2","l");
//gle1g55->AddEntry(fitfun2,"Fit function", "l");
totrightleg2->Draw("same");
totright2->Print("TOTsmearingSBright2.pdf");
totright2->Print("TOTsmearingSBright2.png");


auto totright3 = new TCanvas("TOTsmearingSBright", "TOTsmearingSBright", 0, 0, 650, 550);
totright3->SetRightMargin(right);
totright3->SetTopMargin(top);

dummyHist->Draw();
//dummyHist->GetYaxis()->SetRangeUser(0.85, 1.05);
//dummyHist->GetXaxis()->SetNdivisions(504);
SBrightc->SetTitle("p-K#rightarrow p-(K^{+}K^{-}) SB right 1.027<M_{K^{+}K^{-}}<1.1");
SBrightc->Draw("same");
//smearedfitSBright->Draw("L3same");
smearedfitSBrightpol->Draw("L3same");
//smearedSBright->Draw("L3same");

//fitfun2->Draw("L3same");
auto totrightleg3 = new TLegend(0.65, 0.2, 0.65 + 0.25, 0.45);
totrightleg3->SetBorderSize(0);
totrightleg3->SetTextFont(42);
totrightleg3->SetTextSize(gStyle->GetTextSize() * 0.7);
totrightleg3->SetFillStyle(0);
totrightleg3->AddEntry(SBrightc,"SB right", "pez");
totrightleg3->AddEntry(smearedfitSBrightpol,"pol5","l");
//gle1g55->AddEntry(fitfun2,"Fit function", "l");
totrightleg3->Draw("same");
totright3->Print("TOTsmearingSBright3.pdf");
totright3->Print("TOTsmearingSBright3.png");


auto totright4 = new TCanvas("TOTsmearingSBright", "TOTsmearingSBright", 0, 0, 650, 550);
totright4->SetRightMargin(right);
totright4->SetTopMargin(top);

dummyHist->Draw();
//dummyHist->GetYaxis()->SetRangeUser(0.85, 1.05);
//dummyHist->GetXaxis()->SetNdivisions(504);
SBrightc->SetTitle("p-K#rightarrow p-(K^{+}K^{-}) SB right 1.027<M_{K^{+}K^{-}}<1.1");
SBrightc->Draw("same");
//smearedfitSBright->Draw("L3same");
//smearedfitSBrightpol->Draw("L3same");
smearedSBright->Draw("L3same");

//fitfun2->Draw("L3same");
auto totrightleg4 = new TLegend(0.65, 0.2, 0.65 + 0.25, 0.45);
totrightleg4->SetBorderSize(0);
totrightleg4->SetTextFont(42);
totrightleg4->SetTextSize(gStyle->GetTextSize() * 0.7);
totrightleg4->SetFillStyle(0);
totrightleg4->AddEntry(SBrightc,"SB right", "pez");
totrightleg4->AddEntry(smearedSBright,"smeared data","l");
//gle1g55->AddEntry(fitfun2,"Fit function", "l");
totrightleg4->Draw("same");
totright4->Print("TOTsmearingSBright4.pdf");
totright4->Print("TOTsmearingSBright4.png");


DreamPlot::SetStyleGraph(smearedfitSBright, 20, kGray, 0.7);
smearedfitSBright->SetLineWidth(lineWidth);

auto totrightmove = new TCanvas("TOTsmearingSBrightmove", "TOTsmearingSBrightmove", 0, 0, 650, 550);
totrightmove->SetRightMargin(right);
totrightmove->SetTopMargin(top);

dummyHist->Draw();
//dummyHist->GetYaxis()->SetRangeUser(0.85, 1.05);
//dummyHist->GetXaxis()->SetNdivisions(504);
SBrightc->SetTitle("p-K#rightarrow p-(K^{+}K^{-}) SB right 1.027<M_{K^{+}K^{-}}<1.1");
SBrightc->Draw("same");
smearedfitSBright->Draw("L3same");
movedright->Draw("L3same");
//smearedSBright->Draw("L3same");

//fitfun2->Draw("L3same");
auto totrightlegmove = new TLegend(0.65, 0.2, 0.65 + 0.25, 0.45);
totrightlegmove->SetBorderSize(0);
totrightlegmove->SetTextFont(42);
totrightlegmove->SetTextSize(gStyle->GetTextSize() * 0.7);
totrightlegmove->SetFillStyle(0);
totrightlegmove->AddEntry(SBrightc,"SB right", "pez");
totrightlegmove->AddEntry(smearedfitSBright,"gaus+pol2","l");
totrightlegmove->AddEntry(movedright,"moved CF","l");

//gle1g55->AddEntry(fitfun2,"Fit function", "l");
totrightlegmove->Draw("same");
totrightmove->Print("TOTsmearingSBrightmove.pdf");
totrightmove->Print("TOTsmearingSBrightmove.png");


DreamPlot::SetStyleGraph(smearedSBright, 20, kGray, 0.7);
smearedSBright->SetLineWidth(lineWidth);

auto totrightmove2 = new TCanvas("TOTsmearingSBrightmove2", "TOTsmearingSBrightmove2", 0, 0, 650, 550);
totrightmove2->SetRightMargin(right);
totrightmove2->SetTopMargin(top);

dummyHist->Draw();
//dummyHist->GetYaxis()->SetRangeUser(0.85, 1.05);
//dummyHist->GetXaxis()->SetNdivisions(504);
SBrightc->SetTitle("p-K#rightarrow p-(K^{+}K^{-}) SB right 1.027<M_{K^{+}K^{-}}<1.1");
SBrightc->Draw("same");
smearedSBright->Draw("L3same");
movedright2->Draw("L3same");
//smearedSBright->Draw("L3same");

//fitfun2->Draw("L3same");
auto totrightlegmove2 = new TLegend(0.65, 0.2, 0.65 + 0.25, 0.45);
totrightlegmove2->SetBorderSize(0);
totrightlegmove2->SetTextFont(42);
totrightlegmove2->SetTextSize(gStyle->GetTextSize() * 0.7);
totrightlegmove2->SetFillStyle(0);
totrightlegmove2->AddEntry(SBrightc,"SB right", "pez");
totrightlegmove2->AddEntry(smearedSBright,"smeared CF","l");
totrightlegmove2->AddEntry(movedright2,"moved CF","l");

//gle1g55->AddEntry(fitfun2,"Fit function", "l");
totrightlegmove2->Draw("same");
totrightmove2->Print("TOTsmearingSBrightmove2.pdf");
totrightmove2->Print("TOTsmearingSBrightmove2.png");




auto totrightnorm = new TCanvas("TOTsmearingSBrightnorm", "TOTsmearingSBrightnorm", 0, 0, 650, 550);
totrightnorm->SetRightMargin(right);
totrightnorm->SetTopMargin(top);

dummyHist->Draw();
//dummyHist->GetYaxis()->SetRangeUser(0.85, 1.05);
//dummyHist->GetXaxis()->SetNdivisions(504);
SBrightc->SetTitle("p-K#rightarrow p-(K^{+}K^{-}) SB right 1.027<M_{K^{+}K^{-}}<1.1");
SBrightc->Draw("same");
smearedfitSBright->Draw("L3same");
normSBright->Draw("L3same");
//smearedSBright->Draw("L3same");

//fitfun2->Draw("L3same");
auto totrightlegnorm = new TLegend(0.65, 0.2, 0.65 + 0.25, 0.45);
totrightlegnorm->SetBorderSize(0);
totrightlegnorm->SetTextFont(42);
totrightlegnorm->SetTextSize(gStyle->GetTextSize() * 0.7);
totrightlegnorm->SetFillStyle(0);
totrightlegnorm->AddEntry(SBrightc,"SB right", "pez");
totrightlegnorm->AddEntry(smearedfitSBright,"gaus+pol2","l");
totrightlegnorm->AddEntry(normSBright,"norm CF","l");

//gle1g55->AddEntry(fitfun2,"Fit function", "l");
totrightlegnorm->Draw("same");
totrightnorm->Print("TOTsmearingSBrightnorm.pdf");
totrightnorm->Print("TOTsmearingSBrightnorm.png");


DreamPlot::SetStyleGraph(smearedSBright, 20, kGray, 0.7);
smearedSBright->SetLineWidth(lineWidth);

auto totrightnorm2 = new TCanvas("TOTsmearingSBrightnorm2", "TOTsmearingSBrightnorm2", 0, 0, 650, 550);
totrightnorm2->SetRightMargin(right);
totrightnorm2->SetTopMargin(top);

dummyHist->Draw();
//dummyHist->GetYaxis()->SetRangeUser(0.85, 1.05);
//dummyHist->GetXaxis()->SetNdivisions(504);
SBrightc->SetTitle("p-K#rightarrow p-(K^{+}K^{-}) SB right 1.027<M_{K^{+}K^{-}}<1.1");
SBrightc->Draw("same");
smearedSBright->Draw("L3same");
normSBright2->Draw("L3same");
//smearedSBright->Draw("L3same");

//fitfun2->Draw("L3same");
auto totrightlegnorm2 = new TLegend(0.65, 0.2, 0.65 + 0.25, 0.45);
totrightlegnorm2->SetBorderSize(0);
totrightlegnorm2->SetTextFont(42);
totrightlegnorm2->SetTextSize(gStyle->GetTextSize() * 0.7);
totrightlegnorm2->SetFillStyle(0);
totrightlegnorm2->AddEntry(SBrightc,"SB right", "pez");
totrightlegnorm2->AddEntry(smearedSBright,"smeared CF","l");
totrightlegnorm2->AddEntry(normSBright2,"norm CF","l");

//gle1g55->AddEntry(fitfun2,"Fit function", "l");
totrightlegnorm2->Draw("same");
totrightnorm2->Print("TOTsmearingSBrightnorm2.pdf");
totrightnorm2->Print("TOTsmearingSBrightnorm2.png");










dummyHist->GetYaxis()->SetRangeUser(0.9, 1.15);


auto tot12 = new TCanvas("TOTsmearingSB12", "TOTsmearingSB12", 0, 0, 650, 550);
tot12->SetRightMargin(right);
tot12->SetTopMargin(top);

dummyHist->Draw();
//dummyHist->GetYaxis()->SetRangeUser(0.85, 1.05);
//dummyHist->GetXaxis()->SetNdivisions(504);
SB12c->SetTitle("p-K#rightarrow p-(K^{+}K^{-}) SB 1.1<M_{K^{+}K^{-}}<1.2");
SB12c->Draw("same");
smearedfitSB12->Draw("L3same");
smearedfitSB12pol->Draw("L3same");
smearedSB12->Draw("L3same");

//fitfun2->Draw("L3same");
auto tot12leg = new TLegend(0.65, 0.2, 0.65 + 0.25, 0.45);
tot12leg->SetBorderSize(0);
tot12leg->SetTextFont(42);
tot12leg->SetTextSize(gStyle->GetTextSize() * 0.7);
tot12leg->SetFillStyle(0);
tot12leg->AddEntry(SB12c,"SB CF", "pez");
tot12leg->AddEntry(smearedfitSB12,"gaus+pol2","l");
tot12leg->AddEntry(smearedfitSB12pol,"pol5","l");
tot12leg->AddEntry(smearedSB12,"smeared data","l");
//gle1g55->AddEntry(fitfun2,"Fit function", "l");
tot12leg->Draw("same");
tot12->Print("TOTsmearingSB12.pdf");
tot12->Print("TOTsmearingSB12.png");




auto tot23 = new TCanvas("TOTsmearingSB23", "TOTsmearingSB23", 0, 0, 650, 550);
tot23->SetRightMargin(right);
tot23->SetTopMargin(top);

dummyHist->Draw();
//dummyHist->GetYaxis()->SetRangeUser(0.85, 1.05);
//dummyHist->GetXaxis()->SetNdivisions(504);
SB23c->SetTitle("p-K#rightarrow p-(K^{+}K^{-}) SB 1.2<M_{K^{+}K^{-}}<1.3");
SB23c->Draw("same");
smearedfitSB23->Draw("L3same");
smearedfitSB23pol->Draw("L3same");
smearedSB23->Draw("L3same");

//fitfun2->Draw("L3same");
auto tot23leg = new TLegend(0.65, 0.2, 0.65 + 0.25, 0.45);
tot23leg->SetBorderSize(0);
tot23leg->SetTextFont(42);
tot23leg->SetTextSize(gStyle->GetTextSize() * 0.7);
tot23leg->SetFillStyle(0);
tot23leg->AddEntry(SB23c,"SB CF", "pez");
tot23leg->AddEntry(smearedfitSB23,"gaus+pol2","l");
tot23leg->AddEntry(smearedfitSB23pol,"pol5","l");
tot23leg->AddEntry(smearedSB23,"smeared data","l");
//gle1g55->AddEntry(fitfun2,"Fit function", "l");
tot23leg->Draw("same");
tot23->Print("TOTsmearingSB23.pdf");
tot23->Print("TOTsmearingSB23.png");




auto tot34 = new TCanvas("TOTsmearingSB34", "TOTsmearingSB34", 0, 0, 650, 550);
tot34->SetRightMargin(right);
tot34->SetTopMargin(top);

dummyHist->Draw();
//dummyHist->GetYaxis()->SetRangeUser(0.85, 1.05);
//dummyHist->GetXaxis()->SetNdivisions(504);
SB34c->SetTitle("p-K#rightarrow p-(K^{+}K^{-}) SB 1.3<M_{K^{+}K^{-}}<1.4");
SB34c->Draw("same");
smearedfitSB34->Draw("L3same");
smearedfitSB34pol->Draw("L3same");
smearedSB34->Draw("L3same");

//fitfun2->Draw("L3same");
auto tot34leg = new TLegend(0.65, 0.2, 0.65 + 0.25, 0.45);
tot34leg->SetBorderSize(0);
tot34leg->SetTextFont(42);
tot34leg->SetTextSize(gStyle->GetTextSize() * 0.7);
tot34leg->SetFillStyle(0);
tot34leg->AddEntry(SB34c,"SB CF", "pez");
tot34leg->AddEntry(smearedfitSB34,"gaus+pol2","l");
tot34leg->AddEntry(smearedfitSB34pol,"pol5","l");
tot34leg->AddEntry(smearedSB34,"smeared data","l");
//gle1g55->AddEntry(fitfun2,"Fit function", "l");
tot34leg->Draw("same");
tot34->Print("TOTsmearingSB34.pdf");
tot34->Print("TOTsmearingSB34.png");



auto tot45 = new TCanvas("TOTsmearingSB45", "TOTsmearingSB45", 0, 0, 650, 550);
tot45->SetRightMargin(right);
tot45->SetTopMargin(top);

dummyHist->Draw();
//dummyHist->GetYaxis()->SetRangeUser(0.85, 1.05);
//dummyHist->GetXaxis()->SetNdivisions(504);
SB45c->SetTitle("p-K#rightarrow p-(K^{+}K^{-}) SB 1.4<M_{K^{+}K^{-}}<1.5");
SB45c->Draw("same");
smearedfitSB45->Draw("L3same");
smearedfitSB45pol->Draw("L3same");
smearedSB45->Draw("L3same");

//fitfun2->Draw("L3same");
auto tot45leg = new TLegend(0.65, 0.2, 0.65 + 0.25, 0.45);
tot45leg->SetBorderSize(0);
tot45leg->SetTextFont(42);
tot45leg->SetTextSize(gStyle->GetTextSize() * 0.7);
tot45leg->SetFillStyle(0);
tot45leg->AddEntry(SB45c,"SB CF", "pez");
tot45leg->AddEntry(smearedfitSB45,"gaus+pol2","l");
tot45leg->AddEntry(smearedfitSB45pol,"pol5","l");
tot45leg->AddEntry(smearedSB45,"smeared data","l");
//gle1g55->AddEntry(fitfun2,"Fit function", "l");
tot45leg->Draw("same");
tot45->Print("TOTsmearingSB45.pdf");
tot45->Print("TOTsmearingSB45.png");



auto tot56 = new TCanvas("TOTsmearingSB56", "TOTsmearingSB56", 0, 0, 650, 550);
tot56->SetRightMargin(right);
tot56->SetTopMargin(top);

dummyHist->Draw();
//dummyHist->GetYaxis()->SetRangeUser(0.85, 1.05);
//dummyHist->GetXaxis()->SetNdivisions(504);
SB56c->SetTitle("p-K#rightarrow p-(K^{+}K^{-}) SB 1.5<M_{K^{+}K^{-}}<1.6");
SB56c->Draw("same");
smearedfitSB56->Draw("L3same");
smearedfitSB56pol->Draw("L3same");
smearedSB56->Draw("L3same");

//fitfun2->Draw("L3same");
auto tot56leg = new TLegend(0.65, 0.2, 0.65 + 0.25, 0.45);
tot56leg->SetBorderSize(0);
tot56leg->SetTextFont(42);
tot56leg->SetTextSize(gStyle->GetTextSize() * 0.7);
tot56leg->SetFillStyle(0);
tot56leg->AddEntry(SB56c,"SB CF", "pez");
tot56leg->AddEntry(smearedfitSB56,"gaus+pol2","l");
tot56leg->AddEntry(smearedfitSB56pol,"pol5","l");
tot56leg->AddEntry(smearedSB56,"smeared data","l");
//gle1g55->AddEntry(fitfun2,"Fit function", "l");
tot56leg->Draw("same");
tot56->Print("TOTsmearingSB56.pdf");
tot56->Print("TOTsmearingSB56.png");



auto tot67 = new TCanvas("TOTsmearingSB67", "TOTsmearingSB67", 0, 0, 650, 550);
tot67->SetRightMargin(right);
tot67->SetTopMargin(top);

dummyHist->Draw();
//dummyHist->GetYaxis()->SetRangeUser(0.85, 1.05);
//dummyHist->GetXaxis()->SetNdivisions(504);
SB67c->SetTitle("p-K#rightarrow p-(K^{+}K^{-}) SB 1.6<M_{K^{+}K^{-}}<1.7");
SB67c->Draw("same");
smearedfitSB67->Draw("L3same");
smearedfitSB67pol->Draw("L3same");
smearedSB67->Draw("L3same");

//fitfun2->Draw("L3same");
auto tot67leg = new TLegend(0.65, 0.2, 0.65 + 0.25, 0.45);
tot67leg->SetBorderSize(0);
tot67leg->SetTextFont(42);
tot67leg->SetTextSize(gStyle->GetTextSize() * 0.7);
tot67leg->SetFillStyle(0);
tot67leg->AddEntry(SB67c,"SB CF", "pez");
tot67leg->AddEntry(smearedfitSB67,"gaus+pol2","l");
tot67leg->AddEntry(smearedfitSB67pol,"pol5","l");
tot67leg->AddEntry(smearedSB67,"smeared data","l");
//gle1g55->AddEntry(fitfun2,"Fit function", "l");
tot67leg->Draw("same");
tot67->Print("TOTsmearingSB67.pdf");
tot67->Print("TOTsmearingSB67.png");



auto tot78 = new TCanvas("TOTsmearingSB78", "TOTsmearingSB78", 0, 0, 650, 550);
tot78->SetRightMargin(right);
tot78->SetTopMargin(top);

dummyHist->Draw();
//dummyHist->GetYaxis()->SetRangeUser(0.85, 1.05);
//dummyHist->GetXaxis()->SetNdivisions(504);
SB78c->SetTitle("p-K#rightarrow p-(K^{+}K^{-}) SB 1.7<M_{K^{+}K^{-}}<1.8");
SB78c->Draw("same");
smearedfitSB78->Draw("L3same");
smearedfitSB78pol->Draw("L3same");
smearedSB78->Draw("L3same");

//fitfun2->Draw("L3same");
auto tot78leg = new TLegend(0.65, 0.2, 0.65 + 0.25, 0.45);
tot78leg->SetBorderSize(0);
tot78leg->SetTextFont(42);
tot78leg->SetTextSize(gStyle->GetTextSize() * 0.7);
tot78leg->SetFillStyle(0);
tot78leg->AddEntry(SB78c,"SB CF", "pez");
tot78leg->AddEntry(smearedfitSB78,"gaus+pol2","l");
tot78leg->AddEntry(smearedfitSB78pol,"pol5","l");
tot78leg->AddEntry(smearedSB78,"smeared data","l");
//gle1g55->AddEntry(fitfun2,"Fit function", "l");
tot78leg->Draw("same");
tot78->Print("TOTsmearingSB78.pdf");
tot78->Print("TOTsmearingSB78.png");



auto tot89 = new TCanvas("TOTsmearingSB89", "TOTsmearingSB89", 0, 0, 650, 550);
tot89->SetRightMargin(right);
tot89->SetTopMargin(top);

dummyHist->Draw();
//dummyHist->GetYaxis()->SetRangeUser(0.85, 1.05);
//dummyHist->GetXaxis()->SetNdivisions(504);
SB89c->SetTitle("p-K#rightarrow p-(K^{+}K^{-}) SB 1.8<M_{K^{+}K^{-}}<1.9");
SB89c->Draw("same");
smearedfitSB89->Draw("L3same");
smearedfitSB89pol->Draw("L3same");
smearedSB89->Draw("L3same");

//fitfun2->Draw("L3same");
auto tot89leg = new TLegend(0.65, 0.2, 0.65 + 0.25, 0.45);
tot89leg->SetBorderSize(0);
tot89leg->SetTextFont(42);
tot89leg->SetTextSize(gStyle->GetTextSize() * 0.7);
tot89leg->SetFillStyle(0);
tot89leg->AddEntry(SB89c,"SB CF", "pez");
tot89leg->AddEntry(smearedfitSB89,"gaus+pol2","l");
tot89leg->AddEntry(smearedfitSB89pol,"pol5","l");
tot89leg->AddEntry(smearedSB89,"smeared data","l");
//gle1g55->AddEntry(fitfun2,"Fit function", "l");
tot89leg->Draw("same");
tot89->Print("TOTsmearingSB89.pdf");
tot89->Print("TOTsmearingSB89.png");



auto tot90 = new TCanvas("TOTsmearingSB90", "TOTsmearingSB90", 0, 0, 650, 550);
tot90->SetRightMargin(right);
tot90->SetTopMargin(top);

dummyHist->Draw();
//dummyHist->GetYaxis()->SetRangeUser(0.85, 1.05);
//dummyHist->GetXaxis()->SetNdivisions(504);
SB90c->SetTitle("p-K#rightarrow p-(K^{+}K^{-}) SB 1.9<M_{K^{+}K^{-}}<2.0");
SB90c->Draw("same");
smearedfitSB90->Draw("L3same");
smearedfitSB90pol->Draw("L3same");
smearedSB90->Draw("L3same");

//fitfun2->Draw("L3same");
auto tot90leg = new TLegend(0.65, 0.2, 0.65 + 0.25, 0.45);
tot90leg->SetBorderSize(0);
tot90leg->SetTextFont(42);
tot90leg->SetTextSize(gStyle->GetTextSize() * 0.7);
tot90leg->SetFillStyle(0);
tot90leg->AddEntry(SB90c,"SB CF", "pez");
tot90leg->AddEntry(smearedfitSB90,"gaus+pol2","l");
tot90leg->AddEntry(smearedfitSB90pol,"pol5","l");
tot90leg->AddEntry(smearedSB90,"smeared data","l");
//gle1g55->AddEntry(fitfun2,"Fit function", "l");
tot90leg->Draw("same");
tot90->Print("TOTsmearingSB90.pdf");
tot90->Print("TOTsmearingSB90.png");








smearedfitCF->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedfitCF->SetMarkerColor(2);
smearedfitCF->SetMarkerSize(1.5);
smearedfitCF->SetMarkerStyle(21);

smearedfitSBleft->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedfitSBleft->SetMarkerColor(2);
smearedfitSBleft->SetMarkerSize(1.5);
smearedfitSBleft->SetMarkerStyle(21);

smearedfitSBright->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedfitSBright->SetMarkerColor(2);
smearedfitSBright->SetMarkerSize(1.5);
smearedfitSBright->SetMarkerStyle(21);

smearedfitSB12->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedfitSB12->SetMarkerColor(2);
smearedfitSB12->SetMarkerSize(1.5);
smearedfitSB12->SetMarkerStyle(21);

smearedfitSB23->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedfitSB23->SetMarkerColor(2);
smearedfitSB23->SetMarkerSize(1.5);
smearedfitSB23->SetMarkerStyle(21);

smearedfitSB34->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedfitSB34->SetMarkerColor(2);
smearedfitSB34->SetMarkerSize(1.5);
smearedfitSB34->SetMarkerStyle(21);


smearedfitSB45->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedfitSB45->SetMarkerColor(2);
smearedfitSB45->SetMarkerSize(1.5);
smearedfitSB45->SetMarkerStyle(21);


smearedfitSB56->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedfitSB56->SetMarkerColor(2);
smearedfitSB56->SetMarkerSize(1.5);
smearedfitSB56->SetMarkerStyle(21);


smearedfitSB67->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedfitSB67->SetMarkerColor(2);
smearedfitSB67->SetMarkerSize(1.5);
smearedfitSB67->SetMarkerStyle(21);


smearedfitSB78->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedfitSB78->SetMarkerColor(2);
smearedfitSB78->SetMarkerSize(1.5);
smearedfitSB78->SetMarkerStyle(21);


smearedfitSB89->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedfitSB89->SetMarkerColor(2);
smearedfitSB89->SetMarkerSize(1.5);
smearedfitSB89->SetMarkerStyle(21);


smearedfitSB90->SetLineColor(2);
//smearedCF->SetLineWidth(4);
smearedfitSB90->SetMarkerColor(2);
smearedfitSB90->SetMarkerSize(1.5);
smearedfitSB90->SetMarkerStyle(21);




  outfile->cd();

  histSmearTOTSBleft->Write("histSmearTOTSBleft");
  histSmearTOTPeak->Write("histSmearTOTPeak");
  histSmearTOTSBright->Write("histSmearTOTSBright");
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


  projection0->Write("projection0");


  pKCF->Write("pKCFgaus");
  pKCF66->Write("pKCFgausver2");

  pKBG->Write("pKBG");
  pKCF2->Write("pKCFpol");
  pKCFDimi->Write("pKCFpolDIMI");
  pKCFDimipol5->Write("pKCFpolDIMIpol5");
  grint->Write("confint");

  ConfintUP->Write("high");
  ConfintLOW->Write("low");

  smearedconfSBleftlow->Write("leftlow");
  smearedconfSBleftup->Write("leftup");
  smearedconfSBrightlow->Write("rightlow");
  smearedconfSBrightup->Write("rightup");
  normSBleftLOW->Write("leftlownorm");
  normSBleftUP->Write("leftupnorm");
  normSBrightLOW->Write("rightlownorm");
  normSBrightUP->Write("rightupnorm");



  smearedCF->Write("smearedCF");
  smearedSBleft->Write("smearedSBleft");
  smearedSBright->Write("smearedSBright");
  smearedSB12->Write("smearedSB12");
  smearedSB23->Write("smearedSB23");
  smearedSB34->Write("smearedSB34");
  smearedSB45->Write("smearedSB45");
  smearedSB56->Write("smearedSB56");
  smearedSB67->Write("smearedSB67");
  smearedSB78->Write("smearedSB78");
  smearedSB89->Write("smearedSB89");
  smearedSB90->Write("smearedSB90");




  //  smearedCFpol->Write("smearedCFpol");
  //  smearedSBleftpol->Write("smearedSBleftpol");
  //  smearedSBrightpol->Write("smearedSBrightpol");
  //  smearedSB12pol->Write("smearedSB12pol");
  //  smearedSB23pol->Write("smearedSB23pol");
  //  smearedSB34pol->Write("smearedSB34pol");
  //  smearedSB45pol->Write("smearedSB45pol");
  //  smearedSB56pol->Write("smearedSB56pol");
  //  smearedSB67pol->Write("smearedSB67pol");
  //  smearedSB78pol->Write("smearedSB78pol");
  //  smearedSB89pol->Write("smearedSB89pol");
  //  smearedSB90pol->Write("smearedSB90pol");

//  total->Write("totfit");
  fitfunpol->Write("fitfunpol");
  fitfungaus->Write("fitfungaus");

//  smearedfitCF->Write("smearedfitCF");
//  smearedfitSBleft->Write("smearedfitSBleft");
//  smearedfitSBright->Write("smearedfitSBright");
//  smearedfitSB12->Write("smearedfitSB12");
//  smearedfitSB23->Write("smearedfitSB23");
//  smearedfitSB34->Write("smearedfitSB34");


  smearedfitCF->Write("smearedfitCF");
  smearedfitSBleft->Write("smearedfitSBleft");
  smearedfitSBright->Write("smearedfitSBright");
  smearedfitSB12->Write("smearedfitSB12");
  smearedfitSB23->Write("smearedfitSB23");
  smearedfitSB34->Write("smearedfitSB34");
  smearedfitSB45->Write("smearedfitSB45");
  smearedfitSB56->Write("smearedfitSB56");
  smearedfitSB67->Write("smearedfitSB67");
  smearedfitSB78->Write("smearedfitSB78");
  smearedfitSB89->Write("smearedfitSB89");
  smearedfitSB90->Write("smearedfitSB90");


  smearedfitCFpol->Write("smearedfitCFpol");
  smearedfitSBleftpol->Write("smearedfitSBleftpol");
  smearedfitSBrightpol->Write("smearedfitSBrightpol");
  smearedfitSB12pol->Write("smearedfitSB12pol");
  smearedfitSB23pol->Write("smearedfitSB23pol");
  smearedfitSB34pol->Write("smearedfitSB34pol");
  smearedfitSB45pol->Write("smearedfitSB45pol");
  smearedfitSB56pol->Write("smearedfitSB56pol");
  smearedfitSB67pol->Write("smearedfitSB67pol");
  smearedfitSB78pol->Write("smearedfitSB78pol");
  smearedfitSB89pol->Write("smearedfitSB89pol");
  smearedfitSB90pol->Write("smearedfitSB90pol");



  movedCF->Write("movedCF");
  movedleft->Write("movedleft");
  movedright->Write("movedright");

  movedCF2->Write("movedCF2");
  movedleft2->Write("movedleft2");
  movedright2->Write("movedright2");


  normCF->Write("normCF");
  normSBleft->Write("normSBleft");
  normSBright->Write("normSBright");
//  moved12->Write("moved12");
//  moved23->Write("moved23");
//  moved34->Write("moved34");
//  moved45->Write("moved45");
//  moved56->Write("moved56");
//  moved67->Write("moved67");
//  moved78->Write("moved78");
//  moved89->Write("moved89");
//  moved90->Write("moved90");

  outfile->Close();

}

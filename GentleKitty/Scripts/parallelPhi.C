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
#include "TVirtualFitter.h"
#include "TList.h"
#include "TMinuit.h"
#include "TGraph2DErrors.h"
#include "TFitter.h"

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



double gausspol2(double x, double *par) {

    //cout<<par[0]<<"  "<<par[1]<<"  "<<par[2]<<"  "<<par[3]<<"  "<<par[4]<<"  "<<par[5]<<"  "<<endl;
    double y = par[0] * TMath::Gaus(x, par[1], par[2], false) + par[3] + par[4] * x +
            par[5] * x * x;
    return y;
}


// double integrate(TGraph* CF, double min, double max){
// double sum=0;
// int diff=max-min;
// int a=0;
// for (int i=min;i<=max;i++){
// double ycf= CF->Eval(i);
// cout<<"value: "<< i<< " ----> "<<ycf<<endl;
// sum = sum+ycf;
// a++;
//}
// cout<<sum<<endl;
// return sum/a;
// cout<<sum/a<<endl;

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
//      Form("/home/emma/FemtoPhiHM_checks2/800-900/CFOutput_pPhi_%s_%s.root",
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


TGraphAsymmErrors* GetpKeff(const char* path) {
  const char* name =
      Form("%s/SidebandInspect.root",
           path);
  auto file = TFile::Open(name);
  // file->ls();
  TGraphAsymmErrors* pK = (TGraphAsymmErrors*)file->FindObjectAny("pKCFgaus");
  //file->TFile::Close();
  return pK;
}


TH2F* GetMatix(const char* path, const char* filename) {
  const char* name =
      Form("%s/SmearSideband.root",
           path);
  auto file = TFile::Open(name);
  // file->ls();
  TH2F* matrix = (TH2F*)file->FindObjectAny(Form("%s", filename));
 // file->TFile::Close();
  return matrix;
}


TH1F* GetSB(const char* path, const char* filename) {
  const char* name =
      Form("%s/SidebandInspect.root",
           path);
  auto file = TFile::Open(name);
 // file->ls();
  TH1F* SB = (TH1F*)file->FindObjectAny(Form("%s", filename));
 // file->TFile::Close();
  return SB;
}

TH1F* GetSB() {
  const char* name ="/home/emma/FemtoPhiHM_allSB/parallelsmer/SidebandInspect.root";
  auto file = TFile::Open(name);
  //file->ls();
  TH1F* SB = (TH1F*)file->FindObjectAny("SBrightc");
  //file->TFile::Close();
  return SB;
}

TGraphAsymmErrors* GetBG(const char* path) {
  const char* name =
      Form("%s/SidebandInspect.root",
           path);
  auto file = TFile::Open(name);
  // file->ls();
  TGraphAsymmErrors* pK = (TGraphAsymmErrors*)file->FindObjectAny("pKBG");
  //file->TFile::Close();
  return pK;
}


//TGraphAsymmErrors* GetpKeff(const char* path) {
//  const char* name =
//      Form("%s/SidebandInspect.root",
//           path);
//  auto file = TFile::Open(name);
//  // file->ls();
//  TGraphAsymmErrors* pK = (TGraphAsymmErrors*)file->FindObjectAny("pKCFgaus");
//  return pK;
//}


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

double GetProjWidth(TH1* proj) {
  const Int_t nbins_transformed = proj->GetNbinsX();
  double bin_width = proj->GetXaxis()->GetBinWidth(1);
  int count = 0;
  for (int momTrans = 0; momTrans < nbins_transformed; momTrans++) {
    Double_t matrixvalues = 0.;
    matrixvalues = proj->GetBinContent(momTrans + 1);
    if (matrixvalues != 0) {
      count++;
    }
  }
  double width = count * bin_width;
  return width;
}

double GetExright(TH1* proj, double mean) {
  const Int_t nbins_transformed = proj->GetNbinsX();
  double bin_width = proj->GetXaxis()->GetBinWidth(1);
  int count = 0;
  int maxbin = 0;
  for (int momTrans = 0; momTrans < nbins_transformed; momTrans++) {
    Double_t matrixvalues = 0.;
    matrixvalues = proj->GetBinContent(momTrans + 1);
    if (matrixvalues != 0) {
      count++;
      maxbin = momTrans;
    }
  }
  double width = count * bin_width;
  double max = proj->GetBinCenter(maxbin) + 1;
  // cout<<"wifht: "<< width<<" maxbin: "<<maxbin<<"maxvaue: "<<max<<"mean: "<<
  // mean<<endl;

  double er = max - mean;

  return er;
}

double GetExleft(TH1* proj, double mean) {
  const Int_t nbins_transformed = proj->GetNbinsX();
  double bin_width = proj->GetXaxis()->GetBinWidth(1);
  int count = 0;
  int maxbin = 0;
  for (int momTrans = 0; momTrans < nbins_transformed; momTrans++) {
    Double_t matrixvalues = 0.;
    matrixvalues = proj->GetBinContent(momTrans + 1);
    if (matrixvalues != 0) {
      count++;
      maxbin = momTrans;
    }
  }
  double width = count * bin_width;

  double max = proj->GetBinCenter(maxbin) + 1;
  // cout<<"wifht: "<< width<<" maxbin: "<<maxbin<<"maxvaue: "<<max<<"mean: "<<
  // mean<<endl;

  double el = mean - (max - width);

  return el;
}




//float getchisquared(TH1F* hist, TGraph* func) {
//  int nbins = hist->GetNbinsX();
//  double chi2 = 0;
//  double tmp = 0;
//  cout << "nbins=" << nbins << endl;
//  double momentum = 0;
//  for (int i = 0; i < nbins; ++i) {
//    momentum = hist->GetBinCenter(i + 1);
//    if (momentum <= 1000) {
//      double funval = func->Eval(momentum);
//      double histval = hist->GetBinContent(i + 1);
//      double error = hist->GetBinError(i + 1);
//      tmp = (histval - funval) / error;
//      chi2 += tmp * tmp;

//      // cout<< "momentum: "<< momentum<< "   funval: "<<funval<<"   histval:
//      // "<< histval<< "   error: "<< error<<endl;
//    }
//  }
//  cout << "chi2: ..." << chi2 << endl;
//  return chi2;
//}


double myfunction1(double *x, double *par){
   double lim = par[0];

   double f = par[1]+par[2]*x[0]+par[3]*x[0]*x[0];
   double df = par[2]+2*par[3]*lim;
   double valflim = par[1]+par[2]*lim+par[3]*lim*lim;
   double d=valflim-df*lim;

  // cout<< lim<<" valf"<< valflim<< " d"<< d<<endl;

  if(x[0] < lim) return f;

  else {
       double f2 =  d+df*x[0];
       return f2;

  }

}


Double_t myGaus1(Double_t *x, Double_t *par){
   double lim = par[0];

   Double_t f = par[1]+par[2]*x[0]+par[3]*x[0]*x[0]+par[4] * TMath::Gaus(x[0], par[5], par[6], false);
   double df = par[2]+2*par[3]*lim;
   double valflim = par[1]+par[2]*lim+par[3]*lim*lim;
   double d=valflim-df*lim;

  // cout<< lim<<" valf"<< valflim<< " d"<< d<<endl;

  if(x[0] < lim) return f;

  else {
       Double_t f2 =  d+df*x[0]+par[4] * TMath::Gaus(x[0], par[5], par[6], false);
       return f2;

  }

}


Double_t myfunction2(Double_t *x, Double_t *par){
   Double_t f = par[0]+par[1]*x[0]+par[2]*x[0]*x[0];
   double df = par[1]+2*par[2]*x[0];
   Double_t x0=-par[1]/(2*par[2]);
   double valflx0 = par[0]+par[1]*x0+par[2]*x0*x0;
  if(x[0] < x0) return f;

  else {
       Double_t f2 =  valflx0;
       return f2;

  }

}


//const int nMJPars = 6;
const int nMJPars = 7;

auto MJFit = [ ] (double *x, double *p) {
  return p[0] + p[1] * x[0] + p[2] * x[0] * x[0] + p[3] * x[0] * x[0] * x[0] + p[4] * x[0] * x[0] * x[0] * x[0] + p[5] * x[0] * x[0] * x[0] * x[0] * x[0]+ p[6] * x[0] * x[0] * x[0] * x[0] * x[0]* x[0];
//  return p[0] + p[1] * x[0] + p[2] * x[0] * x[0] + p[3] * x[0] * x[0] * x[0] + p[4] * x[0] * x[0] * x[0] * x[0] ;

};


TH1F* GetMJ( const char* prefix,  int a) {
  const char* addon=Form("%d",a);
  const char* name=Form("/home/emma/FemtoPhiHM_ROTMC/phitruth10001200/CFOutput_pPhi_%s_%s.root",prefix, addon);
  auto file = TFile::Open(name);
  //file->ls();
  TH1F* corr = (TH1F*)file->FindObjectAny("hCk_ReweightedMeV_2");
  return corr;

}


Double_t myGaus2(Double_t *x, Double_t *par){
   Double_t f = par[0]+par[1]*x[0]+par[2]*x[0]*x[0]+par[3] * TMath::Gaus(x[0], par[4], par[5], false);
   double df = par[1]+2*par[2]*x[0];
   Double_t x0=-par[1]/(2*par[2]);
   double valflx0 = par[0]+par[1]*x0+par[2]*x0*x0;
  if(x[0] < x0) return f;

  else {
       Double_t f2 =  valflx0+par[3] * TMath::Gaus(x[0], par[4], par[5], false);
       return f2;

  }

}



float getchisquared(TH1F* hist, TH2F* smermatrix , Double_t *p) {
  int nbins = hist->GetNbinsX();
  double chi2 = 0;
  double tmp = 0;
  //cout << "nbins=" << nbins << endl;

  const int binwidth = 1;
  const int NumMomBins = 1000;
  //double momentum1 = 0;
  Double_t x[NumMomBins], y[NumMomBins], *momentum1[NumMomBins] ;
  for (Int_t i=0; i<NumMomBins; i++) {
    // momentum1 = hist->GetBinCenter(i + 1);
     x[i] = i*binwidth;
     momentum1[i]=&x[i];
     //y[i] = gausspol2(x[i] , p);


     y[i] = myGaus1(momentum1[i] , p);
     //y[i] = myGaus2(momentum1[i] , p);


//     if(x[i]<200){
//         cout<<x[i]<<"  yigaus;"<<y[i]<<endl;
//     }
   }

  TGraph *fitfungaus = new TGraph (NumMomBins, x, y);


  auto smearedfit=GetSmearedCF(fitfungaus,smermatrix);

  double IntegralCF =integrate(smearedfit,800,900);

  Double_t scale = 1/IntegralCF;

  auto normCF=GetSmearedCFscaled(fitfungaus, smermatrix, scale);

  double momentum = 0;
  for (int i = 0; i < nbins; ++i) {
    momentum = hist->GetBinCenter(i + 1);
    if (momentum <= 1000) {
      double funval = normCF->Eval(momentum);
      double histval = hist->GetBinContent(i + 1);
      //cout<<momentum<< "  fitval: "<<funval<<"  SB:"<<histval<<endl;
      double error = hist->GetBinError(i + 1);
      tmp = (histval - funval) / error;
      chi2 += tmp * tmp;

      // cout<< "momentum: "<< momentum<< "   funval: "<<funval<<"   histval:
      // "<< histval<< "   error: "<< error<<endl;
    }
  }
  //cout << "chi2: " << chi2 << endl;
  return chi2;
}


//float getchisquared(TF1* fit, TGraphAsymmErrors* pKeff) {
//  int ngraph = pKeff->GetN();
//  double chi2 = 0;
//  double tmp = 0;
//  cout << "ngraph=" << ngraph << endl;
//  double momentum = 0;
//  for (int i = 0; i < ngraph; ++i) {
//    Double_t x,y;
//    pKeff->GetPoint(i,x,y);
//    momentum =x;
//    if (momentum <= 1000) {
//      double fitval = func->Eval(momentum);
//      double funval = pKeff->Eval(momentum);
//      double errorxL = pKeff->GetErrorXlow(i);
//      double errorxH = pKeff->GetErrorXhigh(i);
//      double erroryL = pKeff->GetErrorYlow(i);
//      double erroryH = pKeff->GetErrorYhigh(i);

//      double meanX= (errorxL+errorxH)/2;
//      double meanY= (erroryL+erroryH)/2;

//      tmp = (funval - fitval) / meanY;
//      chi2 += tmp * tmp;

//      // cout<< "momentum: "<< momentum<< "   funval: "<<funval<<"   histval:
//      // "<< histval<< "   error: "<< error<<endl;
//    }
//  }
//  cout << "chi2: ..." << chi2 << endl;
//  return chi2;
//}


float getchisquared(TGraphAsymmErrors* pKeff, Double_t *p) {
  int ngraph = pKeff->GetN();
  double chi2 = 0;
  double tmp = 0;
  //cout << "ngraph=" << ngraph << endl;
  double *momentum = 0;
  for (int i = 0; i < ngraph; ++i) {
    Double_t x,y;
    pKeff->GetPoint(i,x,y);
   // cout<<"x: "<<x<<"   y: "<<y<<endl;
    momentum =&x;
    if (x <= 900) {
      double funval = pKeff->Eval(x);
      double erroryL = pKeff->GetErrorYlow(i);
      double erroryH = pKeff->GetErrorYhigh(i);
      double meanY= (erroryL+erroryH)/2;
      double fitval=myGaus1(momentum, p);
      //double fitval=myGaus2(momentum, p);

     // cout<<"gausval: "<<fitval<<"   error:"<<meanY<<endl;
      tmp = (funval - fitval) / meanY;
      chi2 += tmp * tmp;

      // cout<< "momentum: "<< momentum<< "   funval: "<<funval<<"   histval:
      // "<< histval<< "   error: "<< error<<endl;
    }
  }
  //cout << "chi2:" << chi2 << endl;
  return chi2;
}




//double getChisq(TF1* fit, TH1F* SBhist, TGraph* smerfunc) {
//  double chi2 = 0;
//  chi2 += fit->GetChisquare();
//  chi2 += getchisquared(SBhist, smerfunc);
//  return chi2;
//}


auto smearingMatrixright= GetMatix("/home/emma/FemtoPhiHM_allSB/parallelsmer", "histSmearTOTSBright");
auto smearingMatrixleft= GetMatix("/home/emma/FemtoPhiHM_allSB/parallelsmer", "histSmearTOTSBleft");
auto smearingMatrixPeak= GetMatix("/home/emma/FemtoPhiHM_allSB/parallelsmer", "histSmearTOTPeak");

auto pKeffective= GetpKeff("/home/emma/FemtoPhiHM_allSB/parallelsmer");
auto sidebandright= GetSB("/home/emma/FemtoPhiHM_allSB/parallelsmer", "SBrightc");
auto sidebandleft= GetSB("/home/emma/FemtoPhiHM_allSB/parallelsmer", "SBleftc");
auto pKBackground= GetBG("/home/emma/FemtoPhiHM_allSB/parallelsmer");


auto smearingMatrixrightBIG= GetMatix("/home/emma/FemtoPhiHM_allSB", "histSmearTOTSBright");
auto smearingMatrixleftBIG= GetMatix("/home/emma/FemtoPhiHM_allSB", "histSmearTOTSBleft");
auto smearingMatrixPeakBIG= GetMatix("/home/emma/FemtoPhiHM_allSB", "histSmearTOTPeak");

void getChisq(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */  ) {
  double chi2 = 0;
 // cout<<"tot chi2 calc:"<<endl;

  chi2 += getchisquared(pKeffective, p);;
  chi2 += getchisquared(sidebandright, smearingMatrixrightBIG, p);
  chi2 += getchisquared(sidebandleft, smearingMatrixleftBIG, p);

  fval = chi2;

  cout<<"tot chi2 calc:"<<chi2<<endl;


}


void getChipkeff(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */  ) {
  double chi2 = 0;
  //cout<<"tot chi2 calc:"<<endl;

  chi2 += getchisquared(pKeffective, p);;
  fval = chi2;

}




//void TFitter::GetConfidenceIntervals(Int_t n, Int_t ndim, const Double_t *x, Double_t *ci, Double_t cl)
// {
//    TF1 *f = (TF1*)fUserFunc;
//    Int_t npar = f->GetNumberFreeParameters();
//    Int_t npar_real = f->GetNpar();
//    Double_t *grad = new Double_t[npar_real];
//    Double_t *sum_vector = new Double_t[npar];
//    Bool_t *fixed=0;
//    Double_t al, bl;
//    if (npar_real != npar){
//       fixed = new Bool_t[npar_real];
//       memset(fixed,0,npar_real*sizeof(Bool_t));

//       for (Int_t ipar=0; ipar<npar_real; ipar++){
//          fixed[ipar]=0;
//          f->GetParLimits(ipar,al,bl);
//          if (al*bl != 0 && al >= bl) {
//             //this parameter is fixed
//             fixed[ipar]=1;
//          }
//       }
//    }
//    Double_t c=0;

//    for (Int_t ipar=0; ipar<npar_real; ipar++){

//  cout<<ipar<<" fixed: "<<fixed[ipar]<<endl;
//    }
//   // Double_t *matr = GetCovarianceMatrix();

//     Double_t *matr =TVirtualFitter::GetFitter()->GetCovarianceMatrix();
//    if (!matr){
//       delete [] grad;
//       delete [] sum_vector;
//       if (fixed)
//          delete [] fixed;
//       return;
//    }

//    Double_t t = TMath::StudentQuantile(0.5 + cl/2, f->GetNDF());
//    cout<< "t: "<<t<<endl;
//    Double_t chidf = TMath::Sqrt(f->GetChisquare()/f->GetNDF());
//    cout<< "chidf: "<<chidf<<endl;

//    Int_t igrad, ifree=0;
//    for (Int_t ipoint=0; ipoint<n; ipoint++){
//       c=0;
//       f->GradientPar(x+ndim*ipoint, grad);

//       cout<<"gradpar: "<<*x+ndim*ipoint<<" dim: "<<ndim<<" point: "<<ipoint<<" x: "<<*x<<" " <<x<<endl;
//       //multiply the covariance matrix by gradient
//       for (Int_t irow=0; irow<npar; irow++){
//          sum_vector[irow]=0;
//          igrad = 0;
//          for (Int_t icol=0; icol<npar; icol++){
//             igrad = 0;
//             ifree=0;
//             if (fixed) {
//                //find the free parameter #icol
//                while (ifree<icol+1){
//                   if (fixed[igrad]==0) ifree++;
//                   igrad++;
//                }
//                igrad--;
//                //now the [igrad] element of gradient corresponds to [icol] element of cov.matrix
//             } else {
//                igrad = icol;
//             }
//             sum_vector[irow]+=matr[irow*npar_real+icol]*grad[igrad];
//             cout<<"row "<<irow<< " col: "<<icol<<" vecsum: "<< sum_vector[irow]<<" grad"<<igrad<<": " <<grad[igrad]<<" matrix:"<<matr[irow*npar_real+icol]<<endl;

//          }
//          cout<<"row TOT "<<irow<< " vecsum: "<< sum_vector[irow]<<endl;
//       }
//       igrad = 0;
//       for (Int_t i=0; i<npar; i++){
//          igrad = 0; ifree=0;
//          if (fixed) {
//             //find the free parameter #icol
//             while (ifree<i+1){
//                if (fixed[igrad]==0) ifree++;
//                igrad++;
//             }
//             igrad--;
//          } else {
//             igrad = i;
//          }
//          c+=grad[igrad]*sum_vector[i];
//          cout<<i<<"c: "<<c<<" sumvec:"<<sum_vector[i]<<" igrad:"<<grad[igrad]<<endl;
//       }

//       c=TMath::Sqrt(c);
//       ci[ipoint]=c*t*chidf;
//       cout<<"point "<<ipoint<<c<<": "<<ci[ipoint]<<endl;
//    }

//    delete [] grad;
//    delete [] sum_vector;
//    if (fixed)
//       delete [] fixed;
// }





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

  auto filename = TString::Format("%s/SidebandParallel.root", InputDir.Data());
  auto outfile = new TFile(filename, "RECREATE");




  //auto smearingMatrix= GetMatix("/home/emma/FemtoPhiHM_allSB/parallelsmer", "histSmearTOTSBright");
  auto pKCF= GetpKeff("/home/emma/FemtoPhiHM_allSB/parallelsmer");
  auto sidebandRIGHT= GetSB("/home/emma/FemtoPhiHM_allSB/parallelsmer", "SBrightc");
  auto sidebandLEFT= GetSB("/home/emma/FemtoPhiHM_allSB/parallelsmer", "SBleftc");
  auto pphiCF= GetSB("/home/emma/FemtoPhiHM_allSB/parallelsmer", "Peakc");
  auto pKBG= GetBG("/home/emma/FemtoPhiHM_allSB/parallelsmer");



  TF1 *f13 = new TF1("fgpf13","gaus(0)+pol2(3)",0,650);

  f13->SetParameter(0,0.214);
  f13->SetParLimits(1,175,195); //188.117
  f13->SetParLimits(2,120,160); //141.988
  f13->SetParameter(3,0.075);
  f13->SetParLimits(3,0.6,0.80); //0.072,0.95
  f13->SetParameter(4,0.00069999);
  f13->SetParLimits(4,0.0007,0.00085);
  f13->SetParLimits(5,-0.0000006,-0.0000005);

double iniParams[6] = { 0.247383, 175.134, 160.0, 0.711071, 0.00085, -0.000000563462};

//double iniParams[10] = { 0.231929, 170.508, 160.0, 0.734362, 0.000759403, -0.000000483643};


TF1* fun = new TF1("fun2","pol2",0,800);
pKBG->Fit(fun);

int limit=600;
//--------------------------------FIT1-------------------------------------

  TF1* fun1 = new TF1("fun1",myfunction1,0,900,4);
  fun1->FixParameter(0,limit);
  fun1->SetParameter(1,fun->GetParameter(0));
  fun1->SetParameter(2,fun->GetParameter(1));
  fun1->SetParameter(3,fun->GetParameter(2));
  pKBG->Fit(fun1);


TF1* fungaus1 = new TF1("fun2",myGaus1,0,900,7);
fungaus1->FixParameter(0,fun1->GetParameter(0));
fungaus1->SetParameter(1,fun1->GetParameter(1));
fungaus1->SetParameter(2,fun1->GetParameter(2));
fungaus1->SetParameter(3,fun1->GetParameter(3));
fungaus1->SetParameter(4,0.214);
fungaus1->SetParLimits(5,175,195); //188.117
fungaus1->SetParLimits(6,120,160); //141.988

pKeffective->Fit(fungaus1);


//--------------------------------FIT2------------------------------------------

//TF1* fun1 = new TF1("fun1",myfunction2,0,800,3);
//fun1->SetParameter(0,fun->GetParameter(0));
//fun1->SetParameter(1,fun->GetParameter(1));
//fun1->SetParameter(2,fun->GetParameter(2));
//pKBG->Fit(fun1);


//TF1* fungaus1 = new TF1("fun2",myGaus2,0,800,6);
//fungaus1->SetParameter(0,fun1->GetParameter(0));
//fungaus1->SetParameter(1,fun1->GetParameter(1));
//fungaus1->SetParameter(2,fun1->GetParameter(2));
//fungaus1->SetParameter(3,0.214);
//fungaus1->SetParLimits(4,150,200); //188.117
//fungaus1->SetParLimits(5,120,160); //141.988

//pKeffective->Fit(fungaus1);


//----------------------------------------------------------------------------------------------------------------





int npo=25;

//pKCF->Fit(f14);

TGraphErrors *GraphConf1 = new TGraphErrors(npo);
GraphConf1->SetTitle("Fit Polynomial with .95 confidence;x;f(x)");
for (int i=0; i<npo; i++){
     GraphConf1->SetPoint(i, pKeffective->GetX()[i], 0);
     //cout<<"i:"<<i<<" "<< pKCF->GetX()[i] <<endl;

}
//(TVirtualFitter::GetFitter())->GetConfidenceIntervals(GraphConf1,0.68);
Double_t grad[7];
for (int i=0; i<npo; i++){
fungaus1->GradientPar(pKeffective->GetX()+1*npo, grad);
cout<<"x"<<i<<": "<<pKeffective->GetX()[i]<<endl;
for (int j=0; j<7; j++){
    cout<<"grad"<<i<<": "<<grad[j]<<endl;
}
}

(TVirtualFitter::GetFitter())->GetConfidenceIntervals(GraphConf1->GetN(),1,GraphConf1->GetX(), GraphConf1->GetEY(), 0.68);

for (Int_t i=0; i<GraphConf1->GetN(); i++)
{
   GraphConf1->SetPoint(i, GraphConf1->GetX()[i], ((TF1*)(fungaus1))->Eval(GraphConf1->GetX()[i]));

    cout<<i<<" x:"<<GraphConf1->GetX()[i]<< " y:"<<((TF1*)(fungaus1))->Eval(GraphConf1->GetX()[i])<< " ey"<<GraphConf1->GetEY()[i]<<endl;

}

//--------------------------------FIT1-------------------------------------

const int binwidth1 = 1;
const int NumMomBins1 = 1000;
float p10=fungaus1->GetParameter(0);
float p11=fungaus1->GetParameter(1);
float p12=fungaus1->GetParameter(2);
float p13=fungaus1->GetParameter(3);
float p14=fungaus1->GetParameter(4);
float p15=fungaus1->GetParameter(5);
float p16=fungaus1->GetParameter(6);
Int_t N1 =NumMomBins1;

Double_t x1[N1], y1[N1];

for (Int_t i = 0; i < N1; i++) {
  x1[i] = i * binwidth1;
  //y3[i] = p1+p2 * x3[i] +p3*x3[i] *x3[i] +p4* TMath::Gaus(x3[i], p5, p6, false);
  double lim = p10;

  Double_t f = p11+p12*x1[i]+p13*x1[i]*x1[i]+p14 * TMath::Gaus(x1[i], p15, p16, false);
  double df = p12+2*p13*lim;
  double valflim =p11+p12*lim+p13*lim*lim;
  double d=valflim-df*lim;

 // cout<< lim<<" valf"<< valflim<< " d"<< d<<endl;

 if(x1[i] < lim) y1[i]=f;

 else {
      y1[i] =  d+df*x1[i]+p14 *  TMath::Gaus(x1[i], p15, p16, false);
 }
  // cout<< x2[i]<< "    y= "<<y2[i]<<endl;
}



Int_t N2 =2000;

Double_t x2[N2], y2[N2];

for (Int_t i = 0; i < N2; i++) {
  x2[i] = i * binwidth1;
  //y3[i] = p1+p2 * x3[i] +p3*x3[i] *x3[i] +p4* TMath::Gaus(x3[i], p5, p6, false);
  double lim = p10;

  Double_t f = p11+p12*x2[i]+p13*x2[i]*x2[i]+p14 * TMath::Gaus(x2[i], p15, p16, false);
  double df = p12+2*p13*lim;
  double valflim =p11+p12*lim+p13*lim*lim;
  double d=valflim-df*lim;

 // cout<< lim<<" valf"<< valflim<< " d"<< d<<endl;

 if(x2[i] < lim) y2[i]=f;

 else {
      y2[i] =  d+df*x2[i]+p14 *  TMath::Gaus(x2[i], p15, p16, false);
 }
  // cout<< x2[i]<< "    y= "<<y2[i]<<endl;
}



//--------------------------------FIT2-------------------------------------

//const int binwidth1 = 1;
//const int NumMomBins1 = 1000;
//float p10=fungaus1->GetParameter(0);
//float p11=fungaus1->GetParameter(1);
//float p12=fungaus1->GetParameter(2);
//float p13=fungaus1->GetParameter(3);
//float p14=fungaus1->GetParameter(4);
//float p15=fungaus1->GetParameter(5);
//Int_t N1 =NumMomBins1;

//Double_t x1[N1], y1[N1];

//for (Int_t i = 0; i < N1; i++) {
//    x1[i] = i * binwidth1;
//    Double_t f = p10+p11*x1[i]+p12*x1[i]*x1[i]+p13 * TMath::Gaus(x1[i], p14, p15, false);
//    double df = p11+2*p12*x1[i];
//    Double_t x0=-p11/(2*p12);
//    double valflx0 = p10+p11*x0+p12*x0*x0;
//  //cout<< par[0]<< par[1]<< par[2]<<endl;
//    //cout<< x0 <<" valf"<< valflx0<< endl;

//   if(x1[i] < x0) {
//       y1[i]=f;
//   }
//   else {
//        y1[i]=valflx0+p13 * TMath::Gaus(x1[i], p14, p15, false);
//   }
//}

//Int_t N2 =2000;

//Double_t x2[N2], y2[N2];

//for (Int_t i = 0; i < N2; i++) {
//    x2[i] = i * binwidth1;
//    Double_t f = p10+p11*x2[i]+p12*x2[i]*x2[i]+p13 * TMath::Gaus(x2[i], p14, p15, false);
//    double df = p11+2*p12*x2[i];
//    Double_t x0=-p11/(2*p12);
//    double valflx0 = p10+p11*x0+p12*x0*x0;
//  //cout<< par[0]<< par[1]<< par[2]<<endl;
//    //cout<< x0 <<" valf"<< valflx0<< endl;

//   if(x2[i] < x0) {
//       y2[i]=f;
//   }
//   else {
//        y2[i]=valflx0+p13 * TMath::Gaus(x2[i], p14, p15, false);
//   }
//}




//----------------------------------------------------------------------------------------------------------------
//  TGraph* fitfunpol = new TGraph(nn, x2, y2);
//  TGraph* fitfungaus = new TGraph(nn, x3, y3);


TGraph *fitfungausold = new TGraph (NumMomBins1, x1, y1);


auto smearedfitrightold=GetSmearedCF(fitfungausold,smearingMatrixright);

double IntegralCFrold =integrate(smearedfitrightold,800,900);

Double_t scalerold = 1/IntegralCFrold;

auto normCFrightold=GetSmearedCFscaled(fitfungausold, smearingMatrixright, scalerold);



auto smearedfitleftold=GetSmearedCF(fitfungausold,smearingMatrixleft);

double IntegralCFlold =integrate(smearedfitleftold,800,900);

Double_t scalelold = 1/IntegralCFlold;

auto normCFleftold=GetSmearedCFscaled(fitfungausold, smearingMatrixleft, scalelold);



auto smearedfitpeakold=GetSmearedCF(fitfungausold,smearingMatrixPeak);

double IntegralCFold =integrate(smearedfitpeakold,800,900);

Double_t scaleold = 1/IntegralCFold;

auto normCFold=GetSmearedCFscaled(fitfungausold, smearingMatrixPeak, scaleold);


cout<<"P11"<<endl;




TGraph *fitfungausold2 = new TGraph (N2, x2, y2);


auto smearedfitrightold2=GetSmearedCF(fitfungausold2,smearingMatrixrightBIG);

double IntegralCFrold2 =integrate(smearedfitrightold2,800,900);

Double_t scalerold2 = 1/IntegralCFrold2;

auto normCFrightold2=GetSmearedCFscaled(fitfungausold2, smearingMatrixrightBIG, scalerold2);



auto smearedfitleftold2=GetSmearedCF(fitfungausold2,smearingMatrixleftBIG);

double IntegralCFlold2 =integrate(smearedfitleftold2,800,900);

Double_t scalelold2 = 1/IntegralCFlold2;

auto normCFleftold2=GetSmearedCFscaled(fitfungausold2, smearingMatrixleftBIG, scalelold2);



auto smearedfitpeakold2=GetSmearedCF(fitfungausold2,smearingMatrixPeakBIG);

double IntegralCFold2 =integrate(smearedfitpeakold2,800,900);

Double_t scaleold2 = 1/IntegralCFold2;

auto normCFold2=GetSmearedCFscaled(fitfungausold2, smearingMatrixPeakBIG, scaleold2);


//--------------------------------FIT1-------------------------------------

  TVirtualFitter::SetDefaultFitter("Minuit");
  TVirtualFitter * minuit = TVirtualFitter::Fitter(0,7);

  minuit->SetParameter(0, fungaus1->GetParName(0), fungaus1->GetParameter(0), 0, 0,0);
  minuit->SetParameter(1, fungaus1->GetParName(1), fungaus1->GetParameter(1), fungaus1->GetParameter(1)/100, 0,0);
  minuit->SetParameter(2, fungaus1->GetParName(2), fungaus1->GetParameter(2), fungaus1->GetParameter(2)/100, 0,0);
  minuit->SetParameter(3, fungaus1->GetParName(3), fungaus1->GetParameter(3), fungaus1->GetParameter(3)/100, 0,0);
  minuit->SetParameter(4, fungaus1->GetParName(4), fungaus1->GetParameter(4), fungaus1->GetParameter(4)/100, 0,0);
  minuit->SetParameter(5, fungaus1->GetParName(5), fungaus1->GetParameter(5), fungaus1->GetParameter(5)/100, 0,0);
  minuit->SetParameter(6, fungaus1->GetParName(6), fungaus1->GetParameter(6), fungaus1->GetParameter(6)/100, 0,0);

  minuit->SetFCN(getChisq);

  double arglist[100];
  arglist[0] = 0;
  minuit->ExecuteCommand("SET PRINT",arglist,2);


  arglist[0] = 1.0;
  minuit->ExecuteCommand("SET ERRordef",arglist,1);

  cout<<"-----------MIGRAD"<<endl;

  arglist[0] = 10000; // number of function callsn 1000
  arglist[1] = 0.05; // tolerance 0.1

  minuit->ExecuteCommand("MIGRAD",arglist,2);

  cout<<"-----------HESSE"<<endl;

  arglist[0] = 10000; // number of function callsn 1000

  minuit->ExecuteCommand("HESSE",arglist,1);

//  get result


  double minParams[7];
  double parErrors[7];
  for (int i = 0; i < 7; ++i) {
    minParams[i] = minuit->GetParameter(i);
    parErrors[i] = minuit->GetParError(i);
    cout<<i<<". minpar: "<<minParams[i]<<"   mparerror: "<<parErrors[i] <<endl;
  }
  double chi2, edm, errdef;
  int nvpar, nparx;
  minuit->GetStats(chi2,edm,errdef,nvpar,nparx);
  fungaus1->SetParameters(minParams);
  fungaus1->SetParErrors(parErrors);
  fungaus1->SetChisquare(chi2);
//  int ndf = coords.size()-nvpar;
  int ndf =25+25+25;
  fungaus1->SetNDF(ndf);
  std::cout << "Chi2 Fit = " << chi2 << " ndf = " << ndf << "  " << fungaus1->GetNDF() << std::endl;
  // add to list of functions
  fungaus1->FixParameter(3,fungaus1->GetParameter(3));

  pKeffective->GetListOfFunctions()->Add(fungaus1);
  //sideband->GetListOfFunctions()->Add(f13);

  int paramnum=(TVirtualFitter::GetFitter())->GetNumberFreeParameters ();
  int paratot=fungaus1->GetNpar();
  int praf1=fungaus1->GetNumberFreeParameters();



cout<<"parafree fitter: "<<paramnum<< " vs funtion"<<praf1<<" par3 value:"<<   fun1->GetParameter(3)<<" paratot:"<<paratot<<endl;





cout<<"lalalal"<<endl;
  const int binwidth = 1;
  const int NumMomBins = 1000;
  float p0=fungaus1->GetParameter(0);
  float p1=fungaus1->GetParameter(1);
  float p2=fungaus1->GetParameter(2);
  float p3=fungaus1->GetParameter(3);
  float p4=fungaus1->GetParameter(4);
  float p5=fungaus1->GetParameter(5);
  float p6=fungaus1->GetParameter(6);
  Int_t nn =NumMomBins;


  cout<<"P88"<<endl;




  Double_t x3[nn], y3[nn];
  Double_t *momentum[nn];
  for (Int_t i = 0; i < nn; i++) {
    x3[i] = i * binwidth;
    //y3[i] = p1+p2 * x3[i] +p3*x3[i] *x3[i] +p4* TMath::Gaus(x3[i], p5, p6, false);
    double lim = p0;

    Double_t f = p1+p2*x3[i]+p3*x3[i]*x3[i]+p4 * TMath::Gaus(x3[i], p5, p6, false);
    double df = p2+2*p3*lim;
    double valflim =p1+p2*lim+p3*lim*lim;
    double d=valflim-df*lim;

   // cout<< lim<<" valf"<< valflim<< " d"<< d<<endl;

   if(x3[i] < lim) y3[i]=f;

   else {
        y3[i] =  d+df*x3[i]+p4 *  TMath::Gaus(x3[i], p5, p6, false);
   }
    // cout<< x2[i]<< "    y= "<<y2[i]<<endl;
  }



  TGraph *fitfungaus = new TGraph (NumMomBins, x3, y3);


  Int_t nn2=2000;
  Double_t x4[nn2];
  Double_t y4[nn2];
  for (Int_t i = 0; i < nn2; i++) {
    x4[i] = i * binwidth;
    //y3[i] = p1+p2 * x4[i] +p3*x4[i] *x4[i] +p4* TMath::Gaus(x4[i], p5, p6, false);
    double lim = p0;

    Double_t f = p1+p2*x4[i]+p3*x4[i]*x4[i]+p4 * TMath::Gaus(x4[i], p5, p6, false);
    double df = p2+2*p3*lim;
    double valflim =p1+p2*lim+p3*lim*lim;
    double d=valflim-df*lim;

   // cout<< lim<<" valf"<< valflim<< " d"<< d<<endl;

   if(x4[i] < lim) y4[i]=f;

   else {
        y4[i] =  d+df*x4[i]+p4 *  TMath::Gaus(x4[i], p5, p6, false);
   }
    // cout<< x2[i]<< "    y= "<<y2[i]<<endl;
  }



  TGraph *fitfungaus2 = new TGraph (nn2, x4, y4);



  //--------------------------------FIT2-------------------------------------




//TVirtualFitter::SetDefaultFitter("Minuit");
//TVirtualFitter * minuit = TVirtualFitter::Fitter(0,6);

////minuit->SetParameter(0, fungaus1->GetParName(0), fungaus1->GetParameter(0), fungaus1->GetParameter(0)/100, 0.0,1.0);
//minuit->SetParameter(0, fungaus1->GetParName(0), fungaus1->GetParameter(0), fungaus1->GetParameter(0)/100, 0,0);
//minuit->SetParameter(1, fungaus1->GetParName(1), fungaus1->GetParameter(1), fungaus1->GetParameter(1)/100, 0,0);
//minuit->SetParameter(2, fungaus1->GetParName(2), fungaus1->GetParameter(2), fungaus1->GetParameter(2)/100, 0,0);
//minuit->SetParameter(3, fungaus1->GetParName(3), fungaus1->GetParameter(3), fungaus1->GetParameter(3)/100, 0,0);
//minuit->SetParameter(4, fungaus1->GetParName(4), fungaus1->GetParameter(4), fungaus1->GetParameter(4)/100, 0,0);
//minuit->SetParameter(5, fungaus1->GetParName(5), fungaus1->GetParameter(5), fungaus1->GetParameter(5)/100, 0,0);

//minuit->SetFCN(getChisq);

//double arglist[100];
//arglist[0] = 0;
//minuit->ExecuteCommand("SET PRINT",arglist,2);


//arglist[0] = 1.0;
//minuit->ExecuteCommand("SET ERRordef",arglist,1);

//cout<<"-----------MIGRAD"<<endl;

//arglist[0] = 10000; // number of function callsn 1000
//arglist[1] = 0.05; // tolerance 0.1

//minuit->ExecuteCommand("MIGRAD",arglist,2);

//cout<<"-----------HESSE"<<endl;

//arglist[0] = 10000; // number of function callsn 1000

//minuit->ExecuteCommand("HESSE",arglist,1);

////get result


//double minParams[6];
//double parErrors[6];
//for (int i = 0; i < 6; ++i) {
//  minParams[i] = minuit->GetParameter(i);
//  parErrors[i] = minuit->GetParError(i);
//  cout<<i<<". minpar: "<<minParams[i]<<"   mparerror: "<<parErrors[i] <<endl;
//}
//double chi2, edm, errdef;
//int nvpar, nparx;
//minuit->GetStats(chi2,edm,errdef,nvpar,nparx);
//fungaus1->SetParameters(minParams);
//fungaus1->SetParErrors(parErrors);
//fungaus1->SetChisquare(chi2);
////  int ndf = coords.size()-nvpar;
//int ndf =25+25+25;
//fungaus1->SetNDF(ndf);
//std::cout << "Chi2 Fit = " << chi2 << " ndf = " << ndf << "  " << fungaus1->GetNDF() << std::endl;
//// add to list of functions
//fungaus1->FixParameter(2,fungaus1->GetParameter(2));

//pKeffective->GetListOfFunctions()->Add(fungaus1);
////sideband->GetListOfFunctions()->Add(f13);

//int paramnum=(TVirtualFitter::GetFitter())->GetNumberFreeParameters ();
//int paratot=fungaus1->GetNpar();
//int praf1=fungaus1->GetNumberFreeParameters();



//cout<<"parafree fitter: "<<paramnum<< " vs funtion"<<praf1<<" par3 value:"<<   fun1->GetParameter(3)<<" paratot:"<<paratot<<endl;





//cout<<"lalalal"<<endl;
//const int binwidth = 1;
//const int NumMomBins = 1000;
//float p0=fungaus1->GetParameter(0);
//float p1=fungaus1->GetParameter(1);
//float p2=fungaus1->GetParameter(2);
//float p3=fungaus1->GetParameter(3);
//float p4=fungaus1->GetParameter(4);
//float p5=fungaus1->GetParameter(5);
//Int_t nn =NumMomBins;


//cout<<p0<<" "<< p1<< " "<< p2<<" "<< p3<<" "<< p4<<" "<< p5<<endl;




//Double_t x3[nn], y3[nn];
//Double_t *momentum[nn];
//for (Int_t i = 0; i < nn; i++) {
//    x3[i] = i * binwidth;
//    Double_t f = p0+p1*x3[i]+p2*x3[i]*x3[i]+p3 * TMath::Gaus(x3[i], p4, p5, false);
//    double df = p1+2*p2*x3[i];
//    Double_t x0=-p1/(2*p2);
//    double valflx0 = p0+p1*x0+p2*x0*x0;
//  //cout<< par[0]<< par[1]<< par[2]<<endl;
//   //cout<< x0 <<" valf"<< valflx0<< endl;

//   if(x3[i] < x0) {
//       y3[i]=f;
//   }
//   else {
//        y3[i]=valflx0+p3 * TMath::Gaus(x3[i], p4, p5, false);
//   }
//}



//TGraph *fitfungaus = new TGraph (NumMomBins, x3, y3);

//Int_t nn2=2000;
//Double_t x4[nn2], y4[nn2];
//for (Int_t i = 0; i < nn2; i++) {
//    x4[i] = i * binwidth;
//    Double_t f = p0+p1*x4[i]+p2*x4[i]*x4[i]+p3 * TMath::Gaus(x4[i], p4, p5, false);
//    double df = p1+2*p2*x4[i];
//    Double_t x0=-p1/(2*p2);
//    double valflx0 = p0+p1*x0+p2*x0*x0;
//  //cout<< par[0]<< par[1]<< par[2]<<endl;
//   //cout<< x0 <<" valf"<< valflx0<< endl;

//   if(x4[i] < x0) {
//       y4[i]=f;
//   }
//   else {
//        y4[i]=valflx0+p3 * TMath::Gaus(x4[i], p4, p5, false);
//   }
//}



//TGraph *fitfungaus2 = new TGraph (nn2, x4, y4);

//----------------------------------------------------------------------------------------------------------------

  auto smearedfitright=GetSmearedCF(fitfungaus,smearingMatrixright);

  double IntegralCFr =integrate(smearedfitright,800,900);

  Double_t scaler = 1/IntegralCFr;

  auto normCFright=GetSmearedCFscaled(fitfungaus, smearingMatrixright, scaler);

  cout<<"P10"<<endl;


  auto smearedfitleft=GetSmearedCF(fitfungaus,smearingMatrixleft);

  double IntegralCFl =integrate(smearedfitleft,800,900);

  Double_t scalel = 1/IntegralCFl;

  auto normCFleft=GetSmearedCFscaled(fitfungaus, smearingMatrixleft, scalel);

  cout<<"P11"<<endl;



  auto smearedfitpeak=GetSmearedCF(fitfungaus,smearingMatrixPeak);

  double IntegralCF =integrate(smearedfitpeak,800,900);

  Double_t scale = 1/IntegralCF;

  auto normCF=GetSmearedCFscaled(fitfungaus, smearingMatrixPeak, scale);






  cout<<"lalalal"<<endl;



  auto smearedfitrightnew2=GetSmearedCF(fitfungaus2,smearingMatrixrightBIG);

  double IntegralCFrnew2 =integrate(smearedfitrightnew2,800,900);

  Double_t scalernew2 = 1/IntegralCFrnew2;

  auto normCFrightnew2=GetSmearedCFscaled(fitfungaus2, smearingMatrixrightBIG, scalernew2);



  auto smearedfitleftnew2=GetSmearedCF(fitfungaus2,smearingMatrixleftBIG);

  double IntegralCFlnew2 =integrate(smearedfitleftnew2,800,900);

  Double_t scalelnew2 = 1/IntegralCFlnew2;

  auto normCFleftnew2=GetSmearedCFscaled(fitfungaus2, smearingMatrixleftBIG, scalelnew2);



  auto smearedfitpeaknew2=GetSmearedCF(fitfungaus2,smearingMatrixPeakBIG);

  double IntegralCFnew2 =integrate(smearedfitpeaknew2,800,900);

  Double_t scalenew2 = 1/IntegralCFnew2;

  auto normCFnew2=GetSmearedCFscaled(fitfungaus2, smearingMatrixPeakBIG, scalenew2);


//---------------------------------------------------  auto scaledCF=normCF*0.3;

cout<< "minijet FITTTTTTTTT-------------------------"<<endl;
  int nP=1000;

  TGraph *Graph=new TGraph(nP);
   for (int i=0; i<nP; i++){
         Graph->SetPoint(i, normCF->GetX()[i], 1+((normCF->GetY()[i]-1)*0.45));

}


   TH1* MJ = GetMJ("HMPhi",0);

 //  TF1 *MJFitFct = new TF1("sidebandFit", "pol4", 0, 1200);
   TF1 *MJFitFct = new TF1("sidebandFit", MJFit, 0, 2000,
                                 nMJPars);  // TODO check that the ranges make sense
   MJ->Fit(MJFitFct, "R");



   const int binwidthmj = 1;
   const int NumMomBinsmj = 1000;
   float mjp1=MJFitFct->GetParameter(0);
   float mjp2=MJFitFct->GetParameter(1);
   float mjp3=MJFitFct->GetParameter(2);
   float mjp4=MJFitFct->GetParameter(3);
   float mjp5=MJFitFct->GetParameter(4);
   float mjp6=MJFitFct->GetParameter(5);
   float mjp7=MJFitFct->GetParameter(6);

   Int_t Nmj =NumMomBinsmj;


cout<<mjp1<<" "<<mjp2<<" "<<mjp3<<" "<<mjp4<<" "<<mjp5<<" "<<endl;
   Double_t x11mj[Nmj], y11mj[Nmj];

   for (int i = 0; i < Nmj; ++i) {
     Double_t x,y;
     normCFnew2->GetPoint(i,x,y);
     x11mj[i]=x;
     y11mj[i]=y;
     cout<< x11mj[i]<<" : "<<  y11mj[i]<<endl;
   }

   TGraph *SBGraph = new TGraph (NumMomBinsmj, x11mj, y11mj);


   Double_t x1mj[Nmj], y1mj[Nmj];

   for (Int_t i = 0; i < Nmj; i++) {
    x1mj[i] = x11mj[i];
    y1mj[i] = mjp1+mjp2*x1mj[i]+mjp3*x1mj[i]*x1mj[i] +mjp4*x1mj[i]*x1mj[i]*x1mj[i] +mjp5*x1mj[i]*x1mj[i]*x1mj[i]*x1mj[i]+mjp6*x1mj[i]*x1mj[i]*x1mj[i]*x1mj[i]*x1mj[i]+mjp7*x1mj[i]*x1mj[i]*x1mj[i]*x1mj[i]*x1mj[i]*x1mj[i];
    //1mj[i] = mjp1+mjp2*x1mj[i]+mjp3*x1mj[i]*x1mj[i] +mjp4*x1mj[i]*x1mj[i]*x1mj[i] +mjp5*x1mj[i]*x1mj[i]*x1mj[i]*x1mj[i];

     cout<< x1mj[i]<<" : "<<  y1mj[i]<<endl;
   }


   TGraph *MJGraph = new TGraph (NumMomBinsmj, x1mj, y1mj);


   Double_t x2mj[Nmj], y2mj[Nmj];


   for (Int_t i = 0; i < Nmj; i++) {
     x2mj[i] = x11mj[i];
     y2mj[i] = (y11mj[i])/(y1mj[i]) ;
     cout<< x2mj[i]<<" : "<<  y2mj[i]<<endl;

   }

   TGraph *Diff = new TGraph (NumMomBinsmj, x2mj, y2mj);







   Double_t grad2[7];
   for (int i=0; i<npo; i++){
   fungaus1->GradientPar(pKeffective->GetX()+1*npo, grad);
   cout<<"x"<<i<<": "<<pKeffective->GetX()[i]<<endl;
   for (int j=0; j<7; j++){
       cout<<"grad"<<i<<": "<<grad2[j]<<endl;
   }
   }



   TVirtualFitter::GetFitter()->SetUserFunc(fungaus1);
   (TVirtualFitter::GetFitter())->SetObjectFit(pKeffective);


   int nPci=25;

   //pKCF->Fit(f14);

  TGraphErrors *GraphConf = new TGraphErrors(nPci);
  GraphConf->SetTitle("Fit Polynomial with .95 confidence;x;f(x)");
  for (int i=0; i<nPci; i++){
        GraphConf->SetPoint(i, pKeffective->GetX()[i], 0);
        //cout<<"i:"<<i<<" "<< pKCF->GetX()[i] <<endl;

  }
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(GraphConf,0.68);



  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(GraphConf->GetN(),1,GraphConf->GetX(), GraphConf->GetEY(), 0.68);

   for (Int_t i=0; i<GraphConf->GetN(); i++)
   {
      GraphConf->SetPoint(i, GraphConf->GetX()[i], ((TF1*)(fungaus1))->Eval(GraphConf->GetX()[i]));

       cout<<i<<" x:"<<GraphConf->GetX()[i]<< " y:"<<((TF1*)(fungaus1))->Eval(GraphConf->GetX()[i])<< " ey"<<GraphConf->GetEY()[i]<<endl;

   }


  for (int i=0;i<6;i++){
      for (int j=0;j<6;j++){

          cout<<i<<j<<" fittr: "<<TVirtualFitter::GetFitter()->GetCovarianceMatrixElement(i,j)<<"minuit: "<<minuit->GetCovarianceMatrixElement(i,j)<<endl;

  }
  }




  const float right = 0.04;
  const float top = 0.025;

  auto dummyHist = new TH1F("dummyHist",
                            ";#it{k}* (MeV/#it{c}); #it{C}(#it{k}*)", 100, 4,
                            1000);
  DreamPlot::SetStyleHisto(dummyHist, 20, kGreen + 2);

  auto gc3336 = new TCanvas("pKBG", "pKBG", 0, 0, 650, 550);
  gc3336->SetRightMargin(right);
  gc3336->SetTopMargin(top);
  dummyHist->Draw();
  //dummyHist->GetYaxis()->SetRangeUser(0.8, 1.1);
  //dummyHist->GetXaxis()->SetRangeUser(0.0, 200);

  //dummyHist->GetXaxis()->SetNdivisions(504);
  GraphConf->SetLineColor(kRed);
  GraphConf->Draw("ap");
  pKeffective->SetMarkerStyle(5);
  pKeffective->SetMarkerSize(0.7);
  pKeffective->Draw("ALPsame");
  gc3336->Print("confintss.pdf");
  gc3336->Print("confintss.png");

//  ConfintUP->SetLineColor(kRed);
//  ConfintLOW->SetLineColor(kRed);



  outfile->cd();

//gr01->Write("confint");
Graph->Write("dd");
  cout<<"ppqp"<<endl;

  pKeffective->Write("pKCFgaus");
  pKBG->Write("pKBGpol1");


  fungaus1->Write("fit");

  f13->Write("orifit");
//  GraphConf2->Write("ddd");
 // f15->Write("refit");

 GraphConf->Write("Confint");
 GraphConf1->Write("Confint1");

fitfungausold->Write("GRAPH gauss old");
fitfungaus->Write("GRAPH gauss new");


  sidebandRIGHT->Write("rightsb");

  normCFright->Write("rightsbSMEAR");

  sidebandLEFT->Write("leftsb");

  normCFleft->Write("leftsbSMEAR");

  pphiCF->Write("CF");

  normCF->Write("cfSMEAR");




  normCFrightold->Write("rightsbSMEARold");

  normCFleftold->Write("leftsbSMEARold");

  normCFold->Write("cfSMEARold");





  normCFrightnew2->Write("rightsbSMEAR2");

  normCFleftnew2->Write("leftsbSMEAR2");

  normCFnew2->Write("cfSMEAR2");




  normCFrightold2->Write("rightsbSMEARold2");

  normCFleftold2->Write("leftsbSMEARold2");

  normCFold2->Write("cfSMEARold2");


  MJ->Write("MinijetBG/MCTruth");
  MJGraph->Write("MJFit");
  Diff->Write("Ratio");

  SBGraph->Write("SBgrapg");
  outfile->Close();
}

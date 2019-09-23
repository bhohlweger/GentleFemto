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
std::vector<int> fFillColors = {kGray+1, kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9, kCyan-8, kYellow-7};
std::vector<int> fColors     = {kBlack, kRed+1 , kBlue+2, kGreen+3, kMagenta+1, kOrange-1, kCyan+2, kYellow+2};
std::vector<int> fMarkers    = {kFullCircle, kFullSquare, kOpenCircle, kOpenSquare, kOpenDiamond, kOpenCross, kFullCross, kFullDiamond, kFullStar, kOpenStar};

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// void SetStyle(Bool_t graypalette) {
//   gStyle->Reset("Plain");
//   gStyle->SetOptTitle(0);
//   gStyle->SetOptStat(0);
//   if(graypalette) gStyle->SetPalette(8,0);
//   else gStyle->SetPalette(1);
//   gStyle->SetCanvasColor(10);
//   gStyle->SetCanvasBorderMode(0);
//   gStyle->SetFrameLineWidth(1);
//   gStyle->SetFrameFillColor(kWhite);
//   gStyle->SetPadColor(10);
//   gStyle->SetPadTickX(1);
//   gStyle->SetPadTickY(1);
//   gStyle->SetPadBottomMargin(0.15);
//   gStyle->SetPadLeftMargin(0.15);
//   gStyle->SetHistLineWidth(1);
//   gStyle->SetHistLineColor(kRed);
//   gStyle->SetFuncWidth(2);
//   gStyle->SetFuncColor(kGreen);
//   gStyle->SetLineWidth(2);
//   gStyle->SetLabelSize(0.045,"xyz");
//   gStyle->SetLabelOffset(0.01,"y");
//   gStyle->SetLabelOffset(0.01,"x");
//   gStyle->SetLabelColor(kBlack,"xyz");
//   gStyle->SetTitleSize(0.05,"xyz");
//   gStyle->SetTitleOffset(1.25,"y");
//   gStyle->SetTitleOffset(1.2,"x");
//   gStyle->SetTitleFillColor(kWhite);
//   gStyle->SetTextSizePixels(26);
//   gStyle->SetTextFont(42);
//   //  gStyle->SetTickLength(0.04,"X");  gStyle->SetTickLength(0.04,"Y");
//
//   gStyle->SetLegendBorderSize(0);
//   gStyle->SetLegendFillColor(kWhite);
//   //  gStyle->SetFillColor(kWhite);
//   gStyle->SetLegendFont(42);
//
//
// }
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// a and c are the weighting entities, i.e. wMean = (weightA*A + weightB*d) / (a+weightB)
float weightedMean(float weightA, float A, float weightB, float B)
{
  return (weightA*A + weightB*B)/(weightA+weightB);
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// a and weightB are the weighting entities, i.e. wMean = (weightA*A + weightB*B) / (a+weightB)
float weightedMeanError(float weightA, float A, float weightB, float B, float weightAErr, float AErr, float weightBErr, float BErr)
{
  return std::sqrt(weightAErr*weightAErr*std::pow(((weightB*(A - B))/((weightA + weightB)*(weightA + weightB))), 2) + AErr*AErr* std::pow(weightA/(weightA + weightB), 2)  + weightBErr*weightBErr*std::pow((weightA* (-A + B))/((weightA + weightB)*(weightA + weightB)) , 2) + BErr*BErr* std::pow(weightB/(weightA + weightB), 2));
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TH1F *getSignalHisto(TF1 *function, TH1F *histo, float rangeLow, float rangeHigh, const char *name)
{
  const int firstBin = histo->FindBin(rangeLow);
  const int lastBin  = histo->FindBin(rangeHigh);
  TH1F *result = new TH1F(Form("result_%f_%f_%s", rangeLow, rangeHigh, name), "", histo->GetNbinsX(), histo->GetXaxis()->GetXmin(), histo->GetXaxis()->GetXmax());
  for(int i = firstBin; i<lastBin; ++i) {
    float weight = histo->GetBinContent(i) - function->Eval(histo->GetBinCenter(i));
    result->Fill(histo->GetBinCenter(i), weight);
    result->SetBinError(i, histo->GetBinError(i));
  }
  result->SetFillColor(fFillColors[0]);
  result->SetLineColor(fFillColors[0]);
  result->Sumw2();
  return result;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void FitLambda(TH1F* histo, float &signal, float &signalErr, float &background, float &backgroundErr, float lowerBound, float upperBound)
{
  histo->Sumw2();
  // Fit Background with second order polynomial, excluding Mlambda +/- 10 MeV, f("f", obj, xmin, xmax , npar)
  TF1 *fBackground = new TF1("fBackground", [&](double *x, double *p) { if (x[0] > 1.1075 && x[0] < 1.1235) {TF1::RejectPoint(); return (double)0; } return p[0] + p[1]*x[0] + p[2]*x[0]*x[0]; }, 1.095, 1.15, 3);
  TFitResultPtr backgroundR = histo->Fit("fBackground", "SRQ0", "", 1.095, 1.15);

  // parse then to proper TF1
  TF1 *fBackground2 = new TF1("fBackground2","pol2", 0, 1.5);
  fBackground2->SetParameter(0, fBackground->GetParameter(0));
  fBackground2->SetParameter(1, fBackground->GetParameter(1));
  fBackground2->SetParameter(2, fBackground->GetParameter(2));

  // remove background from signal
  TH1F *signalOnly = getSignalHisto(fBackground2, histo, 1.0, 1.3, Form("%s_signal_only", histo->GetName()));
  signalOnly->Sumw2();
  signalOnly->Draw("same");

  // fit signal only
  TF1 *fSignalSingleGauss = new TF1("fSignalSingleGauss", "gaus", 1.095, 1.15);
  //  fSignalSingleGauss->SetParameter(1, 1.115);
  signalOnly->Fit("fSignalSingleGauss", "SRQ0", "", 1.1075, 1.1235);

  TF1 *fSignalGauss = new TF1("fSignalGauss", "gaus(0) + gaus(3)", 1.1, 1.3);
  fSignalGauss->SetParameter(0, 0.05 * histo->GetMaximum());
  fSignalGauss->SetParameter(1, fSignalSingleGauss->GetParameter(1));
  fSignalGauss->SetParLimits(1, 1.115-0.01, 1.115+0.01);
  fSignalGauss->SetParameter(2, 5.f*fSignalSingleGauss->GetParameter(2));
  fSignalGauss->SetParameter(3, 0.95 * histo->GetMaximum());
  fSignalGauss->SetParameter(4, fSignalSingleGauss->GetParameter(1));
  fSignalGauss->SetParLimits(4, 1.115-0.01, 1.115+0.01);
  fSignalGauss->SetParameter(5, fSignalSingleGauss->GetParameter(2));
  TFitResultPtr r = signalOnly->Fit("fSignalGauss", "SRQ0", "", 1.1075, 1.1235);

  // Extract signal as integral
  signal = fSignalGauss->Integral(lowerBound, upperBound) /double(histo->GetBinWidth(1));
  signalErr = fSignalGauss->IntegralError(lowerBound, upperBound, r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray()) /double(histo->GetBinWidth(1));

  TF1 *fLambda = new TF1("fLambda", "fBackground2 + fSignalGauss", 1.1, 1.13);
  fLambda->SetNpx(1000);
  fLambda->SetParameter(3, 0.75 * histo->GetMaximum());
  fLambda->SetParameter(4, fSignalGauss->GetParameter((1)));
  fLambda->SetParameter(5, fSignalGauss->GetParameter((2)));
  fLambda->SetParameter(6, 0.2 * histo->GetMaximum());
  fLambda->SetParameter(7, fSignalGauss->GetParameter((4)));
  fLambda->SetParameter(8, fSignalGauss->GetParameter((5)));
  fLambda->SetLineColor(fColors[1]);
  histo->Fit("fLambda", "SRQ", "", 1.095, 1.15);

  TF1 *fLambda_background = new TF1("fLambda_background", "pol2(0)", 1.05, 1.25);
  fLambda_background->SetParameter(0, fLambda->GetParameter(0));
  fLambda_background->SetParameter(1, fLambda->GetParameter(1));
  fLambda_background->SetParameter(2, fLambda->GetParameter(2));
  fLambda_background->SetLineStyle(3);
  fLambda_background->SetLineColor(fColors[1]);

  background = fLambda_background->Integral(lowerBound, upperBound) /double(histo->GetBinWidth(1));
  backgroundErr = fLambda_background->IntegralError(lowerBound, upperBound, backgroundR->GetParams(), backgroundR->GetCovarianceMatrix().GetMatrixArray()) /double(histo->GetBinWidth(1));

  histo->GetListOfFunctions()->Add(fLambda_background);

  delete signalOnly;
  delete fSignalGauss;
  delete fSignalSingleGauss;
}
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void BbarB_QAplots(const char* filename, const char* prefix,
                     const char* addon = "", const char* date = "") {

std::cout<<"----------------------"<<std::endl;

std::cout<<date<<std::endl;

DreamPlot::SetStyle();

TString add1="1";
TString add2="2";
TString add3="3";
TString add4="4";
TString add5="5";


//Accessing all the directory in the root file
  // const char filepath = "../../../BBbar/";
  // const char fileresults = strcpy(filepath,filename);
  TFile* _file0=TFile::Open(filename);
  TDirectoryFile *dirQA=(TDirectoryFile*)(_file0->FindObjectAny(Form("%sQA%s", prefix, addon)));
  TDirectoryFile *dirEvtCuts=(TDirectoryFile*)(_file0->FindObjectAny(Form("%sEvtCuts%s", prefix, addon)));

  TDirectoryFile *dirTrackCuts=(TDirectoryFile*)(_file0->FindObjectAny(Form("%sTrackCuts%s", prefix, addon)));
  TDirectoryFile *dirv0Cuts=(TDirectoryFile*)(_file0->FindObjectAny(Form("%sv0Cuts%s", prefix, addon)));

  TDirectoryFile *dirAntiTrackCuts=(TDirectoryFile*)(_file0->FindObjectAny(Form("%sAntiTrackCuts%s", prefix, addon)));
  TDirectoryFile *dirAntiv0Cuts=(TDirectoryFile*)(_file0->FindObjectAny(Form("%sAntiv0Cuts%s", prefix, addon)));

  TDirectoryFile *dirResults=(TDirectoryFile*)(_file0->FindObjectAny(Form("%sResults%s", prefix, addon)));

  //Accessing all the Lists in the root file

  TList *listQA=0;
  TList *listEvtCuts=0;
  TList *listTrackCuts=0;
  TList *listv0Cuts=0;
  TList *listAntiTrackCuts=0;
  TList *listAntiv0Cuts=0;
  TList *listResults=0;

// QA plots
std::cout << "------------------------EVENTS QA------------------------" << std::endl;
std::cout << std::endl;

  dirQA->GetObject(Form("%sQA%s",prefix,addon),listQA);
  TList* listAliEvtCuts=(TList*)listQA->FindObject("AliEventCuts");
  TList* listPairCleaner=(TList*)listQA->FindObject("PairCleaner");
  // importing the histos
  auto* hCutStats = (TH1F*)listAliEvtCuts->FindObject("fCutStats");
  auto* cCutStats = new TCanvas("cCutStats","cCutStats",0,0,1500,1000);
  DreamPlot::SetStyleHisto(hCutStats, 8, 2, 20);
  hCutStats->Draw();

  double nTriggers = hCutStats->GetBinContent(1);
  double nEvt = hCutStats->GetBinContent(17);
// 0 flag proton from Lambda (1 antiproton from antiLambda)
  auto* hSharTracks0 = (TH1F*)listPairCleaner->FindObject("DaugthersSharedTracks_0");
  auto* hSharTracks1 = (TH1F*)listPairCleaner->FindObject("DaugthersSharedTracks_1");
  auto* hSharDaughs0 = (TH1F*)listPairCleaner->FindObject("DaugthersSharedDaughters_0");
  auto* hSharDaughs1 = (TH1F*)listPairCleaner->FindObject("DaugthersSharedDaughters_1");
  auto* cSharTracks = new TCanvas();
  cSharTracks->Divide(2,2);
// need to define a Pad in order to set a log scale
  auto* p1 = cSharTracks->cd(1);
  p1->SetLogy();
  hSharTracks0->SetTitle("; # p-V0 with shared tracks / event; Entries");
  hSharTracks0->Draw();
  auto* p2 = cSharTracks->cd(2);
  p2->SetLogy();
  hSharTracks1->SetTitle("; # #bar{p}-#bar{V0} with shared tracks/ event; Entries");
  hSharTracks1->Draw();

  auto* p3 = cSharTracks->cd(3);
  p3->SetLogy();
  hSharDaughs0->SetTitle("; # V0 pairs with shared tracks / event; Entries");
  hSharDaughs0->Draw();
  auto* p4 = cSharTracks->cd(4);
  p4->SetLogy();
  hSharDaughs1->SetTitle("; # #bar{V0} pairs with shared tracks /event; Entries");
  hSharDaughs1->Draw();

//  std::cout<< "# check (4)= " << hCutStats->GetBinContent(4)<<std::endl;
//  std::cout<< "# check= " << hCutStats->GetBinContent(5)<<std::endl;
//  std::cout<< "# check (6)= " << hCutStats->GetBinContent(6)<<std::endl;

  std::cout<< "# events before cuts = " << nTriggers <<std::endl;
  std::cout<< "# events after cuts = " << nEvt <<std::endl;

//  cCutStats->WriteOutput(Form("%sCutStats.pdf", foldername.Data()));


// Event Cuts plots
  dirEvtCuts->GetObject(Form("%sEvtCuts%s",prefix,addon),listEvtCuts);
  TList* listEvtCutsBefore=(TList*)listEvtCuts->FindObject("before");
  TList* listEvtCutsAfter=(TList*)listEvtCuts->FindObject("after");
//  std::cout << "List " << listEvtCutsBefore->GetName()<< std::endl;
  auto* hEvtCounter = (TH1F*)listEvtCuts->FindObject("EventCounter");
  auto* hCutConfig = (TH1F*)listEvtCuts->FindObject("CutConfig");
  auto* cEvtCounter = new TCanvas();
  cEvtCounter->Divide(2,1);
  cEvtCounter->cd(1);
  DreamPlot::SetStyleHisto(hEvtCounter, 8, 2, 8);
  hEvtCounter->Draw();
  cEvtCounter->cd(2);
   hCutConfig->GetXaxis()->SetRangeUser(0,13);
   hCutConfig->GetYaxis()->SetRangeUser(0.,1.1);
   hCutConfig->GetYaxis()->SetNdivisions(11,4,0,kFALSE);
  DreamPlot::SetStyleHisto(hCutConfig, 2, 2, 8);
  hCutConfig->Draw();
  auto* hSpherBefore = (TH1F*)listEvtCutsBefore->FindObject("Sphericity_before");
  auto* hSpherAfter = (TH1F*)listEvtCutsAfter->FindObject("Sphericity_after");
  auto* cSpher = new TCanvas();
  cSpher->Divide(2,1);
  DreamPlot::SetStyleHisto(hSpherBefore, 8, 2);
  DreamPlot::SetStyleHisto(hSpherAfter, 8, 2);
  cSpher->cd(1);
  hSpherBefore->Draw();
  cSpher->cd(2);
  hSpherAfter->Draw();

//Extracting the number of events w/o sphericity and with sphericity cuts
    Double_t EvtNoSpher = hSpherBefore->GetEntries();
    Double_t EvtSpher = hSpherAfter->GetEntries();
    std::cout << "# evts no sphericity cuts = " << EvtNoSpher << std::endl;
    std::cout << "# evts WITH sphericity cuts = " << EvtSpher << std::endl;


// Track Cuts plots
std::cout << "------------------------Protons------------------------" << std::endl;
std::cout << std::endl;

 dirTrackCuts->GetObject(Form("%sTrackCuts%s",prefix,addon),listTrackCuts);
  TList* listTrackCutsBefore=(TList*)listTrackCuts->FindObject("before");
  TList* listTrackCutsAfter=(TList*)listTrackCuts->FindObject("after");

  auto* hprotonDCAxyBefore = (TH2F*)listTrackCutsBefore->FindObject("DCAXY_before");
  TH1F *hprotonDCAxyBefore_projectionY = (TH1F*)hprotonDCAxyBefore->ProjectionY();

// getting min and max of x=pT and y= DCAxy values
  int binminxBefore = hprotonDCAxyBefore->FindFirstBinAbove(0,1);
  int binmaxxBefore = hprotonDCAxyBefore->FindLastBinAbove(3,1);
  int binminyBefore = hprotonDCAxyBefore->FindFirstBinAbove(-6,2);
  int binmaxyBefore = hprotonDCAxyBefore->FindLastBinAbove(5,2);
  double minxBefore = hprotonDCAxyBefore->GetXaxis()->GetXmin();
  double maxxBefore = hprotonDCAxyBefore->GetXaxis()->GetXmax();
  double minyBefore = hprotonDCAxyBefore->GetYaxis()->GetXmin();
  double maxyBefore = hprotonDCAxyBefore->GetYaxis()->GetXmin();
  int nbinsy = hprotonDCAxyBefore->GetYaxis()->GetNbins();
  int nbinsx = hprotonDCAxyBefore->GetXaxis()->GetNbins();
  std::cout << "xmin = " << minxBefore << "---" << "xmax = " << maxxBefore << std::endl;
  std::cout << "ymin = " << minyBefore << "---" << "ymax = " << maxyBefore << std::endl;
  std::cout << "number of bins in x " << nbinsx << "---" << "number of bins in y " << nbinsy << std::endl;


  auto* hprotonDCAxyCuts = (TH2F*)listTrackCuts->FindObject("DCAXYPtBinningTot");

  auto* hprotonDCAxyAfter = (TH2F*)listTrackCutsAfter->FindObject("DCAXY_after");
  TH1F *hprotonDCAxyAfter_projectionY = (TH1F*) hprotonDCAxyAfter->ProjectionY();
  auto* hprotonDCAz = (TH2F*)listTrackCutsAfter->FindObject("DCAZ_after");
  auto* hprotonPt = (TH1F*)listTrackCutsAfter->FindObject("pTDist_after");
  Double_t nProton = hprotonPt->GetEntries();
  std::cout << "# protons = " << nProton << std::endl;
  auto* hprotonPhi = (TH1F*)listTrackCutsAfter->FindObject("phiDist_after");
  auto* hprotonEta = (TH1F*)listTrackCutsAfter->FindObject("EtaDist_after");

  auto* cProtons = new TCanvas("Can_p","Can_p",0,0,6000,6000);
  cProtons->Divide(2,3);
  hprotonDCAxyBefore_projectionY->GetXaxis()->SetRangeUser(-3,3);
  hprotonDCAxyAfter_projectionY->GetXaxis()->SetRangeUser(-0.2,0.2);
  hprotonDCAxyBefore_projectionY->SetLineColor(kBlack);
  hprotonDCAxyBefore_projectionY->SetLineWidth(2);
  hprotonDCAxyAfter_projectionY->SetLineColor(kRed);
  hprotonDCAxyAfter_projectionY->SetLineWidth(2);
  hprotonDCAz->GetYaxis()->SetRangeUser(-0.5,0.5);
  hprotonEta->GetXaxis()->SetRangeUser(-1,1);
  hprotonPhi->GetYaxis()->SetRangeUser(0,1.1*hprotonPhi->GetMaximum());


    cProtons->cd(1);
   hprotonDCAxyBefore_projectionY->Draw();
   cProtons->cd(2);
   hprotonDCAxyAfter_projectionY->Draw();
  cProtons->cd(3);
  hprotonPt->Draw();
  cProtons->cd(4);
  hprotonPhi->Draw();
  cProtons->cd(6);
  hprotonEta->Draw();
  // cProtons->cd(7);
  // hpTIntprotonDCAxyBefore->Draw();

  auto* cProtonspT = new TCanvas("cProtonspT","cProtonspT",0,0,5000,5000);
  hprotonPt->GetYaxis()->SetTitle("Counts");
  hprotonPt->Draw();

  auto* cProtonsphi = new TCanvas("cProtonsphi","cProtonsphi",0,0,5000,5000);
  hprotonPhi->GetYaxis()->SetTitle("Counts");
  hprotonPhi->Draw();

  auto* cProtonsEta = new TCanvas("cProtonsEta","cProtonsEta",0,0,5000,5000);
  hprotonEta->GetYaxis()->SetTitle("Counts");
  hprotonEta->Draw();

  auto* cProtonsDCAxy = new TCanvas("cProtonsDCAxy","cProtonsDCAxy",0,0,5000,5000);
  hprotonDCAxyAfter_projectionY->GetYaxis()->SetTitle("Counts");
  hprotonDCAxyAfter_projectionY->Draw();


  std::cout << std::endl;
  std::cout << "------------------------------------------------" << std::endl;
  // Anti-Tracks plots
  std::cout << std::endl;
  std::cout << "------------------------Anti-Protons------------------------" << std::endl;

  dirAntiTrackCuts->GetObject(Form("%sAntiTrackCuts%s",prefix,addon),listAntiTrackCuts);
  TList* listAntiTrackCutsBefore=(TList*)listAntiTrackCuts->FindObject("before");
  auto* hAntiprotonDCAxyBefore = (TH2F*)listAntiTrackCutsBefore->FindObject("DCAXY_before");
  TH1F *hAntiprotonDCAxyBefore_projectionY = (TH1F*)hAntiprotonDCAxyBefore->ProjectionY();
  TList* listAntiTrackCutsAfter=(TList*)listAntiTrackCuts->FindObject("after");
  auto* hAntiprotonDCAxyAfter = (TH2F*)listAntiTrackCutsAfter->FindObject("DCAXY_after");
  TH1F *hAntiprotonDCAxyAfter_projectionY = (TH1F*)hAntiprotonDCAxyAfter->ProjectionY();

  auto* hAntiprotonDCAz = (TH2F*)listAntiTrackCutsAfter->FindObject("DCAZ_after");
  auto* hAntiprotonPt = (TH1F*)listAntiTrackCutsAfter->FindObject("pTDist_after");
  Double_t nAntiProton = hAntiprotonPt->GetEntries();
  std::cout << "# protons = " << nAntiProton << std::endl;
  auto* hAntiprotonPhi = (TH1F*)listAntiTrackCutsAfter->FindObject("phiDist_after");
  auto* hAntiprotonEta = (TH1F*)listAntiTrackCutsAfter->FindObject("EtaDist_after");

  auto* cAntiProtons = new TCanvas("Can_Antip","Can_Antip",0,0,5000,3100);
  cAntiProtons->Divide(2,3);
  hAntiprotonDCAxyAfter->GetYaxis()->SetRangeUser(-0.2,0.2);
  hAntiprotonDCAz->GetYaxis()->SetRangeUser(-0.5,0.5);
  hAntiprotonEta->GetXaxis()->SetRangeUser(-1,1);
  hAntiprotonPhi->GetYaxis()->SetRangeUser(0,1.1*hAntiprotonPhi->GetMaximum());

   // DreamPlot::SetStyleHisto(hAntiprotonDCAxy, 8, 2);
   // DreamPlot::SetStyleHisto(hAntiprotonDCAz, 8, 2);
   // DreamPlot::SetStyleHisto(hAntiprotonPt, 8, 2);
   // DreamPlot::SetStyleHisto(hAntiprotonPhi, 8, 2);
   // DreamPlot::SetStyleHisto(hAntiprotonEta, 8, 2);
   cAntiProtons->cd(1);
   hAntiprotonDCAxyBefore_projectionY->Draw();
   cAntiProtons->cd(2);
   hAntiprotonDCAxyAfter_projectionY->Draw();
  cAntiProtons->cd(3);
  hAntiprotonPt->Draw();
  cAntiProtons->cd(4);
  hAntiprotonPhi->Draw();
  cAntiProtons->cd(5);
  hAntiprotonEta->Draw();

  auto* cAntiProtonspT = new TCanvas("cAntiProtonspT","cAntiProtonspT",0,0,5000,5000);
  hAntiprotonPt->GetYaxis()->SetTitle("Counts");
  hAntiprotonPt->Draw();

  auto* cAntiProtonsphi = new TCanvas("cAntiProtonsphi","cAntiProtonsphi",0,0,5000,5000);
  hAntiprotonPhi->GetYaxis()->SetTitle("Counts");
  hAntiprotonPhi->Draw();

  auto* cAntiProtonsEta = new TCanvas("cAntiProtonsEta","cAntiProtonsEta",0,0,5000,5000);
  hAntiprotonEta->GetYaxis()->SetTitle("Counts");
  hAntiprotonEta->Draw();

  auto* cAntiProtonsDCAxy = new TCanvas("cAntiProtonsDCAxy","cAntiProtonsDCAxy",0,0,5000,5000);
  hAntiprotonDCAxyAfter_projectionY->GetYaxis()->SetTitle("Counts");
  hAntiprotonDCAxyAfter_projectionY->Draw();

//Comparison p/antip
  auto* cAntiProtonspTComp = new TCanvas("cAntiProtonspTComp","cAntiProtonspTComp",0,0,5000,5000);
  hAntiprotonPt->GetYaxis()->SetTitle("Counts");
  hAntiprotonPt->SetLineColor(kRed);
  hprotonPt->SetLineColor(kBlue);
  TLegend *leg2 = new TLegend( 0.48, 0.7, 0.88, 0.8);//x1,y1,x2,y2
leg2->SetFillStyle(0);
gStyle->SetLegendTextSize(0.02);
leg2->AddEntry(hprotonPt,"p");
leg2->AddEntry(hAntiprotonPt,"#bar{p}");
leg2->SetLineColor(0);
  hprotonPt->Draw();
  leg2->Draw("same");
  hAntiprotonPt->Draw("same");


  auto* cAntiProtonsphiComp = new TCanvas("cAntiProtonsphiComp","cAntiProtonsphiComp",0,0,5000,5000);
  hAntiprotonPhi->GetYaxis()->SetTitle("Counts");
  hAntiprotonPhi->SetLineColor(kRed);
  hprotonPhi->SetLineColor(kBlue);
  hprotonPhi->Draw();
  hAntiprotonPhi->Draw("same");

  auto* cAntiProtonsEtaComp = new TCanvas("cAntiProtonsEtaComp","cAntiProtonsEtaComp",0,0,5000,5000);
  hAntiprotonEta->GetYaxis()->SetTitle("Counts");
  hAntiprotonEta->SetLineColor(kRed);
  hprotonEta->SetLineColor(kBlue);
  hprotonEta->Draw();
  hAntiprotonEta->Draw("same");


  auto* cAntiProtonsDCAxyComp = new TCanvas("cAntiProtonsDCAxyComp","cAntiProtonsDCAxyComp",0,0,5000,5000);
  hAntiprotonDCAxyAfter_projectionY->GetYaxis()->SetTitle("Counts");
  hAntiprotonDCAxyAfter_projectionY->SetLineColor(kRed);
  hprotonDCAxyAfter_projectionY->SetLineColor(kBlue);
  hprotonDCAxyAfter_projectionY->Draw();
  hAntiprotonDCAxyAfter_projectionY->Draw("same");


std::cout << std::endl;
std::cout << "------------------------------------------------" << std::endl;
// V0 Cuts plots
std::cout << std::endl;
std::cout << "------------------------LAMBDAS------------------------" << std::endl;
std::cout << std::endl;
  dirv0Cuts->GetObject(Form("%sv0Cuts%s",prefix,addon),listv0Cuts);
  TList* listv0CutsCutsfolder=(TList*)listv0Cuts->FindObject("v0Cuts");
  TList* listv0CutsBefore=(TList*)listv0CutsCutsfolder->FindObject("before");
  TList* listv0CutsAfter=(TList*)listv0CutsCutsfolder->FindObject("after");
  TList* listv0PosCutsAfter=(TList*)listv0Cuts->FindObject("PosCuts/after");
  TList* listv0NegCuts=(TList*)listv0Cuts->FindObject("NegCuts");

  auto* hlambdaInvMass = (TH1F*)listv0CutsCutsfolder->FindObject("InvMassBefK0Rej");
  auto* hlambdaK0InvMass = (TH1F*)listv0CutsCutsfolder->FindObject("InvMassKaon");
  auto* hlambdaK0InvMassCut = (TH1F*)listv0CutsCutsfolder->FindObject("InvMasswithCuts");

   auto* hlambdaPt = (TH1F*)listv0CutsAfter->FindObject("pTDist_after");
   const float nLambdaEntries = hlambdaPt->GetEntries();
   std::cout << "# lambdas = " << nLambdaEntries << std::endl;
   auto* hlambdaPhi = (TH1F*)listv0CutsAfter->FindObject("PhiDist_after");
   auto* hlambdaEta = (TH1F*)listv0CutsAfter->FindObject("EtaDist_after");
   auto* hlambdaDecayVtxXPV = (TH1F*)listv0CutsAfter->FindObject("DecayVtxXPV_after");
   auto* hlambdaDecayVtxYPV = (TH1F*)listv0CutsAfter->FindObject("DecayVtxYPV_after");
   auto* hlambdaDecayVtxZPV = (TH1F*)listv0CutsAfter->FindObject("DecayVtxZPV_after");
   auto* hlambdaTransRadius = (TH1F*)listv0CutsAfter->FindObject("TransverseRadius_after");
   auto* hlambdaDauPDCAToPV = (TH1F*)listv0CutsAfter->FindObject("DCADauPToPV_after");
   auto* hlambdaDauNDCAToPV = (TH1F*)listv0CutsAfter->FindObject("DCADauNToPV_after");
   auto* hlambdaDauDCAToSecVtx = (TH1F*)listv0CutsAfter->FindObject("DCADauToVtx_after");
   auto* hlambdaCPA = (TH1F*)listv0CutsAfter->FindObject("PointingAngle_after");

   auto* cLambdaCuts = new TCanvas();
   cLambdaCuts->Divide(2,2);

    hlambdaPt->GetXaxis()->SetRangeUser(0.1,10);
    hlambdaPt->GetXaxis()->SetNdivisions(510,kTRUE);
    hlambdaPt->SetTitleOffset(1.5, "y");
    hlambdaPt->SetTitle("; p_{T} [GeV/c]; Entries");
    hlambdaEta->GetXaxis()->SetRangeUser(-1,1);
    hlambdaEta->SetTitle("; #eta; Entries");
    hlambdaEta->SetTitleOffset(1.5, "y");
    hlambdaPhi->SetTitle("; #phi; Entries");
    hlambdaPhi->SetMinimum(0.);
    hlambdaPhi->GetYaxis()->SetRangeUser(0,1.1*hlambdaPhi->GetMaximum());
    hlambdaPhi->SetTitleOffset(1.5, "y");

    hlambdaInvMass->GetXaxis()->SetRangeUser(1.1,1.13);
    hlambdaInvMass->SetTitle("; M_{p#pi^{-}} [GeV/c^{2}]; ");
    hlambdaK0InvMassCut->GetXaxis()->SetRangeUser(1.1,1.13);
    hlambdaK0InvMassCut->SetTitle("; M_{p#pi^{-}} [GeV/c^{2}]; ");
    hlambdaK0InvMass->GetXaxis()->SetRangeUser(0.44,0.56);
    hlambdaK0InvMass->SetTitle("; M_{#pi^{+}#pi^{-}} [GeV/c^{2}]; ");

    auto* pL01 = cLambdaCuts->cd(1);
    pL01->SetLogy();
    hlambdaPt->Draw();
    cLambdaCuts->cd(2);
    hlambdaPhi->Draw();
    cLambdaCuts->cd(3);
    hlambdaEta->Draw();
   //
   auto* cLambdaCutsMass = new TCanvas();
   cLambdaCutsMass->Divide(2,2);
   cLambdaCutsMass->cd(1);
   hlambdaK0InvMass->Draw();
   cLambdaCutsMass->cd(2);
   hlambdaInvMass->Draw();
   cLambdaCutsMass->cd(3);
   hlambdaK0InvMassCut->Draw();
   cLambdaCutsMass->cd(4);
   DreamPlot::SetStyleHisto(hlambdaK0InvMass, 8, 2, 10);
   DreamPlot::SetStyleHisto(hlambdaInvMass, 8, 4, 10);
   DreamPlot::SetStyleHisto(hlambdaK0InvMassCut, 8, 2, 10);
   hlambdaInvMass->Draw();
   hlambdaK0InvMassCut->Draw("same");
   //leg(x1,y1,x2,y2)
     auto* leg= new TLegend(0.25,0.65,0.45,0.8);
     leg->AddEntry(hlambdaInvMass, "before cuts", "l");//"l" sets the legend as lines
     leg->AddEntry(hlambdaK0InvMassCut, "K^{0} rej. and with cuts", "l");
     leg->Draw("same");


  auto* cLambdaCutsDau = new TCanvas();
  cLambdaCutsDau->Divide(2,3);
  auto* pL8 = cLambdaCutsDau->cd(1);
  pL8->SetLogy();
  hlambdaDecayVtxXPV->SetTitle("; 0.2 < r_{xy} < 100 [cm]; ");
  DreamPlot::SetStyleHisto(hlambdaDecayVtxXPV, 8, 2, 10);
  DreamPlot::SetStyleHisto(hlambdaDecayVtxYPV, 8, 4, 10);
  DreamPlot::SetStyleHisto(hlambdaDecayVtxZPV, 8, 8, 10);
  hlambdaDecayVtxXPV->Draw();
  hlambdaDecayVtxYPV->Draw("same");
  hlambdaDecayVtxZPV->Draw("same");
  auto* legL1= new TLegend(0.65,0.65,0.85,0.8);
  legL1->AddEntry(hlambdaDecayVtxXPV, "x", "l");//"l" sets the legend as lines
  legL1->AddEntry(hlambdaDecayVtxYPV, "y", "l");
  legL1->AddEntry(hlambdaDecayVtxZPV, "z", "l");
  legL1->Draw("same");
  auto* pL9 = cLambdaCutsDau->cd(2);
  pL9->SetLogy();
  hlambdaTransRadius->Draw();
  auto* pL10 = cLambdaCutsDau->cd(3);
  pL10->SetLogy();
  DreamPlot::SetStyleHisto(hlambdaDauPDCAToPV, 8, 2, 10);
  DreamPlot::SetStyleHisto(hlambdaDauNDCAToPV, 8, 4, 10);
  hlambdaDauPDCAToPV->SetTitle("; DCA to PV > 0.05 [cm]");
  hlambdaDauPDCAToPV->Draw();
  hlambdaDauNDCAToPV->Draw("same");
  auto* legL2= new TLegend(0.65,0.65,0.85,0.8);
  legL2->AddEntry(hlambdaDauPDCAToPV, "Pos. Daughter", "l");//"l" sets the legend as lines
  legL2->AddEntry(hlambdaDauNDCAToPV, "Neg. Daughter", "l");
  legL2->Draw("same");
  auto* pL11 = cLambdaCutsDau->cd(4);
  pL11->SetLogy();
  hlambdaDauDCAToSecVtx->GetXaxis()->SetRangeUser(0.,3.);
  hlambdaDauDCAToSecVtx->Draw();
  auto* pL12 = cLambdaCutsDau->cd(5);
  pL12->SetLogy();
  hlambdaCPA->GetXaxis()->SetRangeUser(0.96,1.1);
  hlambdaCPA->Draw();



  std::cout << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << std::endl;
  std::cout<<"-------------------Extracting S/B for Lambda-----------------"<<std::endl;
  std::cout << std::endl;


  const float marginLambda = 0.004;
  const float massLambda = 1.115;
  float nLambda = 0.;
  float purLambda = 0.;


  auto* cLambda = new TCanvas();
  auto* hRecoLambdaM = (TH1F*)hlambdaK0InvMassCut->Clone();
  hRecoLambdaM->Draw("PE");
  float lambdaSignalAll, lambdaSignalAllErr, lambdaBackgroundAll, lambdaBackgroundAllErr;
  TList *funListLambda;
  FitLambda(hRecoLambdaM, lambdaSignalAll, lambdaSignalAllErr, lambdaBackgroundAll, lambdaBackgroundAllErr, massLambda-marginLambda, massLambda+marginLambda);
  std::cout << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  std::cout << "Lambda \n";
  std::cout << "Signal " << lambdaSignalAll << " Background " << lambdaBackgroundAll << " S/B " << lambdaSignalAll/lambdaBackgroundAll << " Purity " << lambdaSignalAll/(lambdaSignalAll+lambdaBackgroundAll)*100.f << "\n";
  std::cout << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  nLambda = lambdaSignalAll;
  purLambda = lambdaSignalAll/(lambdaSignalAll+lambdaBackgroundAll)*100.f;
  hRecoLambdaM->GetXaxis()->SetRangeUser(1.1, 1.13);
  funListLambda = hRecoLambdaM->GetListOfFunctions();
  TF1 *fLambdaTotal = (TF1*)funListLambda->FindObject("fLambda");

  float amp1 = fLambdaTotal->GetParameter(3);
  float amp2 = fLambdaTotal->GetParameter(6);
  float mean1 = fLambdaTotal->GetParameter(4);
  float mean2 = fLambdaTotal->GetParameter(7);
  float width1 = fLambdaTotal->GetParameter(5);
  float width2 = fLambdaTotal->GetParameter(8);

  float meanMass = weightedMean(amp1, mean1, amp2, mean2);
  float meanWidth = weightedMean(amp1, width1, amp2, width2);
  TLatex LambdaLabel;
  LambdaLabel.SetNDC(kTRUE);
  LambdaLabel.SetTextSize(gStyle->GetTextSize()*0.8);
  LambdaLabel.DrawLatex(gPad->GetUxmax()-0.8, gPad->GetUymax()-0.39,
  Form("#splitline{#splitline{#splitline{S/B: %.1f}{m_{#Lambda} = %.1f (MeV/#it{c}^{2})}}{#sigma_{#Lambda} = %.1f (MeV/#it{c}^{2})}}{Purity = %.1f %%}", lambdaSignalAll/lambdaBackgroundAll, meanMass*1000.f, meanWidth*1000.f, lambdaSignalAll/(lambdaSignalAll+lambdaBackgroundAll)*100.f));
  TLatex BeamTextLambda;
  BeamTextLambda.SetNDC(kTRUE);
  BeamTextLambda.DrawLatex(gPad->GetUxmax()-0.8, gPad->GetUymax()-0.2, "pp #sqrt{#it{s}} = 13 TeV");

  std::cout << std::endl;
  std::cout << "------------------------------------------------" << std::endl;
  // Anti-V0 Cuts plots
  std::cout << std::endl;
  std::cout << "------------------------Anti-Lambdas------------------------" << std::endl;
  std::cout << std::endl;
    dirAntiv0Cuts->GetObject(Form("%sAntiv0Cuts%s",prefix,addon),listAntiv0Cuts);
    TList* listAntiv0CutsCutsfolder=(TList*)listAntiv0Cuts->FindObject("v0Cuts");
    TList* listAntiv0CutsBefore=(TList*)listAntiv0CutsCutsfolder->FindObject("before");
    TList* listAntiv0CutsAfter=(TList*)listAntiv0CutsCutsfolder->FindObject("after");
    TList* listAntiv0PosCutsAfter=(TList*)listAntiv0Cuts->FindObject("PosCuts/after");
    TList* listAntiv0NegCuts=(TList*)listAntiv0Cuts->FindObject("NegCuts/after");

    auto* hAntilambdaInvMass = (TH1F*)listAntiv0CutsCutsfolder->FindObject("InvMassBefK0Rej");
    auto* hAntilambdaK0InvMass = (TH1F*)listAntiv0CutsCutsfolder->FindObject("InvMassKaon");
    auto* hAntilambdaK0InvMassCut = (TH1F*)listAntiv0CutsCutsfolder->FindObject("InvMasswithCuts");

     auto* hAntilambdaPt = (TH1F*)listAntiv0CutsAfter->FindObject("pTDist_after");
     const float nAntiLambdaEntries = hAntilambdaPt->GetEntries();
     std::cout << "# Antilambdas = " << nAntiLambdaEntries << std::endl;
     auto* hAntilambdaPhi = (TH1F*)listAntiv0CutsAfter->FindObject("PhiDist_after");
     auto* hAntilambdaEta = (TH1F*)listAntiv0CutsAfter->FindObject("EtaDist_after");
     auto* hAntilambdaDecayVtxXPV = (TH1F*)listAntiv0CutsAfter->FindObject("DecayVtxXPV_after");
     auto* hAntilambdaDecayVtxYPV = (TH1F*)listAntiv0CutsAfter->FindObject("DecayVtxYPV_after");
     auto* hAntilambdaDecayVtxZPV = (TH1F*)listAntiv0CutsAfter->FindObject("DecayVtxZPV_after");
     auto* hAntilambdaTransRadius = (TH1F*)listAntiv0CutsAfter->FindObject("TransverseRadius_after");
     auto* hAntilambdaDauPDCAToPV = (TH1F*)listAntiv0CutsAfter->FindObject("DCADauPToPV_after");
     auto* hAntilambdaDauNDCAToPV = (TH1F*)listAntiv0CutsAfter->FindObject("DCADauNToPV_after");
     auto* hAntilambdaDauDCAToSecVtx = (TH1F*)listAntiv0CutsAfter->FindObject("DCADauToVtx_after");
     auto* hAntilambdaCPA = (TH1F*)listAntiv0CutsAfter->FindObject("PointingAngle_after");

     auto* cAntiLambdaCuts = new TCanvas();

     cAntiLambdaCuts->Divide(2,2);
      hAntilambdaPt->GetXaxis()->SetRangeUser(0.1,6.9);
      hAntilambdaPt->GetXaxis()->SetNdivisions(510,kTRUE);
      hAntilambdaPt->SetTitleOffset(1.5, "y");
      hAntilambdaPt->SetTitle("; p_{T} [GeV/c]; Entries");
      hAntilambdaEta->GetXaxis()->SetRangeUser(-1,1);
      hAntilambdaEta->SetTitle("; #eta; Entries");
      hAntilambdaEta->SetTitleOffset(1.5, "y");
      hAntilambdaPhi->SetTitle("; #phi; Entries");
      hAntilambdaPhi->SetMinimum(0.);
      hAntilambdaPhi->GetYaxis()->SetRangeUser(0,1.1*hAntilambdaPhi->GetMaximum());
      hAntilambdaPhi->SetTitleOffset(1.5, "y");

      hAntilambdaInvMass->GetXaxis()->SetRangeUser(1.1,1.13);
      hAntilambdaInvMass->SetTitle("; M_{#bar{p}#pi^{+}} [GeV/c^{2}]; ");
      hAntilambdaK0InvMassCut->GetXaxis()->SetRangeUser(1.1,1.13);
      hAntilambdaK0InvMassCut->SetTitle("; M_{#bar{p}#pi^{+}} [GeV/c^{2}]; ");
      hAntilambdaK0InvMass->GetXaxis()->SetRangeUser(0.44,0.56);
      hAntilambdaK0InvMass->SetTitle("; M_{#pi^{-}#pi^{+}} [GeV/c^{2}]; ");

      auto* pAL01 = cAntiLambdaCuts->cd(1);
      pAL01->SetLogy();
      cAntiLambdaCuts->cd(1);
      hAntilambdaPt->Draw();
      cAntiLambdaCuts->cd(2);
      hAntilambdaPhi->Draw();
      cAntiLambdaCuts->cd(3);
      hAntilambdaEta->Draw();
     //
     auto* cAntiLambdaCutsMass = new TCanvas();
     cAntiLambdaCutsMass->Divide(2,2);
     cAntiLambdaCutsMass->cd(1);
     hAntilambdaK0InvMass->Draw();
     cAntiLambdaCutsMass->cd(2);
     hAntilambdaInvMass->Draw();
     cAntiLambdaCutsMass->cd(3);
     hAntilambdaK0InvMassCut->Draw();
     cAntiLambdaCutsMass->cd(4);
     DreamPlot::SetStyleHisto(hAntilambdaInvMass, 8, 2, 10);
     DreamPlot::SetStyleHisto(hAntilambdaK0InvMassCut, 8, 4, 10);
     DreamPlot::SetStyleHisto(hAntilambdaK0InvMass, 8, 4, 10);
     hAntilambdaInvMass->Draw();
     hAntilambdaK0InvMassCut->Draw("same");
     //leg(x1,y1,x2,y2)
       auto* Antileg= new TLegend(0.25,0.65,0.45,0.8);
       Antileg->AddEntry(hAntilambdaInvMass, "before cuts", "l");//"l" sets the legend as lines
       Antileg->AddEntry(hAntilambdaK0InvMassCut, "#bar{K^{0}} rej. and with cuts", "l");
       Antileg->Draw("same");


    auto* cAntiLambdaCutsDau = new TCanvas();
    cAntiLambdaCutsDau->Divide(2,3);
    auto* pAntiL8 = cAntiLambdaCutsDau->cd(1);
    pAntiL8->SetLogy();
    hAntilambdaDecayVtxXPV->SetTitle("; 0.2 < r_{xy} < 100 [cm]; ");
    DreamPlot::SetStyleHisto(hAntilambdaDecayVtxXPV, 8, 2, 10);
    DreamPlot::SetStyleHisto(hAntilambdaDecayVtxYPV, 8, 4, 10);
    DreamPlot::SetStyleHisto(hAntilambdaDecayVtxZPV, 8, 8, 10);

    hAntilambdaDecayVtxXPV->Draw();
    hAntilambdaDecayVtxYPV->Draw("same");
    hAntilambdaDecayVtxZPV->Draw("same");
    auto* legAntiL1= new TLegend(0.65,0.65,0.85,0.8);
    legAntiL1->AddEntry(hAntilambdaDecayVtxXPV, "x", "l");//"l" sets the legend as lines
    legAntiL1->AddEntry(hAntilambdaDecayVtxYPV, "y", "l");
    legAntiL1->AddEntry(hAntilambdaDecayVtxZPV, "z", "l");
    legAntiL1->Draw("same");
    auto* pAntiL9 = cAntiLambdaCutsDau->cd(2);
    pAntiL9->SetLogy();
    hAntilambdaTransRadius->Draw();
    auto* pAntiL10 = cAntiLambdaCutsDau->cd(3);
    pAntiL10->SetLogy();
    DreamPlot::SetStyleHisto(hAntilambdaDauPDCAToPV, 8, 2, 10);
    DreamPlot::SetStyleHisto(hAntilambdaDauNDCAToPV, 8, 4, 10);
    hAntilambdaDauPDCAToPV->SetTitle("; DCA to PV > 0.05 [cm]");
    hAntilambdaDauPDCAToPV->Draw();
    hAntilambdaDauNDCAToPV->Draw("same");
    auto* legAntiL2= new TLegend(0.65,0.65,0.85,0.8);
    legAntiL2->AddEntry(hAntilambdaDauPDCAToPV, "Pos. Daughter", "l");//"l" sets the legend as lines
    legAntiL2->AddEntry(hAntilambdaDauNDCAToPV, "Neg. Daughter", "l");
    legAntiL2->Draw("same");
    auto* pAntiL11 = cAntiLambdaCutsDau->cd(4);
    pAntiL11->SetLogy();
    hAntilambdaDauDCAToSecVtx->GetXaxis()->SetRangeUser(0.,3.);
    hAntilambdaDauDCAToSecVtx->Draw();
    auto* pAntiL12 = cAntiLambdaCutsDau->cd(5);
    pAntiL12->SetLogy();
    hAntilambdaCPA->GetXaxis()->SetRangeUser(0.96,1.1);
    hAntilambdaCPA->Draw();

    std::cout << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    std::cout << std::endl;
    std::cout<<"-------------------Extracting S/B for AntiLambda-----------------"<<std::endl;
    std::cout << std::endl;

    float nAntiLambda = 0.;
    float purAntiLambda = 0.;


    auto* cAntiLambda = new TCanvas();
    auto* hRecoAntiLambdaM = (TH1F*)hAntilambdaK0InvMassCut->Clone();
    hRecoAntiLambdaM->Draw("PE");
    FitLambda(hRecoAntiLambdaM, lambdaSignalAll, lambdaSignalAllErr, lambdaBackgroundAll, lambdaBackgroundAllErr, massLambda-marginLambda, massLambda+marginLambda);
    std::cout << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    std::cout << "AntiLambda \n";
    std::cout << "Signal " << lambdaSignalAll << " Background " << lambdaBackgroundAll << " S/B " << lambdaSignalAll/lambdaBackgroundAll << " Purity " << lambdaSignalAll/(lambdaSignalAll+lambdaBackgroundAll)*100.f << "\n";
    std::cout << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    nAntiLambda = lambdaSignalAll;
    purAntiLambda = lambdaSignalAll/(lambdaSignalAll+lambdaBackgroundAll)*100.f;
    hRecoAntiLambdaM->GetXaxis()->SetRangeUser(1.1, 1.13);
    funListLambda = hRecoAntiLambdaM->GetListOfFunctions();
    fLambdaTotal = (TF1*)funListLambda->FindObject("fLambda");

    amp1 = fLambdaTotal->GetParameter(3);
    amp2 = fLambdaTotal->GetParameter(6);
    mean1 = fLambdaTotal->GetParameter(4);
    mean2 = fLambdaTotal->GetParameter(7);
    width1 = fLambdaTotal->GetParameter(5);
    width2 = fLambdaTotal->GetParameter(8);

    meanMass = weightedMean(amp1, mean1, amp2, mean2);
    meanWidth = weightedMean(amp1, width1, amp2, width2);
    LambdaLabel.DrawLatex(gPad->GetUxmax()-0.8, gPad->GetUymax()-0.39,
    Form("#splitline{#splitline{#splitline{S/B: %.1f}{m_{#bar{#Lambda}} = %.1f (MeV/#it{c}^{2})}}{#sigma_{#bar{#Lambda}} = %.1f (MeV/#it{c}^{2})}}{Purity = %.1f %%}", lambdaSignalAll/lambdaBackgroundAll, meanMass*1000.f, meanWidth*1000.f, lambdaSignalAll/(lambdaSignalAll+lambdaBackgroundAll)*100.f));
    BeamTextLambda.DrawLatex(gPad->GetUxmax()-0.8, gPad->GetUymax()-0.2, "pp #sqrt{#it{s}} = 13 TeV");
  std::cout << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;


  std::cout << std::endl;
  std::cout << "------------------------------------------------------" << std::endl;
  auto* LambdapTComp = new TCanvas();
  hlambdaPt->SetLineColor(kBlue);
  hAntilambdaPt->SetLineColor(kRed);
  hlambdaPt->Draw();
  hAntilambdaPt->Draw("same");
  TLegend *leg2L = new TLegend( 0.48, 0.7, 0.88, 0.8);//x1,y1,x2,y2
leg2L->SetFillStyle(0);
gStyle->SetLegendTextSize(0.02);
leg2L->AddEntry(hlambdaPt,"#Lambda");
leg2L->AddEntry(hAntilambdaPt,"#bar{#Lambda}");
leg2L->SetLineColor(0);
leg2L->Draw("same");


auto* LambdaPhiComp = new TCanvas();
hlambdaPhi->SetLineColor(kBlue);
hAntilambdaPhi->SetLineColor(kRed);
hlambdaPhi->Draw();
hAntilambdaPhi->Draw("same");

auto* LambdaEtaComp = new TCanvas();
hlambdaEta->SetLineColor(kBlue);
hAntilambdaEta->SetLineColor(kRed);
hlambdaEta->Draw();
hAntilambdaEta->Draw("same");

  std::cout << std::endl;
  std::cout<<"The system only dreams in total darkness"<<std::endl;
  gSystem->Exec(TString::Format("mkdir /Users/Valentina/cernbox/Analysis/BBbar/GentleFemto_Output/QA/QA_%s_%s/",prefix,date));
  TString foldername = TString::Format("/Users/Valentina/cernbox/Analysis/BBbar/GentleFemto_Output/QA/QA_%s_%s/",prefix,date);


  if(strcmp(addon, add1)==0){
    cCutStats->SaveAs(foldername+"CutStats_st1.pdf");
    cSharTracks->SaveAs(foldername+"DaughSharTracks_st1.pdf");
    cEvtCounter->SaveAs(foldername+"EventCuts_st1.pdf");
    cSpher->SaveAs(foldername+"SphericityDist_st1.pdf");
    cProtons->SaveAs(foldername+"ProtonsQA_st1.pdf");
    cAntiProtons->SaveAs(foldername+"AntiProtonsQA_st1.pdf");
    cLambdaCuts->SaveAs(foldername+"LambdaQASel_st1.pdf");
    cLambdaCutsMass->SaveAs(foldername+"LambdaQAMass_st1.pdf");
    cLambdaCutsDau->SaveAs(foldername+"LambdaQADau_st1.pdf");
    cLambda->SaveAs(foldername+"SignalLambda_st1.pdf");
    cAntiLambdaCuts->SaveAs(foldername+"AntiLambdaQASel_st1.pdf");
    cAntiLambdaCutsMass->SaveAs(foldername+"AntiLambdaQAMass_st1.pdf");
    cAntiLambdaCutsDau->SaveAs(foldername+"AntiLambdaQADau_st1.pdf");
    cAntiLambda->SaveAs(foldername+"AntiSignalLambda_st1.pdf");
}
else if(strcmp(addon, add2)==0){
  cCutStats->SaveAs(foldername+"CutStats_st2.pdf");
  cSharTracks->SaveAs(foldername+"DaughSharTracks_st2.pdf");
  cEvtCounter->SaveAs(foldername+"EventCuts_st2.pdf");
  cSpher->SaveAs(foldername+"SphericityDist_st2.pdf");
  cProtons->SaveAs(foldername+"ProtonsQA_st2.pdf");
  cAntiProtons->SaveAs(foldername+"AntiProtonsQA_st2.pdf");
  cLambdaCuts->SaveAs(foldername+"LambdaQASel_st2.pdf");
  cLambdaCutsMass->SaveAs(foldername+"LambdaQAMass_st2.pdf");
  cLambdaCutsDau->SaveAs(foldername+"LambdaQADau_st2.pdf");
  cLambda->SaveAs(foldername+"SignalLambda_st2.pdf");
  cAntiLambdaCuts->SaveAs(foldername+"AntiLambdaQASel_st2.pdf");
  cAntiLambdaCutsMass->SaveAs(foldername+"AntiLambdaQAMass_st2.pdf");
  cAntiLambdaCutsDau->SaveAs(foldername+"AntiLambdaQADau_st2.pdf");
  cAntiLambda->SaveAs(foldername+"AntiSignalLambda_st2.pdf");
}
else if(strcmp(addon, add3)==0){
  cCutStats->SaveAs(foldername+"CutStats_st3.pdf");
  cSharTracks->SaveAs(foldername+"DaughSharTracks_st3.pdf");
  cEvtCounter->SaveAs(foldername+"EventCuts_st3.pdf");
  cSpher->SaveAs(foldername+"SphericityDist_st3.pdf");
  cProtons->SaveAs(foldername+"ProtonsQA_st3.pdf");
  cAntiProtons->SaveAs(foldername+"AntiProtonsQA_st3.pdf");
  cLambdaCuts->SaveAs(foldername+"LambdaQASel_st3.pdf");
  cLambdaCutsMass->SaveAs(foldername+"LambdaQAMass_st3.pdf");
  cLambdaCutsDau->SaveAs(foldername+"LambdaQADau_st3.pdf");
  cLambda->SaveAs(foldername+"SignalLambda_st3.pdf");
  cAntiLambdaCuts->SaveAs(foldername+"AntiLambdaQASel_st3.pdf");
  cAntiLambdaCutsMass->SaveAs(foldername+"AntiLambdaQAMass_st3.pdf");
  cAntiLambdaCutsDau->SaveAs(foldername+"AntiLambdaQADau_st3.pdf");
  cAntiLambda->SaveAs(foldername+"AntiSignalLambda_st3.pdf");
}
else if(strcmp(addon, add4)==0){
  cCutStats->SaveAs(foldername+"CutStats_full.pdf");
  cSharTracks->SaveAs(foldername+"DaughSharTracks_full.pdf");
  cEvtCounter->SaveAs(foldername+"EventCuts_full.pdf");
  cSpher->SaveAs(foldername+"SphericityDist_full.pdf");
  cProtons->SaveAs(foldername+"ProtonsQA_full.pdf");
  cProtonspT->SaveAs(foldername+"ProtonspT_full.pdf");
  cProtonsphi->SaveAs(foldername+"ProtonsPhi_full.pdf");
  cProtonsEta->SaveAs(foldername+"ProtonsEta_full.pdf");
  cProtonsDCAxy->SaveAs(foldername+"ProtonsDCAxy_full.pdf");
  cAntiProtons->SaveAs(foldername+"AntiProtonsQA_full.pdf");
  cAntiProtonspT->SaveAs(foldername+"AntiProtonspT_full.pdf");
  cAntiProtonsphi->SaveAs(foldername+"AntiProtonsPhi_full.pdf");
  cAntiProtonsEta->SaveAs(foldername+"AntiProtonsEta_full.pdf");
  cAntiProtonsDCAxy->SaveAs(foldername+"AntiProtonsDCAxy_full.pdf");

  cAntiProtonspTComp->SaveAs(foldername+"Comp_pT_full.pdf");
  cAntiProtonsphiComp->SaveAs(foldername+"Comp_Phi_full.pdf");
  cAntiProtonsEtaComp->SaveAs(foldername+"Comp_Eta_full.pdf");
  cAntiProtonsDCAxyComp->SaveAs(foldername+"Comp_DCAxy_full.pdf");

  cLambdaCuts->SaveAs(foldername+"LambdaQASel_full.pdf");
  cLambdaCutsMass->SaveAs(foldername+"LambdaQAMass_full.pdf");
  cLambdaCutsDau->SaveAs(foldername+"LambdaQADau_full.pdf");
  cLambda->SaveAs(foldername+"SignalLambda_full.pdf");
  cAntiLambdaCuts->SaveAs(foldername+"AntiLambdaQASel_full.pdf");
  cAntiLambdaCutsMass->SaveAs(foldername+"AntiLambdaQAMass_full.pdf");
  cAntiLambdaCutsDau->SaveAs(foldername+"AntiLambdaQADau_full.pdf");
  cAntiLambda->SaveAs(foldername+"AntiSignalLambda_full.pdf");

  LambdapTComp->SaveAs(foldername+"LambdaComp_pT_full.pdf");
  LambdaPhiComp->SaveAs(foldername+"LambdaComp_Phi_full.pdf");
  LambdaEtaComp->SaveAs(foldername+"LambdaComp_Eta_full.pdf");
}

//  TList *Results;
//  dirResults->GetObject(Form("%sResults%s", prefix, addon),Results);
}

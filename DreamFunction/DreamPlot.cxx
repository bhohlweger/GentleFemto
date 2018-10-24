/*
 * DreamPlot.cxx
 *
 *  Created on: 29 Aug 2018
 *      Author: bernhardhohlweger
 */

#include "DreamPlot.h"
#include "TStyle.h"
#include "TFile.h"
#include "TLatex.h"
DreamPlot::DreamPlot()
    : fProtonProton(nullptr),
      fProtonLambda(nullptr),
      fLambdaLambda(nullptr),
      fProtonXi(nullptr),
      fRadius(0),
      fRadiusStat(0),
      fRadiusSysUp(0),
      fRadiusSysLow(0),
      fEnergy(0),
      fCollisionSystem(""),
      fMonteCarloGen("") {
  fProtonProton = new DreamData("ProtonProton");
  fProtonLambda = new DreamData("ProtonLambda");
  fLambdaLambda = new DreamData("LambdaLambda");
  fProtonXi = new DreamData("ProtonXi");
}

DreamPlot::~DreamPlot() {
  // TODO Auto-generated destructor stub
}

void DreamPlot::ReadData(const char* PathToDataFolder,
                         const char* PathToSysFolder, int binWidth) {
  TString HistName = "hCk_Reweighted_";
  if (binWidth == 16) {
    HistName += "0";
  } else if (binWidth == 20) {
    HistName += "1";
  } else {
    std::cout << "Unknown bin width " << binWidth << std::endl;
  }
  TFile* CFFile_pp = TFile::Open(Form("%sCFOutput_pp.root", PathToDataFolder));
  TFile* CFFile_ppSys = TFile::Open(
      Form("%sC2totalsysPP.root", PathToSysFolder));
  fProtonProton->SetCorrelationFunction(
      (TH1F*) CFFile_pp->Get("hCk_Reweighted_0"));
  fProtonProton->SetSystematics((TF1*) CFFile_ppSys->Get("RelSysPPUnbinned"), 1,
                                0.002);

  TFile* CFFile_pL = TFile::Open(Form("%sCFOutput_pL.root", PathToDataFolder));
  TFile* CFFile_pLSys = TFile::Open(
      Form("%sC2totalsysPL.root", PathToSysFolder));
  fProtonLambda->SetCorrelationFunction(
      (TH1F*) CFFile_pL->Get(HistName.Data()));
  fProtonLambda->SetSystematics((TF1*) CFFile_pLSys->Get("RelSysPLUnbinned"), 1,
                                0.002);

  TFile* CFFile_LL = TFile::Open(Form("%sCFOutput_LL.root", PathToDataFolder));
  TFile* CFFile_LLSys = TFile::Open(
      Form("%sC2totalsysLL.root", PathToSysFolder));
  fLambdaLambda->SetCorrelationFunction(
      (TH1F*) CFFile_LL->Get(HistName.Data()));
  fLambdaLambda->SetSystematics((TF1*) CFFile_LLSys->Get("RelSysLLUnbinned"), 1,
                                0.002);

  TFile* CFFile_pXi = TFile::Open(
      Form("%sCFOutput_pXi.root", PathToDataFolder));
  TFile* CFFile_pXiSys = TFile::Open(
      Form("%sC2totalsysPXi.root", PathToSysFolder));
  fProtonXi->SetCorrelationFunction((TH1F*) CFFile_pXi->Get(HistName.Data()));
  fProtonXi->SetSystematics((TF1*) CFFile_pXiSys->Get("RelSysPXiUnbinned"), 1,
                            0.002);

  return;
}

void DreamPlot::ReadSimulation(const char* PathToSimFolder, int binWidth) {
  TString HistName = "hCk_Reweighted_";
  if (binWidth == 16) {
    HistName += "0";
  } else if (binWidth == 20) {
    HistName += "1";
  } else {
    std::cout << "Unknown bin width " << binWidth << std::endl;
  }
  TFile* CFFile_pp = TFile::Open(Form("%sCFOutput_pp.root", PathToSimFolder));
  fProtonProton->SetCorrelationFunctionSimulation(
      (TH1F*) CFFile_pp->Get("hCk_Reweighted_0"));

  TFile* CFFile_pL = TFile::Open(Form("%sCFOutput_pL.root", PathToSimFolder));
  fProtonLambda->SetCorrelationFunctionSimulation(
      (TH1F*) CFFile_pL->Get(HistName.Data()));

  TFile* CFFile_LL = TFile::Open(Form("%sCFOutput_LL.root", PathToSimFolder));
  fLambdaLambda->SetCorrelationFunctionSimulation(
      (TH1F*) CFFile_LL->Get(HistName.Data()));

  TFile* CFFile_pXi = TFile::Open(Form("%sCFOutput_pXi.root", PathToSimFolder));
  fProtonXi->SetCorrelationFunctionSimulation(
      (TH1F*) CFFile_pXi->Get(HistName.Data()));
}

void DreamPlot::ReadFit(const char* fitPath) {
  TString PathToFile = Form("%sSYSTEMATICS_CutVarAdd_Global_Radius_Normal.root",
                            fitPath);
  TFile *systFit = TFile::Open(PathToFile.Data());
  if (systFit) {
    TGraph* grpp_default = (TGraph*) systFit->Get("ppGraphDefault");
    TGraph* grppLow = (TGraph*) systFit->Get("ppGraphLowerLim");
    TGraph* grppUp = (TGraph*) systFit->Get("ppGraphUpperLim");
    if (!grpp_default) {
      std::cout << "no pp Default file \n";
    } else if (!grppLow) {
      std::cout << "no pp lower file \n";
    } else if (!grppUp) {
      std::cout << "no pp upper file \n";
    } else {
      fProtonProton->FemtoModelFitBands(grpp_default, grppLow, grppUp, 1000, 2);
    }
    TGraph* gr_pLNLO_default = (TGraph*) systFit->Get("pLamGraphDefault_NLO");
    TGraph* gr_pLNLO_low = (TGraph*) systFit->Get("Copy_pLamGraphLowerLim_NLO");
    TGraph* gr_pLNLO_up = (TGraph*) systFit->Get("pLamGraphUpperLim_NLO");
    if (!gr_pLNLO_default) {
      std::cout << "no pL NLO Default file \n";
    } else if (!gr_pLNLO_low) {
      std::cout << "no pL NLO lower file \n";
    } else if (!gr_pLNLO_up) {
      std::cout << "no pL NLO upper file \n";
    } else {
      fProtonLambda->FemtoModelFitBands(gr_pLNLO_default, gr_pLNLO_low,
                                        gr_pLNLO_up, 1000, 1);
    }
    TGraph* gr_pLLO_default = (TGraph*) systFit->Get("pLamGraphDefault");
    TGraph* grpLLowLO = (TGraph*) systFit->Get("pLamGraphLowerLim");
    TGraph* grpLUpLO = (TGraph*) systFit->Get("pLamGraphUpperLim");
    if (!gr_pLLO_default) {
      std::cout << "no pL LO Default file \n";
    } else if (!grpLLowLO) {
      std::cout << "no pL LO lower file \n";
    } else if (!grpLUpLO) {
      std::cout << "no pL LO upper file \n";
    } else {
      fProtonLambda->FemtoModelFitBands(gr_pLLO_default, grpLLowLO, grpLUpLO,
                                        1000, 3);
    }
    TGraph* grLLDefault = (TGraph*) systFit->Get("LamLamGraphDefault");
    TGraph* grLLLow = (TGraph*) systFit->Get("LamLamGraphLowerLim");
    TGraph* grLLUp = (TGraph*) systFit->Get("LamLamGraphUpperLim");
    if (!grLLDefault) {
      std::cout << "no LL Default file \n";
    } else if (!grLLLow) {
      std::cout << "no LL lower file \n";
    } else if (!grLLUp) {
      std::cout << "no LL upper file \n";
    } else {
      fLambdaLambda->FemtoModelFitBands(grLLDefault, grLLLow, grLLUp, 1000, 5);
    }
    TGraph* grpXiDefault = (TGraph*) systFit->Get("pXimGraphDefault");
    TGraph* grpXiLower = (TGraph*) systFit->Get("pXimGraphLowerLim");
    TGraph* grpXiUpper = (TGraph*) systFit->Get("pXimGraphUpperLim");
    if (!grpXiDefault) {
      std::cout << "no pXi Default file \n";
    } else if (!grpXiLower) {
      std::cout << "no pXi lower file \n";
    } else if (!grpXiUpper) {
      std::cout << "no pXi upper file \n";
    } else {
      fProtonXi->FemtoModelFitBands(grpXiDefault, grpXiLower, grpXiUpper, 1000,
                                    6);
    }
    TGraph* grpXiDefaultCoulomb = (TGraph*) systFit->Get(
        "pXimGraphDefault_COULOMB");
    TGraph* grpXiLowerCoulomb = (TGraph*) systFit->Get(
        "pXimGraphLowerLim_COULOMB");
    TGraph* grpXiUpperCoulomb = (TGraph*) systFit->Get(
        "pXimGraphUpperLim_COULOMB");
    if (!grpXiDefaultCoulomb) {
      std::cout << "no pXi Coulomb Default file \n";
    } else if (!grpXiLowerCoulomb) {
      std::cout << "no pXi Coulomb lower file \n";
    } else if (!grpXiUpperCoulomb) {
      std::cout << "no pXi Coulomb upper file \n";
    } else {
      fProtonXi->FemtoModelFitBands(grpXiDefaultCoulomb, grpXiLowerCoulomb,
                                    grpXiUpperCoulomb, 1000, 7);
    }
  } else {
    std::cout << "No Cats file!  \n";
  }
  return;
}

void DreamPlot::SetStyle(bool graypalette, bool title) {
  const int NCont = 255;
  gStyle->Reset("Plain");
  gStyle->SetNumberContours(NCont);
  gStyle->SetOptTitle(title);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetOptStat(0);
  if (graypalette)
    gStyle->SetPalette(8, 0);
  else
    gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(2);
  gStyle->SetLabelSize(0.045, "xyz");
  gStyle->SetLabelOffset(0.01, "y");
  gStyle->SetLabelOffset(0.01, "x");
  gStyle->SetLabelColor(kBlack, "xyz");
  gStyle->SetTitleSize(0.05, "xyz");
  gStyle->SetTitleOffset(1.25, "y");
  gStyle->SetTitleOffset(1.2, "x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetLegendFont(42);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetErrorX(0.001);
  const int NRGBs = 6;
  Double_t stops[NRGBs];
  for (int i = 0; i < NRGBs; ++i)
    stops[i] = float(i) / (NRGBs - 1);

  Double_t red[NRGBs] = { 1., 29. / 255., 25. / 255., 27. / 255., 32. / 255.,
      24. / 255. };
  Double_t green[NRGBs] = { 1., 221. / 255., 160. / 255., 113. / 255., 74.
      / 255., 37. / 255. };
  Double_t blue[NRGBs] = { 1., 221. / 255., 184. / 255., 154. / 255., 129.
      / 255., 98. / 255. };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
}

void DreamPlot::SetStyleHisto(TH1 *histo, int marker, int color) {
  histo->GetXaxis()->SetLabelSize(0.045);
  histo->GetXaxis()->SetTitleSize(0.05);
  histo->GetXaxis()->SetLabelOffset(0.01);
  histo->GetXaxis()->SetTitleOffset(1.2);
  histo->GetXaxis()->SetLabelFont(42);
  histo->GetYaxis()->SetLabelSize(0.045);
  histo->GetYaxis()->SetTitleSize(0.05);
  histo->GetYaxis()->SetLabelOffset(0.01);
  histo->GetYaxis()->SetTitleOffset(1.25);
  histo->SetMarkerStyle(marker);
  histo->SetMarkerColor(color);
  histo->SetLineColor(color);
}


void DreamPlot::DrawCorrelationFunctions() {
  SetStyle();
  const float right = 0.025;
  const float top = 0.025;
  TLatex ref;
  ref.SetTextSize(gStyle->GetTextSize() * 0.7);
  ref.SetNDC(kTRUE);
  TCanvas* c_PP = new TCanvas("CFpp", "CFpp", 0, 0, 650, 550);
  c_PP->SetRightMargin(right);
  c_PP->SetTopMargin(top);
  fProtonProton->SetLegendName("p-p #oplus #bar{p}-#bar{p} pairs");
  fProtonProton->SetLegendName("Femtoscopic fit");
  fProtonProton->SetRangePlotting(0, 0.125, 0.5, 3.5);
  fProtonProton->DrawCorrelationPlot(c_PP);
  DrawSystemInfo(c_PP);
  c_PP->SaveAs("CF_pp_Gauss_prelim.pdf");

  TCanvas* c_PL = new TCanvas("CFpL", "CFpL", 0, 0, 650, 550);
  c_PL->SetRightMargin(right);
  c_PL->SetTopMargin(top);
  fProtonLambda->SetLegendName("p-#Lambda #oplus #bar{p}-#bar{#Lambda} pairs");
  fProtonLambda->SetLegendName("Femtoscopic fit (#chiEFT NLO)");
  fProtonLambda->SetLegendName("Femtoscopic fit (#chiEFT LO)");
  fProtonLambda->SetNDivisions(505);
  fProtonLambda->SetRangePlotting(0, 0.2, 0.8, 2);
  fProtonLambda->DrawCorrelationPlot(c_PL);
  DrawSystemInfo(c_PL);
  ref.DrawLatex(0.562, 0.465, "Nucl. Phys. A915 (2013) 24");
  c_PL->SaveAs("CF_pL_Gauss_prelim.pdf");

  TCanvas* c_LL = new TCanvas("CFLL", "CFLL", 0, 0, 650, 550);
  c_LL->SetRightMargin(right);
  c_LL->SetTopMargin(top);
  fLambdaLambda->SetLegendName("#Lambda-#Lambda #oplus #bar{#Lambda}-#bar{#Lambda} pairs");
  fLambdaLambda->SetLegendName("Femtoscopic fit ");
  fLambdaLambda->SetNDivisions(505);
  fLambdaLambda->SetRangePlotting(0, 0.2, 0.35, 2.);
  fLambdaLambda->DrawCorrelationPlot(c_LL);
  DrawSystemInfo(c_LL);
  c_LL->SaveAs("CF_LL_Gauss_prelim.pdf");

  TCanvas* c_pXi = new TCanvas("CFpXi", "CFpXi", 0, 0, 650, 550);
  c_pXi->SetRightMargin(right);
  c_pXi->SetTopMargin(top);
  fProtonXi->SetLegendName("p-#Xi^{-} #oplus #bar{p}-#bar{#Xi}^{+} pairs");
  fProtonXi->SetLegendName("Femtoscopic fit (HAL-QCD)");
  fProtonXi->SetLegendName("Femtoscopic fit (Coulomb)");
  fProtonXi->SetNDivisions(505);
  fProtonXi->SetRangePlotting(0, 0.2, 0.75, 4.);
  fProtonXi->DrawCorrelationPlot(c_pXi);
  DrawSystemInfo(c_pXi);
  ref.DrawLatex(0.562, 0.465, "Nucl. Phys. A967 (2017) 856");
  ref.DrawLatex(0.562, 0.415, "PoS LATTICE2016 (2017) 116");
  c_pXi->SaveAs("CF_pXi_Gauss_prelim.pdf");
}

void DreamPlot::DrawSystemInfo(TCanvas* c) {
  c->cd();
  TLatex BeamText;
  TLatex text;
  BeamText.SetTextSize(gStyle->GetTextSize() * 0.85);
  BeamText.SetNDC(kTRUE);
  BeamText.DrawLatex(0.5, 0.875, "ALICE Preliminary");
  TString CollisionSystem=Form("%s",fCollisionSystem);
  std::cout << CollisionSystem.Data() << std::endl;
  if (CollisionSystem.Index("Pb") > 0)
    BeamText.DrawLatex(
        0.5, 0.825,
        Form("%s #sqrt{#it{s}_{NN}} = %.2f TeV", fCollisionSystem, fEnergy));
  else
    BeamText.DrawLatex(
        0.5, 0.825,
        Form("%s #sqrt{#it{s}} = %i TeV", fCollisionSystem, (int) fEnergy));
  text.SetTextSize(gStyle->GetTextSize() * 0.75);
  text.SetNDC();
  text.SetTextColor(1);
  text.DrawLatex(
      0.5,
      0.775,
      Form("#it{r_{0}} = %.3f #pm %.3f ^{+%.3f}_{-%.3f} fm", fRadius,
           fRadiusStat, fRadiusSysUp, fRadiusSysLow));
}

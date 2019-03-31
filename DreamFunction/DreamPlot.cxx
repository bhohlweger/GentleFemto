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
      fProtonSigma(nullptr),
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
  fProtonSigma = new DreamData("ProtonSigma0");
}

DreamPlot::~DreamPlot() {
  // TODO Auto-generated destructor stub
}

void DreamPlot::ReadData(const char* PathToDataFolder,
                         const char* PathToSysFolder, int binWidth,
                         int UnitConvData) {
  fProtonProton->SetUnitConversionData(UnitConvData);
  fProtonLambda->SetUnitConversionData(UnitConvData);
  fLambdaLambda->SetUnitConversionData(UnitConvData);
  fProtonXi->SetUnitConversionData(UnitConvData);
  TString HistName = "hCk_ReweightedMeV_";
  if (binWidth == 16) {
    HistName += "0";
  } else if (binWidth == 20) {
    HistName += "1";
  } else {
    std::cout << "Unknown bin width " << binWidth << std::endl;
  }
  double sysWidth = 5;
  TFile* CFFile_pp = TFile::Open(Form("%s/CFOutput_pp.root", PathToDataFolder));
  TFile* CFFile_ppSys = TFile::Open(
      Form("%s/Systematics_pp.root", PathToSysFolder));
  fProtonProton->SetCorrelationFunction(
      (TH1F*) CFFile_pp->Get("hCk_ReweightedMeV_0"));

//  TF1 *pPbSys = new TF1("pPbSyst", [&](double *x, double *p) {
//    if (x[0] > 0 && x[0] < 0.04) {return 0.051;}
//    else {return p[0] + p[1]*x[0] + p[2]*x[0]*x[0];}},
//                        0, 3.0, 3);
//  TF1 *systematics = (TF1*) CFFile_ppSys->Get("RelSysPPUnbinned");
//  pPbSys->SetParameter(0, systematics->GetParameter(0));
//  pPbSys->SetParameter(1, systematics->GetParameter(1));
//  pPbSys->SetParameter(2, systematics->GetParameter(2));

  fProtonProton->SetSystematics((TF1*) CFFile_ppSys->Get("SystError"),
                                2);

  TFile* CFFile_pL = TFile::Open(Form("%s/CFOutput_pL.root", PathToDataFolder));
  TFile* CFFile_pLSys = TFile::Open(
      Form("%s/Systematics_pL.root", PathToSysFolder));
  fProtonLambda->SetCorrelationFunction(
      (TH1F*) CFFile_pL->Get(HistName.Data()));
  fProtonLambda->SetSystematics((TF1*) CFFile_pLSys->Get("SystError"),
                                sysWidth);

  TFile* CFFile_LL = TFile::Open(Form("%s/CFOutput_LL.root", PathToDataFolder));
  TFile* CFFile_LLSys = TFile::Open(
      Form("%s/Systematics_LL.root", PathToSysFolder));
  fLambdaLambda->SetCorrelationFunction(
      (TH1F*) CFFile_LL->Get(HistName.Data()));
  fLambdaLambda->SetSystematics((TF1*) CFFile_LLSys->Get("SystError"),
                                sysWidth);

  TFile* CFFile_pXi = TFile::Open(
      Form("%s/CFOutput_pXi.root", PathToDataFolder));
  TFile* CFFile_pXiSys = TFile::Open(
      Form("%s/Systematics_pXi.root", PathToSysFolder));
  fProtonXi->SetCorrelationFunction((TH1F*) CFFile_pXi->Get(HistName.Data()));
  fProtonXi->SetSystematics((TF1*) CFFile_pXiSys->Get("SystError"),
                            sysWidth);

  return;
}

void DreamPlot::ReadDataSigma(const char* PathToDataFolder,
                              const char* PathToSysFolder) {
  auto CFFile = TFile::Open(Form("%s/CFOutput_pSigma.root", PathToDataFolder));
  auto CFFile_Sys = TFile::Open(
      Form("%s/Systematics_pSigma0.root", PathToSysFolder));
  fProtonSigma->SetCorrelationFunction(
      (TH1F*) CFFile->Get("hCk_ReweightedMeV_3"));
  fProtonSigma->SetSystematics(
      (TF1*) CFFile_Sys->Get("SystError"),
      fProtonSigma->GetCorrelationFunction()->GetBinWidth(1) / 6.f);
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
  TFile* CFFile_pp = TFile::Open(Form("%s/CFOutput_pp.root", PathToSimFolder));
  fProtonProton->SetCorrelationFunctionSimulation(
      (TH1F*) CFFile_pp->Get("hCk_Reweighted_0"));

  TFile* CFFile_pL = TFile::Open(Form("%s/CFOutput_pL.root", PathToSimFolder));
  fProtonLambda->SetCorrelationFunctionSimulation(
      (TH1F*) CFFile_pL->Get(HistName.Data()));

  TFile* CFFile_LL = TFile::Open(Form("%s/CFOutput_LL.root", PathToSimFolder));
  fLambdaLambda->SetCorrelationFunctionSimulation(
      (TH1F*) CFFile_LL->Get(HistName.Data()));

  TFile* CFFile_pXi = TFile::Open(Form("%s/CFOutput_pXi.root", PathToSimFolder));
  fProtonXi->SetCorrelationFunctionSimulation(
      (TH1F*) CFFile_pXi->Get(HistName.Data()));
}

void DreamPlot::ReadFit(const char* fitPath, int UnitConvCATS) {
  fProtonProton->SetUnitConversionCATS(UnitConvCATS);
  fProtonLambda->SetUnitConversionCATS(UnitConvCATS);
  fLambdaLambda->SetUnitConversionCATS(UnitConvCATS);
  fProtonXi->SetUnitConversionCATS(UnitConvCATS);
  TString PathToFile = Form("%s/SYSTEMATICS_CutVarAdd_Global_Radius_Normal.root",
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
      fProtonProton->FemtoModelFitBands(grpp_default, grppLow, grppUp, 2, 1, 3,
                                        -3000);
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
                                        gr_pLNLO_up, 1, 1, 3, 3000);
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
      fProtonLambda->FemtoModelFitBands(gr_pLLO_default, grpLLowLO, grpLUpLO, 3,
                                        7, 3, 3244);
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
      fLambdaLambda->FemtoModelFitBands(grLLDefault, grLLLow, grLLUp, 5, 1, 3,
                                        3000);
    }
    TGraph* grpXiDefault = (TGraph*) systFit->Get("HALCoulombpXimGraphDefault");
    TGraph* grpXiLower = (TGraph*) systFit->Get("HALCoulombpXimGraphLowerLim");
    TGraph* grpXiUpper = (TGraph*) systFit->Get("HALCoulombpXimGraphUpperLim");
    if (!grpXiDefault) {
      std::cout << "no pXi Default file \n";
    } else if (!grpXiLower) {
      std::cout << "no pXi lower file \n";
    } else if (!grpXiUpper) {
      std::cout << "no pXi upper file \n";
    } else {
      fProtonXi->FemtoModelFitBands(grpXiDefault, grpXiLower, grpXiUpper, 10,
                                    10, 0, 3252);
    }
    TGraph* grpXiDefaultCoulomb = (TGraph*) systFit->Get(
        "CoulombpXimGraphDefault");
    TGraph* grpXiLowerCoulomb = (TGraph*) systFit->Get(
        "CoulombpXimGraphLowerLim");
    TGraph* grpXiUpperCoulomb = (TGraph*) systFit->Get(
        "CoulombpXimGraphUpperLim");
    if (!grpXiDefaultCoulomb) {
      std::cout << "no pXi Coulomb Default file \n";
    } else if (!grpXiLowerCoulomb) {
      std::cout << "no pXi Coulomb lower file \n";
    } else if (!grpXiUpperCoulomb) {
      std::cout << "no pXi Coulomb upper file \n";
    } else {
      fProtonXi->FemtoModelFitBands(grpXiDefaultCoulomb, grpXiLowerCoulomb,
                                    grpXiUpperCoulomb, 12, 1, 3, 4000);
    }
    TGraph* grpXiDefaultSideband= (TGraph*) systFit->Get(
        "CoulombpXimGraphSidebandDefault");
    TGraph* grpXiLowerSideband = (TGraph*) systFit->Get(
        "CoulombpXimGraphSidebandDown");
    TGraph* grpXiUpperSideband = (TGraph*) systFit->Get(
        "CoulombpXimGraphSidebandUp");
    if (!grpXiDefaultSideband) {
      std::cout << "no pXi Sideband Default file \n";
    } else if (!grpXiLowerSideband) {
      std::cout << "no pXi Sideband lower file \n";
    } else if (!grpXiUpperSideband) {
      std::cout << "no pXi Sideband upper file \n";
    } else {
      fProtonXi->FemtoModelFitBands(grpXiDefaultSideband, grpXiLowerSideband,
                                    grpXiUpperSideband, 4, 2, 3, 3003,true);
    }
  } else {
    std::cout << "No Cats file!  \n";
  }
  return;
}

void DreamPlot::ReadFitSigma(const char* fitPath) {
  auto ledniFile = TFile::Open(Form("%s/Param_pSigma0_2.root", fitPath));
  if (ledniFile) {
    auto ledniband = (TGraphErrors*) ledniFile->Get("CF_fit");
    if (!ledniband) {
      std::cout << "No coupled Lednicky \n";
    } else {
      fProtonSigma->FemtoModelFitBands(ledniband, kRed + 2, 0, 0, 3252, true);
    }
  } else {
    std::cout << "No Lednicky file!  \n";
  }

  auto haidenbauerFile = TFile::Open(Form("%s/Param_pSigma0_3.root", fitPath));
  if (haidenbauerFile) {
    auto haidenbauerband = (TGraphErrors*) haidenbauerFile->Get("CF_fit");
    auto sideband = (TGraphErrors*) haidenbauerFile->Get("CF_sidebands");
    if (!haidenbauerband) {
      std::cout << "No coupled Lednicky \n";
    } else if (!sideband) {
      std::cout << "No sideband \n";
    } else {
      fProtonSigma->FemtoModelFitBands(haidenbauerband, kAzure - 3, 0, 0, 3225,
                                       true);
      fProtonSigma->FemtoModelFitBands(sideband, kCyan + 2, 0.5, true);
    }
  } else {
    std::cout << "No Haidenbauer file!  \n";
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
  gStyle->SetErrorX(0.005);
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
  ref.SetTextSize(gStyle->GetTextSize() * 0.4);
  ref.SetNDC(kTRUE);
  TLatex Numbering;
  Numbering.SetTextSize(gStyle->GetTextSize() * 1.3);
  Numbering.SetNDC(kTRUE);
  TCanvas* c_PP = new TCanvas("CFpp", "CFpp", 0, 0, 650, 550);
  c_PP->SetRightMargin(right);
  c_PP->SetTopMargin(top);
  fProtonProton->SetLegendName("p-p #oplus #bar{p}-#bar{p}", "fpe");
  fProtonProton->SetLegendName("Coulomb + Argonne #nu_{18} (fit)", "l");
  fProtonProton->SetRangePlotting(0, 200, 0.7, 3.5);
  fProtonProton->SetInletRangePlotting(50,200,0.94,1.06);
  fProtonProton->SetInletCoordinates(0.35, 0.27, 0.95, 0.61);
  fProtonProton->SetNDivisions(505);
  fProtonProton->SetLegendCoordinates(
      0.37, 0.71 - 0.09 * fProtonProton->GetNumberOfModels(), 0.7, 0.8);
  fProtonProton->DrawCorrelationPlot(c_PP);
  DrawSystemInfo(c_PP, true, 0.39);
  c_PP->cd();
  Numbering.DrawLatex( 0.3,
                       0.9,"#bf{a)}");
  c_PP->SaveAs("CF_pp_Gauss_prelim.pdf");

  TCanvas* c_PL = new TCanvas("CFpL", "CFpL", 0, 0, 650, 550);
  c_PL->SetRightMargin(right);
  c_PL->SetTopMargin(top);
  fProtonLambda->SetLegendName("p-#Lambda #oplus #bar{p}-#bar{#Lambda}",
                               "fpe");
  fProtonLambda->SetLegendName("Femtoscopic fit (#chiEFT NLO)", "l");
  fProtonLambda->SetLegendName("Femtoscopic fit (#chiEFT LO)", "l");
  fProtonLambda->SetNDivisions(505);
  fProtonLambda->SetRangePlotting(0, 200, 0.8, 2);
  fProtonLambda->SetLegendCoordinates(
      0.5, 0.785 - 0.09 * fProtonLambda->GetNumberOfModels(), 0.7, 0.875);
  fProtonLambda->DrawCorrelationPlot(c_PL);
  DrawSystemInfo(c_PL, false);
  if (fProtonLambda->GetNumberOfModels() > 0) {
    ref.DrawLatex(0.5075, 0.765 - 0.09 * fProtonLambda->GetNumberOfModels(),
                  "Nucl. Phys. A915 (2013) 24");
  }
  c_PL->SaveAs("CF_pL_Gauss_prelim.pdf");

  TCanvas* c_LL = new TCanvas("CFLL", "CFLL", 0, 0, 650, 550);
  c_LL->SetRightMargin(right);
  c_LL->SetTopMargin(top);
  fLambdaLambda->SetLegendName(
      "#Lambda-#Lambda #oplus #bar{#Lambda}-#bar{#Lambda}", "fpe");
  fLambdaLambda->SetLegendName("Femtoscopic fit ", "l");
  fLambdaLambda->SetNDivisions(505);
  fLambdaLambda->SetRangePlotting(0, 200, 0.35, 2.);
  fLambdaLambda->SetLegendCoordinates(
      0.5, 0.785 - 0.09 * fLambdaLambda->GetNumberOfModels(), 0.7, 0.875);
  fLambdaLambda->DrawCorrelationPlot(c_LL);
  DrawSystemInfo(c_LL, false);
  c_LL->SaveAs("CF_LL_Gauss_prelim.pdf");

  TCanvas* c_pXi = new TCanvas("CFpXi", "CFpXi", 0, 0, 650, 550);
  c_pXi->SetRightMargin(right);
  c_pXi->SetTopMargin(top);
  fProtonXi->SetLegendName("p-#Xi^{-} #oplus #bar{p}-#bar{#Xi}^{+}",
                           "fpe");
  fProtonXi->SetLegendName("Coulomb + HAL-QCD", "fl");
  fProtonXi->SetLegendName("Coulomb", "l");
  fProtonXi->SetLegendName("p-#Xi^{-} sideband background", "l");
  fProtonXi->SetNDivisions(505);
  fProtonXi->SetRangePlotting(0, 300, 0.8, 2.5);
  fProtonXi->SetLegendCoordinates(0.35,
                                  0.785 - 0.09 * fProtonXi->GetNumberOfModels(),
                                  0.7, 0.875);
  fProtonXi->DrawCorrelationPlot(c_pXi);
  DrawSystemInfo(c_pXi, false, 0.37);
  Numbering.DrawLatex( 0.3,
                       0.9,"#bf{b)}");
  c_pXi->SaveAs("CF_pXi_Gauss_prelim.pdf");
}

void DreamPlot::DrawCorrelationFunctionSigma(const char* fitPath) {
  SetStyle();
  const float right = 0.025;
  const float top = 0.025;
  TLatex ref;
  ref.SetTextSize(gStyle->GetTextSize() * 0.4);
  ref.SetNDC(kTRUE);
  TLatex Numbering;
  Numbering.SetTextSize(gStyle->GetTextSize() * 1.3);
  Numbering.SetNDC(kTRUE);
  auto c = new TCanvas("CFpSigma", "CFpSigma", 0, 0, 650, 550);
  c->SetRightMargin(right);
  c->SetTopMargin(top);
  fProtonSigma->SetLegendName("p#minus#Sigma^{0} #oplus #bar{p}#minus#bar{#Sigma^{0}}", "fpe");
  fProtonSigma->SetLegendName("Lednicky coupled channel", "fl");
  fProtonSigma->SetLegendName("#chiEFT (NLO)", "fl");
  fProtonSigma->SetLegendName("p-#Sigma^{0} sideband background", "l");
  fProtonSigma->SetRangePlotting(0, 450, 0.9, 1.6);
  fProtonSigma->SetNDivisions(505);
  fProtonSigma->SetLegendCoordinates(
      0.45, 0.71 - 0.09 * fProtonSigma->GetNumberOfModels(), 0.7, 0.8);
  // Necessary fix to get the right unit on the axes
  fProtonSigma->SetUnitConversionData(2);
  fProtonSigma->DrawCorrelationPlot(c, 13);
  DrawSystemInfo(c, false, 0.46, true);
  c->cd();
  c->SaveAs(Form("%s/CF_pSigma_prelim.pdf", fitPath));
}

void DreamPlot::DrawSystemInfo(TCanvas* c, bool plotRadius, float xMin, bool isPreliminary) {
  c->cd();
  TLatex BeamText;
  TLatex text;
  BeamText.SetTextSize(gStyle->GetTextSize() * 0.95);
  BeamText.SetNDC(kTRUE);
//  BeamText.DrawLatex(0.5, 0.875, "ALICE");
  TString CollisionSystem = Form("%s", fCollisionSystem);
  if (CollisionSystem.Index("Pb") > 0) {
    if (isPreliminary && !plotRadius) {
      BeamText.DrawLatex(xMin, 0.9, "ALICE Preliminary");
      BeamText.DrawLatex(
          xMin, 0.825,
          Form("%s #sqrt{#it{s}_{NN}} = %.2f TeV", fCollisionSystem, fEnergy));
    } else {
      BeamText.DrawLatex(
          xMin,
          0.9,
          Form("ALICE %s #sqrt{#it{s}_{NN}} = %.2f TeV", fCollisionSystem,
               fEnergy));
    }
  } else {
    if (isPreliminary && !plotRadius) {
      BeamText.DrawLatex(xMin, 0.9, "ALICE Preliminary");
      BeamText.DrawLatex(
          xMin, 0.825,
          Form("%s #sqrt{#it{s}} = %i TeV", fCollisionSystem, (int) fEnergy));
    } else {
      BeamText.DrawLatex(
          xMin,
          0.9,
          Form("ALICE %s #sqrt{#it{s}} = %i TeV", fCollisionSystem,
               (int) fEnergy));
    }
  }
  if (plotRadius) {
    text.SetNDC();
    text.SetTextColor(1);
    text.SetTextSize(gStyle->GetTextSize() * 0.83);
    text.DrawLatex(
        xMin,
        0.825,
        Form("#it{r}_{0} = %.3f #pm %.3f (stat.) ^{+%.3f}_{-%.3f} (syst.) fm", fRadius,
             fRadiusStat, fRadiusSysUp, fRadiusSysLow));
  }
}

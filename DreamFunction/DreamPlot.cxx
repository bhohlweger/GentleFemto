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
      fProtonSigmaSideband(nullptr),
      fProtonAntiProton(nullptr),
      fProtonAntiLambda(nullptr),
      fLambdaAntiLambda(nullptr),
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
  fProtonSigmaSideband = new DreamData("ProtonSigma0Sidebands");
  fProtonAntiProton = new DreamData("ProtonAntiProton");
  fProtonAntiLambda = new DreamData("ProtonAntiLambda");
  fLambdaAntiLambda = new DreamData("LambdaAntiLambda");
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
  if (!CFFile_ppSys) {
    CFFile_ppSys = TFile::Open(Form("%s/C2totalsysPP.root", PathToSysFolder));
  }
  fProtonProton->SetCorrelationFunction(
      (TH1F*) CFFile_pp->Get("hCk_ReweightedMeV_0"));
  TF1* sysParam = (TF1*) CFFile_ppSys->Get("SystError");
  if (sysParam) {
    fProtonProton->SetSystematics(sysParam, 2);
  } else {
    fProtonProton->SetSystematics((TF1*) CFFile_ppSys->Get("RelSysPPUnbinned"), 2);
  }
  TFile* CFFile_pL = TFile::Open(Form("%s/CFOutput_pL.root", PathToDataFolder));
  TFile* CFFile_pLSys = TFile::Open(
      Form("%s/Systematics_pL.root", PathToSysFolder));
  if (!CFFile_pLSys) {
    CFFile_pLSys = TFile::Open(Form("%s/C2totalsysPL.root", PathToSysFolder));
  }
  fProtonLambda->SetCorrelationFunction(
      (TH1F*) CFFile_pL->Get(HistName.Data()));
  if (CFFile_pLSys->Get("SystError")) {
    fProtonLambda->SetSystematics((TF1*) CFFile_pLSys->Get("SystError"),
                                  sysWidth);
  } else {
    fProtonLambda->SetSystematics((TF1*) CFFile_pLSys->Get("RelSysPLUnbinned"),
                                  sysWidth);
  }
  TFile* CFFile_LL = TFile::Open(Form("%s/CFOutput_LL.root", PathToDataFolder));
  TFile* CFFile_LLSys = TFile::Open(
      Form("%s/Systematics_LL.root", PathToSysFolder));
  if (!CFFile_LLSys) {
    CFFile_LLSys = TFile::Open(Form("%s/C2totalsysLL.root", PathToSysFolder));
  }
  fLambdaLambda->SetCorrelationFunction(
      (TH1F*) CFFile_LL->Get(HistName.Data()));
  if (CFFile_LLSys->Get("SystError")) {
    fLambdaLambda->SetSystematics((TF1*) CFFile_LLSys->Get("SystError"),
                                  sysWidth);
  } else {
    fLambdaLambda->SetSystematics((TF1*) CFFile_LLSys->Get("RelSysLLUnbinned"),
                                  sysWidth);
  }

  TFile* CFFile_pXi = TFile::Open(
      Form("%s/CFOutput_pXi.root", PathToDataFolder));
  TFile* CFFile_pXiSys = TFile::Open(
      Form("%s/Systematics_pXi.root", PathToSysFolder));
  if (!CFFile_pXiSys) {
    CFFile_pXiSys = TFile::Open(Form("%s/C2totalsysPXi.root", PathToSysFolder));
  }
  fProtonXi->SetCorrelationFunction((TH1F*) CFFile_pXi->Get(HistName.Data()));
  if (CFFile_pXiSys->Get("SystError")) {
    fProtonXi->SetSystematics((TF1*) CFFile_pXiSys->Get("SystError"), sysWidth);
  } else {
    fProtonXi->SetSystematics((TF1*) CFFile_pXiSys->Get("RelSysPXiUnbinned"),
                              sysWidth);
  }
  return;
}

void DreamPlot::ReadDataSigma(const char* PathToDataFolder,
                              const char* PathToSysFolder) {
  auto CFFile = TFile::Open(Form("%s/CFOutput_pSigma.root", PathToDataFolder));
  fProtonSigma->SetCorrelationGraph(
      (TGraphAsymmErrors*) CFFile->Get("Graph_from_hCk_Reweighted_3MeV"));
  auto CFFile_Sys = TFile::Open(
      Form("%s/Systematics_pSigma0.root", PathToSysFolder));
  fProtonSigma->SetSystematics(
      (TF1*) CFFile_Sys->Get("SystError"), 6.5);
}

void DreamPlot::ReadSidebandSigma(const char* PathToFitFolder,
                                  const char* PathToSysFolder) {
  auto fitFile = TFile::Open(Form("%s/Param_pSigma0_2.root", PathToFitFolder));
  auto histSideband = (TGraphAsymmErrors*) fitFile->Get("SidebandMerged_0");
  histSideband->SetName(Form("%sMeV", histSideband->GetName()));
  fProtonSigmaSideband->SetCorrelationGraph(histSideband);

  auto CFFile_Sys = TFile::Open(
      Form("%s/SystematicsSidebands_pSigma0.root", PathToSysFolder));
  fProtonSigmaSideband->SetSystematics((TF1*) CFFile_Sys->Get("SystError"), 6.5);
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

  TFile* CFFile_pXi = TFile::Open(
      Form("%s/CFOutput_pXi.root", PathToSimFolder));
  fProtonXi->SetCorrelationFunctionSimulation(
      (TH1F*) CFFile_pXi->Get(HistName.Data()));
}

void DreamPlot::ReadFit(const char* fitPath, int UnitConvCATS) {
  fProtonProton->SetUnitConversionCATS(UnitConvCATS);
  fProtonLambda->SetUnitConversionCATS(UnitConvCATS);
  fLambdaLambda->SetUnitConversionCATS(UnitConvCATS);
  fProtonXi->SetUnitConversionCATS(UnitConvCATS);
  TString PathToFile = Form(
      "%s/SYSTEMATICS_CutVarAdd_Global_Radius_Normal.root", fitPath);
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
//    TGraph* grpXiHALDefault = (TGraph*) systFit->Get("pXimGraphDefault");
//    TGraph* grpXiHALLower = (TGraph*) systFit->Get("pXimGraphLowerLim");
//    TGraph* grpXiHALUpper = (TGraph*) systFit->Get("pXimGraphUpperLim");
    TGraph* grpXiHALDefault = (TGraph*) systFit->Get("HALQCDpXimGraphDefault");
    TGraph* grpXiHALLower = (TGraph*) systFit->Get("HALQCDpXimGraphLowerLim");
    TGraph* grpXiHALUpper = (TGraph*) systFit->Get("HALQCDpXimGraphUpperLim");
    if (!grpXiHALDefault) {
      std::cout << "no pXi Default file \n";
    } else if (!grpXiHALLower) {
      std::cout << "no pXi lower file \n";
    } else if (!grpXiHALUpper) {
      std::cout << "no pXi upper file \n";
    } else {
      fProtonXi->FemtoModelFitBands(grpXiHALDefault, grpXiHALLower,
                                    grpXiHALUpper, 10, 10, 0, 3252);
    }
    TGraph* grpXiTomDefault = (TGraph*) systFit->Get("RikkenpXimGraphDefault");
    TGraph* grpXiTomLower = (TGraph*) systFit->Get("RikkenpXimGraphLowerLim");
    TGraph* grpXiTomUpper = (TGraph*) systFit->Get("RikkenpXimGraphUpperLim");
    if (!grpXiTomDefault) {
      std::cout << "no pXi Default file \n";
    } else if (!grpXiTomLower) {
      std::cout << "no pXi lower file \n";
    } else if (!grpXiTomUpper) {
      std::cout << "no pXi upper file \n";
    } else {
      fProtonXi->FemtoModelFitBands(grpXiTomDefault, grpXiTomLower,
                                    grpXiTomUpper, 11, 8, 0, 3225);
    }
//    TGraph* grpXiDefaultCoulomb = (TGraph*) systFit->Get(
//        "pXimGraphDefault_COULOMB");
//    TGraph* grpXiLowerCoulomb = (TGraph*) systFit->Get(
//        "pXimGraphLowerLim_COULOMB");
//    TGraph* grpXiUpperCoulomb = (TGraph*) systFit->Get(
//        "pXimGraphUpperLim_COULOMB");
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
//    TGraph* grpXiDefaultSideband = (TGraph*) systFit->Get("pXimGraphSidebandDefault");
//    TGraph* grpXiLowerSideband = (TGraph*) systFit->Get("pXimGraphSidebandDown");
//    TGraph* grpXiUpperSideband = (TGraph*) systFit->Get("pXimGraphSidebandUp");
    TGraph* grpXiDefaultSideband = (TGraph*) systFit->Get("pXiSidebandDefault");
    TGraph* grpXiLowerSideband = (TGraph*) systFit->Get("pXiSidebandDown");
    TGraph* grpXiUpperSideband = (TGraph*) systFit->Get("pXiSidebandUp");
    if (!grpXiDefaultSideband) {
      std::cout << "no pXi Sideband Default file \n";
    } else if (!grpXiLowerSideband) {
      std::cout << "no pXi Sideband lower file \n";
    } else if (!grpXiUpperSideband) {
      std::cout << "no pXi Sideband upper file \n";
    } else {
      fProtonXi->FemtoModelFitBands(grpXiDefaultSideband, grpXiLowerSideband,
                                    grpXiUpperSideband, 4, 2, 3, 3003, true);
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
    auto lednibandFlat = (TGraphErrors*) ledniFile->Get("CF_fit_deviation");
    if (!ledniband) {
      std::cout << "No coupled Lednicky \n";
    } else {
      fProtonSigma->FemtoModelFitBands(ledniband, kRed + 2 , 0.5, true);
      fProtonSigma->FemtoModelDeviations(lednibandFlat, kRed + 1, 0, 0, 3385,
                                         false);
    }
  } else {
    std::cout << "No Lednicky file!  \n";
  }

  auto haidenbauerFile = TFile::Open(Form("%s/Param_pSigma0_3.root", fitPath));
  if (haidenbauerFile) {
    auto haidenbauerband = (TGraphErrors*) haidenbauerFile->Get("CF_fit");
    auto haidenbauerbandFlat = (TGraphErrors*) ledniFile->Get("CF_fit_deviation");
    if (!haidenbauerband) {
      std::cout << "No coupled Lednicky \n";
    } else {
      fProtonSigma->FemtoModelFitBands(haidenbauerband, kAzure, 0.5, true);
      fProtonSigma->FemtoModelDeviations(haidenbauerbandFlat, kAzure, 0, 0,
                                         3325, false);
    }
  } else {
    std::cout << "No Haidenbauer file!  \n";
  }

  auto ESC16File = TFile::Open(Form("%s/Param_pSigma0_4.root", fitPath));
  if (ESC16File) {
    auto esc16band = (TGraphErrors*) ESC16File->Get("CF_fit");
    auto esc16bandFlat = (TGraphErrors*) ESC16File->Get("CF_fit_deviation");
    if (!esc16band) {
      std::cout << "No ESC16 \n";
    } else {
      fProtonSigma->FemtoModelFitBands(esc16band, kGreen + 2, 0.6, true);
      fProtonSigma->FemtoModelDeviations(esc16bandFlat, kGreen + 2, 0, 0, 3352, false);
    }
  } else {
    std::cout << "No ESC16 file!  \n";
  }

  auto NSC97fFile = TFile::Open(Form("%s/Param_pSigma0_5.root", fitPath));
  if (NSC97fFile) {
    auto NSC97fband = (TGraphErrors*) NSC97fFile->Get("CF_fit");
    auto NSC97fbandFlat = (TGraphErrors*) NSC97fFile->Get("CF_fit_deviation");
    auto sideband = (TGraphErrors*) NSC97fFile->Get("CF_sidebands");
    auto sidebandUnscaled = (TGraphErrors*) NSC97fFile->Get("CF_genuineSidebands");
    if (!NSC97fband) {
      std::cout << "No NSC97f \n";
    } else if (!sideband) {
      std::cout << "No sideband \n";
    } else if (!sidebandUnscaled) {
      std::cout << "No sideband (unscaled) \n";
    } else {
      fProtonSigma->FemtoModelFitBands(NSC97fband, kOrange - 3, 0.5, true);
      fProtonSigma->FemtoModelDeviations(NSC97fbandFlat, kOrange - 3, 0, 0, 3358, false);
      fProtonSigmaSideband->FemtoModelFitBands(sidebandUnscaled, kGray + 1, 0.5, true);
    }
  } else {
    std::cout << "No NSC97f file! \n";
  }

  auto sidebandFile = TFile::Open(Form("%s/Param_pSigma0_6.root", fitPath));
  if (sidebandFile) {
    auto sideband = (TGraphErrors*) sidebandFile->Get("CF_fit");
    auto sidebandFlat = (TGraphErrors*) sidebandFile->Get("CF_fit_deviation");
    auto correlatedErrorSB = (TGraphErrors*) sidebandFile->Get("CF_correlatedError");
    if (!sideband) {
      std::cout << "No Sideband \n";
    } else {
      fProtonSigma->FemtoModelFitBands(sideband, kGray + 1, 0.9, true);
      fProtonSigma->FemtoModelDeviations(sidebandFlat, kGray, 0, 0, 0, false);
      fProtonSigma->SetCorrelatedError(correlatedErrorSB, kGray + 1, 3352, false);
    }
  } else {
    std::cout << "No Sideband file! \n";
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
  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadTopMargin(0.025);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(2);
  gStyle->SetLabelFont(43, "xyz");
  gStyle->SetTitleFont(43, "xyz");
  gStyle->SetLabelSize(28, "xyz");
  gStyle->SetTitleSize(28, "xyz");
  gStyle->SetLabelOffset(0.01, "xy");
  gStyle->SetLabelColor(kBlack, "xyz");
  gStyle->SetTitleOffset(1.25, "y");
  gStyle->SetTitleOffset(1.25, "x");
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetLegendFont(42);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetErrorX(0.005);
  gStyle->SetHatchesLineWidth(1);
  gStyle->SetHatchesSpacing(0.5);
  gStyle->SetPalette(kCividis);
}

void DreamPlot::SetStyleHisto(TH1 *histo, int marker, int color, float alpha) {
  histo->GetXaxis()->SetLabelSize(28);
  histo->GetXaxis()->SetTitleSize(28);
  histo->GetXaxis()->SetLabelOffset(0.01);
  histo->GetXaxis()->SetTitleOffset(1.2);
  histo->GetXaxis()->SetLabelFont(43);
  histo->GetXaxis()->SetTitleFont(43);
  histo->GetYaxis()->SetLabelSize(28);
  histo->GetYaxis()->SetTitleSize(28);
  histo->GetYaxis()->SetLabelOffset(0.01);
  histo->GetYaxis()->SetTitleOffset(1.25);
  histo->GetYaxis()->SetLabelFont(43);
  histo->GetYaxis()->SetTitleFont(43);
  histo->SetMarkerStyle(marker);
  if (alpha > 0) {
    histo->SetMarkerColorAlpha(color, alpha);
    histo->SetLineColorAlpha(color, alpha);
  } else {
    histo->SetMarkerColor(color);
    histo->SetLineColor(color);
  }
}

void DreamPlot::SetStyleHistoCF(TH1 *histo, int marker, int color, int labelsize) {
 histo->GetXaxis()->SetLabelSize(labelsize);
 histo->GetXaxis()->SetTitleSize(0.05);
 histo->GetXaxis()->SetLabelOffset(0.01);
 histo->GetXaxis()->SetTitleOffset(1.2);
 histo->GetXaxis()->SetLabelFont(43);
 histo->GetXaxis()->SetTitle("k* [MeV/c]");
  // histo->GetYaxis()->SetLabelSize(15);
 histo->GetYaxis()->SetTitleSize(0.05);
 histo->GetYaxis()->SetLabelOffset(0.01);
 histo->GetYaxis()->SetTitleOffset(1.25);
 histo->GetYaxis()->SetTitle("C(k^{*})");
 histo->SetMarkerSize(0.6);
 histo->SetMarkerStyle(marker);
 histo->SetMarkerColor(color);
 histo->SetLineColor(color);
}

void DreamPlot::SetStyleGraph(TGraph *histo, int marker, int color, float alpha) {
  histo->GetXaxis()->SetLabelSize(28);
  histo->GetXaxis()->SetTitleSize(28);
  histo->GetXaxis()->SetLabelOffset(0.01);
  histo->GetXaxis()->SetTitleOffset(1.2);
  histo->GetXaxis()->SetLabelFont(43);
  histo->GetYaxis()->SetLabelSize(28);
  histo->GetYaxis()->SetTitleSize(28);
  histo->GetYaxis()->SetLabelOffset(0.01);
  histo->GetYaxis()->SetTitleOffset(1.25);
  histo->SetMarkerStyle(marker);
  histo->SetMarkerColorAlpha(color, alpha);
  histo->SetLineColorAlpha(color, alpha);
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
  fProtonProton->SetRangePlotting(0, 200, 0.8, 3.5);
  fProtonProton->SetInletRangePlotting(50, 375, 0.94, 1.06);
  fProtonProton->SetInletCoordinates(0.27, 0.22, 0.95, 0.61);
  fProtonProton->SetNDivisions(505);
  fProtonProton->SetLegendCoordinates(
      0.30, 0.71 - 0.09 * fProtonProton->GetNumberOfModels(), 0.7, 0.8);
  fProtonProton->DrawCorrelationPlot(c_PP);
  DrawSystemInfo(c_PP, true, 0.32);
  c_PP->cd();
//  Numbering.DrawLatex(0.19, 0.9, "#bf{a)}");
  c_PP->SaveAs("CF_pp_Gauss_prelim.pdf");

  TCanvas* c_PL = new TCanvas("CFpL", "CFpL", 0, 0, 650, 550);
  c_PL->SetRightMargin(right);
  c_PL->SetTopMargin(top);
  fProtonLambda->SetLegendName("p-#Lambda #oplus #bar{p}-#bar{#Lambda}", "fpe");
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
  fProtonXi->SetLegendName("p-#Xi^{-} #oplus #bar{p}-#bar{#Xi}^{+}", "fpe");
  fProtonXi->SetLegendName("Coulomb + HAL-QCD", "fl");
  fProtonXi->SetLegendName("Coulomb + ESC 16", "fl");
  fProtonXi->SetLegendName("Coulomb", "l");
  fProtonXi->SetLegendName("p-#Xi^{-} sideband background", "l");
  fProtonXi->SetNDivisions(505);
  fProtonXi->SetRangePlotting(0, 300, 0.8, 2.6);
  fProtonXi->SetLegendCoordinates(0.4,
                                  0.75 - 0.09 * fProtonXi->GetNumberOfModels(),
                                  0.7, 0.75);
  fProtonXi->DrawCorrelationPlot(c_pXi);
  DrawSystemInfo(c_pXi, false, 0.42, 2);
//  Numbering.DrawLatex(0.19, 0.9, "#bf{b)}");
  c_pXi->SaveAs("CF_pXi_Gauss_prelim.pdf");
}

void DreamPlot::DrawCorrelationFunctionSigma(const char* fitPath) {
  SetStyle();
  const float right = 0.01;
  const float top = 0.02;
  auto c = new TCanvas("CFpSigma", "CFpSigma", 0, 0, 650, 1000);
  TPad *p1 = new TPad("p1", "p1", 0., 0.45, 1., 1.);
  p1->SetBottomMargin(0.0);
  p1->SetRightMargin(right);
  p1->SetTopMargin(top);
  c->SetRightMargin(right);
  c->SetTopMargin(top);
  fProtonSigma->SetLegendName(
      "p#minus#kern[-0.95]{ }#Sigma^{0} #oplus #bar{p}#minus#kern[-0.85]{ }#bar{#Sigma^{0}}", "fpe");
  fProtonSigma->SetLegendName("fss2", "l");
  fProtonSigma->SetLegendName("#chiEFT (NLO)", "l");
  fProtonSigma->SetLegendName("ESC16", "l");
  fProtonSigma->SetLegendName("NSC97f", "l");
  fProtonSigma->SetLegendName("p#minus#kern[-1.]{ }(#Lambda#gamma) baseline", "l");
  fProtonSigma->SetRangePlotting(0, 365, 0.85, 1.7);
  fProtonSigma->SetNDivisions(504);
  fProtonSigma->SetForceAxisRanges(true);
  const float leftX = 0.475;
  const float upperY = 0.82;
  fProtonSigma->SetLegendCoordinates(
      leftX, upperY - 0.07 * (fProtonSigma->GetNumberOfModels() + 1), 0.7, upperY);
  // Necessary fix to get the right unit on the axes
  fProtonSigma->SetUnitConversionData(2);
  fProtonSigma->DrawCorrelationPlot(p1, 13, kBlue + 3, 0.9);

  DrawSystemInfo(p1, false, leftX + 0.01, 0);
  c->cd();
  p1->Draw();
  fProtonSigma->DrawDeviationPerBin(c, 0, .45, 3);
  c->SaveAs(Form("%s/CF_pSigma_deviation.pdf", fitPath));
  c->SaveAs(Form("%s/CF_pSigma_deviation.root", fitPath));

  auto d = new TCanvas("CFpSigmaSideband", "CFpSigmaSideband", 0, 0, 650, 550);
  d->SetRightMargin(right);
  d->SetTopMargin(top);
  fProtonSigmaSideband->SetLegendName(
      "p#minus#kern[-0.65]{ }(#Lambda#gamma) #oplus #bar{p}#minus#kern[-0.4]{ }(#bar{#Lambda}#gamma)", "fpe");
  fProtonSigmaSideband->SetLegendName("Parametrization", "l");
  fProtonSigmaSideband->SetRangePlotting(0, 365, 0.9, 1.7);
  fProtonSigmaSideband->SetNDivisions(504);
  fProtonSigmaSideband->SetForceAxisRanges(true);
  fProtonSigmaSideband->SetLegendCoordinates(
      leftX, upperY - 0.075 * (fProtonSigmaSideband->GetNumberOfModels() + 1), 0.7, upperY);
  // Necessary fix to get the right unit on the axes
  fProtonSigmaSideband->SetUnitConversionData(2);
  fProtonSigmaSideband->DrawCorrelationPlot(d, 13, kBlue + 3, 0.9);
  DrawSystemInfo(d, false, leftX + 0.01, 0);
  d->cd();
  d->SaveAs(Form("%s/CF_pSigma_sideband.pdf", fitPath));
  d->SaveAs(Form("%s/CF_pSigma_sideband.root", fitPath));

  auto e = new TCanvas("CFpSigmaSingle", "CFpSigmaSingle", 0, 0, 650, 550);
  e->SetRightMargin(right);
  e->SetTopMargin(top);
  fProtonSigma->SetRangePlotting(0, 365, 0.8, 1.8);
  fProtonSigma->DrawCorrelationPlot(e, 13, kBlue + 3, 0.9);
  DrawSystemInfo(e, false, leftX + 0.01, 0);
  e->SaveAs(Form("%s/CF_pSigma.pdf", fitPath));
  e->SaveAs(Form("%s/CF_pSigma.root", fitPath));
}


void DreamPlot::DrawCorrelationFunctionProtonProton(const char* path) {
  SetStyle();
  gStyle->SetHatchesSpacing(0.5);
  gStyle->SetHatchesLineWidth(1);

  TFile* CFFile_pp = TFile::Open(Form("%s/CFOutput_pp.root", path));
  fProtonProton->SetCorrelationFunction(
      (TH1F*) CFFile_pp->Get("hCk_ReweightedMeV_1"));
  TFile* CFFile_ppSys = TFile::Open("~/cernbox/SystematicsAndCalib/ppRun2_HM/Systematics_pp.root");
  TF1* sysParam = (TF1*) CFFile_ppSys->Get("SystError");
  fProtonProton->SetSystematics(sysParam, 2);

  auto fitFile = TFile::Open(Form("%s/tmp_0.root", path));
  if (fitFile) {
    auto ppFit = (TGraphErrors*) fitFile->Get("Model");
      fProtonProton->FemtoModelFitBands(ppFit, kTeal + 2, 0.8, true);
  }

  const float leftX = 0.45;
  const float upperY = 0.82;
  const float right = 0.01;
  const float top = 0.02;
  auto c = new TCanvas("CFpp", "CFpp", 0, 0, 650, 550);
  c->SetRightMargin(right);
  c->SetTopMargin(top);
  fProtonProton->SetLegendName(
      "p#minus#kern[-0.95]{ }p #oplus #bar{p}#minus#kern[-0.85]{ }#bar{p}", "fpe");
  fProtonProton->SetLegendName("Coulomb + Argonne #nu_{18} (fit)", "l");
  fProtonProton->SetRangePlotting(0, 225, 0.8, 3.5);
  fProtonProton->SetInletRangePlotting(55, 345, 0.94, 1.06);
  fProtonProton->SetInletCoordinates(0.26, 0.25, 0.975, 0.665);
  fProtonProton->SetTextSizeLegend(20);
  fProtonProton->SetAxisOffsetInlet(2.7, 1.1);

  fProtonProton->SetNDivisions(505);
  fProtonProton->SetForceAxisRanges(true);
  fProtonProton->SetLegendCoordinates(
      leftX, upperY - 0.07 * (fProtonProton->GetNumberOfModels() + 1), 0.7, upperY);
  // Necessary fix to get the right unit on the axes
  fProtonProton->SetUnitConversionData(2);
  fProtonProton->DrawCorrelationPlot(c, 13, kBlue + 3, 0.9, 1.);
  DrawSystemInfo(c, false, leftX + 0.01, 0);
  c->SaveAs(Form("%s/CF_pp.pdf", path));
  c->SaveAs(Form("%s/CF_pp.root", path));
}

void DreamPlot::DrawCorrelationFunctionsBBar(int pAp_model) {
  SetStyle();
  const float right = 0.025;
  const float top = 0.025;
  TLatex ref;
  ref.SetTextSize(gStyle->GetTextSize() * 0.4);
  ref.SetNDC(kTRUE);
  TLatex Numbering;
  Numbering.SetTextSize(gStyle->GetTextSize() * 1.3);
  Numbering.SetNDC(kTRUE);
  TCanvas* c_PAP = new TCanvas("CFpAp", "CFpAp", 0, 0, 650, 550);
  c_PAP->SetRightMargin(right);
  c_PAP->SetTopMargin(top);
  fProtonAntiProton->SetLegendName("p-#bar{p}", "fpe");
  if(pAp_model=0)fProtonAntiProton->SetLegendName("#chiEFT + Coulomb", "l");
  if(pAp_model=1)fProtonAntiProton->SetLegendName("Lednicky-Lyuboshits + Coulomb", "l");
  if(pAp_model=2)fProtonAntiProton->SetLegendName("Coulomb", "l");
  fProtonAntiProton->SetRangePlotting(0, 400, 0.8, 3.5);
  fProtonAntiProton->SetNDivisions(505);
  fProtonAntiProton->SetLegendCoordinates(
      0.30, 0.71 - 0.09 * fProtonAntiProton->GetNumberOfModels(), 0.7, 0.8);
  fProtonAntiProton->DrawCorrelationPlot(c_PAP);
  DrawSystemInfo(c_PAP, true, 0.32);
  c_PAP->cd();
  c_PAP->SaveAs("CF_pAp_%s_prelim.pdf");
}

void DreamPlot::DrawSystemInfo(TPad* c, bool plotRadius, float xMin,
                               int isPreliminary) {
  c->cd();
  TLatex BeamText;
  TLatex text;
  BeamText.SetTextSize(gStyle->GetTextSize() * .96);
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
    if (isPreliminary == 1 && !plotRadius) {
//      BeamText.DrawLatex(xMin, 0.9, "ALICE Preliminary");
      BeamText.DrawLatex(
          xMin, 0.9,
          Form("ALICE Preliminary %s #sqrt{#it{s}} = %i TeV", fCollisionSystem, (int) fEnergy));
    } else if (isPreliminary == 2 && !plotRadius) {
      BeamText.SetTextSize(gStyle->GetTextSize() * .9);
      BeamText.DrawLatex(xMin, 0.9, "#bf{ALICE Preliminary}");
      BeamText.DrawLatex(
          xMin,
          0.84,
          Form("%s #sqrt{#it{s}} = %i TeV", fCollisionSystem, (int) fEnergy));
      BeamText.DrawLatex(
          xMin,
          0.78, "High Mult. (0-0.072% INEL)");
    } else if (isPreliminary ==0 && !plotRadius) {
      BeamText.SetTextSize(gStyle->GetTextSize() * .9);
            BeamText.DrawLatex(xMin, 0.9, Form("ALICE %s #sqrt{#it{s}} = %i TeV", fCollisionSystem, (int) fEnergy));      BeamText.DrawLatex(
                xMin,
                0.84, "High-mult. (0#kern[-0.65]{ }#minus#kern[-0.65]{ }0.072#kern[-0.9]{ }% INEL#kern[-0.5]{ }>#kern[-0.5]{ }0)");
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
    text.SetTextSize(gStyle->GetTextSize() * 0.96);
    text.DrawLatex(
        xMin,
        0.825,
        Form("#it{r}_{0} = %.3f #pm %.3f (stat.) ^{+%.3f}_{-%.3f} (syst.) fm",
             fRadius, fRadiusStat, fRadiusSysUp, fRadiusSysLow));
  }
}

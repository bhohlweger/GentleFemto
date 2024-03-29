#include "TSystem.h"
#include "DLM_Source.h"
#include "DLM_Potentials.h"
#include "DLM_CkModels.h"
#include "DLM_WfModel.h"
#include "CATS.h"
#include "DLM_CkDecomposition.h"
#include "DLM_Fitters.h"
#include "TidyCats.h"
#include "CATSInputSigma0.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TF1.h"
#include "TCanvas.h"
#include <iostream>
#include "CATSLambdaParam.h"
#include "SidebandSigma.h"
#include "TMinuit.h"
#include "TMath.h"
#include "TPaveText.h"
#include "DreamPlot.h"
#include "TNtuple.h"
#include "DreamSystematics.h"
#include "TVirtualFitter.h"
#include "TSpline.h"

#include <fstream>
#include <iostream>
#include <string>

#ifndef GENTLEKITTY_SCRIPTS_PDMESON_DMESONTOOLS_H_
#define GENTLEKITTY_SCRIPTS_PDMESON_DMESONTOOLS_H_

/// =====================================================================================
TGraphAsymmErrors GetCorrelationGraph(TString filename, TString appendix,
                                      TString suffix, TString graphName,
                                      const double normLower,
                                      const double normUpper, const int rebin,
                                      MomentumGami* folderFolder = nullptr) {
  TGraphAsymmErrors outputGraph;
  ReadDreamFile *DreamFile = new ReadDreamFile(4, 4);
  DreamFile->SetAnalysisFile(filename.Data(), appendix.Data(), suffix.Data());

  DreamCF *CF_pDminus = new DreamCF();
  DreamPair *pDminus = new DreamPair("Part", normLower, normUpper);
  DreamPair *apDplus = new DreamPair("AntiPart", normLower, normUpper);

  pDminus->SetPair(DreamFile->GetPairDistributions(0, 3, ""));
  apDplus->SetPair(DreamFile->GetPairDistributions(1, 2, ""));
  pDminus->ShiftForEmpty(pDminus->GetPair());
  apDplus->ShiftForEmpty(apDplus->GetPair());
  pDminus->FixShift(pDminus->GetPairShiftedEmpty(0),
                    apDplus->GetPairShiftedEmpty(0), apDplus->GetFirstBin());
  apDplus->FixShift(apDplus->GetPairShiftedEmpty(0),
                    pDminus->GetPairShiftedEmpty(0), pDminus->GetFirstBin());
  pDminus->Rebin(pDminus->GetPairFixShifted(0), rebin, true);
  apDplus->Rebin(apDplus->GetPairFixShifted(0), rebin, true);
  pDminus->ReweightMixedEvent(pDminus->GetPairRebinned(0), 0.2, 0.9,
                              pDminus->GetPair());
  apDplus->ReweightMixedEvent(apDplus->GetPairRebinned(0), 0.2, 0.9,
                              apDplus->GetPair());
  if (folderFolder) {
    auto meUnscaled = pDminus->GetPair()->GetMEDist();
    auto meScaled = pDminus->GetPairReweighted(0)->GetMEDist();
    double ScalingFactor_pDminus = meUnscaled->Integral()
        / (double) meScaled->Integral();

    meUnscaled = apDplus->GetPair()->GetMEDist();
    meScaled = apDplus->GetPairReweighted(0)->GetMEDist();
    double ScalingFactor_ApDplus = meUnscaled->Integral()
        / (double) meScaled->Integral();

    pDminus->UnfoldMomentum(pDminus->GetPair(), folderFolder,
                            ScalingFactor_pDminus);
    apDplus->UnfoldMomentum(apDplus->GetPair(), folderFolder,
                            ScalingFactor_ApDplus);
    pDminus->ShiftForEmpty(pDminus->GetPairUnfolded(0));
    apDplus->ShiftForEmpty(apDplus->GetPairUnfolded(0));
    pDminus->FixShift(pDminus->GetPairShiftedEmpty(1),
                      apDplus->GetPairShiftedEmpty(1), apDplus->GetFirstBin());
    apDplus->FixShift(apDplus->GetPairShiftedEmpty(1),
                      pDminus->GetPairShiftedEmpty(1), pDminus->GetFirstBin());
    pDminus->Rebin(pDminus->GetPairFixShifted(1), rebin, true);
    apDplus->Rebin(apDplus->GetPairFixShifted(1), rebin, true);
//    pDminus->ReweightMixedEvent(pDminus->GetPairRebinned(1), 0.2, 0.9,
//                                pDminus->GetPair());
//    apDplus->ReweightMixedEvent(apDplus->GetPairRebinned(1), 0.2, 0.9,
//                                apDplus->GetPair());
  }

  CF_pDminus->SetPairs(pDminus, apDplus);
  CF_pDminus->GetCorrelations();

  for (auto it : CF_pDminus->GetCorrelationFunctionGraphs()) {
    TString itName = it->GetName();
    std::cout << it->GetName() << std::endl;
    if (graphName == itName) {
      std::cout << it->GetName() << std::endl;
      outputGraph = *it;
    }
  }

  delete apDplus;
  delete pDminus;
  //delete CF_pDminus;
  delete DreamFile;
  return outputGraph;
}

/// =====================================================================================
TH1F GetMEDist(TString filename, TString appendix, TString suffix,
               TString graphName, const double normLower,
               const double normUpper, const int rebin) {
  ReadDreamFile *DreamFile = new ReadDreamFile(4, 4);
  DreamFile->SetAnalysisFile(filename.Data(), appendix.Data(), suffix.Data());

  DreamPair *pDminus = new DreamPair("Part", normLower, normUpper);
  DreamPair *apDplus = new DreamPair("AntiPart", normLower, normUpper);

  pDminus->SetPair(DreamFile->GetPairDistributions(0, 3, ""));
  apDplus->SetPair(DreamFile->GetPairDistributions(1, 2, ""));
  pDminus->Rebin(pDminus->GetPair(), rebin, true);
  apDplus->Rebin(apDplus->GetPair(), rebin, true);
  pDminus->ReweightMixedEvent(pDminus->GetPairRebinned(0), 0.2, 0.9,
                              pDminus->GetPair());
  apDplus->ReweightMixedEvent(apDplus->GetPairRebinned(0), 0.2, 0.9,
                              apDplus->GetPair());
  auto meDist = *(pDminus->GetPairReweighted(0)->GetMEDist());
  meDist.Add(apDplus->GetPairReweighted(0)->GetMEDist());

  delete apDplus;
  delete pDminus;
  delete DreamFile;
  return meDist;
}

/// =====================================================================================
TGraphErrors* WeightedMean(TGraphErrors *gr1, TGraphErrors *gr2,
                           const double weight1) {
  if (gr1->GetN() != gr2->GetN()) {
    std::cout << "Unequal number of data points - aborting!\n";
    return nullptr;
  }
  auto grOut = (TGraphErrors*) gr1->Clone(
      Form("weightedMean_%i", int(gRandom->Uniform() * 10000.f)));
  static double x1, x2, y1, y2;
  for (int i = 0; i < gr1->GetN(); ++i) {
    gr1->GetPoint(i, x1, y1);
    gr2->GetPoint(i, x2, y2);
    if (std::abs(x1 - x2) > 1e-3) {
      std::cout << "x values are not equal - aborting!\n";
      return nullptr;
    }
    grOut->SetPoint(i, x1, weight1 * y1 + (1. - weight1) * y2);
    grOut->SetPointError(
        i,
        gr1->GetErrorX(i),
        std::sqrt(
            weight1 * weight1 * gr1->GetErrorY(i) * gr1->GetErrorY(i)
                + (1. - weight1) * (1. - weight1) * gr2->GetErrorY(i)
                    * gr2->GetErrorY(i)));
  }
  return grOut;
}

/// =====================================================================================
DLM_Ck* getDLMCk(TGraph *gr) {
  const double *x = gr->GetX();
  DLM_Ck *out = new DLM_Ck(gr->GetN() - 1, x[0], x[gr->GetN() - 1]);
  static double xVal, yVal;
  for (int i = 0; i < gr->GetN(); ++i) {
    gr->GetPoint(i, xVal, yVal);
    out->SetBinContent(i, yVal);
  }
  out->Update();
  return out;
}

/// =====================================================================================
DLM_Ck* getDLMCk(TGraph *gr, int nbins, double kmin, double kmax) {
  DLM_Ck *out = new DLM_Ck(nbins, kmin, kmax);
  for (int i = 0; i < nbins; ++i) {
    out->SetBinContent(i, gr->Eval(out->GetBinCenter(0, i)));
  }
  out->Update();
  return out;
}

/// =====================================================================================
TGraphErrors* getCkFromYuki(int potential, double rad = 0.9) {
  TGraphErrors* ingraph;
  double x, y, x2, y2;
  TString HomeDir = gSystem->GetHomeDirectory().c_str();
  if (potential == 3) {
    ingraph = new TGraphErrors(
        TString::Format(
            "%s/CERNHome/D-mesons/Analysis/Models/corr_model1_%.2f_fm_wC.dat",
            HomeDir.Data(), rad));
    for (int i = 0; i < ingraph->GetN(); ++i) {
      ingraph->GetPoint(i, x, y);
      ingraph->SetPoint(i, x, y);
      ingraph->SetPointError(i, 0., 0.);
    }
  } else if (potential == 4) {
    ingraph = new TGraphErrors(
        TString::Format(
            "%s/CERNHome/D-mesons/Analysis/Models/corr_model3_%.2f_fm_wC.dat",
            HomeDir.Data(), rad));
    for (int i = 0; i < ingraph->GetN(); ++i) {
      ingraph->GetPoint(i, x, y);
      ingraph->SetPoint(i, x, y);
      ingraph->SetPointError(i, 0., 0.);
    }
  } else if (potential == 5) {
    ingraph = new TGraphErrors(
        TString::Format(
            "%s/CERNHome/D-mesons/Analysis/Models/corr_model4_1_%.2f_fm_wC.dat",
            HomeDir.Data(), rad));
    for (int i = 0; i < ingraph->GetN(); ++i) {
      ingraph->GetPoint(i, x, y);
      ingraph->SetPoint(i, x, y);
      ingraph->SetPointError(i, 0., 0.);
    }
  } else if (potential == 6) {
    ingraph = new TGraphErrors(
        TString::Format(
            "%s/CERNHome/D-mesons/Analysis/Models/corr_model4_2_%.2f_fm_wC.dat",
            HomeDir.Data(), rad));
    for (int i = 0; i < ingraph->GetN(); ++i) {
      ingraph->GetPoint(i, x, y);
      ingraph->SetPoint(i, x, y);
      ingraph->SetPointError(i, 0., 0.);
    }
  } else {
    std::cout << "ERROR: getCkFromYuki - potential not available\n";
    return nullptr;
  }
  return ingraph;
}

/// =====================================================================================
TGraph* getCkPotential(int potVal, double rad) {
  TString HomeDir = gSystem->GetHomeDirectory().c_str();
  TString filename = HomeDir.Data();
  filename += "/CERNHome/D-mesons/Analysis/Models/corr_VI0_";
  filename += TString::Format("%.2f/corr", rad);
  if (potVal > 0) {
    filename += TString::Format("%04dMeV_", potVal);
  } else {
    filename += TString::Format("-%04dMeV_", std::abs(potVal));
  }
  filename += TString::Format("%.2ffm_wC.dat", rad);
  auto grOut = new TGraph();
  int count = 0;
  std::ifstream potFile;
  potFile.open(filename.Data());
  std::string line;
  double p, kstar, c1, c2;
  while (!potFile.eof()) {
    getline(potFile, line);
    std::istringstream is(line);
    while (is >> p >> kstar >> c1 >> c2) {
      grOut->SetPoint(count++, kstar, 1.f + c1 + c2);
    }
  }
  return grOut;
}

/// =====================================================================================
TH2F* CutRange(const TH2F *input) {
  auto HistExtendes = new TH2F(
      TString::Format("%s_extended", input->GetName()).Data(),
      input->GetTitle(), 595, 0.025, 3, 600, 0, 3);
  for (int iBinX = 1; iBinX <= input->GetNbinsX(); iBinX++) {
    for (int iBinY = 1; iBinY <= input->GetNbinsY(); iBinY++) {
      HistExtendes->Fill(input->GetXaxis()->GetBinCenter(iBinX),
                         input->GetYaxis()->GetBinCenter(iBinY),
                         input->GetBinContent(iBinX, iBinY));
    }
  }
  return HistExtendes;
}

/// =====================================================================================
TH2F* TransformToMeV(const TH2F *input) {
  auto histMeV = new TH2F(
      TString::Format("%s_MeV", input->GetName()).Data(), input->GetTitle(),
      input->GetNbinsX(), 1000. * input->GetXaxis()->GetBinLowEdge(1),
      1000. * input->GetXaxis()->GetBinUpEdge(input->GetNbinsX()),
      input->GetNbinsY(), 1000. * input->GetYaxis()->GetBinLowEdge(1),
      1000. * input->GetYaxis()->GetBinUpEdge(input->GetNbinsY()));
  for (int iBinX = 1; iBinX <= input->GetNbinsX(); iBinX++) {
    for (int iBinY = 1; iBinY <= input->GetNbinsY(); iBinY++) {
      histMeV->SetBinContent(iBinX, iBinY, input->GetBinContent(iBinX, iBinY));
    }
  }
  return histMeV;
}

/// =====================================================================================
template<typename T>
T getBootstrapFromVec(const std::vector<T> &vec) {
  return vec.at(gRandom->Uniform() * vec.size());
}

/// =====================================================================================
TGraphAsymmErrors* getBootstrapGraph(TGraphAsymmErrors *gr) {
  auto grOut = (TGraphAsymmErrors*) gr->Clone(
      Form("bootstrap_%s_%i", gr->GetName(),
           int(gRandom->Uniform() * 10000.f)));
  static double xVal, yVal;
  for (int i = 0; i < gr->GetN(); ++i) {
    gr->GetPoint(i, xVal, yVal);
    grOut->SetPoint(i, xVal, gRandom->Gaus(yVal, gr->GetErrorY(i)));
  }
  return grOut;
}

/// =====================================================================================
TGraphErrors* getBootstrapGraph(TGraphErrors *gr) {
  auto grOut = (TGraphErrors*) gr->Clone(
      Form("bootstrap_%s_%i", gr->GetName(),
           int(gRandom->Uniform() * 10000.f)));
  static double xVal, yVal;
  for (int i = 0; i < gr->GetN(); ++i) {
    gr->GetPoint(i, xVal, yVal);
    grOut->SetPoint(i, xVal, gRandom->Gaus(yVal, gr->GetErrorY(i)));
  }
  return grOut;
}

/// =====================================================================================
void getImprovedStartParamsPol3(TGraphAsymmErrors *gr, const double range,
                                std::vector<double> &params) {
  auto pol1 = new TF1("pol1", "pol1", 0, range);
  gr->Fit(pol1, "RQ");
  auto pol2 = new TF1("pol2", "pol2", 0, range);
  pol2->SetParameter(0, pol1->GetParameter(0));
  pol2->SetParameter(1, pol1->GetParameter(1));
  gr->Fit(pol2, "RQ");
  auto pol3 = new TF1("pol3", "pol3", 0, range);
  pol3->FixParameter(0, pol2->GetParameter(0));
  pol3->FixParameter(1, pol2->GetParameter(1));
  pol3->FixParameter(2, pol2->GetParameter(2));
  gr->Fit(pol3, "RQ");
  pol3->ReleaseParameter(2);
  gr->Fit(pol3, "RQ");
  pol3->ReleaseParameter(1);
  gr->Fit(pol3, "RQ");
  pol3->ReleaseParameter(0);
  gr->Fit(pol3, "RQ");
  params = { {pol3->GetParameter(0), pol3->GetParameter(1), pol3->GetParameter(
          2), pol3->GetParameter(3)}};
  delete pol3;
  delete pol2;
  delete pol1;
}

/// =====================================================================================
TGraphAsymmErrors* EvalBootstrap(TNtuple *tuple, TGraphAsymmErrors* grCF,
                                 TList* debug, TString OutputDir,
                                 TString potName) {
  auto grOut = new TGraphAsymmErrors(Form("gr_%s", tuple->GetName()));
  int count = 0;
  double kstar, yVal;
  for (int i = 0; i < grCF->GetN(); ++i) {
    grCF->GetPoint(i, kstar, yVal);
    tuple->Draw("cf >> h(1000, 0, 10)",
                Form("TMath::Abs(kstar - %.3f) < 4.99", kstar), "N");  // the bins are shifted for the syst. variations
    auto hist = (TH1F*) gROOT->FindObject("h");
    if (hist->GetEntries() == 0) {
      continue;
    }
    double error = hist->GetRMS();

    tuple->Draw("cf >> h2",
                Form("TMath::Abs(kstar - %.3f) < 4.99 && BootID == 0", kstar));  // the bins are shifted for the syst. variations
    auto hist2 = (TH1F*) gROOT->FindObject("h2");
    if (hist2->GetEntries() == 0) {
      continue;
    }
    double cfDefault = hist2->GetMean();
    grOut->SetPoint(count, kstar, cfDefault);
    grOut->SetPointError(count, grCF->GetErrorX(i), grCF->GetErrorX(i), error,
                         error);
    if (debug) {
      auto c = new TCanvas();
      DreamPlot::SetStyleHisto(hist);
      hist->SetTitle(
          Form("#it{k}* = %.1f MeV/#it{c} ;#Delta C(#it{k}*); Entries", kstar));
      hist->Draw();
      hist->GetXaxis()->SetRangeUser(hist->GetMean() - 10 * hist->GetRMS(),
                                     hist->GetMean() + 10 * hist->GetRMS());
      auto gr = new TGraphErrors();
      gr->SetMarkerStyle(20);
      gr->SetMarkerSize(1);
      gr->SetMarkerColor(kRed + 2);
      gr->SetLineColor(kRed + 2);
      gr->SetPoint(0, cfDefault, 0.5 * hist->GetMaximum());
      gr->SetPointError(0, error, 0);
      gr->Draw("pe same");
      c->Print(
          Form("%s/Debug/Debug_%s_%s_%i.pdf", OutputDir.Data(),
               tuple->GetName(), potName.Data(), int(kstar)));
      delete gr;
      delete c;
      debug->Add(hist);
    }

    delete hist;
    delete hist2;

    count++;
  }
  return grOut;
}

/// =====================================================================================
TGraphErrors* EvalBootstrap(TNtuple *tuple, TList* debug, TString OutputDir,
                            TString potName, const double kmin,
                            const double kmax, const double binWidth,
                            const bool useRMS = true) {
  auto grOut = new TGraphErrors(Form("gr_%s", tuple->GetName()));
  int count = 0;
  for (double kstar = kmin; kstar < kmax;) {
    tuple->Draw("cf >> h(10000, 0, 100)",
                Form("TMath::Abs(kstar - %.3f) < 0.01", kstar), "N");
    auto hist = (TH1F*) gROOT->FindObject("h");
    if (hist->GetEntries() == 0) {
      kstar += binWidth;
      continue;
    }
    double mean = 0;
    double error = 0;
    if (useRMS) {
      mean = hist->GetMean();
      error = hist->GetRMS();
    } else {
      const float binLow = hist->GetXaxis()->GetBinLowEdge(
          hist->FindFirstBinAbove(0.1, 1));
      const float binUp = hist->GetXaxis()->GetBinUpEdge(
          hist->FindLastBinAbove(0.1, 1));
      mean = (binUp + binLow) / 2.;
      error = std::abs(binLow - binUp) / std::sqrt(12);
    }

    grOut->SetPoint(count, kstar, mean);
    grOut->SetPointError(count, 0, error);
    if (debug) {
      auto c = new TCanvas();
      DreamPlot::SetStyleHisto(hist);
      hist->SetTitle(
          Form("#it{k}* = %.1f MeV/#it{c} ;#Delta C(#it{k}*); Entries", kstar));
      hist->Draw();
      hist->GetXaxis()->SetRangeUser(hist->GetMean() - 10 * hist->GetRMS(),
                                     hist->GetMean() + 10 * hist->GetRMS());
      auto gr = new TGraphErrors();
      gr->SetMarkerStyle(20);
      gr->SetMarkerSize(1);
      gr->SetMarkerColor(kRed + 2);
      gr->SetLineColor(kRed + 2);
      gr->SetPoint(0, mean, 0.5 * hist->GetMaximum());
      gr->SetPointError(0, error, 0);
      gr->Draw("pe same");
      c->Print(
          Form("%s/Debug/Debug_%s_%i.pdf", OutputDir.Data(), potName.Data(),
               int(kstar)));
      delete gr;
      delete c;
      debug->Add(hist);
    }

    delete hist;
    count++;
    kstar += binWidth;
  }
  return grOut;
}

/// =====================================================================================
TNtuple* EvalPotentials(TNtuple *tuple, std::vector<float> &pots,
                        std::vector<double> &sourceSizes, int nSystVars) {
  auto tupleOut = new TNtuple("potentialsChi2", "potentialsChi2",
                              "potVal:chi2:femtoRad:systID");
  int count = 0;
  for (const auto &potIt : pots) {
    for (const auto &radIt : sourceSizes) {
      for (int systIt = 0; systIt <= nSystVars; ++systIt) {
        tuple->Draw(
            "chi2 >> htemp(5000,0,50)",
            Form(
                "TMath::Abs(potVal - %f) < 1e-3 && TMath::Abs(femtoRad - %f) < 1e-3 && systID == %d",
                potIt, radIt, systIt),
            "N");
        auto hist = (TH1F*) gROOT->FindObject("htemp");
        if (hist->GetEntries() == 0) {
          continue;
        }
        tupleOut->Fill(potIt, hist->GetMean(), radIt, systIt);
      }
    }
  }
  return tupleOut;
}

#endif /* GENTLEKITTY_SCRIPTS_PDMESON_DMESONTOOLS_H_ */

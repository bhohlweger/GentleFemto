#include "TSystem.h"
#include "TROOT.h"
#include "TKey.h"
#include "TClass.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TString.h"
#include "TNtuple.h"
#include "CATSInputSigma0.h"
#include "DreamPlot.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TLegend.h"
#include <iostream>

const int tupleLength = 14;

void ComputeChi2(TH1F* dataHist, TGraphErrors *grFit, double &bestChi2,
                 double &defaultChi2, double &worstChi2, int &nBins) {
  double theoryX, theoryY, theoryErr;
  double dataY, dataErr, dataErrSq;
  double chi2Down = 0.;
  double chi2Def = 0.;
  double chi2Up = 0.;
  int maxkStarBinChi2 = dataHist->FindBin(250);
  for (unsigned uBin = 1; uBin <= maxkStarBinChi2; uBin++) {
    double mom = dataHist->GetBinCenter(uBin);
    grFit->GetPoint(uBin - 1, theoryX, theoryY);
    theoryErr = grFit->GetErrorY(uBin - 1);

    if (mom != theoryX) {
      std::cerr << "PROBLEM Sigma0 " << mom << '\t' << theoryX << std::endl;
    }
    dataY = dataHist->GetBinContent(uBin);
    dataErr = dataHist->GetBinError(uBin);
    dataErrSq = (dataErr * dataErr);

    chi2Up += (dataY - (theoryY + theoryErr)) * (dataY - (theoryY + theoryErr))
        / dataErrSq;
    chi2Def += (dataY - theoryY) * (dataY - theoryY) / dataErrSq;
    chi2Down += (dataY - (theoryY - theoryErr))
        * (dataY - (theoryY - theoryErr)) / dataErrSq;
    ++nBins;
  }
  worstChi2 = (chi2Up > chi2Down) ? chi2Up : chi2Down;
  defaultChi2 = chi2Def;
  bestChi2 = (chi2Up > chi2Down) ? chi2Down : chi2Up;

}

// =========================================
// Compute per k* bin the 68% variations
void EvalError(TNtuple *tuple, const int iBranches, TH1F* histCF,
               TGraphErrors *grOut, bool debugPlot, TString varFolder) {
  float kVal = histCF->GetBinCenter(iBranches);
  float CkExp = histCF->GetBinContent(iBranches);

  tuple->Draw(Form("delta%u >> h", iBranches));
  TH1F* hist = (TH1F*) gROOT->FindObject("h");

  double binLow = hist->GetXaxis()->GetBinLowEdge(
      hist->FindFirstBinAbove(1, 1));
  double binUp = hist->GetXaxis()->GetBinUpEdge(hist->FindLastBinAbove(1, 1));
  double DeltaCoulomb = (binLow - binUp) / TMath::Sqrt(12);
  double CoulombdefVal = (binUp + binLow) / 2.;

  if (debugPlot) {
    auto gr = new TGraphErrors();
    DreamPlot::SetStyleGraph(gr, 20, kRed + 2);
    gr->SetLineWidth(2);
    gr->SetPoint(0, CoulombdefVal, 1);
    gr->SetPointError(0, std::abs(DeltaCoulomb), 0);
    DreamPlot::SetStyleHisto(hist);
    hist->SetTitle(
        Form("#it{k}* = %.1f MeV/#it{c} ;#Delta C(#it{k}*); Entries", kVal));
    auto c = new TCanvas();
    hist->Draw();
    gr->Draw("pez same");
    c->SaveAs(
        Form("%s/Delta_%s_%i.pdf", varFolder.Data(), tuple->GetName(),
             iBranches));
    delete gr;
    delete c;
  }

  //std::cout << CkExp - binLow << " " << CkExp - binUp << " " << CoulombdefVal << "\n";
  grOut->SetPoint(iBranches - 1, kVal, CkExp - CoulombdefVal);
  grOut->SetPointError(iBranches - 1, 0, DeltaCoulomb);

  delete hist;
}

// =========================================
// Draw all systematic variations available
void DrawSigma(const unsigned& NumIter, TString varFolder, const int& potential,
               std::vector<double> params) {
  bool debugPlots = false;
  bool fancyPlot = false;
  double d0, REf0inv, IMf0inv, deltap0, deltap1, deltap2, etap0, etap1, etap2;

  DreamPlot::SetStyle(false, true);
  TString dataHistName = "hCk_ReweightedpSigma0MeV_0";
  TH1F* CF_Histo;

  TString graphfilename;
  if (potential == 0) {
    if (params.size() != 3) {
      std::cout << "ERROR: Wrong number of scattering parameters\n";
      return;
    }
    d0 = params[0];
    REf0inv = params[1];
    IMf0inv = params[2];

    graphfilename = TString::Format("%s/Param_pSigma0_%i_%.3f_%.3f_%.3f.root",
                                    varFolder.Data(), potential, d0, REf0inv,
                                    IMf0inv);
  } else if (potential == 1) {
    if (params.size() != 6) {
      std::cout << "ERROR: Wrong number of parameters for delta/eta\n";
      return;
    }
    deltap0 = params[0];
    deltap1 = params[1];
    deltap2 = params[2];
    etap0 = params[3];
    etap1 = params[4];
    etap2 = params[5];

    graphfilename = TString::Format(
        "%s/Param_pSigma0_%i_%.1f_%.4f_%.7f_%.2f_%.5f_%.8f.root",
        varFolder.Data(), potential, deltap0, deltap1, deltap2, etap0, etap1,
        etap2);

  } else {
    graphfilename = TString::Format("%s/Param_pSigma0_%i.root",
                                    varFolder.Data(), potential);
  }

  auto file = TFile::Open(graphfilename.Data(), "update");
  float tupler[tupleLength];
  auto fit = new TNtuple("fit", "fit",
                         "delta1:delta2:delta3:delta4:delta5:delta6:"
                         "delta7:delta8:delta9:delta10:delta11:"
                         "delta12:delta13:delta14");
  auto sideband = new TNtuple("sideband", "sideband",
                              "delta1:delta2:delta3:delta4:delta5:delta6:"
                              "delta7:delta8:delta9:delta10:delta11:"
                              "delta12:delta13:delta14");

  auto c1 = new TCanvas("CF_var", "CF_var");
  TIter next(file->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*) next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TDirectoryFile")) {
      continue;
    }
    auto dirFile = (TDirectoryFile*) key->ReadObj();
    TIter nextnext(dirFile->GetListOfKeys());
    TKey *nextkey;
    while ((nextkey = (TKey*) nextnext())) {

      auto hist = (TH1F*) nextkey->ReadObj();
      TString histName = Form("%s", hist->GetName());
      if (histName.Contains(dataHistName)) {
        CF_Histo = hist;
        c1->cd();
        CF_Histo->GetXaxis()->SetRangeUser(0, 600);
        CF_Histo->GetYaxis()->SetRangeUser(0.8, 1.6);
        CF_Histo->SetTitle("; #it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
        DreamPlot::SetStyleHisto(hist, 24, kBlack);
        CF_Histo->Draw();
      }
    }

    TIter nextnext2(dirFile->GetListOfKeys());
    while ((nextkey = (TKey*) nextnext2())) {

      auto Graph = (TGraph*) nextkey->ReadObj();
      TString graphName = Form("%s", Graph->GetName());

      if (graphName.Contains("Sigma0Sideband")) {
        c1->cd();
        Graph->SetLineColor(kGray + 2);
        Graph->SetLineStyle(0);
        Graph->Draw("L3same");
        for (int iPnt = 0; iPnt < tupleLength; ++iPnt) {
          tupler[iPnt] = CF_Histo->GetBinContent(iPnt + 1)
              - Graph->Eval(CF_Histo->GetBinCenter(iPnt + 1));
        }
        sideband->Fill(tupler);
      }

      if (graphName.Contains("pSigma0Graph")) {
        c1->cd();
        Graph->SetLineColor(kRed + 2);
        Graph->Draw("L3same");
        for (int iPnt = 0; iPnt < tupleLength; ++iPnt) {
          tupler[iPnt] = CF_Histo->GetBinContent(iPnt + 1)
              - Graph->Eval(CF_Histo->GetBinCenter(iPnt + 1));
        }
        fit->Fill(tupler);
      }
    }
  }

  if (debugPlots) {
    if (potential == 0) {
      c1->Print(
          Form("%s/CF_pSigma_model_%.3f_%.3f_%.3f.pdf", varFolder.Data(), d0,
               REf0inv, IMf0inv));
    } else if (potential == 1) {
      c1->Print(
          Form("%s/CF_pSigma_model_%.1f_%.4f_%.7f_%.2f_%.5f_%.8f.pdf",
               varFolder.Data(), deltap0, deltap1, deltap2, etap0, etap1,
               etap2));
    } else {
      c1->Print(Form("%s/CF_pSigma_model_%i.pdf", varFolder.Data(), potential));
    }
    CF_Histo->GetYaxis()->SetRangeUser(0.9, 1.1);
    if (potential == 0) {
      c1->Print(
          Form("%s/CF_pSigma_model_zoom_%.3f_%.3f_%.3f.pdf", varFolder.Data(),
               d0, REf0inv, IMf0inv));
    } else if (potential == 1) {
      c1->Print(
          Form("%s/CF_pSigma_model_zoom_%.1f_%.4f_%.7f_%.2f_%.5f_%.8f.pdf",
               varFolder.Data(), deltap0, deltap1, deltap2, etap0, etap1,
               etap2));
    } else {
      c1->Print(
          Form("%s/CF_pSigma_model_zoom_%i.pdf", varFolder.Data(), potential));
    }
  }

  file->cd();
  c1->Write();
  fit->Write();
  sideband->Write();

  auto grCF = new TGraphErrors();
  grCF->SetName("CF_fit");
  auto grSidebands = new TGraphErrors();
  grSidebands->SetName("CF_sidebands");

  for (int iBranches = 1; iBranches < tupleLength + 1; ++iBranches) {
    EvalError(fit, iBranches, CF_Histo, grCF, debugPlots, varFolder);
    EvalError(sideband, iBranches, CF_Histo, grSidebands, debugPlots,
              varFolder);
  }

  // in case we're running the exclusion task, we're interested in the best/worst chi2 for a given set of scattering parameters
  double bestChi2, defaultChi2, worstChi2;
  int nBins = 0;
  ComputeChi2(CF_Histo, grCF, bestChi2, defaultChi2, worstChi2, nBins);
  double pvalBest = TMath::Prob(bestChi2, round(nBins));
  double nSigmaBest = TMath::Sqrt(2) * TMath::ErfcInverse(pvalBest);
  double pvalDefault= TMath::Prob(defaultChi2, round(nBins));
  double nSigmaDefault = TMath::Sqrt(2) * TMath::ErfcInverse(pvalDefault);
  double pvalWorst = TMath::Prob(worstChi2, round(nBins));
  double nSigmaWorst = TMath::Sqrt(2) * TMath::ErfcInverse(pvalWorst);

  std::cout << "=============\n";
  std::cout << "nSigma maker for p-Sigma0 speaks \n";
  std::cout << "You have a nice potential, Sir!\n"
  std::cout << "  best nSigma  " << nSigmaBest << "\n";
  std::cout << "  def nSigma   " << nSigmaDefault << "\n";
  std::cout << "  worst nSigma " << nSigmaWorst << "\n";
  std::cout << "=============\n";

  double bestChi2Sideband, defaultChi2Sideband, worstChi2Sideband;
  int nBinsSideband = 0;
  ComputeChi2(CF_Histo, grSidebands, bestChi2Sideband, defaultChi2Sideband, worstChi2Sideband, nBinsSideband);
  double pvalBestSideband = TMath::Prob(bestChi2Sideband, round(nBinsSideband));
  double nSigmaBestSideband = TMath::Sqrt(2) * TMath::ErfcInverse(pvalBestSideband);
  double pvalDefaultSideband = TMath::Prob(defaultChi2Sideband, round(nBinsSideband));
  double nSigmaDefaultSideband = TMath::Sqrt(2) * TMath::ErfcInverse(pvalDefaultSideband);
  double pvalWorstSideband = TMath::Prob(worstChi2Sideband, round(nBinsSideband));
  double nSigmaWorstSideband = TMath::Sqrt(2) * TMath::ErfcInverse(pvalWorstSideband);

  std::cout << "Your sideband fits even better! \n";
  std::cout << "  best nSigma  " << nSigmaBestSideband << "\n";
  std::cout << "  def nSigma   " << nSigmaDefaultSideband << "\n";
  std::cout << "  worst nSigma " << nSigmaWorstSideband << "\n";
  std::cout << "Thanks for your input - well done!\n";
  std::cout << "=============\n";


  if (potential == 0) {
    auto fitTuple = (TNtuple*) file->Get("fitResult");

    TString resultfilename = TString::Format("%s/Result_%.3f_%.3f_%.3f.root",
                                             varFolder.Data(), d0, REf0inv,
                                             IMf0inv);
    auto resultFile = TFile::Open(resultfilename.Data(), "recreate");
    TNtuple* ntResult = new TNtuple(
        "exclusion", "exclusion",
        "CFneg:d0:REf0inv:IMf0inv:bestChi2:defChi2:worstChi2");

    Float_t ntBuffer[7];
    fitTuple->Draw("CFneg >> h");
    TH1F* hist = (TH1F*) gROOT->FindObject("h");
    ntBuffer[0] = hist->GetMean();
    ntBuffer[1] = d0;
    ntBuffer[2] = REf0inv;
    ntBuffer[3] = IMf0inv;
    ntBuffer[4] = bestChi2;
    ntBuffer[5] = defaultChi2;
    ntBuffer[6] = worstChi2;

    ntResult->Fill(ntBuffer);
    ntResult->Write();
    resultFile->Close();
  } else if (potential == 1) {
    auto fitTuple = (TNtuple*) file->Get("fitResult");

    TString resultfilename = TString::Format(
        "%s/Result_%.1f_%.4f_%.7f_%.2f_%.5f_%.8f.root", varFolder.Data(),
        deltap0, deltap1, deltap2, etap0, etap1, etap2);
    auto resultFile = TFile::Open(resultfilename.Data(), "recreate");
    TNtuple* ntResult =
        new TNtuple(
            "exclusion",
            "exclusion",
            "CFneg:deltap0:deltap1:deltap2:etap0:etap1:etap2:bestChi2:defChi2:worstChi2");

    Float_t ntBuffer[10];
    fitTuple->Draw("CFneg >> h");
    TH1F* hist = (TH1F*) gROOT->FindObject("h");
    ntBuffer[0] = hist->GetMean();
    ntBuffer[1] = deltap0;
    ntBuffer[2] = deltap1;
    ntBuffer[3] = deltap2;
    ntBuffer[4] = etap0;
    ntBuffer[5] = etap1;
    ntBuffer[6] = etap2;
    ntBuffer[7] = bestChi2;
    ntBuffer[8] = defaultChi2;
    ntBuffer[9] = worstChi2;

    ntResult->Fill(ntBuffer);
    ntResult->Write();
    resultFile->Close();
  }

  file->cd();
  auto c = new TCanvas("pSigma0Correlation");
  auto histDummy = new TH1F("histDummy",
                            "; #it{k}* (MeV/#it{c}); #it{C}(#it{k}*)", 100, 0,
                            450);
  histDummy->GetYaxis()->SetRangeUser(0.8, 1.6);
  histDummy->Draw();
  CF_Histo->Draw("same");
  grCF->Draw("L3same");
  grCF->SetLineColor(kRed + 2);
  grCF->SetFillColor(kRed + 2);
  grSidebands->Draw("l3 same");
  grSidebands->SetFillColorAlpha(kBlack, 0.5);
  grSidebands->SetLineColorAlpha(kBlack, 0.0);
  TLatex BeamText;
  BeamText.SetNDC(kTRUE);
  if (!fancyPlot) {
    BeamText.SetTextSize(0.8 * gStyle->GetTextSize());
    if (potential == 0) {
      BeamText.DrawLatex(0.45, 0.8, TString::Format("d_{0} = %.3f fm", d0));
      BeamText.DrawLatex(
          0.45, 0.73,
          TString::Format("#Rgothic(f_{0}^{-1}) = %.3f fm^{-1}", REf0inv));
      BeamText.DrawLatex(
          0.45, 0.66,
          TString::Format("#Jgothic(f_{0}^{-1}) = %.3f fm^{-1}", IMf0inv));
    } else if (potential == 1) {
      BeamText.DrawLatex(0.3, 0.8, "Phase shift");
      BeamText.DrawLatex(
          0.3, 0.73,
          TString::Format("%.2f, %.4f, %.6f", deltap0, deltap1, deltap2));
      BeamText.DrawLatex(0.3, 0.66, "Elasticity");
      BeamText.DrawLatex(
          0.3, 0.59, TString::Format("%.2f, %.4f, %.6f", etap0, etap1, etap2));
    }
    BeamText.DrawLatex(0.7, 0.8,
                       TString::Format("#chi^{2}_{best} = %.3f", bestChi2));
    BeamText.DrawLatex(0.7, 0.73,
                       TString::Format("#chi^{2}_{def} = %.3f", defaultChi2));
    BeamText.DrawLatex(0.7, 0.66,
                       TString::Format("#chi^{2}_{worst} = %.3f", worstChi2));
  } else {
    BeamText.SetTextSize(gStyle->GetTextSize());
    BeamText.DrawLatex(0.525, 0.8, "ALICE ");
    BeamText.DrawLatex(0.525, 0.75, "pp #sqrt{s} = 13 TeV HM");
    auto leg = new TLegend(0.51, 0.53, 0.85, 0.73);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(gStyle->GetTextSize());
    leg->AddEntry(CF_Histo, "p#Sigma^{0} #oplus #bar{p}#bar{#Sigma}^{0}", "pe");
    leg->AddEntry(grSidebands, "Background", "f");
    leg->AddEntry(grCF, "Femtoscopic fit", "f");
    leg->Draw("same");
  }

  if (debugPlots) {
    if (potential == 0) {
      c->Print(
          Form("%s/CF_pSigma_fit_%.3f_%.3f_%.3f.pdf", varFolder.Data(), d0,
               REf0inv, IMf0inv));
    } else if (potential == 1) {
      c->Print(
          Form("%s/CF_pSigma_fit_%.1f_%.4f_%.7f_%.2f_%.5f_%.8f.pdf",
               varFolder.Data(), deltap0, deltap1, deltap2, etap0, etap1,
               etap2));
    } else {
      c->Print(Form("%s/CF_pSigma_fit_%i.pdf", varFolder.Data(), potential));
    }
  }
  c->Write();
  grCF->Write();
  grSidebands->Write();

  delete histDummy;
  delete c;
  delete grSidebands;
  delete grCF;
  delete c1;
  delete sideband;
  delete fit;
  file->Close();
}

// =========================================
// Draw all systematic variations available
void DrawSigma(char *argv[]) {
  const unsigned& NumIter = atoi(argv[1]);
  TString varFolder = argv[5];
  const int potential = atoi(argv[6]);
  std::vector<double> params;
  if (potential == 0) {
    if (!argv[7] || !argv[8] || !argv[9]) {
      std::cout << "ERROR: Missing the scattering parameters\n";
      return;
    }
    params.push_back(atof(argv[7]));  // d0
    params.push_back(atof(argv[8]));  // REf0inv
    params.push_back(atof(argv[9]));  // IMf0inv
  } else if (potential == 1) {
    if (!argv[7] || !argv[8] || !argv[9] || !argv[10] || !argv[11]
        || !argv[12]) {
      std::cout << "ERROR: Missing the parameters for delta/eta\n";
      return;
    }
    params.push_back(atof(argv[7]));   // deltap0
    params.push_back(atof(argv[8]));   // deltap1
    params.push_back(atof(argv[9]));   // deltap2
    params.push_back(atof(argv[10]));  // etap0
    params.push_back(atof(argv[11]));  // etap1
    params.push_back(atof(argv[12]));  // etap2
  }
  DrawSigma(NumIter, varFolder, potential, params);
}

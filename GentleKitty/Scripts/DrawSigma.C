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

void nSigmaMaker(TNtuple* resultTuple, int upperRange, bool debugPlot,
                 TString varFolder, int potential,
                 std::array<double, 3> &nSigma) {

  resultTuple->Draw(Form("nSigma%i >> hist%i", upperRange, upperRange));
  TH1F* hist = (TH1F*) gROOT->FindObject(Form("hist%i", upperRange));

  const float binLow = hist->GetXaxis()->GetBinLowEdge(
      hist->FindFirstBinAbove(0.1, 1));
  const float binUp = hist->GetXaxis()->GetBinUpEdge(
      hist->FindLastBinAbove(0.1, 1));
  const float Delta = std::abs((binLow - binUp) / std::sqrt(12));
  const float nSigmaDefault = (binUp + binLow) / 2.;
  const float nSigmaBest = nSigmaDefault - Delta;
  const float nSigmaWorst = nSigmaDefault + Delta;
  std::cout << " - " << upperRange << " MeV/c cut-off: " << nSigmaBest << " / "
            << nSigmaDefault << " / " << nSigmaWorst << "\n";

  if (debugPlot) {
    auto c = new TCanvas();
    hist->Draw();
    hist->SetTitle(
        Form("; #it{n}_{#sigma} (#it{k}* < %i MeV/#it{c}); Entries",
             upperRange));

    TLatex nSigmaText;
    nSigmaText.SetNDC(kTRUE);
    nSigmaText.SetTextSize(0.8 * gStyle->GetTextSize());
    nSigmaText.DrawLatex(
        0.7, 0.8, TString::Format("#it{n}_{#sigma, best} = %.1f", nSigmaBest));
    nSigmaText.DrawLatex(
        0.7, 0.75,
        TString::Format("#it{n}_{#sigma, default} = %.1f", nSigmaDefault));
    nSigmaText.DrawLatex(
        0.7, 0.7,
        TString::Format("#it{n}_{#sigma, worst} = %.1f", nSigmaWorst));
    c->Print(
        Form("%s/nSigma%i_%i.pdf", varFolder.Data(), upperRange, potential));
    delete c;
  }
  nSigma[0] = nSigmaBest;
  nSigma[1] = nSigmaDefault;
  nSigma[2] = nSigmaWorst;
}

// =========================================
// Compute per k* bin the 68% variations
void EvalError(TNtuple *tuple, const int ikstar, TGraph* histCF,
               TGraphErrors *grOut, bool debugPlot, TString varFolder) {
  double kVal, CFval;
  histCF->GetPoint(ikstar, kVal, CFval);

  tuple->Draw("modelVal >> h", Form("std::abs(kstar - %.3f) < 1e-3", kVal));
  TH1F* hist = (TH1F*) gROOT->FindObject("h");

  double binLow = hist->GetXaxis()->GetBinLowEdge(
      hist->FindFirstBinAbove(0.1, 1));
  double binUp = hist->GetXaxis()->GetBinUpEdge(hist->FindLastBinAbove(0.1, 1));
  double DeltaCoulomb = (binLow - binUp) / TMath::Sqrt(12);
  double DefaultVal = (binUp + binLow) / 2.;

  if (debugPlot) {
    auto gr = new TGraphErrors();
    DreamPlot::SetStyleGraph(gr, 20, kRed + 2);
    gr->SetLineWidth(2);
    gr->SetPoint(0, DefaultVal, 1);
    gr->SetPointError(0, std::abs(DeltaCoulomb), 0);
    DreamPlot::SetStyleHisto(hist);
    hist->SetTitle(
        Form("#it{k}* = %.1f MeV/#it{c} ;#Delta C(#it{k}*); Entries", kVal));
    auto c = new TCanvas();
    hist->Draw();
    gr->Draw("pez same");
    c->SaveAs(
        Form("%s/Delta_%s_%i.pdf", varFolder.Data(), tuple->GetName(), ikstar));
    delete gr;
    delete c;
  }

  //std::cout << ikstar << " " << kVal << " " << CFval << " " << DefaultVal << "\n";
  grOut->SetPoint(ikstar, kVal, DefaultVal);
  grOut->SetPointError(ikstar, 0, DeltaCoulomb);

  delete hist;
}

// =========================================
// Draw all systematic variations available
void DrawSigma(const unsigned& NumIter, TString varFolder, const int& potential,
               std::vector<double> params) {
  bool debugPlots = false;
  double d0, REf0inv, IMf0inv, deltap0, deltap1, deltap2, etap0, etap1, etap2;

  DreamPlot::SetStyle(false, true);
  TString dataHistName = "hCk_ReweightedpSigma0MeV_0";
  TH1F* CF_Histo;
  TGraph* CF_Model;
  TGraph* CF_Sidebands;

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
  auto fit = new TNtuple("fit", "fit", "kstar:modelVal");
  auto sideband = new TNtuple("sideband", "sideband", "kstar:modelVal");

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
        if (!CF_Histo) {
          std::cout << "ERROR: No default histogram found!\n";
          return;
        }
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

      if (graphName.Contains("SBExtrapolate_")) {
        CF_Sidebands = Graph;
        c1->cd();
        Graph->SetLineColor(kGray + 2);
        Graph->SetLineStyle(0);
        Graph->Draw("L3same");
        double x, y;
        for (int iPnt = 0; iPnt < Graph->GetN(); ++iPnt) {
          Graph->GetPoint(iPnt, x, y);
          sideband->Fill(x, y);
        }
      }

      if (graphName.Contains("S0Extrapolate_")) {
        CF_Model = Graph;
        c1->cd();
        Graph->SetLineColor(kRed + 2);
        Graph->Draw("L3same");
        double x, y;
        for (int iPnt = 0; iPnt < Graph->GetN(); ++iPnt) {
          Graph->GetPoint(iPnt, x, y);
          fit->Fill(x, y);
        }
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

  auto resultTuple = (TNtuple*) file->Get("fitResult");
  std::array<double, 3> nSigma250;
  std::array<double, 3> nSigma200;
  std::array<double, 3> nSigma150;

  std::cout << "=============\n";
  std::cout << "This is nSigma maker for p-Sigma0 speaking\n";
  std::cout << "I gracefully acknowledge your input\n";
  std::cout << "You have a nice potential, Sir!\n";
  nSigmaMaker(resultTuple, 250, debugPlots, varFolder, potential, nSigma250);
  nSigmaMaker(resultTuple, 200, debugPlots, varFolder, potential, nSigma200);
  nSigmaMaker(resultTuple, 150, debugPlots, varFolder, potential, nSigma150);
  std::cout << "Thanks for your interest in potential " << potential << "\n";
  std::cout << "=============\n";

  auto grCF = new TGraphErrors();
  grCF->SetName("CF_fit");
  auto grSidebands = new TGraphErrors();
  grSidebands->SetName("CF_sidebands");

  double x, y;
  for (int ikstar = 0; ikstar < CF_Model->GetN(); ++ikstar) {
    CF_Model->GetPoint(ikstar, x, y);
    EvalError(fit, ikstar, CF_Model, grCF, debugPlots, varFolder);
    EvalError(sideband, ikstar, CF_Model, grSidebands, debugPlots, varFolder);
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
  BeamText.SetTextSize(0.8 * gStyle->GetTextSize());
  BeamText.DrawLatex(
      0.6,
      0.8,
      TString::Format("n#sigma_{250} = %.1f / %.1f / %.1f", nSigma250[0],
                      nSigma250[1], nSigma250[2]));
  BeamText.DrawLatex(
      0.6,
      0.73,
      TString::Format("n#sigma_{200} = %.1f / %.1f / %.1f", nSigma200[0],
                      nSigma200[1], nSigma200[2]));
  BeamText.DrawLatex(
      0.6,
      0.66,
      TString::Format("n#sigma_{150} = %.1f / %.1f / %.1f", nSigma150[0],
                      nSigma150[1], nSigma150[2]));
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
  return;
}

// =========================================
// Draw all systematic variations available
void DrawSigma(char *argv[]) {
  const unsigned& NumIter = atoi(argv[1]);
  TString varFolder = argv[6];
  const int potential = atoi(argv[7]);
  std::vector<double> params;
  if (potential == 0) {
    if (!argv[8] || !argv[9] || !argv[10]) {
      std::cout << "ERROR: Missing the scattering parameters\n";
      return;
    }
    params.push_back(atof(argv[8]));  // d0
    params.push_back(atof(argv[9]));  // REf0inv
    params.push_back(atof(argv[10]));  // IMf0inv
  } else if (potential == 1) {
    if (!argv[8] || !argv[9] || !argv[10] || !argv[11] || !argv[12]
        || !argv[13]) {
      std::cout << "ERROR: Missing the parameters for delta/eta\n";
      return;
    }
    params.push_back(atof(argv[8]));   // deltap0
    params.push_back(atof(argv[9]));   // deltap1
    params.push_back(atof(argv[10]));   // deltap2
    params.push_back(atof(argv[11]));  // etap0
    params.push_back(atof(argv[12]));  // etap1
    params.push_back(atof(argv[13]));  // etap2
  }
  DrawSigma(NumIter, varFolder, potential, params);
}

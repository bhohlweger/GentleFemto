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

  tuple->Draw("modelVal >> h", Form("TMath::Abs(kstar - %.3f) < 1e-1", kVal));
  TH1F* hist = (TH1F*) gROOT->FindObject("h");
  if(hist->GetEntries() == 0) return;

  double binLow = hist->GetXaxis()->GetBinLowEdge(
      hist->FindFirstBinAbove(0.1, 1));
  double binUp = hist->GetXaxis()->GetBinUpEdge(hist->FindLastBinAbove(0.1, 1));
  double DeltaCoulomb = (binLow - binUp) / TMath::Sqrt(12);
  double DefaultVal = (binUp + binLow) / 2.;

  if (debugPlot && kVal < 350) {
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

  grOut->SetPoint(ikstar, kVal, DefaultVal);
  grOut->SetPointError(ikstar, 0, DeltaCoulomb);

  delete hist;
}

// =========================================
// Draw all systematic variations available
void DrawSigma(TString varFolder, const int& potential) {
  bool debugPlots = false;
  double d0, REf0inv, IMf0inv, deltap0, deltap1, deltap2, etap0, etap1, etap2;

  DreamPlot::SetStyle(false, true);
  TString dataHistName = "Graph_from_hCk_ReweightedpSigma0_0MeV";
  TString sidebandHistName = "SidebandMerged_0";
  TGraphAsymmErrors* CF_Histo;
  TGraphAsymmErrors* CF_Sideband;
  TGraph* CF_Model;
  TGraph* CF_Sidebands;

  TString graphfilename = TString::Format("%s/Param_pSigma0_%i.root",
                                    varFolder.Data(), potential);
  auto file = TFile::Open(graphfilename.Data(), "update");
  auto fit = new TNtuple("fit", "fit", "kstar:modelVal");
  auto sideband = new TNtuple("sideband", "sideband", "kstar:modelVal");
  auto genuineSideband = new TNtuple("genuineSideband", "genuineSideband", "kstar:modelVal");

  auto c1 = new TCanvas("CF_var", "CF_var");
  auto c2 = new TCanvas("CF_sideband_var", "CFsideband_var");
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

      auto gr = (TGraphAsymmErrors*) nextkey->ReadObj();
      auto hist = (TGraphAsymmErrors*) nextkey->ReadObj();
      TString histName = Form("%s", gr->GetName());
      if (histName.Contains(dataHistName)) {
        CF_Histo = gr;
        if (!CF_Histo) {
          std::cout << "ERROR: No default histogram found!\n";
          return;
        }
        c1->cd();
        CF_Histo->GetXaxis()->SetRangeUser(0, 600);
        CF_Histo->GetYaxis()->SetRangeUser(0.8, 1.7);
        CF_Histo->SetTitle("; #it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
        DreamPlot::SetStyleGraph(CF_Histo, 24, kBlue + 3);
        CF_Histo->Draw("APEZ");
      } else if (histName.Contains(sidebandHistName)) {
        c2->cd();
        CF_Sideband = hist;
        CF_Sideband->GetXaxis()->SetRangeUser(0, 600);
        CF_Sideband->GetYaxis()->SetRangeUser(0.8, 1.7);
        CF_Sideband->SetTitle("; #it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
        DreamPlot::SetStyleGraph(CF_Sideband, 24, kBlue + 3);
        CF_Sideband->Draw("APEZ");
        if (!CF_Sideband) {
          std::cout << "ERROR: No default sideband histogram found!\n";
          return;
        }
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

      if (graphName.Contains("GenuineSideBand_")) {
        double x, y;
        for (int iPnt = 0; iPnt < Graph->GetN(); ++iPnt) {
          c2->cd();
          Graph->GetPoint(iPnt, x, y);
          Graph->SetLineColor(kGray + 2);
          Graph->SetLineStyle(0);
          Graph->Draw("L3same");
          genuineSideband->Fill(x, y);
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

  c1->cd();
  CF_Histo->Draw("pezsame");
  c2->cd();
  CF_Sideband->Draw("pezsame");

  if (debugPlots) {
    c1->Print(Form("%s/CF_pSigma_model_%i.pdf", varFolder.Data(), potential));
    c2->Print(Form("%s/CF_Sideband_model_%i.pdf", varFolder.Data(), potential));
    CF_Histo->GetYaxis()->SetRangeUser(0.9, 1.1);
    CF_Sideband->GetYaxis()->SetRangeUser(0.9, 1.1);
    c1->Print(
        Form("%s/CF_pSigma_model_zoom_%i.pdf", varFolder.Data(), potential));
    c2->Print(
        Form("%s/CF_Sideband_model_zoom_%i.pdf", varFolder.Data(), potential));
  }

  file->cd();
  c1->Write();
  fit->Write();
  sideband->Write();
  genuineSideband->Write();

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
  auto grGenuineSidebands = new TGraphErrors();
  grGenuineSidebands->SetName("CF_genuineSidebands");

  double x, y;
  for (int ikstar = 0; ikstar < CF_Model->GetN(); ++ikstar) {
    CF_Model->GetPoint(ikstar, x, y);
    EvalError(fit, ikstar, CF_Model, grCF, debugPlots, varFolder);
    EvalError(sideband, ikstar, CF_Model, grSidebands, debugPlots, varFolder);
    EvalError(genuineSideband, ikstar, CF_Model, grGenuineSidebands, debugPlots, varFolder);
  }

  file->cd();
  auto c = new TCanvas("pSigma0Correlation");
  auto histDummy = new TH1F("histDummy",
                            "; #it{k}* (MeV/#it{c}); #it{C}(#it{k}*)", 100, 0,
                            350);
  histDummy->GetYaxis()->SetRangeUser(0.8, 1.6);
  histDummy->Draw();
  grCF->Draw("L3same");
  grCF->SetLineColor(kRed + 2);
  grCF->SetFillColor(kRed + 2);
  grSidebands->Draw("l3 same");
  grSidebands->SetFillColorAlpha(kBlack, 0.5);
  grSidebands->SetLineColorAlpha(kBlack, 0.0);
  CF_Histo->Draw("pezsame");
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
  c->Print(Form("%s/CF_pSigma_fit_%i.pdf", varFolder.Data(), potential));

  auto d = new TCanvas("pSigma0Sideband");
  histDummy->Draw();
  grGenuineSidebands->Draw("L3same");
  grGenuineSidebands->SetLineColor(kRed + 2);
  grGenuineSidebands->SetFillColor(kRed + 2);
  CF_Sideband->Draw("pezsame");
  d->Print(Form("%s/CF_pSigma_sideband_%i.pdf", varFolder.Data(), potential));
  c->Write();
  grCF->Write();
  grSidebands->Write();
  grGenuineSidebands->Write();
  CF_Sideband->Write();

  delete histDummy;
  delete c;
  delete d;
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
  TString OutputDir = argv[4];
  const int potential = atoi(argv[5]);
  DrawSigma(OutputDir, potential);
}

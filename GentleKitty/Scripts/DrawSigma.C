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
                 double &defaultChi2, double &worstChi2) {
  double theoryX, theoryY, theoryErr;
  double dataY, dataErr, dataErrSq;
  double chi2Best, chi2Worst;
  double chi2Down, chi2Def, chi2Up;
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
  }
  worstChi2 = (chi2Up > chi2Down) ? chi2Up : chi2Down;
  defaultChi2 = chi2Def;
  bestChi2 = (chi2Up > chi2Down) ? chi2Down : chi2Up;

}

// =========================================
// Compute per k* bin the 68% variations
void EvalError(TNtuple *tuple, const int iBranches, TH1F* histCF,
               TGraphErrors *grOut) {
  float kVal = histCF->GetBinCenter(iBranches);
  float CkExp = histCF->GetBinContent(iBranches);

  tuple->Draw(Form("delta%u >> h", iBranches));
  TH1F* hist = (TH1F*) gROOT->FindObject("h");

  double binLow = hist->GetXaxis()->GetBinLowEdge(
      hist->FindFirstBinAbove(1, 1));
  double binUp = hist->GetXaxis()->GetBinUpEdge(hist->FindLastBinAbove(1, 1));
  double DeltaCoulomb = (binLow - binUp) / TMath::Sqrt(12);
  double CoulombdefVal = (binUp + binLow) / 2.;

  //std::cout << CkExp - binLow << " " << CkExp - binUp << " " << CoulombdefVal << "\n";
  grOut->SetPoint(iBranches - 1, kVal, CkExp - CoulombdefVal);
  grOut->SetPointError(iBranches - 1, 0, DeltaCoulomb);

  delete hist;
}

// =========================================
// Draw all systematic variations available
void DrawSigma(const unsigned& NumIter, TString InputDir, TString appendix,
               TString varFolder, const bool& isExclusion, const float d0,
               const float f0inv) {
  bool batchmode = true;

  DreamPlot::SetStyle();
  CATSInputSigma0 *CATSinput = new CATSInputSigma0();
  CATSinput->ReadSigma0CorrelationFile(InputDir.Data(), appendix.Data());
  CATSinput->ObtainCFs(10, 340, 440);
  TString dataHistName = "hCk_ReweightedpSigma0MeV_0";
  auto CF_Histo = CATSinput->GetCF("pSigma0", dataHistName.Data());
  CF_Histo->GetXaxis()->SetRangeUser(0, 600);
  CF_Histo->GetYaxis()->SetRangeUser(0.8, 1.6);
  CF_Histo->SetTitle("; #it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  DreamPlot::SetStyleHisto(CF_Histo, 24, kBlack);

  TString graphfilename;
  if (isExclusion) {
    graphfilename = TString::Format("%s/Param_pSigma0_%i_%.3f_%.3f.root",
                                    varFolder.Data(), NumIter, d0, f0inv);
  } else {
    graphfilename = TString::Format("%s/Param_pSigma0_%i.root",
                                    varFolder.Data(), NumIter);

  }

  float tupler[tupleLength];
  auto fit = new TNtuple("fit", "fit",
                         "delta1:delta2:delta3:delta4:delta5:delta6:"
                         "delta7:delta8:delta9:delta10:delta11:"
                         "delta12:delta13:delta14");
  auto sideband = new TNtuple("sideband", "sideband",
                              "delta1:delta2:delta3:delta4:delta5:delta6:"
                              "delta7:delta8:delta9:delta10:delta11:"
                              "delta12:delta13:delta14");
  std::vector<TGraph*> Graph;
  auto file = TFile::Open(graphfilename.Data(), "update");
  auto c1 = new TCanvas("CF_var", "CF_var");
  CF_Histo->Draw();
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

  if (!batchmode) {
    if (isExclusion) {
      c1->Print(
          Form("%s/CF_pSigma_model_%.3f_%.3f.pdf", varFolder.Data(), d0,
               f0inv));
    } else {
      c1->Print(Form("%s/CF_pSigma_model.pdf", varFolder.Data()));
    }
    CF_Histo->GetYaxis()->SetRangeUser(0.9, 1.1);
    if (isExclusion) {
      c1->Print(
          Form("%s/CF_pSigma_model_zoom_%.3f_%.3f.pdf", varFolder.Data(), d0,
               f0inv));
    } else {
      c1->Print(Form("%s/CF_pSigma_model_zoom.pdf", varFolder.Data()));
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
    EvalError(fit, iBranches, CF_Histo, grCF);
    EvalError(sideband, iBranches, CF_Histo, grSidebands);
  }

  // in case we're running the exclusion task, we're interested in the best/worst chi2 for a given set of scattering parameters
  if (isExclusion) {
    auto fitTuple = (TNtuple*) file->Get("fitResult");

    TString resultfilename = TString::Format("%s/Result_%.3f_%.3f.root",
                                             varFolder.Data(), d0, f0inv);
    auto resultFile = TFile::Open(resultfilename.Data(), "recreate");
    TNtuple* ntResult = new TNtuple(
        "exclusion", "exclusion", "CFneg:d0:f0inv:bestChi2:defChi2:worstChi2");
    double bestChi2, defaultChi2, worstChi2;
    ComputeChi2(CF_Histo, grCF, bestChi2, defaultChi2, worstChi2);

    Float_t ntBuffer[6];
    fitTuple->Draw("CFneg >> h");
    TH1F* hist = (TH1F*) gROOT->FindObject("h");
    ntBuffer[0] = hist->GetMean();
    ntBuffer[1] = d0;
    ntBuffer[2] = f0inv;
    ntBuffer[3] = bestChi2;
    ntBuffer[4] = defaultChi2;
    ntBuffer[5] = worstChi2;

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
  BeamText.SetTextSize(gStyle->GetTextSize());
  BeamText.SetNDC(kTRUE);
  BeamText.DrawLatex(0.525, 0.8, "ALICE ");
  BeamText.DrawLatex(0.525, 0.75, "pp #sqrt{s} = 13 TeV HM");
  auto leg = new TLegend(0.51, 0.53, 0.85, 0.73);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(gStyle->GetTextSize());
  leg->AddEntry(CF_Histo, "p#Sigma^{0} #oplus #bar{p}#bar{#Sigma}^{0}", "pe");
  leg->AddEntry(grCF, "Femtoscopic fit", "f");
  leg->AddEntry(grSidebands, "Background", "f");
  leg->Draw("same");

  if (!batchmode) {
    if (isExclusion) {
      c->Print(
          Form("%s/CF_pSigma_fit_%.3f_%.3f.pdf", varFolder.Data(), d0, f0inv));
    } else {
      c->Print(Form("%s/CF_pSigma_fit.pdf", varFolder.Data()));
    }
  }
  c->Write();
  grCF->Write();
  grSidebands->Write();
}

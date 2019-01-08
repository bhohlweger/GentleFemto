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
#include <iostream>

const int tupleLength = 14;

// =========================================
// Draw all systematic variations available
void DrawSigma(TString InputDir, TString appendix, TString varFolder) {
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
  const char* beg = "Graph";
  const char* ext = ".root";

  char* dir = gSystem->ExpandPathName(varFolder.Data());
  void* dirp = gSystem->OpenDirectory(dir);
  TString MygraphName = "";

  float tupler[tupleLength];
  auto fit = new TNtuple("fit", "fit",
                         "delta1:delta2:delta3:delta4:delta5:delta6:"
                         "delta7:delta8:delta9:delta10:delta11:"
                         "delta12:delta13:delta14");
  auto sideband = new TNtuple("sideband", "sideband",
                              "delta1:delta2:delta3:delta4:delta5:delta6:"
                              "delta7:delta8:delta9:delta10:delta11:"
                              "delta12:delta13:delta14");

  const char* entry;
  std::vector<TGraph*> Graph;
  Int_t n = 0;
  TString str;
  auto c1 = new TCanvas();
  CF_Histo->Draw();
  TNtuple *outCoulomb = new TNtuple("fitResults", "coulomb", "delta");
  while ((entry = (char*) gSystem->GetDirEntry(dirp))) {
    str = entry;
    if (str.EndsWith(ext) && str.BeginsWith(beg)) {
      auto file = TFile::Open(gSystem->ConcatFileName(dir, entry));
      TIter next(file->GetListOfKeys());
      TKey *key;
      while ((key = (TKey*) next())) {
        TClass *cl = gROOT->GetClass(key->GetClassName());
        if (!cl->InheritsFrom("TGraph")) {
          continue;
        }
        TGraph *Graph = (TGraph*) key->ReadObj();
        TString graphName = Form("%s", Graph->GetName());

        if (graphName.Contains("Sigma0Sideband")) {
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
          Graph->SetLineColor(kRed + 2);
          Graph->Draw("L3same");
          for (int iPnt = 0; iPnt < tupleLength; ++iPnt) {
            tupler[iPnt] = CF_Histo->GetBinContent(iPnt + 1)
                - Graph->Eval(CF_Histo->GetBinCenter(iPnt + 1));
          }
          fit->Fill(tupler);
        }
      }

      file->Close();
    }
  }

  c1->SaveAs(Form("%s/CF_pSigma_model.pdf", varFolder.Data()));
  CF_Histo->GetYaxis()->SetRangeUser(0.9, 1.1);
  c1->SaveAs(Form("%s/CF_pSigma_model_zoom.pdf", varFolder.Data()));
  auto output = TFile::Open(Form("%s/outfile.root", varFolder.Data()),
                            "RECREATE");
  output->cd();
  c1->Write();
  fit->Write();
  sideband->Write();
  output->Close();
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
// Obtain per k* bin the 68% variations for signal and sidebands
void EvalError(TString InputDir, TString appendix, TString varFolder) {
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
  auto output = TFile::Open(Form("%s/outfile.root", varFolder.Data()),
                            "update");
  auto fit = (TNtuple*) output->Get("fit");
  auto sideband = (TNtuple*) output->Get("sideband");

  auto grCF = new TGraphErrors();
  grCF->SetName("CF_fit");
  auto grSidebands = new TGraphErrors();
  grSidebands->SetName("CF_sidebands");

  for (int iBranches = 1; iBranches < tupleLength + 1; ++iBranches) {
    EvalError(fit, iBranches, CF_Histo, grCF);
    EvalError(sideband, iBranches, CF_Histo, grSidebands);
  }

  auto c = new TCanvas();
  auto histDummy = new TH1F("histDummy", "; #it{k}* (MeV/#it{c}); #it{C}(#it{k}*)", 100, 0, 450);
  histDummy->GetYaxis()->SetRangeUser(0.8, 1.6);
  histDummy->Draw();
  CF_Histo->Draw("same");
  grCF->Draw("L3same");
  grCF->SetLineColor(kRed + 2);
  grCF->SetFillColor(kRed + 2);
  grSidebands->Draw("l3 same");
  grSidebands->SetFillColorAlpha(kBlack, 0.5);
  grSidebands->SetLineColorAlpha(kBlack, 0.0);
  c->SaveAs(Form("%s/CF_pSigma_fit.pdf", varFolder.Data()));
  output->cd();
  c->Write();
  grCF->Write();
  grSidebands->Write();
}

// =========================================
// main
int main(int argc, char *argv[]) {
  DrawSigma(argv[1], argv[2], argv[3]);
  EvalError(argv[1], argv[2], argv[3]);
  return 0;
}

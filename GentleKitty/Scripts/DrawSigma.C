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
#include <iostream>

void DrawSigma(TString InputDir, TString appendix, TString varFolder) {
  DreamPlot::SetStyle();
  CATSInputSigma0 *CATSinput = new CATSInputSigma0();
  CATSinput->ReadSigma0CorrelationFile(InputDir.Data(), appendix.Data());
  CATSinput->ObtainCFs(10, 340, 440);
  TString dataHistName = "hCk_ReweightedpSigma0MeV_0";
  auto CF_Histo = CATSinput->GetCF("pSigma0", dataHistName.Data());
  CF_Histo->GetXaxis()->SetRangeUser(0, 700);
  CF_Histo->GetYaxis()->SetRangeUser(0.8, 1.6);
  CF_Histo->SetTitle("; #it{k}* (MeV/#it{c}); #it{C}(#it{k}*)");
  DreamPlot::SetStyleHisto(CF_Histo, 24, kBlack);
  const char* beg = "Graph";
  const char* ext = ".root";

  char* dir = gSystem->ExpandPathName(varFolder.Data());
  void* dirp = gSystem->OpenDirectory(dir);
  TString MygraphName = "";

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

        if (graphName.Contains("pSigma0Graph")) {
          Graph->SetLineColor(kRed + 2);
          Graph->Draw("L3same");
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
  output->Close();
}

int main(int argc, char *argv[]) {
  DrawSigma(argv[1], argv[2], argv[3]);
  return 0;
}

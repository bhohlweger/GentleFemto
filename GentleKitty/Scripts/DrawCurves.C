#include "TSystem.h"
#include "TROOT.h"
#include "TKey.h"
#include "TClass.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TString.h"
#include "CATSInput.h"
#include "DreamPlot.h"
#include <iostream>

void DrawCurves(const char* pair, const char* cfpath, const char* varFolder) {
  CATSInput *CATSinput = new CATSInput();
  CATSinput->ReadCorrelationFile(cfpath);
  CATSinput->ObtainCFs(5, 240, 340);
  TString HistName = Form("hCk_Reweighted%sMeV_0", pair);
  TH1F* CF_Histo = CATSinput->GetCF(pair, HistName);
  CF_Histo->GetXaxis()->SetRangeUser(0, 500);
  CF_Histo->SetTitle("; #it{k}* (GeV/#it{c}); #it{C}(#it{k}*)");
  CF_Histo->SetStats(false);
  DreamPlot::SetStyle();
  DreamPlot::SetStyleHisto(CF_Histo);
  const char* beg = "GraphFile";
  const char* ext = ".root";

  char* dir = gSystem->ExpandPathName(varFolder);
  void* dirp = gSystem->OpenDirectory(dir);
  TString MygraphName = "";
  if (TString(pair) == "pp") {
    MygraphName += "Fit";
  } else if (TString(pair) == "pXi") {
    MygraphName += "pXimGraph";
  }
  const char* entry;
  std::vector<TGraph*> Graph;
  Int_t n = 0;
  TString str;
  TCanvas* c1 = new TCanvas("c1", "c1");
  c1->cd(0);
  c1->SetCanvasSize(1920, 1280);
  c1->SetMargin(0.15, 0.05, 0.2, 0.05);
  CF_Histo->Draw();
  while ((entry = (char*) gSystem->GetDirEntry(dirp))) {
    str = entry;
//    std::cout << str << std::endl;
    if (str.EndsWith(ext) && str.BeginsWith(beg)) {
      TFile* file = TFile::Open(gSystem->ConcatFileName(dir, entry));
      TIter next(file->GetListOfKeys());
      TKey *key;
      while ((key = (TKey*) next())) {
        TClass *cl = gROOT->GetClass(key->GetClassName());
        if (!cl->InheritsFrom("TGraph"))
          continue;
        TGraph *Graph = (TGraph*) key->ReadObj();
        TString graphName = Form("%s", Graph->GetName());
        if (graphName.Contains(MygraphName.Data())) {

          if (TString(pair) == "pXi") {
            if (graphName.Contains("COULOMB")) {
              Graph->SetLineColor(3);
            } else {
              Graph->SetLineColor(2);
            }
          } else {
            Graph->SetLineColor(2);
          }
          Graph->Draw("CPsame");
        }
      }
      file->Close();
    }
  }
  c1->SaveAs(Form("%s/c1.png", varFolder));
  TFile* output = TFile::Open(Form("%s/outfile.root", varFolder), "RECREATE");
  output->cd();
  c1->Write();
  output->Close();
}

int main(int argc, char *argv[]) {
  DrawCurves(argv[1], argv[2], argv[3]);
  return 0;
}
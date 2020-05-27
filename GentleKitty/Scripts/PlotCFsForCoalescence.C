#include "TROOT.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TFile.h"
#include "DreamPlot.h"

int main(int argc, char *argv[]) {
  if (!argv[1]) {
    std::cout << "pp RadFile missing\n";
    return -1;
  }
  if (!argv[2]) {
    std::cout << "pp RadFile missing\n";
    return -1;
  }
  if (!argv[3]) {
    std::cout << "Source Name\n";
    return -1;
  }
  const char *ppFile = argv[1];
  const char *pLNLOFile = argv[2];
// const char* pLLOFile = argv[3];
  const char *sourceName = argv[3];
  DreamPlot::SetStyle();
  gStyle->SetHatchesSpacing(0.5);

  TFile *CF_ppFile1 = TFile::Open(ppFile, "read");
  TGraphErrors *GraphCF1 = (TGraphErrors*) CF_ppFile1->Get(
      "Graph_from_hCk_Rebinned_0MeV");
  std::cout << "reached1\n";
  TFile *CF_ppFile2 = TFile::Open(pLNLOFile, "read");
  TGraphErrors *GraphCF2 = (TGraphErrors*) CF_ppFile2->Get(
      "Graph_from_hCk_Rebinned_0MeV");

  if (!GraphCF1) {
    std::cout << "Graph_from_hCk_Rebinned_MeV0 from first file missing\n";
    return -1;
  }
  if (!GraphCF2) {
    std::cout << "Graph_from_hCk_Rebinned_MeV0 from second file missing\n";
    return -1;
  }

  TFile *out = TFile::Open(Form("%s.root", sourceName), "recreate");
  out->cd();
  auto c4 = new TCanvas("c8", "c8", 1200, 800);
  c4->cd();
  TLegend *leg = new TLegend(0.4, 0.5, 0.73, 0.78);
  leg->SetFillStyle(0);
  leg->SetTextFont(43);
//  leg->SetNColumns(2);
  leg->SetTextSizePixels(40);

  GraphCF1->SetTitle("; #it{k}*  (MeV); #it{C}(k*)");
  GraphCF1->GetXaxis()->SetTitleSize(40);
  GraphCF1->GetYaxis()->SetTitleSize(40);
  GraphCF1->GetXaxis()->SetTitleOffset(1.35);
  GraphCF1->GetYaxis()->SetTitleOffset(1.4);
  GraphCF1->GetXaxis()->SetLabelSize(40);
  GraphCF1->GetYaxis()->SetLabelSize(40);
  GraphCF1->GetXaxis()->SetLabelOffset(.02);
  GraphCF1->GetYaxis()->SetLabelOffset(.02);
  GraphCF1->SetMarkerColor(kRed + 1);
  GraphCF1->SetLineColor(kRed + 1);
  GraphCF1->SetFillColorAlpha(kRed + 1, 0.7);
  GraphCF1->SetMarkerStyle(21);
  GraphCF1->SetMarkerSize(1.0);
  GraphCF1->GetXaxis()->SetLimits(0.0, 400);
  GraphCF1->GetYaxis()->SetRangeUser(0.5, 4.5);
  GraphCF1->Draw("AP");

  GraphCF2->SetMarkerColor(kBlue + 2);
  GraphCF2->SetLineColor(kBlue + 2);
  GraphCF2->SetLineWidth(1);
  GraphCF2->SetMarkerStyle(21);
  GraphCF2->SetMarkerSize(1.0);
  GraphCF2->Draw("pez same");

  leg->AddEntry(GraphCF1,
                "p#minus#kern[0.4]{p} Avg Ref. Mult_{#eta < 0.8} ~30.1", "pef");
  leg->AddEntry(GraphCF2, "p#minus#kern[0.4]{p} Avg Ref. Mult_{#eta < 0.8} ~23",
                "pef");

  leg->Draw("same");
  c4->SaveAs(Form("%s/CFPlots.pdf", gSystem->pwd()));
  c4->Write();
  out->Write();
  out->Close();
  CF_ppFile1->Close();
  CF_ppFile2->Close();
}

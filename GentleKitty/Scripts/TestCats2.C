#include "TFile.h"
#include "TGraph.h"
#include "TF1.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TCanvas.h"

void testCats() {
  TFile* myfile =
      TFile::Open(
          "~/cernbox/HM13TeV/AnalysisData/ClosePairRej/SelectedPairs/fmTOutputProper.root",
          "read");
  TGraph* ppSource = (TGraph*) myfile->Get("hRad_pp");
  ppSource->SetTitle("Source using a Gaussian core plus resonances");
  ppSource->SetLineWidth(5);
  ppSource->SetLineColor(kBlue);
  ppSource->GetXaxis()->SetRangeUser(0., 8.);
  ppSource->GetXaxis()->SetTitle("r (fm)");
  ppSource->GetXaxis()->SetTitleSize(0.05);
  ppSource->GetXaxis()->SetTitleOffset(0.95);
  ppSource->GetXaxis()->SetLabelSize(0.05);
//  ppSource->GetXaxis()->SetLabelOffset(0.95);

  ppSource->GetXaxis()->CenterTitle(true);
  ppSource->GetYaxis()->SetRangeUser(0., 0.55);
  ppSource->GetYaxis()->SetTitle("S(r) (fm^{-1})");
  ppSource->GetYaxis()->CenterTitle(true);
  ppSource->GetYaxis()->SetTitleSize(0.05);
  ppSource->GetYaxis()->SetTitleOffset(0.95);
  ppSource->GetYaxis()->SetLabelSize(0.05);
//  ppSource->GetYaxis()->SetLabelOffset(0.95);

  TGraph* pLSource = (TGraph*) myfile->Get("hRad_pL");
  pLSource->SetLineWidth(5);
  pLSource->SetLineColor(30);
  TGraph* pS0Source = (TGraph*) myfile->Get("hRad_pS0");
  pS0Source->SetLineWidth(5);
  pS0Source->SetLineColor(38);
  TGraph* pXiSource = (TGraph*) myfile->Get("hRad_pXim");
  pXiSource->SetLineWidth(5);
  pXiSource->SetLineColor(46);
  TGraph* pOmSource = (TGraph*) myfile->Get("hRad_pOmega");
  pOmSource->SetLineWidth(5);
  pOmSource->SetLineColor(93);

  TF1* ppGauss = (TF1*) myfile->Get("pp");
  ppGauss->SetLineWidth(5);
  ppGauss->SetLineStyle(5);
  ppGauss->SetLineColor(kBlue);
  ppGauss->SetNpx(500);
  TF1* pLGauss = (TF1*) myfile->Get("pL");
  pLGauss->SetLineWidth(5);
  pLGauss->SetLineStyle(5);
  pLGauss->SetLineColor(30);
  pLGauss->SetNpx(500);
  TF1* pS0Gauss = (TF1*) myfile->Get("pS0");
  pS0Gauss->SetLineWidth(5);
  pS0Gauss->SetLineStyle(5);
  pS0Gauss->SetLineColor(38);
  pS0Gauss->SetNpx(500);
  TF1* pXimGauss = (TF1*) myfile->Get("pXi");
  pXimGauss->SetLineWidth(5);
  pXimGauss->SetLineStyle(5);
  pXimGauss->SetLineColor(46);
  pXimGauss->SetNpx(500);
  TF1* pOmGauss = (TF1*) myfile->Get("pOm");
  pOmGauss->SetLineWidth(5);
  pOmGauss->SetLineStyle(5);
  pOmGauss->SetLineColor(93);
  pOmGauss->SetNpx(500);

  auto* c1 = new TCanvas("IamBeautiful", "Iambautiful");
  c1->cd();
  ppSource->Draw("APL");
  pLSource->Draw("same");
  pS0Source->Draw("same");
  pXiSource->Draw("same");
  pOmSource->Draw("same");

  ppGauss->Draw("same");
  pLGauss->Draw("same");
  pS0Gauss->Draw("same");
  pXimGauss->Draw("same");
  pOmGauss->Draw("same");

  TLegend* leg = new TLegend(0.4, 0.4, 0.87, 0.87);
  leg->SetBorderSize(0);
  leg->SetFillStyle(4000);
  leg->AddEntry(
      ppSource,
      TString::Format(
          "p#minusp #LT#it{m}_{T}#GT =1.35 GeV/#it{c}^{2} (R_{G,#it{eff}} = %.2f fm)",
          ppGauss->GetParameter(0)),
      "l");
  leg->AddEntry(
      pLSource,
      TString::Format(
          "p#minus#Lambda #LT#it{m}_{T}#GT=1.55 GeV/#it{c}^{2} (R_{G,#it{eff}} = %.2f fm)",
          pLGauss->GetParameter(0)),
      "l");
  leg->AddEntry(
      pS0Source,
      TString::Format(
          "p#minus#Sigma^{0} #LT#it{m}_{T}#GT=2.07 GeV/#it{c}^{2} (R_{G,#it{eff}} = %.2f fm)",
          pS0Gauss->GetParameter(0)),
      "l");
  leg->AddEntry(
      pXiSource,
      TString::Format(
          "p#minus#Xi^{#minus} #LT#it{m}_{T}#GT=1.85 GeV/#it{c}^{2} (R_{G,#it{eff}} = %.2f fm)",
          pXimGauss->GetParameter(0)),
      "l");
  leg->AddEntry(
      pOmSource,
      TString::Format(
          "p#minus#Omega^{#minus} #LT#it{m}_{T}#GT=2.17 GeV/#it{c}^{2} (R_{G,#it{eff}} = %.2f fm)",
          pOmGauss->GetParameter(0)),
      "l");
  leg->Draw("same");

  TFile* output = new TFile(
      "~/cernbox/HM13TeV/AnalysisData/ClosePairRej/SelectedPairs/bautiful.root",
      "RECREATE");
  c1->SaveAs("Beautiful.png");
  c1->SaveAs("Beautiful.pdf");
  output->cd();
  c1->Write();
  output->Close();
  return;
}

int main(int argc, char *argv[]) {
  testCats();
  return 0;
}

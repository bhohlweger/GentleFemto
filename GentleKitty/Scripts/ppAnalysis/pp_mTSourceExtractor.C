#include "TFile.h"
#include "TGraphErrors.h"
#include "TNtuple.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <iostream>
#include "TLegend.h"
void mTRadiusExtractor(const char *inputFile, float &radius, float &radErrSys, float &radErrStat) {
  TFile* inFile = TFile::Open(inputFile, "READ");
  if(!inFile) {
    std::cout << "No input found for " << inputFile << ", exiting \n";
    return;
  }
  TNtuple* mTSysVar = (TNtuple*) inFile->Get("ntResult");
  auto histRad = new TH1D("hRad", "hRad", 200, 0.4, 2.0);
  auto histErrRad = new TH1D("hRadErr", "hRadErr", 50, 0., 0.1);
  mTSysVar->Draw("Radius_pp>>hRad");
  mTSysVar->Draw("RadiusErr_pp>>hRadErr");
  radius = histRad->GetMean();
  radErrSys = histRad->GetRMS();
  radErrStat = histErrRad->GetMean();
  delete histErrRad;
  delete histRad;
}

void mTRadiusPlot(const char* inputFolder, const char* avgmTFile) {
  TFile* mTFile = TFile::Open(avgmTFile, "READ");
  TGraphErrors* avgmT = (TGraphErrors*) mTFile->Get("AveragemT");
  TGraphErrors* mTRadiusSyst = new TGraphErrors();
  TGraphErrors* mTRadiusStat = new TGraphErrors();
  for (int imT = 0; imT < 7; ++imT) {
    double mT, dummy;
    float radius, radErrSyst, radErrStat;
    TString inputFile = TString::Format("%s/OutFileVarpp_%u.root",
                                        inputFolder, imT);
    avgmT->GetPoint(imT, dummy, mT);
//    std::cout << avgmT->GetErrorY(imT) << std::endl;
    mTRadiusExtractor(inputFile.Data(), radius, radErrSyst, radErrStat);
    std::cout << "radius: " << radius << " radErrSyst: " << radErrSyst << " radErrStat: " << radErrStat << std::endl;
    mTRadiusSyst->SetPoint(imT, mT, radius);
    mTRadiusStat->SetPoint(imT, mT, radius);
    mTRadiusSyst->SetPointError(imT, avgmT->GetErrorY(imT), radErrSyst);
    mTRadiusStat->SetPointError(imT, avgmT->GetErrorY(imT), radErrStat);
  }
  auto c1 = new TCanvas("c2","c2");
  mTRadiusSyst->SetLineColor(kBlack);
  mTRadiusSyst->SetTitle("; < m_{T} >  (MeV/#it{c}^{2}); r_{Core} (fm)");

  mTRadiusSyst->GetXaxis()->SetTitleSize(0.05);
  mTRadiusSyst->GetYaxis()->SetTitleSize(0.05);
  mTRadiusSyst->GetXaxis()->SetTitleOffset(0.9);
  mTRadiusSyst->GetYaxis()->SetTitleOffset(0.9);

  mTRadiusSyst->GetXaxis()->SetLabelSize(0.05);
  mTRadiusSyst->GetYaxis()->SetLabelSize(0.05);
  mTRadiusSyst->GetXaxis()->SetLabelOffset(0.004);
  mTRadiusSyst->GetYaxis()->SetLabelOffset(0.004);

  mTRadiusSyst->SetMarkerColor(kBlack);
  mTRadiusSyst->SetLineWidth(3);
  mTRadiusSyst->Draw("Ap");
  mTRadiusSyst->SetFillColorAlpha(kBlack, 0.4);
  mTRadiusSyst->Draw("2 same");

  TLegend* leg = new TLegend(0.6,0.6,0.9,0.9);
  leg->AddEntry(mTRadiusSyst,"p#minus p (AV18)","lef");

  mTRadiusStat->SetMarkerColor(kBlack);
  mTRadiusStat->SetLineWidth(3);
  mTRadiusStat->Draw("pe same");
  leg->Draw("same");
  c1->SaveAs("mTvsRad.pdf");
  TFile* output = TFile::Open("mTRad.root", "RECREATE");
  output->cd();
  c1->Write();
  mTRadiusSyst->SetName("mTRadiusSyst");
  mTRadiusSyst->Write();
  mTRadiusStat->SetName("mTRadiusStat");
  mTRadiusStat->Write();
  output->Close();
}

int main(int argc, char *argv[]) {
  mTRadiusPlot(argv[1], argv[2]);
  return 0;
}

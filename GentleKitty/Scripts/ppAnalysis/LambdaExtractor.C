#include <iostream>
#include "TH1F.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TMath.h"

int main(int argc, char* argv[]) {
  const char* baseDir = argv[1];
  const int nMTBins = atoi(argv[2]);
  const char* pair = argv[3];
  const char* dirSuffix = argv[4] == "" ? "" : argv[4];
  TGraphErrors lam_pp;
  TGraphErrors lam_ppL;
  TGraphErrors lam_pL;
  TGraphErrors lam_pLS0;
  TGraphErrors lam_pLXim;
  for (int imT = 0; imT < nMTBins; ++imT) {
    TFile* output = TFile::Open(
        TString::Format("%s/mTBin_%u/%s/OutFileVar%s.root", baseDir, imT + 1,
                        dirSuffix, pair).Data(),
        "read");
    if (!output) {
      std::cout
          << "File "
          << TString::Format("%s/mTBin_%u/%s/OutFileVar%s.root", baseDir,
                             imT + 1, dirSuffix, pair).Data()
          << " not found \n";
      return -1;
    }
    TTree* ppTree = (TTree*) output->FindObjectAny("ppTree");
    if (!ppTree) {
      std::cout <<"Tree not found \n";
      return -1;
    }
    ppTree->Draw("mTValue>>mTHist");
    double mTValue = ((TH1F*) gROOT->FindObject("mTHist"))->GetMean(1);
    std::cout << "mT:" << mTValue << std::endl;
    ppTree->Draw("lam_pp>>hlam_pp");
    TH1F* hlam_pp = (TH1F*) gROOT->FindObject("hlam_pp");

    double binLow = hlam_pp->GetXaxis()->GetBinLowEdge(
        hlam_pp->FindFirstBinAbove(0.1, 1));
    double binUp = hlam_pp->GetXaxis()->GetBinUpEdge(
        hlam_pp->FindLastBinAbove(0.1, 1));
    double Delta = TMath::Abs((binLow - binUp)) / TMath::Sqrt(12);
    double DefaultVal = (binUp + binLow) / 2.;
    lam_pp.SetPoint(imT, mTValue, DefaultVal);
    lam_pp.SetPointError(imT, 0, Delta);

    ppTree->Draw("lam_ppL>>hlam_ppL");
    TH1F* hlam_ppL = (TH1F*) gROOT->FindObject("hlam_ppL");

    binLow = hlam_ppL->GetXaxis()->GetBinLowEdge(
        hlam_ppL->FindFirstBinAbove(0.1, 1));
    binUp = hlam_ppL->GetXaxis()->GetBinUpEdge(
        hlam_ppL->FindLastBinAbove(0.1, 1));
    Delta = TMath::Abs((binLow - binUp)) / TMath::Sqrt(12);
    DefaultVal = (binUp + binLow) / 2.;
    lam_ppL.SetPoint(imT, mTValue, DefaultVal);
    lam_ppL.SetPointError(imT, 0, Delta);

    ppTree->Draw("lam_pL>>hlam_pL");
    TH1F* hlam_pL = (TH1F*) gROOT->FindObject("hlam_pL");

    binLow = hlam_pL->GetXaxis()->GetBinLowEdge(
        hlam_pL->FindFirstBinAbove(0.1, 1));
    binUp = hlam_pL->GetXaxis()->GetBinUpEdge(
        hlam_pL->FindLastBinAbove(0.1, 1));
    Delta = TMath::Abs((binLow - binUp)) / TMath::Sqrt(12);
    DefaultVal = (binUp + binLow) / 2.;
    lam_pL.SetPoint(imT, mTValue, DefaultVal);
    lam_pL.SetPointError(imT, 0, Delta);

    ppTree->Draw("lam_pLS0>>hlam_pLS0");
    TH1F* hlam_pLS0 = (TH1F*) gROOT->FindObject("hlam_pLS0");

    binLow = hlam_pLS0->GetXaxis()->GetBinLowEdge(
        hlam_pLS0->FindFirstBinAbove(0.1, 1));
    binUp = hlam_pLS0->GetXaxis()->GetBinUpEdge(
        hlam_pLS0->FindLastBinAbove(0.1, 1));
    Delta = TMath::Abs((binLow - binUp)) / TMath::Sqrt(12);
    DefaultVal = (binUp + binLow) / 2.;
    lam_pLS0.SetPoint(imT, mTValue, DefaultVal);
    lam_pLS0.SetPointError(imT, 0, Delta);

    ppTree->Draw("lam_pLXim>>hlam_pLXim");
    TH1F* hlam_pLXim = (TH1F*) gROOT->FindObject("hlam_pLXim");

    binLow = hlam_pLXim->GetXaxis()->GetBinLowEdge(
        hlam_pLXim->FindFirstBinAbove(0.1, 1));
    binUp = hlam_pLXim->GetXaxis()->GetBinUpEdge(
        hlam_pLXim->FindLastBinAbove(0.1, 1));
    Delta = TMath::Abs((binLow - binUp)) / TMath::Sqrt(12);
    DefaultVal = (binUp + binLow) / 2.;
    lam_pLXim.SetPoint(imT, mTValue, DefaultVal);
    lam_pLXim.SetPointError(imT, 0, Delta);

    output->Close();

  }
  lam_pp.SetName("Graph_pp");
  lam_pp.SetMarkerColor(kRed);
  lam_pp.SetLineColor(kRed);

  lam_ppL.SetName("Graph_ppL");
  lam_ppL.SetMarkerColor(kRed);
  lam_ppL.SetLineColor(kRed);

  lam_pL.SetName("Graph_pL");
  lam_pL.SetMarkerColor(kRed);
  lam_pL.SetLineColor(kRed);

  lam_pLS0.SetName("Graph_pLS0");
  lam_pLS0.SetMarkerColor(kRed);
  lam_pLS0.SetLineColor(kRed);

  lam_pLXim.SetName("Graph_pLXim");
  lam_pLXim.SetMarkerColor(kRed);
  lam_pLXim.SetLineColor(kRed);

  TFile* out = TFile::Open(TString::Format("out_%s.root", pair), "recreate");
  out->cd();
  lam_pp.Write("Graph_pp",0);
  lam_ppL.Write("Graph_ppL",0);
  lam_pL.Write("Graph_pL",0);
  lam_pLS0.Write("Graph_pLS0",0);
  lam_pLXim.Write("Graph_pLXim",0);
  out->Close();
}

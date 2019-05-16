#include "TFile.h"
#include "TGraphErrors.h"
#include "TNtuple.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <iostream>
#include "TLegend.h"
#include <vector>
void mTRadiusExtractor(int iBin, float &radius, float &radErrSys, float &radErrStat) {
  // Resonances
  // mT Bin 1: Radius Mean: 1.13497 Radius stat Err: 0.00996232 Radius SystErr Down: 0.0198902 Radius Syst Err Up: 0.0225297
  // mT Bin 2: Radius Mean: 1.07765 Radius stat Err: 0.0116435 Radius SystErr Down: 0.0146171 Radius Syst Err Up: 0.019332
  // mT Bin 3: Radius Mean: 1.04481 Radius stat Err: 0.0228825 Radius SystErr Down: 0.015478 Radius Syst Err Up: 0.0150418
  // mT Bin 4: Radius Mean: 0.990884 Radius stat Err: 0.0147561 Radius SystErr Down: 0.00988509 Radius Syst Err Up: 0.0143956
  // mT Bin 5: Radius Mean: 0.924195 Radius stat Err: 0.0139078 Radius SystErr Down: 0.0117057 Radius Syst Err Up: 0.0160959
  // mT Bin 6: Radius Mean: 0.826865 Radius stat Err: 0.00853009 Radius SystErr Down: 0.0126212 Radius Syst Err Up: 0.0163956
  // mT Bin 7: Radius Mean: 0.731138 Radius stat Err: 0.0136254 Radius SystErr Down: 0.0125193 Radius Syst Err Up: 0.0147957
//    std::vector<float> mean =     {1.13497,   1.07765,  1.04481,  0.990884,   0.924195,   0.826865,   0.731138};
//    std::vector<float> statErr =  {0.00996232,0.0116435,0.0228825,0.0147561,  0.0139078,  0.00853009, 0.0136254};
//    std::vector<float> systDown = {0.0198902, 0.0146171,0.015478, 0.00948221, 0.00988509, 0.0126212,  0.0125193};
//    std::vector<float> systUp=    {0.0225297, 0.019332, 0.0150418,0.0145662,  0.0143956,  0.0163956,  0.0147957};
  // Gauss Only
  // mT Bin 1: Radius Mean: 1.48545 Radius stat Err: 0.0164265 Radius SystErr Down: 0.0247313 Radius Syst Err Up: 0.0303012
  // mT Bin 2: Radius Mean: 1.42313 Radius stat Err: 0.0233733 Radius SystErr Down: 0.0181501 Radius Syst Err Up: 0.0265942
  // mT Bin 3: Radius Mean: 1.38728 Radius stat Err: 0.0284505 Radius SystErr Down: 0.0150145 Radius Syst Err Up: 0.0245082
  // mT Bin 4: Radius Mean: 1.32564 Radius stat Err: 0.026964 Radius SystErr  Down: 0.0117462 Radius Syst Err Up: 0.014972
  // mT Bin 5: Radius Mean: 1.25229 Radius stat Err: 0.0246837 Radius SystErr Down: 0.0130806 Radius Syst Err Up: 0.0160792
  // mT Bin 6: Radius Mean: 1.14325 Radius stat Err: 0.0183403 Radius SystErr Down: 0.01267 Radius Syst Err Up: 0.0135278
  ///mT Bin 7: Radius Mean: 1.03737 Radius stat Err: 0.0101932 Radius SystErr Down: 0.0144771 Radius Syst Err Up: 0.0148329
    std::vector<float> mean =     {1.48545,   1.42313,  1.38728,  1.32564,    1.25229,     1.14325,   1.03737};
    std::vector<float> statErr =  {0.0164265, 0.0233733,0.0284505,0.026964,   0.0246837,  0.0183403,  0.0101932};
    std::vector<float> systDown = {0.0247313, 0.0181501,0.0150145,0.0117462,  0.0130806,  0.01267,    0.0144771};
    std::vector<float> systUp=    {0.0303012, 0.0265942,0.0245082,0.014972,   0.0160792,  0.0135278,  0.0148329};
    radius = mean.at(iBin);
    radErrStat = statErr.at(iBin);
    radErrSys = (systDown.at(iBin)+systUp.at(iBin))/2.;
}

void mTRadiusPlot(const char* inputFolder, const char* avgmTFile) {
  TFile* mTFile = TFile::Open(avgmTFile, "READ");
  if (!mTFile) {
    std::cout << "No mT File \n";
    return;
  }
  TGraphErrors* avgmT = (TGraphErrors*) mTFile->Get("AveragemT_ppVar0");
  if(! avgmT) {
    std::cout << "no average mT file " << std::endl;
    return;
  }

  TGraphErrors* mTRadiusSyst = new TGraphErrors();
  TGraphErrors* mTRadiusStat = new TGraphErrors();
  for (int imT = 0; imT < 7; ++imT) {
    double mT, dummy;
    float radius, radErrSyst, radErrStat;
    TString inputFile = TString::Format("%s/OutFileVarpp_%u.root",
                                        inputFolder, imT);
    avgmT->GetPoint(imT, dummy, mT);
//    std::cout << avgmT->GetErrorY(imT) << std::endl;
    mTRadiusExtractor(imT, radius, radErrSyst, radErrStat);
    std::cout << "radius: " << radius << " radErrSyst: " << radErrSyst << " radErrStat: " << radErrStat << std::endl;
    if (radius < 0) {
      continue;
    }
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

//  mTRadiusSyst->GetXaxis()->SetRangeUser(0.95,2.7);
//  mTRadiusSyst->GetYaxis()->SetRangeUser(0.65,1.2);
  mTRadiusSyst->GetXaxis()->SetRangeUser(0.95,2.7);
  mTRadiusSyst->GetYaxis()->SetRangeUser(0.95,1.55);

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



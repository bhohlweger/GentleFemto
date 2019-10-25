/*
 * LambdaGami.cxx
 *
 *  Created on: Oct 24, 2019
 *      Author: schmollweger
 */

#include "LambdaGami.h"
#include "TError.h"
#include <iostream>

LambdaGami::LambdaGami() {

}

LambdaGami::~LambdaGami() {
}

TH1F* LambdaGami::UnfoldResidual(TH1F* cf, TH1F* res, double lamRes) {
  if (!cf) {
    Warning("LambdaGami::UnfoldResidual", "No Input hist");
    return nullptr;
  }
  const unsigned int nBins = cf->GetNbinsX();
  TString outName = TString::Format("%s_ResGami", cf->GetName());
  TH1F* outHist = new TH1F(outName.Data(), outName.Data(), nBins,
                           cf->GetXaxis()->GetXmin(),
                           cf->GetXaxis()->GetXmax());
  for (int iBin = 1; iBin < nBins + 1; ++iBin) {
    if (TMath::Abs(cf->GetBinCenter(iBin) - res->GetBinCenter(iBin)) > 1e-3) {
      Error(
          "LambdaGami::UnfoldResidual",
          "Difference between binning of genuine and residual contribution \n");
      std::cout << "iBin: " << iBin << " cf->GetBinCenter(iBin): "
                << cf->GetBinCenter(iBin) << " res->GetBinCenter(iBin): "
                << res->GetBinCenter(iBin) << std::endl;
      return nullptr;
    }
    double content = cf->GetBinContent(iBin);
    double residual = res->GetBinContent(iBin) - 1;
    residual /= lamRes;
    outHist->SetBinContent(iBin, content - residual);
  }
  return outHist;
}

TH1F* LambdaGami::UnfoldGenuine(TH1F* cf, double lamGen) {
  if (!cf) {
    Warning("LambdaGami::UnfoldGenuine", "No Input hist");
    return nullptr;
  }
  const unsigned int nBins = cf->GetNbinsX();
  TString outName = TString::Format("%s_GenuineGami", cf->GetName());
  TH1F* outHist = new TH1F(outName.Data(), outName.Data(), nBins,
                           cf->GetXaxis()->GetXmin(),
                           cf->GetXaxis()->GetXmax());
  for (int iBin = 1; iBin < nBins + 1; ++iBin) {
    double content = cf->GetBinContent(iBin) - 1;
    content /= lamGen;
    content += 1;
    outHist->SetBinContent(iBin, content);
  }
  return outHist;
}

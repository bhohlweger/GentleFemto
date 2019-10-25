/*
 * LambdaGami.cxx
 *
 *  Created on: Oct 24, 2019
 *      Author: schmollweger
 */

#include "LambdaGami.h"
#include "TError.h"

LambdaGami::LambdaGami() {

}

LambdaGami::~LambdaGami() {
}

TH1F* LambdaGami::UnfoldResidual(TH1F* cf, TH1F* res, double lamRes) {
  if (!cf) {
    Warning("LambdaGami::UnfoldLambda", "No Input hist");
    return nullptr;
  }
  const unsigned int nBins = cf->GetNbinsX() + 1;
  TString outName = TString::Format("%s_ResGami", cf->GetName());
  TH1F* outHist = new TH1F(outName.Data(), outName.Data(), nBins,
                           cf->GetXaxis()->GetXbins()->GetArray());
  for (int iBin = 1; iBin < nBins; ++iBin) {
    if (TMath::Abs(cf->GetBinCenter(iBin) - res->GetBinCenter(iBin)) > 1e-3) {
      Error(
          "LambdaGami::UnfoldLambda",
          "Difference between binning of genuine and residual contribution \n");
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
    Warning("LambdaGami::UnfoldLambda", "No Input hist");
    return nullptr;
  }
  const unsigned int nBins = cf->GetNbinsX() + 1;
  TString outName = TString::Format("%s_GenuineGami", cf->GetName());
  TH1F* outHist = new TH1F(outName.Data(), outName.Data(), nBins,
                           cf->GetXaxis()->GetXbins()->GetArray());
  for (int iBin = 1; iBin < nBins; ++iBin) {
    double content = cf->GetBinContent(iBin) - 1;
    content /= lamGen;
    content += 1;
    outHist->SetBinContent(iBin, content);
  }
  return outHist;
}

/*
 * LambdaGami.cxx
 *
 *  Created on: Oct 24, 2019
 *      Author: schmollweger
 */

#include "LambdaGami.h"
#include "TError.h"

#include <iostream>

LambdaGami::LambdaGami()
    : fRelError(nullptr) {
}

LambdaGami::~LambdaGami() {
}

TH1F* LambdaGami::UnfoldResidual(TH1F* cf, TH1F* res, double lamRes) {
  if (!cf) {
    Error("LambdaGami::UnfoldResidual", "No Input hist");
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
    residual *= lamRes;
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

void LambdaGami::StoreStatErr(TH1F* cfMeasured) {
  if (fRelError) {
    delete fRelError;
  }
  TString Name = TString::Format("%sRelErr", cfMeasured->GetName());
  fRelError = (TH1F*) cfMeasured->Clone(Name.Data());
  for (int iBin = 1; iBin < fRelError->GetNbinsX() + 1; ++iBin) {
    fRelError->SetBinContent(
        iBin, fRelError->GetBinError(iBin) / fRelError->GetBinContent(iBin));
    fRelError->SetBinError(iBin, 0);
  }
}

void LambdaGami::AddStatErr(TH1F* cfMeasured) {
  if (!fRelError) {
    Error("LambdaGami::AddStatErr","No rel error stored!\n");
    return;
  }
  for (int iBin = 1; iBin < fRelError->GetNbinsX() + 1; ++iBin) {
    cfMeasured->SetBinError(
        iBin, cfMeasured->GetBinContent(iBin) * fRelError->GetBinContent(iBin));
  }
}

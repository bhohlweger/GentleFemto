/*
 * MomentumGami.cxx
 *
 *  Created on: Oct 24, 2019
 *      Author: schmollweger
 */

#include "MomentumGami.h"
#include "TError.h"
#include "TF1.h"
#include "TRandom.h"

MomentumGami::MomentumGami(float maxkStar)
    : fResolution(nullptr),
      fResProjection(),
      fToUnfold(nullptr),
      fMaxkStar(maxkStar) {
  // TODO Auto-generated constructor stub
  gRandom->SetSeed(0);
}

MomentumGami::~MomentumGami() {
}

void MomentumGami::SetResolution(TH2F* resoMatrix, float UnitConversion) {
  //User needs to make sure that the target binning is set properly!
  TString hname;
  fResolution = resoMatrix;
  std::cout << "resoMatrix->GetXaxis()->FindBin(fMaxkStar * UnitConversion): "
            << resoMatrix->GetXaxis()->FindBin(fMaxkStar * UnitConversion)
            << std::endl;
  for (int i = 0;
      i < resoMatrix->GetXaxis()->FindBin(fMaxkStar * UnitConversion); i++) {
    hname = Form("hprojY%i", i);
    fResProjection.push_back(
        (TH1F*) resoMatrix->ProjectionY(hname, i + 1, i + 1));
    //normalize
    fResProjection[i]->Scale(
        1. / fResProjection[i]->Integral(1, fResProjection[i]->GetNbinsX()));
  }
}

void MomentumGami::Unfold(TH1F* InputDist) {
  if (InputDist->FindBin(fMaxkStar) > InputDist->GetNbinsX()) {
    Error("MomentumGami::Unfold",
          "Distribution has smaller kStar range than unfolding");
    return;
  }

  fToUnfold = (TH1F*) InputDist->Clone("ToUnfold");
  TF1 * momSmearing = new TF1("momSmearing", this, &MomentumGami::Eval, 0,
                              fMaxkStar, fToUnfold->FindBin(fMaxkStar) - 5,
                              "momSmearing", "momSmearing");  // create TF1 class.
  momSmearing->SetParameter(0, 1.5);
  momSmearing->SetParLimits(0, 0.5, 1.9);
  momSmearing->SetParameter(1, 1.1);
  momSmearing->SetParLimits(1, 0.5, 1.9);
  for (int iPar = 2; iPar < fToUnfold->FindBin(fMaxkStar) - 5; ++iPar) {
    momSmearing->SetParameter(iPar, 1);
    momSmearing->SetParLimits(iPar, 0.9, 1.1);
  }
  fToUnfold->Fit("momSmearing", "QR");

  for (int iBim = 1; iBim < InputDist->FindBin(fMaxkStar) - 4; ++iBim) {
    int ParNmb = iBim - 1;
    InputDist->SetBinContent(
        iBim,
        momSmearing->GetParameter(ParNmb) * InputDist->GetBinContent(iBim));
    //keep the same relative error
    InputDist->SetBinError(
        iBim,
        momSmearing->GetParameter(ParNmb) * InputDist->GetBinError(iBim));
  }
  delete momSmearing;
  delete fToUnfold;
  return;
}

double MomentumGami::Eval(double *x, double *p) {
  if (fToUnfold->FindBin(x[0]) == 0
      || fToUnfold->FindBin(x[0]) > fToUnfold->GetNbinsX()) {
    Error("MomentumGami::operator()",
          "Trying to call with x-Value out of range");
    return -999;
  }
  //corrected "imaginary histo"
  const int nbinsProj = fToUnfold->FindBin(fMaxkStar);
  std::vector<float> zz;
  for (int i = 0; i < nbinsProj - 3; i++) {
    zz.push_back(p[i] * fToUnfold->GetBinContent(i + 1));
  }

//now build uncorrected: take zz "imaginary corrected" and smear it:
  std::vector<float> uncorr;

  for (int iproj = 0; iproj < nbinsProj; iproj++) {
    float binCont = 0;
    for (int ibin = 1; ibin <= nbinsProj - 5; ibin++) {
      binCont += fResProjection[iproj]->GetBinContent(ibin) * zz[ibin - 1];
    }
    uncorr.push_back(binCont);
  }

  return uncorr[fToUnfold->FindBin(x[0]) - 1];
}


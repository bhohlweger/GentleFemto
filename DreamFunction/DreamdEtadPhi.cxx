/*
 * DreamdEtadPhi.cxx
 *
 *  Created on: Feb 4, 2019
 *      Author: schmollweger
 */

#include "DreamdEtadPhi.h"
#include "TMath.h"
#include <iostream>

DreamdEtadPhi::DreamdEtadPhi()
    : fSEdEtadPhi(nullptr),
      fMEdEtadPhi(nullptr),
      fdEtadPhi(nullptr),
      fProjectionY(nullptr) {
  // TODO Auto-generated constructor stub

}

DreamdEtadPhi::~DreamdEtadPhi() {
  // TODO Auto-generated destructor stub
}

void DreamdEtadPhi::ShiftAbovePhi() {
  //Shifts everything above pi*2 - 2 to negative values
  TString SEHistName = Form("%sShifted", fSEdEtadPhi->GetName());

  TH2F* SEReplacement = new TH2F(SEHistName.Data(), SEHistName.Data(),
                                 fSEdEtadPhi->GetNbinsX(),
                                 fSEdEtadPhi->GetXaxis()->GetXmin(),
                                 fSEdEtadPhi->GetXaxis()->GetXmax(),
                                 fSEdEtadPhi->GetYaxis()->GetNbins(),
                                 fSEdEtadPhi->GetYaxis()->GetXmin(),
                                 fSEdEtadPhi->GetYaxis()->GetXmax());

  TString MEHistName = Form("%sShifted", fMEdEtadPhi->GetName());
  TH2F* MEReplacement = new TH2F(MEHistName.Data(), MEHistName.Data(),
                                 fMEdEtadPhi->GetNbinsX(),
                                 fMEdEtadPhi->GetXaxis()->GetXmin(),
                                 fMEdEtadPhi->GetXaxis()->GetXmax(),
                                 fSEdEtadPhi->GetYaxis()->GetNbins(),
                                 fSEdEtadPhi->GetYaxis()->GetXmin(),
                                 fSEdEtadPhi->GetYaxis()->GetXmax());

  for (int iXBin = 1; iXBin <= fSEdEtadPhi->GetNbinsX(); ++iXBin) {
    for (int iYBin = fSEdEtadPhi->GetYaxis()->FindBin(0.);
        iYBin <= fSEdEtadPhi->GetNbinsY(); ++iYBin) {
      float yValNew = fSEdEtadPhi->GetYaxis()->GetBinCenter(iYBin);
      if (yValNew > 2 * TMath::Pi() - 1.9) {
        yValNew += -2 * TMath::Pi();
      }
      int iYBinNew = fSEdEtadPhi->GetYaxis()->FindBin(yValNew);
      SEReplacement->SetBinContent(iXBin, iYBinNew,
                                   fSEdEtadPhi->GetBinContent(iXBin, iYBin));
      SEReplacement->GetBinError(iXBin, iYBinNew,
                                 fSEdEtadPhi->GetBinError(iXBin, iYBin));
      MEReplacement->SetBinContent(iXBin, iYBinNew,
                                   fMEdEtadPhi->GetBinContent(iXBin, iYBin));
      MEReplacement->GetBinError(iXBin, iYBinNew,
                                 fMEdEtadPhi->GetBinError(iXBin, iYBin));
      iYBinNew++;
    }
  }
  delete fSEdEtadPhi;
  fSEdEtadPhi = SEReplacement;
  delete fMEdEtadPhi;
  fMEdEtadPhi = MEReplacement;
  return;
}

void DreamdEtadPhi::DivideSEandME() {
  fdEtadPhi = (TH2F*) fSEdEtadPhi->Clone(
      Form("%sDividedME", fSEdEtadPhi->GetName()));
  std::cout << "GetIntegral SE 2D: " << fdEtadPhi->GetIntegral() << std::endl;
  std::cout << "Integral SE 2D: " << fdEtadPhi->Integral() << std::endl;
  fdEtadPhi->Scale(1. / fdEtadPhi->Integral());
  TH2F* tmpME = (TH2F*) fMEdEtadPhi->Clone("tmp");
  tmpME->Scale(1. / tmpME->Integral());
  fdEtadPhi->Divide(tmpME);
  delete tmpME;
  return;
}

void DreamdEtadPhi::ProjectionY() {
  TString HistName = Form("%sProjection", fSEdEtadPhi->GetName());
  fProjectionY = fSEdEtadPhi->ProjectionY(HistName.Data());
  fProjectionY->Scale(1. / fProjectionY->Integral());
  TH1D* tmpME = fMEdEtadPhi->ProjectionY("TmpME");
  tmpME->Scale(1./tmpME->Integral());
  fProjectionY->Divide(tmpME);
  delete tmpME;
}

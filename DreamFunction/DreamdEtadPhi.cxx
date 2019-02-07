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

void DreamdEtadPhi::DivideSEandME(int rebin) {
  fdEtadPhi = (TH2F*) fSEdEtadPhi->Clone(
      Form("%sDividedME", fSEdEtadPhi->GetName()));
  fdEtadPhi->Rebin2D(rebin,rebin);
  fdEtadPhi->Scale(1. / fdEtadPhi->Integral());
  TH2F* tmpME = (TH2F*) fMEdEtadPhi->Clone("tmp");
  tmpME->Rebin2D(rebin,rebin);
  tmpME->Scale(1. / tmpME->Integral());
  fdEtadPhi->Divide(tmpME);
  fdEtadPhi->SetStats(0);

  delete tmpME;
  return;
}

void DreamdEtadPhi::ProjectionY() {
  TString HistName = Form("%sProjection", fSEdEtadPhi->GetName());
  fProjectionY = fSEdEtadPhi->ProjectionY(HistName.Data());
  fProjectionY->Scale(1. / fProjectionY->Integral());
  TH1D* tmpME = fMEdEtadPhi->ProjectionY("TmpME");
  tmpME->Scale(1. / tmpME->Integral());
  fProjectionY->Divide(tmpME);
  delete tmpME;
}

void DreamdEtadPhi::Draw2D(TPad* p, float Rad) {
  if (fdEtadPhi) {
    p->cd();
    p->SetTopMargin(0.08);
    p->SetBottomMargin(0.1);
    p->SetLeftMargin(0.1);
    p->SetRightMargin(0.01);
    fdEtadPhi->SetStats(0);
    if (Rad != 0) {
      fdEtadPhi->SetTitle(Form("TPC Radius %2.0f cm", Rad));
      fdEtadPhi->SetTitleSize(0.1);
    }
    fdEtadPhi->GetXaxis()->SetTitle("#Delta#eta");
    fdEtadPhi->GetXaxis()->SetTitleSize(0.06);
    fdEtadPhi->GetXaxis()->SetLabelSize(0.06);
    fdEtadPhi->GetXaxis()->SetTitleOffset(0.89);
//    fdEtadPhi->GetXaxis()->SetLabelOffset(0.00);
    fdEtadPhi->GetYaxis()->SetTitle("#Delta#varphi");
    fdEtadPhi->GetYaxis()->SetTitleSize(0.06);
    fdEtadPhi->GetYaxis()->SetLabelSize(0.06);
    fdEtadPhi->GetYaxis()->SetTitleOffset(0.89);
//    fdEtadPhi->GetYaxis()->SetLabelOffset(0.9);
    fdEtadPhi->Draw("COLZ");
  } else {
    std::cout << "No fdEtadPhi for " << Rad << std::endl;
  }
}

void DreamdEtadPhi::DrawProjectionY(TPad* p, float Rad) {
  if (fProjectionY) {
    p->cd();
    p->SetTopMargin(0.08);
    p->SetBottomMargin(0.1);
    p->SetLeftMargin(0.1);
    p->SetRightMargin(0.01);
    fProjectionY->SetStats(0);
    if (Rad != 0) {
      fProjectionY->SetTitle(Form("TPC Radius %2.0f cm", Rad));
      fProjectionY->SetTitleSize(0.1);
    }
    fProjectionY->GetXaxis()->SetTitle("#Delta#varphi");
    fProjectionY->GetXaxis()->SetTitleSize(0.06);
    fProjectionY->GetXaxis()->SetLabelSize(0.06);
    fProjectionY->GetXaxis()->SetTitleOffset(0.89);
    fProjectionY->GetYaxis()->SetTitle("C(#Delta#varphi)");
    fProjectionY->GetYaxis()->SetTitleSize(0.06);
    fProjectionY->GetYaxis()->SetLabelSize(0.06);
    fProjectionY->GetYaxis()->SetTitleOffset(0.89);
    fProjectionY->Draw();
  } else {
    std::cout << "No Projection Y for " << Rad << std::endl;
  }
}

/*
 * DecayQA.cxx
 *
 *  Created on: Feb 13, 2019
 *      Author: schmollweger
 */

#include "DecayQA.h"
#include "TCanvas.h"
#include "TString.h"
DecayQA::DecayQA()
    : fReader(),
      fHairyPlotter(new MakeHistosGreat()),
      fFitter(new ForgivingFitter()),
      fDecayCuts(nullptr),
      fAntiDecayCuts(nullptr),
      fDivCanX(0),
      fDivCanY(0) {
  // TODO Auto-generated constructor stub

}

DecayQA::~DecayQA() {
  // TODO Auto-generated destructor stub
}

void DecayQA::InvariantMassLambda(float CutMin, float CutMax) {
  auto invMassPart = (TH2F*) fReader->Get2DHistInList(
      fReader->GetListInList(fDecayCuts, { "v0Cuts" }), "InvMassPt");
  FitInvariantMass(invMassPart, CutMin, CutMax, "Lambda");
//  auto invMassPart = (TH2F*) fReader->Get2DHistInList(
//      fReader->GetListInList(fDecayCuts, { "v0Cuts" }), "InvMassPt");
//  FitInvariantMass(invMassPart,CutMin,CutMax,"Lambda");
}

void DecayQA::FitInvariantMass(TH2F* invMasspT, float CutMin, float CutMax,
                               const char* outname) {
  //First project the whole thing into one bin
  auto* invMass = (TH1F*) invMasspT->ProjectionY(Form("InvMass%s", outname));
  invMass->GetXaxis()->SetRangeUser(1.107, 1.16);
  fFitter->FitInvariantMass(invMass, CutMin, CutMax);
  fHairyPlotter->FormatHistogram(invMass, 1, 1);
  fHairyPlotter->DrawAndStore( { invMass }, outname);

  auto* cMassBins = new TCanvas(Form("c%s", outname), Form("c%s", outname));
  cMassBins->Divide(fDivCanX, fDivCanY);
  if (invMasspT->GetXaxis()->GetNbins() > fDivCanX * fDivCanY) {
    std::cerr << "FitInvariantMass: Number of divisions not sufficient"
              " to plot all pT bins: \n"
              << "pT Bins: " << invMasspT->GetXaxis()->GetNbins() << '\n';
  }
  for (int ipT = 1; ipT < invMasspT->GetXaxis()->GetNbins(); ++ipT) {
    TPad* CurrentPad = (TPad*) cMassBins->cd(ipT);
    auto invMasspTBin = (TH1F*) invMasspT->ProjectionY(
        Form("%sInvMasspT%u", outname, ipT), ipT, ipT);
    fFitter->FitInvariantMass(invMasspTBin, CutMin, CutMax);
    fHairyPlotter->FormatHistogram(invMasspTBin, 1, 1);
    fHairyPlotter->DrawOnPad( { invMasspTBin }, CurrentPad);
  }
  cMassBins->SaveAs(Form("InvMasspT_%s.pdf", outname));
}

void DecayQA::SetRangesFitting(float signalMin, float signalMax, float bkgMin,
                               float bkgMax) {
  if (fFitter) {
    fFitter->SetRanges(signalMin, signalMax, bkgMin, bkgMax);
  }
}

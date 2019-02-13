/*
 * DecayQA.cxx
 *
 *  Created on: Feb 13, 2019
 *      Author: schmollweger
 */

#include "DecayQA.h"
#include "TCanvas.h"
DecayQA::DecayQA()
    : fReader(),
      fHairyPlotter(new MakeHistosGreat()),
      fFitter(new ForgivingFitter()),
      fDecayCuts(nullptr),
      fAntiDecayCuts(nullptr) {
  // TODO Auto-generated constructor stub

}

DecayQA::~DecayQA() {
  // TODO Auto-generated destructor stub
}

void DecayQA::InvariantMassLambda(float CutMin, float CutMax) {
  auto invMassPart = (TH2F*) fReader->Get2DHistInList(
      fReader->GetListInList(fDecayCuts, { "v0Cuts" }), "InvMassPt");
  invMassPart->GetXaxis()->SetRangeUser(1.107,1.16);
  FitInvariantMass(invMassPart,CutMin,CutMax,"Lambda");
//  auto invMassPart = (TH2F*) fReader->Get2DHistInList(
//      fReader->GetListInList(fDecayCuts, { "v0Cuts" }), "InvMassPt");
//  FitInvariantMass(invMassPart,CutMin,CutMax,"Lambda");
}

void DecayQA::FitInvariantMass(TH2F* invMasspT, float CutMin, float CutMax,
                               const char* outname) {
  //First project the whole thing into one bin
  auto* c = new TCanvas(outname,outname);
  auto* invMass = (TH1F*) invMasspT->ProjectionY(Form("InvMass%s", outname));
  fFitter->FitInvariantMass(invMass, CutMin, CutMax);
//  invMass->Draw();
//  for (int ipT = 0; ipT < 9; ++ipT) {
//  }
  c->cd();
  invMass->Draw();
  TList *funListLambda = invMass->GetListOfFunctions();
  TF1 *fLambdaTotal = (TF1*) funListLambda->FindObject("fLambda");
  fLambdaTotal->Draw("SAME");
  c->SaveAs(Form("%sFit.pdf",outname));
}

void DecayQA::SetRangesFitting(float signalMin, float signalMax, float bkgMin,
                               float bkgMax) {
  if (fFitter) {
    fFitter->SetRanges(signalMin, signalMax, bkgMin, bkgMax);
  }
}

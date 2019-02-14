/*
 * DecayQA.cxx
 *
 *  Created on: Feb 13, 2019
 *      Author: schmollweger
 */

#include "DecayQA.h"
#include "TCanvas.h"
#include "TString.h"
#include "TLine.h"
DecayQA::DecayQA(const char* partLatex)
    : fReader(),
      fHairyPlotter(new MakeHistosGreat()),
      fFitter(new ForgivingFitter()),
      fDecayCuts(nullptr),
      fAntiDecayCuts(nullptr),
      fDivCanX(0),
      fDivCanY(0),
      fPartLatex(partLatex) {
  // TODO Auto-generated constructor stub

}

DecayQA::~DecayQA() {
  // TODO Auto-generated destructor stub
}

void DecayQA::InvariantMassLambda(float CutMin, float CutMax) {
  auto invMassPart = (TH2F*) fReader->Get2DHistInList(
      fReader->GetListInList(fDecayCuts, { "v0Cuts" }), "InvMassPt");
  FitInvariantMass(invMassPart, CutMin, CutMax, "Lambda");
  PlotKaonRejection(
      (TH1F*) fReader->Get1DHistInList(fReader->GetListInList(fDecayCuts, {
                                                                  "v0Cuts" }),
                                       "InvMassKaon"),
      "Lambda");

  auto invMassAntiPart = (TH2F*) fReader->Get2DHistInList(
      fReader->GetListInList(fAntiDecayCuts, { "v0Cuts" }), "InvMassPt");
  FitInvariantMass(invMassAntiPart, CutMin, CutMax, "AntiLambda");
  PlotKaonRejection(
      (TH1F*) fReader->Get1DHistInList(fReader->GetListInList(fAntiDecayCuts, {
                                                                  "v0Cuts" }),
                                       "InvMassKaon"),
      "AntiLambda");
}

void DecayQA::FitInvariantMass(TH2F* invMasspT, float CutMin, float CutMax,
                               const char* outname) {
  //First project the whole thing into one bin
  auto* cMassIntegrated = new TCanvas(Form("cInt%s", outname),
                                      Form("cInt%s", outname));
  auto* invMass = (TH1F*) invMasspT->ProjectionY(Form("InvMass%s", outname), 0,
                                                 -1, "e");
  TPad* intPad = (TPad*) cMassIntegrated->cd();
  fFitter->FitInvariantMass(invMass, CutMin, CutMax);
  invMass->GetXaxis()->SetRangeUser(0.99 * CutMin, 1.01 * CutMax);
  invMass->GetYaxis()->SetRangeUser(0, invMass->GetMaximum() * 1.8);
  invMass->GetYaxis()->SetTitle("d#it{N}/d#it{M} [(GeV/#it{c}^{2})^{-1})]");
  invMass->GetXaxis()->SetTitle("#it{M}_{p#pi} (GeV/#it{c}^{2})");
  fHairyPlotter->FormatHistogram(invMass, 0, 0, 0.8);
  fHairyPlotter->DrawOnPad( { invMass }, intPad, "P");
  fHairyPlotter->DrawLatexLabel(invMasspT->GetXaxis()->GetXmin(),
                                invMasspT->GetXaxis()->GetXmax(), fFitter,
                                intPad, fPartLatex, 0.8, 0.45);
  fHairyPlotter->DrawLine(intPad, 1.112, 1.112, 0,
                          invMass->GetMaximum() * 0.5);
  fHairyPlotter->DrawLine(intPad, 1.120, 1.120, 0,
                          invMass->GetMaximum() * 0.5);
  cMassIntegrated->SaveAs(Form("InvInt%s.pdf", outname));

  auto* cMassBins = new TCanvas(Form("c%s", outname), Form("c%s", outname));
  cMassBins->Divide(fDivCanX, fDivCanY);
  if (invMasspT->GetXaxis()->GetNbins() > fDivCanX * fDivCanY) {
    std::cerr << "FitInvariantMass: Number of divisions not sufficient"
              " to plot all pT bins: \n"
              << "pT Bins: " << invMasspT->GetXaxis()->GetNbins() << '\n';
  }
  auto* Purity = new TH1F(Form("%sPurity", outname), Form("%sPurity", outname),
                          invMasspT->GetXaxis()->GetNbins(),
                          invMasspT->GetXaxis()->GetXmin(),
                          invMasspT->GetXaxis()->GetXmax());
  Purity->GetXaxis()->SetTitle("p_{T} (GeV/#it{c})");
  Purity->GetYaxis()->SetTitle(
      Form("Purity / %.1f (GeV/#it{c})^{-1}", Purity->GetBinWidth(1)));
  for (int ipT = 1; ipT < invMasspT->GetXaxis()->GetNbins(); ++ipT) {
    TPad* CurrentPad = (TPad*) cMassBins->cd(ipT);
    CurrentPad->SetTopMargin(0.08);
    CurrentPad->SetRightMargin(0.03);
    auto invMasspTBin = (TH1F*) invMasspT->ProjectionY(
        Form("%sInvMasspT%u", outname, ipT), ipT, ipT, "e");
    fFitter->FitInvariantMass(invMasspTBin, CutMin, CutMax);
    invMasspTBin->GetXaxis()->SetRangeUser(0.99 * CutMin, 1.01 * CutMax);
    invMasspTBin->GetYaxis()->SetRangeUser(0,
                                           invMasspTBin->GetMaximum() * 1.75);
    invMasspTBin->GetXaxis()->SetTitle("#it{M}_{p#pi} (GeV/#it{c}^{2})");
    invMasspTBin->GetXaxis()->SetMaxDigits(2);
    invMasspTBin->GetXaxis()->SetNdivisions(420);
    invMasspTBin->GetYaxis()->SetMaxDigits(4);
    invMasspTBin->GetYaxis()->SetNdivisions(210);
    fHairyPlotter->FormatSmallHistogram(invMasspTBin, 0, 0, 0.7);
    fHairyPlotter->DrawOnPad( { invMasspTBin }, CurrentPad, "");
    fHairyPlotter->DrawLatexLabel(invMasspT->GetXaxis()->GetBinLowEdge(ipT),
                                  invMasspT->GetXaxis()->GetBinUpEdge(ipT),
                                  fFitter, CurrentPad, fPartLatex, 0.8, 0.35);
    fHairyPlotter->DrawLine(CurrentPad, 1.112, 1.112, 0,
                            invMasspTBin->GetMaximum() * 0.5);
    fHairyPlotter->DrawLine(CurrentPad, 1.120, 1.120, 0,
                            invMasspTBin->GetMaximum() * 0.5);
    float signal = (float) fFitter->GetSignalCounts();
    float background = (float) fFitter->GetBackgroundCounts();
    Purity->SetBinContent(ipT, signal / (signal + background));
    Purity->SetBinError(ipT, 0.03 * signal / (signal + background));
  }
  cMassBins->SaveAs(Form("InvMasspT_%s.pdf", outname));
  fHairyPlotter->FormatHistogram(Purity, 0, 0, 1.0);
  fHairyPlotter->DrawAndStore( { Purity }, Form("Purity%s", outname), "PE1");
}
void DecayQA::PlotKaonRejection(TH1F* invMassKaon, const char* outname) {
  invMassKaon->SetName(Form("%s%s", outname, invMassKaon->GetName()));
  invMassKaon->GetXaxis()->SetRangeUser(0.44, 0.56);
  invMassKaon->GetYaxis()->SetRangeUser(0, 1.8 * invMassKaon->GetMaximum());
  invMassKaon->GetYaxis()->SetTitle("d#it{N}/d#it{M} [(GeV/#it{c}^{2})^{-1})]");
  invMassKaon->GetXaxis()->SetTitle("#it{M}_{#pi#pi} (GeV/#it{c}^{2})");
  ForgivingFitter *kaonFit = new ForgivingFitter();
  kaonFit->SetRanges(0.485, 0.515, 0.44, 0.56);
  kaonFit->FitInvariantMass(invMassKaon, 0.49, 0.51);
  fHairyPlotter->FormatHistogram(invMassKaon, 0, 0, 1.0);
  auto* canKaon = new TCanvas(Form("CanKaon%s", outname),
                              Form("CanKaon%s", outname));
  TPad* tmppad = (TPad*) canKaon->cd();
  fHairyPlotter->DrawOnPad( { invMassKaon }, tmppad, "PE");
  fHairyPlotter->DrawLatexLabel(0.3, 4.3, kaonFit, tmppad, "K^{0}", 0.8, 0.45);
  fHairyPlotter->DrawLine(tmppad, 0.48, 0.48, 0,
                          invMassKaon->GetMaximum() * 0.5);
  fHairyPlotter->DrawLine(tmppad, 0.515, 0.515, 0,
                          invMassKaon->GetMaximum() * 0.5);
  canKaon->SaveAs(Form("%sKaon.pdf", outname));
  delete kaonFit;
  delete canKaon;
}
void DecayQA::SetRangesFitting(float signalMin, float signalMax, float bkgMin,
                               float bkgMax) {
  if (fFitter) {
    fFitter->SetRanges(signalMin, signalMax, bkgMin, bkgMax);
  }
}

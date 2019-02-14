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
DecayQA::DecayQA(const char* partLatex, const char* latexProducts)
    : fReader(),
      fHairyPlotter(new MakeHistosGreat()),
      fFitter(new ForgivingFitter()),
      fDecayCuts(nullptr),
      fAntiDecayCuts(nullptr),
      fDivCanX(0),
      fDivCanY(0),
      fInvMassPtStartBin(1),
      fPartLatex(partLatex),
      fScaleMax(0),
      fTexOffX(0),
      fTexOffY(0),
      fDecChannel(latexProducts) {
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

void DecayQA::InvariantMassXi(float CutMin, float CutMax) {
  auto invMassPart = (TH2F*) fReader->Get2DHistInList(
      fReader->GetListInList(fDecayCuts, { "Cascade" }), "InvMassXi");
  FitInvariantMass(invMassPart, CutMin, CutMax, "Xi");

  auto invMassAntiPart = (TH2F*) fReader->Get2DHistInList(
      fReader->GetListInList(fAntiDecayCuts, { "Cascade" }), "InvMassXi");
  FitInvariantMass(invMassAntiPart, CutMin, CutMax, "AntiXi");
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
  invMass->GetYaxis()->SetMaxDigits(3);
  invMass->GetYaxis()->SetTitle("d#it{N}/d#it{M} [(GeV/#it{c}^{2})^{-1})]");
  invMass->GetXaxis()->SetTitle(Form("#it{M}_{%s} (GeV/#it{c}^{2})",fDecChannel));
  fHairyPlotter->FormatHistogram(invMass, 0, 0, 0.8);
  fHairyPlotter->DrawOnPad( { invMass }, intPad, "P");
  fHairyPlotter->DrawLatexLabel(
      invMasspT->GetXaxis()->GetBinLowEdge(fInvMassPtStartBin),
      invMasspT->GetXaxis()->GetXmax(), fFitter, intPad, fPartLatex, 0.8, 0.45);
  fHairyPlotter->DrawLine(intPad, CutMin, CutMin, 0,
                          invMass->GetMaximum() * 0.5);
  fHairyPlotter->DrawLine(intPad, CutMax, CutMax, 0,
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
  for (int ipT = fInvMassPtStartBin;
      ipT < invMasspT->GetXaxis()->GetNbins() + 1; ++ipT) {
    unsigned int iPad = 0;
    if (fInvMassPtStartBin != 1) {
      iPad = ipT - fInvMassPtStartBin + 1;
    } else {
      iPad = ipT;
    }
    TPad* CurrentPad = (TPad*) cMassBins->cd(iPad);
    CurrentPad->SetTopMargin(0.08);
    CurrentPad->SetRightMargin(0.03);
    auto invMasspTBin = (TH1F*) invMasspT->ProjectionY(
        Form("%sInvMasspT%u", outname, ipT), ipT, ipT, "e");
    fFitter->FitInvariantMass(invMasspTBin, CutMin, CutMax);
    invMasspTBin->GetXaxis()->SetRangeUser(0.99 * CutMin, 1.01 * CutMax);
    invMasspTBin->GetYaxis()->SetRangeUser(
        0, invMasspTBin->GetMaximum() * fScaleMax);
    invMasspTBin->GetXaxis()->SetTitle(Form("#it{M}_{%s} (GeV/#it{c}^{2})",fDecChannel));
    invMasspTBin->GetXaxis()->SetMaxDigits(2);
    invMasspTBin->GetXaxis()->SetNdivisions(420);
    invMasspTBin->GetYaxis()->SetMaxDigits(3);
    invMasspTBin->GetYaxis()->SetNdivisions(210);
    fHairyPlotter->FormatSmallHistogram(invMasspTBin, 0, 0, 0.7);
    fHairyPlotter->DrawOnPad( { invMasspTBin }, CurrentPad, "");
    fHairyPlotter->DrawLatexLabel(invMasspT->GetXaxis()->GetBinLowEdge(ipT),
                                  invMasspT->GetXaxis()->GetBinUpEdge(ipT),
                                  fFitter, CurrentPad, fPartLatex, fTexOffX,
                                  fTexOffY);
    fHairyPlotter->DrawLine(CurrentPad, CutMin, CutMin, 0,
                            invMasspTBin->GetMaximum() * 0.5);
    fHairyPlotter->DrawLine(CurrentPad, CutMax, CutMax, 0,
                            invMasspTBin->GetMaximum() * 0.5);
    float signal = (float) fFitter->GetSignalCounts();
    float background = (float) fFitter->GetBackgroundCounts();
    Purity->SetBinContent(ipT, signal / (signal + background));
//    Purity->SetBinError(ipT, 0.03 * signal / (signal + background));
  }
  Purity->GetYaxis()->SetRangeUser(0.7, 1.1);
  cMassBins->SaveAs(Form("InvMasspT_%s.pdf", outname));
  fHairyPlotter->FormatHistogram(Purity, 0, 0, 1.5);
  fHairyPlotter->DrawAndStore( { Purity }, Form("Purity%s", outname), "P");
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

void DecayQA::PlotQATopologyLambda() {
  PlotQATopologyLambda(fDecayCuts, "Lambda");
  PlotQATopologyLambda(fAntiDecayCuts, "AntiLambda");
}
void DecayQA::PlotQATopologyLambda(TList *v0Cuts, const char* outname) {
  auto dcaDaugVtx = (TH1F*) fReader->Get1DHistInList(
      fReader->GetListInList(v0Cuts, { "v0Cuts", "after" }),
      "DCADauToVtx_after");
  if (!dcaDaugVtx) {
    std::cerr << "dcaDaugter To Vtx is missing for " << outname << std::endl;
  }
  dcaDaugVtx->GetXaxis()->SetRangeUser(0, 3.);
  dcaDaugVtx->GetXaxis()->SetTitle("DCA(p,#pi) (cm)");
  dcaDaugVtx->GetYaxis()->SetTitle("Entries");
  fHairyPlotter->FormatHistogram(dcaDaugVtx, 0, 1);
  fHairyPlotter->DrawLogYAndStore( { dcaDaugVtx },
                                  Form("%sdcaDaugDecVtx", outname));

  auto v0TRad = (TH1F*) fReader->Get1DHistInList(
      fReader->GetListInList(v0Cuts, { "v0Cuts", "after" }),
      "TransverseRadius_after");
  if (!v0TRad) {
    std::cerr << "Transverse Radius Distribution is missing for " << outname
              << std::endl;
  }
  v0TRad->GetXaxis()->SetRangeUser(0, 120.);
  v0TRad->GetXaxis()->SetTitle("r_{xy} (cm)");
  v0TRad->GetYaxis()->SetTitle("Entries");
  fHairyPlotter->FormatHistogram(v0TRad, 0, 1);
  fHairyPlotter->DrawAndStore( { v0TRad }, Form("%sTransRad", outname));

  auto DCAPos = (TH1F*) fReader->Get1DHistInList(
      fReader->GetListInList(v0Cuts, { "v0Cuts", "after" }),
      "DCADauPToPV_after");
  if (!DCAPos) {
    std::cerr << "DCA Pos Daug Distribution is missing for " << outname
              << std::endl;
  }
  DCAPos->GetXaxis()->SetRangeUser(0, 100);
  DCAPos->GetXaxis()->SetTitle("DCA(Pos Daug,PV) (cm)");
  DCAPos->GetYaxis()->SetTitle("Entries");
  fHairyPlotter->FormatHistogram(DCAPos, 0, 1);
  fHairyPlotter->DrawLogYAndStore( { DCAPos }, Form("%sDCAPVPosDaug", outname));

  auto DCANeg = (TH1F*) fReader->Get1DHistInList(
      fReader->GetListInList(v0Cuts, { "v0Cuts", "after" }),
      "DCADauNToPV_after");
  if (!DCANeg) {
    std::cerr << "DCA Pos Daug Distribution is missing for " << outname
              << std::endl;
  }
  DCANeg->GetXaxis()->SetRangeUser(0, 100.);
  DCANeg->GetXaxis()->SetTitle("DCA(Neg Daug,PV) (cm)");
  DCANeg->GetYaxis()->SetTitle("Entries");
  fHairyPlotter->FormatHistogram(DCANeg, 0, 1);
  fHairyPlotter->DrawLogYAndStore( { DCANeg }, Form("%sDCAPVNegDaug", outname));

  auto* pTDist = (TH1F*) fReader->Get1DHistInList(
      fReader->GetListInList(v0Cuts, { "v0Cuts", "after" }), "pTDist_after");
  if (!pTDist) {
    std::cerr << "pT Distribution missing for " << outname << std::endl;
  }
  pTDist->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  pTDist->GetYaxis()->SetTitle(
      Form("Entries/ %.1f GeV/#it{c}", pTDist->GetBinWidth(1)));
  fHairyPlotter->FormatHistogram(pTDist, 0, 1);
  fHairyPlotter->DrawLogYAndStore( { pTDist }, Form("%spT", outname), "PE");

  auto* phiDist = (TH1F*) fReader->Get1DHistInList(
      fReader->GetListInList(v0Cuts, { "v0Cuts", "after" }), "PhiDist_after");
  if (!phiDist) {
    std::cerr << "phi Distribution missing for " << outname << std::endl;
  }
  phiDist->GetXaxis()->SetTitle("#varphi (rad)");
  phiDist->GetYaxis()->SetTitle(
      Form("Entries/ %.3f rad", phiDist->GetBinWidth(1)));
  fHairyPlotter->FormatHistogram(phiDist, 0, 1);
  fHairyPlotter->DrawAndStore( { phiDist }, Form("%s_phi", outname));

  auto* etaDist = (TH1F*) fReader->Get1DHistInList(
      fReader->GetListInList(v0Cuts, { "v0Cuts", "after" }), "EtaDist_after");
  if (!etaDist) {
    std::cerr << "eta Distribution missing for " << outname << std::endl;
  }
  etaDist->GetXaxis()->SetTitle("#eta");
  etaDist->GetXaxis()->SetRangeUser(-1., 1.);
  etaDist->GetYaxis()->SetTitle(
      Form("Entries/ %.2f ", etaDist->GetBinWidth(1)));
  fHairyPlotter->FormatHistogram(etaDist, 0, 1);
  fHairyPlotter->DrawAndStore( { etaDist }, Form("%s_eta", outname), "hist");

}

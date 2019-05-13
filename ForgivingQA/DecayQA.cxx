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

void DecayQA::InvariantMassLambdaSigma0(float CutMin, float CutMax) {
  auto invMassPart = (TH2F*) fReader->Get2DHistInList(fDecayCuts, "InvMassPt");
  invMassPart->RebinX(10);
  FitInvariantMass(invMassPart, CutMin, CutMax, "Lambda");
  PlotKaonRejection(
      (TH1F*) fReader->Get1DHistInList(fDecayCuts, "fHistK0MassAfter"),
      "Lambda");

  auto invMassAntiPart = (TH2F*) fReader->Get2DHistInList(fAntiDecayCuts,
                                                          "InvMassPt");
  invMassAntiPart->RebinX(10);
  FitInvariantMass(invMassAntiPart, CutMin, CutMax, "AntiLambda");
  PlotKaonRejection(
      (TH1F*) fReader->Get1DHistInList(fAntiDecayCuts, "fHistK0MassAfter"),
      "AntiLambda");
}

void DecayQA::InvariantMassSigma0(float massCuts, const char* name, bool isSum) {
  TH2F* invMassPart;
  if (!isSum) {
    if (fDecayCuts) {
      invMassPart = (TH2F*) fReader->Get2DHistInList(fDecayCuts,
                                                     "fHistInvMassPtRaw");
    } else if (fAntiDecayCuts) {
      invMassPart = (TH2F*) fReader->Get2DHistInList(fAntiDecayCuts,
                                                     "fHistInvMassPtRaw");

    } else {
      std::cout << "No cuts set\n";
      return;
    }
  } else {
    invMassPart = (TH2F*) fReader->Get2DHistInList(fDecayCuts,
                                                   "fHistInvMassPtRaw");
    auto invMassAntiPart = (TH2F*) fReader->Get2DHistInList(
        fAntiDecayCuts, "fHistInvMassPtRaw");
    invMassPart->Add(invMassAntiPart);
  }
  invMassPart->RebinX(20);
  FitInvariantMassSigma0(invMassPart, massCuts, name);
}

void DecayQA::GetPeriodQA(float CutMin, float CutMax,
                          std::vector<const char*> pathToList,
                          const char* histname) {
  auto invMassPart = (TH2F*) fReader->Get2DHistInList(
      fReader->GetListInList(fDecayCuts, pathToList), histname);
  auto* invMass = (TH1F*) invMassPart->ProjectionY("InvMassClone", 0, -1, "e");
  fFitter->FitInvariantMass(invMass, CutMin, CutMax);
}

void DecayQA::GetPeriodQASigma(float CutMin, float CutMax, const char* period) {
  auto invMasspT = (TH2F*) fReader->Get2DHistInList(fDecayCuts, "InvMassPt");
  invMasspT->RebinX(10);
  auto* invMass = (TH1F*) invMasspT->ProjectionY("InvMassClone", 0, -1, "e");
  fFitter->FitInvariantMass(invMass, CutMin, CutMax);
}

void DecayQA::GetPeriodQASigma0(float massCuts, const char* period) {
  auto invMasspT = (TH2F*) fReader->Get2DHistInList(fDecayCuts,
                                                    "fHistInvMassPtRaw");
  auto invMassAntiPart = (TH2F*) fReader->Get2DHistInList(fAntiDecayCuts,
                                                          "fHistInvMassPtRaw");
  invMasspT->Add(invMassAntiPart);
  auto invMass = (TH1F*) invMasspT->ProjectionY("InvMassClone", 0, -1, "e");
  fFitter->FitInvariantMassSigma(invMass, massCuts);
}

void DecayQA::InvariantMassXi(float CutMin, float CutMax) {
  auto invMassPart = (TH2F*) fReader->Get2DHistInList(
      fReader->GetListInList(fDecayCuts, { "Cascade" }), "InvMassXi");
  FitInvariantMass(invMassPart, CutMin, CutMax, "Xi");

  auto invMassAntiPart = (TH2F*) fReader->Get2DHistInList(
      fReader->GetListInList(fAntiDecayCuts, { "Cascade" }), "InvMassXi");
  FitInvariantMass(invMassAntiPart, CutMin, CutMax, "AntiXi");
}

void DecayQA::IvariantMassXiLambda() {
  auto invMassXiLa = (TH2F*) fReader->Get2DHistInList(
      fReader->GetListInList(fDecayCuts, { "Cascade" }), "InvMassv0Pt");
  auto invMassLa = (TH1F*) invMassXiLa->ProjectionY("invMassXiLas", 0, -1, "e");
  fHairyPlotter->FormatHistogram(invMassLa, 0, 0);
  auto cInvLa = new TCanvas("invMassXiLa", "invMassXiLa");
  fFitter->ShittyInvariantMass(invMassLa, cInvLa, 0.2, 6.3, "#Lambda");
  cInvLa->SaveAs("InvariantMassXiLa.pdf");

  auto invMassAntiXiLa = (TH2F*) fReader->Get2DHistInList(
      fReader->GetListInList(fAntiDecayCuts, { "Cascade" }), "InvMassv0Pt");
  auto invMassAntiLa = (TH1F*) invMassAntiXiLa->ProjectionY("invMassXiALas", 0,
                                                            -1, "e");
  fHairyPlotter->FormatHistogram(invMassAntiLa, 0, 0);
  auto cInvALa = new TCanvas("invMassAXiALa", "invMassAXiALa");
  fFitter->ShittyInvariantMass(invMassAntiLa, cInvALa, 0.2, 6.3, "#Lambda");
  cInvALa->SaveAs("InvariantMassAXiALa.pdf");

  auto cInvLapT = new TCanvas("InvMassXiLaPT", "InvMassXiLaPT", 0, 0, 2000,
                              1500);
  cInvLapT->Divide(4, 3);
  auto cInvALapT = new TCanvas("InvMassAXiALaPT", "InvMassAXiALaPT", 0, 0, 2000,
                               1500);
  cInvALapT->Divide(4, 3);
  for (int iPt = 1; iPt <= invMassXiLa->GetXaxis()->GetNbins(); ++iPt) {
    TPad* tmppadLa = (TPad*) cInvLapT->cd(iPt);
    auto invMassLapT = (TH1F*) invMassXiLa->ProjectionY(
        Form("invMassLaPt_%u", iPt + 1), iPt + 1, iPt + 1, "e");
    fHairyPlotter->FormatSmallHistogram(invMassLapT, 0, 0);
    fFitter->ShittyInvariantMass(
        invMassLapT, tmppadLa, invMassXiLa->GetXaxis()->GetBinLowEdge(iPt + 1),
        invMassXiLa->GetXaxis()->GetBinUpEdge(iPt + 1), "#Lambda");

    TPad* tmppadALa = (TPad*) cInvALapT->cd(iPt);
    auto invMassALapT = (TH1F*) invMassAntiXiLa->ProjectionY(
        Form("invMassALaPt_%u", iPt + 1), iPt + 1, iPt + 1, "e");
    fHairyPlotter->FormatSmallHistogram(invMassALapT, 0, 0);
    fFitter->ShittyInvariantMass(
        invMassALapT, tmppadALa,
        invMassAntiXiLa->GetXaxis()->GetBinLowEdge(iPt + 1),
        invMassAntiXiLa->GetXaxis()->GetBinUpEdge(iPt + 1), "#Lambda");
  }
  cInvLapT->SaveAs("InvMassXiLaPt.pdf");
  cInvALapT->SaveAs("InvMassAXiALaPt.pdf");
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
  //invMass->GetYaxis()->SetMaxDigits(3);
  invMass->GetYaxis()->SetTitle("d#it{N}/d#it{M} [(GeV/#it{c}^{2})^{-1})]");
  invMass->GetXaxis()->SetTitle(
      Form("#it{M}_{%s} (GeV/#it{c}^{2})", fDecChannel));
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
    invMasspTBin->GetXaxis()->SetTitle(
        Form("#it{M}_{%s} (GeV/#it{c}^{2})", fDecChannel));
    //invMasspTBin->GetXaxis()->SetMaxDigits(2);
    invMasspTBin->GetXaxis()->SetNdivisions(505);
    //invMasspTBin->GetYaxis()->SetMaxDigits(3);
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

void DecayQA::FitInvariantMassSigma0(TH2F* invMasspT, float massCuts,
                                     const char* outname) {
  //First project the whole thing into one bin
  auto* cMassIntegrated = new TCanvas(Form("cInt%s", outname),
                                      Form("cInt%s", outname), 0, 0, 650, 550);
  cMassIntegrated->SetTopMargin(0.05);
  cMassIntegrated->SetRightMargin(0.025);
  auto* invMass = (TH1F*) invMasspT->ProjectionY(Form("InvMass%s", outname), 0,
                                                 -1, "e");
  TPad* intPad = (TPad*) cMassIntegrated->cd();
  fFitter->FitInvariantMassSigma(invMass, massCuts, kRed + 1);
  const double CutMin = fFitter->GetMeanMass() - massCuts;
  const double CutMax = fFitter->GetMeanMass() + massCuts;
  invMass->GetXaxis()->SetRangeUser(1.172, 1.212);
  invMass->GetYaxis()->SetRangeUser(0, invMass->GetMaximum() * 1.8);
  invMass->GetXaxis()->SetMaxDigits(2);
  invMass->GetXaxis()->SetNdivisions(505);
  invMass->GetYaxis()->SetMaxDigits(1);
  invMass->GetYaxis()->SetNdivisions(505);
  invMass->GetYaxis()->SetTitle("d#it{N}/d#it{M} [(GeV/#it{c}^{2})^{-1})]");
  invMass->GetXaxis()->SetTitle(
      Form("#it{M}_{%s} (GeV/#it{c}^{2})", fDecChannel));
  fHairyPlotter->FormatHistogram(invMass, 2, 8, 1.1);
  fHairyPlotter->DrawOnPad( { invMass }, intPad, "P");
  fHairyPlotter->DrawPerformance(fFitter, intPad, fPartLatex, 0.205, 0.87, 0.4, 10);
  fHairyPlotter->DrawLine(intPad, CutMin, CutMin, 0,
                          invMass->GetMaximum() * 0.5, kTeal + 3);
  fHairyPlotter->DrawLine(intPad, CutMax, CutMax, 0,
                          invMass->GetMaximum() * 0.5, kTeal + 3);
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
    fFitter->FitInvariantMassSigma(invMasspTBin, massCuts);
    invMasspTBin->GetXaxis()->SetRangeUser(0.99 * CutMin, 1.01 * CutMax);
    invMasspTBin->GetYaxis()->SetRangeUser(
        0, invMasspTBin->GetMaximum() * fScaleMax);
    invMasspTBin->GetXaxis()->SetTitle(
        Form("#it{M}_{%s} (GeV/#it{c}^{2})", fDecChannel));
    invMasspTBin->GetXaxis()->SetNdivisions(420);
    invMasspTBin->GetXaxis()->SetMaxDigits(2);
    invMasspTBin->GetYaxis()->SetMaxDigits(2);
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
    Purity->SetBinContent(ipT, fFitter->GetPurity());
    Purity->SetBinError(ipT, fFitter->GetPurityErr());
  }
  Purity->GetYaxis()->SetRangeUser(0.1, 0.7);
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

void DecayQA::PlotQATopologySigma0Daughter(TList* v0Cuts, const char* outname) {
  // DCA daughters at the decay vertex
  auto dcaDaugVtx = (TH1F*) (fReader->Get2DHistInList(
      fReader->GetListInList(v0Cuts, { "After" }), "fHistDCADaughtersAfter"))
      ->ProjectionY();
  if (!dcaDaugVtx) {
    std::cerr << "dcaDaugter To Vtx is missing for " << outname << std::endl;
  }
  dcaDaugVtx->GetXaxis()->SetRangeUser(0, 1.5);
  dcaDaugVtx->GetXaxis()->SetTitle("DCA(daughters at decay vtx) (cm)");
  dcaDaugVtx->GetYaxis()->SetTitle("Entries");
  fHairyPlotter->FormatHistogram(dcaDaugVtx, 0, 1);
  fHairyPlotter->DrawAndStore( { dcaDaugVtx },
                              Form("%sdcaDaugDecVtx", outname));

  // Transverse radius
  auto v0TRad =
      (TH1F*) (fReader->Get2DHistInList(
          fReader->GetListInList(v0Cuts, { "After" }),
          "fHistTransverseRadiusAfter"))->ProjectionY();
  if (!v0TRad) {
    std::cerr << "Transverse Radius Distribution is missing for " << outname
              << std::endl;
  }
  v0TRad->GetXaxis()->SetRangeUser(0, 120.);
  v0TRad->GetXaxis()->SetTitle("r_{xy} (cm)");
  v0TRad->GetYaxis()->SetTitle("Entries");
  fHairyPlotter->FormatHistogram(v0TRad, 0, 1);
  fHairyPlotter->DrawAndStore( { v0TRad }, Form("%sTransRad", outname));

  // DCA pos daughter to PV
  auto DCAPos = (TH1F*) (fReader->Get2DHistInList(
      fReader->GetListInList(v0Cuts, { "V0_PosDaughter" }),
      "fHistSingleParticleDCAtoPVAfter_pos"))->ProjectionY();
  if (!DCAPos) {
    std::cerr << "DCA Pos Daug Distribution is missing for " << outname
              << std::endl;
  }
  DCAPos->GetXaxis()->SetRangeUser(0, 7.5);
  DCAPos->GetXaxis()->SetTitle("DCA(Pos Daug,PV) (cm)");
  DCAPos->GetYaxis()->SetTitle("Entries");
  fHairyPlotter->FormatHistogram(DCAPos, 0, 1);
  fHairyPlotter->DrawLogYAndStore( { DCAPos }, Form("%sDCAPVPosDaug", outname));

  // DCA neg daughter to PV
  auto DCANeg = (TH1F*) (fReader->Get2DHistInList(
      fReader->GetListInList(v0Cuts, { "V0_NegDaughter" }),
      "fHistSingleParticleDCAtoPVAfter_neg"))->ProjectionY();
  if (!DCANeg) {
    std::cerr << "DCA Pos Daug Distribution is missing for " << outname
              << std::endl;
  }
  DCANeg->GetXaxis()->SetRangeUser(0, 7.5);
  DCANeg->GetXaxis()->SetTitle("DCA(Neg Daug,PV) (cm)");
  DCANeg->GetYaxis()->SetTitle("Entries");
  fHairyPlotter->FormatHistogram(DCANeg, 0, 1);
  fHairyPlotter->DrawLogYAndStore( { DCANeg }, Form("%sDCAPVNegDaug", outname));

  // pT
  auto pTDist = (TH1F*) (fReader->Get2DHistInList(v0Cuts, "InvMassPt"))
      ->ProjectionX();
  if (!pTDist) {
    std::cerr << "pT Distribution missing for " << outname << std::endl;
  }
  pTDist->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  pTDist->GetYaxis()->SetTitle(
      Form("Entries/ %.1f GeV/#it{c}", pTDist->GetBinWidth(1)));
  fHairyPlotter->FormatHistogram(pTDist, 0, 1);
  fHairyPlotter->DrawLogYAndStore( { pTDist }, Form("%spT", outname), "PE");

  // Phi
  auto* phiDist = (TH1F*) (fReader->Get2DHistInList(v0Cuts, "fHistEtaPhi"))
      ->ProjectionY();
  if (!phiDist) {
    std::cerr << "phi Distribution missing for " << outname << std::endl;
  }
  phiDist->SetMinimum(0);
  phiDist->GetXaxis()->SetTitle("#varphi (rad)");
  phiDist->GetYaxis()->SetTitle(
      Form("Entries/ %.3f rad", phiDist->GetBinWidth(1)));
  fHairyPlotter->FormatHistogram(phiDist, 0, 1);
  fHairyPlotter->DrawAndStore( { phiDist }, Form("%s_phi", outname));

  // Eta
  auto* etaDist = (TH1F*) (fReader->Get2DHistInList(v0Cuts, "fHistEtaPhi"))
      ->ProjectionX();
  if (!etaDist) {
    std::cerr << "eta Distribution missing for " << outname << std::endl;
  }
  etaDist->GetXaxis()->SetTitle("#eta");
  etaDist->SetMinimum(0);
  etaDist->GetXaxis()->SetRangeUser(-1., 1.);
  etaDist->GetYaxis()->SetTitle(
      Form("Entries/ %.2f ", etaDist->GetBinWidth(1)));
  fHairyPlotter->FormatHistogram(etaDist, 0, 1);
  fHairyPlotter->DrawAndStore( { etaDist }, Form("%s_eta", outname), "hist");

  // Cosine pointing angle
  auto* cpa = (TH1F*) (fReader->Get2DHistInList(
      fReader->GetListInList(v0Cuts, { "After" }), "fHistCosPAAfter"))
      ->ProjectionY();
  if (!cpa) {
    std::cerr << "CPA is missing for " << outname << std::endl;
  }
  cpa->GetXaxis()->SetTitle("cos(#alpha)");
  cpa->GetXaxis()->SetRangeUser(0.99, 1.);
  cpa->GetYaxis()->SetTitle("Entries");
  fHairyPlotter->FormatHistogram(cpa, 0, 1);
  fHairyPlotter->DrawAndStore( { cpa }, Form("%sCPA", outname));
}

void DecayQA::PlotPIDSigma0Daughter(TList* v0Cuts, const char* outname) {
  // Armenteros
  TH2F* armenteros = (TH2F*) fReader->Get2DHistInList(
      fReader->GetListInList(v0Cuts, { "After" }), "fHistArmenterosAfter");
  if (!armenteros) {
    std::cerr << "Armenteros is missing for " << outname << std::endl;
  }
  armenteros->SetTitle("Armenteros-Podolandski");
  fHairyPlotter->FormatHistogram(armenteros, 0, 1);
  std::vector<TH2*> drawVecTPC = { armenteros };
  fHairyPlotter->DrawLogZAndStore(drawVecTPC, Form("%s_ArmenterosPID", outname),
                                  "colz");

  // dE/dc pos daughter
  TH2F* posPID = (TH2F*) fReader->Get2DHistInList(
      fReader->GetListInList(v0Cuts, { "V0_PosDaughter" }),
      "fHistSingleParticlePID_pos");
  if (!posPID) {
    std::cerr << "Pos PID is missing for " << outname << std::endl;
  }
  fHairyPlotter->FormatHistogram(posPID, 0, 1);
  drawVecTPC = {posPID};
  fHairyPlotter->DrawLogZAndStore(drawVecTPC, Form("%s_PosPID", outname),
                                  "colz");

  // Armenteros
  TH2F* negPID = (TH2F*) fReader->Get2DHistInList(
      fReader->GetListInList(v0Cuts, { "V0_NegDaughter" }),
      "fHistSingleParticlePID_neg");
  if (!negPID) {
    std::cerr << "Neg PID is missing for " << outname << std::endl;
  }
  fHairyPlotter->FormatHistogram(negPID, 0, 1);
  drawVecTPC = {negPID};
  fHairyPlotter->DrawLogZAndStore(drawVecTPC, Form("%s_NegPID", outname),
                                  "colz");
}

void DecayQA::PlotQATopologySigma0(TList* v0Cuts, const char* outname) {

  // pT
  auto pTDist = (TH1F*) (fReader->Get2DHistInList(v0Cuts, "fHistInvMassPt"))
      ->ProjectionX();
  if (!pTDist) {
    std::cerr << "pT Distribution missing for " << outname << std::endl;
  }
  pTDist->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  pTDist->GetYaxis()->SetTitle(
      Form("Entries/ %.1f GeV/#it{c}", pTDist->GetBinWidth(1)));
  fHairyPlotter->FormatHistogram(pTDist, 0, 1);
  fHairyPlotter->DrawLogYAndStore( { pTDist }, Form("%spT", outname), "PE");

  // Phi
  auto* phiDist = (TH1F*) (fReader->Get2DHistInList(v0Cuts, "fHistEtaPhi"))
      ->ProjectionY();
  if (!phiDist) {
    std::cerr << "phi Distribution missing for " << outname << std::endl;
  }
  phiDist->GetXaxis()->SetTitle("#varphi (rad)");
  phiDist->SetMinimum(0);
  phiDist->GetYaxis()->SetTitle(
      Form("Entries/ %.3f rad", phiDist->GetBinWidth(1)));
  fHairyPlotter->FormatHistogram(phiDist, 0, 1);
  fHairyPlotter->DrawAndStore( { phiDist }, Form("%s_phi", outname));

  // Eta
  auto* etaDist = (TH1F*) (fReader->Get2DHistInList(v0Cuts, "fHistEtaPhi"))
      ->ProjectionX();
  if (!etaDist) {
    std::cerr << "eta Distribution missing for " << outname << std::endl;
  }
  etaDist->SetMinimum(0);
  etaDist->GetXaxis()->SetTitle("#eta");
  etaDist->GetXaxis()->SetRangeUser(-1., 1.);
  etaDist->GetYaxis()->SetTitle(
      Form("Entries/ %.2f ", etaDist->GetBinWidth(1)));
  fHairyPlotter->FormatHistogram(etaDist, 0, 1);
  fHairyPlotter->DrawAndStore( { etaDist }, Form("%s_eta", outname), "hist");

  // Number of Sigma Candidates
  auto* nSigma = (fReader->Get1DHistInList(v0Cuts, "fHistNSigma"));
  if (!nSigma) {
    std::cerr << "Number of Sigma candidates distribution missing for "
              << outname << std::endl;
  }
  nSigma->GetXaxis()->SetTitle(Form("Number of %s candidates", fPartLatex));
  nSigma->GetYaxis()->SetTitle("Entries");
  fHairyPlotter->FormatHistogram(nSigma, 0, 1);
  fHairyPlotter->DrawLogYAndStore( { nSigma }, Form("%s_nSigma", outname),
                                  "hist");

  // Track cleaner Lambda
  auto* nLambdaLabel = (fReader->Get1DHistInList(v0Cuts, "fHistNLambdaLabel"));
  if (!nLambdaLabel) {
    std::cerr << "Number of Lambda candidates cleanup distribution missing for "
              << outname << std::endl;
  }
  nLambdaLabel->GetYaxis()->SetTitle("Entries");
  fHairyPlotter->FormatHistogram(nLambdaLabel, 0, 1);
  fHairyPlotter->DrawLogYAndStore( { nLambdaLabel },
                                  Form("%s_nLambdaLabel", outname), "hist");

  // Track cleaner photon
  auto* nPhotonLabel = (fReader->Get1DHistInList(v0Cuts, "fHistNPhotonLabel"));
  if (!nPhotonLabel) {
    std::cerr << "Number of photon candidates cleanup distribution missing for "
              << outname << std::endl;
  }
  nPhotonLabel->GetYaxis()->SetTitle("Entries");
  fHairyPlotter->FormatHistogram(nPhotonLabel, 0, 1);
  fHairyPlotter->DrawLogYAndStore( { nPhotonLabel },
                                  Form("%s_nPhotonLabel", outname), "hist");
}

/*
 * TrackQA.cxx
 *
 *  Created on: Feb 8, 2019
 *      Author: schmollweger
 */

#include "TrackQA.h"
#include <stdlib.h>
#include <iostream>

TrackQA::TrackQA()
    : fReader(),
      fHairyPlotter(new MakeHistosGreat()),
      fTrackCuts(nullptr),
      fAntiTrackCuts(nullptr) {
}

TrackQA::TrackQA(const char* outname)
    : fReader(),
      fHairyPlotter(new MakeHistosGreat(outname)),
      fTrackCuts(nullptr),
      fAntiTrackCuts(nullptr) {
}


TrackQA::~TrackQA() {
  delete fHairyPlotter;
}

void TrackQA::PlotKinematic() {
  PlotKinematic(fTrackCuts, "Proton");
  PlotKinematic(fAntiTrackCuts, "AntiProton");
}

void TrackQA::PlotKinematicTracks() {
  PlotKinematic(fTrackCuts, "Proton");
}

void TrackQA::PlotKinematicAntiTracks() {
  PlotKinematic(fAntiTrackCuts, "AntiProton");
}

int TrackQA::GetNumberOfTracks() const {
  TH1F* counter = nullptr;
  if (fTrackCuts && fAntiTrackCuts) {
    std::cerr
        << "ERROR: You set both track and anti-track cuts - not defined (yet ;) )\n";
    return -1;
  } else if (!fTrackCuts && !fAntiTrackCuts) {
    std::cerr << "ERROR: No Track List set - can't count\n";
    return -1;
  }

  if (fTrackCuts) {
    auto* pTDist = (TH1F*) fReader->Get1DHistInList(
        fReader->GetListInList(fTrackCuts, { "after" }), "pTDist_after");
    counter = pTDist;
  } else if (fAntiTrackCuts) {
    auto* pTDist = (TH1F*) fReader->Get1DHistInList(
        fReader->GetListInList(fAntiTrackCuts, { "after" }), "pTDist_after");
    counter = pTDist;
  }
  return counter->GetEntries();
}

void TrackQA::PlotKinematic(TList* cuts, const char* outname) {
  auto* pTDist = (TH1F*) fReader->Get1DHistInList(
      fReader->GetListInList(cuts, { "after" }), "pTDist_after");
  if (!pTDist) {
    std::cerr << "pT Distribution missing for " << outname << std::endl;
  }
  pTDist->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  pTDist->GetYaxis()->SetTitle(
      Form("Entries/ %.1f GeV/#it{c}", pTDist->GetBinWidth(1)));
  pTDist->SetMinimum(pTDist->GetBinContent(pTDist->FindLastBinAbove(1000)) * 0.5);
  fHairyPlotter->FormatHistogram(pTDist, fStyler);
  fHairyPlotter->DrawLogYAndStore( { pTDist }, Form("%s_pT", outname));

  auto* phiDist = (TH1F*) fReader->Get1DHistInList(
      fReader->GetListInList(cuts, { "after" }), "phiDist_after");
  if (!phiDist) {
    std::cerr << "phi Distribution missing for " << outname << std::endl;
  }
  phiDist->SetMinimum(0);
  phiDist->GetXaxis()->SetTitle("#phi (rad)");
  phiDist->GetYaxis()->SetTitle(
      Form("Entries/ %.3f rad", phiDist->GetBinWidth(1)));
  fHairyPlotter->FormatHistogram(phiDist, fStyler);
  fHairyPlotter->DrawAndStore( { phiDist }, Form("%s_phi", outname));

  auto* etaDist = (TH1F*) fReader->Get1DHistInList(
      fReader->GetListInList(cuts, { "after" }), "EtaDist_after");
  if (!etaDist) {
    std::cerr << "eta Distribution missing for " << outname << std::endl;
  }
  etaDist->GetXaxis()->SetTitle("#eta");
  etaDist->GetXaxis()->SetRangeUser(-1., 1.);
  etaDist->GetYaxis()->SetTitle(
      Form("Entries/ %.2f ", etaDist->GetBinWidth(1)));
  fHairyPlotter->FormatHistogram(etaDist, fStyler);
  fHairyPlotter->DrawAndStore( { etaDist }, Form("%s_eta", outname), "hist");

  auto* dcaXY2DDistBef = (TH2F*) fReader->Get2DHistInList(cuts,
                                                          "DCAXYPtBinningTot");
  TH1F* dcaXYDistBef = nullptr;
  if (!dcaXY2DDistBef) {
    std::cerr << "DCAXY Before Distribution missing for " << outname
              << std::endl;
  } else {
    dcaXYDistBef = (TH1F*) dcaXY2DDistBef->ProjectionY(
        Form("dcaXYProjBef_%s", outname));
    fHairyPlotter->FormatHistogram(dcaXYDistBef, fStyler);
  }
  auto* dcaXY2DDistAfter = (TH2F*) fReader->Get2DHistInList(
      fReader->GetListInList(cuts, { "after" }), "DCAXYProp_after");
  if (!dcaXY2DDistAfter) {
    std::cerr << "DCAXY After Distribution missing for " << outname
              << std::endl;
  }
  auto dcaXYDistAfter = (TH1F*) dcaXY2DDistAfter->ProjectionY(
      Form("dcaXYProjAfter_%s", outname));


  dcaXYDistAfter->GetXaxis()->SetTitle("DCA_{#it{xy}} (cm)");

  dcaXYDistAfter->GetXaxis()->SetRangeUser(-4., 4.);
  dcaXYDistAfter->GetYaxis()->SetRangeUser(1.e3, 1.e8);

  dcaXYDistAfter->GetYaxis()->SetTitle(
      Form("Entries/ %.2f cm ", dcaXYDistAfter->GetBinWidth(1)));
  fHairyPlotter->FormatHistogram(dcaXYDistAfter, fStyler.drawMarker, fStyler.drawColor - 1 );
  std::vector<TH1*> dcaPlot = { dcaXYDistAfter };
  if (dcaXYDistBef) {
    dcaPlot.push_back(dcaXYDistBef);
  }
  fHairyPlotter->DrawLogYAndStore(dcaPlot, Form("%s_dca", outname));
}

void TrackQA::PlotPID() {
  PlotPID(fTrackCuts, "Proton");
  PlotPID(fAntiTrackCuts, "AntiProton");
}

void TrackQA::PlotPIDTracks() {
  PlotPID(fTrackCuts, "Proton");
}

void TrackQA::PlotPIDAntiTracks() {
  PlotPID(fAntiTrackCuts, "AntiProton");
}


void TrackQA::PlotPID(TList* cuts, const char* outname) {
  auto* TPCdEdx = (TH2F*) fReader->Get2DHistInList(
      fReader->GetListInList(cuts, { "after" }), "NSigTPC_after");
  if (!TPCdEdx) {
    std::cerr << "TPC dEdx missing for " << outname << '\n';
  }
  fHairyPlotter->FormatHistogram(TPCdEdx);
  TPCdEdx->GetYaxis()->SetRangeUser(-5, 5);
  TPCdEdx->GetXaxis()->SetTitle("#it{p}_{TPC} (GeV/#it{c})");
  TPCdEdx->GetYaxis()->SetTitle("#it{n}_{#sigma, TPC}");
  std::vector<TH2*> drawVecTPC = { TPCdEdx };
  fHairyPlotter->DrawLogZAndStore(drawVecTPC, Form("%s_LognSigTPC", outname), "colz");
  fHairyPlotter->DrawAndStore(drawVecTPC, Form("%s_nSigTPC", outname), "COLZ");

  auto* TOF = (TH2F*) fReader->Get2DHistInList(
      fReader->GetListInList(cuts, { "after" }), "NSigTOF_after");
  if (!TOF) {
    std::cerr << "TOF missing for " << outname << '\n';
  }
  fHairyPlotter->FormatHistogram(TOF);
  TOF->GetYaxis()->SetRangeUser(-5, 5);
  TOF->GetXaxis()->SetTitle("#it{p}_{TPC} (GeV/#it{c})");
  TOF->GetYaxis()->SetTitle("#it{n}_{#sigma, TOF}");
  std::vector<TH2*> drawVecTOF = { TOF };
  fHairyPlotter->DrawLogZAndStore(drawVecTOF, Form("%s_LognSigTOF", outname), "COLZ");
  fHairyPlotter->DrawAndStore(drawVecTOF, Form("%s_nSigTOF", outname), "COLZ");

}

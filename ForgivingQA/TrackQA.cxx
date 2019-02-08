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
  // TODO Auto-generated constructor stub

}

TrackQA::~TrackQA() {
  // TODO Auto-generated destructor stub
}

void TrackQA::PlotKinematic() {
  PlotKinematic(fTrackCuts, "Proton");
  PlotKinematic(fAntiTrackCuts, "AntiProton");
}

void TrackQA::PlotKinematic(TList* cuts, const char* outname) {
  auto* pTDist = (TH1F*) fReader->Get1DHistInList(
      fReader->GetListInList(cuts, { { "after" } }), "pTDist_after");
  if (!pTDist) {
    std::cerr << "pT Distribution missing for " << outname << std::endl;
  }
  pTDist->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  pTDist->GetYaxis()->SetTitle(
      Form("Entries/ %.1f GeV/#it{c}", pTDist->GetBinWidth(1)));
  fHairyPlotter->FormatHistogram(pTDist, 0, 1);
  fHairyPlotter->DrawLogYAndStore(pTDist, Form("%s_pT",outname));

  auto* phiDist = (TH1F*) fReader->Get1DHistInList(
      fReader->GetListInList(cuts, { { "after" } }), "phiDist_after");
  if (!phiDist) {
    std::cerr << "phi Distribution missing for " << outname << std::endl;
  }
  phiDist->GetXaxis()->SetTitle("#varphi (rad)");
  phiDist->GetYaxis()->SetTitle(
      Form("Entries/ %.3f rad", phiDist->GetBinWidth(1)));
  fHairyPlotter->FormatHistogram(phiDist, 0, 1);
  fHairyPlotter->DrawAndStore(phiDist, Form("%s_phi",outname));
}

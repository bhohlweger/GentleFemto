/*
 * EventQA.cxx
 *
 *  Created on: Feb 8, 2019
 *      Author: schmollweger
 */

#include "EventQA.h"
#include "TCanvas.h"
#include <stdlib.h>
#include <iostream>
EventQA::EventQA()
    : fReader(),
      fQA(nullptr),
      fEventCuts(nullptr),
      fHairyPlotter(new MakeHistosGreat()),
      fNEvts(-1) {
}

EventQA::EventQA(const char* outname)
    : fReader(),
      fQA(nullptr),
      fEventCuts(nullptr),
      fHairyPlotter(new MakeHistosGreat(outname)),
      fNEvts(-1) {
}

EventQA::~EventQA() {
  delete fHairyPlotter;
}

void EventQA::MakeEventQA() {
  PlotCutCounter();
  PlotEventProperties(200);
}

void EventQA::PlotCutCounter() {
  auto* cutStats = (TH1D*) fReader->Get1DHistInList(
      fReader->GetListInList(fQA, { { "AliEventCuts" } }), "fCutStats");
  if (!cutStats) {
    std::cerr << "PlotCutCounter: No cutStat Hist\n";
  }
  fHairyPlotter->FormatHistogram(cutStats, fStyler);
  fHairyPlotter->DrawAndStore( { cutStats }, "Evt_CutStats");
}

void EventQA::PlotEventProperties(unsigned int multMax) {
  //Method also Sets the Number of Events!

  auto* multiplicity = (TH1F*) fReader->Get1DHistInList(
      fReader->GetListInList(fEventCuts, { { "after" } }),
      "MultiplicityRef08_after");
  if (!multiplicity) {
    std::cerr << "PlotEventProperties: Missing Multiplicity Histogram! \n";
  }
  multiplicity->GetXaxis()->SetRangeUser(0, multMax);
  multiplicity->GetYaxis()->SetTitle(Form("#it{N}_{events}"));
  fHairyPlotter->FormatHistogram(multiplicity, fStyler);
  fHairyPlotter->DrawLogYAndStore( { multiplicity }, "EvtProp_eventMult");

  auto* zVtx = (TH1F*) fReader->Get1DHistInList(
      fReader->GetListInList(fEventCuts, { { "after" } }), "VtxZ_after");
  if (!zVtx) {
    std::cerr << "PlotEventProperties: Missing zVtx Histogram! \n";
  }
  fNEvts = zVtx->GetEntries();
  zVtx->GetXaxis()->SetTitle("Vertex #it{z} (cm)");
  zVtx->GetXaxis()->SetRangeUser(-12, 12);
  zVtx->GetYaxis()->SetTitle(Form("#it{N}_{events}/%.1f cm", zVtx->GetBinWidth(1)));
  fHairyPlotter->FormatHistogram(zVtx, fStyler);
  std::vector<TH1*> drawVec = { { zVtx } };
  fHairyPlotter->DrawAndStore(drawVec, "EvtProp_zVtx");
}

void EventQA::PlotPileUpRejection() {
  auto TkltsVsClusterBefore = (TH2F*) fReader->Get2DHistInList(
      fReader->GetListInList(fEventCuts, { "before" }),
      "SPDTrackletsVsClusterL01Sum_before");
  auto* TkltsVsClusterAfter= (TH2F*) fReader->Get2DHistInList(
      fReader->GetListInList(fEventCuts, { "after" }),
      "SPDTrackletsVsClusterL01Sum_after");
  if (TkltsVsClusterBefore && TkltsVsClusterAfter) {
    fHairyPlotter->FormatHistogram(TkltsVsClusterBefore);
    TkltsVsClusterBefore->GetXaxis()->SetTitle("#it{N}_{SPD Tracklets}");
    TkltsVsClusterBefore->GetYaxis()->SetTitle("#it{N}_{SPD Cluster}");
    std::vector<TH2*> drawVecBefore = { { TkltsVsClusterBefore} };
    fHairyPlotter->DrawLogZAndStore(drawVecBefore,"TrackletsVsClustBefore", "COLZ");

    fHairyPlotter->FormatHistogram(TkltsVsClusterAfter);
    TkltsVsClusterAfter->GetXaxis()->SetTitle("#it{N}_{SPD Tracklets}");
    TkltsVsClusterAfter->GetYaxis()->SetTitle("#it{N}_{SPD Cluster}");
    std::vector<TH2*> drawVecAfter = { { TkltsVsClusterAfter } };
    fHairyPlotter->DrawLogZAndStore(drawVecAfter,"TrackletsVsClustAfter", "COLZ");

  } else {
   std::cerr << "No Tracklets vs. Clust Plot \n";
  }
  return;
}

void EventQA::PlotStatsTrackCleaner(std::vector<const char*> TrackDecay,
                                    std::vector<const char*> DecayDecay,
                                    unsigned int xMax) {
  int itrDec = 0;
  for (auto it : TrackDecay) {
    auto* trackDecay = (TH1F*) fReader->Get1DHistInList(
        fReader->GetListInList(fQA, { { "PairCleaner" } }),
        Form("DaugthersSharedTracks_%u", itrDec));
    if (!trackDecay) {
      std::cerr << "PlotStatsTracksCleaner: "
                << Form("DaugthersSharedTracks_%u", itrDec) << " missing \n";
    }
    TString nameAxis = Form("%s pairs with shared tracks/event", it);
    trackDecay->GetXaxis()->SetTitle(nameAxis.Data());
    trackDecay->GetXaxis()->SetRangeUser(0, xMax);
    trackDecay->GetYaxis()->SetTitle("Entries");
    fHairyPlotter->FormatHistogram(trackDecay, fStyler);
    fHairyPlotter->DrawLogYAndStore( { trackDecay },
                                    Form("DecTrack_%u", itrDec));

    ++itrDec;
  }
  int iDecDec = 0;
  for (auto it : DecayDecay) {
    auto* DecayDecay = (TH1F*) fReader->Get1DHistInList(
        fReader->GetListInList(fQA, { { "PairCleaner" } }),
        Form("DaugthersSharedDaughters_%u", iDecDec));
    if (!DecayDecay) {
      std::cerr << "PlotStatsTracksCleaner: "
                << Form("DaugthersSharedDaughters_%u", iDecDec)
                << " missing \n";
    }
    TString nameAxis = Form("%s pairs with shared tracks/event", it);
    DecayDecay->GetXaxis()->SetTitle(nameAxis.Data());
    DecayDecay->GetYaxis()->SetTitle("Entries");
    DecayDecay->GetXaxis()->SetRangeUser(0, xMax);
    fHairyPlotter->FormatHistogram(DecayDecay, fStyler);
    fHairyPlotter->DrawLogYAndStore( { DecayDecay },
                                    Form("Decdec_%u", iDecDec));

    ++iDecDec;
  }
}

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
      fNEvts(0) {
  // TODO Auto-generated constructor stub
}

EventQA::~EventQA() {
  // TODO Auto-generated destructor stub
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
  fHairyPlotter->FormatHistogram(cutStats, 0, 1);
  fHairyPlotter->DrawAndStore({cutStats}, "Evt_CutStats");
}

void EventQA::PlotEventProperties(unsigned int multMax) {
  //Method also Sets the Number of Events!

  auto* multiplicity = (TH2F*) fReader->Get1DHistInList(
      fReader->GetListInList(fEventCuts, { { "after" } }),
      "MultiplicityRef08_after");
  if (!multiplicity) {
    std::cerr << "PlotEventProperties: Missing Multiplicity Histogram! \n";
  }
  multiplicity->GetXaxis()->SetRangeUser(0, multMax);
  multiplicity->GetYaxis()->SetTitle(Form("N_{Events}"));
  fHairyPlotter->FormatHistogram(multiplicity, 0, 1);
  fHairyPlotter->DrawLogYAndStore({multiplicity}, "EvtProp_eventMult");

  auto* zVtx = (TH2F*) fReader->Get1DHistInList(
      fReader->GetListInList(fEventCuts, { { "after" } }), "VtxZ_after");
  if (!zVtx) {
    std::cerr << "PlotEventProperties: Missing zVtx Histogram! \n";
  }
  fNEvts = zVtx->GetEntries();
  zVtx->GetXaxis()->SetTitle("v_{z} (cm)");
  zVtx->GetXaxis()->SetRangeUser(-11, 11);
  zVtx->GetYaxis()->SetTitle(Form("N_{Events}/%.1f cm", zVtx->GetBinWidth(1)));
  fHairyPlotter->FormatHistogram(zVtx, 0, 1);
  std::vector<TH2*> drawVec = {{zVtx}};
  fHairyPlotter->DrawAndStore(drawVec, "EvtProp_zVtx");
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
    TString nameAxis = Form("# %s pairs with shared tracks/event", it);
    trackDecay->GetXaxis()->SetTitle(nameAxis.Data());
    trackDecay->GetXaxis()->SetRangeUser(0, xMax);
    fHairyPlotter->FormatHistogram(trackDecay, 0, 1);
    fHairyPlotter->DrawLogYAndStore({trackDecay}, Form("DecTrack_%u", itrDec));

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
    TString nameAxis = Form("# %s pairs with shared tracks/event", it);
    DecayDecay->GetXaxis()->SetTitle(nameAxis.Data());
    DecayDecay->GetXaxis()->SetRangeUser(0, xMax);
    fHairyPlotter->FormatHistogram(DecayDecay, 0, 1);
    fHairyPlotter->DrawLogYAndStore({DecayDecay}, Form("Decdec_%u", iDecDec));

    ++iDecDec;
  }
}

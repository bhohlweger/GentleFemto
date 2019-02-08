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
      fHairyPlotter(new MakeHistosGreat()) {
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
  fHairyPlotter->DrawAndStore(cutStats, "Evt_CutStats");
}

void EventQA::PlotEventProperties(unsigned int multMax) {
  auto* multiplicity = (TH2F*) fReader->Get1DHistInList(
      fReader->GetListInList(fEventCuts, { {"after"} }),
      "MultiplicityRef08_after");
  multiplicity->GetXaxis()->SetRangeUser(0, multMax);
  fHairyPlotter->FormatHistogram(multiplicity, 0, 1);
  fHairyPlotter->DrawLogYAndStore(multiplicity, "EvtProp_eventMult");

  auto* zVtx = (TH2F*) fReader->Get1DHistInList(
      fReader->GetListInList(fEventCuts, { {"after"} }),
      "VtxZ_after");
  zVtx->GetXaxis()->SetTitle("v_{z} (cm)");
  zVtx->GetXaxis()->SetRangeUser(-11, 11);
  fHairyPlotter->FormatHistogram(zVtx, 0, 1);
  fHairyPlotter->DrawAndStore(zVtx, "EvtProp_zVtx");

}


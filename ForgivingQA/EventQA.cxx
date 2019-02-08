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
      fPainter(new MakeHistosGreat()) {
  // TODO Auto-generated constructor stub
}

EventQA::~EventQA() {
  // TODO Auto-generated destructor stub
}

void EventQA::MakeEventQA() {
  PlotCutCounter();
}

void EventQA::PlotCutCounter() {
  TH1D* cutStats = (TH1D*) fReader->Get1DHistInList(
      fReader->GetListInList(fQA, { { "AliEventCuts" } }), "fCutStats");
  if (!cutStats) {
    std::cerr << "PlotCutCounter: No cutStat Hist\n";
  }
  fPainter->FormatHistogram(cutStats,0,1);
  fPainter->DrawAndStore(cutStats,"CutStats");
}

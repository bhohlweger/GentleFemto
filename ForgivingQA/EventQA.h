/*
 * EventQA.h
 *
 *  Created on: Feb 8, 2019
 *      Author: schmollweger
 */

#ifndef FORGIVINGQA_EVENTQA_H_
#define FORGIVINGQA_EVENTQA_H_
#include "TList.h"
#include "ForgivingReader.h"
#include "MakeHistosGreat.h"
#include <vector>

class EventQA {
 public:
  EventQA();
  virtual ~EventQA();
  void GetCutStats();
  void SetQAList(TList* QA) {
    fQA = QA;
  }
  ;
  void SetEventCuts(TList* EvtCuts) {
    fEventCuts = EvtCuts;
  }
  ;
  void MakeEventQA();
  void PlotCutCounter();
  void SetTightMargin() {
    fHairyPlotter->SetTightMargin(true);
  }
  ;
  void SetLooseMargin() {
    fHairyPlotter->SetTightMargin(false);
  }
  ;
  void PlotEventProperties(unsigned int multMax);
  void PlotPileUpRejection();
  void PlotStatsTrackCleaner(std::vector<const char*> TrackDecay,
                             std::vector<const char*> DecayDecay,
                             unsigned int xMax);
  unsigned int GetNumberOfEvents();
 private:
  ForgivingReader* fReader;
  TList* fQA;
  TList* fEventCuts;
  MakeHistosGreat* fHairyPlotter;
  int fNEvts;
};

inline
unsigned int EventQA::GetNumberOfEvents() {
  if (fNEvts > 0) {
    return fNEvts;
  } else {
    auto* zVtx = (TH2F*) fReader->Get1DHistInList(
        fReader->GetListInList(fEventCuts, { { "after" } }), "VtxZ_after");
    if (!zVtx) {
      std::cerr << "GetNumberOfEvents: Missing zVtx Histogram! \n";
      return 0;
    }
    return zVtx->GetEntries();
  }
}

#endif /* FORGIVINGQA_EVENTQA_H_ */

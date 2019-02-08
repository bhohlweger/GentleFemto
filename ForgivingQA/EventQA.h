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
  void PlotStatsTrackCleaner(std::vector<const char*> TrackDecay,
                             std::vector<const char*> DecayDecay,
                             unsigned int xMax);
  unsigned int GetNumberOfEvents() {
    return fNEvts;  // only works if the zVtx hist is set & was plotted!
  }
  ;
 private:
  ForgivingReader* fReader;
  TList* fQA;
  TList* fEventCuts;
  MakeHistosGreat* fHairyPlotter;
  unsigned int fNEvts;
};

#endif /* FORGIVINGQA_EVENTQA_H_ */

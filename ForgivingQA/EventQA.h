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
 private:
  ForgivingReader* fReader;
  TList* fQA;
  TList* fEventCuts;
  MakeHistosGreat* fPainter;
};

#endif /* FORGIVINGQA_EVENTQA_H_ */

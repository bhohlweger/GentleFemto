/*
 * TrackQA.h
 *
 *  Created on: Feb 8, 2019
 *      Author: schmollweger
 */

#ifndef FORGIVINGQA_TRACKQA_H_
#define FORGIVINGQA_TRACKQA_H_
#include "TList.h"
#include "ForgivingReader.h"
#include "MakeHistosGreat.h"

class TrackQA {
 public:
  TrackQA();
  TrackQA(const char* outname);
  virtual ~TrackQA();
  void SetStyler(DrawStyle styler) { fStyler = styler; }
  void PlotKinematic();
  void PlotKinematic(TList *cuts, const char* outname);
  void PlotPID();
  void PlotPID(TList* cuts, const char* outname);
  void PlotKinematicTracks();
  void PlotPIDTracks();
  void PlotKinematicAntiTracks();
  void PlotPIDAntiTracks();
  int GetNumberOfTracks() const;
  void SetTrackCuts(TList* trkCuts) {
    fTrackCuts = trkCuts;
  }
  ;
  void SetAntiTrackCuts(TList* trkCuts) {
    fAntiTrackCuts = trkCuts;
  }
  ;
 private:
  ForgivingReader* fReader;
  MakeHistosGreat* fHairyPlotter;
  DrawStyle fStyler;
  TList *fTrackCuts;
  TList *fAntiTrackCuts;
};

#endif /* FORGIVINGQA_TRACKQA_H_ */

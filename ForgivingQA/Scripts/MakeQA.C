#include "ForgivingReader.h"
#include "MakeHistosGreat.h"
#include "EventQA.h"
#include "TrackQA.h"
#include "DecayQA.h"

int main(int argc, char* argv[]) {
  const char* filename = argv[1];
  const char* prefix = argv[2];
  const char* addon = (argv[3]) ? argv[3] : "";
  MakeHistosGreat::SetStyle(false);
  ForgivingReader* reader = new ForgivingReader(filename, prefix, addon);
  EventQA* evtQA = new EventQA();
  evtQA->SetLooseMargin();
  evtQA->SetQAList(reader->GetQA());
  evtQA->SetEventCuts(reader->GetEventCuts());

  evtQA->PlotCutCounter();
  evtQA->PlotEventProperties(200);
  evtQA->PlotPileUpRejection();
  evtQA->SetTightMargin();
  evtQA->PlotStatsTrackCleaner( { "p-#Lambda", "#bar{p}-#bar{#Lambda}", "p-#Xi",
                                   "#bar{p}-#bar{#Xi}" },
                               { "#Lambda-#Lambda",
                                   "#bar{#Lambda}-#bar{#Lambda}", "#Xi-#Xi",
                                   "#bar{#Xi}-#bar{#Xi}" },
                               6);

  TrackQA* trkQA = new TrackQA();
  trkQA->SetTrackCuts(reader->GetTrackCuts());
  trkQA->SetAntiTrackCuts(reader->GetAntiTrackCuts());
//  trkQA->PlotKinematic();
  trkQA->PlotPID();

  DecayQA* v0QA = new DecayQA("#Lambda","p#pi");
  v0QA->SetCanvasDivisions(4, 2);
  v0QA->SetDecayCuts(reader->Getv0Cuts());
  v0QA->SetIMHistoScale(1.75,0.8,0.35);
  v0QA->SetAntiDecayCuts(reader->GetAntiv0Cuts());
  v0QA->SetRangesFitting(1.1075, 1.1235, 1.09, 1.15);
  v0QA->InvariantMassLambda(1.112, 1.120);
  v0QA->PlotQATopologyLambda();

  DecayQA* cascQA = new DecayQA("#Xi^{-}","#pi#Lambda");
  cascQA->SetCanvasDivisions(4, 3);
  cascQA->SetInvMasspTStartBin(2);
  cascQA->SetIMHistoScale(2.5,0.8,0.45);
  cascQA->SetDecayCuts(reader->GetCascadeCuts());
  cascQA->SetAntiDecayCuts(reader->GetAntiCascadeCuts());
  cascQA->SetRangesFitting(1.31, 1.33, 1.3, 1.35);
  cascQA->InvariantMassXi(1.317, 1.327);
  cascQA->IvariantMassXiLambda();
  return 0;
}


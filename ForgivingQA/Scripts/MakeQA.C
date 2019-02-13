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
  evtQA->SetTightMargin();
  evtQA->PlotStatsTrackCleaner( { "p-#Lambda", "#bar{p}-#bar{#Lambda}", "p-#Xi",
                                   "#bar{p}-#bar{#Xi}" },
                               { "#Lambda-#Lambda",
                                   "#bar{#Lambda}-#bar{#Lambda}", "#Xi-#Xi",
                                   "#bar{#Xi}-#bar{#Xi}" }, 6);

  TrackQA* trkQA = new TrackQA();
  trkQA->SetTrackCuts(reader->GetTrackCuts());
  trkQA->SetAntiTrackCuts(reader->GetAntiTrackCuts());

  trkQA->PlotKinematic();
  trkQA->PlotPID();

  DecayQA* v0QA = new DecayQA();
  v0QA->SetDecayCuts(reader->Getv0Cuts());
  v0QA->SetAntiDecayCuts(reader->GetAntiv0Cuts());
  v0QA->SetRangesFitting();
  v0QA->InvariantMassLambda(1.111,1.119);
  return 0;
}


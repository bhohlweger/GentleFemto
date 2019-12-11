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

//  evtQA->PlotCutCounter();
  evtQA->PlotEventProperties(200);
  evtQA->PlotPileUpRejection();
  evtQA->SetTightMargin();
  evtQA->PlotStatsTrackCleaner( {"p-#bar{#Lambda}", "#bar{p}-#Lambda"
                                },{"#Lambda-#bar{#Lambda}" }, 6);

  TrackQA* trkQA = new TrackQA();
  trkQA->SetTrackCuts(reader->GetTrackCuts());
  trkQA->PlotKinematicTracks();
  trkQA->PlotPIDTracks();

  TrackQA* antitrkQA = new TrackQA();
  antitrkQA->SetAntiTrackCuts(reader->GetAntiTrackCuts());
  antitrkQA->PlotKinematicAntiTracks();
  antitrkQA->PlotPIDAntiTracks();

  DecayQA* v0QA = new DecayQA("#Lambda","p#pi");
  v0QA->SetCanvasDivisions(4, 2);
  v0QA->SetDecayCuts(reader->Getv0Cuts());
  v0QA->SetIMHistoScale(1.75,0.8,0.35);
  v0QA->SetAntiDecayCuts(reader->GetAntiv0Cuts());
  v0QA->SetRangesFitting(1.1075, 1.1235, 1.09, 1.15);
  v0QA->InvariantMassPartLambda(1.112, 1.120);
  v0QA->PlotQATopologyPartLambda();

  DecayQA* antiv0QA = new DecayQA("#bar{#Lambda}","#bar{p}#pi");
  antiv0QA->SetCanvasDivisions(4, 2);
  antiv0QA->SetIMHistoScale(1.75,0.8,0.35);
  antiv0QA->SetAntiDecayCuts(reader->GetAntiv0Cuts());
  antiv0QA->SetRangesFitting(1.1075, 1.1235, 1.09, 1.15);
  antiv0QA->InvariantMassAntiPartLambda(1.112, 1.120);
  antiv0QA->PlotQATopologyAntiPartLambda();


  return 0;
}

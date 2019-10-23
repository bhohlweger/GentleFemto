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
  std::cout<<"debug 1"<<std::endl;
  evtQA->SetQAList(reader->GetQA());
  std::cout<<"debug 2"<<std::endl;
  evtQA->SetEventCuts(reader->GetEventCuts());
  std::cout<<"debug 3"<<std::endl;

//  evtQA->PlotCutCounter();
  std::cout<<"debug 4"<<std::endl;
  evtQA->PlotEventProperties(200);
  std::cout<<"debug 5"<<std::endl;
  evtQA->PlotPileUpRejection();
  std::cout<<"debug 6"<<std::endl;
  evtQA->SetTightMargin();
  std::cout<<"debug 7"<<std::endl;
  evtQA->PlotStatsTrackCleaner( {"p-#bar{#Lambda}", "#bar{p}-#Lambda",
                                "#Lambda-#bar{#Lambda}"},{ }, 6);
  std::cout<<"debug 8"<<std::endl;

  TrackQA* trkQA = new TrackQA();
  trkQA->SetTrackCuts(reader->GetTrackCuts());
  std::cout<<"debug 9"<<std::endl;
  trkQA->PlotKinematic();
  std::cout<<"debug 10"<<std::endl;
  trkQA->PlotPID();
  std::cout<<"debug 11"<<std::endl;


  TrackQA* antitrkQA = new TrackQA();
  antitrkQA->SetAntiTrackCuts(reader->GetAntiTrackCuts());
  antitrkQA->PlotKinematic();
  antitrkQA->PlotPID();

  DecayQA* v0QA = new DecayQA("#Lambda","p#pi");
  v0QA->SetCanvasDivisions(4, 2);
  v0QA->SetDecayCuts(reader->Getv0Cuts());
  v0QA->SetIMHistoScale(1.75,0.8,0.35);
  v0QA->SetAntiDecayCuts(reader->GetAntiv0Cuts());
  v0QA->SetRangesFitting(1.1075, 1.1235, 1.09, 1.15);
  v0QA->InvariantMassLambda(1.112, 1.120);
  v0QA->PlotQATopologyLambda();

  DecayQA* antiv0QA = new DecayQA("#bar{#Lambda}","#bar{p}#pi");
  antiv0QA->SetCanvasDivisions(4, 2);
  antiv0QA->SetIMHistoScale(1.75,0.8,0.35);
  antiv0QA->SetAntiDecayCuts(reader->GetAntiv0Cuts());
  antiv0QA->SetRangesFitting(1.1075, 1.1235, 1.09, 1.15);
  antiv0QA->InvariantMassLambda(1.112, 1.120);
  antiv0QA->PlotQATopologyLambda();


  return 0;
}

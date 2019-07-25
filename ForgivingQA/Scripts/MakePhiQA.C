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
  auto file = reader->GetFile();
  auto dir = file->GetDirectory(Form("%sResults%s", prefix, addon));
  TList* list;
  dir->GetObject(Form("%sResults%s", prefix, addon), list);

  EventQA* evtQA = new EventQA();
  evtQA->SetLooseMargin();
  evtQA->SetQAList(list);
  evtQA->SetEventCuts((TList*) list->FindObject("Event Cuts"));

  evtQA->PlotCutCounter();
  evtQA->PlotEventProperties(200);
  evtQA->PlotPileUpRejection();
  evtQA->SetTightMargin();
  evtQA->PlotStatsTrackCleaner( { "p-#bar{p}", "p-#varphi", "#bar{p}-#varphi" },
                               { "#varphi-#varphi" }, 7);

  // Protons
  TrackQA* protonQA = new TrackQA();
  protonQA->SetTrackCuts((TList*) list->FindObject("Proton"));
  protonQA->SetAntiTrackCuts((TList*) list->FindObject("AntiProton"));
  protonQA->PlotKinematic();
  protonQA->PlotPID();

  // analogous for the Kaons
  TrackQA* kaonQA = new TrackQA();
  kaonQA->SetTrackCuts((TList*) list->FindObject("Particle1"));
  kaonQA->SetAntiTrackCuts((TList*) list->FindObject("Particle2"));
  kaonQA->PlotKinematic(((TList*) list->FindObject("Particle2")), "AntiKaon");
  kaonQA->PlotPID(((TList*) list->FindObject("Particle2")), "AntiKaon");
  kaonQA->PlotKinematic(((TList*) list->FindObject("Particle1")), "Kaon");
  kaonQA->PlotPID(((TList*) list->FindObject("Particle1")), "Kaon");

  DecayQA* v0QA = new DecayQA("#varphi","K^{-}K^{+}");
  v0QA->SetCanvasDivisions(4, 2);
  v0QA->PlotQATopologyLambda((TList*)list->FindObject("Phi"), "Phi");
}

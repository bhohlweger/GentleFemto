#include "DecayQA.h"
#include "EventQA.h"
#include "ForgivingReader.h"
#include "MakeHistosGreat.h"
#include "TrackQA.h"

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
  evtQA->PlotStatsTrackCleaner( { "p-#Sigma^{0}", "#bar{p}-#bar{#Sigma^{0}}",
                                   "p-#Lambda#gamma (up)",
                                   "#bar{p}-#bar{#Lambda}#gamma (up)",
                                   "p-#Lambda#gamma (down)",
                                   "#bar{p}-#bar{#Lambda}#gamma (down)" },
                               { }, 6);

  TrackQA* trkQA = new TrackQA();
  trkQA->SetTrackCuts(reader->GetTrackCuts());
  trkQA->SetAntiTrackCuts(reader->GetAntiTrackCuts());
  trkQA->PlotKinematic();
  trkQA->PlotPID();

  DecayQA* v0QA = new DecayQA("#Lambda", "p#pi");
  v0QA->SetCanvasDivisions(5, 2);
  v0QA->SetDecayCuts(reader->Getv0Cuts());
  v0QA->SetIMHistoScale(1.75, 0.8, 0.35);
  v0QA->SetAntiDecayCuts(reader->GetAntiv0Cuts());
  v0QA->SetRangesFitting(1.109, 1.121, 1.09, 1.15);
  v0QA->InvariantMassLambdaSigma0(1.112, 1.120);
  v0QA->PlotQATopologySigma0Daughter(reader->Getv0Cuts(), "Lambda");
  v0QA->PlotPIDSigma0Daughter(reader->Getv0Cuts(), "Lambda");
  v0QA->PlotQATopologySigma0Daughter(reader->GetAntiv0Cuts(), "AntiLambda");
  v0QA->PlotPIDSigma0Daughter(reader->GetAntiv0Cuts(), "AntiLambda");

  DecayQA* gammaQA = new DecayQA("#gamma", "e^{+}e^{-}");
  gammaQA->PlotQATopologySigma0Daughter(reader->GetOtherCuts("PhotonCuts"),
                                        "Photon");
  gammaQA->PlotPIDSigma0Daughter(reader->GetOtherCuts("PhotonCuts"), "Photon");

  DecayQA* sigma0QA = new DecayQA("#Sigma^{0}", "#Lambda#gamma");
  sigma0QA->SetDecayCuts(reader->GetOtherCuts("Sigma0Cuts"));
  sigma0QA->SetCanvasDivisions(3, 2);
  sigma0QA->SetIMHistoScale(1.75, 0.8, 0.35);
  sigma0QA->PlotQATopologySigma0(reader->GetOtherCuts("Sigma0Cuts"),
                                 "Sigma0part");
  sigma0QA->SetRangesFitting(1.19, 1.196, 1.167, 1.217);
  sigma0QA->InvariantMassSigma0(0.003, "Sigma0part", false);
  delete sigma0QA;

  DecayQA* antiSigma0QA = new DecayQA("#bar{#Sigma^{0}}",
                                      "#bar{#Lambda}#gamma");
  antiSigma0QA->SetAntiDecayCuts(reader->GetOtherCuts("AntiSigma0Cuts"));
  antiSigma0QA->SetCanvasDivisions(3, 2);
  antiSigma0QA->SetIMHistoScale(1.75, 0.8, 0.35);
  antiSigma0QA->PlotQATopologySigma0(reader->GetOtherCuts("AntiSigma0Cuts"),
                                     "Sigma0antiPart");
  antiSigma0QA->SetRangesFitting(1.19, 1.196, 1.167, 1.217);
  antiSigma0QA->InvariantMassSigma0(0.003, "Sigma0antiPart", false);
  delete antiSigma0QA;

  DecayQA* sigmaSumQA = new DecayQA("#Sigma^{0} + #bar{#Sigma^{0}}",
                                    "#Lambda#gamma + #bar{#Lambda}#gamma");
  sigmaSumQA->SetDecayCuts(reader->GetOtherCuts("Sigma0Cuts"));
  sigmaSumQA->SetAntiDecayCuts(reader->GetOtherCuts("AntiSigma0Cuts"));
  sigmaSumQA->SetCanvasDivisions(3, 2);
  sigmaSumQA->SetIMHistoScale(1.75, 0.8, 0.35);
  sigmaSumQA->PlotQATopologySigma0(reader->GetOtherCuts("Sigma0Cuts"),
                                   "SigmaSum");
  sigmaSumQA->SetRangesFitting(1.19, 1.196, 1.167, 1.217);
  sigmaSumQA->InvariantMassSigma0(0.003);
  delete sigmaSumQA;
}

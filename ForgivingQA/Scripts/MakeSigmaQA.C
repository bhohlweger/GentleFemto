#include "DecayQA.h"
#include "EventQA.h"
#include "ForgivingReader.h"
#include "MakeHistosGreat.h"
#include "TrackQA.h"

int main(int argc, char* argv[]) {
  gROOT->ProcessLine("gErrorIgnoreLevel = 3001");
  const char* filename = argv[1];
  const char* prefix = argv[2];
  const char* addon = (argv[3]) ? argv[3] : "";
  MakeHistosGreat::SetStyle(false);
  ForgivingReader* reader = new ForgivingReader(filename, prefix, addon);
  TString suffix = TString::Format("%s", addon);

  DrawStyle styler;
  styler.drawMarker = 2;
  styler.drawColor = 8;
  styler.drawSize = 1.1;
  styler.drawLineColor = kTeal + 3;
  styler.drawSignalFitColor = kBlue + 4;
  styler.drawBackgroundFitColor = kGreen + 2;

  EventQA* evtQA = new EventQA("Event");
  evtQA->SetStyler(styler);
  evtQA->SetLooseMargin();
  evtQA->SetQAList(reader->GetQA());
  evtQA->SetEventCuts(reader->GetEventCuts());
  if (suffix == "0") {
    evtQA->PlotCutCounter();
    evtQA->PlotEventProperties(200);
    evtQA->PlotPileUpRejection();
    evtQA->SetTightMargin();
    evtQA->PlotStatsTrackCleaner(
        { "p#minus#kern[-0.95]{ }#Sigma^{0}",
            "#bar{p}#minus#kern[-0.85]{ }#bar{#Sigma^{0}}",
            "p#minus#kern[-0.65]{ }(#Lambda#gamma) (up)",
            "#bar{p}#minus#kern[-0.4]{ }(#bar{#Lambda}#gamma) (up)",
            "p#minus#kern[-0.65]{ }(#Lambda#gamma) (down)",
            "#bar{p}#minus#kern[-0.4]{ }(#bar{#Lambda}#gamma) (down)" },
        { "#Lambda#minus#kern[-0.4]{ }#Lambda", "#bar{#Lambda}#minus#kern[-0.4]{ }#bar{#Lambda}",
            "#gamma#minus#kern[-0.4]{ }#gamma", "#Lambda#minus#kern[-0.4]{ }#bar{#Lambda}",
            "#Lambda#minus#kern[-0.4]{ }#gamma", "#bar{#Lambda}#minus#kern[-0.4]{ }#gamma" },
        6);
  }
  delete evtQA;

  TrackQA* trkQA = new TrackQA("Proton");
  trkQA->SetStyler(styler);
  trkQA->SetTrackCuts(reader->GetTrackCuts());
  trkQA->SetAntiTrackCuts(reader->GetAntiTrackCuts());
  if (suffix == "0") {
    trkQA->PlotKinematic();
    trkQA->PlotPID();
  }
  delete trkQA;

  DecayQA* v0QA = new DecayQA("#Lambda", "p#pi", "Lambda");
  v0QA->SetStyler(styler);
  v0QA->SetCanvasDivisions(4, 2);
  v0QA->SetDecayCuts(reader->Getv0Cuts());
  v0QA->SetAntiDecayCuts(reader->GetAntiv0Cuts());
  v0QA->SetIMHistoScale(1.75, 0.8, 0.35);
  v0QA->SetRangesFitting(1.109, 1.121, 1.09, 1.15);
  v0QA->InvariantMassLambda(1.112, 1.120, (suffix != "0"), 0.493, 0.504);
  if (suffix == "0") {
    v0QA->PlotQATopologyLambda();
    v0QA->PlotPIDLambda();
  }
  delete v0QA;

  if (suffix == "0") {
    DecayQA* gammaQA = new DecayQA("#gamma", "e^{+}e^{-}", "Photon");
    gammaQA->SetStyler(styler);
    gammaQA->PlotQATopologySigma0Daughter(reader->GetOtherCuts("PhotonCuts"),
                                          "Photon");
    gammaQA->PlotPIDSigma0Daughter(reader->GetOtherCuts("PhotonCuts"),
                                   "Photon");
    delete gammaQA;
  }

  DecayQA* sigma0QA = new DecayQA("#Sigma^{0}", "#Lambda#gamma", "Sigma0");
  sigma0QA->SetStyler(styler);
  sigma0QA->SetDecayCuts(reader->GetOtherCuts("Sigma0Cuts"));
  sigma0QA->SetCanvasDivisions(3, 3);
  sigma0QA->SetInvMasspTStartBin(3);
  sigma0QA->SetIMHistoScale(1.75, 0.8, 0.35);
  if (suffix == "0") {
    sigma0QA->PlotQATopologySigma0(reader->GetOtherCuts("Sigma0Cuts"),
                                   "Sigma0part");
  }
  sigma0QA->SetRangesFitting(1.19, 1.196, 1.167, 1.217);
  sigma0QA->InvariantMassSigma0(0.003, "Sigma0part", false);
  delete sigma0QA;

  DecayQA* antiSigma0QA = new DecayQA("#bar{#Sigma^{0}}",
                                      "#bar{#Lambda}#gamma", "AntiSigma0");
  antiSigma0QA->SetStyler(styler);
  antiSigma0QA->SetAntiDecayCuts(reader->GetOtherCuts("AntiSigma0Cuts"));
  antiSigma0QA->SetCanvasDivisions(3, 3);
  antiSigma0QA->SetInvMasspTStartBin(3);
  antiSigma0QA->SetIMHistoScale(1.75, 0.8, 0.35);
  if (suffix == "0") {
    antiSigma0QA->PlotQATopologySigma0(reader->GetOtherCuts("AntiSigma0Cuts"),
                                       "Sigma0antiPart");
  }
  antiSigma0QA->SetRangesFitting(1.19, 1.196, 1.167, 1.217);
  antiSigma0QA->InvariantMassSigma0(0.003, "Sigma0antiPart", false);
  delete antiSigma0QA;

  DecayQA* sigmaSumQA = new DecayQA("#Sigma^{0} + #bar{#Sigma^{0}}",
                                    "#Lambda#gamma + #bar{#Lambda}#gamma", "Sigma0Combined");
  sigmaSumQA->SetStyler(styler);
  sigmaSumQA->SetDecayCuts(reader->GetOtherCuts("Sigma0Cuts"));
  sigmaSumQA->SetAntiDecayCuts(reader->GetOtherCuts("AntiSigma0Cuts"));
  sigmaSumQA->SetCanvasDivisions(3, 3);
  sigmaSumQA->SetInvMasspTStartBin(3);
  sigmaSumQA->SetIMHistoScale(1.75, 0.8, 0.35);
  if (suffix == "0") {
    sigmaSumQA->PlotQATopologySigma0(reader->GetOtherCuts("Sigma0Cuts"),
                                     "SigmaSum");
  }
  sigmaSumQA->SetRangesFitting(1.19, 1.196, 1.167, 1.217);
  sigmaSumQA->InvariantMassSigma0(0.003);
  delete sigmaSumQA;
}

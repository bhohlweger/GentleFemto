#include "PeriodQA.h"
#include "DecayQA.h"
#include "EventQA.h"
#include "TrackQA.h"
#include "TSystemDirectory.h"

PeriodQA::PeriodQA()
    : fDirectory(),
      fPeriods() {
}

void PeriodQA::SetDirectory(const char *dir) {
  fDirectory = dir;

  TSystemDirectory *workdir = new TSystemDirectory("workdir", dir);
  TList *list = workdir->GetListOfFiles();

  TIter next(list);
  TObject *obj = nullptr;
  while (obj = next()) {
    TString objName = obj->GetName();
    if (!objName.Contains("AnalysisResults"))
      continue;

    TString LHCperiod = objName;
    LHCperiod.ReplaceAll("AnalysisResults_", "");
    LHCperiod.ReplaceAll(".root", "");
    fPeriods.push_back(LHCperiod);
  }
  std::sort(fPeriods.begin(), fPeriods.end(),
            [](TString a, TString b) {return a<b;});
}

void PeriodQA::ProcessQA(const char* prefix, const char* addon) {
  auto histnEvents = PeriodQAHist("histQAnEvts", "Number of Events");

  auto histPurityLambda = PeriodQAHist("histQALambda", "Purity #Lambda (%)");
  auto histPurityAntiLambda = PeriodQAHist("histQAAntiLambda",
                                           "Purity #bar{#Lambda} (%)");
  auto histPurityXi = PeriodQAHist("histQAXi", "Purity #Xi (%)");
  auto histPurityAntiXi = PeriodQAHist("histQAAntiXi", "Purity #bar{#Xi} (%)");

  int i = 0;
  for (auto it : fPeriods) {
    TString filename = fDirectory.Data();
    filename += "/AnalysisResults_";
    filename += it;
    filename += ".root";
    ForgivingReader* reader = new ForgivingReader(filename.Data(), prefix,
                                                  addon);
    auto* multiplicity = (TH2F*) reader->Get1DHistInList(
        reader->GetListInList(reader->GetEventCuts(), { { "after" } }),
        "MultiplicityRef08_after");
    EventQA* evtQA = new EventQA();
    evtQA->SetEventCuts(reader->GetEventCuts());
    float nEvents = evtQA->GetNumberOfEvents();
    delete evtQA;
    if (nEvents < 1E-6) {
      delete reader;
      ++i;
      continue;
    }
    nEvents*=multiplicity->GetMean();
    histnEvents->SetBinContent(
        i + 1, float(multiplicity->GetMean()));
    histnEvents->SetBinError(
        i + 1, float(multiplicity->GetMeanError()));
    TrackQA* trkQA = new TrackQA();
    trkQA->SetTrackCuts(reader->GetTrackCuts());

    TrackQA* antiTrkQA = new TrackQA();
    antiTrkQA->SetAntiTrackCuts(reader->GetAntiTrackCuts());

    DecayQA* v0QA = new DecayQA("#Lambda", "p#pi");
    v0QA->SetDecayCuts(reader->Getv0Cuts());
    v0QA->SetCanvasDivisions(5, 2);
    v0QA->SetIMHistoScale(1.75, 0.8, 0.35);
    v0QA->SetRangesFitting(1.109, 1.121, 1.09, 1.15);
    v0QA->GetPeriodQA(1.112, 1.120, { "v0Cuts" }, "InvMassPt");
    float purity = v0QA->GetPurity();
    if (purity > 1E-6) {
      histPurityLambda->SetBinContent(i + 1, purity * 100.f);
    }

    DecayQA* antiv0QA = new DecayQA("#bar{#Lambda}", "p#pi");
    antiv0QA->SetDecayCuts(reader->GetAntiv0Cuts());
    antiv0QA->SetCanvasDivisions(5, 2);
    antiv0QA->SetIMHistoScale(1.75, 0.8, 0.35);
    antiv0QA->SetRangesFitting(1.109, 1.121, 1.09, 1.15);
    antiv0QA->GetPeriodQA(1.112, 1.120, { "v0Cuts" }, "InvMassPt");
    purity = antiv0QA->GetPurity();
    if (purity > 1E-6) {
      histPurityAntiLambda->SetBinContent(i + 1, purity * 100.f);
    }

    DecayQA* cascQA = new DecayQA("#Xi^{-}", "#pi#Lambda");
    cascQA->SetDecayCuts(reader->GetCascadeCuts());
    cascQA->SetCanvasDivisions(4, 3);
    cascQA->SetIMHistoScale(2.5, 0.8, 0.45);
    cascQA->SetRangesFitting(1.31, 1.33, 1.3, 1.35);
    cascQA->GetPeriodQA(1.317, 1.327, { "Cascade" }, "InvMassXi");
    purity = cascQA->GetPurity();
    if (purity > 1E-6) {
      histPurityXi->SetBinContent(i + 1, purity * 100.f);
    }

    DecayQA* anticascQA = new DecayQA("#Xi^{-}", "#pi#Lambda");
    anticascQA->SetDecayCuts(reader->GetAntiCascadeCuts());
    anticascQA->SetCanvasDivisions(4, 3);
    anticascQA->SetIMHistoScale(2.5, 0.8, 0.45);
    anticascQA->SetRangesFitting(1.31, 1.33, 1.3, 1.35);
    anticascQA->GetPeriodQA(1.317, 1.327, { "Cascade" }, "InvMassXi");
    purity = anticascQA->GetPurity();
    if (purity > 1E-6) {
      histPurityAntiXi->SetBinContent(i + 1, purity * 100.f);
    }
    delete trkQA;
    delete antiTrkQA;
    delete v0QA;
    delete antiv0QA;
    delete cascQA;
    delete anticascQA;
    delete reader;
    ++i;
  }
  auto aEvt = new TCanvas();
  histnEvents->Draw();
  aEvt->Print("PeriodMeanMult.pdf");
  auto c = new TCanvas();
  histPurityLambda->Draw();
  c->Print("PeriodQALambda.pdf");
  auto d = new TCanvas();
  histPurityAntiLambda->Draw();
  d->Print("PeriodQAAntiLambda.pdf");
  auto e = new TCanvas();
  histPurityXi->Draw();
  e->Print("PeriodQAXi.pdf");
  auto f = new TCanvas();
  histPurityAntiXi->Draw();
  f->Print("PeriodQAAntiXi.pdf");
}

void PeriodQA::ProcessSigmaQA(const char* prefix, const char* addon) {
  auto histPurityLambda = PeriodQAHist("histQALambda", "Purity #Lambda (%)");
  auto histPurityAntiLambda = PeriodQAHist("histQAAntiLambda",
                                           "Purity #bar{#Lambda} (%)");
  auto histPuritySigma = PeriodQAHist(
      "histQASigma", "Purity #Sigma^{0} #oplus #bar{#Sigma^{0}} (%)");

  auto histNProton = PeriodQAHist("histNProton", "p/Event");
  auto histNAntiProton = PeriodQAHist("histNAntiProton", "#bar{p}/Event");
  auto histNLambda = PeriodQAHist("histNLambda", "#Lambda/Event");
  auto histNAntiLambda = PeriodQAHist("histNAntiLambda", "#bar{#Lambda}/Event");
  auto histNSigma = PeriodQAHist("histNSigma",
                                 "(#Sigma^{0} #oplus #bar{#Sigma^{0}})/Event");

  int i = 0;
  for (auto it : fPeriods) {
    TString filename = fDirectory.Data();
    filename += "AnalysisResults_";
    filename += it;
    filename += ".root";
    ForgivingReader* reader = new ForgivingReader(filename.Data(), prefix,
                                                  addon);
    EventQA* evtQA = new EventQA();
    evtQA->SetEventCuts(reader->GetEventCuts());
    float nEvents = evtQA->GetNumberOfEvents();
    delete evtQA;
    if (nEvents < 1E-6) {
      delete reader;
      ++i;
      continue;
    }

    TrackQA* trkQA = new TrackQA();
    trkQA->SetTrackCuts(reader->GetTrackCuts());
    histNProton->SetBinContent(i + 1,
                               float(trkQA->GetNumberOfTracks()) / nEvents);
    delete trkQA;

    TrackQA* antiTrkQA = new TrackQA();
    antiTrkQA->SetAntiTrackCuts(reader->GetAntiTrackCuts());
    histNAntiProton->SetBinContent(
        i + 1, float(antiTrkQA->GetNumberOfTracks()) / nEvents);
    delete antiTrkQA;

    DecayQA* v0QA = new DecayQA("#Lambda", "p#pi");
    v0QA->SetDecayCuts(reader->Getv0Cuts());
    v0QA->SetCanvasDivisions(5, 2);
    v0QA->SetIMHistoScale(1.75, 0.8, 0.35);
    v0QA->SetRangesFitting(1.109, 1.121, 1.09, 1.15);
    v0QA->GetPeriodQASigma(1.112, 1.120, it.Data());
    float purity = v0QA->GetPurity();
    float purityErr = v0QA->GetPurityErr();
    if (purity > 1E-6) {
      histPurityLambda->SetBinContent(i + 1, purity * 100.f);
      histPurityLambda->SetBinError(i + 1, purityErr * 100.f);
    }
    histNLambda->SetBinContent(i + 1, v0QA->GetSignalCounts() / nEvents);
    delete v0QA;

    DecayQA* antiv0QA = new DecayQA("#bar{#Lambda}", "p#pi");
    antiv0QA->SetDecayCuts(reader->GetAntiv0Cuts());
    antiv0QA->SetCanvasDivisions(5, 2);
    antiv0QA->SetIMHistoScale(1.75, 0.8, 0.35);
    antiv0QA->SetRangesFitting(1.109, 1.121, 1.09, 1.15);
    antiv0QA->GetPeriodQASigma(1.112, 1.120, it.Data());
    purity = antiv0QA->GetPurity();
    purityErr = antiv0QA->GetPurityErr();
    if (purity > 1E-6) {
      histPurityAntiLambda->SetBinContent(i + 1, purity * 100.f);
      histPurityAntiLambda->SetBinError(i + 1, purityErr * 100.f);
    }
    histNAntiLambda->SetBinContent(i + 1,
                                   antiv0QA->GetSignalCounts() / nEvents);
    delete antiv0QA;

    DecayQA* sigma0QA = new DecayQA("#Sigma", "#Lambda#gamma");
    sigma0QA->SetDecayCuts(reader->GetOtherCuts("Sigma0Cuts"));
    sigma0QA->SetAntiDecayCuts(reader->GetOtherCuts("AntiSigma0Cuts"));
    sigma0QA->SetRangesFitting(1.187, 1.199, 1.167, 1.217);
    sigma0QA->GetPeriodQASigma0(0.003, it.Data());
    purity = sigma0QA->GetPurity();
    purityErr = sigma0QA->GetPurityErr();
    if (purity > 1E-6) {
      histPuritySigma->SetBinContent(i + 1, purity * 100.f);
      histPuritySigma->SetBinError(i + 1, purityErr * 100.f);
    }
    histNSigma->SetBinContent(i + 1, sigma0QA->GetSignalCounts() / nEvents);
    histNSigma->SetBinError(i + 1, sigma0QA->GetSignalCountsErr() / nEvents);

    delete sigma0QA;
    delete reader;
    ++i;
  }
  auto c = new TCanvas();
  histPurityLambda->Draw("PE");
  histPurityLambda->Fit("pol0");
  c->Print("PeriodQALambda.pdf");
  auto d = new TCanvas();
  histPurityAntiLambda->Draw("PE");
  histPurityAntiLambda->Fit("pol0");
  d->Print("PeriodQAAntiLambda.pdf");
  auto e = new TCanvas();
  histPuritySigma->Draw("PE");
  histPuritySigma->Fit("pol0");
  e->Print("PeriodQASigma.pdf");

  auto a1 = new TCanvas();
  histNProton->Draw();
  a1->Print("PeriodQANProton.pdf");
  auto b1 = new TCanvas();
  histNAntiProton->Draw();
  b1->Print("PeriodQANAntiProton.pdf");
  auto c1 = new TCanvas();
  histNLambda->Draw();
  c1->Print("PeriodQANLambda.pdf");
  auto d1 = new TCanvas();
  histNAntiLambda->Draw();
  d1->Print("PeriodQANAntiLambda.pdf");
  auto e1 = new TCanvas();
  histNSigma->Draw("PE");
  e1->Print("PeriodQANSigma.pdf");
}

TH1F* PeriodQA::PeriodQAHist(const char* name, const char* title) {
  auto histnew = new TH1F(name, Form("; ; %s", title), fPeriods.size(), 0,
                          fPeriods.size());
  histnew->SetMarkerStyle(20);
  int i = 1;
  for (auto it : fPeriods) {
    histnew->GetXaxis()->SetBinLabel(i++, it.Data());
  }
  histnew->LabelsOption("v", "X");
  return histnew;
}

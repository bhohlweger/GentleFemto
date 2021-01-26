#include "SidebandSigma.h"
#include "DreamPlot.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TROOT.h"

void SigmaEvalSidebands(TString InputDir, TString trigger) {
  gROOT->ProcessLine("gErrorIgnoreLevel = 3001");
  DreamPlot::SetStyle(false, true);

  std::vector<TString> varVec = { { "3 - 25", "3 - 50", "3 - 75", "3 - 100",
      "5 - 25", "5 - 50", "5 - 75", "5 - 100", "10 - 25", "10 - 50", "10 - 75",
      "10 - 100", "25 - 50", "25 - 75", "25 - 100" } };

  std::vector<int> colorVec = { { kBlue - 9, kBlue - 6, kBlue - 2, kBlue + 3,
      kRed - 9, kRed - 6, kRed - 2, kRed + 3, kGreen - 9, kGreen - 6, kGreen
          - 2, kGreen + 3, kOrange + 6, kOrange +1, kOrange - 6 } };

  std::vector<int> markerVec = { { 24, 24, 24, 24, 25, 25, 25, 25, 26, 26, 26,
      26, 27, 27, 27 } };

  auto histDummy = new TH1F("histDummy", "; #it{k}* (MeV/#it{c}); C(#it{k}*)",
                            100, 0, 500);

  auto c = new TCanvas("sideband", "sideband", 0, 0, 650, 550);
  histDummy->Draw();
  histDummy->GetXaxis()->SetNdivisions(505);
  histDummy->SetMinimum(0.85);
  histDummy->SetMaximum(1.7);

  const float leftX = 0.45;
  auto leg = new TLegend(leftX, 0.45, 0.925, 0.76);
  leg->SetTextSize(gStyle->GetTextSize() * .9);
  leg->SetNColumns(3);

  for (int i = 1; i < 15; ++i) {
    if ( i == 4|| i == 8 || i ==12) continue;
    auto side = new SidebandSigma();
    side->SetRebin(10);
    side->SetSideBandFile(InputDir.Data(), trigger.Data(), Form("%i", i));
    side->SetNormalizationRange(250, 400);
    side->SideBandCFs();
    auto SBmerge = side->GetSideBandGraph(5);
    DreamPlot::SetStyleGraph(SBmerge, markerVec[i - 1], colorVec[i - 1], 0.75);
    SBmerge->Draw("pez");
    leg->AddEntry(SBmerge, varVec[i - 1].Data(), "PE");
  }

  TLatex BeamText;
  BeamText.SetTextSize(gStyle->GetTextSize() * .9);
  BeamText.SetNDC(kTRUE);
  BeamText.DrawLatex(leftX + 0.01, 0.9, "ALICE this thesis");
  BeamText.DrawLatex(leftX + 0.01, 0.84,
                     "pp #sqrt{#it{s}} = 13 TeV (High-mult.)");
  BeamText.DrawLatex(
      leftX + 0.01,
      0.78,
      "p#minus#kern[-0.65]{ }(#Lambda#gamma) #oplus #bar{p}#minus#kern[-0.4]{ }(#bar{#Lambda}#gamma)");

  leg->Draw("same");

  c->Print("sideband.pdf");
}

int main(int argc, char* argv[]) {
  SigmaEvalSidebands(argv[1], argv[2]);

  return 1;
}

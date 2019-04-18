#include "SidebandSigma.h"
#include "DreamPlot.h"
#include "TCanvas.h"
#include "TLegend.h"

void SigmaEvalSidebands(TString InputDir, TString trigger) {
  gROOT->ProcessLine("gErrorIgnoreLevel = 3001");
  DreamPlot::SetStyle(false, true);
  auto c = new TCanvas();
  std::vector<TH1F*> histVec;

  std::vector<TString> varVec = { { "[10, 75]", "[10, 100]", "[25, 50]",
      "[25, 75]", "[25, 100]", "[5, 25]", "[5, 50]", "[5, 75]", "[5, 100]" } };

  for (int i = 33; i < 42; ++i) {
    auto side = new SidebandSigma();
    side->SetRebin(10);
    side->SetSideBandFile(InputDir.Data(), trigger.Data(), Form("%i", i));
    side->SetNormalizationRange(300, 500);
    side->SideBandCFs();
    auto SBmerge = side->GetSideBands(5);
    DreamPlot::SetStyleHisto(SBmerge, 20, i);
    if (i == 33) {
      SBmerge->Draw();
      SBmerge->SetTitle("; #it{k}* (MeV/#it{c}); C(#it{k}*)");
    } else {
      SBmerge->Draw("same");
    }
    histVec.push_back(SBmerge);
    delete side;
  }
  auto leg = new TLegend(0.6, 0.6, 0.85, 0.85);
  leg->SetNColumns(2);
  int i = 0;
  for(auto it : histVec) {
    leg->AddEntry(it, varVec[i++].Data(), "PE");
  }
  leg->Draw("same");

  c->Print("sideband.pdf");
}

int main(int argc, char* argv[]) {
  SigmaEvalSidebands(argv[1], argv[2]);

  return 1;
}

#include "SidebandSigma.h"
#include "DreamPlot.h"
#include "TCanvas.h"
#include "TLegend.h"

void SigmaEvalSidebands(TString InputDir, TString trigger) {
  gROOT->ProcessLine("gErrorIgnoreLevel = 3001");
  DreamPlot::SetStyle(false, true);
  auto c = new TCanvas();
  std::vector<TH1F*> histVec;

  std::vector<TString> varVec = { { "[5, 10]", "[3, 10]", "[3, 50]",
      "[3, 100]", "[2, 10]", "[2, 50]", "[2, 100]", "[5, 10]", "[5, 50]",
      "[5, 100]", "[10, 10]", "[10, 50]", "[25, 10]", "[25, 50]", "[25, 100]",
      "[50, 10]", "[50, 50]", "[50, 100]" } };

  for (int i = 0; i < 18; ++i) {
    auto side = new SidebandSigma();
    side->SetRebin(10);
    side->SetSideBandFile(InputDir.Data(), trigger.Data(), Form("%i", i));
    side->SetNormalizationRange(250, 400);
    side->SideBandCFs();
    auto SBmerge = side->GetSideBands(5);
    DreamPlot::SetStyleHisto(SBmerge, 20, i+1);
    if (i == 0) {
      SBmerge->Draw();
      SBmerge->SetTitle("; #it{k}* (MeV/#it{c}); C(#it{k}*)");
      SBmerge->SetMinimum(0.8);
      SBmerge->SetMaximum(2);
      SBmerge->GetXaxis()->SetRangeUser(0, 500);
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

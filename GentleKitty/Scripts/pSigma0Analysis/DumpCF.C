#include "CATSInputSigma0.h"
#include "SidebandSigma.h"

int main(int argc, char *argv[]) {
  TString InputDir = argv[1];
  TString trigger = argv[2];

  TString dataHistName = "hCk_ReweightedpSigma0MeV_0";
  TString dataGrName = "Graph_from_hCk_ReweightedpSigma0_0MeV";
  const int nSidebandHist = 5;
  const int rebinSideband = 10;

  auto CATSinput = new CATSInputSigma0();
  auto appendix = TString::Format("%i", 0);
  CATSinput->ReadSigma0CorrelationFile(InputDir.Data(), trigger.Data(),
                                       appendix.Data());
  CATSinput->ObtainCFs(10, 250, 400);
  auto dataGr = CATSinput->GetCFGr("pSigma0", dataGrName.Data());

  auto side = new SidebandSigma();
  side->SetRebin(rebinSideband);
  side->SetSideBandFile(InputDir.Data(), trigger.Data(), appendix.Data());
  side->SetNormalizationRange(250, 400);
  side->SideBandCFs();
  auto sbGr = side->GetSideBandGraph(nSidebandHist);

  auto file = new TFile("cfoutput.root", "RECREATE");
  dataGr->Write("pSigma0");
  sbGr->Write("pSB");
  file->Close();
}

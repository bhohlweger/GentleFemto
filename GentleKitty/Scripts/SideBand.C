#include "SideBandFit.h"
#include "TF1.h"
#include "TCanvas.h"

int main(int argc, char *argv[]) {
  TString WorkDir = TString::Format("%s", argv[1]);
  const char* prefix = argv[2];
  const char* suffixUp = argv[3];
  const char* suffixDown = argv[4];
  SideBandFit* side = new SideBandFit();
  side->SetRebin(5);
  side->SetSideBandFile(WorkDir.Data(), prefix, suffixUp, suffixDown);

  side->SetNormalizationRange(400, 600);
  side->SideBandCFs(true);
  TH1F* fitme = side->GetSideBands(5);
  double SideBandPars[4];
  side->FitSideBands(fitme, SideBandPars);
//	side->WriteOutput("/home/hohlweger/cernbox/pPb/Sidebands");
  TF1* sideFit = (TF1*) fitme->FindObject("SideBandFit");
  TCanvas* c1 = new TCanvas("c1", "c1");
  fitme->Draw();
  sideFit->Draw("SAME");
  c1->SaveAs(
      TString::Format("%s/fitsideband%s_%s.png", WorkDir.Data(), suffixUp,
                      suffixDown));
  TFile* outfile = TFile::Open(TString::Format("%s/out.root", WorkDir.Data()),
                               "update");
  outfile->cd();
  c1->Write();
  outfile->Close();
  side->WriteOutput(
      TString::Format("%s/Var%s_%s.root", WorkDir.Data(), suffixUp,
                      suffixDown));
//    }
  return 0;
}

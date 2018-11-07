#include "SideBandFit.h"
#include "TF1.h"
#include "TCanvas.h"

int main(int argc, char *argv[]) {
	SideBandFit* side = new SideBandFit();
	side->SetRebin(5);
	side->SetSideBandFile("/home/hohlweger/cernbox/pPb/Sidebands","42","43");

	side->SetNormalizationRange(400,600);
	side->SideBandCFs();
	TH1F* fitme = side->GetSideBands(5);
	double SideBandPars[4];
	side->FitSideBands(fitme,SideBandPars);
//	side->WriteOutput("/home/hohlweger/cernbox/pPb/Sidebands");
    TF1* sideFit = (TF1*)fitme->FindObject("SideBandFit");
    TCanvas* c1 = new TCanvas("c1","c1");
    fitme->Draw();
    sideFit->Draw("SAME");
    c1->SaveAs("/home/hohlweger/cernbox/pPb/Sidebands/fitsideband.png");
    TFile* outfile = TFile::Open("/home/hohlweger/cernbox/pPb/Sidebands/out.root","RECREATE");
    outfile->cd();
    c1->Write();
    outfile->Close();
//    for (auto it = 0; it< 4; ++it) {
//    	std::cout << "SideBand Par " << it << " value " << SideBandPars[it] << std::endl;
//    }
	return 0;
}

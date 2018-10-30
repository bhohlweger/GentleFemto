#include "DLM_Source.h"
#include "DLM_Potentials.h"
#include "DLM_CkModels.h"
#include "CATS.h"
#include "DLM_CkDecomposition.h"
#include "DLM_Fitters.h"
#include "TRandom3.h"
#include "TH2F.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TGraph.h"
#include "TF1.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TMath.h"
#include "TCanvas.h"
#include <iostream>
#include "stdlib.h"

void DrawSomeThings(const char* outfile, int radius) {
	double PotPars1S0[10] = { 0, 0, NN_AV18, v18_Coupled3P2, 1, 1, 1, 0, 0, 0 };
	double PotPars3P0[10] = { 0, 0, NN_AV18, v18_Coupled3P2, 1, 1, 1, 1, 1, 0 };
	double PotPars3P1[10] = { 0, 0, NN_AV18, v18_Coupled3P2, 1, 1, 1, 1, 1, 1 };
	double PotPars3P2[10] = { 0, 0, NN_AV18, v18_Coupled3P2, 1, 1, 1, 1, 1, 2 };

	const double Weight1S0 = 3. / 12.;
	const double Weight3P0 = 1. / 12.;
	const double Weight3P1 = 3. / 12.;
	const double Weight3P2 = 5. / 12.;

	const double Mass_p = 938.272;
	const double Mass_L = 1115.683;

	CATS AB_pp;
	double Pars_pp[6] = { 0, 0, 0, radius * 1.2, radius / 1.2,
			0.5 };
	AB_pp.SetAnaSource(GaussSource, Pars_pp);
	AB_pp.SetUseAnalyticSource(true);
	AB_pp.SetThetaDependentSource(false);

	AB_pp.SetExcludeFailedBins(false);
	AB_pp.SetMomBins(100, 0, 400);

	AB_pp.SetNumChannels(4);
	AB_pp.SetNumPW(0, 2);
	AB_pp.SetNumPW(1, 2);
	AB_pp.SetNumPW(2, 2);
	AB_pp.SetNumPW(3, 2);
	AB_pp.SetSpin(0, 0);
	AB_pp.SetSpin(1, 1);
	AB_pp.SetSpin(2, 1);
	AB_pp.SetSpin(3, 1);
	AB_pp.SetChannelWeight(0, Weight1S0);
	AB_pp.SetChannelWeight(1, Weight3P0);
	AB_pp.SetChannelWeight(2, Weight3P1);
	AB_pp.SetChannelWeight(3, Weight3P2);

	AB_pp.SetQ1Q2(1);
	AB_pp.SetPdgId(2212, 2212);
	AB_pp.SetRedMass(0.5 * Mass_p);

	AB_pp.SetShortRangePotential(0, 0, fDlmPot, PotPars1S0);
	AB_pp.SetShortRangePotential(1, 1, fDlmPot, PotPars3P0);
	AB_pp.SetShortRangePotential(2, 1, fDlmPot, PotPars3P1);
	AB_pp.SetShortRangePotential(3, 1, fDlmPot, PotPars3P2);

	AB_pp.KillTheCat();

	TGraph  cats_pp;
	for(int iBin = 0; iBin < AB_pp.GetNumMomBins(); ++iBin) {
		cats_pp.SetPoint(iBin+1,AB_pp.GetMomentum(iBin),AB_pp.GetCorrFun(iBin));
	}
	cats_pp.SetName("CatsCF");
	TCanvas* cCats = new TCanvas("cCATS","cCATS",1);
	cCats->cd();
	cats_pp.Draw("AP");

	TFile* Output = TFile::Open(Form("%smyCats.root",outfile),"RECREATE");
	cats_pp.Write();
	Output->Close();

	//	DLM_Ck* Ck_pp = new DLM_Ck(NumSourcePars, 0, AB_pp);
//
//	DLM_CkDecomposition CkDec_pp("pp", 0, *Ck_pp, hSigma_pp_MeV);

}


int main(int argc, char *argv[]) {
	DrawSomeThings(argv[1],atof(argv[2]));
	return 0;
}

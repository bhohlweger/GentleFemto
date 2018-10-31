#include "ForBernie.h"
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

void GetXiForRadius(TString InputDir, TString tmpOutputDir,
		double GaussSourceSize, const int tOut, TFile *file,
		const char* GraphName, bool saveCoulomb) {
	std::cout << GaussSourceSize << std::endl;
	const int binwidth = 20;
	const unsigned NumMomBins_pXim = 13;
	double kMinXiP = 8.;
	const double kMin_pXim = kMinXiP;
	const double kMax_pXim = kMin_pXim + binwidth * NumMomBins_pXim;

	std::cout << "kMinXiP: " << kMin_pXim << std::endl;
	std::cout << "kMax_pXim: " << kMax_pXim << std::endl;
	std::cout << "binwidth: " << binwidth << std::endl;
	std::cout << "NumMomBins_pXim: " << NumMomBins_pXim << std::endl;

	double FemtoRegion_pXim[3][2];
	FemtoRegion_pXim[0][0] = kMin_pXim;
	FemtoRegion_pXim[0][1] = 180;
	FemtoRegion_pXim[1][0] = kMin_pXim;
	FemtoRegion_pXim[1][1] = 220;
	FemtoRegion_pXim[2][0] = kMin_pXim;
	FemtoRegion_pXim[2][1] = 260;

	double BlRegion[3][2];
	BlRegion[0][0] = 320;
	BlRegion[0][1] = 480;
	BlRegion[1][0] = 300;
	BlRegion[1][1] = 500;
	BlRegion[2][0] = 300;
	BlRegion[2][1] = 540;

	double PurityProton = 0.984265; //pPb 5 TeV
	double PurityXi = 0.88; //new cuts

	double pp_f0 = 0.862814;
	double pp_f1 = 0.09603;
	double pL_f0 = 0.521433; //fraction of primary Lambdas
	double pL_f1 = 0.173811; //fraction of Sigma0
	double pL_f2 = 0.152378; //fractions of Xi0/m

	double ProtonPrim = pp_f0;
	double arrayPercLamProton[3] = { pp_f1 / (1. - pp_f0) * 0.8, pp_f1
			/ (1. - pp_f0), pp_f1 / (1. - pp_f0) * 1.2 }; //+/- 20%

	const unsigned NumChannels_p = 4;
	double** Purities_p = new double*[3];
	double** Fraction_p = new double*[3];
	for (unsigned uVar = 0; uVar < 3; uVar++) {
		Purities_p[uVar] = new double[NumChannels_p];
		Fraction_p[uVar] = new double[NumChannels_p];

		Purities_p[uVar][0] = PurityProton;
		Purities_p[uVar][1] = PurityProton;
		Purities_p[uVar][2] = PurityProton;
		Purities_p[uVar][3] = 1. - PurityProton;

		Fraction_p[uVar][0] = ProtonPrim;
		Fraction_p[uVar][1] = (1. - ProtonPrim) * (arrayPercLamProton[uVar]);
		Fraction_p[uVar][2] = (1. - ProtonPrim)
				* (1. - arrayPercLamProton[uVar]);
		Fraction_p[uVar][3] = 1.;
	}

	//ratio Xi-(1530) to Xi-
	const double Xim1530_to_Xim = 0.32 * (1. / 3.);
	//ratio Xi0(1530) to Xi0 (n=neutral)
	const double Xin1530_to_Xim = 0.32 * (2. / 3.);
	const double Omegam_to_Xim = 0.1;
	const double OmegamXim_BR = 0.086;

	const unsigned NumChannels_Xim = 5;
	double** Purities_Xim = new double*[3];
	double** Fraction_Xim = new double*[3];
	for (unsigned uVar = 0; uVar < 3; uVar++) {
		Purities_Xim[uVar] = new double[NumChannels_Xim];
		Fraction_Xim[uVar] = new double[NumChannels_Xim];

		Purities_Xim[uVar][0] = PurityXi;
		Purities_Xim[uVar][1] = PurityXi;
		Purities_Xim[uVar][2] = PurityXi;
		Purities_Xim[uVar][3] = PurityXi;
		Purities_Xim[uVar][4] = 1. - PurityXi;

		//the ratios that we have for Xis are referred to the total number of Xi particles (which already include all contributions)
		//hence Xi1530_to_Xi indeed is simply the number of Xis that stem from a Xi1530
		Fraction_Xim[uVar][0] = 1. - Xim1530_to_Xim - Xin1530_to_Xim
				- Omegam_to_Xim * OmegamXim_BR;
		Fraction_Xim[uVar][1] = Xim1530_to_Xim;
		Fraction_Xim[uVar][2] = Xin1530_to_Xim;
		Fraction_Xim[uVar][3] = Omegam_to_Xim * OmegamXim_BR;
		Fraction_Xim[uVar][4] = 1.;
	}

	TString CalibBaseDir = "~/cernbox/SystematicsAndCalib/pPbRun2_MB";
	TString ResMatrixFileName = TString::Format(
			"%s/run2_decay_matrices_old.root", CalibBaseDir.Data());
	TString SigmaMatrixFileName = TString::Format("%s/Sample3_MeV_compact.root",
			CalibBaseDir.Data());

	CATSInput *CATSinput = new CATSInput();
	CATSinput->SetCalibBaseDir(CalibBaseDir.Data());
	CATSinput->SetMomResFileName("run2_decay_matrices_old.root");
	CATSinput->ReadResFile();
	CATSinput->SetSigmaFileName("Sample3_MeV_compact.root");
	CATSinput->ReadSigmaFile();

	CATSinput->ReadCorrelationFile(InputDir.Data());
	//  int vCutID;//which data file (cut combination) should you take. 0 = default
	//int vSource;//which source we use, see above
	int vFemReg_pXim;  //which femto region we use for pXim (1 = default)
	int vBlReg;  //which baseline region to use (1 = default)

	vFemReg_pXim = 1;
	vBlReg = 1;
	if (tOut < 11 && tOut > 13) {
		std::cout << tOut << " not supported \n";
		return;
	}
	double Pars_pXi[6] = { 0, 0, 0, GaussSourceSize * 1.2, GaussSourceSize
					/ 1.2, 0.5 };
	TidyCats* tidy = new TidyCats();
	CATS AB_pXim;
	tidy->GetCatsProtonXiMinus(&AB_pXim, GaussSourceSize, Pars_pXi,
			NumMomBins_pXim, kMin_pXim, kMax_pXim, true, 13);
	AB_pXim.KillTheCat();

	double Pars_pXim1530[6] = { 0, 0, 0, GaussSourceSize * 1.2, GaussSourceSize
			/ 1.2, 0.5 };
	CATS AB_pXim1530;
	tidy->GetCatsProtonXiMinus1530(&AB_pXim1530, GaussSourceSize, Pars_pXim1530,
			NumMomBins_pXim, kMin_pXim, kMax_pXim);
	AB_pXim1530.KillTheCat();

	TString HistpXiName = "hCk_ReweightedpXiMeV_0";
	TH1F* OliHisto_pXim = CATSinput->GetCF("pXi", HistpXiName.Data());
	if (!OliHisto_pXim)
		std::cout << HistpXiName.Data() << " pXi Missing" << std::endl;

	TH1F *OliHisto_pXimFornSigma = nullptr;
	OliHisto_pXimFornSigma = (TH1F*) OliHisto_pXim->Clone("pXiForNSigma");
	CATSinput->AddSystematics("C2totalsysPXi.root", OliHisto_pXim);

	const unsigned NumSourcePars = 1;
	DLM_Ck* Ck_pXim = new DLM_Ck(NumSourcePars, 0, AB_pXim);
	DLM_Ck* Ck_pXim1530 = new DLM_Ck(NumSourcePars, 0, AB_pXim1530);
	Ck_pXim->Update();
	Ck_pXim1530->Update();
	if (!CATSinput->GetSigmaFile(3)) {
		std::cout << "No Sigma file 3 \n";
		return;
	}
	DLM_CkDecomposition CkDec_pXim("pXim", 3, *Ck_pXim,
			CATSinput->GetSigmaFile(3));
	DLM_CkDecomposition CkDec_pXim1530("pXim1530", 0, *Ck_pXim1530, NULL);
	int vFrac_pp_pL = 1;
	const double lam_pXim = Purities_p[vFrac_pp_pL][0]
			* Fraction_p[vFrac_pp_pL][0] * Purities_Xim[0][0]
			* Fraction_Xim[0][0];
	const double lam_pXim_pXim1530 = Purities_p[vFrac_pp_pL][0]
			* Fraction_p[vFrac_pp_pL][0] * Purities_Xim[0][1]
			* Fraction_Xim[0][1];
	const double lam_pXim_fake = Purities_p[vFrac_pp_pL][3] * Purities_Xim[0][0]
			+ Purities_p[vFrac_pp_pL][0] * Purities_Xim[0][4]
			+ Purities_p[vFrac_pp_pL][3] * Purities_Xim[0][4];

	printf("lam_pXim = %.3f\n", lam_pXim);
	printf("lam_pXim_pXim1530 = %.3f\n", lam_pXim_pXim1530);
	printf("lam_pXim_fake = %.3f\n", lam_pXim_fake);
	printf("\n");

	CkDec_pXim.AddContribution(0, lam_pXim_pXim1530,
			DLM_CkDecomposition::cFeedDown, &CkDec_pXim1530,
			CATSinput->GetResFile(3));		//from Xi-(1530)
	CkDec_pXim.AddContribution(1,
			1. - lam_pXim - lam_pXim_pXim1530 - lam_pXim_fake,
			DLM_CkDecomposition::cFeedDown);		//other feed-down (flat)
	CkDec_pXim.AddContribution(2, lam_pXim_fake, DLM_CkDecomposition::cFake);
	CkDec_pXim.Update();
	DLM_Fitter1* fitter = new DLM_Fitter1(1);
//	fitter->SetOutputDir(tmpOutputDir.Data());
	fitter->SetSystem(0, *OliHisto_pXim, 1, CkDec_pXim,
			FemtoRegion_pXim[vFemReg_pXim][0],
			FemtoRegion_pXim[vFemReg_pXim][1], BlRegion[vBlReg][0],
			BlRegion[vBlReg][1]);
//	fitter->SetSeparateBL(0, true);	 	//Separate BL
	fitter->SetSeparateBL(0, false); 	//Simultaneous BL
	fitter->AddSameSource("pXim1530", "pXim", 1);
	//Global Fit default
	//baseline
//	fitter->FixParameter("pXim", DLM_Fitter1::p_a, 1.06826);
//	fitter->FixParameter("pXim", DLM_Fitter1::p_b, 1.2624e-13);
//	fitter->FixParameter("pXim", DLM_Fitter1::p_c, 0);
//	//gaussian radius
//	fitter->FixParameter("pXim", DLM_Fitter1::p_sor0, GaussSourceSize);
//	//normalization
//	fitter->FixParameter("pXim", DLM_Fitter1::p_Cl, -0.933815);

	//Standalone Fit default
	//baseline
//	fitter->FixParameter("pXim", DLM_Fitter1::p_a, 1.05447);
//	fitter->FixParameter("pXim", DLM_Fitter1::p_b, 2.10185e-07);
//	fitter->FixParameter("pXim", DLM_Fitter1::p_c, 0);
//	//gaussian radius
//	fitter->FixParameter("pXim", DLM_Fitter1::p_sor0, GaussSourceSize);
//	//normalization
//	fitter->FixParameter("pXim", DLM_Fitter1::p_Cl, -0.933603);

	//Fit BL & Normalization
	fitter->FixParameter("pXim", DLM_Fitter1::p_c, 0);

	fitter->FixParameter("pXim", DLM_Fitter1::p_a, 1.);
	fitter->FixParameter("pXim", DLM_Fitter1::p_b, 0);
//	fitter->SetParameter("pXim", DLM_Fitter1::p_a, 1.0, 0.7, 1.3);
//	fitter->SetParameter("pXim", DLM_Fitter1::p_b, 1e-4, 0, 2e-3);

//	fitter->SetParameter("pXim", DLM_Fitter1::p_Cl, -0.9, -1.2, -0.8);
	fitter->FixParameter("pXim", DLM_Fitter1::p_Cl, -1.);

	fitter->FixParameter("pXim", DLM_Fitter1::p_sor0, GaussSourceSize);

	double StartPar = 1;
	fitter->GoBabyGo();

	double p_a_strong = fitter->GetParameter("pXim", DLM_Fitter1::p_a);
	double p_b_strong = fitter->GetParameter("pXim", DLM_Fitter1::p_b);
	double Cl_strong = fitter->GetParameter("pXim", DLM_Fitter1::p_c);

	std::cout << "AFTER THE FIT \n";
	std::cout << "p_a: " << fitter->GetParameter("pXim", DLM_Fitter1::p_a)
			<< std::endl;
	std::cout << "p_b: " << fitter->GetParameter("pXim", DLM_Fitter1::p_b)
			<< std::endl;
	std::cout << "p_c: " << fitter->GetParameter("pXim", DLM_Fitter1::p_c)
			<< std::endl;
	std::cout << "p_sor0:" << fitter->GetParameter("pXim", DLM_Fitter1::p_sor0)
			<< std::endl;
	std::cout << "p_Cl: " << fitter->GetParameter("pXim", DLM_Fitter1::p_Cl)
			<< std::endl;

	TGraph FitResult_pXim;
	FitResult_pXim.SetName(TString::Format("pXimGraph%s", GraphName));
	fitter->GetFitGraph(0, FitResult_pXim);

//	CoulombStrong.SetNa 1221 total, 1053 done, 99 error, 69 active, 0 waiting me(
//			TString::Format("FitResult_pXim_%.2f", GaussSourceSize));
//	fitter->GetFitGraph(0, CoulombStrong);

	double Chi2 = fitter->GetChi2();
	unsigned NDF = fitter->GetNdf();

	double Chi2_pXim = 0;
	unsigned EffNumBins_pXim = 0;
	double Chi2_pXim_exclPeak = 0;
	unsigned EffNumBins_pXim_exclPeak = 0;
	double Chi2_pXim_kSm60 = 0;
	unsigned EffNumBins_pXim_kSm60 = 0;
	int maxkStarBin = OliHisto_pXimFornSigma->FindBin(200);
	for (unsigned uBin = 0; uBin < maxkStarBin; uBin++) {

		double mom = AB_pXim.GetMomentum(uBin);
		//double dataX;
		double dataY;
		double dataErr;
		double theoryX;
		double theoryY;

		if (mom > FemtoRegion_pXim[vFemReg_pXim][1])
			continue;

		FitResult_pXim.GetPoint(uBin, theoryX, theoryY);
		if (mom != theoryX) {
			std::cout << mom << '\t' << theoryX << std::endl;
			printf("  PROBLEM pXi!\n");
		}
		dataY = OliHisto_pXimFornSigma->GetBinContent(uBin + 1);
		dataErr = OliHisto_pXimFornSigma->GetBinError(uBin + 1);
		Chi2_pXim += (dataY - theoryY) * (dataY - theoryY)
				/ (dataErr * dataErr);
		EffNumBins_pXim++;
		if (mom < 60) {
			Chi2_pXim_kSm60 += (dataY - theoryY) * (dataY - theoryY)
					/ (dataErr * dataErr);
			EffNumBins_pXim_kSm60++;
		}
		if (!(mom > 60 && mom < 100)) {
			Chi2_pXim_exclPeak += (dataY - theoryY) * (dataY - theoryY)
					/ (dataErr * dataErr);
			EffNumBins_pXim_exclPeak++;
		}
	}
	EffNumBins_pXim--;
	EffNumBins_pXim_kSm60--;
	EffNumBins_pXim_exclPeak--;

	printf("AB_pXim[0] = %.2f\n", AB_pXim.GetCorrFun(0));
	AB_pXim.RemoveShortRangePotential(0, 0);
	AB_pXim.RemoveShortRangePotential(1, 0);
	AB_pXim.RemoveShortRangePotential(2, 0);
	AB_pXim.RemoveShortRangePotential(3, 0);
	AB_pXim.KillTheCat(CATS::kPotentialChanged);
	//AB_pXim.KillTheCat(CATS::kAllChanged);
	printf("NEW AB_pXim[0] = %.2f\n", AB_pXim.GetCorrFun(0));

	CkDec_pXim.Update(true);
	fitter->GoBabyGo();
	TGraph FitResult_pXim_COULOMB;
	FitResult_pXim_COULOMB.SetName(
			TString::Format("pXimGraph%s_COULOMB", GraphName));
	fitter->GetFitGraph(0, FitResult_pXim_COULOMB);

//	Coulomb.SetName(
//			TString::Format("FitResult_pXim_COULOMB_%.2f", GaussSourceSize));
//	fitter->GetFitGraph(0, Coulomb);

	double p_a_coulomb = fitter->GetParameter("pXim", DLM_Fitter1::p_a);
	double p_b_coulomb = fitter->GetParameter("pXim", DLM_Fitter1::p_b);
	double Cl_coulomb = fitter->GetParameter("pXim", DLM_Fitter1::p_c);

	double Chi2_pXim_COULOMB = 0;
	unsigned EffNumBins_pXim_COULOMB = 0;
	double Chi2_pXim_COULOMB_exclPeak = 0;
	unsigned EffNumBins_pXim_COULOMB_exclPeak = 0;
	double Chi2_pXim_COULOMB_kSm60 = 0;
	unsigned EffNumBins_pXim_COULOMB_kSm60 = 0;

	for (unsigned uBin = 0; uBin < maxkStarBin; uBin++) {
		//    for(unsigned uBin=0; uBin<8; uBin++){

		double mom = AB_pXim.GetMomentum(uBin);
		double dataY;
		double dataErr;
		double theoryX;
		double theoryY;

		if (mom > FemtoRegion_pXim[vFemReg_pXim][1])
			continue;
		FitResult_pXim_COULOMB.GetPoint(uBin, theoryX, theoryY);
		if (mom != theoryX) {
			std::cout << mom << '\t' << theoryX << std::endl;
			printf("  PROBLEM pXi!\n");
		}
		dataY = OliHisto_pXimFornSigma->GetBinContent(uBin + 1);
		dataErr = OliHisto_pXimFornSigma->GetBinError(uBin + 1);

		Chi2_pXim_COULOMB += (dataY - theoryY) * (dataY - theoryY)
				/ (dataErr * dataErr);
		EffNumBins_pXim_COULOMB++;
		if (mom < 60) {
			Chi2_pXim_COULOMB_kSm60 += (dataY - theoryY) * (dataY - theoryY)
					/ (dataErr * dataErr);
			EffNumBins_pXim_COULOMB_kSm60++;
		}
		if (!(mom > 60 && mom < 100)) {
			Chi2_pXim_COULOMB_exclPeak += (dataY - theoryY) * (dataY - theoryY)
					/ (dataErr * dataErr);
			EffNumBins_pXim_COULOMB_exclPeak++;
		}
	}
	EffNumBins_pXim_COULOMB--;
	EffNumBins_pXim_COULOMB_kSm60--;
	EffNumBins_pXim_COULOMB_exclPeak--;

	double pvalXi = TMath::Prob(Chi2_pXim, round(EffNumBins_pXim));
	double nSigmaXi = TMath::Sqrt(2) * TMath::ErfcInverse(pvalXi);

	double pvalXi_kSm60 = TMath::Prob(Chi2_pXim_kSm60,
			round(EffNumBins_pXim_kSm60));
	double nSigmaXi_kSm60 = TMath::Sqrt(2) * TMath::ErfcInverse(pvalXi_kSm60);

	double pvalXi_exclPeak = TMath::Prob(Chi2_pXim_exclPeak,
			round(EffNumBins_pXim_exclPeak));
	double nSigmaXi_exclPeak = TMath::Sqrt(2)
			* TMath::ErfcInverse(pvalXi_exclPeak);

	double pvalXiCoulomb = TMath::Prob(Chi2_pXim_COULOMB,
			round(EffNumBins_pXim_COULOMB));
	double nSigmaXiCoulomb = TMath::Sqrt(2) * TMath::ErfcInverse(pvalXiCoulomb);

	double pvalXiCoulomb_kSm60 = TMath::Prob(Chi2_pXim_COULOMB_kSm60,
			round(EffNumBins_pXim_COULOMB_kSm60));
	double nSigmaXiCoulomb_kSm60 = TMath::Sqrt(2)
			* TMath::ErfcInverse(pvalXiCoulomb_kSm60);

	double pvalXiCoulomb_exclPeak = TMath::Prob(Chi2_pXim_COULOMB_exclPeak,
			round(EffNumBins_pXim_COULOMB_exclPeak));
	double nSigmaXiCoulomb_exclPeak = TMath::Sqrt(2)
			* TMath::ErfcInverse(pvalXiCoulomb_exclPeak);

	std::cout << "Default \n";
	std::cout << "Coulomb + Strong - pValue: " << pvalXi << "  sigma:  "
			<< nSigmaXi << "  NDF:  " << EffNumBins_pXim << "  chi2:  "
			<< Chi2_pXim << std::endl;
	std::cout << "Coulomb - pValue: " << pvalXiCoulomb << "  sigma:  "
			<< nSigmaXiCoulomb << "  NDF:  " << EffNumBins_pXim_COULOMB
			<< "  chi2:  " << Chi2_pXim_COULOMB << std::endl;

	std::cout << "k < 60 MeV \n";
	std::cout << "Coulomb + Strong - pValue: " << pvalXi_kSm60 << "  sigma:  "
			<< nSigmaXi_kSm60 << "  NDF:  " << EffNumBins_pXim_kSm60
			<< "  chi2:  " << Chi2_pXim_kSm60 << std::endl;
	std::cout << "Coulomb - pValue: " << pvalXiCoulomb_kSm60 << "  sigma:  "
			<< nSigmaXiCoulomb_kSm60 << "  NDF:  "
			<< EffNumBins_pXim_COULOMB_kSm60 << "  chi2:  "
			<< Chi2_pXim_COULOMB_kSm60 << std::endl;

	std::cout << "Excluding the Peak \n";
	std::cout << "Coulomb + Strong - pValue: " << pvalXi_exclPeak
			<< "  sigma:  " << nSigmaXi_exclPeak << "  NDF:  "
			<< EffNumBins_pXim_exclPeak << "  chi2:  " << Chi2_pXim_exclPeak
			<< std::endl;
	std::cout << "Coulomb - pValue: " << pvalXiCoulomb_exclPeak << "  sigma:  "
			<< nSigmaXiCoulomb_exclPeak << "  NDF:  "
			<< EffNumBins_pXim_COULOMB_exclPeak << "  chi2:  "
			<< Chi2_pXim_COULOMB_exclPeak << std::endl;

	TPaveText* info4 = new TPaveText(0.2, 0.68, 0.9, 0.95, "blNDC"); //lbrt
	info4->SetName("info4");
	info4->SetBorderSize(1);
	info4->SetTextSize(0.04);
	info4->SetFillColor(kWhite);
	info4->SetTextFont(22);
	TString SOURCE_NAME = "Gauss";
	double Yoffset = 1.2;

	info4->AddText(
			TString::Format("R(%s)=%.3f#pm%.3f", SOURCE_NAME.Data(),
					fitter->GetParameter("pXim", DLM_Fitter1::p_sor0),
					fitter->GetParError("pXim", DLM_Fitter1::p_sor0)));
	info4->AddText(
			TString::Format("C_{l}=%.3f#pm%.3f",
					fitter->GetParameter("pXim", DLM_Fitter1::p_Cl),
					fitter->GetParError("pXim", DLM_Fitter1::p_Cl)));
	info4->AddText(
			TString::Format(
					"Local #chi^{2}/ndf=%.0f/%u=%.2f, p_{val}=%.5f, n#sigma=%.3f",
					Chi2_pXim, EffNumBins_pXim,
					Chi2_pXim / double(EffNumBins_pXim), pvalXi, nSigmaXi));
	info4->AddText(
			TString::Format(
					"Coulomb #chi^{2}/ndf=%.0f/%u=%.2f, p_{val}=%.5f, n#sigma=%.3f",
					Chi2_pXim_COULOMB, EffNumBins_pXim,
					Chi2_pXim_COULOMB / double(EffNumBins_pXim), pvalXiCoulomb,
					nSigmaXiCoulomb));
	TH1F* hAxis_pXim = new TH1F("hAxis_pXim", "hAxis_pXim", 600, 0, 600);
	hAxis_pXim->SetStats(false);
	hAxis_pXim->SetTitle("");
	hAxis_pXim->GetXaxis()->SetLabelSize(0.065);
	hAxis_pXim->GetXaxis()->CenterTitle();
	hAxis_pXim->GetXaxis()->SetTitleOffset(1.35);
	hAxis_pXim->GetXaxis()->SetLabelOffset(0.02);
	hAxis_pXim->GetXaxis()->SetTitleSize(0.075);
	hAxis_pXim->GetYaxis()->SetLabelSize(0.065);
	hAxis_pXim->GetYaxis()->CenterTitle();
	hAxis_pXim->GetYaxis()->SetTitleOffset(Yoffset);
	hAxis_pXim->GetYaxis()->SetTitleSize(0.075);
	hAxis_pXim->GetXaxis()->SetRangeUser(0, 500);
	hAxis_pXim->GetYaxis()->SetRangeUser(0.7, 4.5);    //pPb

	TCanvas* cfast = new TCanvas(
			TString::Format("cfast_r%.2f", GaussSourceSize),
			TString::Format("cfast_r%.2f", GaussSourceSize), 1);
	cfast->cd(0);
	cfast->SetCanvasSize(1920, 1280);
	cfast->SetMargin(0.15, 0.05, 0.2, 0.05);    //lrbt

	OliHisto_pXim->SetTitle("p#Xi^{#minus}");
	OliHisto_pXim->SetLineWidth(2);
	OliHisto_pXim->SetLineColor(kBlack);
	FitResult_pXim.SetLineWidth(2);
	FitResult_pXim.SetLineColor(kRed);
	FitResult_pXim.SetMarkerStyle(24);
	FitResult_pXim.SetMarkerColor(kRed);
	FitResult_pXim.SetMarkerSize(1);

	FitResult_pXim_COULOMB.SetLineWidth(2);
	FitResult_pXim_COULOMB.SetLineColor(kGreen);
	FitResult_pXim_COULOMB.SetMarkerStyle(24);
	FitResult_pXim_COULOMB.SetMarkerColor(kGreen);
	FitResult_pXim_COULOMB.SetMarkerSize(1);
	//                FitResult_pXim_COULOMB
	hAxis_pXim->Draw("axis");
	//                OliHisto_woSys_pXim->Draw("same");
	OliHisto_pXim->Draw("same");
	FitResult_pXim.Draw("CP,same");
	FitResult_pXim_COULOMB.Draw("CP,same");
	info4->Draw("same");
	cfast->SaveAs(
			TString::Format("%scfast_r%.2f.png", tmpOutputDir.Data(),
					GaussSourceSize));
	file->cd();
	if (saveCoulomb)
		FitResult_pXim_COULOMB.Write();
	FitResult_pXim.Write();
	return;
}

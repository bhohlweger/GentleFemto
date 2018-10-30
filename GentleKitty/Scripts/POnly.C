#include "ForBernie.h"
#include "DLM_Source.h"
#include "DLM_Potentials.h"
#include "DLM_CkModels.h"
#include "CATS.h"
#include "DLM_CkDecomposition.h"
#include "DLM_Fitters.h"
#include "TidyCats.h"
#include "CATSInput.h"
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
void FitPPVariations(const unsigned& NumIter, const unsigned& NumJobs,
		const unsigned& JobID, int system, TString InputDir, TString OutputDir);

int main(int argc, char *argv[]) {
	FitPPVariations(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]),
			argv[5], argv[6]);
	return 0;
}

void FitPPVariations(const unsigned& NumIter, const unsigned& NumJobs,
		const unsigned& JobID, int system, TString InputDir,
		TString OutputDir) {
	TRandom3 rangen(1 + JobID);
	TString HistName = "hCk_ReweightedMeV_";
	TString HistppName = HistName + "0";
	std::cout << HistppName.Data() << std::endl;
	TString InputFilePrefix = "CFOutput_";
	TString HistRestName = HistName;
	int Rebin = 5;
	if (Rebin == 4) {
		HistRestName += "0";
	} else if (Rebin == 5) {
		HistRestName += "1";
	} else {
		std::cout << "Non supported Rebinning (Rebin =" << Rebin
				<< "), please change your input file! \n";
		return;
	}
	//!SETTING THAT YOU MIGHT NEED
	//What source to use: 0 = Gauss; 1=Cauchy; 2=DoubleGauss
	//int vSource = rangen.Integer(2); if(vSource==1) vSource=3; //chose Gauss/EPOS at random

	const bool FAST_PLOT = true;
	const bool FULL_PLOT = false;
	const bool ExcludeXiSysError = true;
	const bool shiftedBinning = true;
	TString CalibBaseDir = "";
	std::cout << "SYSTEM: " << system << std::endl;
	if (system == 0) { // pPb MB
		CalibBaseDir += "~/cernbox/SystematicsAndCalib/pPbRun2_MB/";
	} else if (system == 1) { // pp MB
		CalibBaseDir += "~/cernbox/SystematicsAndCalib/ppRun2_MB/";
	} else if (system == 2) { // pp HM
		CalibBaseDir += "~/cernbox/SystematicsAndCalib/ppRun2_HM/";
	} else if (system == 11) { // pPb at the Grid
		CalibBaseDir +=
				"/home/gu74req/Analysis/CATS_Input/SystematicsAndCalib/pPbRun2_MB/";
		system = 0;
	}

	//just some temp folder for temp files. Really not important, I will try to get rid of this soon
	//  const TString OutputDir = "/Users/bernhardhohlweger/CATSOutput/";

	//perform the systematics by using the different cut variations as further permutations, do NOT add the systematic errors in quadrature
	//TString SystematicsType = "CutVarIterate";
	//perform the systematics by adding the systematics errors of the bin quadratically to the data points
	//  TString SystematicsType = "CutVarAdd";
	//  const bool ExcludeXi = false; //(fit only pp/pL/LL just as in run1)

	//!Binning
	//This is for the CATS objects, make sure it covers the full femto range
	const int binwidth = 4 * Rebin;
	const unsigned NumMomBins_pp = 50;
	const double kMin_pp = 0;
	const double kMax_pp = kMin_pp + 4 * NumMomBins_pp;  //(4 is the bin width)
	unsigned int momBins = 0;
	if (Rebin == 1) {
		momBins = 75;
	} else if (Rebin == 2) {
		momBins = 37;
	} else if (Rebin == 3) {
		momBins = 22;
	} else if (Rebin == 4) {
		momBins = 16;
	} else if (Rebin == 5) {
		momBins = 13;
	}
	const unsigned NumMomBins_pL = momBins;
	const double kMin_pL = 0;
	const double kMax_pL = kMin_pL + binwidth * NumMomBins_pL;

	const unsigned NumMomBins_LL = momBins;
	const double kMin_LL = 0;
	const double kMax_LL = kMin_LL + binwidth * NumMomBins_LL;

	const unsigned NumMomBins_pXim = momBins;
	double kMinXiP = 0;
	if (shiftedBinning) {
		if (system == 0) {
			kMinXiP = 8.;
		} else if (system == 1) {
			kMinXiP = 12.;
		} else if (system == 2) {
			kMinXiP = 0.;
		}
	}
	const double kMin_pXim = kMinXiP;
	const double kMax_pXim = kMin_pXim + binwidth * NumMomBins_pXim;

	std::cout << "kMinXiP: " << kMin_pXim << std::endl;
	std::cout << "kMax_pXim: " << kMax_pXim << std::endl;
	std::cout << "binwidth: " << binwidth << std::endl;
	std::cout << "NumMomBins_pXim: " << NumMomBins_pXim << std::endl;

	//!The Femtoregions for pp/pL/LL/pXim
	//if you modify you may need to change the CATS ranges somewhere below
	double FemtoRegion_pp[3][2];
	FemtoRegion_pp[0][0] = kMin_pp;
	FemtoRegion_pp[0][1] = 120;
	FemtoRegion_pp[1][0] = kMin_pp;
	FemtoRegion_pp[1][1] = 160;
	FemtoRegion_pp[2][0] = kMin_pp;
	FemtoRegion_pp[2][1] = 200;

	double FemtoRegion_pL[3][2];
	FemtoRegion_pL[0][0] = kMin_pL;
	FemtoRegion_pL[0][1] = 180;
	FemtoRegion_pL[1][0] = kMin_pL;
	FemtoRegion_pL[1][1] = 220;
	FemtoRegion_pL[2][0] = kMin_pL;
	FemtoRegion_pL[2][1] = 260;

	double FemtoRegion_LL[3][2];
	FemtoRegion_LL[0][0] = kMin_LL;
	FemtoRegion_LL[0][1] = 180;
	FemtoRegion_LL[1][0] = kMin_LL;
	FemtoRegion_LL[1][1] = 220;
	FemtoRegion_LL[2][0] = kMin_LL;
	FemtoRegion_LL[2][1] = 260;

	double FemtoRegion_pXim[3][2];
	FemtoRegion_pXim[0][0] = kMin_pXim;
	FemtoRegion_pXim[0][1] = 180;
	FemtoRegion_pXim[1][0] = kMin_pXim;
	FemtoRegion_pXim[1][1] = 220;
	FemtoRegion_pXim[2][0] = kMin_pXim;
	FemtoRegion_pXim[2][1] = 260;

	//!The baseline region (the same for pp systems)
	double BlRegion_pp[3][2];
	BlRegion_pp[0][0] = 250;
	BlRegion_pp[0][1] = 450;

	BlRegion_pp[1][0] = 250;
	BlRegion_pp[1][1] = 400;

	BlRegion_pp[2][0] = 350;
	BlRegion_pp[2][1] = 550;

	//!The baseline region (the same for all other systems)
	double BlRegion[3][2];
	BlRegion[0][0] = 320;
	BlRegion[0][1] = 480;
	BlRegion[1][0] = 300;
	BlRegion[1][1] = 500;
	BlRegion[2][0] = 300;
	BlRegion[2][1] = 540;

	double PurityProton;
	double PurityLambda;
	double PurityXi;
	double pp_f0;
	double pp_f1;
	double pL_f0;
	double pL_f1;
	double pL_f2;

	//(single particle quantities)
	//pPb
	if (system == 0) {
		PurityProton = 0.984266;  //pPb 5 TeV
		PurityLambda = 0.937761;
		PurityXi = 0.88;  //new cuts

		pp_f0 = 0.862814;
		pp_f1 = 0.09603;

		pL_f0 = 0.521433;  //fraction of primary Lambdas
		pL_f1 = 0.173811;  //fraction of Sigma0
		pL_f2 = 0.152378;  //fractions of Xi0/m
	} else { // pp MB + HM
		PurityProton = 0.991213;
		PurityLambda = 0.965964;
		PurityXi = 0.956;

		pp_f0 = 0.874808;
		pp_f1 = 0.0876342; //fraction of

		pL_f0 = 0.619493; //fraction of primary Lambdas
		pL_f1 = 0.206498; //fraction of Sigma0
		pL_f2 = 0.0870044; //fractions of Xi0/m
	}
	//ADVANCED***
	double ProtonPrim = pp_f0;
	double arrayPercLamProton[3] = { pp_f1 / (1. - pp_f0) * 0.8, pp_f1
			/ (1. - pp_f0), pp_f1 / (1. - pp_f0) * 1.2 }; //+/- 20%

	const double SigLambdaPrimDir = pL_f0 + pL_f1;
	double arrayPercSigLambda[3] = { 0.8 * pL_f1 / pL_f0, pL_f1 / pL_f0, 1.2
			* pL_f1 / pL_f0 }; //1/3 +/- 20%
	double arrayPercXiLambda[3] = { pL_f2 / (1. - pL_f0 - pL_f1) * 0.8, pL_f2
			/ (1. - pL_f0 - pL_f1), pL_f2 / (1. - pL_f0 - pL_f1) * 1.2 }; //+/- 20%

	//ratio Xi-(1530) to Xi-
	const double Xim1530_to_Xim = 0.32 * (1. / 3.);
	//ratio Xi0(1530) to Xi0 (n=neutral)
	const double Xin1530_to_Xim = 0.32 * (2. / 3.);
	const double Omegam_to_Xim = 0.1;
	const double OmegamXim_BR = 0.086;

	//following my lambda pars with the 3 possible modifications
	//for the proton:
	//0 = primary
	//1 = from Lambda
	//2 = other feeddown (flat)
	//3 = missidentified
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

	//for the Lambda:
	//0 = primary
	//1 = from Sigma0
	//2 = from Xim
	//3 = from Xi0
	//4 = missidentified
	const unsigned NumChannels_L = 5;
	double*** Purities_L = new double**[3];
	double*** Fraction_L = new double**[3];
	for (unsigned uVarSL = 0; uVarSL < 3; uVarSL++) {
		Purities_L[uVarSL] = new double*[3];
		Fraction_L[uVarSL] = new double*[3];
		for (unsigned uVarXi = 0; uVarXi < 3; uVarXi++) {
			Purities_L[uVarSL][uVarXi] = new double[NumChannels_L];
			Fraction_L[uVarSL][uVarXi] = new double[NumChannels_L];

			Purities_L[uVarSL][uVarXi][0] = PurityLambda;
			Purities_L[uVarSL][uVarXi][1] = PurityLambda;
			Purities_L[uVarSL][uVarXi][2] = PurityLambda;
			Purities_L[uVarSL][uVarXi][3] = PurityLambda;
			Purities_L[uVarSL][uVarXi][4] = 1. - PurityLambda;

			//the array is r = S/L, and S+L=1 are the fractions of Sigmas and Lambdas
			double FracOfLambda = 1. / (1. + arrayPercSigLambda[uVarSL]);
			Fraction_L[uVarSL][uVarXi][0] = SigLambdaPrimDir * FracOfLambda;
			Fraction_L[uVarSL][uVarXi][1] = SigLambdaPrimDir
					* (1. - FracOfLambda);
			Fraction_L[uVarSL][uVarXi][2] = (1. - SigLambdaPrimDir)
					* (arrayPercXiLambda[uVarXi]);
			Fraction_L[uVarSL][uVarXi][3] = (1. - SigLambdaPrimDir)
					* (1. - arrayPercXiLambda[uVarXi]);
			Fraction_L[uVarSL][uVarXi][4] = 1.;
		}
	}

	//for the Xi:
	//0 = primary
	//1 = from Xi-(1530)
	//2 = from Xi0(1530)
	//3 = from Omega
	//4 = missidentified
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

	//starting value, do not worry about it too much
	const double GaussSourceSize = 1.2;

	CATSInput *CalibFiles = new CATSInput();
	CalibFiles->SetBaseDir(CalibBaseDir.Data());
	CalibFiles->SetMomResFileName("run2_decay_matrices_old.root");
	CalibFiles->ReadResFile();
	CalibFiles->SetSigmaFileName("Sample3_MeV_compact.root");
	CalibFiles->ReadSigmaFile();

	//  int vCutID;//which data file (cut combination) should you take. 0 = default
	//int vSource;//which source we use, see above
	int vFemReg_pp;  //which femto region we use for pp (1 = default)
	int vFemReg_pL;  //which femto region we use for pL (1 = default)
	int vFemReg_LL;  //which femto region we use for LL (1 = default)
	int vFemReg_pXim;  //which femto region we use for pXim (1 = default)
	int vBlReg;  //which baseline region to use (1 = default)
	int vFrac_pp_pL; //fraction of protons coming from Lambda variation (1 = default)
	int vFrac_pL_pSigma0; //fraction of Lambdas coming from Sigma0 variation (1 = default)
	int vFrac_pL_pXim; //fraction of Lambdas coming from Xim (1 = default)

	//each JOB produces a separate output file
	TFile* OutFile = new TFile(
			TString::Format("%sOutFile_%s_Iter%u_JOBS%u_ID%u.root",
					OutputDir.Data(), "CutVarAdd", NumIter, NumJobs, JobID),
			"recreate");
	//you save a lot of stuff in an NTuple
	TNtuple* ntResult =
			new TNtuple("ntResult", "ntResult",
					"IterID:vFemReg_pp:vFemReg_pL:vFemReg_pXim:vBlReg:vFrac_pp_pL:"
							"vFrac_pL_pSigma0:vFrac_pL_pXim:Radius_pp:RadiusErr_pp:"
							"pa_pp:paErr_pp:pb_pp:pbErr_pp:pCl_pp:pClErr_pp:"
							"Radius_pL:RadiusErr_pL:Radius_pXim:RadiusErr_pXim:Chi2Ndf:pval:SEPARATE_BL:FIX_CL");
	Float_t ntBuffer[31];

	unsigned WhichPart = JobID + 1;
	unsigned NumJobsPart = NumIter;

	unsigned iSplitInto = NumJobs;
	unsigned FirstIter = 0;
	unsigned LastIter = 0;
	Printf("JobID = %i", JobID);
	while (iSplitInto > 0 && WhichPart > 0) {
		FirstIter = LastIter + bool(LastIter);
		LastIter += NumJobsPart / iSplitInto - 1 + bool(LastIter);
		NumJobsPart = NumIter - LastIter - 1;
		iSplitInto--;
		WhichPart--;
	}
	if (NumIter == NumJobs) {
		FirstIter = JobID;
		LastIter = JobID;
	}

	Printf(" FirstIter = %i", FirstIter);
	Printf(" LastIter = %i", LastIter);
	Printf("  Nruns = %i", LastIter - FirstIter + 1);
	//***
	//default -> True, false
	//true => do the BL separately (as an RUN1), false => fit femto and BL region together (RUN2)
	bool SEPARATE_BL = true;
	//if true renormalization is NOT allowed
	bool FIX_CL = false;
	TidyCats* tidy = new TidyCats();
	//the 0 iter is always the default!
	for (unsigned uIter = FirstIter; uIter <= LastIter; uIter++) {
		//    vCutID=(SystematicsType=="CutVarIterate")?rangen.Integer(31):0;
		vFemReg_pp = rangen.Integer(3);
		vFemReg_pL = rangen.Integer(3);
		vFemReg_LL = rangen.Integer(3);
		vFemReg_pXim = rangen.Integer(3);
		vBlReg = rangen.Integer(3);
		vFrac_pp_pL = rangen.Integer(3);
		vFrac_pL_pSigma0 = rangen.Integer(3);
		vFrac_pL_pXim = rangen.Integer(3);

		if (rangen.Integer(2) == 1) {
			SEPARATE_BL = true;
			FIX_CL = false;
		} else {
			SEPARATE_BL = false;
			FIX_CL = true; // variation
		}

		//The defaults
		if (uIter == 0) {
			//      vCutID = 0;
			vFemReg_pp = 1;
			vFemReg_pL = 1;
			vFemReg_LL = 1;
			vFemReg_pXim = 1;
			vBlReg = 1;
			vFrac_pp_pL = 1;
			vFrac_pL_pSigma0 = 1;
			vFrac_pL_pXim = 1;
			SEPARATE_BL = true;
			FIX_CL = false;
		}
//    vMod_pL = 1; // 0 = Usmani (Close to NLO) 1 = NLO Lednicky 2 = LO Lednicky (default)
//    vMod_pXim = 0;
		printf("uIter=%u\n", uIter);
		//    printf("vCutID=%u\n",vCutID);
		printf("vFemReg_pp=%u\n", vFemReg_pp);
		printf("vFemReg_pL=%u\n", vFemReg_pL);
		printf("vFemReg_LL=%u\n", vFemReg_LL);
		printf("vFemReg_pXim=%u\n", vFemReg_pXim);
		printf("vBlReg=%u\n", vBlReg);
		printf("vFrac_pp_pL=%u\n", vFrac_pp_pL);
		printf("vFrac_pL_pSigma0=%u\n", vFrac_pL_pSigma0);
		printf("vFrac_pL_pXim=%u\n", vFrac_pL_pXim);
		printf("SEPERATE_BL=%u\n", SEPARATE_BL);
		printf("FIX_CL=%u\n", FIX_CL);

		//    InputFilePrefix
		//    InputFileSuffix

		//ADVANCED***
		//#,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
		double PotPars1S0[10] = { 0, 0, NN_AV18, v18_Coupled3P2, 1, 1, 1, 0, 0,
				0 };
		double PotPars3P0[10] = { 0, 0, NN_AV18, v18_Coupled3P2, 1, 1, 1, 1, 1,
				0 };
		double PotPars3P1[10] = { 0, 0, NN_AV18, v18_Coupled3P2, 1, 1, 1, 1, 1,
				1 };
		double PotPars3P2[10] = { 0, 0, NN_AV18, v18_Coupled3P2, 1, 1, 1, 1, 1,
				2 };

		const double Weight1S0 = 3. / 12.;
		const double Weight3P0 = 1. / 12.;
		const double Weight3P1 = 3. / 12.;
		const double Weight3P2 = 5. / 12.;

		const double Mass_p = 938.272;
		const double Mass_L = 1115.683;

		double Pars_pp[6] = { 0, 0, 0, GaussSourceSize * 1.2, GaussSourceSize
				/ 1.2, 0.5 };
		CATS AB_pp;
		tidy->GetCatsProtonProton(&AB_pp, GaussSourceSize, Pars_pp,
				NumMomBins_pp, kMin_pp, kMax_pp);
		AB_pp.KillTheCat();

		double Pars_pL[6] = { 0, 0, 0, GaussSourceSize * 1.2, GaussSourceSize
				/ 1.2, 0.5 };
		CATS AB_pL;
		tidy->GetCatsProtonLambda(&AB_pL, GaussSourceSize, Pars_pL,
				NumMomBins_pL, kMin_pL, kMax_pL);
		AB_pL.KillTheCat();
		double Pars_pXi[6] = { 0, 0, 0, GaussSourceSize * 1.2, GaussSourceSize
				/ 1.2, 0.5 };
		CATS AB_pXim;
		tidy->GetCatsProtonXiMinus(&AB_pXim, GaussSourceSize, Pars_pXi,
				NumMomBins_pXim, kMin_pXim, kMax_pXim, true, 13);
		AB_pXim.KillTheCat();

		double Pars_pXim1530[6] = { 0, 0, 0, GaussSourceSize * 1.2,
				GaussSourceSize / 1.2, 0.5 };
		CATS AB_pXim1530;
		tidy->GetCatsProtonXiMinus1530(&AB_pXim1530, GaussSourceSize,
				Pars_pXim1530, NumMomBins_pXim, kMin_pXim, kMax_pXim);
		AB_pXim1530.KillTheCat();

		//***
		std::cout << "Reading Data \n";
		//! DATA FILE
		//    TString OliFileName_pp =
		TString OliFileName_pp = TString::Format("%s%spp.root", InputDir.Data(),
				InputFilePrefix.Data());
		TFile* OliFile_pp =
				OliFileName_pp != "" ? new TFile(OliFileName_pp, "read") : NULL;
		TH1F* OliHisto_pp =
				OliFile_pp ? (TH1F*) OliFile_pp->Get(HistppName.Data()) :
				NULL;
		if (OliHisto_pp)
			std::cout << OliHisto_pp->GetName() << std::endl;
		else {
			std::cout << OliFileName_pp.Data() << "--" << HistppName.Data()
					<< " Missing" << std::endl;
		}
		//!CHANGE PATH HERE
		TString SystErrFileName_pp = TString::Format("%s/C2totalsysPP.root",
				CalibBaseDir.Data());
		TFile* SystErrFile_pp =
				SystErrFileName_pp != "" ?
						new TFile(SystErrFileName_pp, "read") : NULL;
		TH1F* outputParamPP = (TH1F*) SystErrFile_pp->Get("SysParamPP");
		std::cout << "PP" << std::endl;
		std::cout << outputParamPP->GetBinContent(1) << std::endl;
		std::cout << outputParamPP->GetBinContent(2) << std::endl;
		std::cout << outputParamPP->GetBinContent(3) << std::endl;
		TF1 *RelSystPP = new TF1("sysPP", "pol2", 0, 3);
		RelSystPP->SetParameter(0, outputParamPP->GetBinContent(1));
		RelSystPP->SetParameter(1, outputParamPP->GetBinContent(2));
		RelSystPP->SetParameter(2, outputParamPP->GetBinContent(3));

		int NumSEB_pp = RelSystPP == NULL ? 0 : OliHisto_pp->FindBin(500);

		for (int iBin = 0; iBin < NumSEB_pp; iBin++) {
			const float x = OliHisto_pp->GetBinCenter(iBin + 1);
			const float y = OliHisto_pp->GetBinContent(iBin + 1);
			OliHisto_pp->SetBinError(iBin + 1,
					sqrt(
							pow(OliHisto_pp->GetBinError(iBin + 1), 2.)
									+ pow(y * RelSystPP->Eval(x / 1000.), 2.)));
		}

		//!CHANGE PATH HERE
		//    TString OliFileName_pL =
		TString OliFileName_pL = TString::Format("%s%spL.root", InputDir.Data(),
				InputFilePrefix.Data());
		TFile* OliFile_pL;
		if (OliFileName_pp == OliFileName_pL)
			OliFile_pL = OliFile_pp;
		else
			OliFile_pL =
					OliFileName_pL != "" ?
							new TFile(OliFileName_pL, "read") : NULL;
		TH1F* OliHisto_pL =
				OliFile_pL ? (TH1F*) OliFile_pL->Get(HistRestName.Data()) :
				NULL;

		//!CHANGE PATH HERE
		TString SystErrFileName_pL = TString::Format("%s/C2totalsysPL.root",
				CalibBaseDir.Data());
		TFile* SystErrFile_pL =
				SystErrFileName_pL != "" ?
						new TFile(SystErrFileName_pL, "read") : NULL;
		TH1F* outputParamPL = (TH1F*) SystErrFile_pL->Get("SysParamPL");
		std::cout << "PL" << std::endl;
		std::cout << outputParamPL->GetBinContent(1) << std::endl;
		std::cout << outputParamPL->GetBinContent(2) << std::endl;
		std::cout << outputParamPL->GetBinContent(3) << std::endl;
		TF1 *RelSystPL = new TF1("sysPL", "pol2", 0, 3);
		RelSystPL->SetParameter(0, outputParamPL->GetBinContent(1));
		RelSystPL->SetParameter(1, outputParamPL->GetBinContent(2));
		RelSystPL->SetParameter(2, outputParamPL->GetBinContent(3));

		int NumSEB_pL = RelSystPL == NULL ? 0 : OliHisto_pL->FindBin(500);

		for (int iBin = 0; iBin < NumSEB_pL; iBin++) {
			const float x = OliHisto_pL->GetBinCenter(iBin + 1);
			const float y = OliHisto_pL->GetBinContent(iBin + 1);
			OliHisto_pL->SetBinError(iBin + 1,
					sqrt(
							pow(OliHisto_pL->GetBinError(iBin + 1), 2.)
									+ pow(y * RelSystPL->Eval(x / 1000.), 2.)));
		}

		//!CHANGE PATH HERE
		//    TString OliFileName_LL =
		TString OliFileName_LL = TString::Format("%s%sLL.root", InputDir.Data(),
				InputFilePrefix.Data());
		TFile* OliFile_LL;
		if (OliFileName_pp == OliFileName_LL)
			OliFile_LL = OliFile_pp;
		else
			OliFile_LL =
					OliFileName_LL != "" ?
							new TFile(OliFileName_LL, "read") : NULL;
		TH1F* OliHisto_LL =
				OliFile_LL ? (TH1F*) OliFile_LL->Get(HistRestName.Data()) :
				NULL;

		//!CHANGE PATH HERE
		TString SystErrFileName_LL = TString::Format("%s/C2totalsysLL.root",
				CalibBaseDir.Data());
		TFile* SystErrFile_LL =
				SystErrFileName_LL != "" ?
						new TFile(SystErrFileName_LL, "read") : NULL;
		TH1F* outputParamLL = (TH1F*) SystErrFile_LL->Get("SysParamLL");
		std::cout << "LL" << std::endl;
		std::cout << outputParamLL->GetBinContent(1) << std::endl;
		std::cout << outputParamLL->GetBinContent(2) << std::endl;
		std::cout << outputParamLL->GetBinContent(3) << std::endl;

		TF1 *RelSystLL = new TF1("sysLL", "pol2", 0, 3);
		RelSystLL->SetParameter(0, outputParamLL->GetBinContent(1));
		RelSystLL->SetParameter(1, outputParamLL->GetBinContent(2));
		RelSystLL->SetParameter(2, outputParamLL->GetBinContent(3));

		int NumSEB_LL = RelSystLL == NULL ? 0 : OliHisto_LL->FindBin(500);

		for (int iBin = 0; iBin < NumSEB_LL; iBin++) {
			const float x = OliHisto_LL->GetBinCenter(iBin + 1);
			const float y = OliHisto_LL->GetBinContent(iBin + 1);
			OliHisto_LL->SetBinError(iBin + 1,
					sqrt(
							pow(OliHisto_LL->GetBinError(iBin + 1), 2.)
									+ pow(y * RelSystLL->Eval(x / 1000.), 2.)));
		}

		//!CHANGE PATH HERE
		TString OliFileName_pXim = TString::Format("%s%spXi.root",
				InputDir.Data(), InputFilePrefix.Data());
		TFile* OliFile_pXim;
		if (OliFileName_pp == OliFileName_pXim)
			OliFile_pXim = OliFile_pp;
		else
			OliFile_pXim =
					OliFileName_pXim != "" ?
							new TFile(OliFileName_pXim, "read") : NULL;

		TH1F* OliHisto_pXim =
				OliFile_pXim ?
						(TH1F*) OliFile_pXim->Get(
								Form("%s", HistRestName.Data())) :
						NULL;
		std::cout << OliHisto_pXim->GetXaxis()->GetNbins() << '\t'
				<< OliHisto_pXim->GetXaxis()->GetXmin() << '\t'
				<< OliHisto_pXim->GetXaxis()->GetXmax() << '\n';
		//!CHANGE PATH HERE
		TString SystErrFileName_pXim = Form("%s/C2totalsysPXi.root",
				CalibBaseDir.Data());
		TFile* SystErrFile_pXim =
				SystErrFileName_pXim != "" ?
						new TFile(SystErrFileName_pXim, "read") : NULL;
		TH1F* outputParamPXi = (TH1F*) SystErrFile_pXim->Get("SysParamPXi");
		std::cout << "PXi" << std::endl;
		std::cout << outputParamPXi->GetBinContent(1) << std::endl;
		std::cout << outputParamPXi->GetBinContent(2) << std::endl;
		std::cout << outputParamPXi->GetBinContent(3) << std::endl;
		TF1 *RelSystpXi = new TF1("sysPXi", "pol2", 0, 3);
		RelSystpXi->SetParameter(0, outputParamPXi->GetBinContent(1));
		RelSystpXi->SetParameter(1, outputParamPXi->GetBinContent(2));
		RelSystpXi->SetParameter(2, outputParamPXi->GetBinContent(3));

		int NumSEB_pXim = RelSystpXi == NULL ? 0 : OliHisto_pXim->FindBin(500);
		TH1F *OliHisto_pXimFornSigma = nullptr;
		if (!ExcludeXiSysError) {
			for (int iBin = 0; iBin < NumSEB_pXim; iBin++) {
				const float x = OliHisto_pXim->GetBinCenter(iBin + 1);
				const float y = OliHisto_pXim->GetBinContent(iBin + 1);
				OliHisto_pXim->SetBinError(iBin + 1,
						sqrt(
								pow(OliHisto_pXim->GetBinError(iBin + 1), 2.)
										+ pow(y * RelSystpXi->Eval(x / 1000.),
												2.)));
			}
			OliHisto_pXimFornSigma = (TH1F*) OliHisto_pXim->Clone(
					"pXiForNSigma");
		} else {
			OliHisto_pXimFornSigma = (TH1F*) OliHisto_pXim->Clone(
					"pXiForNSigma");
			for (int iBin = 0; iBin < NumSEB_pXim; iBin++) {
				const float x = OliHisto_pXim->GetBinCenter(iBin + 1);
				const float y = OliHisto_pXim->GetBinContent(iBin + 1);
				OliHisto_pXim->SetBinError(iBin + 1,
						sqrt(
								pow(OliHisto_pXim->GetBinError(iBin + 1), 2.)
										+ pow(y * RelSystpXi->Eval(x / 1000.),
												2.)));
			}
		}
		const unsigned NumSourcePars = 1;

		//this way you define a correlation function using a CATS object.
		//needed inputs: num source/pot pars, CATS obj
		DLM_Ck* Ck_pp = new DLM_Ck(NumSourcePars, 0, AB_pp);
		DLM_Ck* Ck_pL = new DLM_Ck(NumSourcePars, 0, AB_pL);
		//this way you define a correlation function using Lednicky.
		//needed inputs: num source/pot pars, mom. binning, pointer to a function which computes C(k)
		DLM_Ck* Ck_LL = new DLM_Ck(1, 2, NumMomBins_LL, kMin_LL, kMax_LL,
				Lednicky_Identical_Singlet);
		DLM_Ck* Ck_pSigma0 = new DLM_Ck(1, 0, NumMomBins_pL, kMin_pL, kMax_pL,
				Lednicky_gauss_Sigma0);
		Ck_pSigma0->SetSourcePar(0, GaussSourceSize);
		//DLM_Ck* Ck_pXiMinus = new DLM_Ck(1,0,NumMomBins_pL,kMin_pL,kMax_pL,pXi_pheno);
		//Ck_pXiMinus->SetSourcePar(0,3.88);
		DLM_Ck* Ck_pXim = new DLM_Ck(NumSourcePars, 0, AB_pXim);
		DLM_Ck* Ck_pXim1530 = new DLM_Ck(NumSourcePars, 0, AB_pXim1530);

		Ck_pL->SetSourcePar(0, GaussSourceSize);

		Ck_LL->SetSourcePar(0, GaussSourceSize);
		Ck_LL->SetPotPar(0, 1.2);
		Ck_LL->SetPotPar(1, 4.5);

		Ck_pp->Update();
		Ck_pL->Update();
		Ck_LL->Update();
		Ck_pSigma0->Update();
		Ck_pXim->Update();
		Ck_pXim1530->Update();
		if (CalibFiles->GetSigmaFile(0)) {
			std::cout << "No Sigma file 0 \n";
			return;
		}
		if (CalibFiles->GetSigmaFile(1)) {
			std::cout << "No Sigma file 1 \n";
			return;
		}
		if (CalibFiles->GetSigmaFile(2)) {
			std::cout << "No Sigma file 2 \n";
			return;
		}
		if (CalibFiles->GetSigmaFile(3)) {
			std::cout << "No Sigma file 3 \n";
			return;
		}
		DLM_CkDecomposition CkDec_pp("pp", 3, *Ck_pp,
				CalibFiles->GetSigmaFile(0));
		DLM_CkDecomposition CkDec_pL("pLambda", 4, *Ck_pL,
				CalibFiles->GetSigmaFile(1));
		DLM_CkDecomposition CkDec_LL("LambdaLambda", 2, *Ck_LL,
				CalibFiles->GetSigmaFile(2));
		DLM_CkDecomposition CkDec_pSigma0("pSigma0", 0, *Ck_pSigma0,
		NULL);
		DLM_CkDecomposition CkDec_pXim("pXim", 3, *Ck_pXim,
				CalibFiles->GetSigmaFile(3));
		DLM_CkDecomposition CkDec_pXim1530("pXim1530", 0, *Ck_pXim1530,
		NULL);

		double lam_pp = Purities_p[vFrac_pp_pL][0] * Fraction_p[vFrac_pp_pL][0]
				* Purities_p[vFrac_pp_pL][0] * Fraction_p[vFrac_pp_pL][0];
		double lam_pp_pL = Purities_p[vFrac_pp_pL][0]
				* Fraction_p[vFrac_pp_pL][0] * Purities_p[vFrac_pp_pL][1]
				* Fraction_p[vFrac_pp_pL][1] * 2;
		double lam_pp_fake = Purities_p[vFrac_pp_pL][3]
				* Purities_p[vFrac_pp_pL][0]
				+ Purities_p[vFrac_pp_pL][0] * Purities_p[vFrac_pp_pL][3]
				+ Purities_p[vFrac_pp_pL][3] * Purities_p[vFrac_pp_pL][3];

		printf("lam_pp = %.3f\n", lam_pp);
		printf("lam_pp_pL = %.3f\n", lam_pp_pL);
		printf("lam_pp_fake = %.3f\n", lam_pp_fake);
		printf("\n");
		if (!CalibFiles->GetResFile(0)) {
			std::cout << "No Calib 0 \n";
			return;
		}
		CkDec_pp.AddContribution(0, lam_pp_pL, DLM_CkDecomposition::cFeedDown,
				&CkDec_pL, CalibFiles->GetResFile(0));
		CkDec_pp.AddContribution(1, 1. - lam_pp - lam_pp_pL - lam_pp_fake,
				DLM_CkDecomposition::cFeedDown);
		CkDec_pp.AddContribution(2, lam_pp_fake, DLM_CkDecomposition::cFake); //0.02

		double lam_pL = Purities_p[vFrac_pp_pL][0] * Fraction_p[vFrac_pp_pL][0]
				* Purities_L[vFrac_pL_pSigma0][vFrac_pL_pXim][0]
				* Fraction_L[vFrac_pL_pSigma0][vFrac_pL_pXim][0];
		double lam_pL_pS0 = Purities_p[vFrac_pp_pL][0]
				* Fraction_p[vFrac_pp_pL][0]
				* Purities_L[vFrac_pL_pSigma0][vFrac_pL_pXim][1]
				* Fraction_L[vFrac_pL_pSigma0][vFrac_pL_pXim][1];
		double lam_pL_pXm = Purities_p[vFrac_pp_pL][0]
				* Fraction_p[vFrac_pp_pL][0]
				* Purities_L[vFrac_pL_pSigma0][vFrac_pL_pXim][2]
				* Fraction_L[vFrac_pL_pSigma0][vFrac_pL_pXim][2];
		double lam_pL_fake = Purities_p[vFrac_pp_pL][3]
				* Purities_L[vFrac_pL_pSigma0][vFrac_pL_pXim][0]
				+ Purities_p[vFrac_pp_pL][0]
						* Purities_L[vFrac_pL_pSigma0][vFrac_pL_pXim][4]
				+ Purities_p[vFrac_pp_pL][3]
						* Purities_L[vFrac_pL_pSigma0][vFrac_pL_pXim][4];

		printf("lam_pL=%.3f\n", lam_pL);
		printf("lam_pL_pS0=%.3f\n", lam_pL_pS0);
		printf("lam_pL_pXm=%.3f\n", lam_pL_pXm);
		printf("lam_pL_fake=%.3f\n", lam_pL_fake);
		printf("\n");
		if (!CalibFiles->GetResFile(1)) {
			std::cout << "No Calib 1 \n";
			return;
		}
		if (!CalibFiles->GetResFile(2)) {
			std::cout << "No Calib 2 \n";
			return;
		}

		CkDec_pL.AddContribution(0, lam_pL_pS0, DLM_CkDecomposition::cFeedDown,
				&CkDec_pSigma0, CalibFiles->GetResFile(1));
		CkDec_pL.AddContribution(1, lam_pL_pXm, DLM_CkDecomposition::cFeedDown,
				&CkDec_pXim, CalibFiles->GetResFile(2));
		CkDec_pL.AddContribution(2,
				1. - lam_pL - lam_pL_pS0 - lam_pL_pXm - lam_pL_fake,
				DLM_CkDecomposition::cFeedDown);
		CkDec_pL.AddContribution(3, lam_pL_fake, DLM_CkDecomposition::cFake); //0.03

		double lam_LL = Purities_L[vFrac_pL_pSigma0][vFrac_pL_pXim][0]
				* Fraction_L[vFrac_pL_pSigma0][vFrac_pL_pXim][0]
				* Purities_L[vFrac_pL_pSigma0][vFrac_pL_pXim][0]
				* Fraction_L[vFrac_pL_pSigma0][vFrac_pL_pXim][0];
		double lam_LL_fake = Purities_L[vFrac_pL_pSigma0][vFrac_pL_pXim][4]
				* Purities_L[vFrac_pL_pSigma0][vFrac_pL_pXim][0]
				+ Purities_L[vFrac_pL_pSigma0][vFrac_pL_pXim][0]
						* Purities_L[vFrac_pL_pSigma0][vFrac_pL_pXim][4]
				+ Purities_L[vFrac_pL_pSigma0][vFrac_pL_pXim][4]
						* Purities_L[vFrac_pL_pSigma0][vFrac_pL_pXim][4];

		printf("lam_LL=%.3f\n", lam_LL);
		printf("lam_LL_fake=%.3f\n", lam_LL_fake);
		printf("\n");

		CkDec_LL.AddContribution(0, 1. - lam_LL - lam_LL_fake,
				DLM_CkDecomposition::cFeedDown);    //0.65
		CkDec_LL.AddContribution(1, lam_LL_fake, DLM_CkDecomposition::cFake); //0.03

		const double lam_pXim = Purities_p[vFrac_pp_pL][0]
				* Fraction_p[vFrac_pp_pL][0] * Purities_Xim[0][0]
				* Fraction_Xim[0][0];
		const double lam_pXim_pXim1530 = Purities_p[vFrac_pp_pL][0]
				* Fraction_p[vFrac_pp_pL][0] * Purities_Xim[0][1]
				* Fraction_Xim[0][1];
		const double lam_pXim_fake = Purities_p[vFrac_pp_pL][3]
				* Purities_Xim[0][0]
				+ Purities_p[vFrac_pp_pL][0] * Purities_Xim[0][4]
				+ Purities_p[vFrac_pp_pL][3] * Purities_Xim[0][4];

		printf("lam_pXim = %.3f\n", lam_pXim);
		printf("lam_pXim_pXim1530 = %.3f\n", lam_pXim_pXim1530);
		printf("lam_pXim_fake = %.3f\n", lam_pXim_fake);
		printf("\n");
		if (!CalibFiles->GetResFile(3)) {
			std::cout << "No Calib 3 \n";
			return;
		}
		CkDec_pXim.AddContribution(0, lam_pXim_pXim1530,
				DLM_CkDecomposition::cFeedDown, &CkDec_pXim1530,
				CalibFiles->GetResFile(3));    //from Xi-(1530)
		CkDec_pXim.AddContribution(1,
				1. - lam_pXim - lam_pXim_pXim1530 - lam_pXim_fake,
				DLM_CkDecomposition::cFeedDown); //other feed-down (flat)
		CkDec_pXim.AddContribution(2, lam_pXim_fake,
				DLM_CkDecomposition::cFake);

		DLM_Fitter1* fitter = new DLM_Fitter1(4);
		fitter->SetOutputDir(OutputDir.Data());

		std::cout << "pp" << std::endl;
		fitter->SetSystem(0, *OliHisto_pp, 1, CkDec_pp,
				FemtoRegion_pp[vFemReg_pp][0], FemtoRegion_pp[vFemReg_pp][1],
				BlRegion_pp[1][0], BlRegion_pp[1][1]);
		std::cout << "pL" << std::endl;
		fitter->SetSystem(1, *OliHisto_pL, 1, CkDec_pL,
				FemtoRegion_pL[vFemReg_pL][0], FemtoRegion_pL[vFemReg_pL][1],
				BlRegion[vBlReg][0], BlRegion[vBlReg][1]);
		std::cout << "LL" << std::endl;
		fitter->SetSystem(2, *OliHisto_LL, 1, CkDec_LL,
				FemtoRegion_LL[vFemReg_LL][0], FemtoRegion_LL[vFemReg_LL][1],
				BlRegion[vBlReg][0], BlRegion[vBlReg][1]);
		std::cout << "pXi" << std::endl;
		fitter->SetSystem(3, *OliHisto_pXim, 1, CkDec_pXim,
				FemtoRegion_pXim[vFemReg_pXim][0],
				FemtoRegion_pXim[vFemReg_pXim][1], BlRegion[vBlReg][0],
				BlRegion[vBlReg][1]);
		if (SEPARATE_BL) {
			fitter->SetSeparateBL(0, true);
			fitter->SetSeparateBL(1, true);
			fitter->SetSeparateBL(2, true);
			fitter->SetSeparateBL(3, true);
			fitter->FixParameter("pp", DLM_Fitter1::p_a, 1.0);
			fitter->FixParameter("pp", DLM_Fitter1::p_b, 0);
		} else {
			fitter->SetSeparateBL(0, false);
			fitter->SetSeparateBL(1, false);
			fitter->SetSeparateBL(2, false);
			fitter->SetSeparateBL(3, false);
			fitter->SetParameter("pp", DLM_Fitter1::p_a, 1.0, 0.7, 1.3);
			fitter->SetParameter("pp", DLM_Fitter1::p_b, 1e-4, 0, 2e-3);
			std::cout << "Fitting ranges for BL set \n";
		}
		fitter->AddSameSource("LambdaLambda", "pp", 1);
		fitter->AddSameSource("pSigma0", "pLambda", 1);
		fitter->AddSameSource("pXim1530", "pXim", 1);

		fitter->FixParameter("pp", DLM_Fitter1::p_c, 0);

		fitter->SetParameter("pp", DLM_Fitter1::p_sor0, 1.2, 0.8, 1.8);
//		fitter->FixParameter("pp", DLM_Fitter1::p_sor0, 1.2);
		if (FIX_CL) {
			fitter->FixParameter("pp", DLM_Fitter1::p_Cl, -1);
			std::cout << "CL Fixed \n";
		} else {
			fitter->SetParameter("pp", DLM_Fitter1::p_Cl, -0.9, -1.2, -0.8);
		}

		fitter->SetParameter("pLambda", DLM_Fitter1::p_a, 1.0, 0.7, 1.3);
		fitter->SetParameter("pLambda", DLM_Fitter1::p_b, 1e-4, 0, 2e-3);
		fitter->FixParameter("pLambda", DLM_Fitter1::p_c, 0);

		fitter->SetParameter("pLambda", DLM_Fitter1::p_sor0, GaussSourceSize,
				0.8, 1.8);

		if (FIX_CL)
			fitter->FixParameter("pLambda", DLM_Fitter1::p_Cl, -1);
		else
			fitter->SetParameter("pLambda", DLM_Fitter1::p_Cl, -0.9, -1.2,
					-0.8);

		fitter->SetParameter("LambdaLambda", DLM_Fitter1::p_a, 1.0, 0.7, 1.3);
		fitter->SetParameter("LambdaLambda", DLM_Fitter1::p_b, 1e-4, 0, 2e-3);
		fitter->FixParameter("LambdaLambda", DLM_Fitter1::p_c, 0);
		//if(FIX_CL) fitter->FixParameter("LambdaLambda",DLM_Fitter1::p_Cl,-1);
		//else fitter->SetParameter("LambdaLambda",DLM_Fitter1::p_Cl,-0.9,-1.2,-0.8);
		//!fix Cl for now, since the LambdaLambda is very shallow end prob. there is very little deviation from unity.
		fitter->FixParameter("LambdaLambda", DLM_Fitter1::p_Cl, -1);
		fitter->SetParameter("LambdaLambda", DLM_Fitter1::p_pot0, 4, -8, 24);
		fitter->SetParameter("LambdaLambda", DLM_Fitter1::p_pot1, 8, 0, 24);

		fitter->SetParameter("pXim", DLM_Fitter1::p_a, 1.0, 0.7, 1.3);
		fitter->SetParameter("pXim", DLM_Fitter1::p_b, 1e-4, 0, 2e-3);
		fitter->FixParameter("pXim", DLM_Fitter1::p_c, 0);

		fitter->SetParameter("pXim", DLM_Fitter1::p_sor0, GaussSourceSize, 0.8,
				1.8);
		if (FIX_CL)
			fitter->FixParameter("pXim", DLM_Fitter1::p_Cl, -1);
		else
			fitter->SetParameter("pXim", DLM_Fitter1::p_Cl, -0.9, -1.2, -0.8);

//		fitter->FixParameter("pp", DLM_Fitter1::p_sor0, 1.383);
		CkDec_pp.Update();
		CkDec_pL.Update();
		CkDec_LL.Update();
		CkDec_pXim.Update();
		fitter->GoBabyGo();

		ntBuffer[0] = uIter;
		ntBuffer[1] = vFemReg_pp;
		ntBuffer[2] = vFemReg_pL;
		ntBuffer[3] = vFemReg_pXim;
		ntBuffer[4] = vBlReg;
		ntBuffer[5] = vFrac_pp_pL;
		ntBuffer[6] = vFrac_pL_pSigma0;
		ntBuffer[7] = vFrac_pL_pXim;
		ntBuffer[8] = fitter->GetParameter("pp", DLM_Fitter1::p_sor0);
		ntBuffer[9] = fitter->GetParError("pp", DLM_Fitter1::p_sor0);
		ntBuffer[10] = fitter->GetParameter("pp", DLM_Fitter1::p_a);
		ntBuffer[11] = fitter->GetParError("pp", DLM_Fitter1::p_a);
		ntBuffer[12] = fitter->GetParameter("pp", DLM_Fitter1::p_b);
		ntBuffer[13] = fitter->GetParError("pp", DLM_Fitter1::p_b);
		ntBuffer[14] = fitter->GetParameter("pp", DLM_Fitter1::p_Cl);
		ntBuffer[15] = fitter->GetParError("pp", DLM_Fitter1::p_Cl);
		ntBuffer[16] = fitter->GetParameter("pLambda", DLM_Fitter1::p_sor0);
		ntBuffer[17] = fitter->GetParError("pLambda", DLM_Fitter1::p_sor0);
		ntBuffer[18] = fitter->GetParameter("pXim", DLM_Fitter1::p_sor0);
		ntBuffer[19] = fitter->GetParError("pXim", DLM_Fitter1::p_sor0);
		ntBuffer[20] = fitter->GetChi2Ndf();
		ntBuffer[21] = fitter->GetPval();
		ntBuffer[22] = (int) SEPARATE_BL;
		ntBuffer[23] = (int) FIX_CL;
		ntResult->Fill(ntBuffer);

		TFile* GraphFile = new TFile(
				TString::Format("%sGraphFile_%s_Iter%u_uIter%u.root",
						OutputDir.Data(), "CutVarAdd", NumIter, uIter),
				"recreate");

		GraphFile->cd();

		TGraph FitResult_pp;
		FitResult_pp.SetName(TString::Format("FitResult_pp_%u", uIter));
		fitter->GetFitGraph(0, FitResult_pp);

		TGraph FitResult_pL;
		FitResult_pL.SetName(TString::Format("FitResult_pL_%u", uIter));
		fitter->GetFitGraph(1, FitResult_pL);

		TGraph FitResult_LL;
		FitResult_LL.SetName(TString::Format("FitResult_LL_%u", uIter));
		fitter->GetFitGraph(2, FitResult_LL);

		TGraph FitResult_pXim;
		FitResult_pXim.SetName(TString::Format("FitResult_pXim_%u", uIter));
		fitter->GetFitGraph(3, FitResult_pXim);

		double Chi2 = fitter->GetChi2();
		unsigned NDF = fitter->GetNdf();
		double RadiusResult = fitter->GetParameter("pp", DLM_Fitter1::p_sor0);
		double RadiusError = fitter->GetParError("pp", DLM_Fitter1::p_sor0);
		if (Chi2 / double(NDF) != fitter->GetChi2Ndf()) {
			printf("Oh boy...\n");
		}

		double Chi2_pp = 0;
		unsigned EffNumBins_pp = 0;
		for (unsigned uBin = 0; uBin < NumMomBins_pp; uBin++) {

			double mom = AB_pp.GetMomentum(uBin);
			double dataY;
			double dataErr;
			double theoryX;
			double theoryY;

			if (mom > FemtoRegion_pp[vFemReg_pp][1])
				continue;

			FitResult_pp.GetPoint(uBin, theoryX, theoryY);
			if (mom != theoryX) {
				std::cout << mom << '\t' << theoryX << std::endl;
				printf("  PROBLEM pp!\n");
			}
			dataY = OliHisto_pp->GetBinContent(uBin + 1);
			dataErr = OliHisto_pp->GetBinError(uBin + 1);

			Chi2_pp += (dataY - theoryY) * (dataY - theoryY)
					/ (dataErr * dataErr);
			EffNumBins_pp++;
		}

		double Chi2_pL = 0;
		unsigned EffNumBins_pL = 0;
		for (unsigned uBin = 0; uBin < NumMomBins_pL; uBin++) {

			double mom = AB_pL.GetMomentum(uBin);
			//double dataX;
			double dataY;
			double dataErr;
			double theoryX;
			double theoryY;

			if (mom > FemtoRegion_pL[vFemReg_pL][1])
				continue;

			FitResult_pL.GetPoint(uBin, theoryX, theoryY);
			if (mom != theoryX) {
				std::cout << mom << '\t' << theoryX << std::endl;
				printf("  PROBLEM pL!\n");
			}
			dataY = OliHisto_pL->GetBinContent(uBin + 1);
			dataErr = OliHisto_pL->GetBinError(uBin + 1);

			Chi2_pL += (dataY - theoryY) * (dataY - theoryY)
					/ (dataErr * dataErr);
			EffNumBins_pL++;
		}

		double Chi2_LL = 0;
		unsigned EffNumBins_LL = 0;
		for (unsigned uBin = 0; uBin < NumMomBins_LL; uBin++) {

			double mom = Ck_LL->GetBinCenter(uBin);
			//double dataX;
			double dataY;
			double dataErr;
			double theoryX;
			double theoryY;

			if (mom > FemtoRegion_LL[vFemReg_LL][1])
				continue;

			FitResult_LL.GetPoint(uBin, theoryX, theoryY);
			if (mom != theoryX) {
				std::cout << mom << '\t' << theoryX << std::endl;
				printf("  PROBLEM LL!\n");
			}
			dataY = OliHisto_LL->GetBinContent(uBin + 1);
			dataErr = OliHisto_LL->GetBinError(uBin + 1);

			Chi2_LL += (dataY - theoryY) * (dataY - theoryY)
					/ (dataErr * dataErr);
			EffNumBins_LL++;
		}
		double a0_LL = fitter->GetParameter("LambdaLambda",
				DLM_Fitter1::p_pot0);
		double a0_err_LL = fitter->GetParError("LambdaLambda",
				DLM_Fitter1::p_pot0);
		double reff_LL = fitter->GetParameter("LambdaLambda",
				DLM_Fitter1::p_pot1);
		double reff_err_LL = fitter->GetParError("LambdaLambda",
				DLM_Fitter1::p_pot1);

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
		printf("Chi2_pp/bins = %.2f/%u = %.2f\n", Chi2_pp, EffNumBins_pp,
				Chi2_pp / double(EffNumBins_pp));
		printf("Chi2_pL/bins = %.2f/%u = %.2f\n", Chi2_pL, EffNumBins_pL,
				Chi2_pL / double(EffNumBins_pL));
		printf("Chi2_LL/bins = %.2f/%u = %.2f\n", Chi2_LL, EffNumBins_LL,
				Chi2_LL / double(EffNumBins_LL));
		printf("Chi2_pXim/bins = %.2f/%u = %.2f\n", Chi2_pXim, EffNumBins_pXim,
				Chi2_pXim / double(EffNumBins_pXim));

		std::cout << "float radius = "
				<< fitter->GetParameter("pp", DLM_Fitter1::p_sor0) << ";"
				<< std::endl;
		std::cout << "float ppBL0 ="
				<< fitter->GetParameter("pp", DLM_Fitter1::p_a) << ";"
				<< std::endl;
		std::cout << "float ppBL1 ="
				<< fitter->GetParameter("pp", DLM_Fitter1::p_b) << ";"
				<< std::endl;
		std::cout << "float ppNorm ="
				<< fitter->GetParameter("pp", DLM_Fitter1::p_Cl) << ";"
				<< std::endl;

		std::cout << "float pLBL0 = "
				<< fitter->GetParameter("pLambda", DLM_Fitter1::p_a) << ";"
				<< std::endl;
		std::cout << "float pLBL1 = "
				<< fitter->GetParameter("pLambda", DLM_Fitter1::p_b) << ";"
				<< std::endl;

		std::cout << "float LLBL0 = "
				<< fitter->GetParameter("LambdaLambda", DLM_Fitter1::p_a) << ";"
				<< std::endl;
		std::cout << "float LLBL1 = "
				<< fitter->GetParameter("LambdaLambda", DLM_Fitter1::p_b) << ";"
				<< std::endl;

		std::cout << "float pXiBL0 = "
				<< fitter->GetParameter("pXim", DLM_Fitter1::p_a) << ";"
				<< std::endl;
		std::cout << "float pXiBL1 = "
				<< fitter->GetParameter("pXim", DLM_Fitter1::p_b) << ";"
				<< std::endl;
		std::cout << "float pXiNorm = "
				<< fitter->GetParameter("pXim", DLM_Fitter1::p_Cl) << ";"
				<< std::endl;

		GraphFile->cd();
		FitResult_pp.Write();
		FitResult_pL.Write();
		FitResult_LL.Write();
		FitResult_pXim.Write();
		fitter->GetFit()->SetName(TString::Format("GlobalFit_%u", uIter));
		fitter->GetFit()->Write();

		if (FAST_PLOT) {
			TPaveText* info1 = new TPaveText(0.45, 0.65, 0.9, 0.95, "blNDC"); //lbrt
			info1->SetName("info1");
			info1->SetBorderSize(1);
			info1->SetTextSize(0.04);
			info1->SetFillColor(kWhite);
			info1->SetTextFont(22);
			TString SOURCE_NAME;
			SOURCE_NAME = "Gauss";
			info1->AddText(
					TString::Format("R(%s)=%.3f#pm%.3f", SOURCE_NAME.Data(),
							RadiusResult, RadiusError));
			info1->AddText(
					TString::Format("C_{l}=%.3f#pm%.3f",
							fitter->GetParameter("pp", DLM_Fitter1::p_Cl),
							fitter->GetParError("pp", DLM_Fitter1::p_Cl)));
			info1->AddText(
					TString::Format("p_a = %.3f p_b = %.5f",
							fitter->GetParameter("pp", DLM_Fitter1::p_a),
							fitter->GetParameter("pp", DLM_Fitter1::p_b)));

			info1->AddText(
					TString::Format(
							"Global #chi^{2}_{ndf}=%.0f/%u=%.2f, p_{val}=%.3f",
							Chi2, NDF, Chi2 / double(NDF),
							TMath::Prob(Chi2, round(NDF))));
			info1->AddText(
					TString::Format(
							"Local #chi^{2}/ndf=%.0f/%u=%.2f, p_{val}=%.3f",
							Chi2_pp, EffNumBins_pp,
							Chi2_pp / double(EffNumBins_pp),
							TMath::Prob(Chi2_pp, round(EffNumBins_pp))));

			TPaveText* info2 = new TPaveText(0.45, 0.63, 0.9, 0.95, "blNDC"); //lbrt
			info2->SetName("info2");
			info2->SetBorderSize(1);
			info2->SetTextSize(0.04);
			info2->SetFillColor(kWhite);
			info2->SetTextFont(22);
			SOURCE_NAME = "Gauss";
			info2->AddText(
					TString::Format("R(%s)=%.3f#pm%.3f", SOURCE_NAME.Data(),
							fitter->GetParameter("pLambda",
									DLM_Fitter1::p_sor0),
							fitter->GetParError("pLambda",
									DLM_Fitter1::p_sor0)));
			info2->AddText(
					TString::Format("C_{l}=%.3f#pm%.3f",
							fitter->GetParameter("pLambda", DLM_Fitter1::p_Cl),
							fitter->GetParError("pLambda", DLM_Fitter1::p_Cl)));
			info2->AddText(
					TString::Format(
							"Local #chi^{2}/ndf=%.0f/%u=%.2f, p_{val}=%.3f",
							Chi2_pL, EffNumBins_pL,
							Chi2_pL / double(EffNumBins_pL),
							TMath::Prob(Chi2_pL, round(EffNumBins_pL))));

			TPaveText* info3 = new TPaveText(0.45, 0.63, 0.9, 0.95, "blNDC"); //lbrt
			info3->SetName("info3");
			info3->SetBorderSize(1);
			info3->SetTextSize(0.04);
			info3->SetFillColor(kWhite);
			info3->SetTextFont(22);
			SOURCE_NAME = "Gauss";
			info3->AddText(
					TString::Format("R(%s)=%.3f#pm%.3f", SOURCE_NAME.Data(),
							fitter->GetParameter("LambdaLambda",
									DLM_Fitter1::p_sor0),
							fitter->GetParError("LambdaLambda",
									DLM_Fitter1::p_sor0)));
			info3->AddText(
					TString::Format("a_{0}=%.2f#pm%.2f", a0_LL, a0_err_LL));
			info3->AddText(
					TString::Format("r_{eff}=%.2f#pm%.2f", reff_LL,
							reff_err_LL));
			info3->AddText(
					TString::Format("C_{l}=%.2f#pm%.2f",
							fitter->GetParameter("LambdaLambda",
									DLM_Fitter1::p_Cl),
							fitter->GetParError("LambdaLambda",
									DLM_Fitter1::p_Cl)));
			info3->AddText(
					TString::Format(
							"Local #chi^{2}/ndf=%.0f/%u=%.2f, p_{val}=%.3f",
							Chi2_LL, EffNumBins_LL,
							Chi2_LL / double(EffNumBins_LL),
							TMath::Prob(Chi2_LL, round(EffNumBins_LL))));

			TPaveText* info4 = new TPaveText(0.2, 0.68, 0.9, 0.95, "blNDC"); //lbrt
			info4->SetName("info4");
			info4->SetBorderSize(1);
			info4->SetTextSize(0.04);
			info4->SetFillColor(kWhite);
			info4->SetTextFont(22);
			SOURCE_NAME = "Gauss";
			info4->AddText(
					TString::Format("R(%s)=%.3f#pm%.3f", SOURCE_NAME.Data(),
							fitter->GetParameter("pXim", DLM_Fitter1::p_sor0),
							fitter->GetParError("pXim", DLM_Fitter1::p_sor0)));
			info4->AddText(
					TString::Format("C_{l}=%.3f#pm%.3f",
							fitter->GetParameter("pXim", DLM_Fitter1::p_Cl),
							fitter->GetParError("pXim", DLM_Fitter1::p_Cl)));

			double Yoffset = 1.2;
			TH1F* hAxis_pp = new TH1F("hAxis_pp", "hAxis_pp", 600, 0, 600);
			hAxis_pp->SetStats(false);
			hAxis_pp->SetTitle("");
			hAxis_pp->GetXaxis()->SetLabelSize(0.065);
			hAxis_pp->GetXaxis()->CenterTitle();
			hAxis_pp->GetXaxis()->SetTitleOffset(1.35);
			hAxis_pp->GetXaxis()->SetLabelOffset(0.02);
			hAxis_pp->GetXaxis()->SetTitleSize(0.075);
			hAxis_pp->GetYaxis()->SetLabelSize(0.065);
			hAxis_pp->GetYaxis()->CenterTitle();
			hAxis_pp->GetYaxis()->SetTitleOffset(Yoffset);
			hAxis_pp->GetYaxis()->SetTitleSize(0.075);
			hAxis_pp->GetXaxis()->SetRangeUser(0, 600);
			hAxis_pp->GetYaxis()->SetRangeUser(0.5, 3);
			TF1* blPP = new TF1("blPP", "pol1", 0, 500);
			blPP->SetParameters(fitter->GetParameter("pp", DLM_Fitter1::p_a),
					fitter->GetParameter("pp", DLM_Fitter1::p_b));
			blPP->SetLineColor(7);
			blPP->SetLineWidth(4);

			TH1F* hAxis_pL = new TH1F("hAxis_pL", "hAxis_pL", 600, 0, 600);
			hAxis_pL->SetStats(false);
			hAxis_pL->SetTitle("");
			hAxis_pL->GetXaxis()->SetLabelSize(0.065);
			hAxis_pL->GetXaxis()->CenterTitle();
			hAxis_pL->GetXaxis()->SetTitleOffset(1.35);
			hAxis_pL->GetXaxis()->SetLabelOffset(0.02);
			hAxis_pL->GetXaxis()->SetTitleSize(0.075);
			hAxis_pL->GetYaxis()->SetLabelSize(0.065);
			hAxis_pL->GetYaxis()->CenterTitle();
			hAxis_pL->GetYaxis()->SetTitleOffset(Yoffset);
			hAxis_pL->GetYaxis()->SetTitleSize(0.075);
			hAxis_pL->GetXaxis()->SetRangeUser(0, kMax_pL);
			hAxis_pL->GetYaxis()->SetRangeUser(0.8, 2.0);    //pPb

			TH1F* hAxis_LL = new TH1F("hAxis_LL", "hAxis_LL", 600, 0, 600);
			hAxis_LL->SetStats(false);
			hAxis_LL->SetTitle("");
			hAxis_LL->GetXaxis()->SetLabelSize(0.065);
			//hAxis_LL->GetXaxis()->SetTitle("k (MeV)");
			hAxis_LL->GetXaxis()->CenterTitle();
			hAxis_LL->GetXaxis()->SetTitleOffset(1.35);
			hAxis_LL->GetXaxis()->SetLabelOffset(0.02);
			hAxis_LL->GetXaxis()->SetTitleSize(0.075);
			hAxis_LL->GetYaxis()->SetLabelSize(0.065);
			//hAxis_LL->GetYaxis()->SetTitle("C(k)");
			hAxis_LL->GetYaxis()->CenterTitle();
			hAxis_LL->GetYaxis()->SetTitleOffset(Yoffset);
			hAxis_LL->GetYaxis()->SetTitleSize(0.075);
			//hAxis_LL->GetXaxis()->SetRangeUser(0,kMax_LL);
			hAxis_LL->GetYaxis()->SetRangeUser(0.5, 1.4);

			TH1F* hAxis_pXim = new TH1F("hAxis_pXim", "hAxis_pXim", 600, 0,
					600);
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
			hAxis_pXim->GetXaxis()->SetRangeUser(0, kMax_pXim);
			hAxis_pXim->GetYaxis()->SetRangeUser(0.7, 4.5);    //pPb
			TString CanName = Form("cfast");
//					TString CanName = Form("cFastBL_%u_%u.root",
//							(int) BlRegion_pp[1][0],
//							(int) BlRegion_pp[1][1]);
			TCanvas* cfast = new TCanvas(CanName.Data(), CanName.Data(), 1);
			cfast->cd(0);
			cfast->SetCanvasSize(1920, 1280);
			cfast->SetMargin(0.15, 0.05, 0.2, 0.05);    //lrbt

			cfast->cd(1);
			OliHisto_pp->SetStats(false);
			OliHisto_pp->SetTitle("pp");
			OliHisto_pp->SetLineWidth(2);
			OliHisto_pp->SetLineColor(kBlack);
			FitResult_pp.SetLineWidth(2);
			FitResult_pp.SetLineColor(kRed);
			FitResult_pp.SetMarkerStyle(24);
			FitResult_pp.SetMarkerColor(kRed);
			FitResult_pp.SetMarkerSize(1);

			hAxis_pp->Draw("axis");
			OliHisto_pp->Draw("same");
			blPP->Draw("same");
			FitResult_pp.Draw("CP,same");
			info1->Draw("same");

			cfast->Write();
			cfast->SaveAs(
					TString::Format("%sIter%u_uIter%u.png", OutputDir.Data(),
							NumIter, uIter));

			delete info1;
			delete info2;
			delete info3;
			delete info4;
			delete cfast;

		}      //FAST_PLOT

		delete fitter;
		delete Ck_pp;
		delete Ck_pL;
		delete Ck_LL;
		delete Ck_pSigma0;
		delete Ck_pXim;
		delete Ck_pXim1530;

		if (OliFile_pp) {
			delete OliFile_pp;
			OliFile_pp = NULL;
		}
		if (OliFile_pL) {
			delete OliFile_pL;
			OliFile_pL = NULL;
		}
		if (OliFile_LL) {
			delete OliFile_LL;
			OliFile_LL = NULL;
		}
		if (OliFile_pXim) {
			delete OliFile_pXim;
			OliFile_pXim = NULL;
		}
		if (SystErrFile_pp) {
			delete SystErrFile_pp;
			SystErrFile_pp = NULL;
		}
		if (SystErrFile_pL) {
			delete SystErrFile_pL;
			SystErrFile_pL = NULL;
		}
		if (SystErrFile_LL) {
			delete SystErrFile_LL;
			SystErrFile_LL = NULL;
		}
		if (SystErrFile_pXim) {
			delete SystErrFile_pXim;
			SystErrFile_pXim = NULL;
		}

		delete GraphFile;

	}      //END OF THE BIG FOR LOOP OVER ITER
	OutFile->cd();
	ntResult->Write();

	delete ntResult;
	delete OutFile;

	for (unsigned uVar = 0; uVar < 3; uVar++) {
		delete[] Purities_p[uVar];
		delete[] Fraction_p[uVar];
		delete[] Purities_L[uVar];
		delete[] Fraction_L[uVar];
	}
	delete[] Purities_p;
	delete[] Fraction_p;
	delete[] Purities_L;
	delete[] Fraction_L;

}


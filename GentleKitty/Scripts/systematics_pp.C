//#include <XiOnly.C>
#include <iostream>

#include "CATStools.h"
#include "CATS.h"
#include "DLM_Potentials.h"
#include "DLM_Source.h"

#include "TGraph.h"
#include "TFile.h"
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

void RUN2_SYSTEMATICS_MEDIAN(const char* InputFolder, const char* OutDirName) {
	const char* FileName = Form("%sOutFile_CutVarAdd_2048.root", InputFolder);
	TFile* inFile = TFile::Open(FileName, "READ");
	TNtuple* sysVarTree = (TNtuple*) inFile->Get("ntResult");
	if (!sysVarTree) {
		std::cout << "no Tree loaded\n";
	}
	auto histRad = new TH1D("hRad", "hRad", 10000, 1.2, 1.6);
	auto histRadIter = new TH2D("hRadIter", "hRadIter", 10000, 1.2, 1.6, 2000,
			0, 2000);
	sysVarTree->Draw("Radius_pp>>hRad","Chi2NdfLocal<6");
	auto mean = histRad->GetMean();

	int n = histRad->GetXaxis()->GetNbins();
	std::vector<double> x(n);
	histRad->GetXaxis()->GetCenter(&x[0]);
	const double * y = histRad->GetArray();
	// exclude underflow/overflows from bin content array y
	auto median = TMath::Median(n, &x[0], &y[1]);
	auto medianBin = histRad->FindBin(median);
	std::cout << "Median: " << median << "\t Mean: " << mean << std::endl;

	auto histRadCumulative = histRad->GetCumulative();
	auto canRad = new TCanvas("cRad", "cRad", 1200, 800);
	canRad->Divide(3, 3);
	canRad->cd(1);
	histRad->Draw();
	canRad->cd(2);
	histRadCumulative->Scale(1 / (double) histRad->GetEntries());
	histRadCumulative->Draw();

	int binMin = 0;
	int binMax = 0;
	for (int iBin = 0; iBin < histRadCumulative->GetNbinsX(); iBin++) {
		if (binMin == 0
				&& (histRadCumulative->GetBinContent(iBin)
						> histRadCumulative->GetBinContent(medianBin) - 0.32)) {
			binMin = iBin;
		}
		if (binMax == 0
				&& (histRadCumulative->GetBinContent(iBin)
						> histRadCumulative->GetBinContent(medianBin) + 0.32)) {
			binMax = iBin;
			break;
		}
	}
	auto radMin = histRadCumulative->GetXaxis()->GetBinCenter(binMin);
	auto radMax = histRadCumulative->GetXaxis()->GetBinCenter(binMax);

	canRad->cd(3);
	sysVarTree->Draw("IterID:Radius_pp>>hRadIter");
	canRad->cd(4);
	auto hIterMean = histRadIter->ProjectionY("hIterMean",
			histRadIter->GetXaxis()->FindBin(mean),
			histRadIter->GetXaxis()->FindBin(mean));
	hIterMean->Draw();
	int uIterMean;
	hIterMean->GetBinWithContent(1., uIterMean);
	std::cout << "uIterMean: " << uIterMean << std::endl;

	canRad->cd(5);
	auto hIterMedian = histRadIter->ProjectionY("hIterMedian",
			histRadIter->GetXaxis()->FindBin(median),
			histRadIter->GetXaxis()->FindBin(median));
	hIterMedian->Draw();
	int uIterMedian;
	hIterMedian->GetBinWithContent(1., uIterMedian);
	std::cout << "uIterMedian: " << uIterMedian << std::endl;

	canRad->cd(6);
	auto hIterUp = histRadIter->ProjectionY("hIterUp",
			histRadIter->GetXaxis()->FindBin(radMax),
			histRadIter->GetXaxis()->FindBin(radMax));
	hIterUp->Draw();
	int uIterUp;
	hIterUp->GetBinWithContent(1., uIterUp);
	std::cout << "uIterUp: " << uIterUp << std::endl;

	canRad->cd(7);
	auto hIterLow = histRadIter->ProjectionY("hIterLow",
			histRadIter->GetXaxis()->FindBin(radMin),
			histRadIter->GetXaxis()->FindBin(radMin));
	hIterLow->Draw();
	int uIterLow;
	hIterLow->GetBinWithContent(1., uIterLow);
	std::cout << "uIterLow: " << uIterLow << std::endl;

	canRad->SaveAs(Form("%scanRad.pdf", OutDirName));
	TFile* outFile = TFile::Open(
			Form("%sSYSTEMATICS_CutVarAdd_Global_Radius_Normal.root",
					OutDirName), "RECREATE");

	TFile* fileMean = TFile::Open(
			Form("%sGraphFile_CutVarAdd_Iter2048_uIter%i.root", InputFolder,
					uIterMean));
	if (!fileMean) {
		std::cout << "Missing file Mean \n";
	}
	TGraph* FileGrDefault_pp = (TGraph*) fileMean->Get(
			Form("FitResult_pp_%u", uIterMean));
	TGraph GrDefault_pp(*FileGrDefault_pp);
	GrDefault_pp.SetName("ppGraphDefault");
	outFile->cd();
	GrDefault_pp.Write();
	fileMean->Close();

	TFile* fileSysLow = TFile::Open(
			Form("%sGraphFile_CutVarAdd_Iter2048_uIter%i.root", InputFolder,
					uIterLow));
	if (!fileSysLow) {
		std::cout << "Missing file Sys Low \n";
	}
	TGraph* FileGrLow_pp = (TGraph*) fileSysLow->Get(
			Form("FitResult_pp_%u", uIterLow));
	TGraph GrLow_pp(*FileGrLow_pp);
	GrLow_pp.SetName("ppGraphLowerLim");
	outFile->cd();
	GrLow_pp.Write();
	fileSysLow->Close();

	TFile* fileSysUp = TFile::Open(
			Form("%sGraphFile_CutVarAdd_Iter2048_uIter%i.root", InputFolder,
					uIterUp));
	if (!fileSysUp) {
		std::cout << "Missing file SysUp \n";
	}
	TGraph* FileGrUp_pp = (TGraph*) fileSysUp->Get(
			Form("FitResult_pp_%u", uIterUp));
	TGraph GrUp_pp(*FileGrUp_pp);
	GrUp_pp.SetName("ppGraphUpperLim");
	outFile->cd();
	GrUp_pp.Write();
	fileSysUp->Close();

	float uIterIDDefault;
	float rDefault_pp;
	float rErr_pp;
	sysVarTree->SetBranchAddress("IterID", &uIterIDDefault);
	sysVarTree->SetBranchAddress("Radius_pp", &rDefault_pp);
	sysVarTree->SetBranchAddress("RadiusErr_pp", &rErr_pp);
	for (int iEntry = 0; iEntry < sysVarTree->GetEntries(); iEntry++) {
		sysVarTree->GetEntry(iEntry);
		if (uIterIDDefault == 0) {
			break;
		}
	}

	auto errLow = rDefault_pp - radMin;
	auto errUp = radMax - rDefault_pp;

	float rLower = rDefault_pp
			- TMath::Sqrt(
					rErr_pp * rErr_pp
							+ (0.2 * rDefault_pp) * (0.2 * rDefault_pp)
							+ errLow * errLow);
	float rUp = rDefault_pp + TMath::Sqrt(rErr_pp * rErr_pp + errLow * errLow);

	TNtuple* outTuple = new TNtuple("outTuple", "outTuple",
			"Rad_pp:RadStat_pp:RadSysLow_pp:RadSysUp_pp:RadLow_pXi:RadUp_pXi");
	float outArray[6];
	outArray[0] = rDefault_pp;
	outArray[1] = rErr_pp;
	outArray[2] = errLow;
	outArray[3] = errUp;
	outArray[4] = rLower;
	outArray[5] = rUp;
	outTuple->Fill(outArray);
	outFile->cd();
	outTuple->Write();

//	GetXiForRadius("~/cernbox/pPb/v0offlineFix/woDetadPhi/200_400/", OutDirName, rLower, 11,
//			outFile, "UpperLim", true);
//	GetXiForRadius("~/cernbox/pPb/v0offlineFix/woDetadPhi/200_400/", OutDirName, rDefault_pp, 12,
//			outFile, "Default", true);
//	GetXiForRadius("~/cernbox/pPb/v0offlineFix/woDetadPhi/200_400/", OutDirName, rUp, 13, outFile,
//			"LowerLim", true);

	outFile->Close();
}

int main(int argc, char *argv[]) {
	RUN2_SYSTEMATICS_MEDIAN(argv[1], argv[2]);
	return 0;
}

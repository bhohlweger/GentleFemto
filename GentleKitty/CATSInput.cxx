/*
 * CATSInput.cxx
 *
 *  Created on: Oct 30, 2018
 *      Author: hohlweger
 */

#include <CATSInput.h>
#include "TFile.h"
#include <iostream>
CATSInput::CATSInput(int Fraction_Res, int Fraction_Sig, double UnitConv_Res,
		double UnitConv_Sig) :
		fNameBasedir(), fNameMomResFile(), fNameSigmaFile(), fFraction_Res(
				Fraction_Res), fFraction_Sig(Fraction_Sig), fUnitConv_Res(
				UnitConv_Res), fUnitConv_Sig(UnitConv_Sig), fRes(), fSigma() {
	// TODO Auto-generated constructor stub

}

CATSInput::~CATSInput() {
	// TODO Auto-generated destructor stub
}

void CATSInput::ReadResFile() {
	TString ResMatrixFileName = "";
	ResMatrixFileName += fNameBasedir;
	ResMatrixFileName += fNameMomResFile;
	auto FileRes = new TFile(ResMatrixFileName, "read");
	if (!FileRes) {
		std::cout << "No Resolution file found in " << ResMatrixFileName
				<< std::endl;
	} else {
		FileRes->cd();
		std::vector<TH2F*> inputHist;
		std::vector<const char*> FileNames = { "hRes_pp_pL",
				"hRes_pL_pSigma0", "hRes_pL_pXim",
				"hRes_pXim_pXim1530" };
		for (auto it : FileNames) {
			inputHist.push_back((TH2F*) FileRes->Get(it));
		}
		for (auto it : inputHist) {
			TString histName = Form("%s_MeV", it->GetName());
			TH2F* tmpResMeV = new TH2F(histName.Data(), histName.Data(),
					it->GetNbinsX() / fFraction_Res,
					it->GetXaxis()->GetBinLowEdge(1) * fUnitConv_Res,
					it->GetXaxis()->GetBinUpEdge(
							it->GetNbinsX() / fFraction_Res) * fUnitConv_Res,
					it->GetNbinsY() / fFraction_Res,
					it->GetYaxis()->GetBinLowEdge(1) * fUnitConv_Res,
					it->GetXaxis()->GetBinUpEdge(
							it->GetNbinsY() / fFraction_Res) * fUnitConv_Res);
			for (int iBinX = 1; iBinX <= it->GetNbinsX() / fFraction_Res;
					iBinX++) {
				for (int iBinY = 1; iBinY <= it->GetNbinsY() / fFraction_Res;
						iBinY++) {
					tmpResMeV->SetBinContent(iBinX, iBinY,
							it->GetBinContent(iBinX, iBinY));
				}
			}
			fRes.push_back(tmpResMeV);
		}
	}
}

void CATSInput::ReadSigmaFile() {
	TString SigmaMatrixFileName = "";
	SigmaMatrixFileName += fNameBasedir;
	SigmaMatrixFileName += fNameSigmaFile;
	auto FileSigma = new TFile(SigmaMatrixFileName, "read");
	if (!FileSigma) {
		std::cout << "No Sigma file found in " << SigmaMatrixFileName
				<< std::endl;
	} else {
		FileSigma->cd();
		std::vector<TH2F*> inputHist;
		std::vector<const char*> FileNames = { "hSigmaMeV_Proton_Proton",
				"hSigmaMeV_Proton_Lambda", "hSigmaMeV_Lambda_Lambda",
				"hSigmaMeV_Proton_Xim" };
		for (auto it : FileNames) {
			inputHist.push_back((TH2F*) FileSigma->Get(it));
		}
		for (auto it : inputHist) {
			TString histName = Form("%s_MeV", it->GetName());
			TH2F* tmpSigmaMeV = new TH2F(histName.Data(), histName.Data(),
					it->GetNbinsX() / fFraction_Sig,
					it->GetXaxis()->GetBinLowEdge(1) * fUnitConv_Sig,
					it->GetXaxis()->GetBinUpEdge(
							it->GetNbinsX() / fFraction_Sig) * fUnitConv_Sig,
					it->GetNbinsY() / fFraction_Sig,
					it->GetYaxis()->GetBinLowEdge(1) * fUnitConv_Sig,
					it->GetXaxis()->GetBinUpEdge(
							it->GetNbinsY() / fFraction_Sig) * fUnitConv_Sig);
			for (int iBinX = 1; iBinX <= it->GetNbinsX() / fFraction_Sig;
					iBinX++) {
				for (int iBinY = 1; iBinY <= it->GetNbinsY() / fFraction_Sig;
						iBinY++) {
					tmpSigmaMeV->SetBinContent(iBinX, iBinY,
							it->GetBinContent(iBinX, iBinY));
				}
			}
			fSigma.push_back(tmpSigmaMeV);
		}
	}
}

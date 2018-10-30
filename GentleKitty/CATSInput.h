/*
 * CATSInput.h
 *
 *  Created on: Oct 30, 2018
 *      Author: hohlweger
 */

#ifndef GENTLEKITTY_CATSINPUT_H_
#define GENTLEKITTY_CATSINPUT_H_
#include "TString.h"
#include "TH2F.h"
#include <vector>
class CATSInput {
public:

	CATSInput(int Fraction_Res = 2, int Fraction_Sig = 1, double UnitConv_Res =
			1, double UnitConv_Sig = 1);
	virtual ~CATSInput();
	void SetBaseDir(const char* path) {
		fNameBasedir.Clear();
		fNameBasedir.Append(path);
	}
	;
	void SetMomResFileName(const char* filename) {
		fNameMomResFile.Clear();
		fNameMomResFile.Append(filename);
	}
	;
	void ReadResFile();
	TH2F* GetResFile(int iPair) const {return fRes.size()>iPair?fRes[iPair]:nullptr;};
	void SetSigmaFileName(const char* filename) {
		fNameSigmaFile.Clear();
		fNameSigmaFile.Append(filename);
	}
	;
	void ReadSigmaFile();
	TH2F* GetSigmaFile(int iPair) const {return fSigma.size()>iPair?fSigma[iPair]:nullptr;};
private:
	TString fNameBasedir;
	TString fNameMomResFile;
	TString fNameSigmaFile;
	const int fFraction_Res;
	const int fFraction_Sig;
	const double fUnitConv_Res;
	const double fUnitConv_Sig;
	std::vector<TH2F*> fRes;
	std::vector<TH2F*> fSigma;
};

#endif /* GENTLEKITTY_CATSINPUT_H_ */

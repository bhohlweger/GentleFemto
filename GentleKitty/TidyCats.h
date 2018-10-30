/*
 * TidyCats.h
 *
 *  Created on: Oct 30, 2018
 *      Author: hohlweger
 */
#include "CATS.h"
#include "DLM_Potentials.h"
#include "DLM_Source.h"
#ifndef GENTLEKITTY_TIDYCATS_H_
#define GENTLEKITTY_TIDYCATS_H_

class TidyCats {
public:
	TidyCats();
	virtual ~TidyCats();
	void GetCatsProtonProton(CATS* AB_pp, double GaussSourceSize, double* pars,
			int momBins, double kMin, double kMax);
	void GetCatsProtonLambda(CATS* AB_pL, double GaussSourceSize, double* pars,
			int momBins, double kMin, double kMax);
	void GetCatsProtonXiMinus(CATS* AB_pXim, double GaussSourceSize,
			double* pars, int momBins, double kMin, double kMax, bool StrongOn,
			double QCDTime);
	void GetCatsProtonXiMinus1530(CATS* AB_pXim1530, double GaussSourceSize,
			double* pars, int momBins, double kMin, double kMax);
private:
};

#endif /* GENTLEKITTY_TIDYCATS_H_ */

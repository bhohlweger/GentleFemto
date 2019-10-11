/*
 * DreamHEP.cxx
 *
 *  Created on: Mar 12, 2019
 *      Author: schmollweger
 */

#include "DreamHEP.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "TString.h"
DreamHEP::DreamHEP() :
  fRootS(7000.0),
  fkStarMax(999.) {
  // TODO Auto-generated constructor stub

}

DreamHEP::~DreamHEP() {
  // TODO Auto-generated destructor stub
}


void DreamHEP::printTH1HEPdata(const TH1* hist, const TGraphErrors* syst, const char* outname) {
  std::ofstream output;
  output.open (Form("%s.yaml", outname));
  output << "dependent_variables:\n";
  output << "- header: {name: C(k*)}\n";
  output << "  qualifiers:\n";
  output << "  - {name: SQRT(S), units: GeV, value: '";
  output << fRootS << "'}\n";
  output << "  values:\n";

  double x,y;
  bool isMeV = false;
  for(int i=0; i<syst->GetN(); ++i) {
    syst->GetPoint(i, x, y);
    if( i == 0 && x > 1) {
      isMeV = true;
    }
    if(x > fkStarMax) continue;
    output << "  - errors:\n";
    output << "    - {label: stat, symerror: " << hist->GetBinError(i+1) << "}\n";
    output << "    - {label: sys, symerror: " << syst->GetErrorY(i) << "}\n";
    output << "    value: " << y << "\n";
  }

  TString unit = isMeV ? "(MeV/c)" : "(GeV/c)";
  output << "independent_variables:\n";
  output << "- header: {name: k* ";
  output << unit.Data();
  output << "}\n";
  output << "  values: \n";

  for(int i=0; i<syst->GetN(); ++i) {
    syst->GetPoint(i, x, y);
    if(x > fkStarMax) continue;
    output << "  - {value: " << x << "}\n";
  }
  output.close();
}

void DreamHEP::printTGAsymmHEPdata(const TGraphAsymmErrors* gr, const TGraphErrors* syst, const char* outname) {
  std::ofstream output;
  output.open (Form("%s.yaml", outname));
  output << "dependent_variables:\n";
  output << "- header: {name: C(k*)}\n";
  output << "  qualifiers:\n";
  output << "  - {name: SQRT(S), units: GeV, value: '";
  output << fRootS << "'}\n";
  output << "  values:\n";

  double x,y;
  bool isMeV = false;
  for(int i=0; i<syst->GetN(); ++i) {
    syst->GetPoint(i, x, y);
    if( i == 0 && x > 1) {
      isMeV = true;
    }
    if(x > fkStarMax) continue;
    output << "  - errors:\n";
    output << "    - {label: stat, symerror: " << gr->GetErrorY(i) << "}\n";
    output << "    - {label: sys, symerror: " << syst->GetErrorY(i) << "}\n";
    output << "    value: " << y << "\n";
  }

  TString unit = isMeV ? "(MeV/c)" : "(GeV/c)";
  output << "independent_variables:\n";
  output << "- header: {name: k* ";
  output << unit.Data();
  output << "}\n";
  output << "  values: \n";

  for(int i=0; i<gr->GetN(); ++i) {
    gr->GetPoint(i, x, y);
    if(x > fkStarMax) continue;
    output << "  - value: " << x << "\n";
  }
  output.close();
}

TGraphErrors* DreamHEP::GetSystErrHist(TH1* hist, TF1* syst) {
  TString HistName = hist->GetName();
  TGraphErrors* HistWithSysErr = new TGraphErrors();
  HistWithSysErr->SetName(Form("%s_syst", HistName.Data()));
  const double convFac = HistName.Contains("MeV") ? 1000. : 1.;
  std::cout << "Conversion factor correct? " << convFac << std::endl;
  int NumSEB = hist->FindBin(0.5 * convFac);
  for (int iBin = 0; iBin < NumSEB; iBin++) {
    const float x = hist->GetBinCenter(iBin + 1);
    const float y = hist->GetBinContent(iBin + 1);
    HistWithSysErr->SetPoint(iBin, x, y);
    HistWithSysErr->SetPointError(iBin, 0, y * syst->Eval(x / convFac));
  }
  return HistWithSysErr;
}

TGraphErrors* DreamHEP::GetSystErrHist(TGraphAsymmErrors* gr, TF1* syst) {
  TString HistName = gr->GetName();
  TGraphErrors* HistWithSysErr = new TGraphErrors();
  HistWithSysErr->SetName(Form("%s_syst", HistName.Data()));
  const double convFac = HistName.Contains("MeV") ? 1000. : 1.;
  std::cout << "Conversion factor correct? " << convFac << std::endl;
  double x, y;
  for (int i = 0; i < gr->GetN(); ++i) {
    gr->GetPoint(i, x, y);
    HistWithSysErr->SetPoint(i, x, y);
    HistWithSysErr->SetPointError(i, 0, y * syst->Eval(x / (float) convFac));
  }
  return HistWithSysErr;
}


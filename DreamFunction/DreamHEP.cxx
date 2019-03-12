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
DreamHEP::DreamHEP() {
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
  output << "  - {name: SQRT(S), units: GeV, value: '7000.0'}\n";
  output << "  values:\n";

  double x,y;
  for(int i=0; i<syst->GetN(); ++i) {
    syst->GetPoint(i, x, y);
    output << "  - errors:\n";
    output << "    - {label: stat, symerror: " << hist->GetBinError(i+1) << "}\n";
    output << "    - {label: sys, symerror: " << syst->GetErrorY(i) << "}\n";
    output << "    value: " << y << "\n";
  }

  output << "independent_variables:\n";
  output << "- header: {name: k* (GeV/c)}\n";
  output << "  values: \n";

  for(int i=0; i<syst->GetN(); ++i) {
    syst->GetPoint(i, x, y);
    output << "  - {value: " << x << "}\n";
  }
  output.close();
}

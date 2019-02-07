/*
 * DreamdEtadPhi.h
 *
 *  Created on: Feb 4, 2019
 *      Author: schmollweger
 */

#ifndef DREAMFUNCTION_DREAMDETADPHI_H_
#define DREAMFUNCTION_DREAMDETADPHI_H_
#include "TH2F.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TList.h"
#include "TPad.h"

class DreamdEtadPhi {
 public:
  DreamdEtadPhi();
  virtual ~DreamdEtadPhi();
  void SetSEDistribution(TH2F* SE, const char* name) {
    fSEdEtadPhi = (TH2F*)SE->Clone(Form("%s%s",SE->GetName(),name));
  }
  void AddSEDistribution(TH2F* SEAdd) {
    fSEdEtadPhi->Add(SEAdd);
  }
  void SetMEDistribution(TH2F* ME, const char* name) {
    fMEdEtadPhi = (TH2F*)ME->Clone(Form("%s%s",ME->GetName(),name));
  }
  void AddMEDistribution(TH2F* MEAdd) {
    fMEdEtadPhi->Add(MEAdd);
  }
  void ShiftAbovePhi();
  void DivideSEandME(int rebin = 2);
  void ProjectionY();
  void Draw2D(TPad *p, float Rad);
  void DrawProjectionY(TPad *p, float Rad);
  void WriteOutput(TList* output, const char *outname) {
    TList *outList = new TList();
    outList->SetName(Form("%s",outname));
    outList->SetOwner();
    output->Add(outList);
    if (fSEdEtadPhi) {
      outList->Add(fSEdEtadPhi);
    }
    if (fMEdEtadPhi) {
      outList->Add(fMEdEtadPhi);
    }
    if (fdEtadPhi) {
      outList->Add(fdEtadPhi);
    }
    if (fProjectionY) {
      outList->Add(fProjectionY);
    }
  }
 private:
  TH2F* fSEdEtadPhi;
  TH2F* fMEdEtadPhi;
  TH2F* fdEtadPhi;
  TH1D* fProjectionY;
};

#endif /* DREAMFUNCTION_DREAMDETADPHI_H_ */

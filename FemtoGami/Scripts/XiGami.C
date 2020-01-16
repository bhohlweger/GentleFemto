#include "AnalyseProXi.h"
#include "CATSInput.h"
#include "CATSLambdaParam.h"
#include "DreamSystematics.h"
#include "LambdaGami.h"
#include "SideBandFit.h"
#include "TidyCats.h"

#include "TF1.h"
#include "TSystem.h"

#include "DLM_CkDecomposition.h"

static double cutOff = 1000;  // at this value the calculation and doing of the cf stops

int main(int argc, char *argv[]) {
  DreamPlot::SetStyle();
  const char* fileName = argv[1];
  const char* prefix = argv[2];
  AnalyseProXi* ana = new AnalyseProXi(cutOff, 0.997);
  ana->LimitCFRange(0.999);

  DreamSystematics protonNorm(DreamSystematics::pXiNorm);
  protonNorm.SetUpperFitRange(800);
  protonNorm.SetBarlowUpperRange(400);

  DreamSystematics protonFeeddown(DreamSystematics::pXiLam);
  protonFeeddown.SetUpperFitRange(800);
  protonFeeddown.SetBarlowUpperRange(400);

  DreamSystematics protonMomRes(DreamSystematics::pXiRes);
  protonMomRes.SetUpperFitRange(800);
  protonMomRes.SetBarlowUpperRange(400);

  ana->SetAnalysisFile(fileName, prefix);
  ana->Default();
  TH1F* DefVar = ana->GetVariation(0, true);
  if (!DefVar) {
    return -1;
  }

  protonNorm.SetDefaultHist((TH1F*) DefVar->Clone("defVarNorm"));
  protonFeeddown.SetDefaultHist((TH1F*) DefVar->Clone("defVarFeedDown"));
  protonMomRes.SetDefaultHist((TH1F*) DefVar->Clone("defVarMom"));

  //create dream sys object, set default for norm/bl vars
  int varCounter = 1;
  for (int iNorm = 0; iNorm < 10; ++iNorm) {
    for (int iBLfnct = 0; iBLfnct < 2; ++iBLfnct) {
      if (iNorm == 2 && iBLfnct == 0) {
        continue;  //default case
      }
      ana->Default();
      ana->SetNormVar(iNorm);
      ana->SetBaselineVar(iBLfnct);
      TH1F* Var = ana->GetVariation(varCounter++);
      protonNorm.SetVarHist(Var);

    }
  }
  protonNorm.EvalSystematics();
  protonNorm.WriteOutput("");

  //create dream sys object, set default for lam/sideband/pxi rad vars
  for (int iSBNorm = 0; iSBNorm < 3; ++iSBNorm) {
    if (iSBNorm == 1) {
      continue;  // default case
    }
    ana->Default();
    ana->SetSideNormVar(iSBNorm);
    TH1F* Var = ana->GetVariation(varCounter++);
    protonFeeddown.SetVarHist(Var);
  }
  for (int iLamPro = 0; iLamPro < 3; ++iLamPro) {
    if (iLamPro == 1) {
      continue;  // default case
    }
    ana->Default();
    ana->SetLambdaVar(iLamPro, 1, 1);
    TH1F* Var = ana->GetVariation(varCounter++);
    protonFeeddown.SetVarHist(Var);
  }
  for (int iLamOmg = 0; iLamOmg < 3; ++iLamOmg) {
    if (iLamOmg == 1) {
      continue;  // default case
    }
    ana->Default();
    ana->SetLambdaVar(1, iLamOmg, 1);
    TH1F* Var = ana->GetVariation(varCounter++);
    protonFeeddown.SetVarHist(Var);
  }
  for (int iLamXim1530 = 0; iLamXim1530 < 3; ++iLamXim1530) {
    if (iLamXim1530 == 1) {
      continue;  // default case
    }
    ana->Default();
    ana->SetLambdaVar(1, 1, iLamXim1530);
    TH1F* Var = ana->GetVariation(varCounter++);
    protonFeeddown.SetVarHist(Var);
  }
  for (int iRadXim1530 = 0; iRadXim1530 < 3; ++iRadXim1530) {
    if (iRadXim1530 == 1) {
      continue;  // default case
    }
    ana->Default();
    ana->SetRadXim1530Var(iRadXim1530);
    TH1F* Var = ana->GetVariation(varCounter++);
    protonFeeddown.SetVarHist(Var);
  }
  protonFeeddown.EvalSystematics();
  protonFeeddown.WriteOutput("");

  for (int iRes = 0; iRes < 3; iRes++) {
    if (iRes == 1) {
      continue;  // default case
    }
    ana->Default();
    ana->SetMomentumResolutionVar(iRes);
    TH1F* Var = ana->GetVariation(varCounter++);
    protonMomRes.SetVarHist(Var);
  }
  for (int iIter = 1; iIter < 4; iIter++) {
    ana->SetMomentumResolutionIter(iIter + 4);
    TH1F* Var = ana->GetVariation(varCounter++);
    protonMomRes.SetVarHist(Var);
  }

  protonMomRes.EvalSystematics();
  protonMomRes.WriteOutput("");

  return 0;
}


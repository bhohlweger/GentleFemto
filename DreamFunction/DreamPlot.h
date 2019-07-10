/*
 * DreamPlot.h
 *
 *  Created on: 29 Aug 2018
 *      Author: bernhardhohlweger
 */

#ifndef DREAMFUNCTION_DREAMPLOT_H_
#define DREAMFUNCTION_DREAMPLOT_H_
#include "DreamData.h"
class DreamPlot {
 public:
  DreamPlot();
  virtual ~DreamPlot();
  void ReadData(const char* PathToDataFolder, const char* PathToSysFolder,
                int binWidth, int UnitConvData);
  void ReadDataSigma(const char* PathToDataFolder, const char* PathToSysFolder);
  void ReadSidebandSigma(const char* PathToFitFolder, const char* PathToSysFolder);
  void ReadSimulation(const char* PathToSimFolder, int binWidth);
  void ReadFit(const char* fitPath, int UnitConvCATS);
  void ReadFitSigma(const char* fitPath);
  void SetProtonProtonBaseLine(float ppBL0, float ppBL1) {
    fProtonProton->SetBaseline(ppBL0, ppBL1);
  }
  ;
  void SetProtonLambdaBaseLine(float ppBL0, float ppBL1) {
    fProtonLambda->SetBaseline(ppBL0, ppBL1);
  }
  ;
  void SetLambdaLambdaBaseLine(float ppBL0, float ppBL1) {
    fLambdaLambda->SetBaseline(ppBL0, ppBL1);
  }
  ;
  void SetProtonXiBaseLine(float ppBL0, float ppBL1) {
    fProtonXi->SetBaseline(ppBL0, ppBL1);
  }
  ;
  void SetRadius(float radius, float stat, float sysUp, float sysLow) {
    fRadius = radius;
    fRadiusStat = stat;
    fRadiusSysUp = sysUp;
    fRadiusSysLow = sysLow;
  }
  ;
  void SetCollisionSystem(float energy, const char* system, const char* mcGen) {
    fEnergy = energy;
    fCollisionSystem = system;
    fMonteCarloGen = mcGen;
  }
  static void SetStyle(bool graypalette = false, bool title = false);
  static void SetStyleHisto(TH1 *histo, int marker = 20, int color = kBlue + 2);
  static void SetStyleGraph(TGraph *histo, int marker = 20, int color = kBlue + 2);
  void DrawCorrelationFunctions();
  void DrawCorrelationFunctionSigma(const char* fitPath);
  void DrawSystemInfo(TCanvas* c, bool plotRadius = true, float xMin = 0.35,
                      int isPreliminary = 0);
  DreamData* fProtonProton;
  DreamData* fProtonLambda;
  DreamData* fLambdaLambda;
  DreamData* fProtonXi;
  DreamData* fProtonSigma;
  DreamData* fProtonSigmaSideband;
  float fRadius;
  float fRadiusStat;
  float fRadiusSysUp;
  float fRadiusSysLow;
  float fEnergy;
  const char* fCollisionSystem;
  const char* fMonteCarloGen;
};

#endif /* DREAMFUNCTION_DREAMPLOT_H_ */

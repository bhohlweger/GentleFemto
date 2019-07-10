#include "DreamPlot.h"
#include "global.h"
#include <iostream>
#include <fstream>
#include <string>

int main(int argc, char* argv[]) {

  const char* expBaseDir = argv[1];
  const char* sysBaseDir = argv[2];
  const char* catsBaseDir = argv[3];

  DreamPlot* PlotMe = new DreamPlot();
  PlotMe->SetCollisionSystem(globalBeamEnergy, globalCollisionSystem,
                             globalEventGenerator);
  PlotMe->ReadDataSigma(expBaseDir, sysBaseDir);
  PlotMe->ReadSidebandSigma(catsBaseDir, sysBaseDir);
  PlotMe->ReadFitSigma(catsBaseDir);
  PlotMe->DrawCorrelationFunctionSigma(catsBaseDir);
  return 0;
}

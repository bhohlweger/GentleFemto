#include "DreamPlot.h"
#include "global.h"
#include <iostream>
#include <fstream>
#include <string>

int main(int argc, char* argv[]) {

  const char* expBaseDir = argv[1];
  const char* sysBaseDir = argv[2];

  DreamPlot* PlotMe = new DreamPlot();
  PlotMe->SetCollisionSystem(globalBeamEnergy, globalCollisionSystem,
                             globalEventGenerator);
  PlotMe->ReadDataSigma(expBaseDir, sysBaseDir);
//  //PlotMe->ReadFit(catsFile, 1);
  PlotMe->DrawCorrelationFunctionSigma();
  return 0;
}

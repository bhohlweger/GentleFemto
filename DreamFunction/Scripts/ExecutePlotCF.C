#include "DreamPlot.h"
#include "global.h"

int main(int argc, char* argv[]) {

  const char* expBaseDir = argv[1];
  const char* sysBaseDir = argv[2];
  const char* simBaseDir = argv[3];
  const char* catsFile = argv[4];

  DreamPlot* PlotMe = new DreamPlot();
  PlotMe->SetProtonProtonBaseLine(globalBaselinePP0, globalBaselinePP1);
  PlotMe->SetProtonLambdaBaseLine(globalBaselinePL0, globalBaselinePL1);
  PlotMe->SetLambdaLambdaBaseLine(globalBaselineLL0, globalBaselineLL1);
  PlotMe->SetProtonXiBaseLine(globalBaselinePXi0, globalBaselinePXi1);

  PlotMe->SetRadius(globalRadius, globalRadiusStat, globalRadiusSysUp,
                    globalRadiusSysLow);
  PlotMe->SetCollisionSystem(globalBeamEnergy, globalCollisionSystem,
                             globalEventGenerator);
  PlotMe->ReadData(expBaseDir, sysBaseDir, globalBinWidth, 1000);
  TString simDir = Form("%s", simBaseDir);
  if (simDir != "") {
    PlotMe->ReadSimulation(simBaseDir, globalBinWidth);
  }
  PlotMe->ReadFit(catsFile, 1);
  PlotMe->DrawCorrelationFunctions();
  return 0;
}

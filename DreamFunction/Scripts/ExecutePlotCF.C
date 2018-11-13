#include "DreamPlot.h"
#include "global.h"
#include <iostream>
#include <fstream>
#include <string>

int main(int argc, char* argv[]) {

  const char* expBaseDir = argv[1];
  const char* sysBaseDir = argv[2];
  const char* simBaseDir = argv[3];
  const char* catsFile = argv[4];

  // Read the radius from the file from pp_Systematics
  double r, rStatErr, rSystErrUp, rSystErrLow;
  std::ifstream radius;
  radius.open(Form("%s/radius.dat", catsFile));
  if (!radius.good()) {
    std::cerr << "ERROR: Radius file not found\n";
  } else {
    while (!radius.eof()) {
      radius >> r >> rStatErr >> rSystErrLow >> rSystErrUp;
    }
  }

  DreamPlot* PlotMe = new DreamPlot();
  PlotMe->SetRadius(r, rStatErr, rSystErrUp, rSystErrLow);
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

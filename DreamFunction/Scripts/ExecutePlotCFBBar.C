#include "DreamPlot.h"
#include "global.h"
#include <iostream>
#include <fstream>
#include <string>

int main(int argc, char* argv[]) {

  const char* expBaseDir = argv[1];

  DreamPlot* PlotMe = new DreamPlot();
  PlotMe->SetCollisionSystem(globalBeamEnergy, globalCollisionSystem,
                             globalEventGenerator);
  PlotMe->DrawCorrelationFunctionsBBar(0);
  return 0;
}

#include "DreamPlot.h"
#include "global.h"

int main(int argc, char* argv[]) {
  DreamPlot* PlotPP = new DreamPlot();
  PlotPP->SetCollisionSystem(globalBeamEnergy, globalCollisionSystem,
                             globalEventGenerator);
  PlotPP->DrawCorrelationFunctionProtonProton(argv[1]);
  return 0;
}

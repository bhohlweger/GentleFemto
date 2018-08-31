static int binWidth = 20; //in MeV
static const char* Plotsystem  = "p#minusPb";
static const char* eventGenerator = "DPMJET";
static float energy = 5.02;
static float radius = 1.382; //fm
static float radiusStat = 0.007;
static float radiusSysLow = 0.009;
static float radiusSysUp = 0.014;
static float ppBL0 = 0.985;
static float ppBL1 = 0.183;
static float pLBL0 = 1.021;
static float pLBL1 = 0.000;
static float LLBL0 = 0.971;
static float LLBL1 = 0.041;
static float pXiBL0 = 1.066;
static float pXiBL1 = 0.000;
void PlotCF_pPb(const char* expBaseDir, const char* sysBaseDir,
            const char* simBaseDir, const char* catsFile) {

  DreamPlot* PlotMe = new DreamPlot();
  PlotMe->SetProtonProtonBaseLine(ppBL0,ppBL1);
  PlotMe->SetProtonLambdaBaseLine(pLBL0,pLBL1);
  PlotMe->SetLambdaLambdaBaseLine(LLBL0,LLBL1);
  PlotMe->SetProtonXiBaseLine(pXiBL0,pXiBL1);

  PlotMe->SetRadius(radius,radiusStat,radiusSysUp,radiusSysLow);
  PlotMe->SetCollisionSystem(energy,Plotsystem,eventGenerator);
  PlotMe->ReadData(expBaseDir,sysBaseDir,binWidth);
  TString simDir=Form("%s",simBaseDir);
  if (simDir!="") {
    PlotMe->ReadSimulation(simBaseDir,binWidth);
  }
  PlotMe->ReadFit(catsFile);
  PlotMe->DrawCorrelationFunctions();
  return;
}

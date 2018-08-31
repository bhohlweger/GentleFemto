static int binWidth = 20; //in MeV
static const char* Plotsystem  = "pp (HM)";
static const char* eventGenerator = "Pythia";
static float energy = 13;
static float radius = 1.321; //fm
static float radiusStat = 0.005;
static float radiusSysLow = 0.026;
static float radiusSysUp = 0.013;
static float ppBL0 = 1.002;
static float ppBL1 = 0.177;
static float pLBL0 = 0.980;
static float pLBL1 = 0.131;
static float LLBL0 = 0.943;
static float LLBL1 = 0.099;
static float pXiBL0 = 1.029;
static float pXiBL1 = 0.082;
void PlotCF_ppHM(const char* expBaseDir, const char* sysBaseDir,
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

#include "PlayWithCats.h"

int main(int argc, char *argv[]) {
  PlayWithCats *catsPlay = new PlayWithCats();
  const char* Data = (argv[1]) ? argv[1] : "";
  const char* Fit = (argv[2]) ? argv[2] : "";
  if (Fit != "")catsPlay->ExtractUncertaintyFit(Fit);
  if (Data != "")catsPlay->ExtractUncertaintyData(Data);
  catsPlay->GenerateDefault();
//  catsPlay->ShiftBinning();
  catsPlay->GenerateCoulombOnly();
//  catsPlay->PlotPotentials();
  catsPlay->PlotPotentialSum();

  catsPlay->CloseFile();
  return 0;
}

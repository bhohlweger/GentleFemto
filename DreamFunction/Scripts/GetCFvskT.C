void GetCFvskT(const char* filename, const char* prefix, const char* addon ="") {
  ReadDreamFile* DreamFile = new ReadDreamFile(6, 6);
  DreamKayTee* kTDists;
  DreamFile->ReadkTHistos(filename, prefix, addon);
  kTDists = DreamFile->GetkTPairDistributions(0,0,1,1);
  std::vector<float> kTBins = { 0.48, 0.69, 1., 1.5 };
  //std::vector<float> kTBins = { 0.48, 0.69, 0.9, 1.2 };
  kTDists->SetKayTeeBins(kTBins);
  kTDists->SetNormalization(0.2,0.4);
  kTDists->ObtainTheCorrelationFunction();

  return;
}

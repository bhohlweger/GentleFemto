void GetCFvskT(const char* filename, const char* prefix) {
  ReadDreamFile* DreamFile = new ReadDreamFile(6, 6);
  DreamKayTee* kTDists;
  DreamFile->ReadkTHistos(filename, prefix);
  kTDists = DreamFile->GetkTPairDistributions(0,0,1,1);
  std::vector<float> kTBins = { 0.48, 0.69, 1., 1.5 };
  kTDists->SetKayTeeBins(kTBins);
  kTDists->SetNormalization(0.2,0.4);
  kTDists->ObtainTheCorrelationFunction();

  return;
}

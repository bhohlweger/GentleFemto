#include "ReadDreamFile.h"

int main(int argc, char* argv[]) {
  const char* filename = argv[1];
  const char* prefix = argv[2];
  const char* addon = (argv[3]) ? argv[3] : "";

  ReadDreamFile* DreamFile = new ReadDreamFile(6, 6);
  DreamFile->ReaddEtadPhiHists(4,filename,prefix,addon);
  return 1;
}

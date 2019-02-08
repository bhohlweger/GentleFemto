#include "ForgivingReader.h"
#include "MakeHistosGreat.h"
#include "EventQA.h"

int main(int argc, char* argv[]) {
  const char* filename = argv[1];
  const char* prefix = argv[2];
  const char* addon = (argv[3]) ? argv[3] : "";
  MakeHistosGreat::SetStyle(false,false);
  ForgivingReader* reader = new ForgivingReader(filename,prefix,addon);
  EventQA* evtQA = new EventQA();
  evtQA->SetQAList(reader->GetQA());
  evtQA->SetEventCuts(reader->GetEventCuts());
  evtQA->MakeEventQA();
  return 0;
}


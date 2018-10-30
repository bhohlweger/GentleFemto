#ifndef FORBERNIE_H
#define FORBERNIE_H
#include "TString.h"
void RUN2_main(const unsigned& NumIter, const unsigned& NumJobs,
               const unsigned& JobID, TString InputDir,
               TString OutputDir);
void CALL_BERNIE_AND_VALE(int argc, char *argv[]);


#endif

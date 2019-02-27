/*
 * CandidateCounter.h
 *
 *  Created on: Feb 27, 2019
 *      Author: schmollweger
 */

#ifndef FORGIVINGQA_CANDIDATECOUNTER_H_
#define FORGIVINGQA_CANDIDATECOUNTER_H_
#include "ForgivingReader.h"
class CandidateCounter {
 public:
  CandidateCounter();
  virtual ~CandidateCounter();
  void SetNumberOfCandidates(ForgivingReader* reader);
};

#endif /* FORGIVINGQA_CANDIDATECOUNTER_H_ */

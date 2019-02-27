/*
 * CandidateCounter.cxx
 *
 *  Created on: Feb 27, 2019
 *      Author: schmollweger
 */

#include "CandidateCounter.h"
#include <iostream>
CandidateCounter::CandidateCounter() {
  // TODO Auto-generated constructor stub

}

CandidateCounter::~CandidateCounter() {
  // TODO Auto-generated destructor stub
}

void CandidateCounter::SetNumberOfCandidates(ForgivingReader* reader) {
  auto pTProton = (TH1F*) reader->Get1DHistInList(reader->GetTrackCuts(), {
                                                      "pTDist_after" });
  if (!pTProton) {
    std::cout
        << "CandidateCounter::SetNumberOfCandidates Proton Counter missing \n";
  } else {

  }
  auto pTAntiProton = (TH1F*) reader->Get1DHistInList(
      reader->GetAntiTrackCuts(), { "pTDist_after" });
  if (!pTAntiProton) {
    std::cout
        << "CandidateCounter::SetNumberOfCandidates Anti-Proton Counter missing \n";
  } else {

  }
  auto pTv0 = (TH2F*) reader->Get2DHistInList(
      reader->GetListInList(reader->GetAntiv0Cuts(), { "MinimalBooking" }), {
          "InvMassPt" });
  if (!pTv0) {
    std::cout
        << "CandidateCounter::SetNumberOfCandidates v0 Counter missing \n";
  } else {

  }
  auto pTAntiv0 = (TH2F*) reader->Get2DHistInList(
      reader->GetListInList(reader->GetAntiv0Cuts(), { "MinimalBooking" }), {
          "InvMassPt" });
  if (!pTAntiv0) {
    std::cout
        << "CandidateCounter::SetNumberOfCandidates Anti-v0 Counter missing \n";
  } else {

  }
  auto pTCasc = (TH2F*) reader->Get2DHistInList(
      reader->GetListInList(reader->GetAntiv0Cuts(), { "MinimalBooking" }), {
          "InvMassPt" });
  if (!pTCasc) {
    std::cout
        << "CandidateCounter::SetNumberOfCandidates Cascade Counter missing \n";
  } else {

  }
  auto pTAntiCasc = (TH2F*) reader->Get2DHistInList(
      reader->GetListInList(reader->GetAntiv0Cuts(), { "MinimalBooking" }), {
          "InvMassPt" });
  if (!pTAntiCasc) {
    std::cout
        << "CandidateCounter::SetNumberOfCandidates AntiCascade Counter missing \n";
  } else {

  }

}

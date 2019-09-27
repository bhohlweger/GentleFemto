/*
 * CandidateCounter.cxx
 *
 *  Created on: Feb 27, 2019
 *      Author: schmollweger
 */

#include "CandidateCounter.h"
#include "TDatabasePDG.h"
#include <iostream>
CandidateCounter::CandidateCounter()
    : fnTracks(0),
      fnv0s(0),
      fnCascades(0) ,
      fnAntiTracks(0),
      fnAntiv0s(0)
{
  // TODO Auto-generated constructor stub

}

CandidateCounter::~CandidateCounter() {
  // TODO Auto-generated destructor stub
}

void CandidateCounter::SetNumberOfCandidates(ForgivingReader* reader) {
  float LambdaMass = TDatabasePDG::Instance()->GetParticle(3122)->Mass();
  float XiMass = TDatabasePDG::Instance()->GetParticle(3312)->Mass();
  auto pTProton = (TH1F*) reader->Get1DHistInList(reader->GetTrackCuts(), {"pTDist_after" });
  if (!pTProton) {
    std::cout
        << "CandidateCounter::SetNumberOfCandidates Proton Counter missing \n";
  } else {
//    std::cout << pTProton->GetName() << std::endl;
//    std::cout << "Number of Protons: " << pTProton->GetEntries() << std::endl;
    fnTracks += pTProton->GetEntries();
  }
  auto pTAntiProton = (TH1F*) reader->Get1DHistInList(reader->GetAntiTrackCuts(), { "pTDist_after" });
  if (!pTAntiProton) {
    std::cout
        << "CandidateCounter::SetNumberOfCandidates Anti-Proton Counter missing \n";
  } else {
//    std::cout << pTAntiProton->GetName() << std::endl;
//    std::cout << "Number of AndiProtons: "<< pTAntiProton->GetEntries() << std::endl;
    fnTracks += pTAntiProton->GetEntries();

  }
  auto pTv0 = (TH2F*) reader->Get2DHistInList(reader->GetListInList(reader->Getv0Cuts(), { "MinimalBooking" }), {"InvMassPt" });
  if (!pTv0) {
    std::cout
        << "CandidateCounter::SetNumberOfCandidates v0 Counter missing \n";
  } else {
    auto pTIntv0 = (TH1F*) pTv0->ProjectionY("pTv0Integrated", 1, -1);
    fnv0s += pTIntv0->Integral(pTIntv0->FindBin(LambdaMass - 0.004),
                               pTIntv0->FindBin(LambdaMass + 0.004));

  }
  auto pTAntiv0 = (TH2F*) reader->Get2DHistInList(
      reader->GetListInList(reader->GetAntiv0Cuts(), { "MinimalBooking" }), {
          "InvMassPt" });
  if (!pTAntiv0) {
    std::cout
        << "CandidateCounter::SetNumberOfCandidates Anti-v0 Counter missing \n";
  } else {
    auto pTIntAv0 = (TH1F*) pTAntiv0->ProjectionY("pTAntiv0Integrated", 1, -1);
    fnv0s += pTIntAv0->Integral(pTIntAv0->FindBin(LambdaMass - 0.004),
                                pTIntAv0->FindBin(LambdaMass + 0.004));

  }
  auto pTCasc = (TH2F*) reader->Get2DHistInList(
      reader->GetListInList(reader->GetCascadeCuts(), { "MinimalBooking" }), {
          "InvMassXiPt" });
  if (!pTCasc) {
    std::cout
        << "CandidateCounter::SetNumberOfCandidates Cascade Counter missing \n";
  } else {
    auto pTIntCasc = (TH1F*) pTCasc->ProjectionY("pTCascIntegrated", 1, -1);
    fnCascades += pTIntCasc->Integral(pTIntCasc->FindBin(XiMass - 0.006),
                                      pTIntCasc->FindBin(XiMass + 0.006));
  }
  auto pTAntiCasc = (TH2F*) reader->Get2DHistInList(
      reader->GetListInList(reader->GetAntiCascadeCuts(), { "MinimalBooking" }), {
          "InvMassXiPt" });
  if (!pTAntiCasc) {
    std::cout
        << "CandidateCounter::SetNumberOfCandidates AntiCascade Counter missing \n";
  } else {
    auto pTIntACasc = (TH1F*) pTAntiCasc->ProjectionY("pTAntiCascIntegrated", 1,
                                                      -1);
    fnCascades += pTIntACasc->Integral(pTIntACasc->FindBin(XiMass - 0.006),
                                       pTIntACasc->FindBin(XiMass + 0.006));
  }

}



void CandidateCounter::SetNumberOfCandidatesBBar(ForgivingReader* reader) {
  float LambdaMass = TDatabasePDG::Instance()->GetParticle(3122)->Mass();
  auto pTProton = (TH1F*) reader->Get1DHistInList(reader->GetTrackCuts(), {"pTDist_after" });
  if (!pTProton) {
    std::cout
        << "CandidateCounter::SetNumberOfCandidates Proton Counter missing \n";
  } else {
//    std::cout << pTProton->GetName() << std::endl;
//    std::cout << "Number of Protons: " << pTProton->GetEntries() << std::endl;
    fnTracks += pTProton->GetEntries();
    // printf("--- debug in CandidateCounter 1 ---\n");
  }
  auto pTAntiProton = (TH1F*) reader->Get1DHistInList(reader->GetAntiTrackCuts(), { "pTDist_after" });
  if (!pTAntiProton) {
    std::cout
        << "CandidateCounter::SetNumberOfCandidates Anti-Proton Counter missing \n";
  } else {
//    std::cout << pTAntiProton->GetName() << std::endl;
//    std::cout << "Number of AntiProtons: "<< pTAntiProton->GetEntries() << std::endl;
    fnAntiTracks += pTAntiProton->GetEntries();
    // printf("--- debug in CandidateCounter 2 ---\n");

  }
  auto pTv0 = (TH2F*) reader->Get2DHistInList(reader->GetListInList(reader->Getv0Cuts(), { "MinimalBooking" }), {"InvMassPt" });
  if (!pTv0) {
    std::cout
        << "CandidateCounter::SetNumberOfCandidates v0 Counter missing \n";
  } else {
    auto pTIntv0 = (TH1F*) pTv0->ProjectionY("pTv0Integrated", 1, -1);
    fnv0s += pTIntv0->Integral(pTIntv0->FindBin(LambdaMass - 0.004),
                               pTIntv0->FindBin(LambdaMass + 0.004));
//    std::cout << "Number of Lambdas: " << fnv0s<<std::scientific<<std::endl;

  }
  auto pTAntiv0 = (TH2F*) reader->Get2DHistInList(
      reader->GetListInList(reader->GetAntiv0Cuts(), { "MinimalBooking" }), {
          "InvMassPt" });
  if (!pTAntiv0) {
    std::cout
        << "CandidateCounter::SetNumberOfCandidates Anti-v0 Counter missing \n";
  } else {
    auto pTIntAv0 = (TH1F*) pTAntiv0->ProjectionY("pTAntiv0Integrated", 1, -1);
    fnAntiv0s += pTIntAv0->Integral(pTIntAv0->FindBin(LambdaMass - 0.004),
                                pTIntAv0->FindBin(LambdaMass + 0.004));
//    std::cout << "Number of AntiLambdas: " << fnAntiv0s<<std::scientific<<std::endl;


  }

}

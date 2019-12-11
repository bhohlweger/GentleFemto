/*
 * DreamCF.cxx
 *
 *  Created on: Aug 21, 2018
 *      Author: hohlweger
 */

#include "DreamPair.h"
#include "TCanvas.h"
#include <iostream>
DreamPair::DreamPair(const char* name, float normleft, float normright)
    : fPair(nullptr),
      fPairShifted(),
      fPairRebinned(),
      fPairReweighted(),
      fPairUnfolded(),
      fFirstBin(-99),
      fNormLeft(normleft),
      fNormRight(normright),
      fName(name) {
}

DreamPair::~DreamPair() {
  if (fPair) {
    delete fPair;
  }
  for (auto it : fPairShifted) {
    delete it;
  }
  for (auto it : fPairFixShifted) {
    delete it;
  }
  for (auto it : fPairRebinned) {
    delete it;
  }
  for (auto it : fPairReweighted) {
    delete it;
  }
  for (auto it : fPairUnfolded) {
    delete it;
  }
}

int DreamPair::GetNDists() {
  //1 for the fPair itself
  return 1 + fPairShifted.size() + fPairRebinned.size() + fPairReweighted.size()
      + fPairUnfolded.size();
}

void DreamPair::ShiftForEmpty(DreamDist* pair) {
  DreamDist* PairShifted = new DreamDist();

  TH1F* SE = pair->GetSEDist();
  TH1F* ME = pair->GetMEDist();

  TH2F* SEMult = pair->GetSEMultDist();
  TH2F* MEMult = pair->GetMEMultDist();
  int FirstBin = SE->FindFirstBinAbove(0);
  //+1 since we start counting bins at 1
  const int nBinsEffSE = SE->GetNbinsX() - FirstBin + 1;
  fFirstBin = SE->GetXaxis()->GetBinLowEdge(FirstBin);
  const float SEkMax = SE->GetXaxis()->GetBinUpEdge(SE->GetNbinsX());
  const int multBins = SEMult->GetNbinsY();
  const int multMax = SEMult->GetYaxis()->GetBinUpEdge(multBins);

  const char* SEHistName = Form("%s_", SE->GetName());
  TH1F* SEShifted = new TH1F(SEHistName, SEHistName, nBinsEffSE, fFirstBin,
                             SEkMax);
  SEShifted->Sumw2();
  const char* MEHistName = Form("%s_", ME->GetName());
  TH1F* MEShifted = new TH1F(MEHistName, MEHistName, nBinsEffSE, fFirstBin,
                             SEkMax);
  MEShifted->Sumw2();
  const char* SEMultHistName = Form("%s_", SEMult->GetName());
  TH2F* SEMultShifted = new TH2F(SEMultHistName, SEMultHistName, nBinsEffSE,
                                 fFirstBin, SEkMax, multBins, 1, multMax);
  SEMultShifted->Sumw2();
  const char* MEMultHistName = Form("%s_", MEMult->GetName());
  TH2F* MEMultShifted = new TH2F(MEMultHistName, MEMultHistName, nBinsEffSE,
                                 fFirstBin, SEkMax, multBins, 1, multMax);
  MEMultShifted->Sumw2();
  int ckBin = 1;
  for (int ikBin = FirstBin; ikBin <= SE->GetNbinsX(); ++ikBin) {
    SEShifted->SetBinContent(ckBin, SE->GetBinContent(ikBin));
    SEShifted->SetBinError(ckBin, SE->GetBinError(ikBin));

    MEShifted->SetBinContent(ckBin, ME->GetBinContent(ikBin));
    MEShifted->SetBinError(ckBin, ME->GetBinError(ikBin));
    for (int iMult = 1; iMult <= multBins; ++iMult) {
      SEMultShifted->SetBinContent(ckBin, iMult,
                                   SEMult->GetBinContent(ikBin, iMult));
      SEMultShifted->SetBinError(ckBin, iMult,
                                 SEMult->GetBinError(ikBin, iMult));

      MEMultShifted->SetBinContent(ckBin, iMult,
                                   MEMult->GetBinContent(ikBin, iMult));
      MEMultShifted->SetBinError(ckBin, iMult,
                                 MEMult->GetBinError(ikBin, iMult));
    }
    ckBin++;
  }
  PairShifted->SetSEDist(SEShifted, "Shifted");
  PairShifted->SetSEMultDist(SEMultShifted, "Shifted");
  PairShifted->SetMEDist(MEShifted, "Shifted");
  PairShifted->SetMEMultDist(MEMultShifted, "Shifted");
  PairShifted->Calculate_CF(fNormLeft, fNormRight);
  fPairShifted.push_back(PairShifted);
  delete SEShifted;
  delete SEMultShifted;
  delete MEShifted;
  delete MEMultShifted;
  //  TH1F* SEProjMult[multBins];
  //  TH1F* SEProjMultShifted[multBins];
  //  TH1F* MEProjMult[multBins];
  //  TH1F* MEProjMultShifted[multBins];
  //  const char* can1Name=Form("c1%s",SE->GetName());
  //  TCanvas *c1 = new TCanvas(can1Name,can1Name,2000,1000);
  //  c1->Divide(2,2);
  //  c1->cd(1);
  //  SEMult->Draw();
  //  c1->cd(2);
  //  SEMultShifted->Draw();
  //  c1->cd(3);
  //  MEMult->Draw();
  //  c1->cd(4);
  //  MEMultShifted->Draw();
  //
  //  const char* can2Name=Form("c2%s",SE->GetName());
  //  TCanvas *c2 = new TCanvas(can2Name,can2Name,2000,1000);
  //  c2->Divide(7,4);

  //  for (int iMult=1;iMult<=multBins;++iMult)
  //  {
  //    TString projSEMultName=Form("%sProj%i",SEMult->GetName(),iMult);
  //    SEProjMult[iMult-1]=(TH1F*)SEMult->ProjectionX(projSEMultName.Data(),iMult,iMult);
  ////    TString projSEMultShiftedName=Form("%sProj%i",SEMultShifted->GetName(),iMult);
  ////    SEProjMultShifted[iMult-1]=(TH1F*)SEMultShifted->ProjectionX("",iMult,iMult);
  //    SEProjMult[iMult-1]->Divide((TH1F*)SEMultShifted->ProjectionX("",iMult,iMult));
  //    c1->cd(iMult);
  //    SEProjMult[iMult-1]->Draw("hist");
  //
  //    TString projMEMultName=Form("%sProj%i",MEMult->GetName(),iMult);
  //    MEProjMult[iMult-1]=(TH1F*)MEMult->ProjectionX(projMEMultName.Data(),iMult,iMult);
  ////    TString projMEMultShiftedName=Form("%sProj%i",MEMultShifted->GetName(),iMult);
  ////    MEProjMultShifted[iMult-1]=(TH1F*)MEMultShifted->ProjectionX("",iMult,iMult);
  //    MEProjMult[iMult-1]->Divide((TH1F*)MEMultShifted->ProjectionX("",iMult,iMult));
  //    c2->cd(iMult);
  //    MEProjMult[iMult-1]->Draw("hist");
  //  }
  return;
}

void DreamPair::FixShift(DreamDist* pair, DreamDist* otherDist, float kMin,
                         const bool fixedShift) {
  //in case of a manual fix shift, the user NEEDS to make sure
  //that kMin is on the boundaries of one bin
  if (!fixedShift && ((fFirstBin == -99) || kMin == -99)) {
    std::cout << "Internal kStar=" << fFirstBin << " and external one kStar="
              << kMin << ". Check if you ran ShiftForEmpty!" << std::endl;
  } else {
    //we only need to do this in case this pair has a different starting bin.
    if (fFirstBin > kMin || fixedShift) {
      TH1F* otherSE = fixedShift ? nullptr : otherDist->GetSEDist();

      TH1F* SE = pair->GetSEDist();
      TH2F* SEMult = pair->GetSEMultDist();
      TH1F* ME = pair->GetMEDist();
      TH2F* MEMult = pair->GetMEMultDist();
      //The epsilon ensures that always the lower bin is picked!
      double epsilon = SE->GetBinWidth(1) * 1e-2;
      int nBins =
          fixedShift ?
              SE->GetXaxis()->GetNbins() - SE->FindBin(kMin - epsilon) :
              otherSE->GetXaxis()->GetNbins();
      float xMin = fixedShift ? kMin : otherSE->GetXaxis()->GetXmin();
      float xMax =
          fixedShift ?
              SE->GetXaxis()->GetXmax() : otherSE->GetXaxis()->GetXmax();

      int multBins = 0;
      int multMax = 0;
      if (SEMult && MEMult) {
        multBins = SEMult->GetYaxis()->GetNbins();
        multMax = SEMult->GetYaxis()->GetXmax();
      }
      TH1F* SEShifted;
      TH1F* MEShifted;
      TH2F* SEMultShifted = nullptr;
      TH2F* MEMultShifted = nullptr;
      const char* SEHistName = Form("%s_", SE->GetName());
      std::cout << "Fix shifed nBins: " << nBins << " xMin: " << xMin
                << " xMax: " << xMax << std::endl;
      SEShifted = new TH1F(SEHistName, SEHistName, nBins, xMin, xMax);
      SEShifted->Sumw2();
      const char* MEHistName = Form("%s_", ME->GetName());
      MEShifted = new TH1F(MEHistName, MEHistName, nBins, xMin, xMax);
      MEShifted->Sumw2();
      if (SEMult && MEMult) {
        const char* SEMultHistName = Form("%s_", SEMult->GetName());
        SEMultShifted = new TH2F(SEMultHistName, SEMultHistName, nBins, xMin,
                                 xMax, multBins, 1, multMax);
        SEMultShifted->Sumw2();
        const char* MEMultHistName = Form("%s_", MEMult->GetName());
        MEMultShifted = new TH2F(MEMultHistName, MEMultHistName, nBins, xMin,
                                 xMax, multBins, 1, multMax);
        MEMultShifted->Sumw2();
      }
      int startBin = fixedShift ? 1 : SEShifted->FindBin(fFirstBin);
      int endBin = SEShifted->GetNbinsX();
      //for the fixedShift == true case the first bin has to be found,
      //while if the other histogram is already shifted it is just the first bin.
      int iOtherBin = fixedShift ? SE->FindBin(kMin) : 1;
      int endOtherBin = SE->GetXaxis()->GetNbins();
      for (int iBin = startBin; iBin <= endBin; ++iBin) {
        if (iOtherBin < endOtherBin) {
          SEShifted->SetBinContent(iBin, SE->GetBinContent(iOtherBin));
          SEShifted->SetBinError(iBin, SE->GetBinError(iOtherBin));
          MEShifted->SetBinContent(iBin, ME->GetBinContent(iOtherBin));
          MEShifted->SetBinError(iBin, ME->GetBinError(iOtherBin));
          if (SEMult && MEMult) {
            for (int iMult = 1; iMult <= multBins; ++iMult) {
              SEMultShifted->SetBinContent(
                  iBin, iMult, SEMult->GetBinContent(iOtherBin, iMult));
              SEMultShifted->SetBinError(iBin, iMult,
                                         SEMult->GetBinError(iOtherBin, iMult));

              MEMultShifted->SetBinContent(
                  iBin, iMult, MEMult->GetBinContent(iOtherBin, iMult));
              MEMultShifted->SetBinError(iBin, iMult,
                                         MEMult->GetBinError(iOtherBin, iMult));
            }
          }
          iOtherBin++;
        } else {
          continue;
        }
      }
      DreamDist* PairFixShifted = new DreamDist();
      PairFixShifted->SetSEDist(SEShifted, "FixShifted");
      PairFixShifted->SetMEDist(MEShifted, "FixShifted");
      if (SEMultShifted)PairFixShifted->SetSEMultDist(SEMultShifted, "FixShifted");
      if (MEMultShifted)PairFixShifted->SetMEMultDist(MEMultShifted, "FixShifted");
      PairFixShifted->Calculate_CF(fNormLeft, fNormRight);
      fPairFixShifted.push_back(PairFixShifted);
      delete SEShifted;
      delete MEShifted;
      if (SEMultShifted)delete SEMultShifted;
      if (MEMultShifted)delete MEMultShifted;
    } else {
      DreamDist* PairFixShifted = new DreamDist(pair, "_FixShifted");
      PairFixShifted->Calculate_CF(fNormLeft, fNormRight);
      fPairFixShifted.push_back(PairFixShifted);
    }
  }
  return;
}

void DreamPair::FixShift(DreamDist* pair, DreamDist* otherPair1,
                         DreamDist* otherPair2, float kMin1, float kMin2) {
  const float kMin = (kMin1 < kMin2) ? kMin1 : kMin2;
  DreamDist *otherDist = (kMin1 < kMin2) ? otherPair1 : otherPair2;
  FixShift(pair, otherDist, kMin, true);
}

void DreamPair::Rebin(DreamDist* pair, int rebin, bool seMean) {
  DreamDist* Rebinned = new DreamDist(pair, Form("_Rebinned_%i", rebin));
  std::cout << "Before Rebinned->GetSEDist()->GetBinWidth(1): "
            << Rebinned->GetSEDist()->GetBinWidth(1) << std::endl;
  Rebinned->GetSEDist()->Rebin(rebin);
  std::cout << "After rebinning by: " << rebin
            << " Rebinned->GetSEDist()->GetBinWidth(1): "
            << Rebinned->GetSEDist()->GetBinWidth(1) << std::endl;
  if (Rebinned->GetSEMultDist())
    Rebinned->GetSEMultDist()->Rebin2D(rebin, 1);
  Rebinned->GetMEDist()->Rebin(rebin);
  if (Rebinned->GetMEMultDist())
    Rebinned->GetMEMultDist()->Rebin2D(rebin, 1);
  if (seMean) {
    Rebinned->Calculate_CF(fNormLeft, fNormRight, pair->GetSEDist());
  } else {
    Rebinned->Calculate_CF(fNormLeft, fNormRight);
  }
  fPairRebinned.push_back(Rebinned);
  return;
}

void DreamPair::ReweightMixedEvent(DreamDist* pair, float kSMin, float kSMax,
                                   DreamDist* pairNotRebinned) {
  DreamDist* PairReweighted = new DreamDist(pair, "_Reweighted");

  TH1F* SE = pair->GetSEDist();
  TH1F* ME = pair->GetMEDist();
  TH2F* SEMult = pair->GetSEMultDist();
  TH2F* MEMult = pair->GetMEMultDist();

  TH1F* MEReweighted = PairReweighted->GetMEDist();
  TH2F* MEMultReweighted = PairReweighted->GetMEMultDist();
  MEReweighted->Reset();  // this one we fill from scratch, we only want to keep the dimensions
  MEMultReweighted->Reset();

  //k* range of the normalization of ME and SE in Multiplicity
  int firstBin = SEMult->GetXaxis()->FindBin(kSMin);
  int lastBin = SEMult->GetXaxis()->FindBin(kSMax);
  TH1D* MultProjSE = SEMult->ProjectionY(Form("%skSProj", SEMult->GetName()),
                                         firstBin, lastBin);
  TH1D* MultProjME = MEMult->ProjectionY(Form("%skSProj", MEMult->GetName()),
                                         firstBin, lastBin);
  TString MultProjRewName = Form("Weighted%s", MultProjME->GetName());
  //  TH1D* MultProjMEReweighted = (TH1D*)MultProjME->Clone(MultProjRewName.Data());
  //  MultProjMEReweighted->Reset();

  int nMultBins = MEMult->GetYaxis()->GetNbins();
  int multMax = MEMult->GetYaxis()->GetBinUpEdge(nMultBins);
  double weight = 0;
  for (int iMult = 1; iMult <= nMultBins; ++iMult) {
    if (MultProjME->GetBinContent(iMult) > 0) {
      weight = MultProjSE->GetBinContent(iMult)
          / MultProjME->GetBinContent(iMult);
    } else {
      weight = 0;
//      std::cout << "Weight = 0 for Mult Bin iMult = " << iMult
//                << " in the case of " << SE->GetName() << std::endl;
    }
    //    MultProjMEReweighted->SetBinContent(iMult,weight*MultProjME->GetBinContent(iMult));
    TString MultBinName = Form("MEBin%i", iMult);
    MEReweighted->Add(MEMult->ProjectionX(MultBinName.Data(), iMult, iMult),
                      weight);

    for (int ikStar = 1; ikStar <= MEMult->GetNbinsX(); ++ikStar) {
      MEMultReweighted->SetBinContent(
          ikStar, iMult, MEMult->GetBinContent(ikStar, iMult) * weight);
      MEMultReweighted->SetBinError(
          ikStar, iMult, MEMult->GetBinError(ikStar, iMult) * weight);
    }
  }
  //  TString CanName = Form("Can%s",SE->GetName());
  //  TCanvas* c1 = new TCanvas(CanName.Data(),CanName.Data(),2000,1000);
  //  c1->cd();
  //  MultProjSE->SetLineColor(2);
  //  MultProjSE->Scale(1./MultProjSE->Integral());
  //  MultProjSE->DrawCopy();
  //  MultProjME->SetLineColor(3);
  //  MultProjME->Scale(1./MultProjME->Integral());
  //  MultProjME->DrawCopy("same");
  //  MultProjMEReweighted->SetLineColor(4);
  //  MultProjMEReweighted->SetLineStyle(4);
  //  MultProjMEReweighted->Scale(1./MultProjMEReweighted->Integral());
  //  MultProjMEReweighted->DrawCopy("same");
  if (pairNotRebinned) {
    PairReweighted->Calculate_CF(fNormLeft, fNormRight,
                                 pairNotRebinned->GetSEDist());
  } else {
    PairReweighted->Calculate_CF(fNormLeft, fNormRight);
  }
  fPairReweighted.push_back(PairReweighted);
  delete MultProjSE;
  delete MultProjME;

  //  const char* can1Name=Form("c1%s",SE->GetName());
  //  TCanvas *c1 = new TCanvas(can1Name,can1Name,2000,1000);
  //  c1->Divide(2,2);
  //  c1->cd(1);
  //  PairReweighted->GetSEMultDist()->Draw("COLZ");
  //  c1->cd(2);
  //  SEMult->Draw("COLZ");
  //  c1->cd(3);
  //  MEMultReweighted->Draw("COLZ");
  //  c1->cd(4);
  //  MEMult->Draw("COLZ");
  //
  //  const char* can2Name=Form("c2%s",SE->GetName());
  //  TCanvas *c2 = new TCanvas(can2Name,can2Name,2000,1000);
  //  c2->Divide(2,2);
  //  c2->cd(1);
  //  SE->Draw();
  //  c2->cd(2);
  //  PairReweighted->GetSEDist()->Draw("same");
  //
  //  c2->cd(3);
  //  ME->Scale(1./ME->Integral());
  //  ME->Draw();
  //  c2->cd(4);
  //
  //  MEReweighted->SetLineColor(2);
  //  MEReweighted->Scale(1./MEReweighted->Integral());
  //  MEReweighted->Draw("same");
  return;
}

void DreamPair::UnfoldMomentum(DreamDist* pair, MomentumGami *mom) {
  if (!pair) {
    std::cout << "No pair set\n";
    return;
  }
  DreamDist* Unfolded = new DreamDist();
  Unfolded->SetSEDist(mom->UnfoldviaRooResp(pair->GetSEDist()),"");
  Unfolded->SetMEDist(mom->UnfoldviaRooResp(pair->GetMEDist()),"");
  Unfolded->Calculate_CF(fNormLeft, fNormRight);
  fPairUnfolded.push_back(Unfolded);

  return;
}

void DreamPair::WriteOutput(TList *Outlist) {
  TList *PairList = new TList();
  PairList->SetName("Pair");
  PairList->SetOwner();
  Outlist->Add(PairList);

  TList *PairShiftedList = new TList();
  PairShiftedList->SetName("PairShifted");
  PairShiftedList->SetOwner();
  Outlist->Add(PairShiftedList);

  TList *PairFixShiftedList = new TList();
  PairFixShiftedList->SetName("PairFixShifted");
  PairFixShiftedList->SetOwner();
  Outlist->Add(PairFixShiftedList);

  TList *PairRebinnedList = new TList();
  PairRebinnedList->SetName("PairRebinned");
  PairRebinnedList->SetOwner();
  Outlist->Add(PairRebinnedList);

  TList *PairReweightedList = new TList();
  PairReweightedList->SetName("PairReweighted");
  PairReweightedList->SetOwner();
  Outlist->Add(PairReweightedList);

  TList *PairUnfoldedList = new TList();
  PairUnfoldedList->SetName("PairUnfolded");
  PairUnfoldedList->SetOwner();
  Outlist->Add(PairUnfoldedList);

  fPair->WriteOutput(PairList);
  for (auto& it : fPairShifted)
    it->WriteOutput(PairShiftedList);
  for (auto& it : fPairFixShifted)
    it->WriteOutput(PairFixShiftedList);
  for (auto& it : fPairRebinned)
    it->WriteOutput(PairRebinnedList);
  for (auto& it : fPairReweighted)
    it->WriteOutput(PairReweightedList);
  for (auto& it : fPairUnfolded)
    it->WriteOutput(PairUnfoldedList);
}

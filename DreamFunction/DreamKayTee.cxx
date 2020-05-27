/*
 * DreamKayTee.cxx
 *
 *  Created on: Sep 10, 2018
 *      Author: hohlweger
 */

#include "DreamKayTee.h"
#include "DreamPlot.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH2F.h"
#include "TSystem.h"
#include <iostream>
#include <algorithm> 
DreamKayTee::DreamKayTee(const int imTMult)
  : fIskT(true),
    nmTBins(imTMult), 
    fFixShift(),
    fFixShiftValue(),
    fKayTeeBins(),
    fRebin(),
    fNKayTeeBins(0),
    fSEmTMult(nullptr),
    fMEmTMult(nullptr),
    fAveragekT(0),
    fCFPart(nullptr),
    fSum(nullptr),
    fNormleft(0),
    fNormright(0),
    fSEMEReweighting(nullptr),
    fSEMEReweightingMeV(nullptr) {
  fSEkT[0] = nullptr;
  fSEkT[1] = nullptr;
  fMEkT[0] = nullptr;
  fMEkT[1] = nullptr;
  fSEmTMult = new TH2F**[imTMult]; 
  fMEmTMult = new TH2F**[imTMult];
  for (int imT = 0; imT < imTMult; ++imT) {
    fSEmTMult[imT] = new TH2F*[2]; //particle and antiparticle pair
    fMEmTMult[imT] = new TH2F*[2]; //particle and antiparticle pair
  }
}

DreamKayTee::~DreamKayTee() {
  // TODO Auto-generated destructor stub
}

void DreamKayTee::ObtainTheCorrelationFunction(const char* outFolder,
                                               const char* prefix,
                                               const char* pair) {
  if (!fSEMEReweighting) {
    std::cout << "No Reweighting set. Call SetSEMEReweightingRatio first! \n";
  } else {
    const char* variable = (fIskT) ? "kT" : "mT";
    const int nBins = (int) fKayTeeBins.size();
    if (fSEkT[0]) {
      fKayTeeBins.push_back(fSEkT[0]->GetYaxis()->GetXmax());
      fCFPart = new DreamPair**[2];
      TString PartName[2] = { "Part", "AntiPart" };
      for (int iPart = 0; iPart < 2; ++iPart) {
        fCFPart[iPart] = new DreamPair*[fNKayTeeBins];
        for (int ikT = 0; ikT < fNKayTeeBins - 1; ++ikT) {
          TString PairName = Form("%s_%s_%i", PartName[iPart].Data(), variable,
                                  ikT);
          fCFPart[iPart][ikT] = new DreamPair(PairName.Data(), fNormleft,
                                              fNormright);
          //multiplications due in order not to hit the border of the bin!
          int kTminBin = fSEkT[iPart]->GetYaxis()->FindBin(
              fKayTeeBins[ikT] * 1.0001);
          int kTmaxBin = fSEkT[iPart]->GetYaxis()->FindBin(
              fKayTeeBins[ikT + 1] * 0.9999);

          DreamDist* kTDist = new DreamDist();
          TString SEkTBinName = Form("SE%s", PairName.Data());
          kTDist->SetSEDist(
              (TH1F*) fSEkT[iPart]->ProjectionX(SEkTBinName.Data(), kTminBin,
                                                kTmaxBin, "e"),
              "");
          TString MEkTBinName = Form("ME%s", PairName.Data());
          kTDist->SetMEDist(
              (TH1F*) fMEkT[iPart]->ProjectionX(MEkTBinName.Data(), kTminBin,
                                                kTmaxBin, "e"),
              "");
          //little hack to make things work.
          TH2F* fake = new TH2F("fakeMult", "fakeMult", 10, 0, 1, 10, 0, 1);
          kTDist->SetSEMultDist(fake, "SE");
          kTDist->SetMEMultDist(fake, "ME");
          fCFPart[iPart][ikT]->SetPair(kTDist);
          if (fFixShift.size() > 0 && fFixShift.at(ikT)) {
            fCFPart[iPart][ikT]->ShiftForEmpty(fCFPart[iPart][ikT]->GetPair());
            fCFPart[iPart][ikT]->FixShift(fCFPart[iPart][ikT]->GetPair(),
                                          nullptr, fFixShiftValue.at(ikT),
                                          true);
            for (auto it : fRebin) {
              fCFPart[iPart][ikT]->Rebin(
                  fCFPart[iPart][ikT]->GetFixShifted().at(0), it);
            }
          } else {
            for (auto it : fRebin) {
              fCFPart[iPart][ikT]->Rebin(fCFPart[iPart][ikT]->GetPair(), it);
            }
          }
        }
      }
      this->AveragekT(pair);
      fSum = new DreamCF*[fNKayTeeBins];
      TString outname = TString::Format("%s/AveragemT.root",
                                        outFolder);
      TFile* allCFsOut = TFile::Open(outname.Data(), "UPDATE");
      allCFsOut->cd();
      if (fAveragekT) {
        TString name = TString::Format("Average%s_%s", variable,pair);
        TGraphErrors* copy = new TGraphErrors(*fAveragekT);
        copy->GetXaxis()->SetTitle("m_{T} Bin");
        copy->GetYaxis()->SetTitle("<m_{T}> (GeV/c^{2})");
        copy->SetName(name.Data());
        copy->Write(name.Data());
        allCFsOut->Write();
      }
      allCFsOut->Close();
      for (int ikT = 0; ikT < fNKayTeeBins - 1; ++ikT) {
        fSum[ikT] = new DreamCF();
        fSum[ikT]->SetPairs(fCFPart[0][ikT], fCFPart[1][ikT]);
        fSum[ikT]->GetCorrelations(pair);
        std::vector<TH1F*> CFs = fSum[ikT]->GetCorrelationFunctions();
        TString outfileName = TString::Format("%s/CFOutput_%s_%s_%s_%u.root",
                                              outFolder, variable, pair, prefix,
                                              ikT);
        for (auto &it : CFs) {
          TString HistName = it->GetName();
          for (int iBin = 1; iBin < it->GetNbinsX(); ++iBin) {
            float binCont = it->GetBinContent(iBin);
            float binErr = it->GetBinError(iBin);
            float binCent = it->GetBinCenter(iBin);
            float corrFactor = 0;
            if (HistName.Contains("MeV")) {
              corrFactor = fSEMEReweightingMeV->GetBinContent(
                  fSEMEReweightingMeV->FindBin(binCent));
            } else {
              corrFactor = fSEMEReweighting->GetBinContent(
                  fSEMEReweighting->FindBin(binCent));
            }
            if (corrFactor != 0) {
              it->SetBinContent(iBin, binCont * corrFactor);
              it->SetBinError(iBin, binErr);
            }
//            else {
//              //dont worry if iBin seems to change, its due to the rebinning!
//              std::cout << "======================================\n";
//              std::cout << it->GetName() << std::endl;
//              std::cout << "!for ikT = " << ikT << " and iBin = " << iBin
//                        << "!\n";
//              std::cout << "Correction factor 0, CFs meaningless!!! \n";
//              std::cout << "======================================\n";
//              std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
//            }
          }
        }
        fSum[ikT]->WriteOutput(outfileName.Data());
      }
    } else {
      std::cout << "No SE " << variable << " histogram with index 0 \n";
    }
  }
  return;
}

void DreamKayTee::ObtainTheCorrelationFunctionBBar(const char* outFolder,
                                               const char* prefix,
                                               const char* pair) {
  if (!fSEMEReweighting) {
    std::cout << "No Reweighting set. Call SetSEMEReweightingRatio first! \n";
  } else {
    const char* variable = (fIskT) ? "kT" : "mT";
    const int nBins = (int) fKayTeeBins.size();
    if (fSEkT[0]) {
      fKayTeeBins.push_back(fSEkT[0]->GetYaxis()->GetXmax());
      fCFPart = new DreamPair**[1];
      TString PartName[2] = { "Part"};
      for (int iPart = 0; iPart < 1; ++iPart) {
        fCFPart[iPart] = new DreamPair*[fNKayTeeBins];
        for (int ikT = 0; ikT < fNKayTeeBins - 1; ++ikT) {
          TString PairName = Form("%s_%s_%i", PartName[iPart].Data(), variable,
                                  ikT);
          fCFPart[iPart][ikT] = new DreamPair(PairName.Data(), fNormleft,
                                              fNormright);
          //multiplications due in order not to hit the border of the bin!
          int kTminBin = fSEkT[iPart]->GetYaxis()->FindBin(
              fKayTeeBins[ikT] * 1.0001);
          int kTmaxBin = fSEkT[iPart]->GetYaxis()->FindBin(
              fKayTeeBins[ikT + 1] * 0.9999);

          DreamDist* kTDist = new DreamDist();
          TString SEkTBinName = Form("SE%s", PairName.Data());
          kTDist->SetSEDist(
              (TH1F*) fSEkT[iPart]->ProjectionX(SEkTBinName.Data(), kTminBin,
                                                kTmaxBin, "e"),
              "");
          TString MEkTBinName = Form("ME%s", PairName.Data());
          kTDist->SetMEDist(
              (TH1F*) fMEkT[iPart]->ProjectionX(MEkTBinName.Data(), kTminBin,
                                                kTmaxBin, "e"),
              "");
          //little hack to make things work.
          TH2F* fake = new TH2F("fakeMult", "fakeMult", 10, 0, 1, 10, 0, 1);
          kTDist->SetSEMultDist(fake, "SE");
          kTDist->SetMEMultDist(fake, "ME");
          fCFPart[iPart][ikT]->SetPair(kTDist);
          if (fFixShift.size() > 0 && fFixShift.at(ikT)) {
            fCFPart[iPart][ikT]->ShiftForEmpty(fCFPart[iPart][ikT]->GetPair());
            fCFPart[iPart][ikT]->FixShift(fCFPart[iPart][ikT]->GetPair(),
                                          nullptr, fFixShiftValue.at(ikT),
                                          true);
            for (auto it : fRebin) {
              fCFPart[iPart][ikT]->Rebin(
                  fCFPart[iPart][ikT]->GetFixShifted().at(0), it);
            }
          } else {
            for (auto it : fRebin) {
              fCFPart[iPart][ikT]->Rebin(fCFPart[iPart][ikT]->GetPair(), it);
            }
          }
        }
      }
      this->AveragekT(pair);
      fSum = new DreamCF*[fNKayTeeBins];
      TString outname = TString::Format("%s/AveragemT.root",
                                        outFolder);
      TFile* allCFsOut = TFile::Open(outname.Data(), "UPDATE");
      allCFsOut->cd();
      if (fAveragekT) {
        TString name = TString::Format("Average%s_%s", variable,pair);
        TGraphErrors* copy = new TGraphErrors(*fAveragekT);
        copy->GetXaxis()->SetTitle("m_{T} Bin");
        copy->GetYaxis()->SetTitle("<m_{T}> (GeV/c^{2})");
        copy->SetName(name.Data());
        copy->Write(name.Data());
        allCFsOut->Write();
      }
      allCFsOut->Close();
      for (int ikT = 0; ikT < fNKayTeeBins - 1; ++ikT) {
        fSum[ikT] = new DreamCF();
        fSum[ikT]->SetPairsBBar(fCFPart[0][ikT]);
        fSum[ikT]->GetCorrelations(pair);
        std::vector<TH1F*> CFs = fSum[ikT]->GetCorrelationFunctions();
        TString outfileName = TString::Format("%s/CFOutput_%s_%s_%s_%u.root",
                                              outFolder, variable, pair, prefix,
                                              ikT);
        for (auto &it : CFs) {
          TString HistName = it->GetName();
          for (int iBin = 1; iBin < it->GetNbinsX(); ++iBin) {
            float binCont = it->GetBinContent(iBin);
            float binErr = it->GetBinError(iBin);
            float binCent = it->GetBinCenter(iBin);
            float corrFactor = 0;
            if (HistName.Contains("MeV")) {
              corrFactor = fSEMEReweightingMeV->GetBinContent(
                  fSEMEReweightingMeV->FindBin(binCent));
            } else {
              corrFactor = fSEMEReweighting->GetBinContent(
                  fSEMEReweighting->FindBin(binCent));
            }
            if (corrFactor != 0) {
              it->SetBinContent(iBin, binCont * corrFactor);
              it->SetBinError(iBin, binErr);
            }
//            else {
//              //dont worry if iBin seems to change, its due to the rebinning!
//              std::cout << "======================================\n";
//              std::cout << it->GetName() << std::endl;
//              std::cout << "!for ikT = " << ikT << " and iBin = " << iBin
//                        << "!\n";
//              std::cout << "Correction factor 0, CFs meaningless!!! \n";
//              std::cout << "======================================\n";
//              std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
//            }
          }
        }
        fSum[ikT]->WriteOutput(outfileName.Data());
      }
    } else {
      std::cout << "No SE " << variable << " histogram with index 0 \n";
    }
  }
  return;
}

void DreamKayTee::AveragekT(const char *pair) {
  DreamPlot::SetStyle();
  const char* variable = (fIskT) ? "kT" : "mT";
  if (fAveragekT) {
    delete fAveragekT;
  }
  fAveragekT = new TGraphErrors();
  TH2F* kTkStar = (TH2F*) fSEkT[0]->Clone(
      Form("%s%skStarForAverage", variable, pair));
  kTkStar->Add(fSEkT[1]);
  TH1F* kTProjection = (TH1F*) kTkStar->ProjectionY(
      Form("%s%sDist", variable, pair), kTkStar->GetXaxis()->FindBin(0.),
      kTkStar->GetXaxis()->FindBin(0.2), "e");
  auto *c1 = new TCanvas(Form("c1%s", pair), Form("c1%s", pair));
  auto *c2 = new TCanvas(Form("c2%s", pair), Form("c2%s", pair));
  c1->cd();
  kTkStar->Draw("COLZ");
  c1->SaveAs(Form("mTvskStar%s.pdf", pair));

  float totalPairs = kTProjection->Integral();
  DreamPlot::SetStyleHisto(kTProjection);
  kTProjection->GetXaxis()->SetTitle("<m_{T}> (GeV/c^{2})");
  c2->cd();
  kTProjection->DrawClone();
  TString outCanName = TString::Format("mTDistribution%s.pdf", pair);
  c2->SaveAs(outCanName.Data());
  for (int ikT = 0; ikT < fNKayTeeBins - 1; ++ikT) {
    int binLow = kTProjection->GetXaxis()->FindBin(
        fKayTeeBins.at(ikT) * 1.0001);
    int binUp = kTProjection->GetXaxis()->FindBin(
        fKayTeeBins.at(ikT + 1) * 0.9999);
    kTProjection->GetXaxis()->SetRange(binLow, binUp);
    float binErr = kTProjection->GetRMS();
    if (binLow == binUp) {
      std::cout << "Low Edge: " << kTProjection->GetBinLowEdge(binLow)
                << " Up Edge: " << kTProjection->GetBinLowEdge(binLow + 1)
                << std::endl;
      binErr = (kTProjection->GetBinLowEdge(binLow + 1)
          - kTProjection->GetBinLowEdge(binLow)) / 2.;
      std::cout << "Percent of total: "
                << kTProjection->GetBinContent(binLow) * 100. / totalPairs
                << std::endl;
    } else {
      std::cout << "Percent of total: "
                << kTProjection->Integral(binLow, binUp) * 100. / totalPairs
                << std::endl;
    }
    std::cout << ikT << '\t' << "Bin Low: " << binLow << "(="
              << fKayTeeBins.at(ikT) << ")" << '\t' << "Bin Up: " << binUp
              << "(=" << fKayTeeBins.at(ikT + 1) << ")" << '\t'
              << kTProjection->GetMean() << '\t' << binErr << std::endl;
    fAveragekT->SetPoint(ikT, ikT + 1, kTProjection->GetMean());
    fAveragekT->SetPointError(ikT, 0, binErr);
  }
}

void DreamKayTee::SetSEMEReweightingRatio(const char* pathToFile,
                                          const char* HistNum,
                                          bool useRebinned) {
  TFile* CalibFile = TFile::Open(pathToFile, "read");
  if (!CalibFile) {
    std::cout << "No Calibration file found in " << pathToFile << std::endl;
  }
  TString HistName = useRebinned ? "Rebinned" : "FixShifted";
  SetSEMEReweightingRatio(
      (TH1F*) CalibFile->Get(
          TString::Format("hCk_%s_%s", HistName.Data(), HistNum)),
      (TH1F*) CalibFile->Get(TString::Format("hCk_Reweighted_%s", HistNum)),
      (TH1F*) CalibFile->Get(
          TString::Format("hCk_%sMeV_%s", HistName.Data(), HistNum)),
      (TH1F*) CalibFile->Get(TString::Format("hCk_ReweightedMeV_%s", HistNum)));
  return;
}

void DreamKayTee::SetSEMEReweightingRatio(TH1F* Rebinned, TH1F* Reweighted,
                                          TH1F* RebinnedMeV,
                                          TH1F* ReweightedMeV) {
  fSEMEReweighting = (TH1F*) Reweighted->Clone(
					       TString::Format("hCk_Reweighting"));
  fSEMEReweighting->Divide(Rebinned);
  fSEMEReweightingMeV = (TH1F*) ReweightedMeV->Clone(
						     TString::Format("hCk_ReweightingMeV"));
  fSEMEReweightingMeV->Divide(RebinnedMeV);
  return;
}


void DreamKayTee::SetSEmTMultDist(int iPart, int imT, TH2F* SEmTMult) {
  fSEmTMult[imT][iPart] = SEmTMult; 
  return; 
}


void DreamKayTee::SetMEmTMultDist(int iPart, int imT, TH2F* MEmTMult) {
  fMEmTMult[imT][iPart] = MEmTMult; 
  return; 
}


std::vector<DreamCF*> DreamKayTee::GetmTMultBinned(int imT, int Varcount) {
  std::vector<DreamCF*> outVec;
  if (!(imT < nmTBins)) {
    std::cout << "imT : " << imT << " too large for nmTBins: " << nmTBins << std::endl;
    return outVec;
  }
  //binning starts at 0, however since we need to do the binning steps starting at +1 this
  //is very find to put 0 here!!1elf!!
  std::vector<int> MultBins = {0, fSEmTMult[imT][0]->GetYaxis()->GetNbins()}; 
  for (auto it : fMultBins) {
    MultBins.push_back(fSEmTMult[imT][0]->GetYaxis()->FindBin(it)); 
  }
  std::sort(MultBins.begin(), MultBins.end());
  for (auto it = MultBins.begin(); it < MultBins.end()-1; ++it) { 
    auto itNextElem = it;
    itNextElem++;
    int binMin = *it + 1;
    int binMax = *itNextElem; 
    int counter = 0;
    
    TString ProjName = TString::Format("MultBin_%i_%i", binMin, binMax); 
    TH1F* PartSEProj = (TH1F*)
      fSEmTMult[imT][0]->ProjectionX(TString::Format("SEPart_%s",ProjName.Data()), binMin, binMax);
    
    TH1F* PartMEProj = (TH1F*)
      fMEmTMult[imT][0]->ProjectionX(TString::Format("MEPart_%s",ProjName.Data()), binMin, binMax);

    
    TH2F* PartSEMultLimited = (TH2F*)
      fSEmTMult[imT][0]->Clone(TString::Format("SEMultPart_%s",ProjName.Data())); 
    
    TH2F* PartMEMultLimited = (TH2F*)
      fMEmTMult[imT][0]->Clone(TString::Format("MEMultPart_%s",ProjName.Data())); 
    
    TH1F* AntiPartSEProj = (TH1F*)
      fSEmTMult[imT][1]->ProjectionX(TString::Format("SEAntiPart_%s",ProjName.Data()), binMin, binMax);
    
    TH1F* AntiPartMEProj = (TH1F*)
      fMEmTMult[imT][1]->ProjectionX(TString::Format("MEAntiPart_%s",ProjName.Data()), binMin, binMax);
    
    TH2F* AntiPartSEMultLimited = (TH2F*)
      fSEmTMult[imT][1]->Clone(TString::Format("SEMultAntiPart_%s",ProjName.Data())); 
    
    TH2F* AntiPartMEMultLimited = (TH2F*)
      fMEmTMult[imT][1]->Clone(TString::Format("MEMultAntiPart_%s",ProjName.Data())); 
    //now blind the whole thing outside of the multiplicity range!
    
    for (int iMult = 1; iMult < fSEmTMult[imT][0]->GetYaxis()->GetNbins()+1; ++iMult) {
      std::cout << iMult << std::endl; 
      if ( binMin <= iMult && iMult <= binMax) {
	continue;
      } else {
	for (int iks = 1; iks <= fSEmTMult[imT][0]->GetXaxis()->GetNbins(); ++iks) {
	  PartSEMultLimited->SetBinContent(iks, iMult, 0);
	  PartSEMultLimited->SetBinError(iks, iMult, 0);

	  PartMEMultLimited->SetBinContent(iks, iMult, 0);
	  PartMEMultLimited->SetBinError(iks, iMult, 0);
	 
	  AntiPartSEMultLimited->SetBinContent(iks, iMult, 0);
	  AntiPartSEMultLimited->SetBinError(iks, iMult, 0);
	 
	  AntiPartMEMultLimited->SetBinContent(iks, iMult, 0);
	  AntiPartMEMultLimited->SetBinError(iks, iMult, 0); 
	}
      }
    }
    
    DreamDist* pair = new DreamDist();
    pair->SetSEDist(PartSEProj,"");
    pair->SetSEMultDist(PartSEMultLimited,"");

    pair->SetMEDist(PartMEProj, "");
    pair->SetMEMultDist(PartMEMultLimited, "");
    
    DreamDist* Antipair = new DreamDist();
    Antipair->SetSEDist(AntiPartSEProj,"");
    Antipair->SetSEMultDist(AntiPartSEMultLimited,"");

    Antipair->SetMEDist(AntiPartMEProj, "");
    Antipair->SetMEMultDist(AntiPartMEMultLimited, "");

    DreamCF* CF_pp = new DreamCF();
    DreamPair* pp = new DreamPair("Part", 0.24, 0.34);
    DreamPair* ApAp = new DreamPair("AntiPart", 0.24, 0.34);
    pp->SetPair(pair); 
    ApAp->SetPair(Antipair);

    pp->ReweightMixedEvent(pp->GetPair(), 0.2, 0.9);
    ApAp->ReweightMixedEvent(ApAp->GetPair(), 0.2, 0.9);
    
  //There was a bug due to "our dear lovely ROOT" in this part which led to different results in Radius... Don't use it without discussing with Bernie

  /*pp->ShiftForEmpty(pp->GetPairReweighted(0));
    ApAp->ShiftForEmpty(ApAp->GetPairReweighted(0));
    pp->ShiftForEmpty(pp->GetPair());
    ApAp->ShiftForEmpty(ApAp->GetPair());
    pp->FixShift(pp->GetPairShiftedEmpty(0), ApAp->GetPairShiftedEmpty(0),
		ApAp->GetFirstBin());
    ApAp->FixShift(ApAp->GetPairShiftedEmpty(0), pp->GetPairShiftedEmpty(0),
		pp->GetFirstBin());
	*/

    pp->FixShift(pp->GetPairReweighted(0), ApAp->GetPairReweighted(0),
                 fFixShiftValue[imT], fFixShift[imT]);
     ApAp->FixShift(ApAp->GetPairReweighted(0), pp->GetPairReweighted(0),
                    fFixShiftValue[imT], fFixShift[imT]);
    
    pp->Rebin(pp->GetPairFixShifted(0), 8);
    ApAp->Rebin(ApAp->GetPairFixShifted(0), 8);

    //pp->Rebin(pp->GetPairFixShifted(0), 20);
    //ApAp->Rebin(ApAp->GetPairFixShifted(0), 20);

    CF_pp->SetPairs(pp,ApAp);
    CF_pp->GetCorrelations();
    CF_pp->WriteOutput(TString::Format("%s/CF_ppVar%u_mTBin_%i_%s.root",gSystem->pwd(),Varcount,
                                       imT,ProjName.Data()));
  }
  

  return outVec; 
}

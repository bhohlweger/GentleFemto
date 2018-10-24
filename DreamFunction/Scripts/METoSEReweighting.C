#include "DreamPlot.h"
#include "TStyle.h"

void GetQADistributions(const char* PairName, DreamDist* PairOrg,
                        DreamDist* PairRew) {
  gStyle->SetOptStat(false);
  DreamPlot::SetStyle();

  TString CanvasName = Form("Can%s", PairName);
  TCanvas* c1 = new TCanvas(CanvasName.Data(), CanvasName.Data(), 2000, 1000);
  c1->Divide(2, 2);

  TH1F* MERaw = PairOrg->GetMEDist();
  TH2F* SEMultRaw = PairOrg->GetSEMultDist();
  TH2F* MEMultRaw = PairOrg->GetMEMultDist();
  TH1F* CFRaw = PairOrg->GetCF();

  TH1F* MERew = PairRew->GetMEDist();
  TH2F* MEMultRew = PairRew->GetMEMultDist();
  TH1F* CFRew = PairRew->GetCF();
  c1->cd(1);
  TLegend leg = TLegend();

  TH1F* SEMultProjRaw = (TH1F*) SEMultRaw->ProjectionY("SEMultProjRaw",
                                                       SEMultRaw->FindBin(0.2),
                                                       SEMultRaw->FindBin(0.4));
  SEMultProjRaw->SetTitle("; Multiplicity Bin; Normalized counts");
  leg.AddEntry(SEMultProjRaw, "SE", "lep");
  SEMultProjRaw->Scale(1. / SEMultProjRaw->Integral());
  TH1F* MEMultProjRaw = (TH1F*) MEMultRaw->ProjectionY("MEMultProjRaw",
                                                       MEMultRaw->FindBin(0.2),
                                                       MEMultRaw->FindBin(0.4));
  MEMultProjRaw->Scale(1. / MEMultProjRaw->Integral());
  leg.AddEntry(MEMultProjRaw, "ME Unweighted", "lep");

  TH1F* MEMultProjRew = (TH1F*) MEMultRew->ProjectionY("MEMultProjRew",
                                                       MEMultRew->FindBin(0.2),
                                                       MEMultRew->FindBin(0.4));
  MEMultProjRew->Scale(1. / MEMultProjRew->Integral());
  leg.AddEntry(MEMultProjRew, "ME Weighted", "lep");

  DreamPlot::SetStyleHisto(SEMultProjRaw, 21, kGreen+2);
  DreamPlot::SetStyleHisto(MEMultProjRaw, 22, kRed+2);
  DreamPlot::SetStyleHisto(MEMultProjRew, 23, kBlue+2);
  SEMultProjRaw->DrawCopy();
  SEMultProjRaw->SetStats(false);
  MEMultProjRaw->DrawCopy("same");
  MEMultProjRew->DrawCopy("Same");
  leg.Draw("same");

  c1->cd(2);
  TH1F* MERawNorm = (TH1F*) MERaw->Clone(Form("%s_Norm", MERaw->GetName()));
  MERawNorm->Scale(1. / MERawNorm->Integral());
  MERawNorm->SetTitle("; #it{k}* (GeV/#it{c}); Normalized counts");
  MERawNorm->Draw();
  MERawNorm->SetStats(false);
  TH1F* MERewNorm = (TH1F*) MERew->Clone(Form("%s_Norm", MERew->GetName()));
  MERewNorm->Scale(1. / MERewNorm->Integral());
  MERewNorm->Draw("same");
  DreamPlot::SetStyleHisto(MERawNorm, 22, kRed+2);
  DreamPlot::SetStyleHisto(MERewNorm, 23, kBlue+2);
  TLegend leg2 = TLegend();
  leg2.AddEntry(MERawNorm, "ME Unweighted", "lep");
  leg2.AddEntry(MERewNorm, "ME Weighted", "lep");
  leg2.Draw("same");


  c1->cd(3);
  CFRaw->SetTitle(";#it{k}* (GeV/#it{c}); C(#it{k}*)");
  CFRaw->Draw();
  CFRaw->SetStats(false);
  CFRew->Draw("same");
  DreamPlot::SetStyleHisto(CFRaw, 22, kRed+2);
  DreamPlot::SetStyleHisto(CFRew, 23, kBlue+2);
  TLegend leg3 = TLegend();
  leg3.AddEntry(CFRaw, "CF Unweighted", "lep");
  leg3.AddEntry(CFRew, "CF Weighted", "lep");
  leg3.Draw("same");

  c1->cd(4);
  TH1F* Ratio = (TH1F*)CFRaw->Clone("Ratio");
  if (Ratio->Divide(CFRew)) {
    Ratio->Draw();
    Ratio->SetStats(false);
    DreamPlot::SetStyleHisto(Ratio, 20, kBlue+2);
    Ratio->SetTitle("; #it{k}* (GeV/#it{c}); C(#it{k}*)_{raw} / C(#it{k}*)_{reweighted}");
  }
  c1->Write();
  delete c1;
  delete SEMultProjRaw;
  delete MEMultProjRaw;
  delete MEMultProjRew;
  delete MERawNorm;
  delete MERewNorm;
  return;
}
void ReweightingQA(TList* PairList) {
  TList* UntouchedPairList = (TList*) PairList->FindObject("PairFixShifted");
  TList* PairRebinnedList = (TList*) PairList->FindObject("PairRebinned");
  TList* PairReweightedList = (TList*) PairList->FindObject("PairReweighted");
  const int nEntries = PairReweightedList->GetEntries() / 5;
  for (int iQA = 0; iQA < nEntries; ++iQA) {
    DreamDist* PairRew = new DreamDist();
    PairRew->SetSEDist((TH1F*) PairReweightedList->At(5 * iQA), "_");
    PairRew->SetSEMultDist((TH2F*) PairReweightedList->At(5 * iQA + 1), "_");
    PairRew->SetMEDist((TH1F*) PairReweightedList->At(5 * iQA + 2), "_");
    PairRew->SetMEMultDist((TH2F*) PairReweightedList->At(5 * iQA + 3), "_");
    PairRew->Calculate_CF(0.2, 0.4);
    //Get the corresponding unmodified pair, which was just reweighted!
    DreamDist* Pair = new DreamDist();
    TString SEName = PairRew->GetSEDist()->GetName();
    TString PartName = SEName(SEName.Index("Particle"), 19);
    TString CanName = PartName;

    TString test = "_Rebinned_";
    int iStart = SEName.Index(test.Data());

    std::cout << PartName.Data() << std::endl;
    if (iStart > 0) {
      iStart += test.Length();
      std::cout << SEName.Data() << std::endl;
      TString RebinName = SEName(iStart, 1);
      CanName += "_Rebinned_";
      CanName += RebinName;
      if (PairRebinnedList->FindObject(
          Form("SEDist_%s_Shifted_FixShifted_Rebinned_%s", PartName.Data(),
               RebinName.Data()))) {
        Pair->SetSEDist(
            (TH1F*) PairRebinnedList->FindObject(
                Form("SEDist_%s_Shifted_FixShifted_Rebinned_%s",
                     PartName.Data(), RebinName.Data())),
            "_");
        Pair->SetSEMultDist(
            (TH2F*) PairRebinnedList->FindObject(
                Form("SEMultDist_%s_Shifted_FixShifted_Rebinned_%s",
                     PartName.Data(), RebinName.Data())),
            "_");
        Pair->SetMEDist(
            (TH1F*) PairRebinnedList->FindObject(
                Form("MEDist_%s_Shifted_FixShifted_Rebinned_%s",
                     PartName.Data(), RebinName.Data())),
            "_");
        Pair->SetMEMultDist(
            (TH2F*) PairRebinnedList->FindObject(
                Form("MEMultDist_%s_Shifted_FixShifted_Rebinned_%s",
                     PartName.Data(), RebinName.Data())),
            "_");
        Pair->Calculate_CF(0.2, 0.4);
      } else {
        std::cout
            << "===========" << '\n' << "==Missing==" << '\n' << "==========="
            << '\n'
            << Form("SEDist_%s_Rebinned_%s", PartName.Data(), RebinName.Data())
            << std::endl;
      }
    } else {
      //else we don't rebinn, e.g. pp case, then take the original CF
      std::cout << Form("SEDist_%s", PartName.Data()) << std::endl;
      if (UntouchedPairList->FindObject(
          Form("SEDist_%s_Shifted_FixShifted", PartName.Data()))) {
        Pair->SetSEDist(
            (TH1F*) UntouchedPairList->FindObject(
                Form("SEDist_%s_Shifted_FixShifted", PartName.Data())),
            "_");
        Pair->SetSEMultDist(
            (TH2F*) UntouchedPairList->FindObject(
                Form("SEMultDist_%s_Shifted_FixShifted", PartName.Data())),
            "_");
        Pair->SetMEDist(
            (TH1F*) UntouchedPairList->FindObject(
                Form("MEDist_%s_Shifted_FixShifted", PartName.Data())),
            "_");
        Pair->SetMEMultDist(
            (TH2F*) UntouchedPairList->FindObject(
                Form("MEMultDist_%s_Shifted_FixShifted", PartName.Data())),
            "_");
        Pair->Calculate_CF(0.2, 0.4);
      } else {
        std::cout << "===========" << '\n' << "==Missing==" << '\n'
                  << "===========" << '\n'
                  << Form("SEDist_%s_Shifted_FixShifted", PartName.Data())
                  << std::endl;
      }
    }
    std::cout << CanName.Data() << std::endl;
    GetQADistributions(CanName.Data(), Pair, PairRew);
    delete Pair;
    delete PairRew;
  }
  return;
}

void METoSEReweighting(const char* foldername) {
  const char* filenames[4] = { "pp", "pXi", "pL", "LL" };
  for (int iFile = 0; iFile < 4; ++iFile) {
    TString FileName = Form("%sCFOutput_%s.root", foldername, filenames[iFile]);
    std::cout << FileName.Data() << std::endl;
    TFile* file = TFile::Open(FileName, "update");
    TList* PairDist = (TList*) file->Get("PairDist");
    if (PairDist) {
      ReweightingQA(PairDist);
    } else {
      file->ls();
    }
    TList* AntiPairDist = (TList*) file->Get("AntiPairDist");
    if (AntiPairDist) {
      ReweightingQA(AntiPairDist);
    }
//    TH1F* CFRaw;
//    TH1F* CFRew;
//    if (iFile == 0) {
//      CFRaw=(TH1F*)file->Get("hCkTotNormWeight");
//      CFRew = (TH1F*)file->Get("hCkReweighted_0");
//    }
//    else{
//      CFRaw=(TH1F*)file->Get("hCkReweighted_3");
//      CFRew = (TH1F*)file->Get("hCkReweighted_1");
//    }
//
//    TH1F* Ratio = (TH1F*)CFRaw->Clone(Form("Ratio_%i",iFile));
//    Ratio->Divide(CFRew);
//    TCanvas *c1=new TCanvas(Form("CanRatio%i",iFile),Form("CanRatio%i",iFile),2000,1000);
//    c1->Divide(2,1);
//    c1->cd(1);
//    CFRaw->SetLineColor(3);
//    CFRaw->DrawCopy();
//    CFRew->SetLineColor(4);
//    CFRew->DrawCopy("same");
//    c1->cd(2);
//    Ratio->DrawCopy();
//    c1->Write();
//    delete Ratio;
  }
}

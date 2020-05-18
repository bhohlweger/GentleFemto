#include <iostream>
#include <fstream>
#include <string>

#include "TFile.h"
#include "TGraphErrors.h"

using namespace std;

int main (int argc, char *argv[]) {
  TFile* out = TFile::Open("/home/schmollweger/cernbox/WaveFunctions/HAL_pXiPaper/Sasaki_pXi_potentials.root", "recreate"); 
  ifstream myfile ("/home/schmollweger/cernbox/WaveFunctions/HAL_pXiPaper/Sasaki_pXi_potentials.txt"); 
  //TFile* out = TFile::Open("/home/schmollweger/cernbox/WaveFunctions/HAL_pXiPaper/Sasaki_pXi_operator.root", "recreate"); 
  //ifstream myfile ("/home/schmollweger/cernbox/WaveFunctions/HAL_pXiPaper/Sasaki_pXi_operator.txt"); 
  int counter = 0;
  TGraphErrors* gr_I0S0_t11 = new TGraphErrors();
  gr_I0S0_t11->SetName("gr_I0S0_t11");

  TGraphErrors* gr_I0S1_t11 = new TGraphErrors();
  gr_I0S1_t11->SetName("gr_I0S1_t11");

  TGraphErrors* gr_I1S0_t11 = new TGraphErrors();
  gr_I1S0_t11->SetName("gr_I1S0_t11");

  TGraphErrors* gr_I1S1_t11 = new TGraphErrors();
  gr_I1S1_t11->SetName("gr_I1S1_t11");
  
  TGraphErrors* gr_I0S0_t12 = new TGraphErrors();
  gr_I0S0_t12->SetName("gr_I0S0_t12");

  TGraphErrors* gr_I0S1_t12 = new TGraphErrors();
  gr_I0S1_t12->SetName("gr_I0S1_t12");

  TGraphErrors* gr_I1S0_t12 = new TGraphErrors();
  gr_I1S0_t12->SetName("gr_I1S0_t12");

  TGraphErrors* gr_I1S1_t12 = new TGraphErrors();
  gr_I1S1_t12->SetName("gr_I1S1_t12");

  TGraphErrors* gr_I0S0_t13 = new TGraphErrors();
  gr_I0S0_t13->SetName("gr_I0S0_t13");

  TGraphErrors* gr_I0S1_t13 = new TGraphErrors();
  gr_I0S1_t13->SetName("gr_I0S1_t13");

  TGraphErrors* gr_I1S0_t13 = new TGraphErrors();
  gr_I1S0_t13->SetName("gr_I1S0_t13");

  TGraphErrors* gr_I1S1_t13 = new TGraphErrors();
  gr_I1S1_t13->SetName("gr_I1S1_t13");


  gr_I0S0_t11->SetLineColor(9);
  gr_I0S1_t11->SetLineColor(9);
  gr_I1S0_t11->SetLineColor(9);
  gr_I1S1_t11->SetLineColor(9); 
  
  gr_I0S0_t12->SetLineColor(2);
  gr_I0S1_t12->SetLineColor(2);
  gr_I1S0_t12->SetLineColor(2);
  gr_I1S1_t12->SetLineColor(2); 
  
  gr_I0S0_t13->SetLineColor(8);
  gr_I0S1_t13->SetLineColor(8);
  gr_I1S0_t13->SetLineColor(8);
  gr_I1S1_t13->SetLineColor(8); 
  
  if (myfile.is_open()) {
    string rad;
    string I0S0t11, dI0S0t11, I0S1t11, dI0S1t11, I1S0t11, dI1S0t11, I1S1t11, dI1S1t11;
    string I0S0t12, dI0S0t12, I0S1t12, dI0S1t12, I1S0t12, dI1S0t12, I1S1t12, dI1S1t12;
    string I0S0t13, dI0S0t13, I0S1t13, dI0S1t13, I1S0t13, dI1S0t13, I1S1t13, dI1S1t13; 
    while ( counter<2.08e4) {
      myfile >> rad >> I0S0t11 >> dI0S0t11 >> I0S1t11 >> dI0S1t11 >> I1S0t11 >> dI1S0t11 >> I1S1t11 >> dI1S1t11 >> I0S0t12 >> dI0S0t12 >> I0S1t12 >> dI0S1t12 >> I1S0t12 >> dI1S0t12 >> I1S1t12 >> dI1S1t12 >> I0S0t13 >> dI0S0t13 >> I0S1t13 >> dI0S1t13 >> I1S0t13 >> dI1S0t13 >> I1S1t13 >> dI1S1t13;
      
      gr_I0S0_t11->SetPoint(counter, stof(rad), stof(I0S0t11)); 
      gr_I0S0_t11->SetPointError(counter, 0, stof(dI0S0t11));

      gr_I0S1_t11->SetPoint(counter, stof(rad), stof(I0S1t11)); 
      gr_I0S1_t11->SetPointError(counter, 0, stof(dI0S1t11));

      gr_I1S0_t11->SetPoint(counter, stof(rad), stof(I1S0t11)); 
      gr_I1S0_t11->SetPointError(counter, 0, stof(dI1S0t11));

      gr_I1S1_t11->SetPoint(counter, stof(rad), stof(I1S1t11)); 
      gr_I1S1_t11->SetPointError(counter, 0, stof(dI1S1t11));

      gr_I0S0_t12->SetPoint(counter, stof(rad), stof(I0S0t12)); 
      gr_I0S0_t12->SetPointError(counter, 0, stof(dI0S0t12));

      gr_I0S1_t12->SetPoint(counter, stof(rad), stof(I0S1t12)); 
      gr_I0S1_t12->SetPointError(counter, 0, stof(dI0S1t12));

      gr_I1S0_t12->SetPoint(counter, stof(rad), stof(I1S0t12)); 
      gr_I1S0_t12->SetPointError(counter, 0, stof(dI1S0t12));

      gr_I1S1_t12->SetPoint(counter, stof(rad), stof(I1S1t12)); 
      gr_I1S1_t12->SetPointError(counter, 0, stof(dI1S1t12));

      gr_I0S0_t13->SetPoint(counter, stof(rad), stof(I0S0t13)); 
      gr_I0S0_t13->SetPointError(counter, 0, stof(dI0S0t13));

      gr_I0S1_t13->SetPoint(counter, stof(rad), stof(I0S1t13)); 
      gr_I0S1_t13->SetPointError(counter, 0, stof(dI0S1t13));

      gr_I1S0_t13->SetPoint(counter, stof(rad), stof(I1S0t13)); 
      gr_I1S0_t13->SetPointError(counter, 0, stof(dI1S0t13));

      gr_I1S1_t13->SetPoint(counter, stof(rad), stof(I1S1t13)); 
      gr_I1S1_t13->SetPointError(counter, 0, stof(dI1S1t13));
      counter++; 
    }
    myfile.close();
  } else {
    cout << "Unable to open file";
  }
  out->cd();
  gr_I0S0_t11->Write("gr_I0S0_t11");
  gr_I0S1_t11->Write("gr_I0S1_t11");
  gr_I1S0_t11->Write("gr_I1S0_t11");
  gr_I1S1_t11->Write("gr_I1S1_t11"); 

  gr_I0S0_t12->Write("gr_I0S0_t12");
  gr_I0S1_t12->Write("gr_I0S1_t12");
  gr_I1S0_t12->Write("gr_I1S0_t12");
  gr_I1S1_t12->Write("gr_I1S1_t12"); 

  gr_I0S0_t13->Write("gr_I0S0_t13");
  gr_I0S1_t13->Write("gr_I0S1_t13");
  gr_I1S0_t13->Write("gr_I1S0_t13");
  gr_I1S1_t13->Write("gr_I1S1_t13"); 

  out->Close(); 
  return 0; 
}


/* dumping ground 
   
      std::cout << "COUNTER: " << counter << "   " << strtof(rad.c_str(), NULL) << "    " << strtof(I0S0t11.c_str(), NULL) << "    " << strtof(dI0S0t11.c_str(), NULL) << "    " << strtof(I0S1t11.c_str(), NULL) << "    " << strtof(dI0S1t11.c_str(), NULL) << "    " << strtof(I1S0t11.c_str(), NULL) << "    " << strtof(dI1S0t11.c_str(), NULL) << "    " << strtof(I1S1t11.c_str(), NULL) << "    " << strtof(dI1S1t11.c_str(), NULL) << "    " << strtof(I0S0t12.c_str(), NULL) << "    " << strtof(dI0S0t12.c_str(), NULL) << "    " << strtof(I0S1t12.c_str(), NULL) << "    " << strtof(dI0S1t12.c_str(), NULL) << "    " << strtof(I1S0t12.c_str(), NULL) << "    " << strtof(dI1S0t12.c_str(), NULL) << "    " << strtof(I1S1t12.c_str(), NULL) << "    " << strtof(dI1S1t12.c_str(), NULL) << "    " << strtof(I0S0t13.c_str(), NULL) << "    " << strtof(dI0S0t13.c_str(), NULL) << "    " << strtof(I0S1t13.c_str(), NULL) << "    " << strtof(dI0S1t13.c_str(), NULL) << "    " << strtof(I1S0t13.c_str(), NULL) << "    " << strtof(dI1S0t13.c_str(), NULL) << "    " << strtof(I1S1t13.c_str(), NULL) << "    " << strtof(dI1S1t13.c_str(), NULL) << std::endl;

      

      std::cout << counter << " \t " << rad << " \t " << I0S0t11 << " \t " << dI0S0t11 << " \t " << I0S1t11 << " \t " << dI0S1t11 << " \t " << I1S0t11 << " \t " << dI1S0t11 << " \t " << I1S1t11 << " \t " << dI1S1t11 << " \t " << I0S0t12 << " \t " << dI0S0t12 << " \t " << I0S1t12 << " \t " << dI0S1t12 << " \t " << I1S0t12 << " \t " << dI1S0t12 << " \t " << I1S1t12 << " \t " << dI1S1t12 << " \t " << I0S0t13 << " \t " << dI0S0t13 << " \t " << I0S1t13 << " \t " << dI0S1t13 << " \t " << I1S0t13 << " \t " << dI1S0t13 << " \t " << I1S1t13 << " \t " << dI1S1t13 << std::endl; 

 */ 

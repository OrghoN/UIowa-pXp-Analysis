#include <iostream>
#include <cstdio>
//
#include <TString.h>
#include <TH1F.h>
#include <TTree.h>
#include <TFile.h>
#include <TROOT.h>
#include <TMath.h>
#include <TH2F.h>
#include <TClass.h>
//
//important
#include <TPad.h>
#include <TStyle.h>
#include <TPave.h>
#include <TPaletteAxis.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TColor.h>
//
#include "TSystem.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include <TDirectory.h>
#include <TSystem.h>
//
using namespace std;

void norm2()
{
    TFile *f = TFile::Open("t0RP110relui4.root");

    //hm2rec2OS_ttbb2varbin

  TH1F *h1 = (TH1F*)(hm2rec2OS_ttbb2varbin->Clone("h1"));
  TH1F *h2 = (TH1F*)(hm2rec2OS_diag2varbin->Clone("h2"));
 
    //Method 1:

  int j;
  for(j=0; j<266; ++j) {

    double num1 = h1->GetBinContent(j);
    double num2 = h2->GetBinContent(j);
    cout << " num1 = " << num1 << "\n";
    cout << " num1 = " << num2 << "\n";

    double den1 = h1->GetBinWidth(j);
    double den2 = h2->GetBinWidth(j);
    cout << "den1 = " << den1 << "\n";
    cout << "den2 = " << den2 << "\n";

 double value1 = 0;
 double value2 = 0;
 
 if (den1!=0)
  {
     value1 = num1/den1;
     cout << "value1 = " << value1 << "\n";
     h1->SetBinContent(j, value1);
  }
 
 if (den2!=0)
  {
     value2 = num2/den2;
     cout << "value2 = " << value2 << "\n";
     h2->SetBinContent(j, value2);
  }
  }
  
  //gPad->SetLogx();
  //gPad->SetLogy();
  h1->GetXaxis()->SetTitle(M_{#pi#pi} (GeV/c^{2}));
  h1->GetYaxis()->SetTitle("# of Events/0.02,0.025,0.05GeV/c^{2}");
  h1->GetXaxis()->CenterTitle();
  h1->GetYaxis()->CenterTitle();
  h1->GetXaxis()->SetTitleOffset(1.2);
  h1->GetYaxis()->SetTitleOffset(1.5);

  //h1->SetFillColor(3);
  //h2->SetFillColor(2);
  h1->SetLineColor(kBlue);
  h2->SetLineColor(kRed);

  h1->Draw("HIST");
  h2->Draw("HISTSAME");

 //the end
}

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
//
#include "TSystem.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include <TDirectory.h>
#include <TSystem.h>
//
using namespace std;

void makeplot()
{
  //...using all 2015 data
  //TFile f("t0RPrelui2015all.root");
  //TFile *f = TFile::Open("t0RPrelui2015all.root");
  //TFile *f = TFile::Open("t0RP2trackrelui.root");
  //TFile *f = TFile::Open("t0RP4trackrelui.root");
  TFile *f = TFile::Open("t0RP4track2relui.root");
  //TFile *f = TFile::Open("t0RP42relui.root");
  //TFile *f = TFile::Open("t0RP25relui.root");
  //TFile *f = TFile::Open("t0RP69relui.root");
  //
  cout << "                                 " << endl;
  cout << " * * *  TOTEM making plots  * * *" << endl;
  cout << "             in root:            " << endl;
  cout << "          .x makeplot.cc         " << endl;
  //cout << "        .L makeplot.cc+        " << endl;
  //cout << "        .L makeplot.cc++       " << endl;
  cout << "                                 " << endl;
  cout << "          Iowa-FNAL Team         " << endl;
  cout << "                                 " << endl;
  //
  cout << " making 1D plot, type: " << endl;
  cout << " plot1d(histoname,\"x-title\",\"binning-with-units\",\"y-lin/log\",\"x-lin/log\")" << endl;
  //
  cout << " making 2D plot, type: " << endl;
  cout << " plot2d(histoname,\"y-title\",\"x-title\",\"y-lin/log\",\"x-lin/log\")" << endl;
  cout << "                                 " << endl;  
  cout << " for the flash mode, type: " << endl;
  cout << " plotall()" << endl;
  cout << " no pop-ups:                     " << endl;
  cout << " root -b makeplot.cc             " << endl;
  cout << "                                 " << endl;
  //
  cout << " 1D-example: plot1d(rp_y_125,\"Y(cm)\",\"0.1cm\",\"lin\",\"lin\")" << endl;
  cout << "     output: rp_y_125.png                   " << endl;
  cout << " 2D-example: plot2d(hdedx,\"dE/dx (MeV/cm)\",\"p (GeV/c)\",\"lin\",\"log\")" << endl;
  cout << "     output: hdedx.png                      " << endl;
  cout << "                                        " << endl;
  cout << "          using all 2015 data           " << endl;
  cout << "    data file: t0RPrelui2015all.root    " << endl;
  cout << "                                        " << endl;
  cout << " important:                             " << endl;
  cout << " change root file at the ::Open statement for new data " << endl;
  cout << "                                        " << endl;
  //
  //TString myplot="makeplot.cc";
  //gROOT->Macro(myplot.Data());
}

void plot1d(TH1F* string1, TString xtitle, TString binning, TString scaley1, TString scalex1)
{
  TCanvas * c1 = new TCanvas("c1","c1",1200,800);
  //TCanvas *c1 = new TCanvas("c1","c1",600,400);
  string1->Draw();
  //stats
  //TPaveStats* st + string1 = (TPaveStats*)string1->FindObject("stats");
  //st + string1 ->SetX2NDC(0.89);
  //st + string1 ->SetX1NDC(0.61);
  string1->GetXaxis()->SetTitle(xtitle);
  string1->GetYaxis()->SetTitle("# of Events/" + binning);
  string1->GetXaxis()->CenterTitle();
  string1->GetYaxis()->CenterTitle();
  string1->GetXaxis()->SetTitleOffset(1.2);
  string1->GetYaxis()->SetTitleOffset(1.5);
  c1->Update();
  if (scalex1 == "log"){
    gPad->SetLogx();
    }
  if (scaley1 == "log"){
    gPad->SetLogy();
    }
  c1->Modified();
  //
  //FULL SCREEN BEFORE SAVING...need to check the resolution...ok
  //
  TString hname1 = string1->GetName();
  c1->SaveAs(hname1 + ".png");
  //back to linear
  //gPad->SetLogx(0);
  //gPad->SetLogy(0);
  //c1->Clear();
  //c1->Close();
}

void plot2d(TH2F* string2, TString ytitle, TString xtitle, TString scaley2, TString scalex2)
{
  TCanvas * c1 = new TCanvas("c1","c1",1200,800);
  Int_t nentries;
  string2->Draw("SCATCONT0Z");
  //stats
  //TPaveStats *st + string1 = (TPaveStats*)string1->FindObject("stats");
  //st + string1 ->SetX2NDC(0.89);
  //st + string1 ->SetX1NDC(0.61);
  //
  string2->GetXaxis()->SetTitle(xtitle);
  string2->GetYaxis()->SetTitle(ytitle);
  string2->GetXaxis()->CenterTitle();
  string2->GetYaxis()->CenterTitle();
  string2->GetXaxis()->SetTitleOffset(1.2);
  string2->GetYaxis()->SetTitleOffset(1.5);
  //
  gStyle->SetPalette(107,0);
  gStyle->SetNumberContours(100);
  //
  c1->Update();
  nentries = string2->GetEntries();
  if(nentries){
  TPaletteAxis *palette = (TPaletteAxis*)string2->GetListOfFunctions()->FindObject("palette");
  c1->Update();
  //palette->SetX1NDC(0.875);        
  //palette->SetY1NDC(0.115);
  //palette->SetX2NDC(0.925);
  //palette->SetY2NDC(0.675);
  palette->SetX1NDC(0.885);   
  palette->SetY1NDC(0.115);
  palette->SetX2NDC(0.915);
  palette->SetY2NDC(0.665);
  }
  string2->GetZaxis()->SetLabelSize(0.02);
  string2->GetZaxis()->SetRangeUser(0, string2->GetMaximum());
  string2->GetZaxis()->SetNdivisions(10);
  c1->Modified();
  if (scalex2 == "log"){
    gPad->SetLogx();
    string2->GetXaxis()->SetMoreLogLabels();
   }
  if (scaley2 == "log"){
    gPad->SetLogy();
    //string2->GetYaxis()->SetMoreLogLabels();    
   }
  c1->Modified();
  //
  //FULL SCREEN BEFORE SAVING...need to check the resolution...ok
  // 
  TString hname2 = string2->GetName();
  c1->SaveAs(hname2 + ".png");
  //gPad->SetLogx(0);
  //gPad->SetLogy(0);
  //c1->Clear();
  //c1->Close();
}

void plotall()
{
  //  
  //TString myplot="makeplot.cc";
  //gROOT->Macro(myplot.Data());
  //
  //TH1* readThis = 0;
  //file->GetObject("hpx", readThis);
  //
  //TFile *f = TFile::Open("t0RPrelui2015all.root");
  //plot1d(f->GetObject("rp_y_125",readThis),"X (cm)","0.1cm","lin","lin");
  //plot2d(f->GetObject("rp_yx_125",readThis),"Y (cm)","X (cm)","lin","lin");
  //
//Histogram variable
TH1F *hnconf; gFile->GetObject("hnconf", hnconf);
TH1F *rp_x_020; gFile->GetObject("rp_x_020", rp_x_020);
TH1F *rp_x_021; gFile->GetObject("rp_x_021", rp_x_021);
TH1F *rp_x_024; gFile->GetObject("rp_x_024", rp_x_024);
TH1F *rp_x_025; gFile->GetObject("rp_x_025", rp_x_025);
TH1F *rp_x_120; gFile->GetObject("rp_x_120", rp_x_120);
TH1F *rp_x_121; gFile->GetObject("rp_x_121", rp_x_121);
TH1F *rp_x_124; gFile->GetObject("rp_x_124", rp_x_124);
TH1F *rp_x_125; gFile->GetObject("rp_x_125", rp_x_125);
TH1F *rp_y_020; gFile->GetObject("rp_y_020", rp_y_020);
TH1F *rp_y_021; gFile->GetObject("rp_y_021", rp_y_021);
TH1F *rp_y_024; gFile->GetObject("rp_y_024", rp_y_024);
TH1F *rp_y_025; gFile->GetObject("rp_y_025", rp_y_025);
TH1F *rp_y_120; gFile->GetObject("rp_y_120", rp_y_120);
TH1F *rp_y_121; gFile->GetObject("rp_y_121", rp_y_121);
TH1F *rp_y_124; gFile->GetObject("rp_y_124", rp_y_124);
TH1F *rp_y_125; gFile->GetObject("rp_y_125", rp_y_125);
TH1F *rp2_x_020; gFile->GetObject("rp2_x_020", rp2_x_020);
TH1F *rp2_x_021; gFile->GetObject("rp2_x_021", rp2_x_021);
TH1F *rp2_x_024; gFile->GetObject("rp2_x_024", rp2_x_024);
TH1F *rp2_x_025; gFile->GetObject("rp2_x_025", rp2_x_025);
TH1F *rp2_x_120; gFile->GetObject("rp2_x_120", rp2_x_120);
TH1F *rp2_x_121; gFile->GetObject("rp2_x_121", rp2_x_121);
TH1F *rp2_x_124; gFile->GetObject("rp2_x_124", rp2_x_124);
TH1F *rp2_x_125; gFile->GetObject("rp2_x_125", rp2_x_125);
TH1F *rp2_y_020; gFile->GetObject("rp2_y_020", rp2_y_020);
TH1F *rp2_y_021; gFile->GetObject("rp2_y_021", rp2_y_021);
TH1F *rp2_y_024; gFile->GetObject("rp2_y_024", rp2_y_024);
TH1F *rp2_y_025; gFile->GetObject("rp2_y_025", rp2_y_025);
TH1F *rp2_y_120; gFile->GetObject("rp2_y_120", rp2_y_120);
TH1F *rp2_y_121; gFile->GetObject("rp2_y_121", rp2_y_121);
TH1F *rp2_y_124; gFile->GetObject("rp2_y_124", rp2_y_124);
TH1F *rp2_y_125; gFile->GetObject("rp2_y_125", rp2_y_125);
TH1F *thyEla; gFile->GetObject("thyEla", thyEla);
TH1F *thxEla; gFile->GetObject("thxEla", thxEla);
TH1F *thyEla_diag; gFile->GetObject("thyEla_diag", thyEla_diag);
TH1F *thxEla_diag; gFile->GetObject("thxEla_diag", thxEla_diag);
TH1F *thyEla_ttbb; gFile->GetObject("thyEla_ttbb", thyEla_ttbb);
TH1F *thxEla_ttbb; gFile->GetObject("thxEla_ttbb", thxEla_ttbb);
TH1F *proton_right_xi; gFile->GetObject("proton_right_xi", proton_right_xi);
TH1F *proton_left_xi; gFile->GetObject("proton_left_xi", proton_left_xi);
TH1F *proton_right_logXi; gFile->GetObject("proton_right_logXi", proton_right_logXi);
TH1F *proton_left_logXi; gFile->GetObject("proton_left_logXi", proton_left_logXi);
TH1F *proton_right_t; gFile->GetObject("proton_right_t", proton_right_t);
TH1F *proton_left_t; gFile->GetObject("proton_left_t", proton_left_t);
TH1F *proton_right_t_diag; gFile->GetObject("proton_right_t_diag", proton_right_t_diag);
TH1F *proton_left_t_diag; gFile->GetObject("proton_left_t_diag", proton_left_t_diag);
TH1F *proton_right_t_ttbb; gFile->GetObject("proton_right_t_ttbb", proton_right_t_ttbb);
TH1F *proton_left_t_ttbb; gFile->GetObject("proton_left_t_ttbb", proton_left_t_ttbb);
TH1F *eHF; gFile->GetObject("eHF", eHF);
TH1F *nHF; gFile->GetObject("nHF", nHF);
TH1F *totem_py; gFile->GetObject("totem_py", totem_py);
TH1F *totem_px; gFile->GetObject("totem_px", totem_px);
TH1F *totem_pyy; gFile->GetObject("totem_pyy", totem_pyy);
TH1F *totem_pxx; gFile->GetObject("totem_pxx", totem_pxx);
TH1F *proton_dx0; gFile->GetObject("proton_dx0", proton_dx0);
TH1F *hLS; gFile->GetObject("hLS", hLS);
TH1F *htopo; gFile->GetObject("htopo", htopo);
TH1F *hthyEla2_diag; gFile->GetObject("hthyEla2_diag", hthyEla2_diag);
TH1F *hthxEla2_diag; gFile->GetObject("hthxEla2_diag", hthxEla2_diag);
TH1F *hthyEla2_ttbb; gFile->GetObject("hthyEla2_ttbb", hthyEla2_ttbb);
TH1F *hthxEla2_ttbb; gFile->GetObject("hthxEla2_ttbb", hthxEla2_ttbb);
//Histogram variable
TH2F *rp_yx_020; gFile->GetObject("rp_yx_020", rp_yx_020);
TH2F *rp_yx_021; gFile->GetObject("rp_yx_021", rp_yx_021);
TH2F *rp_yx_024; gFile->GetObject("rp_yx_024", rp_yx_024);
TH2F *rp_yx_025; gFile->GetObject("rp_yx_025", rp_yx_025);
TH2F *rp_yx_120; gFile->GetObject("rp_yx_120", rp_yx_120);
TH2F *rp_yx_121; gFile->GetObject("rp_yx_121", rp_yx_121);
TH2F *rp_yx_124; gFile->GetObject("rp_yx_124", rp_yx_124);
TH2F *rp_yx_125; gFile->GetObject("rp_yx_125", rp_yx_125);
TH2F *rp2_yx_020; gFile->GetObject("rp2_yx_020", rp2_yx_020);
TH2F *rp2_yx_021; gFile->GetObject("rp2_yx_021", rp2_yx_021);
TH2F *rp2_yx_024; gFile->GetObject("rp2_yx_024", rp2_yx_024);
TH2F *rp2_yx_025; gFile->GetObject("rp2_yx_025", rp2_yx_025);
TH2F *rp2_yx_120; gFile->GetObject("rp2_yx_120", rp2_yx_120);
TH2F *rp2_yx_121; gFile->GetObject("rp2_yx_121", rp2_yx_121);
TH2F *rp2_yx_124; gFile->GetObject("rp2_yx_124", rp2_yx_124);
TH2F *rp2_yx_125; gFile->GetObject("rp2_yx_125", rp2_yx_125);
TH2F *proton_x0_RvsL; gFile->GetObject("proton_x0_RvsL", proton_x0_RvsL);
//Histogram variable
TH1F *hlooper; gFile->GetObject("hlooper", hlooper);
TH1F *hpt; gFile->GetObject("hpt", hpt);
TH1F *heta; gFile->GetObject("heta", heta);
TH1F *hphi; gFile->GetObject("hphi", hphi);
TH1F *hptP; gFile->GetObject("hptP", hptP);
TH1F *hetaP; gFile->GetObject("hetaP", hetaP);
TH1F *hphiP; gFile->GetObject("hphiP", hphiP);
TH1F *hptM; gFile->GetObject("hptM", hptM);
TH1F *hetaM; gFile->GetObject("hetaM", hetaM);
TH1F *hphiM; gFile->GetObject("hphiM", hphiM);
TH1F *hptRes; gFile->GetObject("hptRes", hptRes);
TH1F *hetaRes; gFile->GetObject("hetaRes", hetaRes);
TH1F *hphiRes; gFile->GetObject("hphiRes", hphiRes);
TH1F *hthyEla_diag; gFile->GetObject("hthyEla_diag", hthyEla_diag);
TH1F *hthxEla_diag; gFile->GetObject("hthxEla_diag", hthxEla_diag);
TH1F *hthyEla_ttbb; gFile->GetObject("hthyEla_ttbb", hthyEla_ttbb);
TH1F *hthxEla_ttbb; gFile->GetObject("hthxEla_ttbb", hthxEla_ttbb);
TH1F *hntrk0; gFile->GetObject("hntrk0", hntrk0);
TH1F *hntrk; gFile->GetObject("hntrk", hntrk);
TH1F *hntrkvtx; gFile->GetObject("hntrkvtx", hntrkvtx);
TH1F *hntrkntrkvtx2; gFile->GetObject("hntrkntrkvtx2", hntrkntrkvtx2);
TH1F *hntrk2ntrkvtx; gFile->GetObject("hntrk2ntrkvtx", hntrk2ntrkvtx);
//Histogram variable
TH2F *hntrkntrkvtx; gFile->GetObject("hntrkntrkvtx", hntrkntrkvtx);
//Histogram variable
TH1F *hvtx; gFile->GetObject("hvtx", hvtx);
TH1F *hvtx2; gFile->GetObject("hvtx2", hvtx2);
TH1F *hvtx3; gFile->GetObject("hvtx3", hvtx3);
TH1F *hnvtx; gFile->GetObject("hnvtx", hnvtx);
TH1F *hvtxx; gFile->GetObject("hvtxx", hvtxx);
TH1F *hvtxy; gFile->GetObject("hvtxy", hvtxy);
TH1F *hvtxz; gFile->GetObject("hvtxz", hvtxz);
TH1F *hvtxchi2; gFile->GetObject("hvtxchi2", hvtxchi2);
TH1F *hvtxchi2fin; gFile->GetObject("hvtxchi2fin", hvtxchi2fin);
TH1F *heHF; gFile->GetObject("heHF", heHF);
TH1F *hnHF; gFile->GetObject("hnHF", hnHF);
TH1F *hxiL; gFile->GetObject("hxiL", hxiL);
TH1F *hxiR; gFile->GetObject("hxiR", hxiR);
TH1F *hxiL2; gFile->GetObject("hxiL2", hxiL2);
TH1F *hxiR2; gFile->GetObject("hxiR2", hxiR2);
TH1F *hm; gFile->GetObject("hm", hm);
TH1F *hmxicut; gFile->GetObject("hmxicut", hmxicut);
TH1F *hm2rec; gFile->GetObject("hm2rec", hm2rec);
TH1F *hm2recbis; gFile->GetObject("hm2recbis", hm2recbis);
TH1F *hm2recPP; gFile->GetObject("hm2recPP", hm2recPP);
TH1F *hm2recKK; gFile->GetObject("hm2recKK", hm2recKK);
TH1F *hm2recMM; gFile->GetObject("hm2recMM", hm2recMM);
TH1F *hm2recEE; gFile->GetObject("hm2recEE", hm2recEE);
TH1F *hm2recOS; gFile->GetObject("hm2recOS", hm2recOS);
TH1F *hm2recSS; gFile->GetObject("hm2recSS", hm2recSS);
TH1F *hm2recOS_diag; gFile->GetObject("hm2recOS_diag", hm2recOS_diag);
TH1F *hm2recSS_diag; gFile->GetObject("hm2recSS_diag", hm2recSS_diag);
TH1F *hm2recOS_ttbb; gFile->GetObject("hm2recOS_ttbb", hm2recOS_ttbb);
TH1F *hm2recSS_ttbb; gFile->GetObject("hm2recSS_ttbb", hm2recSS_ttbb);
TH1F *hm2rec2OS; gFile->GetObject("hm2rec2OS", hm2rec2OS);
TH1F *hm2rec2SS; gFile->GetObject("hm2rec2SS", hm2rec2SS);
TH1F *hm2rec2OS_diag; gFile->GetObject("hm2rec2OS_diag", hm2rec2OS_diag);
TH1F *hm2rec2SS_diag; gFile->GetObject("hm2rec2SS_diag", hm2rec2SS_diag);
TH1F *hm2rec2OS_ttbb; gFile->GetObject("hm2rec2OS_ttbb", hm2rec2OS_ttbb);
TH1F *hm2rec2SS_ttbb; gFile->GetObject("hm2rec2SS_ttbb", hm2rec2SS_ttbb);
TH1F *hm2rec2OS_diag_trkP; gFile->GetObject("hm2rec2OS_diag_trkP", hm2rec2OS_diag_trkP);
TH1F *hm2rec2OS_diag_trkM; gFile->GetObject("hm2rec2OS_diag_trkM", hm2rec2OS_diag_trkM);
TH1F *hm2rec2OS_ttbb_trkP; gFile->GetObject("hm2rec2OS_ttbb_trkP", hm2rec2OS_ttbb_trkP);
TH1F *hm2rec2OS_ttbb_trkM; gFile->GetObject("hm2rec2OS_ttbb_trkM", hm2rec2OS_ttbb_trkM);
TH1F *hm2rec2OS_diag_pypxP; gFile->GetObject("hm2rec2OS_diag_pypxP", hm2rec2OS_diag_pypxP);
TH1F *hm2rec2OS_diag_pypxM; gFile->GetObject("hm2rec2OS_diag_pypxM", hm2rec2OS_diag_pypxM);
TH1F *hm2rec2OS_ttbb_pypxP; gFile->GetObject("hm2rec2OS_ttbb_pypxP", hm2rec2OS_ttbb_pypxP);
TH1F *hm2rec2OS_ttbb_pypxM; gFile->GetObject("hm2rec2OS_ttbb_pypxM", hm2rec2OS_ttbb_pypxM);
TH1F *hm2rec3OS; gFile->GetObject("hm2rec3OS", hm2rec3OS);
TH1F *hm2rec3SS; gFile->GetObject("hm2rec3SS", hm2rec3SS);
TH1F *hm2rec3OS_diag; gFile->GetObject("hm2rec3OS_diag", hm2rec3OS_diag);
TH1F *hm2rec3SS_diag; gFile->GetObject("hm2rec3SS_diag", hm2rec3SS_diag);
TH1F *hm2rec3OS_ttbb; gFile->GetObject("hm2rec3OS_ttbb", hm2rec3OS_ttbb);
TH1F *hm2rec3SS_ttbb; gFile->GetObject("hm2rec3SS_ttbb", hm2rec3SS_ttbb);
TH1F *hm2rec3OS_diag_trkP; gFile->GetObject("hm2rec3OS_diag_trkP", hm2rec3OS_diag_trkP);
TH1F *hm2rec3OS_diag_trkM; gFile->GetObject("hm2rec3OS_diag_trkM", hm2rec3OS_diag_trkM);
TH1F *hm2rec3OS_ttbb_trkP; gFile->GetObject("hm2rec3OS_ttbb_trkP", hm2rec3OS_ttbb_trkP);
TH1F *hm2rec3OS_ttbb_trkM; gFile->GetObject("hm2rec3OS_ttbb_trkM", hm2rec3OS_ttbb_trkM);
TH1F *hm2rec3OS_diag_pypxP; gFile->GetObject("hm2rec3OS_diag_pypxP", hm2rec3OS_diag_pypxP);
TH1F *hm2rec3OS_diag_pypxM; gFile->GetObject("hm2rec3OS_diag_pypxM", hm2rec3OS_diag_pypxM);
TH1F *hm2rec3OS_ttbb_pypxP; gFile->GetObject("hm2rec3OS_ttbb_pypxP", hm2rec3OS_ttbb_pypxP);
TH1F *hm2rec3OS_ttbb_pypxM; gFile->GetObject("hm2rec3OS_ttbb_pypxM", hm2rec3OS_ttbb_pypxM);
TH1F *hm2rec4OS; gFile->GetObject("hm2rec4OS", hm2rec4OS);
TH1F *hm2rec4SS; gFile->GetObject("hm2rec4SS", hm2rec4SS);
TH1F *hm2rec4OS_diag; gFile->GetObject("hm2rec4OS_diag", hm2rec4OS_diag);
TH1F *hm2rec4SS_diag; gFile->GetObject("hm2rec4SS_diag", hm2rec4SS_diag);
TH1F *hm2rec4OS_ttbb; gFile->GetObject("hm2rec4OS_ttbb", hm2rec4OS_ttbb);
TH1F *hm2rec4SS_ttbb; gFile->GetObject("hm2rec4SS_ttbb", hm2rec4SS_ttbb);
TH1F *hm2rec4OS_diag_trkP; gFile->GetObject("hm2rec4OS_diag_trkP", hm2rec4OS_diag_trkP);
TH1F *hm2rec4OS_diag_trkM; gFile->GetObject("hm2rec4OS_diag_trkM", hm2rec4OS_diag_trkM);
TH1F *hm2rec4OS_ttbb_trkP; gFile->GetObject("hm2rec4OS_ttbb_trkP", hm2rec4OS_ttbb_trkP);
TH1F *hm2rec4OS_ttbb_trkM; gFile->GetObject("hm2rec4OS_ttbb_trkM", hm2rec4OS_ttbb_trkM);
TH1F *hm2rec4OS_diag_pypxP; gFile->GetObject("hm2rec4OS_diag_pypxP", hm2rec4OS_diag_pypxP);
TH1F *hm2rec4OS_diag_pypxM; gFile->GetObject("hm2rec4OS_diag_pypxM", hm2rec4OS_diag_pypxM);
TH1F *hm2rec4OS_ttbb_pypxP; gFile->GetObject("hm2rec4OS_ttbb_pypxP", hm2rec4OS_ttbb_pypxP);
TH1F *hm2rec4OS_ttbb_pypxM; gFile->GetObject("hm2rec4OS_ttbb_pypxM", hm2rec4OS_ttbb_pypxM);
TH1F *hm2rec5OS; gFile->GetObject("hm2rec5OS", hm2rec5OS);
TH1F *hm2rec5SS; gFile->GetObject("hm2rec5SS", hm2rec5SS);
TH1F *hm2rec5OS_diag; gFile->GetObject("hm2rec5OS_diag", hm2rec5OS_diag);
TH1F *hm2rec5SS_diag; gFile->GetObject("hm2rec5SS_diag", hm2rec5SS_diag);
TH1F *hm2rec5OS_ttbb; gFile->GetObject("hm2rec5OS_ttbb", hm2rec5OS_ttbb);
TH1F *hm2rec5SS_ttbb; gFile->GetObject("hm2rec5SS_ttbb", hm2rec5SS_ttbb);
TH1F *hm2rec5OS_diag_trkP; gFile->GetObject("hm2rec5OS_diag_trkP", hm2rec5OS_diag_trkP);
TH1F *hm2rec5OS_diag_trkM; gFile->GetObject("hm2rec5OS_diag_trkM", hm2rec5OS_diag_trkM);
TH1F *hm2rec5OS_ttbb_trkP; gFile->GetObject("hm2rec5OS_ttbb_trkP", hm2rec5OS_ttbb_trkP);
TH1F *hm2rec5OS_ttbb_trkM; gFile->GetObject("hm2rec5OS_ttbb_trkM", hm2rec5OS_ttbb_trkM);
TH1F *hm2rec5OS_diag_pypxP; gFile->GetObject("hm2rec5OS_diag_pypxP", hm2rec5OS_diag_pypxP);
TH1F *hm2rec5OS_diag_pypxM; gFile->GetObject("hm2rec5OS_diag_pypxM", hm2rec5OS_diag_pypxM);
TH1F *hm2rec5OS_ttbb_pypxP; gFile->GetObject("hm2rec5OS_ttbb_pypxP", hm2rec5OS_ttbb_pypxP);
TH1F *hm2rec5OS_ttbb_pypxM; gFile->GetObject("hm2rec5OS_ttbb_pypxM", hm2rec5OS_ttbb_pypxM);
TH1F *hm2rec6OS; gFile->GetObject("hm2rec6OS", hm2rec6OS);
TH1F *hm2rec6SS; gFile->GetObject("hm2rec6SS", hm2rec6SS);
TH1F *hm2rec6OS_diag; gFile->GetObject("hm2rec6OS_diag", hm2rec6OS_diag);
TH1F *hm2rec6SS_diag; gFile->GetObject("hm2rec6SS_diag", hm2rec6SS_diag);
TH1F *hm2rec6OS_ttbb; gFile->GetObject("hm2rec6OS_ttbb", hm2rec6OS_ttbb);
TH1F *hm2rec6SS_ttbb; gFile->GetObject("hm2rec6SS_ttbb", hm2rec6SS_ttbb);
TH1F *hm2rec6OS_diag_trkP; gFile->GetObject("hm2rec6OS_diag_trkP", hm2rec6OS_diag_trkP);
TH1F *hm2rec6OS_diag_trkM; gFile->GetObject("hm2rec6OS_diag_trkM", hm2rec6OS_diag_trkM);
TH1F *hm2rec6OS_ttbb_trkP; gFile->GetObject("hm2rec6OS_ttbb_trkP", hm2rec6OS_ttbb_trkP);
TH1F *hm2rec6OS_ttbb_trkM; gFile->GetObject("hm2rec6OS_ttbb_trkM", hm2rec6OS_ttbb_trkM);
TH1F *hm2recHFvetoOS; gFile->GetObject("hm2recHFvetoOS", hm2recHFvetoOS);
TH1F *hm2recHFvetoSS; gFile->GetObject("hm2recHFvetoSS", hm2recHFvetoSS);
TH1F *hm2rec45OS; gFile->GetObject("hm2rec45OS", hm2rec45OS);
TH1F *hm2rec45SS; gFile->GetObject("hm2rec45SS", hm2rec45SS);
TH1F *hm2rec4515OS; gFile->GetObject("hm2rec4515OS", hm2rec4515OS);
TH1F *hm2rec4515SS; gFile->GetObject("hm2rec4515SS", hm2rec4515SS);
TH1F *hm2rec9919; gFile->GetObject("hm2rec9919", hm2rec9919);
TH1F *hm2rec9922; gFile->GetObject("hm2rec9922", hm2rec9922);
TH1F *hm2rec9971; gFile->GetObject("hm2rec9971", hm2rec9971);
TH1F *hm2rec9978; gFile->GetObject("hm2rec9978", hm2rec9978);
TH1F *hnclusters; gFile->GetObject("hnclusters", hnclusters);
TH1F *hnclusters2; gFile->GetObject("hnclusters2", hnclusters2);
TH1F *hnclustersOSdiag; gFile->GetObject("hnclustersOSdiag", hnclustersOSdiag);
TH1F *hnclusters2OSdiag; gFile->GetObject("hnclusters2OSdiag", hnclusters2OSdiag);
TH1F *halgo; gFile->GetObject("halgo", halgo);
TH1F *hnhits; gFile->GetObject("hnhits", hnhits);
TH1F *hchi2; gFile->GetObject("hchi2", hchi2);
TH1F *hdz; gFile->GetObject("hdz", hdz);
TH1F *hd0; gFile->GetObject("hd0", hd0);
TH1F *halgov; gFile->GetObject("halgov", halgov);
TH1F *hnhitsv; gFile->GetObject("hnhitsv", hnhitsv);
TH1F *hchi2v; gFile->GetObject("hchi2v", hchi2v);
TH1F *hdzv; gFile->GetObject("hdzv", hdzv);
TH1F *hd0v; gFile->GetObject("hd0v", hd0v);
TH1F *hchi2fin; gFile->GetObject("hchi2fin", hchi2fin);
TH1F *hdzfin; gFile->GetObject("hdzfin", hdzfin);
TH1F *hd0fin; gFile->GetObject("hd0fin", hd0fin);
TH1F *hdeltaR; gFile->GetObject("hdeltaR", hdeltaR);
TH1F *hdeltaR2; gFile->GetObject("hdeltaR2", hdeltaR2);
//Histogram variable
TH2F *h2dimdpyAll; gFile->GetObject("h2dimdpyAll", h2dimdpyAll);
TH2F *h2dimdpy; gFile->GetObject("h2dimdpy", h2dimdpy);
TH2F *h2dimdpy_diag; gFile->GetObject("h2dimdpy_diag", h2dimdpy_diag);
TH2F *h2dimdpy_ttbb; gFile->GetObject("h2dimdpy_ttbb", h2dimdpy_ttbb);
//Histogram variable
TH1F *hdpyAll; gFile->GetObject("hdpyAll", hdpyAll);
TH1F *hdpy; gFile->GetObject("hdpy", hdpy);
TH1F *hdpy_diag; gFile->GetObject("hdpy_diag", hdpy_diag);
TH1F *hdpy_ttbb; gFile->GetObject("hdpy_ttbb", hdpy_ttbb);
//Histogram variable
TH2F *h2dimdpxAll; gFile->GetObject("h2dimdpxAll", h2dimdpxAll);
TH2F *h2dimdpx; gFile->GetObject("h2dimdpx", h2dimdpx);
TH2F *h2dimdpx_diag; gFile->GetObject("h2dimdpx_diag", h2dimdpx_diag);
TH2F *h2dimdpx_ttbb; gFile->GetObject("h2dimdpx_ttbb", h2dimdpx_ttbb);
//Histogram variable
TH1F *hdpxAll; gFile->GetObject("hdpxAll", hdpxAll);
TH1F *hdpx; gFile->GetObject("hdpx", hdpx);
TH1F *hdpx_diag; gFile->GetObject("hdpx_diag", hdpx_diag);
TH1F *hdpx_ttbb; gFile->GetObject("hdpx_ttbb", hdpx_ttbb);
//Histogram variable
TH2F *h2dimxVtxRL; gFile->GetObject("h2dimxVtxRL", h2dimxVtxRL);
TH2F *h2dimxVtxcmsR; gFile->GetObject("h2dimxVtxcmsR", h2dimxVtxcmsR);
TH2F *h2dimxVtxcmsL; gFile->GetObject("h2dimxVtxcmsL", h2dimxVtxcmsL);
TH2F *h2dimxVtxcmsRL; gFile->GetObject("h2dimxVtxcmsRL", h2dimxVtxcmsRL);
TH2F *h2dimxVtxcmsR2; gFile->GetObject("h2dimxVtxcmsR2", h2dimxVtxcmsR2);
TH2F *h2dimxVtxcmsL2; gFile->GetObject("h2dimxVtxcmsL2", h2dimxVtxcmsL2);
TH2F *h2dimxVtxcmsRL2; gFile->GetObject("h2dimxVtxcmsRL2", h2dimxVtxcmsRL2);
TH2F *h2dimxVtx_zVtx_CT; gFile->GetObject("h2dimxVtx_zVtx_CT", h2dimxVtx_zVtx_CT);
TH2F *h2dimxVtx_zVtx_C; gFile->GetObject("h2dimxVtx_zVtx_C", h2dimxVtx_zVtx_C);
TH2F *h2dimxVtx_zVtx_T; gFile->GetObject("h2dimxVtx_zVtx_T", h2dimxVtx_zVtx_T);
//Histogram variable
TH1F *hxVtxRL; gFile->GetObject("hxVtxRL", hxVtxRL);
TH1F *hxVtxcmsR; gFile->GetObject("hxVtxcmsR", hxVtxcmsR);
TH1F *hxVtxcmsL; gFile->GetObject("hxVtxcmsL", hxVtxcmsL);
TH1F *hxVtxcmsRL; gFile->GetObject("hxVtxcmsRL", hxVtxcmsRL);
TH1F *hxVtxRL_diag; gFile->GetObject("hxVtxRL_diag", hxVtxRL_diag);
TH1F *hxVtxcmsR_diag; gFile->GetObject("hxVtxcmsR_diag", hxVtxcmsR_diag);
TH1F *hxVtxcmsL_diag; gFile->GetObject("hxVtxcmsL_diag", hxVtxcmsL_diag);
TH1F *hxVtxcmsRL_diag; gFile->GetObject("hxVtxcmsRL_diag", hxVtxcmsRL_diag);
TH1F *hxVtxRL_ttbb; gFile->GetObject("hxVtxRL_ttbb", hxVtxRL_ttbb);
TH1F *hxVtxcmsR_ttbb; gFile->GetObject("hxVtxcmsR_ttbb", hxVtxcmsR_ttbb);
TH1F *hxVtxcmsL_ttbb; gFile->GetObject("hxVtxcmsL_ttbb", hxVtxcmsL_ttbb);
TH1F *hxVtxcmsRL_ttbb; gFile->GetObject("hxVtxcmsRL_ttbb", hxVtxcmsRL_ttbb);
//Histogram variable
TH2F *hdedx; gFile->GetObject("hdedx", hdedx);
//
TH2F *hlndedx; gFile->GetObject("hlndedx", hlndedx);
TH2F *hl10dedx; gFile->GetObject("hl10dedx", hl10dedx);
//
TH2F *phi_proton_right_t; gFile->GetObject("phi_proton_right_t", phi_proton_right_t);
TH2F *phi_proton_left_t; gFile->GetObject("phi_proton_left_t", phi_proton_left_t);
//
TH2F *phi_proton_right_t_diag; gFile->GetObject("phi_proton_right_t_diag", phi_proton_right_t_diag);
TH2F *phi_proton_left_t_diag; gFile->GetObject("phi_proton_left_t_diag", phi_proton_left_t_diag);
//
TH2F *phi_proton_right_t_ttbb; gFile->GetObject("phi_proton_right_t_ttbb", phi_proton_right_t_ttbb);
TH2F *phi_proton_left_t_ttbb; gFile->GetObject("phi_proton_left_t_ttbb", phi_proton_left_t_ttbb);
//
TH2F *phi_proton_right_t_tt; gFile->GetObject("phi_proton_right_t_tt", phi_proton_right_t_tt);
TH2F *phi_proton_left_t_tt; gFile->GetObject("phi_proton_left_t_tt", phi_proton_left_t_tt);
//
TH2F *phi_proton_right_t_bb; gFile->GetObject("phi_proton_right_t_bb", phi_proton_right_t_bb);
TH2F *phi_proton_left_t_bb; gFile->GetObject("phi_proton_left_t_bb", phi_proton_left_t_bb);
//
TH1F *hrapy; gFile->GetObject("hrapy", hrapy);
TH1F *hrapy2; gFile->GetObject("hrapy2", hrapy2);
//
  //TDirectory *dir = _file0->GetDirectory("analyzeHiMassTau");
  //if (dir) {
  //   dir->cd();
  // }
  // else {
  //  fprintf(stderr,"Missing directory analyzeHiMassTau\n");
  //}
  //
  gSystem->mkdir("plotall", kTRUE);
  // ...
  //gSystem->cd("../");
  //gSystem->Exec("mkdir plotall");
  cout << "                                 " << endl;
  cout << "    ...the subdirectory \"plotall\" has been created. " << endl;
  cout << "                                 " << endl;
  //
  //gSystem->Exec("cd ./plotall");
  gSystem->cd("plotall");
  
//Histogram variable-----------------------------Plot title
  plot1d(hnconf,"n conf","1.0","lin","lin");// Number of configurations (TB or BT or TT or BB)
  plot1d(rp_x_020,"X (cm)","0.1cm","lin","lin");// x RP
  plot1d(rp_x_021,"X (cm)","0.1cm","lin","lin");// x RP
  plot1d(rp_x_024,"X (cm)","0.1cm","lin","lin");// x RP
  plot1d(rp_x_025,"X (cm)","0.1cm","lin","lin");// x RP
  plot1d(rp_x_120,"X (cm)","0.1cm","lin","lin");// x RP
  plot1d(rp_x_121,"X (cm)","0.1cm","lin","lin");// x RP
  plot1d(rp_x_124,"X (cm)","0.1cm","lin","lin");// x RP
  plot1d(rp_x_125,"X (cm)","0.1cm","lin","lin");// x RP
  plot1d(rp_y_020,"Y (cm)","0.2cm","lin","lin");// y RP
  plot1d(rp_y_021,"Y (cm)","0.2cm","lin","lin");// y RP
  plot1d(rp_y_024,"Y (cm)","0.2cm","lin","lin");// y RP
  plot1d(rp_y_025,"Y (cm)","0.2cm","lin","lin");// y RP
  plot1d(rp_y_120,"Y (cm)","0.2cm","lin","lin");// y RP
  plot1d(rp_y_121,"Y (cm)","0.2cm","lin","lin");// y RP
  plot1d(rp_y_124,"Y (cm)","0.2cm","lin","lin");// y RP
  plot1d(rp_y_125,"Y (cm)","0.2cm","lin","lin");// y RP
  plot1d(rp2_x_020,"X (cm)","0.1cm","lin","lin");// x RP
  plot1d(rp2_x_021,"X (cm)","0.1cm","lin","lin");// x RP
  plot1d(rp2_x_024,"X (cm)","0.1cm","lin","lin");// x RP
  plot1d(rp2_x_025,"X (cm)","0.1cm","lin","lin");// x RP
  plot1d(rp2_x_120,"X (cm)","0.1cm","lin","lin");// x RP
  plot1d(rp2_x_121,"X (cm)","0.1cm","lin","lin");// x RP
  plot1d(rp2_x_124,"X (cm)","0.1cm","lin","lin");// x RP
  plot1d(rp2_x_125,"X (cm)","0.1cm","lin","lin");// x RP
  plot1d(rp2_y_020,"Y (cm)","0.2cm","lin","lin");// y RP
  plot1d(rp2_y_021,"Y (cm)","0.2cm","lin","lin");// y RP
  plot1d(rp2_y_024,"Y (cm)","0.2cm","lin","lin");// y RP
  plot1d(rp2_y_025,"Y (cm)","0.2cm","lin","lin");// y RP
  plot1d(rp2_y_120,"Y (cm)","0.2cm","lin","lin");// y RP
  plot1d(rp2_y_121,"Y (cm)","0.2cm","lin","lin");// y RP
  plot1d(rp2_y_124,"Y (cm)","0.2cm","lin","lin");// y RP
  plot1d(rp2_y_125,"Y (cm)","0.2cm","lin","lin");// y RP
  plot1d(thyEla,"#theta_{y}^{LR} (mrad)","2.0E-07mrad","lin","lin");// thyL+thyR
  plot1d(thxEla,"#theta_{x}^{LR} (mrad)","2.0E-07mrad","lin","lin");// thxL+thxR
  plot1d(thyEla_diag,"#theta_{y}^{LR} (mrad)","2.0E-07mrad","lin","lin");// thyL+thyR diagonals
  plot1d(thxEla_diag,"#theta_{x}^{LR} (mrad)","2.0E-07mrad","lin","lin");// thxL+thxR diagonals
  plot1d(thyEla_ttbb,"#theta_{y}^{LR} (mrad)","2.0E-07mrad","lin","lin");// thyL+thyR TT/BB
  plot1d(thxEla_ttbb,"#theta_{x}^{LR} (mrad)","2.0E-07mrad","lin","lin");// thxL+thxR TT/BB
  plot1d(proton_right_xi,"#xi^{R}","0.005","lin","lin");// #xi
  plot1d(proton_left_xi,"#xi^{L}","0.005","lin","lin");// #xi
  plot1d(proton_right_logXi,"log(#xi)^{R}","0.025","lin","lin");// log(#xi)
  plot1d(proton_left_logXi,"log(#xi)^{L}","0.025","lin","lin");// log(#xi)
  plot1d(proton_right_t,"|t| (GeV^{2})","0.005GeV^{2}","log","lin");// -t
  plot1d(proton_left_t,"|t| (GeV^{2})","0.005GeV^{2}","log","lin");// -t
  plot1d(proton_right_t_diag,"|t| (GeV^{2})","0.005GeV^{2}","log","lin");// -t diagonal
  plot1d(proton_left_t_diag,"|t| (GeV^{2})","0.005GeV^{2}","log","lin");// -t diagonal
  plot1d(proton_right_t_ttbb,"|t| (GeV^{2})","0.005GeV^{2}","log","lin");// -t
  plot1d(proton_left_t_ttbb,"|t| (GeV^{2})","0.005GeV^{2}","log","lin");// -t TT/BB
  plot1d(eHF,"energy HF tower (GeV)","0.2GeV","lin","lin");// energy HF tower (GeV)
  plot1d(nHF,"n HF tower","1.0","lin","lin");// n HF tower (eHF>5 GeV)
  plot1d(totem_py,"p_{Y} (GeV/c)","0.02GeV/c","lin","lin");// p_{Y} TOTEM
  plot1d(totem_px,"p_{X} (GeV/c)","0.02GeV/c","lin","lin");// p_{X} TOTEM
  plot1d(totem_pyy,"p_{Y} (GeV/c)","0.002GeV/c","lin","lin");// p_{Y} TOTEM
  plot1d(totem_pxx,"p_{X} (GeV/c)","0.002GeV/c","lin","lin");// p_{X} TOTEM
  plot1d(proton_dx0,"#Deltap_{X0} (GeV/c)","2.0E-06GeV/c","lin","lin");// xVtx_{56}-xVtx_{45}
  plot1d(hLS,"LS","1.0","lin","lin");// LS
  plot1d(htopo,"topology","1.0","lin","lin");// 1=TB 2=BT 3=TT 4=BB topology
  plot1d(hthyEla2_diag,"#theta_{y}^{LR} (mrad)","2.0E-07mrad","lin","lin");// thyL+thyR dig
  plot1d(hthxEla2_diag,"#theta_{x}^{LR} (mrad)","2.0E-07mrad","lin","lin");// thxL+thxR dig
  plot1d(hthyEla2_ttbb,"#theta_{y}^{LR} (mrad)","2.0E-07mrad","lin","lin");// thyL+thyR TTBB
  plot1d(hthxEla2_ttbb,"#theta_{x}^{LR} (mrad)","2.0E-07mrad","lin","lin");// thxL+thxR TTBB
//Histogram variable Plot title
  plot2d(rp_yx_020,"Y (cm)","X (cm)","lin","lin");// y vs x RP
  plot2d(rp_yx_021,"Y (cm)","X (cm)","lin","lin");// y vs x RP
  plot2d(rp_yx_024,"Y (cm)","X (cm)","lin","lin");// y vs x RP
  plot2d(rp_yx_025,"Y (cm)","X (cm)","lin","lin");// y vs x RP
  plot2d(rp_yx_120,"Y (cm)","X (cm)","lin","lin");// y vs x RP
  plot2d(rp_yx_121,"Y (cm)","X (cm)","lin","lin");// y vs x RP
  plot2d(rp_yx_124,"Y (cm)","X (cm)","lin","lin");// y vs x RP
  plot2d(rp_yx_125,"Y (cm)","X (cm)","lin","lin");// y vs x RP
  plot2d(rp2_yx_020,"Y (cm)","X (cm)","lin","lin");// y vs x RP
  plot2d(rp2_yx_021,"Y (cm)","X (cm)","lin","lin");// y vs x RP
  plot2d(rp2_yx_024,"Y (cm)","X (cm)","lin","lin");// y vs x RP
  plot2d(rp2_yx_025,"Y (cm)","X (cm)","lin","lin");// y vs x RP
  plot2d(rp2_yx_120,"Y (cm)","X (cm)","lin","lin");// y vs x RP
  plot2d(rp2_yx_121,"Y (cm)","X (cm)","lin","lin");// y vs x RP
  plot2d(rp2_yx_124,"Y (cm)","X (cm)","lin","lin");// y vs x RP
  plot2d(rp2_yx_125,"Y (cm)","X (cm)","lin","lin");// y vs x RP
  plot2d(proton_x0_RvsL,"X (cm)","X (cm)","lin","lin");// xVtx_{56} vs xVtx_{45}
//Histogram variable Plot title
  plot1d(hlooper,"Looper","1.0","lin","lin");// isLooper
  plot1d(hpt,"p_{T} (GeV/c)","0.05GeV/c","lin","lin");// p_{T}
  plot1d(heta,"#eta","0.1","lin","lin");// #eta
  plot1d(hphi,"#varphi (rad)","0.1rad","lin","lin");// #varphi
  plot1d(hptP,"p_{T} (GeV/c)","0.02GeV/c","lin","lin");// p_{T} #pi+
  plot1d(hetaP,"#eta","0.1","lin","lin");// #eta #pi+
  plot1d(hphiP,"#varphi (rad)","0.1rad","lin","lin");// #varphi #pi+
  plot1d(hptM,"p_{T} (GeV/c)","0.02GeV/c","lin","lin");// p_{T} #pi-
  plot1d(hetaM,"#eta","0.1","lin","lin");// #eta #pi-
  plot1d(hphiM,"#varphi (rad)","0.1rad","lin","lin");// #varphi #pi-
  plot1d(hptRes,"p_{T} (GeV/c)","0.02GeV/c","lin","lin");// p_{T} #pi#pi
  plot1d(hetaRes,"#eta","0.1","lin","lin");// #eta #pi#pi
  plot1d(hphiRes,"#varphi (rad)","0.1rad","lin","lin");// #varphi #pi#pi
  plot1d(hthyEla_diag,"#theta_{y}^{RL,Ela} (mrad)","2.0E-07mrad","lin","lin");// thyL+thyR dig
  plot1d(hthxEla_diag,"#theta_{x}^{RL,Ela} (mrad)","2.0E-07mrad","lin","lin");// thxL+thxR dig
  plot1d(hthyEla_ttbb,"#theta_{y}^{RL,Ela} (mrad)","2.0E-07mrad","lin","lin");// thyL+thyR TTBB
  plot1d(hthxEla_ttbb,"#theta_{x}^{RL,Ela} (mrad)","2.0E-07mrad","lin","lin");// thxL+thxR TTBB
  plot1d(hntrk0,"N trk0","1.0","lin","lin");// Ntrk
  plot1d(hntrk,"N trk","1.0","lin","lin");// Ntrk for nPixelHits>0
  plot1d(hntrkvtx,"N trk vtx","1.0","lin","lin");// Ntrkvtx
  plot1d(hntrkntrkvtx2,"N trk","1.0","lin","lin");// Ntrk for Ntrkvtx==2
  plot1d(hntrk2ntrkvtx,"N trk vtx","1.0","lin","lin");// Ntrkvtx for Ntrk==2
//Histogram variable Plot title
  plot2d(hntrkntrkvtx,"N trk","N trk vtx","lin","lin");// Ntrk vs Ntrkvtx
//Histogram variable Plot title
  plot1d(hvtx,"vtx","1.0","lin","lin");// vtx.isFake()
  plot1d(hvtx2,"vtx2","1.0","lin","lin");// vtx.isFake() 2 tracks
  plot1d(hvtx3,"vtx3","1.0","lin","lin");// vtx.isFake() 2 tracks both |#eta|<2.5 and OS
  plot1d(hnvtx,"N vtx","1.0","lin","lin");// Nvtx
  plot1d(hvtxx,"X vtx (cm)","0.002cm","lin","lin");// X vtx
  plot1d(hvtxy,"Y vtx (cm)","0.002cm","lin","lin");// Y vtx
  plot1d(hvtxz,"Z vtx (cm)","0.2cm","lin","lin");// Z vtx
  plot1d(hvtxchi2,"#chi^{2} vtx","1.0","lin","lin");// chi2 vtx
  plot1d(hvtxchi2fin,"#chi^{2} vtx fin","1.0","lin","lin");// chi2 vtx fin
  plot1d(heHF,"HF tower energy","0.2GeV","lin","lin");// HF tower energy
  plot1d(hnHF,"n HF towers","1.0","lin","lin");// n HF towers (E>5 GeV)
  plot1d(hxiL,"#xi L ","0.002","lin","lin");// #xiL 
  plot1d(hxiR,"#xi R ","0.002","lin","lin");// #xiR 
  plot1d(hxiL2,"#xi L2","0.002","lin","lin");// #xiL 
  plot1d(hxiR2,"#xi R2","0.002","lin","lin");// #xiR 
  plot1d(hm,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} 
  plot1d(hmxicut,"M_{#pi#pi} (GeV/c^{2}) #xi cut","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi}
  plot1d(hm2rec,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} 
  plot1d(hm2recbis,"M_{#pi#pi} (GeV/c^{2}) 2rec bis","0.01GeV/c^{2}","lin","lin");// M_{#pi#pi}
  plot1d(hm2recPP,"M_{#pi#pi} (GeV/c^{2}) 2rec PP","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} 
  plot1d(hm2recKK,"M_{KK} (GeV/c^{2}) 2rec KK","0.02GeV/c^{2}","lin","lin");// M_{KK} 
  plot1d(hm2recMM,"M_{#mu#mu} (GeV/c^{2}) 2rec MM","0.02GeV/c^{2}","lin","lin");// M_{#mu#mu} 
  plot1d(hm2recEE,"M_{ee} (GeV/c^{2}) 2rec EE","0.02GeV/c^{2}","lin","lin");// M_{ee} 
  plot1d(hm2recOS,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} OS
  plot1d(hm2recSS,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} SS
  plot1d(hm2recOS_diag,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TB/BT OS
  plot1d(hm2recSS_diag,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TB/BT SS
  plot1d(hm2recOS_ttbb,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TT/BB OS
  plot1d(hm2recSS_ttbb,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TT/BB SS
  plot1d(hm2rec2OS,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} OS
  plot1d(hm2rec2SS,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} SS
  plot1d(hm2rec2OS_diag,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TB/BT OS
  plot1d(hm2rec2SS_diag,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TB/BT SS
  plot1d(hm2rec2OS_ttbb,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TT/BB OS
  plot1d(hm2rec2SS_ttbb,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TT/BB SS
  plot1d(hm2rec2OS_diag_trkP,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TB/BT OS py_{#pi1}py_{#pi2}>0
  plot1d(hm2rec2OS_diag_trkM,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TB/BT OS py_{#pi1}py_{#pi2}<0
  plot1d(hm2rec2OS_ttbb_trkP,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TT/BB OS py_{#pi1}py_{#pi2}>0
  plot1d(hm2rec2OS_ttbb_trkM,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TT/BB OS py_{#pi1}py_{#pi2}<0
  plot1d(hm2rec2OS_diag_pypxP,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TB/BT OS |py/px|_{#pi#pi} > 1
  plot1d(hm2rec2OS_diag_pypxM,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TB/BT OS |py/px|_{#pi#pi} < 1
  plot1d(hm2rec2OS_ttbb_pypxP,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TT/BB OS |py/px|_{#pi#pi} > 1
  plot1d(hm2rec2OS_ttbb_pypxM,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TT/BB OS |py/px|_{#pi#pi} < 1
  plot1d(hm2rec3OS,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} OS
  plot1d(hm2rec3SS,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} SS
  plot1d(hm2rec3OS_diag,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TB/BT OS
  plot1d(hm2rec3SS_diag,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TB/BT SS
  plot1d(hm2rec3OS_ttbb,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TT/BB OS
  plot1d(hm2rec3SS_ttbb,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TT/BB SS
  plot1d(hm2rec3OS_diag_trkP,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TB/BT OS py_{#pi1}py_{#pi2}>0
  plot1d(hm2rec3OS_diag_trkM,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TB/BT OS py_{#pi1}py_{#pi2}<0
  plot1d(hm2rec3OS_ttbb_trkP,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TT/BB OS py_{#pi1}py_{#pi2}>0
  plot1d(hm2rec3OS_ttbb_trkM,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TT/BB OS py_{#pi1}py_{#pi2}<0
  plot1d(hm2rec3OS_diag_pypxP,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TB/BT OS |py/px|_{#pi#pi} > 1
  plot1d(hm2rec3OS_diag_pypxM,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TB/BT OS |py/px|_{#pi#pi} < 1
  plot1d(hm2rec3OS_ttbb_pypxP,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TT/BB OS |py/px|_{#pi#pi} > 1
  plot1d(hm2rec3OS_ttbb_pypxM,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TT/BB OS |py/px|_{#pi#pi} < 1
  plot1d(hm2rec4OS,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} OS
  plot1d(hm2rec4SS,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} SS
  plot1d(hm2rec4OS_diag,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TB/BT OS
  plot1d(hm2rec4SS_diag,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TB/BT SS
  plot1d(hm2rec4OS_ttbb,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TT/BB OS
  plot1d(hm2rec4SS_ttbb,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TT/BB SS
  plot1d(hm2rec4OS_diag_trkP,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi}TB/BT OS py_{#pi1}py_{#pi2}>0
  plot1d(hm2rec4OS_diag_trkM,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi}TB/BT OS py_{#pi1}py_{#pi2}<0
  plot1d(hm2rec4OS_ttbb_trkP,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi}TT/BB OS py_{#pi1}py_{#pi2}>0
  plot1d(hm2rec4OS_ttbb_trkM,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi}TT/BB OS py_{#pi1}py_{#pi2}<0
  plot1d(hm2rec4OS_diag_pypxP,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TB/BT OS |py/px|_{#pi#pi} > 1
  plot1d(hm2rec4OS_diag_pypxM,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TB/BT OS |py/px|_{#pi#pi} < 1
  plot1d(hm2rec4OS_ttbb_pypxP,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TT/BB OS |py/px|_{#pi#pi} > 1
  plot1d(hm2rec4OS_ttbb_pypxM,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TT/BB OS |py/px|_{#pi#pi} < 1
  plot1d(hm2rec5OS,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} OS
  plot1d(hm2rec5SS,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} SS
  plot1d(hm2rec5OS_diag,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TB/BT OS
  plot1d(hm2rec5SS_diag,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TB/BT SS
  plot1d(hm2rec5OS_ttbb,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TT/BB OS
  plot1d(hm2rec5SS_ttbb,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TT/BB SS
  plot1d(hm2rec5OS_diag_trkP,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TB/BT OS py_{#pi1}py_{#pi2}>0
  plot1d(hm2rec5OS_diag_trkM,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TB/BT OS py_{#pi1}py_{#pi2}<0
  plot1d(hm2rec5OS_ttbb_trkP,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TB/BT OS py_{#pi1}py_{#pi2}>0
  plot1d(hm2rec5OS_ttbb_trkM,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TB/BT OS py_{#pi1}py_{#pi2}<0
  plot1d(hm2rec5OS_diag_pypxP,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TB/BT OS |py/px|_{#pi#pi} > 1
  plot1d(hm2rec5OS_diag_pypxM,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TB/BT OS |py/px|_{#pi#pi} < 1
  plot1d(hm2rec5OS_ttbb_pypxP,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TT/BB OS |py/px|_{#pi#pi} > 1
  plot1d(hm2rec5OS_ttbb_pypxM,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TT/BB OS |py/px|_{#pi#pi} < 1
  plot1d(hm2rec6OS,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} OS
  plot1d(hm2rec6SS,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} SS
  plot1d(hm2rec6OS_diag,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TB/BT OS
  plot1d(hm2rec6SS_diag,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TB/BT SS
  plot1d(hm2rec6OS_ttbb,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TT/BB OS
  plot1d(hm2rec6SS_ttbb,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TT/BB SS
  plot1d(hm2rec6OS_diag_trkP,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TB/BT OS py_{#pi1}py_{#pi2}>0
  plot1d(hm2rec6OS_diag_trkM,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TB/BT OS py_{#pi1}py_{#pi2}<0
  plot1d(hm2rec6OS_ttbb_trkP,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TT/BB OS py_{#pi1}py_{#pi2}>0
  plot1d(hm2rec6OS_ttbb_trkM,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} TT/BB OS py_{#pi1}py_{#pi2}<0
  plot1d(hm2recHFvetoOS,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} HFv OS
  plot1d(hm2recHFvetoSS,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} HFv SS
  plot1d(hm2rec45OS,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} OS
  plot1d(hm2rec45SS,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} SS
  plot1d(hm2rec4515OS,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} OS
  plot1d(hm2rec4515SS,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} SS
  plot1d(hm2rec9919,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} 9919
  plot1d(hm2rec9922,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} 9919 9922
  plot1d(hm2rec9971,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} 9971
  plot1d(hm2rec9978,"M_{#pi#pi} (GeV/c^{2})","0.02GeV/c^{2}","lin","lin");// M_{#pi#pi} 9978
  plot1d(hnclusters,"n Pixel Clusters","1.0","lin","lin");// nPixelClusters
  plot1d(hnclusters2,"n Strip Clusters","1.0","lin","lin");// nStripClusters
  plot1d(hnclustersOSdiag,"n Pixel Clusters","1.0","lin","lin");// nPixelClusters
  plot1d(hnclusters2OSdiag,"n Strip Clusters","1.0","lin","lin");// nStripClusters
  plot1d(halgo,"Algo","1.0","lin","lin");// Algo
  plot1d(hnhits,"nhits pix+strip","1.0","lin","lin");// nhits pix+strip
  plot1d(hchi2,"normalized #chi^{2}","1.0","lin","lin");// normalized #chi^{2}
  plot1d(hdz,"dz","0.1mm","lin","lin");// dz
  plot1d(hd0,"d0","0.1mm","lin","lin");// d0
  plot1d(halgov,"Algo","1.0","lin","lin");// Algo
  plot1d(hnhitsv,"nhits pixel","1.0","lin","lin");// nhits pixel
  plot1d(hchi2v,"normalized #chi^{2} vtx-fitted","1.0","lin","lin");// normalized #chi^{2} vtx-fitted
  plot1d(hdzv,"dz vtx-fitted","0.2mm","lin","lin");// dz vtx-fitted
  plot1d(hd0v,"d0 vtx-fitted","0.02mm","lin","lin");// d0 vtx-fitted
  plot1d(hchi2fin,"normalized #chi^{2} vtx-fitted","1.0","lin","lin");// normalized #chi^{2} vtx-fitted
  plot1d(hdzfin,"dz vtx-fitted","0.2mm","lin","lin");// dz vtx-fitted
  plot1d(hd0fin,"d0 vtx-fitted","0.02mm","lin","lin");// d0 vtx-fitted
  plot1d(hdeltaR,"#DeltaR trk-trk","0.05cm","lin","lin");// #Delta R trk-trk
  plot1d(hdeltaR2,"#DeltaR trk-trk","0.05cm","lin","lin");// #Delta R trk-trk
//Histogram variable Plot title
  plot2d(h2dimdpyAll,"p_{y}^{TOTEM} (GeV/c)","p_{y}^{CMS} (GeV/c)","lin","lin");// p_{y}^{TOTEM} vs p_{y}^{CMS} All
  plot2d(h2dimdpy,"p_{y}^{TOTEM} (GeV/c)","p_{y}^{CMS} (GeV/c)","lin","lin");// p_{y}^{TOTEM} vs p_{y}^{CMS}
  plot2d(h2dimdpy_diag,"p_{y}^{TOTEM} (GeV/c)","p_{y}^{CMS} (GeV/c)","lin","lin");// p_{y}^{TOTEM} vs p_{y}^{CMS} diag
  plot2d(h2dimdpy_ttbb,"p_{y}^{TOTEM} (GeV/c)","p_{y}^{CMS} (GeV/c)","lin","lin");// p_{y}^{TOTEM} vs p_{y}^{CMS} TT/BB
//Histogram variable Plot title
  plot1d(hdpyAll,"#Deltap_{Y}^{CMS-TOTEM} (GeV/c)","0.002GeV/c","lin","lin");// #Delta p_{Y} CMS-TOTEM All
  plot1d(hdpy,"#Deltap_{Y}^{CMS-TOTEM} (GeV/c)","0.002GeV/c","lin","lin");// #Delta p_{Y} CMS-TOTEM
  plot1d(hdpy_diag,"#Deltap_{Y}^{CMS-TOTEM} (GeV/c)","0.002GeV/c","lin","lin");// #Delta p_{Y} CMS-TOTEM TB/BT
  plot1d(hdpy_ttbb,"#Deltap_{Y}^{CMS-TOTEM} (GeV/c)","0.002GeV/c","lin","lin");// #Delta p_{Y} CMS-TOTEM TT/BB
//Histogram variable Plot title
  plot2d(h2dimdpxAll,"p_{x}^{TOTEM} (GeV/c)","p_{x}^{CMS} (GeV/c)","lin","lin");// p_{x}^{TOTEM} vs p_{x}^{CMS} All
  plot2d(h2dimdpx,"p_{x}^{TOTEM} (GeV/c)","p_{x}^{CMS} (GeV/c)","lin","lin");// p_{x}^{TOTEM} vs p_{x}^{CMS}
  plot2d(h2dimdpx_diag,"p_{x}^{TOTEM} (GeV/c)","p_{x}^{CMS} (GeV/c)","lin","lin");// p_{x}^{TOTEM} vs p_{x}^{CMS} diag
  plot2d(h2dimdpx_ttbb,"p_{x}^{TOTEM} (GeV/c)","p_{x}^{CMS} (GeV/c)","lin","lin");// p_{x}^{TOTEM} vs p_{x}^{CMS} TT/BB
//Histogram variable Plot title
  plot1d(hdpxAll,"#Deltap_{X}^{CMS-TOTEM} (GeV/c)","0.002GeV/c","lin","lin");// #Delta p_{X} CMS-TOTEM All
  plot1d(hdpx,"#Deltap_{X}^{CMS-TOTEM} (GeV/c)","0.002GeV/c","lin","lin");// #Delta p_{X} CMS-TOTEM
  plot1d(hdpx_diag,"#Deltap_{X}^{CMS-TOTEM} (GeV/c)","0.002GeV/c","lin","lin");// #Delta p_{X} CMS-TOTEM TB/BT
  plot1d(hdpx_ttbb,"#Deltap_{X}^{CMS-TOTEM} (GeV/c)","0.002GeV/c","lin","lin");// #Delta p_{X} CMS-TOTEM TT/BB
//Histogram variable Plot title
  plot2d(h2dimxVtxRL,"x Vtx L (m)","x Vtx R (m)","lin","lin");// xVtxL vs xVtxR (m)
  plot2d(h2dimxVtxcmsR,"x Vtx CMS (cm)","x Vtx R (cm)","lin","lin");// xVtxCMS vs xVtxR (cm)
  plot2d(h2dimxVtxcmsL,"x Vtx CMS (cm)","x Vtx L (cm)","lin","lin");// xVtxCMS vs xVtxL (cm)
  plot2d(h2dimxVtxcmsRL,"x Vtx CMS (cm)","x Vtx RL (cm)","lin","lin");// xVtxCMS vs xVtxRL (cm)
  plot2d(h2dimxVtxcmsR2,"x Vtx CMS (cm)","x Vtx R (cm)","lin","lin");// xVtxCMS vs xVtxR (cm) (|xVtxL-xVtxR|<3e-5)
  plot2d(h2dimxVtxcmsL2,"x Vtx CMS (cm)","x Vtx L (cm)","lin","lin");// xVtxCMS vs xVtxL (cm) (|xVtxL-xVtxR|<3e-5)
  plot2d(h2dimxVtxcmsRL2,"x Vtx CMS (cm)","x Vtx RL (cm)","lin","lin");// xVtxCMS vs xVtxRL (cm)
  plot2d(h2dimxVtx_zVtx_CT,"x Vtx CMS - x Vtx TOTEM (cm)","z Vtx (cm)","lin","lin");// xVtxCMS-xVtxTOTEM vs zVtx (cm)
  plot2d(h2dimxVtx_zVtx_C,"x Vtx CMS (cm)","z Vtx (cm)","lin","lin");// xVtxCMS vs zVtx (cm)
  plot2d(h2dimxVtx_zVtx_T,"x Vtx TOTEM (cm)","z Vtx (cm)","lin","lin");// xVtxTOTEM vs zVtx (cm)
//Histogram variable Plot title
  plot1d(hxVtxRL,"x Vtx R-x Vtx L (m)","1.0E-06m","lin","lin");// xVtxR-xVtxL (m)
  plot1d(hxVtxcmsR,"x Vtx CMS - x Vtx R (cm)","2.0E-03cm","lin","lin");// xVtxCMS-xVtxR (cm)
  plot1d(hxVtxcmsL,"x Vtx CMS - x Vtx L (cm)","2.0E-03cm","lin","lin");// xVtxCMS-xVtxL (cm)
  plot1d(hxVtxcmsRL,"x Vtx CMS - x Vtx TOTEM (cm)","2.0E-03cm","lin","lin");// xVtxCMS-xVtxTOTEM (cm)
  plot1d(hxVtxRL_diag,"x Vtx R - x Vtx L (m)","1.0E-06m","lin","lin");// xVtxR-xVtxL (m)
  plot1d(hxVtxcmsR_diag,"x Vtx CMS - x Vtx R (cm)","2.0E-03cm","lin","lin");// xVtxCMS-xVtxR (cm)
  plot1d(hxVtxcmsL_diag,"x Vtx CMS - x Vtx L (cm)","2.0E-03cm","lin","lin");// xVtxCMS-xVtxL (cm)
  plot1d(hxVtxcmsRL_diag,"x Vtx CMS - x Vtx TOTEM (cm)","2.0E-03cm","lin","lin");// xVtxCMS-xVtxTOTEM (cm)
  plot1d(hxVtxRL_ttbb,"x Vtx R - x Vtx L (m)","1.0E-06m","lin","lin");// xVtxR-xVtxL (m)
  plot1d(hxVtxcmsR_ttbb,"x Vtx CMS - x Vtx R (cm)","2.0E-03cm","lin","lin");// xVtxCMS-xVtxR (cm)
  plot1d(hxVtxcmsL_ttbb,"x Vtx CMS - x Vtx L (cm)","2.0E-03cm","lin","lin");// xVtxCMS-xVtxL (cm)
  plot1d(hxVtxcmsRL_ttbb,"x Vtx CMS - x Vtx TOTEM (cm)","2.0E-03cm","lin","lin");// xVtxCMS-xVtxTOTEM (cm)
//Histogram variable Plot title
  plot2d(hdedx,"dE/dx (MeV/cm)","p (GeV/c)","lin","lin");// dE/dx vs p
//
  plot2d(hlndedx,"ln(dE/dx/(MeV/cm))","p (GeV/c)","lin","lin");// ln dE/dx vs p
  plot2d(hl10dedx,"log(dE/dx/(MeV/cm))","p (GeV/c)","lin","lin");// log dE/dx vs p
  //
  plot2d(phi_proton_right_t,"#varphi (rad)","|-t| (GeV^{2})","lin","lin");// #varphi vs -t
  plot2d(phi_proton_left_t,"#varphi (rad)","|-t| (GeV^{2})","lin","lin");// #varphi vs -t
  //
  plot2d(phi_proton_right_t_diag,"#varphi (rad)","|-t| (GeV^{2})","lin","lin");// #varphi vs -t
  plot2d(phi_proton_left_t_diag,"#varphi (rad)","|-t| (GeV^{2})","lin","lin");// #varphi vs -t
  //
  plot2d(phi_proton_right_t_ttbb,"#varphi (rad)","|-t| (GeV^{2})","lin","lin");// #varphi vs -t
  plot2d(phi_proton_left_t_ttbb,"#varphi (rad)","|-t| (GeV^{2})","lin","lin");// #varphi vs -t
  //
  plot2d(phi_proton_right_t_tt,"#varphi (rad)","|-t| (GeV^{2})","lin","lin");// #varphi vs -t
  plot2d(phi_proton_left_t_tt,"#varphi (rad)","|-t| (GeV^{2})","lin","lin");// #varphi vs -t
  //
  plot2d(phi_proton_right_t_bb,"#varphi (rad)","|-t| (GeV^{2})","lin","lin");// #varphi vs -t
  plot2d(phi_proton_left_t_bb,"#varphi (rad)","|-t| (GeV^{2})","lin","lin");// #varphi vs -t
  //
  plot1d(hrapy,"rapidity","0.01","lin","lin");// rapidity
  plot1d(hrapy2,"rapidity 2","0.01","lin","lin");// rapidity 2
  //
  //gSystem->Exec("cd ../");
  gSystem->cd("../");
  //
  //cout << " ...coming soon ! " << endl;
  cout << "                                 " << endl;
  cout << "    ...all set !                 " << endl;
  cout << "                                 " << endl;
  //
}

//int main()
//  {
//    makeplot();
//  }

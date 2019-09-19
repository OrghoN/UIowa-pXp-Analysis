TFile *_file0 = TFile::Open("t0RP95relui4.root")

  THStack *hs = new THStack("hs","M_{4#pi} TTBB, TBBT");

hs->Add(hm2rec2OS_ttbb3);
hs->Add(hm2rec2OS_ttbb4);
hs->Add(hm2rec2OS_ttbb5);

hs->Add(hm2rec2OS_diag3);
hs->Add(hm2rec2OS_diag4);
hs->Add(hm2rec2OS_diag5);


  hs->GetXaxis()->SetTitle("M_{4#pi} (GeV/c^{2})");
  hs->GetYaxis()->SetTitle("# of Events/0.02,0.025,0.05/GeV/c^{2}");
  hs->GetXaxis()->CenterTitle();
  hs->GetYaxis()->CenterTitle();
  hs->GetXaxis()->SetTitleOffset(1.2);
  hs->GetYaxis()->SetTitleOffset(1.2);

  gPad->SetLogy();

//auto legend = new TLegend(0.1,0.7,0.48,0.9);
auto legend = new TLegend(0.7,0.7,0.8,0.87);
legend->SetHeader("binning","L")
legend->AddEntry(hs,"125 bins: 0.0 to 2.5 GeV/c^{2}","L");
legend->AddEntry(hs," 60 bins: 2.5 to 4.0 GeV/c^{2}","L");
legend->AddEntry(hs," 80 bins: 4.0 to 8.0 GeV/c^{2}","L");
//legend->SetHeader("binning","L"); // option "C" allows to center the header
//legend->AddEntry(hs,"Histogram filled with random numbers","f");
//legend->AddEntry("f1","Function abs(#frac{sin(x)}{x})","l");
//legend->AddEntry("gr","Graph with error bars","lep");
//legend->Draw();


hs->Draw("NOSTACK");
  c1->Update();
 c1->Modified();

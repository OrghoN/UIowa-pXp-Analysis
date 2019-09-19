

TH1D *hy1 = dphi_proton_mrec->ProjectionY("cutgy1",0,55,"[cutg]")
TH1D *hy2 = dphi_proton_mrec->ProjectionY("cutgy2",56,80,"[cutg]")
TH1D *hy3 = dphi_proton_mrec->ProjectionY("cutgy3",81,100,"[cutg]")
TH1D *hy4 = dphi_proton_mrec->ProjectionY("cutgy4",101,400,"[cutg]")

  hy1->Draw("HIST");
  hy2->Draw("HIST");
  hy3->Draw("HIST");
  hy4->Draw("HIST");

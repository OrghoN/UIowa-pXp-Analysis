//gEnv->GetValue(“Canvas.SavePrecision”, -1)
//gEnv->SetValue(“Canvas.SavePrecision”, 16)

//non-resonant background (4p or rpp) of the form ...see Robert's draft
   // A · (m - B)^C · exp^(D·m)  *** 4 parameters
   Double_t nonresonantbackground(Double_t *x, Double_t *par) {
     return par[0]*pow((x[0]-par[1]),par[2])*TMath::Exp(par[3]*x[0]);
}
//Gauss function: g(x) = p0*exp(-0.5*((x-p1)/p2)^2)
Double_t mygaussfunction(Double_t *x, Double_t *par) {
  return par[0]*TMath::Exp(pow(-0.5*((x[0]-par[1])/par[2]),2));
}
// Quadratic background function
Double_t background(Double_t *x, Double_t *par) {
   return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}

// Lorenzian Peak function
Double_t lorentzianPeak(Double_t *x, Double_t *par) {
  return (0.5*par[0]*par[1]/TMath::Pi()) /
    TMath::Max( 1.e-10,(x[0]-par[2])*(x[0]-par[2])
   + .25*par[1]*par[1]);
}

// Sum of background and peak function
Double_t fitFunction(Double_t *x, Double_t *par) {
  return background(x,par) + lorentzianPeak(x,&par[3]);
}

// Sum of nonresonant background and gaussian peak function - K0sK0s
Double_t fitkoskos(Double_t *x, Double_t *par) {
  return nonresonantbackground(x,par) + mygaussfunction(x,&par[4]);
}

   TCanvas *c1 = new TCanvas("c1","K0sK0s",10,10,700,500);
   c1->SetFillColor(33);
   c1->SetFrameFillColor(41);
   c1->SetGrid();

//gEnv->Print();

//...not obsolete
//gEnv->GetValue(“Canvas.SavePrecision”, -1);
//gEnv->SetValue(“Canvas.SavePrecision”, 16);

   double xlow = 1.5 ;
   double xhigh = 1.6 ;

//   TH1F *hdpy = new TH1F("histo",
//   "Lorentzian Peak on Quadratic Background",500,-0.5,0.5);
   hm2rec2OSm1234->SetMarkerStyle(21);
   hm2rec2OSm1234->SetMarkerSize(0.8);
//   hm2rec2OSm1234->SetStats(0);

   // create a TF1 with the range from 1.0 to 2.0 and 7 parameters
   TF1 *fitK = new TF1("fitK",fitkoskos,xlow,xhigh,7);
   fitK->SetNpx(500);
   fitK->SetLineWidth(4);
   fitK->SetLineColor(kMagenta);

   // first try without starting values for the parameters
   // This defaults to 1 for each param.
   // this results in an ok fit for the polynomial function
   // however the non-linear part (lorenzian) does not
   // respond well.
   fitK->SetParameters(1,1,1,1,1,1);
   hm2rec2OSm1234->Fit("fitK","0");

   // second try: set start values for some parameters
   fitK->SetParameter(1,0.2); // 
   fitK->SetParameter(2,0.7); // 
   fitK->SetParameter(3,0.3); // 
   fitK->SetParameter(4,0.5); // 
   fitK->SetParameter(5,1.55); // width
   fitK->SetParameter(6,1.55); // peak
   fitK->SetParameter(7,0.2); // 
   hm2rec2OSm1234->GetYaxis()->SetRangeUser(0,20);
   hm2rec2OSm1234->Fit("fitK","V+","ep");

   // improve the picture:
   TF1 *backFcn = new TF1("backFcn",nonresonantbackground,xlow,xhigh,4);
   backFcn->SetLineColor(kRed);
   TF1 *signalFcn = new TF1("signalFcn",mygaussfunction,xlow,xhigh,3);
   signalFcn->SetLineColor(kBlue);
   signalFcn->SetNpx(500);
   Double_t par[7];

// writes the fit results into the par array
   fitK->GetParameters(par);

   backFcn->SetParameters(par);
   backFcn->Draw("same");

   signalFcn->SetParameters(&par[4]);
   signalFcn->Draw("same");

   // draw the legend
   TLegend *legend=new TLegend(0.585,0.683,0.879,0.850);
   legend->SetTextFont(72);
   legend->SetTextSize(0.03);
   legend->AddEntry(hm2rec2OSm1234,"Data","lpe");
   legend->AddEntry(backFcn,"Background fit","l");
   legend->AddEntry(signalFcn,"Signal fit","l");
   legend->AddEntry(fitK,"Global fit","l");
   legend->Draw();

//signalFcn->Integral(1.0,2.0);
double r = signalFcn->Integral(1.5,1.6);
std::cout << " r = " << r << std::endl

double rr = backFcn->Integral(1.5,1.6);
std::cout << " rr = " << rr << std::endl

double rrr = fitK->Integral(1.5,1.6);
std::cout << " rrr = " << rrr << std::endl

//double rrrr = fitK->Integral(1.0,3.0);
//std::cout << " rrrr = " << rrrr << std::endl

//main->Fit("func","QWEMR"); // or with likelihood: "LQWEMR"

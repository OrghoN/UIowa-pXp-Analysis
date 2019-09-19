void binwidth()
{
  //TFile *f = TFile::Open("t0RP102relui4.root");

//0-2.5 (125bins), 2.5-4(60bins), 4-8(80bins)
double xmin1 = 0.;
double xmax1 = 2.5;
int nbins1 = 125;

double xmin2 = 2.5;
double xmax2 = 4.;
int nbins2 = 60;

double xmin3 = 4.;
double xmax3 = 8.;
int nbins3 = 80;

double bwidth1 = (xmax1 - xmin1)/nbins1;
double bwidth2 = (xmax2 - xmin2)/nbins2;
double bwidth3 = (xmax3 - xmin3)/nbins3;

int nbinstot = nbins1 + nbins2 + nbins3;
 cout << " nbinstot = " << nbinstot  << "\n" ;
//...Luiz
double edges[nbinstot+1] ;

//nbinstot++;

int nbins=0;

 for( int i=0; i<nbins1; i++){ edges[nbins] = xmin1 + bwidth1 * i; nbins++; cout << " i = " << i << "\n" ;}
 cout << " a nbins = " << nbins  << "\n" ; 
 for( int i=0; i<nbins2; i++){ edges[nbins] = xmin2 + bwidth2 * i; nbins++; cout << " i = " << i << "\n" ;}
 cout << " b nbins = " << nbins  << "\n" ; 
//...Luiz
 for( int i=0; i<=nbins3; i++){ edges[nbins] = xmin3 + bwidth3 * i; nbins++; cout << " i = " << i << "\n" ;}
 cout << " c nbins = " << nbins-1  << "\n" ; 
//h = new TH1D("h1","title",nbinstot,edges);
//or
//h = new TH1D("h1","title",nbins + 1,edges);

 for( int i=0; i<=nbinstot; i++){ cout << " edges["<<i <<"] = " << edges[i] << "\n" ;}
 
 TH1* httbb2 = new TH1D("httbb2","title",nbinstot,edges);
 //TH1* hdiag2 = new TH1D("hdiag2","title",nbinstot,edges);
 //TH1* hm2rec2OS_ttbb2 = new TH1D("httbb2","title",nbinstot,edges);
 //TH1* hm2rec2OS_diag2 = new TH1D("hdiag2","title",nbinstot,edges);
 //httbb2 = hm2rec2OS_ttbb2->Clone(); won't work!

 //httbb2->Draw();
 //hdiag2->Draw("same");
}

//STANDARD ROOT INCLUDES
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TDirectory.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>
#include <TMath.h>

//OUR OWN CLASSES TO READ THE TREE
#include "UATree/UADataFormat/interface/MyEvtId.h"
#include "UATree/UADataFormat/interface/MyHLTrig.h"
#include "UATree/UADataFormat/interface/MyCaloTower.h"
#include "UATree/UADataFormat/interface/MyTracks.h"
#include "UATree/UADataFormat/interface/MyVertex.h"
#include "UATree/UADataFormat/interface/MySiPixelCluster.h"
#include "UATree/UADataFormat/interface/MySiStripCluster.h"

// TOTEM data formats
#include "RPRootDumpReconstructedProton.h"
#include "RPRootDumpReconstructedProtonPair.h"
#include "RPRootDumpTrackInfo.h"
#include "RPRootDumpDigiInfo.h"
#include "RPRootDumpPatternInfo.h"
#include "TriggerData.h"

#include "analysis_tools_Mirko2015.h"
//#include "rp_aperture_config.h"

//STANDARD C++ INCLUDES
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <cstdlib>


using namespace std;

void anaRP(vector<string> const &fileNames, string const &outputFileName, string const &outputFileTOTEM, const Int_t nevt_max);

int main(int argc, char **argv)
{

    if (argc < 5)
    {
        cout << "anaRP fileListName=filenames.txt outputFileName=output.root outputFileTOTEM=totemlist.txt nevt_max=-100" << endl;
        exit(-1);
    }

    string fileListName = argv[1];
    string outputFileName = argv[2];
    string outputFileTOTEM = argv[3];
    int nevt_max = -1;
    stringstream maxEvents_ss;
    maxEvents_ss << argv[4];
    maxEvents_ss >> nevt_max;

    ifstream infile( fileListName.c_str() );

    vector<string> fileNames;
    string file;
    while( infile >> file )
    {
        cout << "Adding " << file << endl;
        fileNames.push_back( file );
    }
    infile.close();

    anaRP(fileNames, outputFileName, outputFileTOTEM, nevt_max);
    return 0;
}

/*
void anaRP(string const& fileListName="filenames.txt", string const& outputFileName = "output.root", const Int_t nevt_max = -100){

   string file;

   ifstream infile( fileListName.c_str() );

   vector<string> fileNames;
   while( infile >> file ) {
      cout << "Adding " << file << endl;
      fileNames.push_back( file );
   }
   infile.close();

   anaRP(fileNames,outputFileName,nevt_max);
}
*/

void anaRP(vector<string> const &fileNames, string const &outputFileName = "output.root", string const &outputFileTOTEM = "totemlist.txt", const Int_t nevt_max = -100)
{

    bool isMC  = false;
    string treeName = (!isMC) ? "cms_totem" : "evt";

    double wei = 1.;

    //==============================
    const Int_t nevt_max_corr = (nevt_max >= 0) ? nevt_max : 99999999;
    cout << "nevt_max_corr = " << nevt_max_corr << endl;

    //  ofstream fout(outputFileTOTEM.c_str());


    // Declaration of histograms
    map<string, TH1F *> histosTH1F;

    vector<string> selections;
    selections.push_back("TOTEM0");
    selections.push_back("2valid");
    selections.push_back("anyTB/BT/TT/BB");
    selections.push_back("exclusiveTB/BT/TT/BB");
    selections.push_back("fiducialXY");
    selections.push_back("notElastic");
    selections.push_back("#xi<0.1");

    int nBinsEventSelection = selections.size();
    histosTH1F["EventSelection"] = new TH1F("EventSelection", " ", nBinsEventSelection, 0, nBinsEventSelection);
    for(size_t k = 0; k < selections.size(); ++k)
        histosTH1F["EventSelection"]->GetXaxis()->SetBinLabel( (k + 1), selections[k].c_str() );

    histosTH1F["hnconf"] = new TH1F("hnconf", "Number of configurations (TB or BT or TT or BB)", 5, 0., 5.);

    histosTH1F["rp_x_020"] = new TH1F("rp_x_020", "x RP", 200, -10., 10.);
    histosTH1F["rp_x_021"] = new TH1F("rp_x_021", "x RP", 200, -10., 10.);
    histosTH1F["rp_x_024"] = new TH1F("rp_x_024", "x RP", 200, -10., 10.);
    histosTH1F["rp_x_025"] = new TH1F("rp_x_025", "x RP", 200, -10., 10.);

    histosTH1F["rp_x_120"] = new TH1F("rp_x_120", "x RP", 200, -10., 10.);
    histosTH1F["rp_x_121"] = new TH1F("rp_x_121", "x RP", 200, -10., 10.);
    histosTH1F["rp_x_124"] = new TH1F("rp_x_124", "x RP", 200, -10., 10.);
    histosTH1F["rp_x_125"] = new TH1F("rp_x_125", "x RP", 200, -10., 10.);

    histosTH1F["rp_y_020"] = new TH1F("rp_y_020", "y RP", 500, -50., 50.);
    histosTH1F["rp_y_021"] = new TH1F("rp_y_021", "y RP", 500, -50., 50.);
    histosTH1F["rp_y_024"] = new TH1F("rp_y_024", "y RP", 500, -50., 50.);
    histosTH1F["rp_y_025"] = new TH1F("rp_y_025", "y RP", 500, -50., 50.);

    histosTH1F["rp_y_120"] = new TH1F("rp_y_120", "y RP", 500, -50., 50.);
    histosTH1F["rp_y_121"] = new TH1F("rp_y_121", "y RP", 500, -50., 50.);
    histosTH1F["rp_y_124"] = new TH1F("rp_y_124", "y RP", 500, -50., 50.);
    histosTH1F["rp_y_125"] = new TH1F("rp_y_125", "y RP", 500, -50., 50.);

    //--- from TT/BB, above from TB/BT
    histosTH1F["rp2_x_020"] = new TH1F("rp2_x_020", "x RP", 200, -10., 10.);
    histosTH1F["rp2_x_021"] = new TH1F("rp2_x_021", "x RP", 200, -10., 10.);
    histosTH1F["rp2_x_024"] = new TH1F("rp2_x_024", "x RP", 200, -10., 10.);
    histosTH1F["rp2_x_025"] = new TH1F("rp2_x_025", "x RP", 200, -10., 10.);

    histosTH1F["rp2_x_120"] = new TH1F("rp2_x_120", "x RP", 200, -10., 10.);
    histosTH1F["rp2_x_121"] = new TH1F("rp2_x_121", "x RP", 200, -10., 10.);
    histosTH1F["rp2_x_124"] = new TH1F("rp2_x_124", "x RP", 200, -10., 10.);
    histosTH1F["rp2_x_125"] = new TH1F("rp2_x_125", "x RP", 200, -10., 10.);

    histosTH1F["rp2_y_020"] = new TH1F("rp2_y_020", "y RP", 500, -50., 50.);
    histosTH1F["rp2_y_021"] = new TH1F("rp2_y_021", "y RP", 500, -50., 50.);
    histosTH1F["rp2_y_024"] = new TH1F("rp2_y_024", "y RP", 500, -50., 50.);
    histosTH1F["rp2_y_025"] = new TH1F("rp2_y_025", "y RP", 500, -50., 50.);

    histosTH1F["rp2_y_120"] = new TH1F("rp2_y_120", "y RP", 500, -50., 50.);
    histosTH1F["rp2_y_121"] = new TH1F("rp2_y_121", "y RP", 500, -50., 50.);
    histosTH1F["rp2_y_124"] = new TH1F("rp2_y_124", "y RP", 500, -50., 50.);
    histosTH1F["rp2_y_125"] = new TH1F("rp2_y_125", "y RP", 500, -50., 50.);

    histosTH1F["thyEla"] = new TH1F("thyEla", "thyL+thyR", 4000, -0.0004, 0.0004);
    histosTH1F["thxEla"] = new TH1F("thxEla", "thxL+thxR", 4000, -0.0004, 0.0004);

    histosTH1F["thyEla_diag"] = new TH1F("thyEla_diag", "thyL+thyR diagonals", 4000, -0.0004, 0.0004);
    histosTH1F["thxEla_diag"] = new TH1F("thxEla_diag", "thxL+thxR diagonals", 4000, -0.0004, 0.0004);

    histosTH1F["thyEla_ttbb"] = new TH1F("thyEla_ttbb", "thyL+thyR TT/BB", 4000, -0.0004, 0.0004);
    histosTH1F["thxEla_ttbb"] = new TH1F("thxEla_ttbb", "thxL+thxR TT/BB", 4000, -0.0004, 0.0004);

    //histosTH1F["proton_right_xi"] = new TH1F("proton_right_xi", "#xi" , 200 , -1. , 1.);
    //...Luiz
    histosTH1F["proton_right_xi"] = new TH1F("proton_right_xi", "#xi", 200, -0.5, 0.5);
    histosTH1F["proton_right_logXi"] = new TH1F("proton_right_logXi", "log(#xi)", 200, -5., 0.);

    //histosTH1F["proton_left_xi"] = new TH1F("proton_left_xi", "#xi" , 200 , -1. , 1.);
    //...Luiz
    histosTH1F["proton_left_xi"] = new TH1F("proton_left_xi", "#xi", 200, -0.5, 0.5);
    histosTH1F["proton_left_logXi"] = new TH1F("proton_left_logXi", "log(#xi)", 200, -5., 0.);

    histosTH1F["proton_right_t"] = new TH1F("proton_right_t", "-t", 1000, 0., 5.);
    histosTH1F["proton_left_t"] = new TH1F("proton_left_t", "-t", 1000, 0., 5.);

    histosTH1F["proton_right_t_diag"] = new TH1F("proton_right_t_diag", "-t diagonal", 1000, 0., 5.);
    histosTH1F["proton_left_t_diag"] = new TH1F("proton_left_t_diag", "-t diagonal", 1000, 0., 5.);

    histosTH1F["proton_right_t_ttbb"] = new TH1F("proton_right_t_ttbb", "-t", 1000, 0., 5.);
    histosTH1F["proton_left_t_ttbb"] = new TH1F("proton_left_t_ttbb", "-t TT/BB", 1000, 0., 5.);

    histosTH1F["eHF"] = new TH1F("eHF", "energy HF tower (GeV)", 500, 0., 100.);
    histosTH1F["nHF"] = new TH1F("nHF", "n HF tower (eHF>5 GeV)", 200, 0., 200.);

    histosTH1F["totem_py"] = new TH1F("totem_py", "p_{Y} TOTEM", 500, -5., 5.);
    histosTH1F["totem_px"] = new TH1F("totem_px", "p_{X} TOTEM", 500, -5., 5.);

    //...Luiz
    histosTH1F["totem_pyy"] = new TH1F("totem_pyy", "p_{Y} TOTEM", 1000, -1., 1.);
    histosTH1F["totem_pxx"] = new TH1F("totem_pxx", "p_{X} TOTEM", 1000, -1., 1.);

    histosTH1F["proton_dx0"]  = new TH1F("proton_dx0", "xVtx_{56}-xVtx_{45}", 300, -0.0003, 0.0003);

    histosTH1F["hLS"] = new TH1F("hLS", "LS", 800, 0., 800.);

    histosTH1F["htopo"] = new TH1F("htopo", "1=TB 2=BT 3=TT 4=BB topology", 5, 0, 5);

    //histosTH1F["hthyEla2_diag"] = new TH1F("hthyEla2_diag", "thyL+thyR dig" , 2000 , -0.0004 , 0.0004);
    //histosTH1F["hthxEla2_diag"] = new TH1F("hthxEla2_diag", "thxL+thxR dig" , 2000 , -0.0004 , 0.0004);
    //histosTH1F["hthyEla2_ttbb"] = new TH1F("hthyEla2_ttbb", "thyL+thyR TTBB" , 2000 , -0.0004 , 0.0004);
    //histosTH1F["hthxEla2_ttbb"] = new TH1F("hthxEla2_ttbb", "thxL+thxR TTBB" , 2000 , -0.0004 , 0.0004);
    //...Luiz
    histosTH1F["hthyEla2_diag"] = new TH1F("hthyEla2_diag", "thyL+thyR dig", 4000, -0.0004, 0.0004);
    histosTH1F["hthxEla2_diag"] = new TH1F("hthxEla2_diag", "thxL+thxR dig", 4000, -0.0004, 0.0004);
    histosTH1F["hthyEla2_ttbb"] = new TH1F("hthyEla2_ttbb", "thyL+thyR TTBB", 4000, -0.0004, 0.0004);
    histosTH1F["hthxEla2_ttbb"] = new TH1F("hthxEla2_ttbb", "thxL+thxR TTBB", 4000, -0.0004, 0.0004);

    map<string, TH2F *> histosTH2F;

    histosTH2F["rp_yx_020"] = new TH2F("rp_yx_020", "y vs x RP", 200, -10., 10., 500, -50., 50.);
    histosTH2F["rp_yx_021"] = new TH2F("rp_yx_021", "y vs x RP", 200, -10., 10., 500, -50., 50.);
    histosTH2F["rp_yx_024"] = new TH2F("rp_yx_024", "y vs x RP", 200, -10., 10., 500, -50., 50.);
    histosTH2F["rp_yx_025"] = new TH2F("rp_yx_025", "y vs x RP", 200, -10., 10., 500, -50., 50.);

    histosTH2F["rp_yx_120"] = new TH2F("rp_yx_120", "y vs x RP", 200, -10., 10., 500, -50., 50.);
    histosTH2F["rp_yx_121"] = new TH2F("rp_yx_121", "y vs x RP", 200, -10., 10., 500, -50., 50.);
    histosTH2F["rp_yx_124"] = new TH2F("rp_yx_124", "y vs x RP", 200, -10., 10., 500, -50., 50.);
    histosTH2F["rp_yx_125"] = new TH2F("rp_yx_125", "y vs x RP", 200, -10., 10., 500, -50., 50.);

    //--- from TT/BB, above from TB/BT
    histosTH2F["rp2_yx_020"] = new TH2F("rp2_yx_020", "y vs x RP", 200, -10., 10., 500, -50., 50.);
    histosTH2F["rp2_yx_021"] = new TH2F("rp2_yx_021", "y vs x RP", 200, -10., 10., 500, -50., 50.);
    histosTH2F["rp2_yx_024"] = new TH2F("rp2_yx_024", "y vs x RP", 200, -10., 10., 500, -50., 50.);
    histosTH2F["rp2_yx_025"] = new TH2F("rp2_yx_025", "y vs x RP", 200, -10., 10., 500, -50., 50.);

    histosTH2F["rp2_yx_120"] = new TH2F("rp2_yx_120", "y vs x RP", 200, -10., 10., 500, -50., 50.);
    histosTH2F["rp2_yx_121"] = new TH2F("rp2_yx_121", "y vs x RP", 200, -10., 10., 500, -50., 50.);
    histosTH2F["rp2_yx_124"] = new TH2F("rp2_yx_124", "y vs x RP", 200, -10., 10., 500, -50., 50.);
    histosTH2F["rp2_yx_125"] = new TH2F("rp2_yx_125", "y vs x RP", 200, -10., 10., 500, -50., 50.);

    histosTH2F["proton_x0_RvsL"]  = new TH2F("proton_x0_RvsL", "xVtx_{56} vs xVtx_{45}", 3000, -0.005, 0.001, 3000, -0.005, 0.001);

    //...Luiz
    histosTH2F["phi_proton_right_t"] = new TH2F("phi_proton_right_t", "#varphi vs |-t|", 1000, 0., 5., 64, -3.2, 3.2);
    histosTH2F["phi_proton_left_t"]  = new TH2F("phi_proton_left_t", "#varphi vs |-t|", 1000, 0., 5., 64, -3.2, 3.2);
    //...Luiz
    histosTH2F["phi_proton_right_t_diag"] = new TH2F("phi_proton_right_t_diag", "#varphi vs |-t|", 1000, 0., 5., 64, -3.2, 3.2);
    histosTH2F["phi_proton_left_t_diag"]  = new TH2F("phi_proton_left_t_diag", "#varphi vs |-t|", 1000, 0., 5., 64, -3.2, 3.2);
    //...Luiz
    histosTH2F["phi_proton_right_t_ttbb"] = new TH2F("phi_proton_right_t_ttbb", "#varphi vs |-t|", 1000, 0., 5., 64, -3.2, 3.2);
    histosTH2F["phi_proton_left_t_ttbb"]  = new TH2F("phi_proton_left_t_ttbb", "#varphi vs |-t|", 1000, 0., 5., 64, -3.2, 3.2);
    //...Luiz
    histosTH2F["phi_proton_right_t_tt"] = new TH2F("phi_proton_right_t_tt", "#varphi vs |-t|", 1000, 0., 5., 64, -3.2, 3.2);
    histosTH2F["phi_proton_left_t_tt"]  = new TH2F("phi_proton_left_t_tt", "#varphi vs |-t|", 1000, 0., 5., 64, -3.2, 3.2);
    //...Luiz
    histosTH2F["phi_proton_right_t_bb"] = new TH2F("phi_proton_right_t_bb", "#varphi vs |-t|", 1000, 0., 5., 64, -3.2, 3.2);
    histosTH2F["phi_proton_left_t_bb"]  = new TH2F("phi_proton_left_t_bb", "#varphi vs |-t|", 1000, 0., 5., 64, -3.2, 3.2);
    //

    //---------------------------------------------------

    int nbins_eta = 80;
    int nbins_pt = 100;
    int nbins_phi = 64;

    histosTH1F["hlooper"] = new TH1F("hlooper", "isLooper", 5, 0, 5);

    histosTH1F["hpt"] = new TH1F("hpt", "p_{T}", nbins_pt, 0, 5);
    histosTH1F["heta"] = new TH1F("heta", "#eta", nbins_eta, -4, 4);
    histosTH1F["hphi"] = new TH1F("hphi", "#varphi", nbins_phi, -3.2, 3.2);

    //histosTH1F["hptP"] = new TH1F("hptP","p_{T} #pi+",nbins_pt,0,3);
    //...Luiz
    histosTH1F["hptP"] = new TH1F("hptP", "p_{T} #pi+", 2.0 * nbins_pt, 0, 4);
    histosTH1F["hetaP"] = new TH1F("hetaP", "#eta #pi+", nbins_eta, -4, 4);
    histosTH1F["hphiP"] = new TH1F("hphiP", "#varphi #pi+", nbins_phi, -3.2, 3.2);

    //histosTH1F["hptM"] = new TH1F("hptM","p_{T} #pi-",nbins_pt,0,3);
    //...Luiz
    histosTH1F["hptM"] = new TH1F("hptM", "p_{T} #pi-", 2.0 * nbins_pt, 0, 4);
    histosTH1F["hetaM"] = new TH1F("hetaM", "#eta #pi-", nbins_eta, -4, 4);
    histosTH1F["hphiM"] = new TH1F("hphiM", "#varphi #pi-", nbins_phi, -3.2, 3.2);

    //histosTH1F["hptRes"] = new TH1F("hptRes","p_{T} #pi#pi",nbins_pt,0,3);
    //...Luiz
    histosTH1F["hptRes"] = new TH1F("hptRes", "p_{T} #pi#pi", 2.0 * nbins_pt, 0, 4);
    histosTH1F["hetaRes"] = new TH1F("hetaRes", "#eta #pi#pi", nbins_eta * 1.5, -6, 6);
    histosTH1F["hphiRes"] = new TH1F("hphiRes", "#varphi #pi#pi", nbins_phi, -3.2, 3.2);

    //  histosTH1F["htopo"] = new TH1F("htopo","1=TB 2=BT 3=TT 4=BB topology",5,0,5);

    //histosTH1F["hthyEla_diag"] = new TH1F("hthyEla_diag", "thyL+thyR dig" , 2000 , -0.0004 , 0.0004);
    //histosTH1F["hthxEla_diag"] = new TH1F("hthxEla_diag", "thxL+thxR dig" , 2000 , -0.0004 , 0.0004);
    //histosTH1F["hthyEla_ttbb"] = new TH1F("hthyEla_ttbb", "thyL+thyR TTBB" , 2000 , -0.0004 , 0.0004);
    //histosTH1F["hthxEla_ttbb"] = new TH1F("hthxEla_ttbb", "thxL+thxR TTBB" , 2000 , -0.0004 , 0.0004);
    //...Luiz
    histosTH1F["hthyEla_diag"] = new TH1F("hthyEla_diag", "thyL+thyR dig", 4000, -0.0004, 0.0004);
    histosTH1F["hthxEla_diag"] = new TH1F("hthxEla_diag", "thxL+thxR dig", 4000, -0.0004, 0.0004);
    histosTH1F["hthyEla_ttbb"] = new TH1F("hthyEla_ttbb", "thyL+thyR TTBB", 4000, -0.0004, 0.0004);
    histosTH1F["hthxEla_ttbb"] = new TH1F("hthxEla_ttbb", "thxL+thxR TTBB", 4000, -0.0004, 0.0004);

    histosTH1F["hntrk0"] = new TH1F("hntrk0", "Ntrk", 150, 0, 150);
    histosTH1F["hntrk"] = new TH1F("hntrk", "Ntrk for nPixelHits>0", 150, 0, 150);
    histosTH1F["hntrkvtx"] = new TH1F("hntrkvtx", "Ntrkvtx", 150, 0, 150);
    histosTH1F["hntrkntrkvtx2"] = new TH1F("hntrkntrkvtx2", "Ntrk for Ntrkvtx==2", 150, 0, 150);
    histosTH1F["hntrk2ntrkvtx"] = new TH1F("hntrk2ntrkvtx", "Ntrkvtx for Ntrk==2", 150, 0, 150);

    histosTH2F["hntrkntrkvtx"] = new TH2F("hntrkntrkvtx", "Ntrk vs Ntrkvtx", 150, 0, 150, 150, 0, 150);

    histosTH1F["hvtx"] = new TH1F("hvtx", "vtx.isFake()", 2, 0, 2);
    histosTH1F["hvtx2"] = new TH1F("hvtx2", "vtx.isFake() 2 tracks", 2, 0, 2);
    histosTH1F["hvtx3"] = new TH1F("hvtx3", "vtx.isFake() 2 tracks both |#eta|<2.5 and OS", 2, 0, 2);

    histosTH1F["hnvtx"] = new TH1F("hnvtx", "Nvtx", 10, 0, 10);
    histosTH1F["hvtxx"] = new TH1F("hvtxx", "X vtx", 1000, -1., 1.);
    histosTH1F["hvtxy"] = new TH1F("hvtxy", "Y vtx", 1000, -1., 1.);
    histosTH1F["hvtxz"] = new TH1F("hvtxz", "Z vtx", 300, -30., 30.);
    //histosTH1F["hvtxchi2"] = new TH1F("hvtxchi2","chi2 vtx",1100,-100.,1000.);
    //histosTH1F["hvtxchi2fin"] = new TH1F("hvtxchi2fin","chi2 vtx",1100,-100.,1000.);
    //...Luiz
    histosTH1F["hvtxchi2"] = new TH1F("hvtxchi2", "#chi^{2} vtx", 1100, -100., 1000.);
    histosTH1F["hvtxchi2fin"] = new TH1F("hvtxchi2fin", "#chi^{2} vtx fin", 1100, -100., 1000.);

    //histosTH1F["heHF"] = new TH1F("heHF","HF tower energy",550,-10,100);
    //...Luiz
    histosTH1F["heHF"] = new TH1F("heHF", "HF tower energy (GeV)", 550, -10, 100);
    histosTH1F["hnHF"] = new TH1F("hnHF", "n HF towers (E>5 GeV)", 200, 0, 200);

    histosTH1F["hxiL"] = new TH1F("hxiL", "#xiL ", 100, -0.1, 0.1);
    histosTH1F["hxiR"] = new TH1F("hxiR", "#xiR ", 100, -0.1, 0.1);
    //...Luiz
    histosTH1F["hrapy"] = new TH1F("hrapy", "rapidity", 2000, -10, 10);
    //
    histosTH1F["hxiL2"] = new TH1F("hxiL2", "#xiL ", 100, -0.1, 0.1);
    histosTH1F["hxiR2"] = new TH1F("hxiR2", "#xiR ", 100, -0.1, 0.1);
    //...Luiz
    histosTH1F["hrapy2"] = new TH1F("hrapy2", "rapidity 2", 2000, -10, 10);
    //

    int massbins = 250;

    histosTH1F["hm"] = new TH1F("hm", "M_{#pi#pi} ", massbins, 0, 5.);
    //...Luiz
    //  histosTH1F["hmxicut"] = new TH1F("hmxicit","M_{#pi#pi} ",massbins,0,5.);
    histosTH1F["hmxicut"] = new TH1F("hmxicut", "M_{#pi#pi} ", massbins, 0, 5.);

    histosTH1F["hm2rec"] = new TH1F("hm2rec", "M_{#pi#pi} ", massbins, 0, 5.);
    histosTH1F["hm2recbis"] = new TH1F("hm2recbis", "M_{#pi#pi}", 2 * massbins, 0, 5.);

    //...Luiz
    histosTH1F["hm2recPPPP"] = new TH1F("hm2recPPPP", "M_{4#pi} ", massbins, 0, 5.);
    //histosTH1F["hm2recPP"] = new TH1F("hm2recPP","M_{#pi#pi} ",massbins,0,5.);
    histosTH1F["hm2recKK"] = new TH1F("hm2recKK", "M_{KK} ", massbins, 0, 5.);
    histosTH1F["hm2recMM"] = new TH1F("hm2recMM", "M_{#mu#mu} ", massbins, 0, 5.);
    histosTH1F["hm2recEE"] = new TH1F("hm2recEE", "M_{ee} ", massbins, 0, 5.);

    histosTH1F["hm2recOS"] = new TH1F("hm2recOS", "M_{#pi#pi} OS", massbins, 0, 5.);
    histosTH1F["hm2recSS"] = new TH1F("hm2recSS", "M_{#pi#pi} SS", massbins, 0, 5.);
    histosTH1F["hm2recOS_diag"] = new TH1F("hm2recOS_diag", "M_{#pi#pi} TB/BT OS", massbins, 0, 5.);
    histosTH1F["hm2recSS_diag"] = new TH1F("hm2recSS_diag", "M_{#pi#pi} TB/BT SS", massbins, 0, 5.);
    histosTH1F["hm2recOS_ttbb"] = new TH1F("hm2recOS_ttbb", "M_{#pi#pi} TT/BB OS", massbins, 0, 5.);
    histosTH1F["hm2recSS_ttbb"] = new TH1F("hm2recSS_ttbb", "M_{#pi#pi} TT/BB SS", massbins, 0, 5.);

    histosTH1F["hm2rec2OS"] = new TH1F("hm2rec2OS", "M_{#pi#pi} OS", massbins, 0, 5.);
    histosTH1F["hm2rec2OS_k1k2"] = new TH1F("hm2recOS_k1k2", "M_{k1k2} OS", 2.0 * massbins, 0, 10.);
    histosTH1F["hm2rec2OS_k3k4"] = new TH1F("hm2recOS_k3k4", "M_{k3k4} OS", 2.0 * massbins, 0, 10.);
    histosTH1F["hm2rec2OS_k1k3"] = new TH1F("hm2recOS_k1k3", "M_{k1k3} OS", 2.0 * massbins, 0, 10.);
    histosTH1F["hm2rec2OS_k2k4"] = new TH1F("hm2recOS_k2k4", "M_{k2k4} OS", 2.0 * massbins, 0, 10.);
    histosTH1F["hm2rec2OS_k1k2v2"] = new TH1F("hm2recOS_k1k2v2", "M_{k1k2} OS", 2.0 * massbins, 0, 10.);
    histosTH1F["hm2rec2OS_k3k4v2"] = new TH1F("hm2recOS_k3k4v2", "M_{k3k4} OS", 2.0 * massbins, 0, 10.);
    histosTH1F["hm2rec2OS_k1k3v2"] = new TH1F("hm2recOS_k1k3v2", "M_{k1k3} OS", 2.0 * massbins, 0, 10.);
    histosTH1F["hm2rec2OS_k2k4v2"] = new TH1F("hm2recOS_k2k4v2", "M_{k2k4} OS", 2.0 * massbins, 0, 10.);
    histosTH1F["hm2rec2SS"] = new TH1F("hm2rec2SS", "M_{#pi#pi} SS", massbins, 0, 5.);
    histosTH1F["hm2rec2OS_diag"] = new TH1F("hm2rec2OS_diag", "M_{#pi#pi} TB/BT OS", massbins, 0, 5.);
    histosTH1F["hm2rec2SS_diag"] = new TH1F("hm2rec2SS_diag", "M_{#pi#pi} TB/BT SS", massbins, 0, 5.);
    histosTH1F["hm2rec2OS_ttbb"] = new TH1F("hm2rec2OS_ttbb", "M_{#pi#pi} TT/BB OS", massbins, 0, 5.);
    histosTH1F["hm2rec2SS_ttbb"] = new TH1F("hm2rec2SS_ttbb", "M_{#pi#pi} TT/BB SS", massbins, 0, 5.);

    //histosTH1F["hm2rec2OS_diag_trkP"] = new TH1F("hm2rec2OS_diag_trkP","M_{#pi#pi} TB/BT OS, py_{#pi1}py_{#pi2}>0",massbins,0,5.);
    //histosTH1F["hm2rec2OS_diag_trkM"] = new TH1F("hm2rec2OS_diag_trkM","M_{#pi#pi} TB/BT OS, py_{#pi1}py_{#pi2}<0",massbins,0,5.);
    //histosTH1F["hm2rec2OS_ttbb_trkP"] = new TH1F("hm2rec2OS_ttbb_trkP","M_{#pi#pi} TT/BB OS, py_{#pi1}py_{#pi2}>0",massbins,0,5.);
    //histosTH1F["hm2rec2OS_ttbb_trkM"] = new TH1F("hm2rec2OS_ttbb_trkM","M_{#pi#pi} TT/BB OS, py_{#pi1}py_{#pi2}<0",massbins,0,5.);
    //...Luiz
    histosTH1F["hm2rec2OS_diag_trkP"] = new TH1F("hm2rec2OS_diag_trkP", "M_{#pi#pi} TB/BT OS, py_{#pi_{1}} py_{#pi_{2}}>0", massbins, 0, 5.);
    histosTH1F["hm2rec2OS_diag_trkM"] = new TH1F("hm2rec2OS_diag_trkM", "M_{#pi#pi} TB/BT OS, py_{#pi_{1}} py_{#pi_{2}}<0", massbins, 0, 5.);
    histosTH1F["hm2rec2OS_ttbb_trkP"] = new TH1F("hm2rec2OS_ttbb_trkP", "M_{#pi#pi} TT/BB OS, py_{#pi_{1}} py_{#pi_{2}}>0", massbins, 0, 5.);
    histosTH1F["hm2rec2OS_ttbb_trkM"] = new TH1F("hm2rec2OS_ttbb_trkM", "M_{#pi#pi} TT/BB OS, py_{#pi_{1}} py_{#pi_{2}}<0", massbins, 0, 5.);

    histosTH1F["hm2rec2OS_diag_pypxP"] = new TH1F("hm2rec2OS_diag_pypxP", "M_{#pi#pi} TB/BT OS, |py/px|_{#pi#pi} > 1", massbins, 0, 5.);
    histosTH1F["hm2rec2OS_diag_pypxM"] = new TH1F("hm2rec2OS_diag_pypxM", "M_{#pi#pi} TB/BT OS, |py/px|_{#pi#pi} < 1", massbins, 0, 5.);
    histosTH1F["hm2rec2OS_ttbb_pypxP"] = new TH1F("hm2rec2OS_ttbb_pypxP", "M_{#pi#pi} TT/BB OS, |py/px|_{#pi#pi} > 1", massbins, 0, 5.);
    histosTH1F["hm2rec2OS_ttbb_pypxM"] = new TH1F("hm2rec2OS_ttbb_pypxM", "M_{#pi#pi} TT/BB OS, |py/px|_{#pi#pi} < 1", massbins, 0, 5.);

    histosTH1F["hm2rec3OS"] = new TH1F("hm2rec3OS", "M_{#pi#pi} OS", massbins, 0, 5.);
    histosTH1F["hm2rec3SS"] = new TH1F("hm2rec3SS", "M_{#pi#pi} SS", massbins, 0, 5.);
    histosTH1F["hm2rec3OS_diag"] = new TH1F("hm2rec3OS_diag", "M_{#pi#pi} TB/BT OS", massbins, 0, 5.);
    histosTH1F["hm2rec3SS_diag"] = new TH1F("hm2rec3SS_diag", "M_{#pi#pi} TB/BT SS", massbins, 0, 5.);
    histosTH1F["hm2rec3OS_ttbb"] = new TH1F("hm2rec3OS_ttbb", "M_{#pi#pi} TT/BB OS", massbins, 0, 5.);
    histosTH1F["hm2rec3SS_ttbb"] = new TH1F("hm2rec3SS_ttbb", "M_{#pi#pi} TT/BB SS", massbins, 0, 5.);

    //histosTH1F["hm2rec3OS_diag_trkP"] = new TH1F("hm2rec3OS_diag_trkP","M_{#pi#pi} TB/BT OS, py_{#pi1}py_{#pi2}>0",massbins,0,5.);
    //histosTH1F["hm2rec3OS_diag_trkM"] = new TH1F("hm2rec3OS_diag_trkM","M_{#pi#pi} TB/BT OS, py_{#pi1}py_{#pi2}<0",massbins,0,5.);
    //histosTH1F["hm2rec3OS_ttbb_trkP"] = new TH1F("hm2rec3OS_ttbb_trkP","M_{#pi#pi} TT/BB OS, py_{#pi1}py_{#pi2}>0",massbins,0,5.);
    //histosTH1F["hm2rec3OS_ttbb_trkM"] = new TH1F("hm2rec3OS_ttbb_trkM","M_{#pi#pi} TT/BB OS, py_{#pi1}py_{#pi2}<0",massbins,0,5.);
    //...Luiz
    histosTH1F["hm2rec3OS_diag_trkP"] = new TH1F("hm2rec3OS_diag_trkP", "M_{#pi#pi} TB/BT OS, py_{#pi_{1}}py_{#pi_{2}}>0", massbins, 0, 5.);
    histosTH1F["hm2rec3OS_diag_trkM"] = new TH1F("hm2rec3OS_diag_trkM", "M_{#pi#pi} TB/BT OS, py_{#pi_{1}}py_{#pi_{2}}<0", massbins, 0, 5.);
    histosTH1F["hm2rec3OS_ttbb_trkP"] = new TH1F("hm2rec3OS_ttbb_trkP", "M_{#pi#pi} TT/BB OS, py_{#pi_{1}}py_{#pi_{2}}>0", massbins, 0, 5.);
    histosTH1F["hm2rec3OS_ttbb_trkM"] = new TH1F("hm2rec3OS_ttbb_trkM", "M_{#pi#pi} TT/BB OS, py_{#pi_{1}}py_{#pi_{2}}<0", massbins, 0, 5.);

    histosTH1F["hm2rec3OS_diag_pypxP"] = new TH1F("hm2rec3OS_diag_pypxP", "M_{#pi#pi} TB/BT OS, |py/px|_{#pi#pi} > 1", massbins, 0, 5.);
    histosTH1F["hm2rec3OS_diag_pypxM"] = new TH1F("hm2rec3OS_diag_pypxM", "M_{#pi#pi} TB/BT OS, |py/px|_{#pi#pi} < 1", massbins, 0, 5.);
    histosTH1F["hm2rec3OS_ttbb_pypxP"] = new TH1F("hm2rec3OS_ttbb_pypxP", "M_{#pi#pi} TT/BB OS, |py/px|_{#pi#pi} > 1", massbins, 0, 5.);
    histosTH1F["hm2rec3OS_ttbb_pypxM"] = new TH1F("hm2rec3OS_ttbb_pypxM", "M_{#pi#pi} TT/BB OS, |py/px|_{#pi#pi} < 1", massbins, 0, 5.);

    histosTH1F["hm2rec4OS"] = new TH1F("hm2rec4OS", "M_{#pi#pi} OS", massbins, 0, 5.);
    histosTH1F["hm2rec4SS"] = new TH1F("hm2rec4SS", "M_{#pi#pi} SS", massbins, 0, 5.);
    histosTH1F["hm2rec4OS_diag"] = new TH1F("hm2rec4OS_diag", "M_{#pi#pi} TB/BT OS", massbins, 0, 5.);
    histosTH1F["hm2rec4SS_diag"] = new TH1F("hm2rec4SS_diag", "M_{#pi#pi} TB/BT SS", massbins, 0, 5.);
    histosTH1F["hm2rec4OS_ttbb"] = new TH1F("hm2rec4OS_ttbb", "M_{#pi#pi} TT/BB OS", massbins, 0, 5.);
    histosTH1F["hm2rec4SS_ttbb"] = new TH1F("hm2rec4SS_ttbb", "M_{#pi#pi} TT/BB SS", massbins, 0, 5.);

    //histosTH1F["hm2rec4OS_diag_trkP"] = new TH1F("hm2rec4OS_diag_trkP","M_{#pi#pi}TB/BT OS, py_{#pi1}py_{#pi2}>0",massbins,0,5.);
    //histosTH1F["hm2rec4OS_diag_trkM"] = new TH1F("hm2rec4OS_diag_trkM","M_{#pi#pi}TB/BT OS, py_{#pi1}py_{#pi2}<0",massbins,0,5.);
    //histosTH1F["hm2rec4OS_ttbb_trkP"] = new TH1F("hm2rec4OS_ttbb_trkP","M_{#pi#pi}TT/BB OS, py_{#pi1}py_{#pi2}>0",massbins,0,5.);
    //histosTH1F["hm2rec4OS_ttbb_trkM"] = new TH1F("hm2rec4OS_ttbb_trkM","M_{#pi#pi}TT/BB OS, py_{#pi1}py_{#pi2}<0",massbins,0,5.);
    //...Luiz
    histosTH1F["hm2rec4OS_diag_trkP"] = new TH1F("hm2rec4OS_diag_trkP", "M_{#pi#pi}TB/BT OS, py_{#pi_{1}}py_{#pi_{2}}>0", massbins, 0, 5.);
    histosTH1F["hm2rec4OS_diag_trkM"] = new TH1F("hm2rec4OS_diag_trkM", "M_{#pi#pi}TB/BT OS, py_{#pi_{1}}py_{#pi_{2}}<0", massbins, 0, 5.);
    histosTH1F["hm2rec4OS_ttbb_trkP"] = new TH1F("hm2rec4OS_ttbb_trkP", "M_{#pi#pi}TT/BB OS, py_{#pi_{1}}py_{#pi_{2}}>0", massbins, 0, 5.);
    histosTH1F["hm2rec4OS_ttbb_trkM"] = new TH1F("hm2rec4OS_ttbb_trkM", "M_{#pi#pi}TT/BB OS, py_{#pi_{1}}py_{#pi_{2}}<0", massbins, 0, 5.);

    histosTH1F["hm2rec4OS_diag_pypxP"] = new TH1F("hm2rec4OS_diag_pypxP", "M_{#pi#pi} TB/BT OS, |py/px|_{#pi#pi} > 1", massbins, 0, 5.);
    histosTH1F["hm2rec4OS_diag_pypxM"] = new TH1F("hm2rec4OS_diag_pypxM", "M_{#pi#pi} TB/BT OS, |py/px|_{#pi#pi} < 1", massbins, 0, 5.);
    histosTH1F["hm2rec4OS_ttbb_pypxP"] = new TH1F("hm2rec4OS_ttbb_pypxP", "M_{#pi#pi} TT/BB OS, |py/px|_{#pi#pi} > 1", massbins, 0, 5.);
    histosTH1F["hm2rec4OS_ttbb_pypxM"] = new TH1F("hm2rec4OS_ttbb_pypxM", "M_{#pi#pi} TT/BB OS, |py/px|_{#pi#pi} < 1", massbins, 0, 5.);

    histosTH1F["hm2rec5OS"] = new TH1F("hm2rec5OS", "M_{#pi#pi} OS", massbins, 0, 5.);
    histosTH1F["hm2rec5SS"] = new TH1F("hm2rec5SS", "M_{#pi#pi} SS", massbins, 0, 5.);
    histosTH1F["hm2rec5OS_diag"] = new TH1F("hm2rec5OS_diag", "M_{#pi#pi} TB/BT OS", massbins, 0, 5.);
    histosTH1F["hm2rec5SS_diag"] = new TH1F("hm2rec5SS_diag", "M_{#pi#pi} TB/BT SS", massbins, 0, 5.);
    histosTH1F["hm2rec5OS_ttbb"] = new TH1F("hm2rec5OS_ttbb", "M_{#pi#pi} TT/BB OS", massbins, 0, 5.);
    histosTH1F["hm2rec5SS_ttbb"] = new TH1F("hm2rec5SS_ttbb", "M_{#pi#pi} TT/BB SS", massbins, 0, 5.);

    //histosTH1F["hm2rec5OS_diag_trkP"] = new TH1F("hm2rec5OS_diag_trkP","M_{#pi#pi} TB/BT OS, py_{#pi1}py_{#pi2}>0",massbins,0,5.);
    //histosTH1F["hm2rec5OS_diag_trkM"] = new TH1F("hm2rec5OS_diag_trkM","M_{#pi#pi} TB/BT OS, py_{#pi1}py_{#pi2}<0",massbins,0,5.);
    //histosTH1F["hm2rec5OS_ttbb_trkP"] = new TH1F("hm2rec5OS_ttbb_trkP","M_{#pi#pi} TB/BT OS, py_{#pi1}py_{#pi2}>0",massbins,0,5.);
    //histosTH1F["hm2rec5OS_ttbb_trkM"] = new TH1F("hm2rec5OS_ttbb_trkM","M_{#pi#pi} TB/BT OS, py_{#pi1}py_{#pi2}<0",massbins,0,5.);
    //...Luiz
    histosTH1F["hm2rec5OS_diag_trkP"] = new TH1F("hm2rec5OS_diag_trkP", "M_{#pi#pi} TB/BT OS, py_{#pi_{1}}py_{#pi_{2}}>0", massbins, 0, 5.);
    histosTH1F["hm2rec5OS_diag_trkM"] = new TH1F("hm2rec5OS_diag_trkM", "M_{#pi#pi} TB/BT OS, py_{#pi_{1}}py_{#pi_{2}}<0", massbins, 0, 5.);
    histosTH1F["hm2rec5OS_ttbb_trkP"] = new TH1F("hm2rec5OS_ttbb_trkP", "M_{#pi#pi} TB/BT OS, py_{#pi_{1}}py_{#pi_{2}}>0", massbins, 0, 5.);
    histosTH1F["hm2rec5OS_ttbb_trkM"] = new TH1F("hm2rec5OS_ttbb_trkM", "M_{#pi#pi} TB/BT OS, py_{#pi_{1}}py_{#pi_{2}}<0", massbins, 0, 5.);

    histosTH1F["hm2rec5OS_diag_pypxP"] = new TH1F("hm2rec5OS_diag_pypxP", "M_{#pi#pi} TB/BT OS, |py/px|_{#pi#pi} > 1", massbins, 0, 5.);
    histosTH1F["hm2rec5OS_diag_pypxM"] = new TH1F("hm2rec5OS_diag_pypxM", "M_{#pi#pi} TB/BT OS, |py/px|_{#pi#pi} < 1", massbins, 0, 5.);
    histosTH1F["hm2rec5OS_ttbb_pypxP"] = new TH1F("hm2rec5OS_ttbb_pypxP", "M_{#pi#pi} TT/BB OS, |py/px|_{#pi#pi} > 1", massbins, 0, 5.);
    histosTH1F["hm2rec5OS_ttbb_pypxM"] = new TH1F("hm2rec5OS_ttbb_pypxM", "M_{#pi#pi} TT/BB OS, |py/px|_{#pi#pi} < 1", massbins, 0, 5.);

    histosTH1F["hm2rec6OS"] = new TH1F("hm2rec6OS", "M_{#pi#pi} OS", massbins, 0, 5.);
    histosTH1F["hm2rec6SS"] = new TH1F("hm2rec6SS", "M_{#pi#pi} SS", massbins, 0, 5.);
    histosTH1F["hm2rec6OS_diag"] = new TH1F("hm2rec6OS_diag", "M_{#pi#pi} TB/BT OS", massbins, 0, 5.);
    histosTH1F["hm2rec6SS_diag"] = new TH1F("hm2rec6SS_diag", "M_{#pi#pi} TB/BT SS", massbins, 0, 5.);
    histosTH1F["hm2rec6OS_ttbb"] = new TH1F("hm2rec6OS_ttbb", "M_{#pi#pi} TT/BB OS", massbins, 0, 5.);
    histosTH1F["hm2rec6SS_ttbb"] = new TH1F("hm2rec6SS_ttbb", "M_{#pi#pi} TT/BB SS", massbins, 0, 5.);

    //histosTH1F["hm2rec6OS_diag_trkP"] = new TH1F("hm2rec6OS_diag_trkP","M_{#pi#pi} TB/BT OS, py_{#pi1}py_{#pi2}>0",massbins,0,5.);
    //histosTH1F["hm2rec6OS_diag_trkM"] = new TH1F("hm2rec6OS_diag_trkM","M_{#pi#pi} TB/BT OS, py_{#pi1}py_{#pi2}<0",massbins,0,5.);
    //histosTH1F["hm2rec6OS_ttbb_trkP"] = new TH1F("hm2rec6OS_ttbb_trkP","M_{#pi#pi} TT/BB OS, py_{#pi1}py_{#pi2}>0",massbins,0,5.);
    //histosTH1F["hm2rec6OS_ttbb_trkM"] = new TH1F("hm2rec6OS_ttbb_trkM","M_{#pi#pi} TT/BB OS, py_{#pi1}py_{#pi2}<0",massbins,0,5.);
    //...Luiz
    histosTH1F["hm2rec6OS_diag_trkP"] = new TH1F("hm2rec6OS_diag_trkP", "M_{#pi#pi} TB/BT OS, py_{#pi_{1}}py_{#pi_{2}}>0", massbins, 0, 5.);
    histosTH1F["hm2rec6OS_diag_trkM"] = new TH1F("hm2rec6OS_diag_trkM", "M_{#pi#pi} TB/BT OS, py_{#pi_{1}}py_{#pi_{2}}<0", massbins, 0, 5.);
    histosTH1F["hm2rec6OS_ttbb_trkP"] = new TH1F("hm2rec6OS_ttbb_trkP", "M_{#pi#pi} TT/BB OS, py_{#pi_{1}}py_{#pi_{2}}>0", massbins, 0, 5.);
    histosTH1F["hm2rec6OS_ttbb_trkM"] = new TH1F("hm2rec6OS_ttbb_trkM", "M_{#pi#pi} TT/BB OS, py_{#pi_{1}}py_{#pi_{2}}<0", massbins, 0, 5.);

    histosTH1F["hm2recHFvetoOS"] = new TH1F("hm2recHFvetoOS", "M_{#pi#pi} HFv OS", massbins, 0, 5.);
    histosTH1F["hm2recHFvetoSS"] = new TH1F("hm2recHFvetoSS", "M_{#pi#pi} HFv SS", massbins, 0, 5.);

    histosTH1F["hm2rec45OS"] = new TH1F("hm2rec45OS", "M_{#pi#pi} OS", massbins, 0, 5.);
    histosTH1F["hm2rec45SS"] = new TH1F("hm2rec45SS", "M_{#pi#pi} SS", massbins, 0, 5.);
    histosTH1F["hm2rec4515OS"] = new TH1F("hm2rec4515OS", "M_{#pi#pi} OS", massbins, 0, 5.);
    histosTH1F["hm2rec4515SS"] = new TH1F("hm2rec4515SS", "M_{#pi#pi} SS", massbins, 0, 5.);

    histosTH1F["hm2rec9919"] = new TH1F("hm2rec9919", "M_{#pi#pi} 9919", massbins, 0, 5.);
    histosTH1F["hm2rec9922"] = new TH1F("hm2rec9922", "M_{#pi#pi} 9919,9922", massbins, 0, 5.);
    histosTH1F["hm2rec9971"] = new TH1F("hm2rec9971", "M_{#pi#pi} 9971", massbins, 0, 5.);
    histosTH1F["hm2rec9978"] = new TH1F("hm2rec9978", "M_{#pi#pi} 9978", massbins, 0, 5.);

    histosTH1F["hnclusters"] = new TH1F("hnclusters", "nPixelClusters", 500, 0, 500.);
    histosTH1F["hnclusters2"] = new TH1F("hnclusters2", "nStripClusters", 500, 0, 500.);
    histosTH1F["hnclustersOSdiag"] = new TH1F("hnclustersOSdiag", "nPixelClusters", 500, 0, 500.);
    histosTH1F["hnclusters2OSdiag"] = new TH1F("hnclusters2OSdiag", "nStripClusters", 500, 0, 500.);

    histosTH1F["halgo"] = new TH1F("halgo", "Algo", 15, 0, 15.);
    histosTH1F["hnhits"] = new TH1F("hnhits", "nhits pix+strip", 40, 0, 40.);
    histosTH1F["hchi2"] = new TH1F("hchi2", "normalized #chi^{2}", 1050, -50, 1000.);
    //histosTH1F["hdz"] = new TH1F("hdz","dz",1000,-200,200.);
    //histosTH1F["hd0"] = new TH1F("hd0","d0",2000,-200,200.);
    //...Luiz
    histosTH1F["hdz"] = new TH1F("hdz", "dz", 2000, -100, 100.);
    histosTH1F["hd0"] = new TH1F("hd0", "d0", 2000, -100, 100.);

    histosTH1F["halgov"] = new TH1F("halgov", "Algo", 15, 0, 15.);
    histosTH1F["hnhitsv"] = new TH1F("hnhitsv", "nhits pixel", 40, 0, 40.);
    histosTH1F["hchi2v"] = new TH1F("hchi2v", "normalized #chi^{2} vtx-fitted", 550, -50, 500.);
    //histosTH1F["hdzv"] = new TH1F("hdzv","dz vtx-fitted",500,-100,100.);
    //...Luiz
    histosTH1F["hdzv"] = new TH1F("hdzv", "dz vtx-fitted", 1000, -100, 100.);
    histosTH1F["hd0v"] = new TH1F("hd0v", "d0 vtx-fitted", 2000, -20, 20.);

    histosTH1F["hchi2fin"] = new TH1F("hchi2fin", "normalized #chi^{2} vtx-fitted", 550, -50, 500.);
    //histosTH1F["hdzfin"] = new TH1F("hdzfin","dz vtx-fitted",500,-100,100.);
    //...Luiz
    histosTH1F["hdzfin"] = new TH1F("hdzfin", "dz vtx-fitted", 1000, -100, 100.);
    histosTH1F["hd0fin"] = new TH1F("hd0fin", "d0 vtx-fitted", 2000, -20, 20.);

    //histosTH1F["hdeltaR"] = new TH1F("hdeltaR","#Delta R trk-trk",200,0,10.);
    //histosTH1F["hdeltaR2"] = new TH1F("hdeltaR2","#Delta R trk-trk",200,0,10.);
    //...Luiz
    histosTH1F["hdeltaR"] = new TH1F("hdeltaR", "#DeltaR trk-trk", 200, 0, 10.);
    histosTH1F["hdeltaR2"] = new TH1F("hdeltaR2", "#DeltaR trk-trk", 200, 0, 10.);

    //-----------------
    histosTH2F["h2dimdpyAll"] = new TH2F("h2dimdpyAll", "p_{y}^{TOTEM} vs p_{y}^{CMS}", 200, -2., 2., 200, -2., 2.);
    histosTH2F["h2dimdpy"] = new TH2F("h2dimdpy", "p_{y}^{TOTEM} vs p_{y}^{CMS}", 200, -2., 2., 200, -2., 2.);
    histosTH2F["h2dimdpy_diag"] = new TH2F("h2dimdpy_diag", "p_{y}^{TOTEM} vs p_{y}^{CMS} diag", 100, -2., 2., 100, -2., 2.);
    histosTH2F["h2dimdpy_ttbb"] = new TH2F("h2dimdpy_ttbb", "p_{y}^{TOTEM} vs p_{y}^{CMS} TT/BB", 100, -2., 2., 100, -2., 2.);

    //histosTH1F["hdpyAll"] = new TH1F("hdpyAll"  ,"#Delta p_{Y} CMS-TOTEM",500,-0.5,0.5);
    //histosTH1F["hdpy"] = new TH1F("hdpy"     ,"#Delta p_{Y} CMS-TOTEM",500,-0.5,0.5);
    //histosTH1F["hdpy_diag"] = new TH1F("hdpy_diag","#Delta p_{Y} CMS-TOTEM TB/BT",500,-0.5,0.5);
    //histosTH1F["hdpy_ttbb"] = new TH1F("hdpy_ttbb","#Delta p_{Y} CMS-TOTEM TT/BB",500,-0.5,0.5);
    //...Luiz
    histosTH1F["hdpyAll"] = new TH1F("hdpyAll", "#Deltap_{Y} CMS-TOTEM", 500, -0.5, 0.5);
    histosTH1F["hdpy"] = new TH1F("hdpy", "#Deltap_{Y} CMS-TOTEM", 500, -0.5, 0.5);
    histosTH1F["hdpy_diag"] = new TH1F("hdpy_diag", "#Deltap_{Y} CMS-TOTEM TB/BT", 500, -0.5, 0.5);
    histosTH1F["hdpy_ttbb"] = new TH1F("hdpy_ttbb", "#Deltap_{Y} CMS-TOTEM TT/BB", 500, -0.5, 0.5);

    histosTH2F["h2dimdpxAll"] = new TH2F("h2dimdpxAll", "p_{x}^{TOTEM} vs p_{x}^{CMS}", 200, -2., 2., 200, -2., 2.);
    histosTH2F["h2dimdpx"] = new TH2F("h2dimdpx", "p_{x}^{TOTEM} vs p_{x}^{CMS}", 200, -2., 2., 200, -2., 2.);
    histosTH2F["h2dimdpx_diag"] = new TH2F("h2dimdpx_diag", "p_{x}^{TOTEM} vs p_{x}^{CMS} diag", 100, -2., 2., 100, -2., 2.);
    histosTH2F["h2dimdpx_ttbb"] = new TH2F("h2dimdpx_ttbb", "p_{x}^{TOTEM} vs p_{x}^{CMS} TT/BB", 100, -2., 2., 100, -2., 2.);

    //histosTH1F["hdpxAll"] = new TH1F("hdpxAll", "#Delta p_{X} CMS-TOTEM",500,-0.5,0.5);
    //histosTH1F["hdpx"] = new TH1F("hdpx", "#Delta p_{X} CMS-TOTEM",500,-0.5,0.5);
    //histosTH1F["hdpx_diag"] = new TH1F("hdpx_diag", "#Delta p_{X} CMS-TOTEM TB/BT",500,-0.5,0.5);
    //histosTH1F["hdpx_ttbb"] = new TH1F("hdpx_ttbb", "#Delta p_{X} CMS-TOTEM TT/BB",500,-0.5,0.5);
    //...Luiz
    histosTH1F["hdpxAll"] = new TH1F("hdpxAll", "#Deltap_{X} CMS-TOTEM", 500, -0.5, 0.5);
    histosTH1F["hdpx"] = new TH1F("hdpx", "#Deltap_{X} CMS-TOTEM", 500, -0.5, 0.5);
    histosTH1F["hdpx_diag"] = new TH1F("hdpx_diag", "#Deltap_{X} CMS-TOTEM TB/BT", 500, -0.5, 0.5);
    histosTH1F["hdpx_ttbb"] = new TH1F("hdpx_ttbb", "#Deltap_{X} CMS-TOTEM TT/BB", 500, -0.5, 0.5);

    //------------------
    histosTH2F["h2dimxVtxRL"] = new TH2F("h2dimxVtxRL", "xVtxL vs xVtxR (m)", 1000, -0.004, 0.001, 1000, -0.004, 0.001);
    histosTH2F["h2dimxVtxcmsR"] = new TH2F("h2dimxVtxcmsR", "xVtxCMS vs xVtxR (cm)", 300, -0.3, 0.3, 400, -0.3, 0.5);
    histosTH2F["h2dimxVtxcmsL"] = new TH2F("h2dimxVtxcmsL", "xVtxCMS vs xVtxL (cm)", 300, -0.3, 0.3, 400, -0.3, 0.5);
    histosTH2F["h2dimxVtxcmsRL"] = new TH2F("h2dimxVtxcmsRL", "xVtxCMS vs xVtxRL (cm)", 300, -0.3, 0.3, 400, -0.3, 0.5);

    histosTH2F["h2dimxVtxcmsR2"] = new TH2F("h2dimxVtxcmsR2", "xVtxCMS vs xVtxR (cm) (|xVtxL-xVtxR|<3e-5)", 300, -0.3, 0.3, 400, -0.3, 0.5);
    histosTH2F["h2dimxVtxcmsL2"] = new TH2F("h2dimxVtxcmsL2", "xVtxCMS vs xVtxL (cm) (|xVtxL-xVtxR|<3e-5)", 300, -0.3, 0.3, 400, -0.3, 0.5);
    histosTH2F["h2dimxVtxcmsRL2"] = new TH2F("h2dimxVtxcmsRL2", "xVtxCMS vs xVtxRL (cm)", 300, -0.3, 0.3, 400, -0.3, 0.5);

    histosTH2F["h2dimxVtx_zVtx_CT"] = new TH2F("h2dimxVtx_zVtx_CT", "xVtxCMS-xVtxTOTEM vs zVtx (cm)", 300, -20., 20., 400, -0.3, 0.5);
    histosTH2F["h2dimxVtx_zVtx_C"] = new TH2F("h2dimxVtx_zVtx_C", "xVtxCMS vs zVtx (cm)", 300, -20., 20., 400, -0.3, 0.5);
    histosTH2F["h2dimxVtx_zVtx_T"] = new TH2F("h2dimxVtx_zVtx_T", "xVtxTOTEM vs zVtx (cm)", 300, -20., 20., 400, -0.3, 0.5);

    histosTH1F["hxVtxRL"] = new TH1F("hxVtxRL", "xVtxR-xVtxL (m)", 300, -0.0003, 0.0003);
    //histosTH1F["hxVtxcmsR"] = new TH1F("hxVtxcmsR","xVtxCMS-xVtxR (cm)",300,-0.5,0.5);
    //histosTH1F["hxVtxcmsL"] = new TH1F("hxVtxcmsL","xVtxCMS-xVtxL (cm)",300,-0.5,0.5);
    //histosTH1F["hxVtxcmsRL"] = new TH1F("hxVtxcmsRL","xVtxCMS-xVtxTOTEM (cm)",300,-0.5,0.5);
    //...Luiz
    histosTH1F["hxVtxcmsR"] = new TH1F("hxVtxcmsR", "xVtxCMS-xVtxR (cm)", 500, -0.5, 0.5);
    histosTH1F["hxVtxcmsL"] = new TH1F("hxVtxcmsL", "xVtxCMS-xVtxL (cm)", 500, -0.5, 0.5);
    histosTH1F["hxVtxcmsRL"] = new TH1F("hxVtxcmsRL", "xVtxCMS-xVtxTOTEM (cm)", 500, -0.5, 0.5);

    histosTH1F["hxVtxRL_diag"] = new TH1F("hxVtxRL_diag", "xVtxR-xVtxL (m)", 300, -0.0003, 0.0003);
    //histosTH1F["hxVtxcmsR_diag"] = new TH1F("hxVtxcmsR_diag","xVtxCMS-xVtxR (cm)",300,-0.5,0.5);
    //histosTH1F["hxVtxcmsL_diag"] = new TH1F("hxVtxcmsL_diag","xVtxCMS-xVtxL (cm)",300,-0.5,0.5);
    //histosTH1F["hxVtxcmsRL_diag"] = new TH1F("hxVtxcmsRL_diag","xVtxCMS-xVtxTOTEM (cm)",300,-0.5,0.5);
    //...Luiz
    histosTH1F["hxVtxcmsR_diag"] = new TH1F("hxVtxcmsR_diag", "xVtxCMS-xVtxR (cm)", 500, -0.5, 0.5);
    histosTH1F["hxVtxcmsL_diag"] = new TH1F("hxVtxcmsL_diag", "xVtxCMS-xVtxL (cm)", 500, -0.5, 0.5);
    histosTH1F["hxVtxcmsRL_diag"] = new TH1F("hxVtxcmsRL_diag", "xVtxCMS-xVtxTOTEM (cm)", 500, -0.5, 0.5);

    histosTH1F["hxVtxRL_ttbb"] = new TH1F("hxVtxRL_ttbb", "xVtxR-xVtxL (m)", 300, -0.0003, 0.0003);
    //histosTH1F["hxVtxcmsR_ttbb"] = new TH1F("hxVtxcmsR_ttbb","xVtxCMS-xVtxR (cm)",300,-0.5,0.5);
    //histosTH1F["hxVtxcmsL_ttbb"] = new TH1F("hxVtxcmsL_ttbb","xVtxCMS-xVtxL (cm)",300,-0.5,0.5);
    //histosTH1F["hxVtxcmsRL_ttbb"] = new TH1F("hxVtxcmsRL_ttbb","xVtxCMS-xVtxTOTEM (cm)",300,-0.5,0.5);
    //...Luiz
    histosTH1F["hxVtxcmsR_ttbb"] = new TH1F("hxVtxcmsR_ttbb", "xVtxCMS-xVtxR (cm)", 500, -0.5, 0.5);
    histosTH1F["hxVtxcmsL_ttbb"] = new TH1F("hxVtxcmsL_ttbb", "xVtxCMS-xVtxL (cm)", 500, -0.5, 0.5);
    histosTH1F["hxVtxcmsRL_ttbb"] = new TH1F("hxVtxcmsRL_ttbb", "xVtxCMS-xVtxTOTEM (cm)", 500, -0.5, 0.5);

    //  histosTH2F["hdedx"] = new TH2F("hdedx","dE/dx vs p", 300, 0.,5.,500, 0.,100.);
    //histosTH2F["hdedx"] = new TH2F("hdedx","dE/dx vs p", 300, 0.,5.,1000, 0.,200.);
    //...Luiz
    histosTH2F["hdedx"] = new TH2F("hdedx", "dE/dx vs p", 500, 0., 5., 1000, 0., 200.);
    //...Luiz
    histosTH2F["hlndedx"]  = new TH2F("hlndedx", "ln dE/dx vs p", 500, 0., 5., 1000, 0., 5.);
    histosTH2F["hl10dedx"] = new TH2F("hl10dedx", "log10 dE/dx vs p", 500, 0., 5., 1000, 0., 5.);

    //---------------------------------------------------

    for(map<string, TH1F *>::const_iterator it = histosTH1F.begin(); it != histosTH1F.end(); ++it)
        it->second->Sumw2();
    for(map<string, TH2F *>::const_iterator it = histosTH2F.begin(); it != histosTH2F.end(); ++it)
        it->second->Sumw2();

    //===================

    //vector<TString>* vfiles = new vector<TString>(1,"merged_reduced_8372_198903_LP_Jets1_1_test_v1.root");
    vector<TString> *vfiles = new vector<TString>;
    for(size_t idx_file = 0; idx_file < fileNames.size(); ++idx_file) vfiles->push_back( fileNames[idx_file] );

    // Declaration of tree and its branches variables
    TTree *tree = NULL;
    MyEvtId           *evtId       = NULL;
    vector<MyCaloTower> *calo_coll = NULL;
    vector<MyTracks>   *track_coll = NULL;
    vector<MyVertex>  *vertex_coll = NULL;
    vector<MySiPixelCluster> *sipixelcluster_coll = NULL;
    MySiStripCluster         *sistripcluster_coll = NULL;


    RPRootDumpReconstructedProton *rec_proton_left  = NULL;
    RPRootDumpReconstructedProton *rec_proton_right = NULL;
    map<unsigned int, RPRootDumpTrackInfo *> rp_track_info;
    //  TriggerData   *trigData       = NULL;

    //===================

    std::map< int, TMatrix> AlltransportMatrixPlus;
    std::map< int, TMatrix> AlltransportMatrixMinus;

    //XRPV.B6R5.B1      B1 is Right , CMS minus
    TMatrix M220M(6, 6);

    M220M(0, 0) = -1.871247999249703e+00  ;
    M220M(0, 1) = 1.733151135160244e-02  ;
    M220M(0, 2) = 0.000000000000000e+00  ;
    M220M(0, 3) = 0.000000000000000e+00  ;
    M220M(0, 4) = -3.821064474332431e-02 ;
    M220M(0, 5) = -3.821064474332431e-02  ;
    M220M(1, 0) = 5.528023408827136e-02 ;
    M220M(1, 1) = -5.349147148886547e-01  ;
    M220M(1, 2) = 0.000000000000000e+00  ;
    M220M(1, 3) = 0.000000000000000e+00  ;
    M220M(1, 4) = 2.332546482011731e-03  ;
    M220M(1, 5) = 2.332546482011731e-03  ;
    M220M(2, 0) = 0.000000000000000e+00  ;
    M220M(2, 1) = 0.000000000000000e+00 ;
    M220M(2, 2) = -2.321378009782771e-08  ;
    M220M(2, 3) = 2.629525462245173e+02  ;
    M220M(2, 4) = 0.000000000000000e+00  ;
    M220M(2, 5) = 0.000000000000000e+00  ;
    M220M(3, 0) = 0.000000000000000e+00  ;
    M220M(3, 1) = 0.000000000000000e+00 ;
    M220M(3, 2) = -3.802967965874805e-03  ;
    M220M(3, 3) = 4.731545364353734e+00  ;
    M220M(3, 4) = 0.000000000000000e+00  ;
    M220M(3, 5) = 0.000000000000000e+00  ;
    M220M(5, 0) = 2.252479551546639e-03  ;
    M220M(5, 1) = 2.039900958275588e-02  ;
    M220M(5, 2) = 0.000000000000000e+00  ;
    M220M(5, 3) = 0.000000000000000e+00  ;
    M220M(5, 4) = 1.000000000000000e+00  ;
    M220M(5, 5) = 9.584144208086515e-05  ;
    M220M(5, 0) = 0.000000000000000e+00  ;
    M220M(5, 1) = 0.000000000000000e+00  ;
    M220M(5, 2) = 0.000000000000000e+00  ;
    M220M(5, 3) = 0.000000000000000e+00  ;
    M220M(5, 4) = 0.000000000000000e+00  ;
    M220M(5, 5) = 1.000000000000000e+00  ;

    //XRPV.B6L5.B2      B2 is Left, CMS plus
    TMatrix M220P(6, 6);
    M220P(0, 0) = -1.897523818078534e+00  ;
    M220P(0, 1) = 1.062411421653394e-01; //1.062411421653394e-01  ;
    M220P(0, 2) = 0.000000000000000e+00   ;
    M220P(0, 3) = 0.000000000000000e+00   ;
    M220P(0, 4) = 5.198622934357949e-02   ; //to cross check
    M220P(0, 5) = 5.198622934357949e-02   ;
    M220P(1, 0) = 5.401504221523073e-02   ;
    M220P(1, 1) = -5.300268751290215e-01  ;
    M220P(1, 2) = 0.000000000000000e+00    ;
    M220P(1, 3) = 0.000000000000000e+00    ;
    M220P(1, 4) = -2.668640114664157e-03  ;
    M220P(1, 5) = -2.668640114664157e-03   ;
    M220P(2, 0) = 0.000000000000000e+00    ;
    M220P(2, 1) = 0.000000000000000e+00   ;
    M220P(2, 2) = -3.186327537105585e-09    ;
    M220P(2, 3) = 2.618610731959413e+02   ;
    M220P(2, 4) = 0.000000000000000e+00  ;
    M220P(2, 5) = 0.000000000000000e+00   ;
    M220P(3, 0) = 0.000000000000000e+00   ;
    M220P(3, 1) = 0.000000000000000e+00   ;
    M220P(3, 2) = -3.818818897730703e-03   ;
    M220P(3, 3) = 4.676450995853369e+00    ;
    M220P(3, 4) = 0.000000000000000e+00    ;
    M220P(3, 5) = 0.000000000000000e+00   ;
    M220P(5, 0) = -2.255769806850952e-03   ;
    M220P(5, 1) = -2.727057931490794e-02   ;
    M220P(5, 2) = 0.000000000000000e+00    ;
    M220P(5, 3) = 0.000000000000000e+00    ;
    M220P(5, 4) = 1.000000000000000e+00   ;
    M220P(5, 5) = 1.107429910138087e-04   ;
    M220P(5, 0) = 0.000000000000000e+00   ;
    M220P(5, 1) = 0.000000000000000e+00   ;
    M220P(5, 2) = 0.000000000000000e+00   ;
    M220P(5, 3) = 0.000000000000000e+00  ;
    M220P(5, 4) = 0.000000000000000e+00  ;
    M220P(5, 5) = 1.000000000000000e+00    ;


    TMatrix M215P(6, 6);

    M215P(0, 0) = -2.187692624858721e+00    ;
    M215P(0, 1) = 2.953545515358119e+00    ;
    M215P(0, 2) = 0.000000000000000e+00    ;
    M215P(0, 3) = 0.000000000000000e+00    ;
    M215P(0, 4) = 6.603743360296731e-02    ;
    M215P(0, 5) = 6.603743360296731e-02    ;
    M215P(1, 0) = 5.401504221523073e-02   ;
    M215P(1, 1) = -5.300268751290215e-01   ;
    M215P(1, 2) = 0.000000000000000e+00   ;
    M215P(1, 3) = 0.000000000000000e+00    ;
    M215P(1, 4) = -2.668640114664157e-03   ;
    M215P(1, 5) = -2.668640114664157e-03   ;
    M215P(2, 0) = 0.000000000000000e+00   ;
    M215P(2, 1) = 0.000000000000000e+00   ;
    M215P(2, 2) = 2.051469193227947e-02    ;
    M215P(2, 3) = 2.367391784462199e+02   ;
    M215P(2, 4) = 0.000000000000000e+00    ;
    M215P(2, 5) = 0.000000000000000e+00    ;
    M215P(3, 0) = 0.000000000000000e+00    ;
    M215P(3, 1) = 0.000000000000000e+00   ;
    M215P(3, 2) = -3.818818897730703e-03    ;
    M215P(3, 3) = 4.676450995853369e+00    ;
    M215P(3, 4) = 0.000000000000000e+00    ;
    M215P(3, 5) = 0.000000000000000e+00   ;
    M215P(4, 0) = -2.271149533403129e-03   ;
    M215P(4, 1) = -2.711966453134992e-02  ;
    M215P(4, 2) = 0.000000000000000e+00    ;
    M215P(4, 3) = 0.000000000000000e+00    ;
    M215P(4, 4) = 1.000000000000000e+00   ;
    M215P(4, 5) = 1.113908988335683e-04   ;
    M215P(5, 0) = 0.000000000000000e+00   ;
    M215P(5, 1) = 0.000000000000000e+00    ;
    M215P(5, 2) = 0.000000000000000e+00   ;
    M215P(5, 3) = 0.000000000000000e+00   ;
    M215P(5, 4) = 0.000000000000000e+00   ;
    M215P(5, 5) = 1.000000000000000e+00  ;

    TMatrix M215M(6, 6);

    M215M(0, 0) = -2.168213416771863e+00  ;
    M215M(0, 1) = 2.890893359733129e+00   ;
    M215M(0, 2) = 0.000000000000000e+00   ;
    M215M(0, 3) = 0.000000000000000e+00   ;
    M215M(0, 4) = -5.074108445998152e-02  ;
    M215M(0, 5) = -5.074108445998152e-02   ;
    M215M(1, 0) = 5.528023408827136e-02  ;
    M215M(1, 1) = -5.349147148886547e-01   ;
    M215M(1, 2) = 0.000000000000000e+00   ;
    M215M(1, 3) = 0.000000000000000e+00   ;
    M215M(1, 4) = 2.332546482011731e-03   ;
    M215M(1, 5) = 2.332546482011731e-03   ;
    M215M(2, 0) = 0.000000000000000e+00   ;
    M215M(2, 1) = 0.000000000000000e+00   ;
    M215M(2, 2) = 2.042952069889703e-02   ;
    M215M(2, 3) = 2.375346845272119e+02   ;
    M215M(2, 4) = 0.000000000000000e+00   ;
    M215M(2, 5) = 0.000000000000000e+00   ;
    M215M(3, 0) = 0.000000000000000e+00   ;
    M215M(3, 1) = 0.000000000000000e+00  ;
    M215M(3, 2) = -3.802967965874805e-03  ;
    M215M(3, 3) = 4.731545364353734e+00   ;
    M215M(3, 4) = 0.000000000000000e+00   ;
    M215M(3, 5) = 0.000000000000000e+00   ;
    M215M(4, 0) = 2.252479550701315e-03   ;
    M215M(4, 1) = 2.039900959093559e-02  ;
    M215M(4, 2) = 0.000000000000000e+00   ;
    M215M(4, 3) = 0.000000000000000e+00   ;
    M215M(4, 4) = 1.000000000000000e+00   ;
    M215M(4, 5) = 9.572950680001605e-05   ;
    M215M(5, 0) = 0.000000000000000e+00   ;
    M215M(5, 1) = 0.000000000000000e+00   ;
    M215M(5, 2) = 0.000000000000000e+00   ;
    M215M(5, 3) = 0.000000000000000e+00   ;
    M215M(5, 4) = 0.000000000000000e+00   ;
    M215M(5, 5) = 1.000000000000000e+00   ;

    AlltransportMatrixPlus.insert(std::make_pair(220, M220P));
    AlltransportMatrixMinus.insert(std::make_pair(220, M220M));

    AlltransportMatrixPlus.insert(std::make_pair(215, M215P));
    AlltransportMatrixMinus.insert(std::make_pair(215, M215M));

    //===================

    int i_tot = 0, nevt_tot = 0;
    //starting Loop over files, stops at end of list of files or when reached nevt_max
    for(vector<TString>::iterator itfiles = vfiles->begin() ; itfiles != vfiles->end() && i_tot < nevt_max_corr ; ++itfiles)
    {

        cout << "Opening file " << *itfiles << endl;

        TFile *file = TFile::Open(*itfiles, "READ");
        if (!file || file->IsZombie())
        {
            cout << "corrupted file - skipping " << endl;
            continue;
        }

        // Access TTree from current file
        tree = (TTree *) file->Get( treeName.c_str() );
        int nev = int(tree->GetEntriesFast());
        nevt_tot += nev;
        //RC
        cout << nev << " entries in " << *itfiles << endl;

        // Add branches to TTree ----------------------------------------------------------------------
        tree->SetBranchAddress("cmsEvtUA", &evtId);
        tree->SetBranchAddress("cmsCaloTowersUA", &calo_coll);
        // tracks
        //    tree->SetBranchAddress("cmsTracksUA",&track_coll);//generalTracks
        tree->SetBranchAddress("cmsTracksPIDUA", &track_coll); // refittedTracks

        tree->SetBranchAddress("cmsVerticesUA", &vertex_coll);
        tree->SetBranchAddress("SiPixelClusters", &sipixelcluster_coll);
        tree->SetBranchAddress("SiStripClusters", &sistripcluster_coll);

        //    tree->SetBranchAddress("trigger_data.",&trigData);
        tree->SetBranchAddress("rec_prot_left.", &rec_proton_left);
        tree->SetBranchAddress("rec_prot_right.", &rec_proton_right);
        std::vector<unsigned int> rp_list;
        rp_list.push_back(20);
        rp_list.push_back(21);
        rp_list.push_back(24);
        rp_list.push_back(25);
        rp_list.push_back(120);
        rp_list.push_back(121);
        rp_list.push_back(124);
        rp_list.push_back(125);
        char br_name[200];
        for (unsigned int a = 0; a < 2; ++a)
        {
            int s = 2;
            for (unsigned int r = 0; r < 6; r++)
            {
                unsigned int id = 100 * a + 10 * s + r;
                if( std::find(rp_list.begin(), rp_list.end(), id) == rp_list.end() ) continue;

                sprintf(br_name, "track_rp_%u.", id);
                //RC
                //	std::cout << br_name << std::endl;
                tree->SetBranchAddress(br_name, &rp_track_info[id]);
            }
        }

        //starting loop over events, stops when reached end of file or nevt_max
        for(int i_evt = 0; i_evt < nev && i_tot < nevt_max_corr; ++i_evt, ++i_tot)
        {

            if( ((i_tot + 1) % 5000) == 0) cout << int(double(i_tot + 1) / 1000) << "k done" << endl;
            tree->GetEntry(i_evt);

            //-------------------------------------------------------------------------------------------------
            // TOTEM RP protons

            histosTH1F["EventSelection"]->Fill( "TOTEM0", wei );

            bool proton_right_valid = rec_proton_right->valid;
            bool proton_left_valid = rec_proton_left->valid;

            if(!(proton_right_valid && proton_left_valid) ) continue;

            histosTH1F["EventSelection"]->Fill( "2valid", wei );


            //--------------------------------------------
            // fiducial cut

            RPRootDumpTrackInfo *rp_020 = rp_track_info[20];
            RPRootDumpTrackInfo *rp_021 = rp_track_info[21];
            RPRootDumpTrackInfo *rp_024 = rp_track_info[24];
            RPRootDumpTrackInfo *rp_025 = rp_track_info[25];

            RPRootDumpTrackInfo *rp_120 = rp_track_info[120];
            RPRootDumpTrackInfo *rp_121 = rp_track_info[121];
            RPRootDumpTrackInfo *rp_124 = rp_track_info[124];
            RPRootDumpTrackInfo *rp_125 = rp_track_info[125];

            bool rp_valid_020 = rp_020->valid;
            bool rp_valid_021 = rp_021->valid;
            bool rp_valid_024 = rp_024->valid;
            bool rp_valid_025 = rp_025->valid;

            bool rp_valid_120 = rp_120->valid;
            bool rp_valid_121 = rp_121->valid;
            bool rp_valid_124 = rp_124->valid;
            bool rp_valid_125 = rp_125->valid;

            //------------------------------------------------
            // -z                    IP                    +z
            //          sec45                   sec56
            //top:  024       020           120       124
            //ver:     023 022                 122 123
            //bot:  025       021           121       125
            //
            //------------------------------------------------

            bool diag_top45_bot56 = rp_valid_020 && rp_valid_024 && rp_valid_121 && rp_valid_125;
            bool diag_bot45_top56 = rp_valid_021 && rp_valid_025 && rp_valid_120 && rp_valid_124;

            bool top45_top56      = rp_valid_020 && rp_valid_024 && rp_valid_120 && rp_valid_124;
            bool bot45_bot56      = rp_valid_021 && rp_valid_025 && rp_valid_121 && rp_valid_125;

            int nconf = 0;
            if(diag_top45_bot56) nconf++;
            if(diag_bot45_top56) nconf++;
            if(top45_top56) nconf++;
            if(bot45_bot56) nconf++;

            //      if(diag_top45_bot56 || diag_bot45_top56 || top45_top56 || bot45_bot56);
            //      else continue;

            if(nconf == 0) continue;

            histosTH1F["EventSelection"]->Fill( "anyTB/BT/TT/BB", wei );
            histosTH1F["hnconf"]->Fill(nconf, wei );

            if(nconf != 1) continue;

            histosTH1F["EventSelection"]->Fill( "exclusiveTB/BT/TT/BB", wei );

            bool fiducialCutTB = true;
            if(diag_top45_bot56)
            {

                double x_020 = rp_020->x;
                double y_020 = rp_020->y;
                histosTH1F["rp_x_020"]->Fill( x_020, wei );
                histosTH1F["rp_y_020"]->Fill( y_020, wei );
                histosTH2F["rp_yx_020"]->Fill( x_020, y_020, wei );

                double x_024 = rp_024->x;
                double y_024 = rp_024->y;
                histosTH1F["rp_x_024"]->Fill( x_024, wei );
                histosTH1F["rp_y_024"]->Fill( y_024, wei );
                histosTH2F["rp_yx_024"]->Fill( x_024, y_024, wei );

                double x_121 = rp_121->x;
                double y_121 = rp_121->y;
                histosTH1F["rp_x_121"]->Fill( x_121, wei );
                histosTH1F["rp_y_121"]->Fill( y_121, wei );
                histosTH2F["rp_yx_121"]->Fill( x_121, y_121, wei );

                double x_125 = rp_125->x;
                double y_125 = rp_125->y;
                histosTH1F["rp_x_125"]->Fill( x_125, wei );
                histosTH1F["rp_y_125"]->Fill( y_125, wei );
                histosTH2F["rp_yx_125"]->Fill( x_125, y_125, wei );

                if( x_020 < -1.5 ) fiducialCutTB = false;
                if( x_024 < -1.5 ) fiducialCutTB = false;
                if( x_121 < -1.5 ) fiducialCutTB = false;
                if( x_125 < -1.5 ) fiducialCutTB = false;

                if( y_020 < 6.0 || y_020 > 26.0) fiducialCutTB = false;
                if( y_024 < 6.7 || y_024 > 28.7) fiducialCutTB = false;
                if( y_121 < -25.8 || y_121 > -6.4) fiducialCutTB = false;
                if( y_125 < -28.6 || y_125 > -7.1) fiducialCutTB = false;

            }

            bool fiducialCutBT = true;
            if(diag_bot45_top56)
            {

                double x_021 = rp_021->x;
                double y_021 = rp_021->y;

                histosTH1F["rp_x_021"]->Fill( x_021, wei );
                histosTH1F["rp_y_021"]->Fill( y_021, wei );
                histosTH2F["rp_yx_021"]->Fill( x_021, y_021, wei );

                double x_025 = rp_025->x;
                double y_025 = rp_025->y;
                histosTH1F["rp_x_025"]->Fill( x_025, wei );
                histosTH1F["rp_y_025"]->Fill( y_025, wei );
                histosTH2F["rp_yx_025"]->Fill( x_025, y_025, wei );

                double x_120 = rp_120->x;
                double y_120 = rp_120->y;
                histosTH1F["rp_x_120"]->Fill( x_120, wei );
                histosTH1F["rp_y_120"]->Fill( y_120, wei );
                histosTH2F["rp_yx_120"]->Fill( x_120, y_120, wei );

                double x_124 = rp_124->x;
                double y_124 = rp_124->y;
                histosTH1F["rp_x_124"]->Fill( x_124, wei );
                histosTH1F["rp_y_124"]->Fill( y_124, wei );
                histosTH2F["rp_yx_124"]->Fill( x_124, y_124, wei );

                if(x_021 < -1.5) fiducialCutBT = false;
                if(x_025 < -1.5) fiducialCutBT = false;
                if(x_120 < -1.5) fiducialCutBT = false;
                if(x_124 < -1.5) fiducialCutBT = false;

                if( y_021 < -26.3 || y_021 > -6.4) fiducialCutBT = false;
                if( y_025 < -29.0 || y_025 > -7.0) fiducialCutBT = false;
                if( y_120 < 7.7 || y_120 > 24.3) fiducialCutBT = false;
                if( y_124 < 8.5 || y_124 > 26.8) fiducialCutBT = false;

            }

            bool fiducialCutTT = true;
            if(top45_top56)
            {

                double x_020 = rp_020->x;
                double y_020 = rp_020->y;
                histosTH1F["rp2_x_020"]->Fill( x_020, wei );
                histosTH1F["rp2_y_020"]->Fill( y_020, wei );
                histosTH2F["rp2_yx_020"]->Fill( x_020, y_020, wei );

                double x_024 = rp_024->x;
                double y_024 = rp_024->y;
                histosTH1F["rp2_x_024"]->Fill( x_024, wei );
                histosTH1F["rp2_y_024"]->Fill( y_024, wei );
                histosTH2F["rp2_yx_024"]->Fill( x_024, y_024, wei );

                double x_120 = rp_120->x;
                double y_120 = rp_120->y;
                histosTH1F["rp2_x_120"]->Fill( x_120, wei );
                histosTH1F["rp2_y_120"]->Fill( y_120, wei );
                histosTH2F["rp2_yx_120"]->Fill( x_120, y_120, wei );

                double x_124 = rp_124->x;
                double y_124 = rp_124->y;
                histosTH1F["rp2_x_124"]->Fill( x_124, wei );
                histosTH1F["rp2_y_124"]->Fill( y_124, wei );
                histosTH2F["rp2_yx_124"]->Fill( x_124, y_124, wei );

                if(x_020 < -1.5 ) fiducialCutTT = false;
                if(x_024 < -1.5 ) fiducialCutTT = false;
                if(x_120 < -1.5) fiducialCutTT = false;
                if(x_124 < -1.5) fiducialCutTT = false;

                if( y_020 < 6.0 || y_020 > 26.0) fiducialCutTT = false;
                if( y_024 < 6.7 || y_024 > 28.7) fiducialCutTT = false;
                if( y_120 < 7.7 || y_120 > 24.3) fiducialCutTT = false;
                if( y_124 < 8.5 || y_124 > 26.8) fiducialCutTT = false;
            }

            bool fiducialCutBB = true;
            if(bot45_bot56)
            {

                double x_021 = rp_021->x;
                double y_021 = rp_021->y;

                histosTH1F["rp2_x_021"]->Fill( x_021, wei );
                histosTH1F["rp2_y_021"]->Fill( y_021, wei );
                histosTH2F["rp2_yx_021"]->Fill( x_021, y_021, wei );

                double x_025 = rp_025->x;
                double y_025 = rp_025->y;
                histosTH1F["rp2_x_025"]->Fill( x_025, wei );
                histosTH1F["rp2_y_025"]->Fill( y_025, wei );
                histosTH2F["rp2_yx_025"]->Fill( x_025, y_025, wei );

                double x_121 = rp_121->x;
                double y_121 = rp_121->y;
                histosTH1F["rp2_x_121"]->Fill( x_121, wei );
                histosTH1F["rp2_y_121"]->Fill( y_121, wei );
                histosTH2F["rp2_yx_121"]->Fill( x_121, y_121, wei );

                double x_125 = rp_125->x;
                double y_125 = rp_125->y;
                histosTH1F["rp2_x_125"]->Fill( x_125, wei );
                histosTH1F["rp2_y_125"]->Fill( y_125, wei );
                histosTH2F["rp2_yx_125"]->Fill( x_125, y_125, wei );

                if(x_021 < -1.5) fiducialCutBB = false;
                if(x_025 < -1.5) fiducialCutBB = false;
                if(x_121 < -1.5) fiducialCutBB = false;
                if(x_125 < -1.5) fiducialCutBB = false;

                if( y_021 < -26.3 || y_021 > -6.4) fiducialCutBB = false;
                if( y_025 < -29.0 || y_025 > -7.0) fiducialCutBB = false;
                if( y_121 < -25.8 || y_121 > -6.4) fiducialCutBB = false;
                if( y_125 < -28.6 || y_125 > -7.1) fiducialCutBB = false;
            }

            int nfidu = 0;
            if(diag_top45_bot56 && fiducialCutTB) nfidu++;
            if(diag_bot45_top56 && fiducialCutBT) nfidu++;
            if(     top45_top56 && fiducialCutTT) nfidu++;
            if(     bot45_bot56 && fiducialCutBB) nfidu++;

            if(nfidu == 0) continue;

            histosTH1F["EventSelection"]->Fill( "fiducialXY", wei );


            //---------------------------------------------
            // here xVtxL and xVtxR, and thxL and thyR
            // elastic approximation

            double ThxR, ThyR, ThxL, ThyL, xVtxL, xVtxR;

            //bool diag_top45_bot56 = rp_valid_020 && rp_valid_024 && rp_valid_121 && rp_valid_125;
            if(diag_top45_bot56) LikeElastic_ThetaLeftThetaRight220FAR(20, 24, 121, 125, rp_track_info, rp_list,
                        AlltransportMatrixPlus, AlltransportMatrixMinus,
                        ThxR, ThyR, ThxL, ThyL, xVtxL, xVtxR) ;

            //bool diag_bot45_top56 = rp_valid_021 && rp_valid_025 && rp_valid_120 && rp_valid_124;
            if(diag_bot45_top56) LikeElastic_ThetaLeftThetaRight220FAR(21, 25, 120, 124, rp_track_info, rp_list,
                        AlltransportMatrixPlus, AlltransportMatrixMinus,
                        ThxR, ThyR, ThxL, ThyL, xVtxL, xVtxR) ;

            //bool top45_top56      = rp_valid_020 && rp_valid_024 && rp_valid_120 && rp_valid_124;
            if(top45_top56) LikeElastic_ThetaLeftThetaRight220FAR(20, 24, 120, 124, rp_track_info, rp_list,
                        AlltransportMatrixPlus, AlltransportMatrixMinus,
                        ThxR, ThyR, ThxL, ThyL, xVtxL, xVtxR) ;

            //bool bot45_bot56      = rp_valid_021 && rp_valid_025 && rp_valid_121 && rp_valid_125;
            if(bot45_bot56) LikeElastic_ThetaLeftThetaRight220FAR(21, 25, 121, 125, rp_track_info, rp_list,
                        AlltransportMatrixPlus, AlltransportMatrixMinus,
                        ThxR, ThyR, ThxL, ThyL, xVtxL, xVtxR) ;


            //notElastic
            // this is average theta_x and thetha_y, both measure the same thing
            // one is negative one is positive, so minus is needed
            //      double thX=0.5*(thx_proton_left-thx_proton_right);
            //      double thY=0.5*(thy_proton_left-thy_proton_right);
            // not needed here

            // diagonal in thxL vs thxR plane, and
            // diagonal in thyL vs thyR plane

            //      histosTH1F["thyEla"]->Fill(thy_proton_left+thy_proton_right, wei);
            //      histosTH1F["thxEla"]->Fill(thx_proton_left+thx_proton_right, wei);
            histosTH1F["thyEla"]->Fill(ThyL + ThyR, wei);
            histosTH1F["thxEla"]->Fill(ThxL + ThxR, wei);

            if(diag_top45_bot56 || diag_bot45_top56)
            {
                //	histosTH1F["thyEla_diag"]->Fill(thy_proton_left+thy_proton_right, wei);
                //	histosTH1F["thxEla_diag"]->Fill(thx_proton_left+thx_proton_right, wei);
                histosTH1F["thyEla_diag"]->Fill(ThyL + ThyR, wei);
                histosTH1F["thxEla_diag"]->Fill(ThxL + ThxR, wei);
            }
            else
            {
                //	histosTH1F["thyEla_ttbb"]->Fill(thy_proton_left+thy_proton_right, wei);
                //	histosTH1F["thxEla_ttbb"]->Fill(thx_proton_left+thx_proton_right, wei);
                histosTH1F["thyEla_ttbb"]->Fill(ThyL + ThyR, wei);
                histosTH1F["thxEla_ttbb"]->Fill(ThxL + ThxR, wei);
            }

            bool isElastic = false;
            //      if(TMath::Abs(thy_proton_left+thy_proton_right)< 8e-6 &&
            //	 TMath::Abs(thx_proton_left+thx_proton_right)<30e-6) isElastic=true;

            if(TMath::Abs(ThyL + ThyR) < 8e-6 &&
                    TMath::Abs(ThxL + ThxR) < 30e-6) isElastic = true;

            if(isElastic) continue;

            histosTH1F["EventSelection"]->Fill( "notElastic", wei );

            //---------------------------------------------
            //xi selection

            double xi_proton_right = rec_proton_right->xi;
            double t_proton_right = rec_proton_right->t;

            double xi_proton_left = rec_proton_left->xi;
            double t_proton_left = rec_proton_left->t;

            //...Luiz
            double phi_proton_right = rec_proton_right->phi;
            double phi_proton_left = rec_proton_left->phi;
            //

            //-----------------------------------
            //from now on xi - positive
            xi_proton_right = -xi_proton_right;
            xi_proton_left = -xi_proton_left;

            histosTH1F["proton_right_xi"]->Fill( xi_proton_right, wei );
            histosTH1F["proton_left_xi"]->Fill( xi_proton_left, wei );

            // Mx_max=130 GeV
            //      bool proton_right_good = xi_proton_right < 0.01;
            //      bool proton_left_good = xi_proton_left < 0.01;

            // Mx_max=1300 GeV
            bool proton_right_good = xi_proton_right < 0.1;
            bool proton_left_good = xi_proton_left < 0.1;

            // Mx_max=2600 GeV, could do but didn't do
            //      bool proton_right_good = xi_proton_right < 0.2;
            //      bool proton_left_good = xi_proton_left < 0.2;

            if(proton_right_good && proton_left_good);
            else continue;

            histosTH1F["EventSelection"]->Fill( "#xi<0.1", wei );

            histosTH1F["proton_right_logXi"]->Fill( log10(xi_proton_right), wei );
            histosTH1F["proton_left_logXi"]->Fill( log10(xi_proton_left), wei );

            histosTH1F["proton_right_t"]->Fill( -t_proton_right, wei );
            histosTH1F["proton_left_t"]->Fill( -t_proton_left, wei );

            //...Luiz
            histosTH2F["phi_proton_right_t"]->Fill( -t_proton_right, phi_proton_right );
            histosTH2F["phi_proton_left_t"]->Fill( -t_proton_left, phi_proton_left );
            //

            if(diag_top45_bot56 || diag_bot45_top56)
            {
                histosTH1F["proton_right_t_diag"]->Fill( -t_proton_right, wei );
                histosTH1F["proton_left_t_diag"]->Fill( -t_proton_left, wei );
                //...Luiz
                histosTH2F["phi_proton_right_t_diag"]->Fill( -t_proton_right, phi_proton_right );
                histosTH2F["phi_proton_left_t_diag"]->Fill( -t_proton_left, phi_proton_left );
                //
            }
            else
            {
                histosTH1F["proton_right_t_ttbb"]->Fill( -t_proton_right, wei );
                histosTH1F["proton_left_t_ttbb"]->Fill( -t_proton_left, wei );
                //...Luiz
                histosTH2F["phi_proton_right_t_ttbb"]->Fill( -t_proton_right, phi_proton_right );
                histosTH2F["phi_proton_left_t_ttbb"]->Fill( -t_proton_left, phi_proton_left );
                //
            }
            //...Luiz
            if(top45_top56)
            {
                histosTH2F["phi_proton_right_t_tt"]->Fill( -t_proton_right, phi_proton_right );
                histosTH2F["phi_proton_left_t_tt"]->Fill( -t_proton_left, phi_proton_left );
            }
            //...Luiz
            if(bot45_bot56)
            {
                histosTH2F["phi_proton_right_t_bb"]->Fill( -t_proton_right, phi_proton_right );
                histosTH2F["phi_proton_left_t_bb"]->Fill( -t_proton_left, phi_proton_left );
            }

            histosTH1F["proton_dx0"]->Fill(xVtxL - xVtxR);
            histosTH2F["proton_x0_RvsL"]->Fill(xVtxL, xVtxR);

            // HF veto
            /*
            int nHF = 0;
            for(vector<MyCaloTower>::iterator it_ct = calo_coll->begin() ; it_ct != calo_coll->end() ; ++it_ct){

            if(it_ct->hasHF){
              double eHF = it_ct->emEnergy + it_ct->hadEnergy;
              histosTH1F["eHF"]->Fill( eHF , wei );
              if(eHF>5.) nHF++;
            }
                 }
                 histosTH1F["nHF"]->Fill( nHF , wei );
                 */


            //comment if not writing to txt file
            //      int HFveto = 0;
            //      if(nHF>0) HFveto = 1;

            //-------------------
            // After selection
            //-------------------
            //

            int run = evtId->Run;
            int evt = evtId->Evt;
            int LS = evtId->LumiSect;
            //      int runTOTEM = trigData->run_num;
            //      int evtTOTEM = trigData->event_num;

            // double mx_TOTEM=13000.*TMath::Sqrt(*xi_proton_left*xi_proton_right);
            // MX_max=13000.*xi_max;
            // xi<0.01 -> m< 130
            // xi<0.1  -> m<1300
            // not only to 0.1 or 0.01 but with the vertex cut it should be limited to ~10e-4
            // 10e-4 = 1e-3
            // xi<0.001 -> m<13
            //      double pyTOTEM= 6500.*(thy_proton_left+thy_proton_right);
            double TOTEMpy = 6500.*(ThyL + ThyR);
            //For p_x it is more delicate and at the moment we do it for low-xi protons only (|xi| < 3. * 0.006). We first reconstruct th_x as for
            //elastic scattering (see attachment) and then do the sum as above:  3*0.006 = 0.018
            double TOTEMpx = -6500 * (ThxL + ThxR);

            histosTH1F["totem_py"]->Fill(TOTEMpy, wei);
            histosTH1F["totem_px"]->Fill(TOTEMpx, wei);

            //...Luiz
            histosTH1F["totem_pxx"]->Fill(TOTEMpx, wei);
            histosTH1F["totem_pyy"]->Fill(TOTEMpy, wei);

            int tb = 0;
            if(diag_top45_bot56) tb = 1;
            if(diag_bot45_top56) tb = 2;
            if(top45_top56) tb = 3;
            if(bot45_bot56) tb = 4;
            int Topol = tb;

            histosTH1F["hLS"]->Fill(LS, wei);
            //      histosTH1F["htopo"]->Fill(Topol, wei);

            bool diag = false;
            if(Topol == 1 || Topol == 2) diag = true;

            //     if(diag){
            //       histosTH1F["hthyEla2_diag"]->Fill(ThyL+ThyR, wei);
            //       histosTH1F["hthxEla2_diag"]->Fill(ThxL+ThxR, wei);
            //     }else{
            //       histosTH1F["hthyEla2_ttbb"]->Fill(ThyL+ThyR, wei);
            //       histosTH1F["hthxEla2_ttbb"]->Fill(ThxL+ThxR, wei);
            //     }

            // old
            //      fout<<run<<" "<<ls<<" "<<evt<<" "<<tb<<" "<<xi_proton_left<<" "<<xi_proton_right<<" "<<pxTOTEM<<" "<<pyTOTEM<<" "<<xVtxL<<" "<<xVtxR<<" "<<HFveto<<endl;
            // latest
            //      fout<<run<<" "<<ls<<" "<<evt<<" "<<tb<<" "<<xi_proton_left<<" "<<xi_proton_right<<" "<<ThxL<<" "<<ThxR<<" "<<ThyL<<" "<<ThyR<<" "<<xVtxL<<" "<<xVtxR<<" "<<HFveto<<endl;


            //-----------------------------------
            ///////////////////////////////////////////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////////////////

            //...Luiz #particle masses, gev/c^2
            double m_pi = 0.13957;
            double m_k = 0.493667;
            //double m_mu = 0.1056583715;
            double m_mu = 0.1056583745;
            double m_e = 0.0005109989461;

            //-------------------------------
            //accept only 9919, 9922
            //   if(run==259237 && ( (LS>=78 && LS <=100) || (LS>=432 && LS <=576) ) );
            //   else continue;
            //-------------------------------
            // remove 9920
            //   if(run==259237 && (LS>=103 && LS <=423)) continue;
            // remove 9940,9950
            //   if(run==259352 && (LS>=6 && LS <=153)) continue;
            //   if(run==259352 && (LS>=248 && LS <=283)) continue;
            // remove 9976
            //   if(run==259388 && (LS>=59 && LS <=360)) continue;
            // remove 9985
            //   if(run==259399 && (LS>=362 && LS <=387)) continue;
            // remove 9998
            //   if(run==259431 && (LS>=43 && LS <=354)) continue;
            //-------------------------------------------------

            double xiL = xi_proton_left;
            //...Luiz
            //       double xiR = xi_proton_left;
            double xiR = xi_proton_right;
            //   int Topol  = totemTopol[itotem];
            //   double ThyL    = totemThyL[itotem];
            //   double ThyR    = totemThyR[itotem];
            //   double ThxL    = totemThxL[itotem];
            //   double ThxR    = totemThxR[itotem];
            //   double xVtxL   = totemxVtxL[itotem];
            //   double xVtxR   = totemxVtxR[itotem];
            //   double TOTEMpx =-6500.*(ThxL+ThxR);
            //   double TOTEMpy = 6500.*(ThyL+ThyR);
            //   bool HFveto = true; // no activity in HF
            //   if(totemHFveto[itotem]>0) HFveto = false;

            //Topol
            //1 - TB, 2 - BT
            //3 - TT, 4 - BB

            //       bool diag=false;
            //       if(Topol==1 || Topol==2) diag = true;

            //   htopo->Fill(Topol);

            //---------------------------------------------------
            //   hthyEla->Fill(ThyL+ThyR);
            //   hthxEla->Fill(ThxL+ThxR);
            if(diag)
            {
                histosTH1F["hthyEla_diag"]->Fill(ThyL + ThyR, wei);
                histosTH1F["hthxEla_diag"]->Fill(ThxL + ThxR, wei);
            }
            else
            {
                histosTH1F["hthyEla_ttbb"]->Fill(ThyL + ThyR, wei);
                histosTH1F["hthxEla_ttbb"]->Fill(ThxL + ThxR, wei);
            }

            //---------------------------------------------------
            // tighter elastic rejection

            bool isElastic2 = false;

            //      if(TMath::Abs(ThyL+ThyR)< 8e-6 &&
            //	 TMath::Abs(ThxL+ThxR)<30e-6) isElastic=true;

            if(TMath::Abs(ThyL + ThyR) < 15e-6 &&
                    TMath::Abs(ThxL + ThxR) < 45e-6) isElastic2 = true;

            if(isElastic2) continue;

            //---------------------------------------------------
            if(diag)
            {
                histosTH1F["hthyEla2_diag"]->Fill(ThyL + ThyR, wei);
                histosTH1F["hthxEla2_diag"]->Fill(ThxL + ThxR, wei);
            }
            else
            {
                histosTH1F["hthyEla2_ttbb"]->Fill(ThyL + ThyR, wei);
                histosTH1F["hthxEla2_ttbb"]->Fill(ThxL + ThxR, wei);
            }

            //---------------------------------------------------
            // after new anti-elastic cut
            histosTH1F["htopo"]->Fill(Topol, wei);

            //---------------------------------------------------

            bool fiducialRegion = false;
            double etaCut = 2.5;
            bool fiducialRegionPt = false;
            //double ptCut= 0.2;
            //...Luiz
            double ptCut = 0.1;

            //tracks in 4track-events (npixelhits>0)
            //...Luiz
            TLorentzVector k1(0., 0., 0., 0.);
            TLorentzVector k2(0., 0., 0., 0.);
            TLorentzVector k3(0., 0., 0., 0.);
            TLorentzVector k4(0., 0., 0., 0.);
            TLorentzVector k1k2Rec(0., 0., 0., 0.);
            TLorentzVector k3k4Rec(0., 0., 0., 0.);
            TLorentzVector k1k3Rec(0., 0., 0., 0.);
            TLorentzVector k2k4Rec(0., 0., 0., 0.);

            //
            //TLorentzVector pipiRec(0.,0.,0.,0.);
            //...Luiz
            TLorentzVector kkkkRec(0., 0., 0., 0.);

            int totcharge = 0;

            TLorentzVector kkRec(0., 0., 0., 0.);
            TLorentzVector mmRec(0., 0., 0., 0.);
            TLorentzVector eeRec(0., 0., 0., 0.);

            //int charray[2]={0,0};
            //double chi2array[2]={0.,0.};
            //double d0array[2]={0.,0.};
            //double dzarray[2]={0.,0.};
            //...Luiz
            int charray[4] = {0, 0, 0, 0};
            double chi2array[4] = {0., 0., 0., 0.};
            double d0array[4] = {0., 0., 0., 0.};
            double dzarray[4] = {0., 0., 0., 0.};
            int pidarray[4] = {0, 0, 0, 0};

            int ntrk0 = 0;
            int ntrk = 0;
            int ntrkvtx = 0;

            //       for(TrackCollection::const_iterator itTrack = tracks->begin();itTrack != tracks->end();++itTrack) {
            for(vector<MyTracks>::iterator itTrack = track_coll->begin() ; itTrack != track_coll->end() ; ++itTrack)
            {

                int looper = itTrack->isLooper;
                double pt = itTrack->pt();
                double pz = itTrack->pz();
                double eta = itTrack->eta();
                double phi = itTrack->phi();
                double charge = itTrack->charge;
                int npixelhits = itTrack->nValidPixelHits;
                int nstriphits = itTrack->nValidStripHits;
                int algo = itTrack->trackAlgo;
                double chi2 = itTrack->chi2n;
                double d0 = itTrack->d0;
                double dz = itTrack->dz;

                histosTH1F["hpt"]->Fill(pt);
                histosTH1F["heta"]->Fill(eta);
                histosTH1F["hphi"]->Fill(phi);
                histosTH1F["halgo"]->Fill(algo);
                histosTH1F["hnhits"]->Fill(npixelhits + nstriphits);

                ntrk0++;

                //inner tracker filter
                if(npixelhits > 0)
                {
                    //    if(npixelhits>0 && TMath::Abs(d0)<1. && TMath::Abs(dz)<20.){

                    histosTH1F["hlooper"]->Fill(looper);
                    histosTH1F["hchi2"]->Fill(chi2);
                    histosTH1F["hd0"]->Fill(d0);
                    histosTH1F["hdz"]->Fill(dz);

                    histosTH2F["hdedx"]->Fill(itTrack->p, itTrack->harmonic2_dEdx);
                    //...Luiz
                    double lndEdx = TMath::Log(itTrack->harmonic2_dEdx);
                    histosTH2F["hlndedx"]->Fill(itTrack->p, lndEdx);
                    double l10dEdx = TMath::Log10(itTrack->harmonic2_dEdx);
                    histosTH2F["hl10dedx"]->Fill(itTrack->p, l10dEdx);

                    //...Luiz
                    totcharge += itTrack->charge;
                    double ene = TMath::Sqrt(pt * pt + pz * pz + m_k * m_k); // energy related to kaon mass
                    TLorentzVector trk_lorentz(itTrack->px(), itTrack->py(), itTrack->pz(), ene);
                    kkkkRec += trk_lorentz;

                    // naming tracks, treating each kaon individually
                    if(ntrk == 0) k1 = trk_lorentz;
                    if(ntrk == 1) k2 = trk_lorentz;
                    if(ntrk == 2) k3 = trk_lorentz;
                    if(ntrk == 3) k4 = trk_lorentz;
                    if(ntrk == 0 || ntrk == 1) k1k2Rec += trk_lorentz;
                    if(ntrk == 2 || ntrk == 3) k3k4Rec += trk_lorentz;
                    if(ntrk == 0 || ntrk == 2) k1k3Rec += trk_lorentz;
                    if(ntrk == 1 || ntrk == 3) k2k4Rec += trk_lorentz;
                    EPID pid2 = GetPIDSafe2(itTrack->p, itTrack->harmonic2_dEdx);

                    if(ntrk == 0)
                    {
                        charray[0] = charge;
                        chi2array[0] = chi2;
                        d0array[0] = d0;
                        dzarray[0] = dz;
                        pidarray[0] = pid2;
                    }
                    if(ntrk == 1)
                    {
                        charray[1] = charge;
                        chi2array[1] = chi2;
                        d0array[1] = d0;
                        dzarray[1] = dz;
                        pidarray[1] = pid2;
                    }
                    if(ntrk == 2)
                    {
                        charray[2] = charge;
                        chi2array[2] = chi2;
                        d0array[2] = d0;
                        dzarray[2] = dz;
                        pidarray[2] = pid2;
                    }
                    if(ntrk == 3)
                    {
                        charray[3] = charge;
                        chi2array[3] = chi2;
                        d0array[3] = d0;
                        dzarray[3] = dz;
                        pidarray[3] = pid2;
                    }

                    ntrk++;

                    //-----------------------
                    double eneK = TMath::Sqrt(pt * pt + pz * pz + m_k * m_k);
                    TLorentzVector trk_lorentzK(itTrack->px(), itTrack->py(), itTrack->pz(), eneK);
                    kkRec += trk_lorentzK;

                    double eneM = TMath::Sqrt(pt * pt + pz * pz + m_mu * m_mu);
                    TLorentzVector trk_lorentzM(itTrack->px(), itTrack->py(), itTrack->pz(), eneM);
                    mmRec += trk_lorentzM;
                    double eneE = TMath::Sqrt(pt * pt + pz * pz + m_e * m_e);
                    TLorentzVector trk_lorentzE(itTrack->px(), itTrack->py(), itTrack->pz(), eneE);
                    eeRec += trk_lorentzE;

                }
            }



            histosTH1F["hntrk0"]->Fill(ntrk0);
            histosTH1F["hntrk"]->Fill(ntrk);

            if(ntrk == 0)
            {
                int nclusters =  sipixelcluster_coll->size();
                int nclusters2 = sistripcluster_coll->nStripClusters;

                histosTH1F["hnclusters"]->Fill(nclusters);
                histosTH1F["hnclusters2"]->Fill(nclusters2);
            }

            int nvtx = 0;
            //       for(VertexCollection::const_iterator itVtx = vertices->begin();itVtx != vertices->end();++itVtx) {
            for(vector<MyVertex>::iterator itVtx = vertex_coll->begin() ; itVtx != vertex_coll->end() ; ++itVtx)  //loop to find vertices
            {
                int vtxisfake = itVtx->fake;
                if(vtxisfake == 0) nvtx++;
                else continue;

                ntrkvtx = itVtx->ntracks;

            }

            histosTH1F["hnvtx"]->Fill(nvtx);
            if(nvtx == 1) histosTH1F["hntrkvtx"]->Fill(ntrkvtx);

            //not yet vertex cut, checking vertex-finding efficiency

            int  isfake = vertex_coll->begin()->fake;
            double xvtx = vertex_coll->begin()->x;
            double yvtx = vertex_coll->begin()->y;
            double zvtx = vertex_coll->begin()->z;

            //       double chi2vtx = vertices->begin()->normalizedChi2();
            // not sure if the same variable.
            // myvertex.chi2        = p->chi2();
            double chi2vtx = vertex_coll->begin()->chi2; //chi^2 fitting
            //       double ndofvtx = vertex_coll->begin()->ndof;

            //       ntrkvtx = vertex_coll->begin()->ntracks;


            //for vertex plots
            //...Luiz  ntrk==4
            //fiducialRegion   = (ntrk==2 && TMath::Abs(pi1.Eta())<etaCut && TMath::Abs(pi2.Eta())<etaCut);
            //fiducialRegionPt = (ntrk==2 && pi1.Pt()>ptCut && pi2.Pt()>ptCut);
            //...Luiz
            // fiducial region is a window decided by detector, eta is invariant variable for theta angle

            //pt, phi, and eta are needed to know where particle is in detector
            //pt cut is 0.1, eta cut is 2.5
            fiducialRegion   = (ntrk == 4 && TMath::Abs(k1.Eta()) < etaCut && TMath::Abs(k2.Eta()) < etaCut &&
                                TMath::Abs(k3.Eta()) < etaCut && TMath::Abs(k4.Eta()) < etaCut);
            fiducialRegionPt = (ntrk == 4 && k1.Pt() > ptCut && k2.Pt() > ptCut &&
                                k3.Pt() > ptCut && k4.Pt() > ptCut);
            ////fiducialRegion   = (ntrk==4);
            ////fiducialRegionPt = (ntrk==4);
            histosTH1F["hvtx"]->Fill( isfake );
            //...Luiz
            if(ntrk == 4)
            {
                histosTH1F["hvtx2"]->Fill( isfake );
                if(fiducialRegion && totcharge == 0) histosTH1F["hvtx3"]->Fill( isfake );
            }

            //...Luiz ?????? cut nvtx=1 or 2
            //if(nvtx!=1) continue;
            //
            // if(nvtx!=2) continue;

            if(nvtx != 1)
            {
                if(nvtx != 2) continue; //select number of vertices we want
            }



            histosTH1F["hvtxx"]->Fill(xvtx);
            histosTH1F["hvtxy"]->Fill(yvtx);
            histosTH1F["hvtxz"]->Fill(zvtx);
            histosTH1F["hvtxchi2"]->Fill(chi2vtx);
            histosTH2F["hntrkntrkvtx"]->Fill(ntrkvtx, ntrk);

            //----------------------------
            //invariant mass
            //...Luiz
            double mrec = kkkkRec.M();
            double mrecKK = kkRec.M();
            double mrecMM = mmRec.M(); //reconstructed masses
            double mrecEE = eeRec.M();

            double mreck1k2 = k1k2Rec.M();
            double mreck3k4 = k3k4Rec.M();
            double mreck1k3 = k1k3Rec.M();
            double mreck2k4 = k2k4Rec.M();

            //----------------------------
            // xi cut
            // Mmax=13000*xi_max
            // 0.1 -> 1300 GeV
            // 0.01 -> 130 GeV
            // 0.001 -> 13 GeV

            //...Luiz...rapidity = 1/2 ln ( xi_proton_2/xi_proton_1 )
            double rapy = 0.5 * TMath::Log(xiR / xiL);
            //

            if(fiducialRegion && fiducialRegionPt)
            {
                histosTH1F["hxiL"]->Fill(xiL);
                histosTH1F["hxiR"]->Fill(xiR);
                histosTH1F["hm"]->Fill(mrec);
                //...Luiz
                histosTH1F["hrapy"]->Fill(rapy);
            }

            // last one, before Simone
            if(TMath::Abs(xiL) < 0.02 && TMath::Abs(xiR) < 0.02);
            //  if(TMath::Abs(xiL)<0.01 && TMath::Abs(xiR)<0.01);
            else continue;

            if(fiducialRegion && fiducialRegionPt)
            {
                histosTH1F["hxiL2"]->Fill(xiL);
                histosTH1F["hxiR2"]->Fill(xiR);
                histosTH1F["hmxicut"]->Fill(mrec);
                //...Luiz
                histosTH1F["hrapy2"]->Fill(rapy);
            }

            //-----------------------------
            // balance cut - py cut


            // was in the first submission for 9919,9922
            //  double CMSpx=pipiRec.Px();
            //  double CMSpy=pipiRec.Py();
            //
            //...Luiz | information for kaon tracks
            double CMSpx = kkkkRec.Px();
            double CMSpy = kkkkRec.Py();

            histosTH2F["h2dimdpyAll"]->Fill(CMSpy, TOTEMpy);
            histosTH1F["hdpyAll"]->Fill(CMSpy + TOTEMpy);

            if(fiducialRegion && fiducialRegionPt)
            {
                histosTH2F["h2dimdpy"]->Fill(CMSpy, TOTEMpy);
                histosTH1F["hdpy"]->Fill(CMSpy + TOTEMpy);

                if(diag)
                {
                    histosTH2F["h2dimdpy_diag"]->Fill(CMSpy, TOTEMpy);
                    histosTH1F["hdpy_diag"]->Fill(CMSpy + TOTEMpy);
                }
                else
                {
                    histosTH2F["h2dimdpy_ttbb"]->Fill(CMSpy, TOTEMpy);
                    histosTH1F["hdpy_ttbb"]->Fill(CMSpy + TOTEMpy);
                }
            }


            // last one, before Simone
            bool CTpycut = TMath::Abs(CMSpy + TOTEMpy) < 0.06;
            //  bool CTpycut = TMath::Abs(CMSpy+TOTEMpy)<0.03;
            //  bool CTpycut = TMath::Abs(CMSpy+TOTEMpy)<0.015; // 1/4

            if(!CTpycut) continue;

            // px for completeness
            histosTH2F["h2dimdpxAll"]->Fill(CMSpx, TOTEMpx);
            histosTH1F["hdpxAll"]->Fill(CMSpx + TOTEMpx);

            if(fiducialRegion && fiducialRegionPt)
            {
                histosTH2F["h2dimdpx"]->Fill(CMSpx, TOTEMpx);
                histosTH1F["hdpx"]->Fill(CMSpx + TOTEMpx);

                if(diag)
                {
                    histosTH2F["h2dimdpx_diag"]->Fill(CMSpx, TOTEMpx);
                    histosTH1F["hdpx_diag"]->Fill(CMSpx + TOTEMpx);
                }
                else
                {
                    histosTH2F["h2dimdpx_ttbb"]->Fill(CMSpx, TOTEMpx);
                    histosTH1F["hdpx_ttbb"]->Fill(CMSpx + TOTEMpx);
                }

            }

            // last one, before Simone
            bool CTpxcut = TMath::Abs(CMSpx + TOTEMpx) < 0.15;
            //  bool CTpxcut = TMath::Abs(CMSpx+TOTEMpx)<0.075;
            //  bool CTpxcut = TMath::Abs(CMSpx+TOTEMpx)<0.0375; // 1/4

            //-----------------------------
            // from now on, only 2 vertex tracks. |eta|<etaCut (=2.5)
            //    if(!fiducialRegion) continue;

            //-----------------------------
            //RP vertex

            //xVtxL,xVtxR in meters, see Fig. 6 of Hubert's PAS
            // so meed to multiply by 100 to get in cm

            // 2012
            //  double vertexResolution = 8.3e-6;
            //  bool RPvertex = abs(xVtxL-xVtxR) < 3*vertexResolution;
            // 2015 - Mirko
            //last one, before Simone
            bool RPvertex = abs(xVtxL - xVtxR) < 3e-5;
            //  bool RPvertex = abs(xVtxL-xVtxR) < 1.5e-5;
            //  bool RPvertex = abs(xVtxL-xVtxR) < 0.75e-5; //1/4

            double xvtxT = (xVtxR + xVtxL) / 2.;
            //last one before Simone
            bool CTvertex = -0.04 < (xvtx - xvtxT * 1e2) && (xvtx - xvtxT * 1e2) < 0.18;
            //  bool CTvertex = 0.015<(xvtx-xvtxT*1e2) && (xvtx-xvtxT*1e2)<0.125;
            //  bool CTvertex = 0.04<(xvtx-xvtxT*1e2) && (xvtx-xvtxT*1e2)<0.1; //1/4

            // my
            //  bool RPvertex = abs(xVtxL-xVtxR) < 5e-5;
            //  bool RPvertex = abs(xVtxL-xVtxR) < 4e-5;

            // reject if no RP vertex
            //0.025mm = 25\mum
            // < 3 * vertexResolution = 0.0000249               // <- mm
            //      	if(TMath::Abs(dx0)>0.025) continue;   // <- mum


            if(fiducialRegion && fiducialRegionPt)
            {

                histosTH2F["h2dimxVtxRL"]->Fill(xVtxL, xVtxR);
                histosTH2F["h2dimxVtxcmsR"]->Fill(xVtxR * 1e2, xvtx);
                histosTH2F["h2dimxVtxcmsL"]->Fill(xVtxL * 1e2, xvtx);

                histosTH2F["h2dimxVtxcmsRL"]->Fill(xvtxT * 1e2, xvtx);

                if(RPvertex)
                {
                    histosTH2F["h2dimxVtxcmsR2"]->Fill(xVtxR * 1e2, xvtx);
                    histosTH2F["h2dimxVtxcmsL2"]->Fill(xVtxL * 1e2, xvtx);
                    histosTH2F["h2dimxVtxcmsRL2"]->Fill(xvtxT * 1e2, xvtx);

                    histosTH2F["h2dimxVtx_zVtx_CT"]->Fill(zvtx, xvtx - xvtxT * 1e2);
                    histosTH2F["h2dimxVtx_zVtx_C"]->Fill(zvtx, xvtx);
                    histosTH2F["h2dimxVtx_zVtx_T"]->Fill(zvtx, xvtxT * 1e2);

                    histosTH1F["hxVtxcmsRL"]->Fill(xvtx - xvtxT * 1.e2);
                    if(diag)  histosTH1F["hxVtxcmsRL_diag"]->Fill(xvtx - xvtxT * 1.e2);
                    else      histosTH1F["hxVtxcmsRL_ttbb"]->Fill(xvtx - xvtxT * 1.e2);
                }

                histosTH1F["hxVtxRL"]->Fill(xVtxR - xVtxL);
                histosTH1F["hxVtxcmsR"]->Fill(xvtx - xVtxR * 1.e2);
                histosTH1F["hxVtxcmsL"]->Fill(xvtx - xVtxL * 1.e2);

                if(diag)
                {
                    histosTH1F["hxVtxRL_diag"]->Fill(xVtxR - xVtxL);
                    histosTH1F["hxVtxcmsR_diag"]->Fill(xvtx - xVtxR * 1e2);
                    histosTH1F["hxVtxcmsL_diag"]->Fill(xvtx - xVtxL * 1e2);
                }
                else
                {
                    histosTH1F["hxVtxRL_ttbb"]->Fill(xVtxR - xVtxL);
                    histosTH1F["hxVtxcmsR_ttbb"]->Fill(xvtx - xVtxR * 1e2);
                    histosTH1F["hxVtxcmsL_ttbb"]->Fill(xvtx - xVtxL * 1e2);
                }

            }

            //...cut 1    nvtx==1 or 2

            //-----------------------------
            ////fR: ntrk==2, nvtx==1, |eta|<etaCut
            //
            //...Luiz
            //fR: ntrk==4, nvtx==1 or 2, |eta|<etaCut

            //  totcharge=totcharge0;

            if(fiducialRegion && fiducialRegionPt)
            {

                // how many tracks with pixel if at vertex 2 tracks
                histosTH1F["hntrkntrkvtx2"]->Fill(ntrk);
                histosTH1F["hntrk2ntrkvtx"]->Fill(ntrkvtx);

                histosTH1F["hm2rec"]->Fill(mrec);
                histosTH1F["hm2recbis"]->Fill(mrec);

                if(totcharge == 0)
                {
                    histosTH1F["hm2recOS"]->Fill(mrec);
                    if(diag) histosTH1F["hm2recOS_diag"]->Fill(mrec);
                    else     histosTH1F["hm2recOS_ttbb"]->Fill(mrec);
                }
                else
                {
                    histosTH1F["hm2recSS"]->Fill(mrec);
                    if(diag) histosTH1F["hm2recSS_diag"]->Fill(mrec);
                    else     histosTH1F["hm2recSS_ttbb"]->Fill(mrec);
                }

                //...cut 2
                // if(CTpxcut)
                //{
                if(pidarray[0] == pidKaon && pidarray[1] == pidKaon && pidarray[2] == pidKaon && pidarray[3] == pidKaon)
                {

                    if(totcharge == 0)
                    {
                        histosTH1F["hm2rec20S"]->Fill(mrec);
                        if(nvtx == 1);
                            {
                        if(charray[0]+charray[1] == 0);
                            histosTH1F["hm2rec2OS_k1k2"]->Fill(mreck1k2);
                        if(charray[2]+charray[3] == 0)
                            histosTH1F["hm2rec2OS_k3k4"]->Fill(mreck3k4);
                        if(charray[0]+charray[2] == 0)
                            histosTH1F["hm2rec2OS_k1k3"]->Fill(mreck1k3);
                        if(charray[1]+charray[3] == 0)
                            histosTH1F["hm2rec2OS_k2k4"]->Fill(mreck2k4);
                    }
                         if(nvtx == 2);
                            {
                        if(charray[0]+charray[1] == 0);
                            histosTH1F["hm2rec2OS_k1k2v2"]->Fill(mreck1k2);
                        if(charray[2]+charray[3] == 0)
                            histosTH1F["hm2rec2OS_k3k4v2"]->Fill(mreck3k4);
                        if(charray[0]+charray[2] == 0)
                            histosTH1F["hm2rec2OS_k1k3v2"]->Fill(mreck1k3);
                        if(charray[1]+charray[3] == 0)
                            histosTH1F["hm2rec2OS_k2k4v2"]->Fill(mreck2k4);
                    }
                        if(diag) histosTH1F["hm2rec2OS_diag"]->Fill(mrec);
                        else     histosTH1F["hm2rec2OS_ttbb"]->Fill(mrec);
                    }
                    else
                    {
                        histosTH1F["hm2rec2SS"]->Fill(mrec);
                        if(diag) histosTH1F["hm2rec2SS_diag"]->Fill(mrec);
                        else     histosTH1F["hm2rec2SS_ttbb"]->Fill(mrec);
                    }
                    //....???????
                    if(totcharge == 0 && diag)
                    {
                        if(k1.Py()*k2.Py()*k3.Py()*k4.Py() > 0) histosTH1F["hm2rec2OS_diag_trkP"]->Fill(mrec);
                        if(k1.Py()*k2.Py()*k3.Py()*k4.Py() < 0) histosTH1F["hm2rec2OS_diag_trkM"]->Fill(mrec);
                    }
                    if(totcharge == 0 && !diag)
                    {
                        if(k1.Py()*k2.Py()*k3.Py()*k4.Py() > 0) histosTH1F["hm2rec2OS_ttbb_trkP"]->Fill(mrec);
                        if(k1.Py()*k2.Py()*k3.Py()*k4.Py() < 0) histosTH1F["hm2rec2OS_ttbb_trkM"]->Fill(mrec);
                    }

                    //...Luiz  ??????
                    if(totcharge == 0 && diag)
                    {
                        if(TMath::Abs(kkkkRec.Py()) > TMath::Abs(kkkkRec.Px())) histosTH1F["hm2rec2OS_diag_pypxP"]->Fill(mrec);
                        else histosTH1F["hm2rec2OS_diag_pypxM"]->Fill(mrec);
                        //...Luiz
                        int pypx = 0;
                        if(TMath::Abs(kkkkRec.Py()) > TMath::Abs(kkkkRec.Px())) pypx = 1;
                        else pypx = 0;
                        //...?????
                        if(mrec >= 1.65 && mrec <= 1.75) cout << "scan2OSdiag: " << run << " " << LS << " " << evt << " " << mrec << " " << pypx << endl;
                    }
                    //...Luiz
                    if(totcharge == 0 && !diag)
                    {
                        if(TMath::Abs(kkkkRec.Py()) > TMath::Abs(kkkkRec.Px())) histosTH1F["hm2rec2OS_ttbb_pypxP"]->Fill(mrec);
                        else  histosTH1F["hm2rec2OS_ttbb_pypxM"]->Fill(mrec);
                    }

                    //--------------------------

                    //...Luiz
                    if(totcharge == 0 && diag)
                    {
                        histosTH1F["hm2recPPPP"]->Fill(mrec);
                        histosTH1F["hm2recKK"]->Fill(mrecKK);
                        histosTH1F["hm2recMM"]->Fill(mrecMM);
                        histosTH1F["hm2recEE"]->Fill(mrecEE);
                    }

                    }
                 //}

                //...cut 3
                //.... OS:totcharge==0 SS:totcharge!=0
                if(RPvertex && CTpxcut)
                {
                    if(totcharge == 0)
                    {
                        histosTH1F["hm2rec3OS"]->Fill(mrec);
                        if(diag) histosTH1F["hm2rec3OS_diag"]->Fill(mrec);
                        else     histosTH1F["hm2rec3OS_ttbb"]->Fill(mrec);
                    }
                    else
                    {
                        histosTH1F["hm2rec3SS"]->Fill(mrec);
                        if(diag) histosTH1F["hm2rec3SS_diag"]->Fill(mrec);
                        else     histosTH1F["hm2rec3SS_ttbb"]->Fill(mrec);
                    }

                    //...Luiz
                    //             OS: pi+pi+ or pi-pi-    ?????
                    if(totcharge == 0 && diag)
                    {
                        if(k1.Py()*k2.Py()*k3.Py()*k4.Py() > 0) histosTH1F["hm2rec3OS_diag_trkP"]->Fill(mrec);
                        if(k1.Py()*k2.Py()*k3.Py()*k4.Py() < 0) histosTH1F["hm2rec3OS_diag_trkM"]->Fill(mrec);
                    }
                    if(totcharge == 0 && !diag)
                    {
                        if(k1.Py()*k2.Py()*k3.Py()*k4.Py() > 0) histosTH1F["hm2rec3OS_ttbb_trkP"]->Fill(mrec);
                        if(k1.Py()*k2.Py()*k3.Py()*k4.Py() < 0) histosTH1F["hm2rec3OS_ttbb_trkM"]->Fill(mrec);
                    }
                    //...Luiz
                    if(totcharge == 0 && diag)
                    {
                        if(TMath::Abs(kkkkRec.Py()) > TMath::Abs(kkkkRec.Px())) histosTH1F["hm2rec3OS_diag_pypxP"]->Fill(mrec);
                        else histosTH1F["hm2rec3OS_diag_pypxM"]->Fill(mrec);
                    }
                    if(totcharge == 0 && !diag)
                    {
                        if(TMath::Abs(kkkkRec.Py()) > TMath::Abs(kkkkRec.Px())) histosTH1F["hm2rec3OS_ttbb_pypxP"]->Fill(mrec);
                        else  histosTH1F["hm2rec3OS_ttbb_pypxM"]->Fill(mrec);
                    }

                }

                //...cut 4
                //.... OS:totcharge==0 SS:totcharge!=0
                //    if(RPvertex && CTpxcut && nvtx==1 && CTvertex){
                if(RPvertex && CTpxcut && CTvertex && TMath::Abs(zvtx) < 5.) // core
                {

                    if(totcharge == 0)
                    {
                        histosTH1F["hm2rec4OS"]->Fill(mrec);
                        if(diag) histosTH1F["hm2rec4OS_diag"]->Fill(mrec);
                        else     histosTH1F["hm2rec4OS_ttbb"]->Fill(mrec);
                    }
                    else
                    {
                        histosTH1F["hm2rec4SS"]->Fill(mrec);
                        if(diag) histosTH1F["hm2rec4SS_diag"]->Fill(mrec);
                        else     histosTH1F["hm2rec4SS_ttbb"]->Fill(mrec);
                    }
                    //...Luiz
                    //             OS: pi+pi+ or pi-pi-    ?????
                    if(totcharge == 0 && diag)
                    {
                        if(k1.Py()*k2.Py()*k3.Py()*k4.Py() > 0) histosTH1F["hm2rec4OS_diag_trkP"]->Fill(mrec);
                        if(k1.Py()*k2.Py()*k3.Py()*k4.Py() < 0) histosTH1F["hm2rec4OS_diag_trkM"]->Fill(mrec);
                    }
                    //  ?????
                    if(totcharge == 0 && !diag)
                    {
                        if(k1.Py()*k2.Py()*k3.Py()*k4.Py() > 0) histosTH1F["hm2rec4OS_ttbb_trkP"]->Fill(mrec);
                        if(k1.Py()*k2.Py()*k3.Py()*k4.Py() < 0) histosTH1F["hm2rec4OS_ttbb_trkM"]->Fill(mrec);
                    }
                    //...Luiz
                    if(totcharge == 0 && diag)
                    {
                        if(TMath::Abs(kkkkRec.Py()) > TMath::Abs(kkkkRec.Px())) histosTH1F["hm2rec4OS_diag_pypxP"]->Fill(mrec);
                        else histosTH1F["hm2rec4OS_diag_pypxM"]->Fill(mrec);
                    }
                    //...Luiz
                    if(totcharge == 0 && !diag)
                    {
                        if(TMath::Abs(kkkkRec.Py()) > TMath::Abs(kkkkRec.Px())) histosTH1F["hm2rec4OS_ttbb_pypxP"]->Fill(mrec);
                        else  histosTH1F["hm2rec4OS_ttbb_pypxM"]->Fill(mrec);
                    }

                }

                //...cut 5   nvtx==1 or 2
                //.... OS:totcharge==0 SS:totcharge!=0
                //    if(RPvertex && CTpxcut && nvtx==1 && CTvertex && TMath::Abs(zvtx)>5.){//tails
                //    if(RPvertex && CTpxcut && nvtx==1 && CTvertex){
                // no dpx cut applied
                if(RPvertex && CTvertex)
                {

                    if(totcharge == 0)
                    {
                        histosTH1F["hm2rec5OS"]->Fill(mrec);
                        if(diag) histosTH1F["hm2rec5OS_diag"]->Fill(mrec);
                        else     histosTH1F["hm2rec5OS_ttbb"]->Fill(mrec);
                    }
                    else
                    {
                        histosTH1F["hm2rec5SS"]->Fill(mrec);
                        if(diag) histosTH1F["hm2rec5SS_diag"]->Fill(mrec);
                        else     histosTH1F["hm2rec5SS_ttbb"]->Fill(mrec);
                    }
                    //...Luiz
                    //             OS: pi+pi+ or pi-pi-    ?????
                    if(totcharge == 0 && diag)
                    {
                        if(k1.Py()*k2.Py()*k3.Py()*k4.Py() > 0) histosTH1F["hm2rec5OS_diag_trkP"]->Fill(mrec);
                        if(k1.Py()*k2.Py()*k3.Py()*k4.Py() < 0) histosTH1F["hm2rec5OS_diag_trkM"]->Fill(mrec);
                    }
                    //   ?????
                    if(totcharge == 0 && !diag)
                    {
                        if(k1.Py()*k2.Py()*k3.Py()*k4.Py() > 0) histosTH1F["hm2rec5OS_ttbb_trkP"]->Fill(mrec);
                        if(k1.Py()*k2.Py()*k3.Py()*k4.Py() < 0) histosTH1F["hm2rec5OS_ttbb_trkM"]->Fill(mrec);
                    }
                    //...Luiz
                    if(totcharge == 0 && diag)
                    {
                        if(TMath::Abs(kkkkRec.Py()) > TMath::Abs(kkkkRec.Px())) histosTH1F["hm2rec5OS_diag_pypxP"]->Fill(mrec);
                        else histosTH1F["hm2rec5OS_diag_pypxM"]->Fill(mrec);
                    }
                    //...Luiz
                    if(totcharge == 0 && !diag)
                    {
                        if(TMath::Abs(kkkkRec.Py()) > TMath::Abs(kkkkRec.Px())) histosTH1F["hm2rec5OS_ttbb_pypxP"]->Fill(mrec);
                        else  histosTH1F["hm2rec5OS_ttbb_pypxM"]->Fill(mrec);
                    }

                }

                //...cut 6    nvtx==1 or 2
                //.... OS:totcharge==0 SS:totcharge!=0
                //    if(RPvertex && CTpxcut && nvtx==1){
                if(RPvertex && CTpxcut && CTvertex)
                {

                    double etaCut2 = 1.5;
                    //...Luiz
                    if(TMath::Abs(k1.Eta()) < etaCut2 && TMath::Abs(k2.Eta()) < etaCut2 &&
                            TMath::Abs(k3.Eta()) < etaCut2 && TMath::Abs(k4.Eta()) < etaCut2 )
                    {
                        if(totcharge == 0)
                        {
                            histosTH1F["hm2rec6OS"]->Fill(mrec);
                            if(diag) histosTH1F["hm2rec6OS_diag"]->Fill(mrec);
                            else     histosTH1F["hm2rec6OS_ttbb"]->Fill(mrec);
                        }
                        else
                        {
                            histosTH1F["hm2rec6SS"]->Fill(mrec);
                            if(diag) histosTH1F["hm2rec6SS_diag"]->Fill(mrec);
                            else     histosTH1F["hm2rec6SS_ttbb"]->Fill(mrec);
                        }
                        //...Luiz           ?????
                        //             OS: pi+pi+ or pi-pi-    ?????
                        if(totcharge == 0 && diag)
                        {
                            if(k1.Py()*k2.Py()*k3.Py()*k4.Py() > 0) histosTH1F["hm2rec6OS_diag_trkP"]->Fill(mrec);
                            if(k1.Py()*k2.Py()*k3.Py()*k4.Py() < 0) histosTH1F["hm2rec6OS_diag_trkM"]->Fill(mrec);
                        }
                        //...Luiz           ?????
                        if(totcharge == 0 && !diag)
                        {
                            if(k1.Py()*k2.Py()*k3.Py()*k4.Py() > 0) histosTH1F["hm2rec6OS_ttbb_trkP"]->Fill(mrec);
                            if(k1.Py()*k2.Py()*k3.Py()*k4.Py() < 0) histosTH1F["hm2rec6OS_ttbb_trkM"]->Fill(mrec);
                        }

                    }
                }

                //...cut 7
                if(diag && RPvertex && CTpxcut)
                {

                    if(totcharge == 0) histosTH1F["hm2recHFvetoOS"]->Fill(mrec);
                    else histosTH1F["hm2recHFvetoSS"]->Fill(mrec);

                    //-------------
                    //...Luiz   ????
                    if(k1.Pt() > 0.45 && k2.Pt() > 0.45 && k3.Pt() > 0.45 && k4.Pt() > 0.45)
                    {
                        if(totcharge == 0) histosTH1F["hm2rec45OS"]->Fill(mrec);
                        else histosTH1F["hm2rec45SS"]->Fill(mrec);

                        double etaCut2 = 1.5;
                        //...Luiz    ??????
                        if(TMath::Abs(k1.Eta()) < etaCut2 && TMath::Abs(k2.Eta()) < etaCut2 &&
                                TMath::Abs(k3.Eta()) < etaCut2 && TMath::Abs(k4.Eta()) < etaCut2)
                        {
                            if(totcharge == 0) histosTH1F["hm2rec4515OS"]->Fill(mrec);
                            else histosTH1F["hm2rec4515SS"]->Fill(mrec);
                        }
                    }

                }

                if(diag && totcharge == 0 && RPvertex && CTpxcut)
                {
                    //      if(run==259237 && LS>=78 && LS<=100) histosTH1F["hm2rec9919"]->Fill(mrec);
                    if(run == 259237 && LS >= 78 && LS <= 100) histosTH1F["hm2rec9922"]->Fill(mrec);
                    if(run == 259237 && LS >= 432 && LS <= 576) histosTH1F["hm2rec9922"]->Fill(mrec);
                    if(run == 259385 && LS >= 253 && LS <= 538) histosTH1F["hm2rec9971"]->Fill(mrec);
                    if(run == 259388 && LS >= 369 && LS <= 747) histosTH1F["hm2rec9978"]->Fill(mrec);
                }

            }

            //-----------------------------
            // track variables
            //...Luiz
            if(ntrk == 4 && CTpycut && CTpxcut && RPvertex)
            {
                if(totcharge == 0 && diag)
                {

                    histosTH1F["hptRes"]->Fill(kkkkRec.Pt());
                    histosTH1F["hetaRes"]->Fill(kkkkRec.Eta());
                    histosTH1F["hphiRes"]->Fill(kkkkRec.Phi());

                    if(charray[0] > 0)
                    {
                        histosTH1F["hptP"]->Fill(k1.Pt());
                        histosTH1F["hetaP"]->Fill(k1.Eta());
                        histosTH1F["hphiP"]->Fill(k1.Phi());
                    }
                    else
                    {
                        histosTH1F["hptM"]->Fill(k1.Pt());
                        histosTH1F["hetaM"]->Fill(k1.Eta());
                        histosTH1F["hphiM"]->Fill(k1.Phi());
                    }

                    if(charray[1] > 0)
                    {
                        histosTH1F["hptP"]->Fill(k2.Pt());
                        histosTH1F["hetaP"]->Fill(k2.Eta());
                        histosTH1F["hphiP"]->Fill(k2.Phi());
                    }
                    else
                    {
                        histosTH1F["hptM"]->Fill(k2.Pt());
                        histosTH1F["hetaM"]->Fill(k2.Eta());
                        histosTH1F["hphiM"]->Fill(k2.Phi());
                    }
                    //...Luiz
                    if(charray[2] > 0)
                    {
                        histosTH1F["hptP"]->Fill(k3.Pt());
                        histosTH1F["hetaP"]->Fill(k3.Eta());
                        histosTH1F["hphiP"]->Fill(k3.Phi());
                    }
                    else
                    {
                        histosTH1F["hptM"]->Fill(k3.Pt());
                        histosTH1F["hetaM"]->Fill(k3.Eta());
                        histosTH1F["hphiM"]->Fill(k3.Phi());
                    }
                    //...Luiz
                    if(charray[3] > 0)
                    {
                        histosTH1F["hptP"]->Fill(k4.Pt());
                        histosTH1F["hetaP"]->Fill(k4.Eta());
                        histosTH1F["hphiP"]->Fill(k4.Phi());
                    }
                    else
                    {
                        histosTH1F["hptM"]->Fill(k4.Pt());
                        histosTH1F["hetaM"]->Fill(k4.Eta());
                        histosTH1F["hphiM"]->Fill(k4.Phi());
                    }

                    histosTH1F["hvtxchi2fin"]->Fill(chi2vtx);

                    //...Luiz
                    histosTH1F["hchi2fin"]->Fill(chi2array[0]);
                    histosTH1F["hchi2fin"]->Fill(chi2array[1]);
                    histosTH1F["hchi2fin"]->Fill(chi2array[2]);
                    histosTH1F["hchi2fin"]->Fill(chi2array[3]);
                    histosTH1F["hd0fin"]->Fill(d0array[0]);
                    histosTH1F["hd0fin"]->Fill(d0array[1]);
                    histosTH1F["hd0fin"]->Fill(d0array[2]);
                    histosTH1F["hd0fin"]->Fill(d0array[3]);
                    histosTH1F["hdzfin"]->Fill(dzarray[0]);
                    histosTH1F["hdzfin"]->Fill(dzarray[1]);
                    histosTH1F["hdzfin"]->Fill(dzarray[2]);
                    histosTH1F["hdzfin"]->Fill(dzarray[3]);

                    int nclustersOSdiag =  sipixelcluster_coll->size();
                    int nclusters2OSdiag =  sistripcluster_coll->nStripClusters;
                    histosTH1F["hnclustersOSdiag"]->Fill(nclustersOSdiag);
                    histosTH1F["hnclusters2OSdiag"]->Fill(nclusters2OSdiag);

                }
            }

        } // End of loop over events in a file


        // Close current file
        file->Close();

    } // End of loop over files


    // Output file
    TFile *output = new TFile(outputFileName.c_str(), "RECREATE");
    output->cd();

    for(map<string, TH1F *>::iterator it_histo = histosTH1F.begin();
            it_histo != histosTH1F.end(); ++it_histo)
        (*it_histo).second->Write();
    for(map<string, TH2F *>::iterator it_histo = histosTH2F.begin();
            it_histo != histosTH2F.end(); ++it_histo)
        (*it_histo).second->Write();

    output->Close();

    //  fout.close();

}

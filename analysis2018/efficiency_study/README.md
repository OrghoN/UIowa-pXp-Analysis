This project is a monte carlo simulation to generate and plot information related to K0sK0s events decaying to positive and negative pions (pi+pi-pi+pi-).

To change data generation parameters such as number of events to be generated without editing the python file, a user can make edits to mymain02.cmnd. See
the information in the Pythia8 documentation on command files for more information on this.

NOTE: THIS PROJECT IS A WORK IN PROGRESS. CERTAIN ELEMENTS OF EVENT EFFICIENCIES AND ERRORS STILL NEED TO BE INCLUDED AND WILL BE ADDED TO OVER TIME. 


To ensure the code runs, please use the following steps:
  1. Install Pythia8 according to documentation found at http://home.thep.lu.se/Pythia/
  2. Download and unzip project
  3. Move extracted folders to the examples folder in the Pythia8 directory
  
******************************************
Functionality
******************************************
  Run data generation: python mymain02.py
  
  Generate hists from generated csv files or existing csv files dowloaded with the project: python sim_hists.py

The project contains csv files that contain data generated from 35 million events in the cmnd file. These contain the following data:

  invmass_eff_err.csv - K0sK0s -> pi+pi-pi+pi- event invariant mass, event efficiency, error on event efficiency
  
  min_max_cut1.csv - min pT, max pT, max abs value of eta for pion daughters of K0sK0s events 0.5 GeV and higher
  
  min_max_cut2.csv - min pT, max pT, max abs value of eta for pion daughters of K0sK0s events 1.0 GeV and higher
  
  min_max_cut3.csv - min pT, max pT, max abs value of eta for pion daughters of K0sK0s events 1.5 GeV and higher

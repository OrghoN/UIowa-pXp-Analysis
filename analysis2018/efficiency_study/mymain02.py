################################################################################
# Author: Thomas McDowell
# Date: December 1, 2020
# 
# A monte carlo simulation to create and plot information on 
#
# K0s ID = 310
# pi+ ID = 211
# pi- ID = -211
#
################################################################################

# Libraries
import sys
import numpy as np
import math
from matplotlib import pyplot as plt

# Particle IDs used within the simulaiton
K0s_ID = 310
PI_POS_ID = 211
PI_NEG_ID = -211
PI0_ID = 111

# Set config options
cfg = open("Makefile.inc")
lib = "../lib"

for line in cfg:
	if line.startswith("PREFIX_LIB="): lib = line[11:-1]; break
sys.path.insert(0, lib)

# Import Pythia, initialize simulation
import pythia8
pythia = pythia8.Pythia()
pythia.readFile("mymain02.cmnd")
pythia.init()

# File Objects to write data to
invmass_eff_err = open("invmass_eff_err.csv", "w")
#invmass_eff_err.write("invariant_mass,event_efficiency,fract_error\n")

# pion_pT_eta_vals is not currently in use but here in case it is needed later on
# pion_pT_eta_vals = open( "pion_pT_eta_vals.txt", "w" )
# pion_pT_eta_vals.write( "prt1_pT,prt_1_eta,prt2_pT,prt_2_eta,prt3_pT,prt_3_eta,prt4_pT,prt_4_eta\n" )

min_max_cut1 = open( "min_max_cut1.csv", "w" )
#min_max_cut1.write( "min_pT,max_pT,max_abs_eta\n" )

min_max_cut2 = open( "min_max_cut2.csv", "w" )
#min_max_cut2.write( "min_pT,max_pT,max_abs_eta\n" )

min_max_cut3 = open( "min_max_cut3.csv", "w" )
#min_max_cut3.write( "min_pT,max_pT,max_abs_eta\n" )

K0s_inv_masses = open( "K0s_inv_masses.csv", "w" )

# A few histogram definitions, may have a better place in the script
eta_hist = pythia8.Hist("pi+ pi- eta", 200, -25, 25)
pT_hist = pythia8.Hist("pi+ pi- pT", 200, -25, 25)
all_pT_hist = pythia8.Hist("pT of all K0s Decay", 200, -25, 25)
id_hist = pythia8.Hist("PID of all K0s Daughters", 1000, 0, 1000)

K0sK0s_daughter_pions = []

event_inv_masses = []

# Set used to not double count K0s in checkK0sK0s
K0s_set = set([])

# Some global variables to store information to be printed
num_K0s = 0
num_K0s_daughters = 0
num_pi0 = 0
num_K0sK0s = 0		

# bin weight list, not currently in use
# bin_weights = []

# Checks if a K0s particle object decays to pi+pi- particles
def checkPiDaughters( K0s ):
	pi_plus = 0
	pi_minus = 0
	for K0s_daughter in K0s.daughterList():
		if (pythia.event[K0s_daughter].id() == PI_POS_ID) and (pythia.event[K0s_daughter].pT() > 0.1) and (abs(pythia.event[K0s_daughter].eta()) < 2.5):
			pi_plus += 1
		if (pythia.event[K0s_daughter].id() == PI_NEG_ID) and (pythia.event[K0s_daughter].pT() > 0.1) and (abs(pythia.event[K0s_daughter].eta()) < 2.5):
			pi_minus += 1
	if (pi_plus == 1) and (pi_minus == 1):
		return True
	else:
		return False


# checkEvent checks the event list for K0s and calls checkK0s() on all that are 
# found while also counting each K0s
def checkEvent(event):
	global num_K0s

	for prt in event:
		if ( prt.id() == K0s_ID ):
			num_K0s += 1
			checkK0sK0s( prt )
			

# checks if a K0s particle comes from a K0sK0s event and if so examines the
# event as long as it has not done so already
def checkK0sK0s( prt ):
        # Global variables used to store values at different cuts and other
        # simulation info

	global invmass_eff_err
	global min_max_cut1
	global min_max_cut2
	global min_max_cut3
	global K0s_inv_masses

	global num_K0sK0s
	global K0s_set
	global K0sK0s_daughter_pions
	global event_inv_masses


	# If prt is a K0s particle object not in K0s_set already look at it's mother particle
	if (prt.id() == K0s_ID) and (prt not in K0s_set) and (checkPiDaughters( prt )):
		K0s_P = math.sqrt( prt.px()**2 + prt.py()**2 + prt.pz()**2 )
		K0s_inv_mass = ( ( prt.e()**2 ) - ( K0s_P**2 ) )
		K0s_inv_masses.write( str( K0s_inv_mass ) + '\n' )
		for mother in prt.motherList():
			# count K0s daughters, looking for two so event is K0sK0s
			K0s_daughter_count = 0
			for daughter in pythia.event[mother].daughterList():
				if ( abs( pythia.event[daughter].id() ) == K0s_ID ) and ( checkPiDaughters( pythia.event[daughter] ) ):
					K0s_daughter_count += 1
			
			# if the mother particle had two K0s daughters, check these daughters for pi+pi-		
			if (K0s_daughter_count == 2):
				# temp list to store a pi+ and pi- from tracks
				temp = []
				for daughter in pythia.event[mother].daughterList():
					if (pythia.event[daughter].id() == K0s_ID) and (pythia.event[daughter] not in K0s_set) and (checkPiDaughters(pythia.event[daughter])):
						# If a valid ( decays to pi+pi- ) K0s not in K0s_set is found, it is added to K0s_set
						K0s_set.add(pythia.event[daughter])
						

						for SecondGenDaughter in pythia.event[daughter].daughterList():
							if ( abs( pythia.event[SecondGenDaughter].id() ) == PI_POS_ID ):
								temp.append(pythia.event[SecondGenDaughter])

						# Make sure we have 4 tracks
						if ( len( temp ) == 4):
                                                        
                                                        event_prt_efficiencies = [ 0, 0, 0, 0 ]
							event_efficiency_errors = [ 0, 0, 0, 0 ]
							event_err = 0
							event_fract_err = 0
							event_momentum_list = [ 0, 0, 0, 0 ]
							event_energy_list = [ 0, 0, 0, 0 ]
							event_prt_inv_masses = [ 0, 0, 0, 0 ]

							# Set efficiency and error for each pi+ and pi- ( or each track ) to
							# values corresponding to those in the analysis note
                                                        for i in range( len( temp ) ):

						        	if temp[i].pT() >= 0.1 and temp[i].pT() < 0.15:
                                                            		event_prt_efficiencies[i] = 0.66
									event_efficiency_errors[i] = 0.06

                                                        	if temp[i].pT() >= 0.15 and temp[i].pT() < 0.2:
                                                                	event_prt_efficiencies[i] = 0.79
									event_efficiency_errors[i] = 0.05

                                                        	if temp[i].pT() >= 0.2 and temp[i].pT() < 0.25:
                                                                	event_prt_efficiencies[i] = 0.88
									event_efficiency_errors[i] = 0.02
                                                              
                                                        	if temp[i].pT() >= 0.25 and temp[i].pT() < 0.3:
                                                                	event_prt_efficiencies[i] = 0.89
									event_efficiency_errors[i] = 0.01

                                                        	if temp[i].pT() >= 0.3 and temp[i].pT() < 0.35:
                                                                	event_prt_efficiencies[i] = 0.88
									event_efficiency_errors[i] = 0.01

                                                        	if temp[i].pT() >= 0.35 and temp[i].pT() < 0.4:
                                                                	event_prt_efficiencies[i] = 0.89
									event_efficiency_errors[i] = 0.01

                                                        	if temp[i].pT() >= 0.4 and temp[i].pT() < 0.45:
                                                                	event_prt_efficiencies[i] = 0.92
									event_efficiency_errors[i] = 0.005

                                                        	if temp[i].pT() >= 0.45 and temp[i].pT() < 0.5:
                                                                	event_prt_efficiencies[i] = 0.92
									event_efficiency_errors[i] = 0.005

                                                        	if temp[i].pT() >= 0.5 and temp[i].pT() < 0.55:
                                                                	event_prt_efficiencies[i] = 0.94
									event_efficiency_errors[i] = 0.005

                                                        	if temp[i].pT() >= 0.55 and temp[i].pT() < 0.7:
                                                                	event_prt_efficiencies[i] = 0.95
									event_efficiency_errors[i] = 0.005

                                                        	if temp[i].pT() >= 0.7 and temp[i].pT() < 0.85:
                                                                	event_prt_efficiencies[i] = 0.96
									event_efficiency_errors[i] = 0.005

                                                        	if temp[i].pT() >= 0.85 and temp[i].pT() < 1.0:
                                                                	event_prt_efficiencies[i] = 0.97
									event_efficiency_errors[i] = 0.005

                                                        	if temp[i].pT() >= 1.0 and temp[i].pT() < 1.35:
                                                                	event_prt_efficiencies[i] = 0.96
									event_efficiency_errors[i] = 0.005

                                                        	if temp[i].pT() >= 1.35 and temp[i].pT() < 1.4:
                                                                	event_prt_efficiencies[i] = 0.99
									event_efficiency_errors[i] = 0.005

                                                        	if temp[i].pT() >= 1.4 and temp[i].pT() < 1.5:
                                                                	event_prt_efficiencies[i] = 0.97
									event_efficiency_errors[i] = 0.01

                                                        	if temp[i].pT() >= 1.5 and temp[i].pT() < 1.55:
                                                                	event_prt_efficiencies[i] = 0.96
									event_efficiency_errors[i] = 0.02

                                                        	if temp[i].pT() >= 1.55 and temp[i].pT() < 1.6:
                                                                	event_prt_efficiencies[i] = 0.91
									event_efficiency_errors[i] = 0.04

                                                        	if temp[i].pT() >= 1.6 and temp[i].pT() < 1.7:
                                                                	event_prt_efficiencies[i] = 0.95
									event_efficiency_errors[i] = 0.04

                                                        	if temp[i].pT() >= 1.7 and temp[i].pT() < 1.75:
                                                                	event_prt_efficiencies[i] = 1.0
									event_efficiency_errors[i] = 0

                                                        	if temp[i].pT() >= 1.75 and temp[i].pT() < 1.8:
                                                                	event_prt_efficiencies[i] = 0.97
									event_efficiency_errors[i] = 0.04

                                                        	if temp[i].pT() >= 1.8:
                                                                	event_prt_efficiencies[i] = 1.0
									event_efficiency_errors[i] = 0
							
								event_momentum_list[i] = math.sqrt( temp[i].px()**2 + temp[i].py()**2 + temp[i].pz()**2 )
								event_energy_list[i] = temp[i].e()
								event_prt_inv_masses[i] = ( ( event_energy_list[i]**2 ) - ( event_momentum_list[i]**2 ) )

							
							# Multiply event_contribution by particle efficiencies for event efficiency
                                                        event_efficiency = 1
                                                        for i in range( len( event_prt_efficiencies ) ):
                                                        	event_efficiency *= event_prt_efficiencies[i]
								event_fract_err += ( event_efficiency_errors[i] / event_prt_efficiencies[i] )**2
							
							

							# Add contribution of this K0sK0s event to K0sK0s event count
                                                        num_K0sK0s += event_efficiency
							# Calculate invariant mass of event
							event_inv_mass = sum( event_prt_inv_masses )
							event_inv_masses.append( event_inv_mass )
							
							# Write some data to sim_values.txt
							invmass_eff_err.write( str( event_inv_mass ) + ',' )
							invmass_eff_err.write( str( event_efficiency ) + ',' )
							invmass_eff_err.write( str( event_fract_err ) + '\n' )

							# Find min and max pT, max abs eta of current event's pions
               						min_pT = temp[0].pT()
               						max_pT = temp[0].pT()
                			 		max_abs_eta = abs(temp[0].eta())

                					for i in range(1, len(temp)):
                        					if temp[i].pT() < min_pT:
                                					min_pT = temp[i].pT()

                        					if temp[i].pT() > max_pT:
                               						max_pT = temp[i].pT()
			
                        					if abs(temp[i].eta()) > max_abs_eta:
                                					max_abs_eta = abs(temp[i].eta())

                					# Add values to files if they pass the pT and eta cuts
							# ADD TO .TXT FILES HERE
                					if ( min_pT >= 0.1 ) and ( max_abs_eta <= 2.5 ):
                        					min_max_cut1.write( str( min_pT ) + ',' )
                       						min_max_cut1.write( str( max_pT ) + ',' )
                        					min_max_cut1.write( str( max_abs_eta ) + '\n' )

							if ( min_pT >= 0.5 ) and ( max_abs_eta <= 2.5 ):
                                                                min_max_cut2.write( str( min_pT ) + ',' )
                                                                min_max_cut2.write( str( max_pT ) + ',' )
                                                                min_max_cut2.write( str( max_abs_eta ) + '\n' )

							if ( min_pT >= 1.0 ) and ( max_abs_eta <= 2.5 ):
                                                        	min_max_cut3.write( str( min_pT ) + ',' )
                                                                min_max_cut3.write( str( max_pT ) + ',' )
                                                                min_max_cut3.write( str( max_abs_eta ) + '\n' )

# Main loop, creates events
nEvent = pythia.mode("Main:numberOfEvents")

for iEvent in range(0, nEvent):
	if not pythia.next(): continue
	checkEvent(pythia.event)


# End of loop statistics
pythia.stat()
sigma_info = pythia.info.sigmaGen();
weightSum = pythia.info.weightSum();

print("Number of K0sK0s events: " + str(num_K0sK0s))
print("Length of K0sK0s_daugter_pions: " + str(len(K0sK0s_daughter_pions)))
print("Estimated cross section: " + str(sigma_info))
print("Weight Sum: " + str(weightSum))
print("(weight sum) * (sigma): " + str(sigma_info * weightSum))

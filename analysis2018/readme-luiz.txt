Analysis code 2018

Create program area:

source /cvmfs/cms.cern.ch/cmsset_default.sh

cmsrel CMSSW_10_6_18
cd CMSSW_10_6_18/src
cmsenv

to solve the python-CMSSW issue with the libraries do:
bash:
export PYTHONHOME=/cvmfs/cms.cern.ch/slc7_amd64_gcc820/external/python/2.7.15/
tcsh (ch):
source PYTHONHOME /cvmfs/cms.cern.ch/slc7_amd64_gcc820/external/python/2.7.15/

mkdir PromptLUAna; cd PromptLUAna 
mkedanlzr PromptAnalyzer
cd PromptAnalyzer/plugins
cp /afs/cern.ch/user/l/lregisem/CMSSW_10_6_18/src/PromptLUAna/PromptAnalyzer/plugins/BuildFile.xml .
cp /afs/cern.ch/user/l/lregisem/CMSSW_10_6_18/src/PromptLUAna/PromptAnalyzer/plugins/PromptAnalyzer.cc .
scram b

# copy my modified files
cd ../python/
cp /afs/cern.ch/user/l/lregisem/CMSSW_10_6_18/src/PromptLUAna/PromptAnalyzer/python/ConfFile_XXX_RECO_withDedx.py .
cp /afs/cern.ch/user/l/lregisem/CMSSW_10_6_18/src/PromptLUAna/PromptAnalyzer/python/submit-condorRECOall.csh .
cp /afs/cern.ch/user/l/lregisem/CMSSW_10_6_18/src/PromptLUAna/PromptAnalyzer/python/t* .
scram b

kinit
voms-proxy-init -voms cms
cp /tmp/x509up_u107606 ~/
The file above is user's unique so you need first to register
with VOMS admin for a membership and then get the GRID credentials.
Warning: the documentation may be obsolete so I will prepare a
step-by-step guide, there are many details in order to get the credentials working.

# change name to your user ID in the script: submit-condorRECOall.csh
# also change the X509 number; mine is u10760.

Ferenc's PID
package: UseCode_EnergyLossPID.tgz
deploy the PID package in /afs/cern.ch/user/l/lregisem/CMSSW_10_6_18/src/
or in your /afs/cern.ch/user/FIRST-LETTER/USERNAME/CMSSW_10_6_18/src/
do this:
cd /afs/cern.ch/user/FIRST-LETTER/USERNAME/CMSSW_10_6_18/src/
mkdir UserCode
cd UserCode/
git cms-init
git clone git@github.com:lgemedia/EnergyLossPID.git .
cd ..
scram b

Please, read this:
------------------------------------------------------------------------------

  For those dE/dx values and probabilities you need to apply a special package
on RECO data, so it is rather complicated. I give the details below:

  You need to add the package UserCode/EnergyLossPID. This package is not part
of official CMSSW, that is why it is under UserCode. I attach a recent tar, you
need to tar -xvzf that in src. It contains necessary code and dE/dx calibration
necessary. It is meant for CMSSW_10_X_X.

  Try to compile that.

  The scripts/ directory contains: addDeDx_conf.py and addDeDx_reco.py

  You would need four favourite json to add ('90m.json'), and modify
addDeDx_reco.py for your needs. It is not expected to work out-of-the-box, be
careful.

  In addDeDx_reco.py the lines
  process.reco = cms.Path(process.MeasurementTrackerEvent
                       * process.siPixelClusterShapeCache)
  take care of preparing the cluster info

  and the line
  process.produceEnergyLoss is the place where the dE/dx calculation is done

  If properly used, it will add *on the fly* vector<reco::DeDxData> structures
with labels "energyLossPixHits", "energyLossStrHits", "energyLossAllHits".

  You need RECO format with pixel/strip cluster info kept, for all that
calculations.

  Ferenc

------------------------------------------------------------------------------

Sometimes is good to clean the scram:
scram b clean; scram b

When fixing code issues it is good to include the debugger: 
scram b clean; scram b USER_CXXFLAGS="-g"

Submission: (t2 and t4 subdirectories have 2,000 job files with
two data sets in each to be run) 
------------
t2 = TB/BT
t4 = TT/BB

for each of the created dir: cd dir ; condor_submit submit_all

./submit-condorRECOall.csh t200 t20.eos_0 
./submit-condorRECOall.csh t201 t20.eos_1 
./submit-condorRECOall.csh t210 t21.eos_0
./submit-condorRECOall.csh t211 t21.eos_1
./submit-condorRECOall.csh t220 t22.eos_0
./submit-condorRECOall.csh t221 t22.eos_1
./submit-condorRECOall.csh t230 t23.eos_0
./submit-condorRECOall.csh t231 t23.eos_1

./submit-condorRECOall.csh t40 t40.eos
./submit-condorRECOall.csh t41 t41.eos
./submit-condorRECOall.csh t42 t42.eos
./submit-condorRECOall.csh t43 t43.eos

if you need to resubmit a particular job, do: condor_submit submit_[jobnr]


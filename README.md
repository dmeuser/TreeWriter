**TreeWriter** to build a ROOT tree from MiniAOD. Designed for ttbar, starting from https://github.com/cms-susy-photon-rwth-1b/TreeWriter

## Building and Running ##
Get CMSSW environment 10_2X

```
cmsrel CMSSW_10_2_11_patch1
cd CMSSW_10_2_11_patch1/src/
cmsenv
git cms-init
git cms-merge-topic cms-egamma:EgammaPostRecoTools
scram b -j7
git clone git@github.com:dmeuser/TreeWriter.git
scram b -j7
cd TreeWriter
```
Run the TreeWriter
- locally
```
voms-proxy-init -voms cms
cmsRun TreeWriter/python/runTreeWriter.py
```
- on the Grid using CRAB3
```
. /cvmfs/cms.cern.ch/crab3/crab.sh
cd crab
```
for a single dataset
```
crab submit -c crabConfig.py
```
for all datasets
```
python2 crabConfig.py
```

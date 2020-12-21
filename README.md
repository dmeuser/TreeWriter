**TreeWriter** to build a ROOT tree from MiniAOD. Designed for ttbar, starting from https://github.com/cms-susy-photon-rwth-1b/TreeWriter

## Building and Running ##
Get CMSSW environment 10_2X

```
cmsrel CMSSW_10_2_18
cd CMSSW_10_2_18/src/
cmsenv
git cms-init
git cms-merge-topic dmeuser:ttbar_tree
git cms-merge-topic yongbinfeng:DeepMET
scram b -j7
git clone git@github.com:cms-egamma/EgammaPostRecoTools.git  EgammaUser/EgammaPostRecoTools
cd  EgammaUser/EgammaPostRecoTools
git checkout master
cd -
scram b -j7
git clone git@github.com:dmeuser/TreeWriter.git
scram b -j7
cd TreeWriter
```
Create Pileup Histograms (Source different CMSSW version before, due to https://hypernews.cern.ch/HyperNews/CMS/get/luminosity/816/1/1.html or use manual fix described there)

```
make -C PUreweighting
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
for all datasets
```
python crabConfig.py
```
Download,resubmitting and killing multiple crab jobs with
```
python multicrab.py
```

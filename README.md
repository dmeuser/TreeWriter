**TreeWriter** to build a ROOT tree from MiniAOD. Designed for ttbar, starting from https://github.com/cms-susy-photon-rwth-1b/TreeWriter

## Building and Running ##
Get CMSSW environment 10_6_X

```
cmsrel CMSSW_10_6_20
cd CMSSW_10_6_20/src/
cmsenv
git cms-init
git cms-merge-topic dmeuser:ttbar_treeUL
scram b -j7
git cms-addpkg RecoEgamma/EgammaTools  ### essentially just checkout the package from CMSSW
git clone https://github.com/cms-egamma/EgammaPostRecoTools.git
mv EgammaPostRecoTools/python/EgammaPostRecoTools.py RecoEgamma/EgammaTools/python/.
git clone -b ULSSfiles_correctScaleSysMC https://github.com/jainshilpi/EgammaAnalysis-ElectronTools.git EgammaAnalysis/ElectronTools/data/
git cms-addpkg EgammaAnalysis/ElectronTools
scram b -j7
cd TopQuarkAnalysis
git clone https://gitlab.cern.ch/CMS-TOPPAG/BFragmentationAnalyzer.git
scram b -j7
cd ..
git clone git@github.com:dmeuser/TreeWriter.git -b UltraLegacy
scram b -j7
cd TreeWriter
```
Create Pileup Histograms

```
make -C PUreweighting
```
Run the TreeWriter
- locally
```
voms-proxy-init -voms cms
cmsRun TreeWriter/python/runTreeWriter2018.py
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

**TreeWriter** to build a ROOT tree from MiniAOD. Designed for ttbar, starting from https://github.com/cms-susy-photon-rwth-1b/TreeWriter

## Building ##
Get CMSSW environment 10_6_X and use the following recipe to setup the TreeWriter

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
## Pileup Histograms ##
To produce nTuples histograms with pileup weights are needed. Usually they are already stored in the repository (e.g. for UL17 and UL18) in `PUreweighting/data/puWeights.root`. Before producing or updating the weights, make sure that the JSONS and pileup distributions in `makefile` are up to date and correctly selected for a given year. In addition the correct pileup mix has to be selected in `createMChist.py`. The the pileup weights can be produced and updated using
```
make -C PUreweighting
```

## Running the TreeWriter locally (testing) ##
The TreeWriter can be tested localy by running
```
voms-proxy-init -voms cms
cmsRun TreeWriter/python/runTreeWriter2018.py
```
This usually takes one single input file defined in `runTreeWriter2018.py` and produces the corresponding nTuple. For testing it is usually recommended to not run on the total input file but set `options.maxEvents=<numberOfEvents>` (keep in mind to reset this number to `-1` when running on crab).

## Running the TreeWriter on crab ##
The actualy production of nTuples is done by running the TreeWriter on crab. All necessary configuration can be found in `crabConfig_<year>.py`. Here the input samples are defined, as well as other crab settings like `splitting` or the output path on dCache. For submitting the jobs denfined in the Config use
```
. /cvmfs/cms.cern.ch/crab3/crab.sh
cd crab
python crabConfig_2018.py
```

## Handling crab jobs ##
Running, killed or finished crab jobs can be handled by the tool found in `crab/crabTools/`. The main command to be used is
```
python multicrab.py
```
By setting different options the following tasks can be performed:
- Checking the status of jobs (would rather prefer to use Grafana for this)
- Resubmitting failed jobs (with adapted `maxjobruntime` or adapted `siteblacklist`)
- Kill individual or all jobs
- Download finished jobs. Here the `maxMergeSize` can be set in order to create only nTuples with a given size (recommended und currently set to 10 GB)
- ...

To see the full list of options use `python multicrab.py --help`.

## Adaptions for new users ##
Since the TreeWriter is currently only used by one person, there are some changes which need to be applied in order for other users to use the TreeWriter. There is no guarantee for the completeness of the following list:

`crabConfig<year>.py`:
- Output path on dCache defined by `config.Data.outLFNDirBase`

`crabInfo.py`:
- User related storage settings in `getOutFileName`

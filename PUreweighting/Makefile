#doku in https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData
# JSON: 36.459/fb https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/2766.html
# see https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData
#2016
ANA_JSON2016=/home/home4/institut_1b/teroerde/cmssw/CMSSW_10_6_20/src/TreeWriter/PUreweighting/JSONS/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt
PU_JSON2016=/home/home4/institut_1b/teroerde/cmssw/CMSSW_10_6_20/src/TreeWriter/PUreweighting/JSONS/pileup_latest_2016.txt
MC_TRUEVERTEX2016=data/nTrueVertexMC_2016.root
DATA_TRUEVERTEX2016=data/PileupHistogram-goldenJSON-13tev-2016-69200ub-99bins.root
DATA_TRUEVERTEX2016_UP=data/PileupHistogram-goldenJSON-13tev-2016-72400ub-99bins.root
DATA_TRUEVERTEX2016_DOWN=data/PileupHistogram-goldenJSON-13tev-2016-66000ub-99bins.root
#2016_PREVFP
ANA_JSON2016_PREVFP=/home/home4/institut_1b/teroerde/cmssw/CMSSW_10_6_20/src/TreeWriter/PUreweighting/JSONS/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt
PU_JSON2016=/home/home4/institut_1b/teroerde/cmssw/CMSSW_10_6_20/src/TreeWriter/PUreweighting/JSONS/pileup_latest_2016_preVFP.txt
MC_TRUEVERTEX2016_PREVFP=data/nTrueVertexMC_2016.root
DATA_TRUEVERTEX2016_PREVFP=data/PileupHistogram-goldenJSON-13tev-2016-preVFP-69200ub-99bins.root
DATA_TRUEVERTEX2016_PREVFP_UP=data/PileupHistogram-goldenJSON-13tev-2016-preVFP-72400ub-99bins.root
DATA_TRUEVERTEX2016_PREVFP_DOWN=data/PileupHistogram-goldenJSON-13tev-2016-preVFP-66000ub-99bins.root
#2016_POSTVFP
ANA_JSON2016_POSTVFP=/home/home4/institut_1b/teroerde/cmssw/CMSSW_10_6_20/src/TreeWriter/PUreweighting/JSONS/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt
PU_JSON2016=/home/home4/institut_1b/teroerde/cmssw/CMSSW_10_6_20/src/TreeWriter/PUreweighting/JSONS/pileup_latest_2016_postVFP.txt
MC_TRUEVERTEX2016_POSTVFP=data/nTrueVertexMC_2016.root
DATA_TRUEVERTEX2016_POSTVFP=data/PileupHistogram-goldenJSON-13tev-2016-postVFP-69200ub-99bins.root
DATA_TRUEVERTEX2016_POSTVFP_UP=data/PileupHistogram-goldenJSON-13tev-2016-postVFP-72400ub-99bins.root
DATA_TRUEVERTEX2016_POSTVFP_DOWN=data/PileupHistogram-goldenJSON-13tev-2016-postVFP-66000ub-99bins.root
#2017
ANA_JSON2017=/home/home4/institut_1b/dmeuser/CMSSW_ttbarTreeUL/CMSSW_10_6_29/src/TreeWriter/PUreweighting/JSONS/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt
PU_JSON2017=/home/home4/institut_1b/dmeuser/CMSSW_ttbarTreeUL/CMSSW_10_6_29/src/TreeWriter/PUreweighting/JSONS/pileup_latest_2917.txt
MC_TRUEVERTEX2017=data/nTrueVertexMC_2917.root
DATA_TRUEVERTEX2017=data/PileupHistogram-goldenJSON-13tev-2017-69200ub-99bins.root
DATA_TRUEVERTEX2017_UP=data/PileupHistogram-goldenJSON-13tev-2017-72400ub-99bins.root
DATA_TRUEVERTEX2017_DOWN=data/PileupHistogram-goldenJSON-13tev-2017-66000ub-99bins.root
#2018
ANA_JSON2018=/home/home4/institut_1b/dmeuser/CMSSW_ttbarTreeUL/CMSSW_10_6_29/src/TreeWriter/PUreweighting/JSONS/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt
PU_JSON2018=/home/home4/institut_1b/dmeuser/CMSSW_ttbarTreeUL/CMSSW_10_6_29/src/TreeWriter/PUreweighting/JSONS/pileup_latest_2018.txt
MC_TRUEVERTEX2018=data/nTrueVertexMC_2018.root
DATA_TRUEVERTEX2018=data/PileupHistogram-goldenJSON-13tev-2018-69200ub-99bins.root
DATA_TRUEVERTEX2018_UP=data/PileupHistogram-goldenJSON-13tev-2018-72400ub-99bins.root
DATA_TRUEVERTEX2018_DOWN=data/PileupHistogram-goldenJSON-13tev-2018-66000ub-99bins.root

#Set Year to process
ANA_JSON=$(ANA_JSON2017)
PU_JSON=$(PU_JSON2017)
MC_TRUEVERTEX=$(MC_TRUEVERTEX2017)
DATA_TRUEVERTEX=$(DATA_TRUEVERTEX2017)
DATA_TRUEVERTEX_UP=$(DATA_TRUEVERTEX2017_UP)
DATA_TRUEVERTEX_DOWN=$(DATA_TRUEVERTEX2017_DOWN)


all: data/puWeights.root

data:
	mkdir -p $@

data/puWeights.root: data createWeights.py $(MC_TRUEVERTEX) $(DATA_TRUEVERTEX) $(DATA_TRUEVERTEX_UP) $(DATA_TRUEVERTEX_DOWN)
	python2.7 createWeights.py --mc $(MC_TRUEVERTEX) --data $(DATA_TRUEVERTEX) --data_up $(DATA_TRUEVERTEX_UP) --data_down $(DATA_TRUEVERTEX_DOWN) -o $@

$(MC_TRUEVERTEX): createMChist.py
	python2.7 $^ $@


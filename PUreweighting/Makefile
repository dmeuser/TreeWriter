#doku in https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData
# JSON: 36.459/fb https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/2766.html
#2016
ANA_JSON=/home/home4/institut_1b/dmeuser/CMSSW_ttbar/CMSSW_10_2_11_patch1/src/TreeWriter/PUreweighting/JSONS/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt
PU_JSON=/home/home4/institut_1b/dmeuser/CMSSW_ttbar/CMSSW_10_2_11_patch1/src/TreeWriter/PUreweighting/JSONS/pileup_latest_2016.txt
MC_TRUEVERTEX=nTrueVertexMC_2016.root
DATA_TRUEVERTEX=nTrueVertexData_2016.root
#2017
# ~ANA_JSON=/home/home4/institut_1b/dmeuser/CMSSW_ttbar/CMSSW_10_2_11_patch1/src/TreeWriter/PUreweighting/JSONS/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt
# ~PU_JSON=/home/home4/institut_1b/dmeuser/CMSSW_ttbar/CMSSW_10_2_11_patch1/src/TreeWriter/PUreweighting/JSONS/pileup_latest_2017.txt
# ~MC_TRUEVERTEX=nTrueVertexMC_2017.root
# ~DATA_TRUEVERTEX=nTrueVertexData_2017.root
#2018
#~ ANA_JSON=/home/home4/institut_1b/dmeuser/CMSSW_ttbar/CMSSW_10_2_11_patch1/src/TreeWriter/PUreweighting/JSONS/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt
#~ PU_JSON=/home/home4/institut_1b/dmeuser/CMSSW_ttbar/CMSSW_10_2_11_patch1/src/TreeWriter/PUreweighting/JSONS/pileup_latest_2018.txt
#~ MC_TRUEVERTEX=nTrueVertexMC_2018.root
#~ DATA_TRUEVERTEX=nTrueVertexData_2018.root



PU_CALC=pileupCalc.py --calcMode true --maxPileupBin 100 --numPileupBins 100 -i $(ANA_JSON) --inputLumiJSON $(PU_JSON)

all: data/puWeights.root

data:
	mkdir -p $@

data/puWeights.root: data createWeights.py $(MC_TRUEVERTEX) $(DATA_TRUEVERTEX)
	python2.7 createWeights.py --mc $(MC_TRUEVERTEX) --data $(DATA_TRUEVERTEX) -o $@

$(MC_TRUEVERTEX): createMChist.py
	python2.7 $^ $@

# For the cross section and its uncertainty, see https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData
$(DATA_TRUEVERTEX): $(ANA_JSON) $(PU_JSON) Makefile
	$(PU_CALC) --minBiasXsec 69200 --pileupHistName pileup tmpHist1.root
	$(PU_CALC) --minBiasXsec 72383 --pileupHistName pileupUp tmpHist2.root
	$(PU_CALC) --minBiasXsec 66017 --pileupHistName pileupDown tmpHist3.root
	hadd -f $@ tmpHist1.root tmpHist2.root tmpHist3.root
	rm tmpHist1.root tmpHist2.root tmpHist3.root


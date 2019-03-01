#doku in https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData
# JSON: 36.459/fb https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/2766.html
ANA_JSON=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt
PU_JSON=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt

PU_CALC=pileupCalc.py --calcMode true --maxPileupBin 75 --numPileupBins 75 -i $(ANA_JSON) --inputLumiJSON $(PU_JSON)

all: data/puWeights.root

data:
	mkdir -p $@

data/puWeights.root: data createWeights.py nTrueVertexMC.root nTrueVertexData.root
	python2.7 createWeights.py --mc nTrueVertexMC.root --data nTrueVertexData.root -o $@

nTrueVertexMC.root: createMChist.py
	python2.7 $^ $@

# For the cross section and its uncertainty, see https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData
nTrueVertexData.root: $(ANA_JSON) $(PU_JSON) Makefile
	$(PU_CALC) --minBiasXsec 69200 --pileupHistName pileup tmpHist1.root
	$(PU_CALC) --minBiasXsec 72383 --pileupHistName pileupUp tmpHist2.root
	$(PU_CALC) --minBiasXsec 66017 --pileupHistName pileupDown tmpHist3.root
	hadd -f $@ tmpHist1.root tmpHist2.root tmpHist3.root
	rm tmpHist1.root tmpHist2.root tmpHist3.root


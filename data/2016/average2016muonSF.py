#Script to average 2016 muons SFs for the different eras
import ROOT


# all int lumis 2016 per run range in fb-1:
B = 5.750490644
C = 2.572903489
D = 4.242291557
E = 4.025228137
F = 3.104509132
G = 7.575824256
H = 8.650628380
BtoF = B + C + D + E + F
GtoH = G + H


def combine2DHistos(hist1, hist2, weight1, weight2):
    outHisto = hist1.Clone()
    for binX in range(0, hist1.GetNbinsX() + 1):
        for binY in range(0, hist1.GetNbinsY() + 1):
            v1 = hist1.GetBinContent(binX, binY)
            v2 = hist2.GetBinContent(binX, binY)
            v1err = hist1.GetBinError(binX, binY)
            v2err = hist2.GetBinError(binX, binY)
            # v1errUp=hist1.GetBinErrorUp(binX,binY)
            # v1errDn=hist1.GetBinErrorUp(binX,binY)
            # v2errUp=hist2.GetBinErrorUp(binX,binY)
            # v2errDn=hist2.GetBinErrorUp(binX,binY)
            outV = (weight1 * v1 + weight2 * v2) / (weight1 + weight2)
            term1 = weight1 / (weight1 + weight2) * v1err
            term2 = weight2 / (weight1 + weight2) * v2err
            outErr = ROOT.TMath.Sqrt(
                (term1 * term1) + (term2 * term2))
            outHisto.SetBinContent(binX, binY, outV)
            outHisto.SetBinError(binX, binY, outErr)
    return outHisto.Clone()

f = ROOT.TFile.Open('RunBCDEF_SF_ID.root', 'read')
f2 = ROOT.TFile.Open('RunGH_SF_ID.root', 'read')

hist_out=combine2DHistos(f.Get("NUM_TightID_DEN_genTracks_eta_pt"),f2.Get("NUM_TightID_DEN_genTracks_eta_pt"),BtoF,GtoH)
hist_out.SaveAs("MuonSF_ID_merged.root")

f = ROOT.TFile.Open('RunBCDEF_SF_ISO.root', 'read')
f2 = ROOT.TFile.Open('RunGH_SF_ISO.root', 'read')

hist_out=combine2DHistos(f.Get("NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt"),f2.Get("NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt"),BtoF,GtoH)
hist_out.SaveAs("MuonSF_ISO_merged.root")

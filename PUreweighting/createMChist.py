import ROOT
import sys

def generateHistoFromList( list, name, title="pileup" ):
    while len(list)<75: list += [0] # same binning as for data needed
    histo = ROOT.TH1F( name, title, len(list), 0, len(list) )
    for bin, content in enumerate( list ):
        histo.SetBinContent( bin+1, content ) # bin0 is underflow!
    return histo

def generateHistoFromMixingModule( mixing ):
    exec( "from SimGeneral.MixingModule.{}_cfi import mix".format( mixing ) )
    puList = list(mix.input.nbPileupEvents.probValue)
    return generateHistoFromList( puList, mixing )

def writeObjectsToFile( list, filename ):
    # list contains objects which will be written to a TFile
    outfile = ROOT.TFile( filename, "recreate" )
    outfile.cd()
    for entry in list:
        entry.Write()
    outfile.Close()


if __name__ == "__main__":
    # nicer to view this way (but errors are not used anyway):
    ROOT.TH1.SetDefaultSumw2(True)

    import sys
    if len(sys.argv) < 2:
        sys.exit('Usage: %s outputFilename.root' % sys.argv[0])
    filename = sys.argv[1]

    hists = []
    hists.append( generateHistoFromMixingModule( "mix_2015_25ns_Startup_PoissonOOTPU" ) )
    hists.append( generateHistoFromMixingModule( "mix_2015_25ns_FallMC_matchData_PoissonOOTPU" ) )
    hists.append( generateHistoFromMixingModule( "mix_2016_25ns_SpringMC_PUScenarioV1_PoissonOOTPU" ) )
    hists.append( generateHistoFromMixingModule( "mix_2016_25ns_Moriond17MC_PoissonOOTPU" ) )

    writeObjectsToFile( hists, filename )



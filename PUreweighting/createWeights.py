import argparse
import ROOT

def readHisto( fileName, histoname, newName="" ):

    f = ROOT.TFile( fileName )
    ROOT.gROOT.cd()

    # This histogram is still associated with the file and will be lost after
    # the end of the function.
    histoInFile = f.Get( histoname ).Clone()
    return histoInFile


def getObjectNames( filename ):
    names = []
    f = ROOT.TFile( filename )
    for k in f.GetListOfKeys():
        if isinstance( k.ReadObj(), ROOT.TH1 ):
            names.append( k.GetName() )
    return names


if __name__ == "__main__":

    arguments = argparse.ArgumentParser( description="Print out filters" )
    arguments.add_argument( "-o", "--outputFilename", default="puWeights.root" )
    arguments.add_argument( "--mc", default="mc.root"  )
    arguments.add_argument( "--data", default="data.root" )
    opts = arguments.parse_args()

    dataHist = readHisto( opts.data, "pileup" )
    dataHistUp = readHisto( opts.data, "pileupUp" )
    dataHistDown = readHisto( opts.data, "pileupDown" )

    for h in [ dataHist, dataHistUp, dataHistDown ]:
        h.Scale( 1./h.Integral() )


    mcHistNames = getObjectNames( opts.mc )

    weightHists = []

    for hName in mcHistNames:
        mcH = readHisto( opts.mc, hName )
        mcH.Scale( 1./mcH.Integral() )

        pileupWeight = dataHist.Clone( "pileupWeight_"+hName )
        pileupWeightUp = dataHistUp.Clone( "pileupWeightUp_"+hName )
        pileupWeightDown = dataHistDown.Clone( "pileupWeightDown_"+hName )

        for h in pileupWeight, pileupWeightUp, pileupWeightDown:
            h.Divide( mcH )

        weightHists.extend( [pileupWeight, pileupWeightUp, pileupWeightDown] )

    from createMChist import writeObjectsToFile
    writeObjectsToFile( weightHists, opts.outputFilename )

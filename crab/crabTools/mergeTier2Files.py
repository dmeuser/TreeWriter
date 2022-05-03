#!/usr/bin/env python2
"""
Script to merge all files of a directory on Tier2.
usage ./mergeTier2Files.py <outputFile>.root <srm-source-path>
- <srm-source-path> has to be of the form:
  "srm://grid-srm.physik.rwth-aachen.de:8443/srm/managerv2?SFN=/pnfs/physik.rwth-aachen.de/cms/store/user/<...>/<date>_<time>/"
  the "<date>_<time>" directory contains numbered directries "0001","0002",... that contain the root files
example usage:
./mergeTier2Files.py out.root srm://grid-srm.physik.rwth-aachen.de:8443/srm/managerv2?SFN=/pnfs/physik.rwth-aachen.de/cms/store/user/jolange/data/Test/GJets_HT-400to600_Tune4C_13TeV-madgraph-tauola/sizeCheck/150504_120017/
"""
import subprocess as sp
import sys
import argparse
import os
import os.path
import multiprocessing
import ROOT
import glob

def openDcacheFile(fname):
    f = ROOT.TFile.Open("dcap://grid-dcap-extern.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms"+fname)

def isCorrupt(fname):
    p = multiprocessing.Process(target=openDcacheFile, name="Foo", args=(fname,))
    p.start()
    p.join(2)
    if p.is_alive():
        p.terminate()
        p.join()
        return True
    return False

def getMergingScheme(localfiles, maxMergeSize):
    scheme = {}
    currentSize = 0
    currentList = []
    for filePath in localfiles:
        fileSize = os.path.getsize(filePath)
        currentSize += fileSize
        if((currentSize*1e-9)<maxMergeSize):
           currentList.append(filePath)
        else:
            scheme[len(scheme)+1]=currentList
            currentList = [filePath]
            currentSize = fileSize
    scheme[len(scheme)+1]=currentList
    return scheme

def mergeFiles(inputFiles,outputFile):
    """
    - uses "hadd" to add all inputFiles
    - prepends xrootd-prefixes for remote access
    - merge result is written to outputFile
    """
    print "corrupt files:"
    ROOT.gErrorIgnoreLevel = ROOT.kError
    uncorruptFiles = []
    for f in inputFiles:
        if isCorrupt(f): print f
        else: uncorruptFiles.append(f)
    ROOT.gErrorIgnoreLevel = ROOT.kInfo

    if len(uncorruptFiles) != len(inputFiles): return False
    # warning: do use with care:
    if True: inputFiles = uncorruptFiles

    # add prefix to access remotely
    gridPrefix="root://xrootd-cms.infn.it//"
    inputFiles = [gridPrefix+f for f in inputFiles]
    print "using hadd to merge",len(inputFiles),"files..."
    if sp.call(["hadd","-f",outputFile]+inputFiles):
        sys.exit(1)
    print "written",outputFile
    return True

def getDirectoryContent(srmDirectoryPath):
    """
    - returns a list of all files and directories contained in a directory
    - srmDirectoryPath has to be a srm-syntax path
    - returned paths start from '/pnfs/'
    """
    # get srmls output
    contents=sp.check_output(["srmls",srmDirectoryPath])
    # separate output at spaces/line breaks
    contents=contents.split()
    # only every second entry is name, the rest are sizes
    # and the very first entry is the directory name
    contents=contents[3::2]
    # ignore the files placed in the subfolders not containing the
    # relevant data files and the folders themselves
    ignores=("/failed/","/log/")
    contents=[f for f in contents if not f.endswith(ignores)]
    return contents

def getFilePaths(srmDirectoryPath):
    """
    - returns a list of all files contained in a directory
    - srmDirectoryPath has to be a srm-syntax path
    - returned paths start from '/store/'
    """
    files=getDirectoryContent(srmDirectoryPath)
    # extract the relevant part of the paths
    files= ["/store/"+f.partition("/cms/store/")[-1] for f in files if f.endswith(".root")]
    return files

def downloadAndMergeFiles(inputFiles, outputFile, maxMergeSize=-1):
    tmpDownloadDir = outputFile.replace(".root", "")
    if not os.path.isdir(tmpDownloadDir): os.mkdir(tmpDownloadDir)
    for ifile, f in enumerate(inputFiles):
        if not os.path.isfile(os.path.join(tmpDownloadDir, os.path.basename(f))):
            print "Downloading {} {}/{}".format(f, ifile+1, len(inputFiles))
            sp.call(["srmcp", "srm://grid-srm.physik.rwth-aachen.de:8443/srm/managerv2?SFN=/pnfs/physik.rwth-aachen.de/cms"+f, "file:///{}/{}".format(tmpDownloadDir, f.split("/")[-1])])
        else:
            print "File {} already in folder".format(f)
    localFiles = [f for f in glob.glob("{}/*root".format(tmpDownloadDir))]
    if len(localFiles) == len(inputFiles):
        if maxMergeSize==-1:    # merge all files
            if sp.call(["hadd","-f",outputFile]+localFiles):
                sys.exit(1)
        else:
            mergingScheme = getMergingScheme(localFiles,maxMergeSize)
            for outputNr in mergingScheme:
                outputFile=outputFile.split(".root")[0]
                if sp.call(["hadd","-f",outputFile+"_"+str(outputNr)+".root"]+mergingScheme[outputNr]):
                    sys.exit(1)
        print "Remove temporary files"
        for f in localFiles:
            os.remove(f)
        os.rmdir(tmpDownloadDir)
        return True
    else:
        print "Do not merge, since not all files downloaded"
        return False

def mergeTier2Files( outputFilePath, inputFilePath, checkDuplicates=False, downloadFirst=False, maxMergeSize=-1 ):
    # get all the subdirectories "/XXXX/" that contain the root files
    dataDirectories=getDirectoryContent(inputFilePath)
    # find all files in these subdirectories
    inputFiles=[]
    for d in dataDirectories:
        inputFiles+=getFilePaths("srm://grid-srm.physik.rwth-aachen.de:8443/srm/managerv2?SFN="+d)
    inputFiles.sort()

    # merge all of them
    if downloadFirst:
        out = downloadAndMergeFiles(inputFiles, outputFilePath, maxMergeSize)
    else:
        out = mergeFiles(inputFiles,outputFilePath)

    if checkDuplicates:
        # check if duplicate events exist
        import DuplicateEventFilter.DuplicateEventFilter as dupFilter
        dupFilter.filterFile(outputFilePath,"TreeWriter/eventTree")
    return out

if __name__=="__main__":

    parser = argparse.ArgumentParser(description="Script to merge all files of a directory on Tier2.")
    parser.add_argument("outFile")
    parser.add_argument("srm_source_path")
    parser.add_argument("-n", "--noDuplicateCheck", action="store_true")
    parser.add_argument("-d", "--downloadFirst", action="store_true", help="First download all files, and then merge")
    args = parser.parse_args()

    mergeTier2Files( args.outFile, args.srm_source_path, not args.noDuplicateCheck, args.downloadFirst )

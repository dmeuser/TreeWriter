#!/usr/bin/env python2

import re
import getpass
import os, subprocess as sp
from CRABClient.ClientUtilities import colors
import mergeTier2Files

def myMatch( regex, string ):
    m = re.match( regex, string )
    return [ m ] if m else []

def modifyDatasetName( dataset ):
    # shorten dataset name
    deletes = [
        "_13TeV",
        "-madgraph",
        "_Tune4C",
        "-tauola",
        "_MSDecaysCKM_central",
        "_TuneCUETP8M1",
        "-amcatnloFXFX",
        "-madspin",
        "-pythia8",
        "MLM",
        "_pythia8",
        "_Pythia8",
        "-powheg",
        "_TuneCUETP8M2T4"
    ]
    if dataset.startswith("WGToLNuG_") or dataset.startswith("TTJets_Tune") or dataset.startswith("WJetsToLNu_") or dataset.startswith("DYJetsToLL_"):
        # to be able to seperate these two datasets
        deletes.remove( "-amcatnloFXFX" )
        deletes.remove( "-madgraph" )
        deletes.remove( "MLM" )
    for d in deletes:
        dataset = dataset.replace( d, "" )
    return dataset

class CrabInfo:

    def __init__( self, infoString ):
        self.user=getpass.getuser()
        if infoString.endswith("/crab.log"):
            self.initFromLog( infoString )
        elif "/pnfs/physik.rwth-aachen.de/cms/store" in infoString:
            self.initFromSrm( infoString )

    def initFromLog( self, logFileName ):
        self.logFileDir = logFileName.replace("/crab.log","")

        with open( logFileName ) as f:
            lines = f.readlines()

        self.validStatusCache=False
        for line in lines:
            for m in myMatch("config.Data.outLFNDirBase = '(.*)'", line ):
                self.outLFNDirBase = m.group(1)
            for m in myMatch( ".*Task name:.*\s([\d_]+):.*_crab_(.*)", line ):
                self.time = m.group(1)
            for m in myMatch( "config.Data.outputDatasetTag = '(.*)'", line ):
                self.outputDatasetTag = m.group(1)
            for m in myMatch( "config.Data.inputDataset = '(.*)'", line ):
                self.inputDataset = m.group(1)
                self.datasetName, self.datasetMiddle, self.datasetType = self.inputDataset.split("/")[1:]
            for m in myMatch( ".*Got information from status cache file: (.*)", line ):
                self.details = eval( m.group(1) )
                self.validStatusCache=True;
            for m in myMatch( ".*Status on the CRAB server:(.*)", line ):
                self.statusCRAB = m.group(1).strip()
            for m in myMatch( ".*Status on the scheduler:(.*)", line ):
                self.statusScheduler = m.group(1).strip()

        if self.validStatusCache==True and self.statusCRAB != "SUBMITFAILED":
           self.nJobs = len(self.details)
           self.jobStates = [x["State"] for x in self.details.values()]
        else:
           self.statusScheduler="Not available since submit failed or not yet finished"

    def initFromSrm( self, srmPath ):
        self.srmPath = srmPath
        srmPathSplitted = srmPath.split("/")
        self.time = srmPathSplitted[-2]
        self.outputDatasetTag = srmPathSplitted[-3]
        self.datasetName = srmPathSplitted[-4]
        self.datasetMiddle="" # cannot easily find this, but crashes when undefined

    def getOutFileName( self ):
        #~ modifiedDatasetName = modifyDatasetName( self.datasetName )
        modifiedDatasetName = self.datasetName
        if "ext1" in self.datasetMiddle: modifiedDatasetName += "_ext"
        elif "ext2" in self.datasetMiddle: modifiedDatasetName += "_ext2"
        elif "ext3" in self.datasetMiddle: modifiedDatasetName += "_ext3"
        elif "ext4" in self.datasetMiddle: modifiedDatasetName += "_ext4"
        elif "ext5" in self.datasetMiddle: modifiedDatasetName += "_ext5"
        elif "backup" in self.datasetMiddle: modifiedDatasetName += "_backup"
        if self.user == "kiesel":
            if hasattr(self, "datasetType"):
                if self.datasetType == "MINIAOD": # data
                    modifiedDatasetName += "_"+self.datasetMiddle
            baseDir = "/user/kiesel/nTuples/"
            baseDir = "/net/scratch_cms1b2/cms/user/kiesel/"
            baseDir = os.path.join(baseDir, self.outputDatasetTag)
            if not os.path.isdir(baseDir): os.mkdir(baseDir)
            return os.path.join(baseDir, modifiedDatasetName+"_nTuple.root")
        elif self.user == "jschulz":
            if "Run2015" in self.datasetMiddle: # distinguish different data runs/recos
                modifiedDatasetName+="_"+self.datasetMiddle
            elif "Run2016" in self.datasetMiddle: # distinguish different data runs/recos
                modifiedDatasetName+="_"+self.datasetMiddle
            elif "T5gg" in self.datasetName: # extract mass
                m=re.search(".*_mGluino-(.*)_mNeutralino-(.*)-.*",self.datasetMiddle)
                if m and len(m.groups())==2:
                    modifiedDatasetName+="_g%s_n%s"%m.groups()
                else:
                    modifiedDatasetName="UNKOWNPATTERN"
            return "/user/jschulz/2016/data/run2/dl/{}.root".format(modifiedDatasetName)
        elif self.user == "swuchterl":
            if hasattr(self, "datasetType"):
                if self.datasetType == "MINIAOD": # data
                    modifiedDatasetName += "_"+self.datasetMiddle
            baseDir = "/user/swuchterl/nTuples/"
            baseDir = "/net/scratch_cms1b2/cms/user/swuchterl/"
            baseDir = os.path.join(baseDir, self.outputDatasetTag)
            if not os.path.isdir(baseDir): os.mkdir(baseDir)
            return os.path.join(baseDir, modifiedDatasetName+"_nTuple.root")
        elif self.user == "dmeuser":
            if "Run2015" in self.datasetMiddle: # distinguish different data runs/recos
                modifiedDatasetName+="_"+self.datasetMiddle
            elif "Run2016" in self.datasetMiddle: # distinguish different data runs/recos
                modifiedDatasetName+="_"+self.datasetMiddle
            elif "T5gg" in self.datasetName: # extract mass
                m=re.search(".*_mGluino-(.*)_mNeutralino-(.*)-.*",self.datasetMiddle)
                if m and len(m.groups())==2:
                    modifiedDatasetName+="_g%s_n%s"%m.groups()
                else:
                    modifiedDatasetName="UNKOWNPATTERN"
            if "RunIISummer16" in self.datasetMiddle or "Run2016" in self.datasetMiddle:
                return "/net/data_cms1b/user/dmeuser/top_analysis/2016/v22/{}.root".format(modifiedDatasetName)
            elif "Run2017" in self.datasetMiddle:
                return "/net/data_cms1b/user/dmeuser/top_analysis/2017/v02/{}.root".format(modifiedDatasetName)
            elif "Run2018" in self.datasetMiddle:
                return "/net/data_cms1b/user/dmeuser/top_analysis/2018/v02/{}.root".format(modifiedDatasetName)
        return "outputFile.root"



    def getSrmPath( self ):
        if not hasattr(self, 'srmPath'):
             self.srmPath = "/pnfs/physik.rwth-aachen.de/cms{}/{}/{}/{}/".format(self.outLFNDirBase,self.datasetName,self.outputDatasetTag,self.time)
        return self.srmPath

    def getSrmPathFull( self ):
        return "srm://grid-srm.physik.rwth-aachen.de:8443/srm/managerv2?SFN="+self.getSrmPath()

    def jobSummary( self ):
        jobsPerStatus = dict([(x,self.jobStates.count(x)) for x in list(set(self.jobStates))])
        for k,v in jobsPerStatus.iteritems():
            if not v: continue
            c = colors.NORMAL
            if k == "failed": c = colors.RED
            if k == "running": c = colors.GRAY
            if k == "transferring": c = colors.BLUE
            if k == "finished": c = colors.GREEN
            print "{}\t{}{}  \t{:.1%} ({}{:3}{}/{:3})".format(c,k,colors.NORMAL,1.*v/self.nJobs,c,v,colors.NORMAL,self.nJobs)

        """# TODO: nice error output
        for j, d in self.details["jobs"].iteritems():
            if d["State"] == "failed":
                for i in d["Error"]:
                    try:
                        if "code" in i: print i
                    except: pass
        """

    def completed(self):
        return self.statusScheduler == "COMPLETED"

    def getMergeCommand( self ):
        outFile = self.getOutFileName()
        srmSrc  = self.getSrmPathFull()
        return "./mergeTier2Files.py {} {}".format(outFile,srmSrc)

    def suggestMergeCommand(self):
        doneDir = self.doneDir()
        print "Suggested merge command:"
        print "  ",self.getMergeCommand()
        print "   mv {} {}".format( self.logFileDir, doneDir)

    def doneDir(self):
        """ Where the directory is moved when files have been downloaded
        """
        doneDir=self.logFileDir.replace("/crab_","/.{}_crab_".format(self.time))
        if self.user=="jschulz":
            origDir = self.logFileDir
            if origDir.endswith('/'): origDir=origDir[:-1]
            doneDir=origDir.split('/')[-1]
            doneDir = doneDir.replace("crab_","{}_crab_".format(self.outputDatasetTag))
            doneDir+="_"+self.time
            doneDir='/'.join(origDir.split('/')[:-1])+"/"+doneDir
        return doneDir
    
    def killedDir(self):
        """ Where the directory is moved when files have been killed
        """
        killedDir=self.logFileDir.replace("/crab_","/.killed/.{}_crab_".format(self.time))
        return killedDir

    def beautifyCrabStatus(self):
        print self.logFileDir
        if self.completed():
            print "{}COMPLETED!{}".format(colors.GREEN+colors.BOLD,colors.NORMAL)
        elif self.statusCRAB== "SUBMITFAILED":
           print colors.BOLD+colors.RED,
           print self.statusCRAB+colors.NORMAL
        elif self.validStatusCache==False:
		   print colors.BOLD+colors.GREEN,
		   print "SUBMIT NOT YET FINISHED"+colors.NORMAL
        else:
            if self.statusScheduler=="RESUBMITFAILED": print colors.BOLD+colors.RED,
            print self.statusScheduler+colors.NORMAL
            self.jobSummary()

    def moveCompleted(self):
        """ move completed directory
        """
        doneDir=self.doneDir()
        print "Moving crab directory to",doneDir
        os.rename(self.logFileDir, doneDir)
    
    def moveKilled(self):
        """ move killed directory
        """
        killedDir=self.killedDir()
        print "Moving crab directory to",killedDir
        os.rename(self.logFileDir, killedDir)

    def download(self, downloadFirst=False):
        """
        automatically download files belonging to the crab directory
        and rename the crab directory
        """
        print "Downloading to",self.getOutFileName()
        # change library path to cmssw default
        cmssw=os.environ['CMSSW_BASE']+"/src"
        cmsenv="eval `scramv1 runtime -sh`"
        crabLibPath=os.environ['LD_LIBRARY_PATH']
        cmsswLibPath=sp.check_output("cd "+cmssw+";"+cmsenv+"; echo $LD_LIBRARY_PATH",shell=True)

        os.environ['LD_LIBRARY_PATH']=cmsswLibPath
        sucess = mergeTier2Files.mergeTier2Files( self.getOutFileName(), self.getSrmPathFull(), downloadFirst=downloadFirst )
        # restore crabs library path
        os.environ['LD_LIBRARY_PATH']=crabLibPath
        if sucess: self.moveCompleted()
        return sucess

if __name__ == '__main__':
    import sys
    if len(sys.argv)<2:
        print "specify crab directory"
        exit(0)
    for dir in sys.argv[1:]:
        info = CrabInfo( dir+"/crab.log" )
        info.beautifyCrabStatus(False)
        print "Tier2 output path:\n  ",info.getSrmPathFull()
        print "Download command:\n  ",info.getMergeCommand()


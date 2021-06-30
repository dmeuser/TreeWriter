#!/usr/bin/env python2
import argparse
import subprocess
import os
import glob
import time
from CRABClient.ClientUtilities import colors

import crabInfo

def crabUpdate( dir ):
    # fist check if proxy existent
    out=""
    try:
        out=subprocess.check_output(["voms-proxy-info","--timeleft"],stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError, e:
        print "Initialize your VOMS proxy!"
        print e.output
        exit(0)
    if int(out)==0:
        print "Your VOMS proxy expired! Please refresh."
        exit(0)
    # try to update
    with open(os.devnull, "w") as FNULL:
        try:
            out=subprocess.check_output(["crab","status",dir],stdin=FNULL,stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            print(e.output)
        if "No credentials found!" in out:
            print "No credentials found! Initialize your VOMS proxy."
            exit(0)

def crabResubmit(directory,maxjobruntime,silent=False):
    with open(os.devnull, "w") as FNULL:
        if maxjobruntime<0:
            out=subprocess.check_output(["crab","resubmit",directory],stdin=FNULL,stderr=subprocess.STDOUT)
        else:
            out=subprocess.check_output(["crab","resubmit",directory,"--maxjobruntime="+str(maxjobruntime)],stdin=FNULL,stderr=subprocess.STDOUT)
        if "No credentials found!" in out:
            print "No credentials found! Initialize your VOMS proxy."
            exit(0)
        if not silent: print out

def crabKill(directory,silent=False):
    with open(os.devnull, "w") as FNULL:
        out=subprocess.check_output(["crab","kill",directory],stdin=FNULL,stderr=subprocess.STDOUT)
        if "No credentials found!" in out:
            print "No credentials found! Initialize your VOMS proxy."
            exit(0)
        if not silent: print out

def multicrab(args):
    dirs = args.dirs or glob.glob(os.environ['CMSSW_BASE']+'/src/TreeWriter/crab/crab_*/')

    iTotal=0
    iComplete=0
    killed=[];
    for dir in dirs:
        print dir
        iTotal+=1
        if not args.noUpdate: crabUpdate( dir )
        info = crabInfo.CrabInfo( dir+"/crab.log" )
        info.beautifyCrabStatus()
        if args.resubmit and info.validStatusCache==True:
			if "failed" in info.jobStates:
				print "Resubmitting..."
				crabResubmit(dir,args.maxjobruntime)
        elif args.kill:
            print "Killing..."
            killed.append(dir.split("/")[-2])
            crabKill(dir)
            info.moveKilled()
        elif args.killSubmitFailed and info.statusCRAB=="SUBMITFAILED":
            print "Killing..."
            killed.append(dir.split("/")[-2])
            info.moveKilled()
        else:
            if args.forceDL: info.download(args.downloadFirst)
            elif info.completed():
                iComplete+=1
                if args.autoDL:
                    info.download(args.downloadFirst)
                else:
                    try:
                        info.suggestMergeCommand()
                    except AttributeError,e:
                        print e
                if args.moveCompleted: info.moveCompleted()
        print

    print "==============================="
    print "Summary: %d/%d tasks completed"%(iComplete,iTotal)
    if args.kill or args.killSubmitFailed:
       print colors.BOLD+colors.RED,
       print "Killed the following samples:"+colors.NORMAL
       for sample in killed:
          print sample

    if args.repeat and iComplete != iTotal:
        print time.strftime("%H-%M-%S:"), "Sleeping for 30 minutes."
        time.sleep(1800)
        multicrab(args)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dirs', nargs='+', default=[] )
    parser.add_argument('--noUpdate', action='store_true' )
    parser.add_argument('--autoDL', action='store_true' )
    parser.add_argument('--forceDL', action='store_true' )
    parser.add_argument('--resubmit', action='store_true' )
    parser.add_argument("--maxjobruntime", type=int, default=-1)
    parser.add_argument('--moveCompleted', action='store_true' )
    parser.add_argument('--repeat', action='store_true' )
    parser.add_argument('--downloadFirst', action='store_true' )
    parser.add_argument('--kill', action='store_true' )
    parser.add_argument('--killSubmitFailed', action='store_true' )
    args = parser.parse_args()

    multicrab(args)

if __name__ == "__main__":
    main()

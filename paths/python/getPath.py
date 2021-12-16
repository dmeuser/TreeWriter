import FWCore.ParameterSet.Config as cms
import getpass
import os
import configparser as ConfigParser

def getPath(pathName, user=None):
    
    if user == None:
        user = getpass.getuser()
    
    if not os.path.isfile(os.path.join(os.path.dirname(__file__),"%s.ini"%user)):
        print "No path config found for user %s"%user
        exit(1)
    
    config = ConfigParser.RawConfigParser()
    config.read(os.path.join(os.path.dirname(__file__),"%s.ini"%user))
    
    if not "paths" in config:
        print "No section [paths] defined in config for user %s"%user
        exit(2)
    if not pathName in config["paths"]:
        print "Name %s not defined in config for user %s"%(pathName, user)
        exit(3)

    return config["paths"][pathName]

# user = None
# modifiedDatasetName = "asdf"
# print os.path.join(getPath("localDirBase", user),getPath("dir2016_preVFP", user),getPath("datasetTag_2016_preVFP", user),"nTuple/{}.root".format(modifiedDatasetName))

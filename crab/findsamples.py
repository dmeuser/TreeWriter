



samples = [
    "/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/MINIAODSIM",
    "/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/MINIAODSIM",
    "/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/MINIAODSIM",
    "/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v3/MINIAODSIM",
    "/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIISummer20UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8-v3/MINIAODSIM",
]

flag = "RunIISummer20UL16MiniAODAPVv2*"


import subprocess as sp

for sample in samples:
    name = sample.split("/")[1]
    out = sp.check_output('dasgoclient -query="dataset dataset=/%s/%s/MINIAOD* status=*"'%(name,flag), shell=True)

    if len(out) > 0:
        for line in out.split("\n"):
            if not len(line) > 0:
                continue
            line = line.strip()
            prodstatus = sp.check_output('dasgoclient -query="dataset dataset=%s status=*" -json'%(line), shell=True)

            if '"status":"PRODUCTION"' in prodstatus:
                print '"'+line+'",', "# production"
            else:
                if not '"status":"INVALID"' in prodstatus:
                    print '"'+line+'",'
    else:
        print "not found:", '/%s/%s/MINIAOD*'%(name,flag)
    

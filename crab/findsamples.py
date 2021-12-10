



samples = [
    "/TTTo2L2Nu_mtop169p5_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAOD-106X_mcRun2_asymptotic_v13-v2/MINIAODSIM",
    "/TTTo2L2Nu_mtop175p5_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAOD-106X_mcRun2_asymptotic_v13-v2/MINIAODSIM",
    "/TTToHadronic_mtop169p5_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAOD-106X_mcRun2_asymptotic_v13-v2/MINIAODSIM",
    "/TTToHadronic_mtop175p5_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAOD-106X_mcRun2_asymptotic_v13-v2/MINIAODSIM",
    "/TTToSemiLeptonic_mtop169p5_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAOD-106X_mcRun2_asymptotic_v13-v1/MINIAODSIM",
    "/TTToSemiLeptonic_mtop175p5_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16MiniAOD-106X_mcRun2_asymptotic_v13-v2/MINIAODSIM",

]

flag = "RunIISummer20UL16MiniAODv2*"


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
                print '\t"'+line+'",', "# production"
            else:
                if not '"status":"INVALID"' in prodstatus:
                    print '\t"'+line+'",'
    else:
        print "not found:", '/%s/%s/MINIAOD*'%(name,flag)
    

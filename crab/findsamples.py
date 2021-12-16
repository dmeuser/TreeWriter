



samples = [
    "/MET",

]

flag = "Run2016*21Feb2020*UL2016*"


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
    

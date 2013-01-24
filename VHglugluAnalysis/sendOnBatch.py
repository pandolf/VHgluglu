#! /usr/bin/env python
import os
import sys
import time
import re
# set parameters to use cmst3 batch 
#######################################
### usage  cmst3_submit_manyfilesperjob.py dataset njobs applicationName queue 
#######################################
if (len(sys.argv) != 3) and (len(sys.argv) != 4) and (len(sys.argv) != 5):
    print "usage sendOnBatch.py dataset filesPerJob analyzerType=\"VHgluglu\" flags=\"\""
    sys.exit(1)
dataset = sys.argv[1]
inputlist = "files_"+dataset+".txt"
#settingfile = "config/RSZZsettings.txt"
# choose among cmt3 8nm 1nh 8nh 1nd 1nw 
#queue = "cmst3"
#queue = "cms8nht3"
queue = "8nh"
#queue = "2nd"
#ijobmax = 40
ijobmax = int(sys.argv[2])

analyzerType = "VHgluglu"
if len(sys.argv) >= 4:
    analyzerType = sys.argv[3]
flags = ""
if len(sys.argv) >= 5:
    flags = sys.argv[4]
application = "do2ndLevel_"+analyzerType

# to write on the cmst3 cluster disks
################################################
castordir = "/castor/cern.ch/user/p/pandolf/NTUPLES/" + dataset
pnfsdir = "/pnfs/roma1.infn.it/data/cms/store/user/pandolf/NTUPLES/" + dataset
afsdir = "/afs/cern.ch/user/p/pandolf/scratch0/NTUPLES/"+dataset
#outputmain = castordir+output
# to write on local disks
################################################
diskoutputdir = "/cmsrm/pc25_2/pandolf/MC/Summer12/"+dataset



match_Summer12 = re.search( r'Summer12', dataset, re.M|re.I)
if match_Summer12:
    diskoutputdir = "/cmsrm/pc25_2/pandolf/MC/Summer12/"+dataset
isData = re.search( r'Run201', dataset, re.M|re.I)
if isData:
    diskoutputdir = "/cmsrm/pc25_2/pandolf/data/"+dataset

diskoutputmain2 = afsdir
diskoutputmain = diskoutputdir
# prepare job to write on the cmst3 cluster disks
################################################
dir = analyzerType + "_" + dataset
os.system("mkdir -p "+dir)
os.system("mkdir -p "+dir+"/log/")
os.system("mkdir -p "+dir+"/input/")
os.system("mkdir -p "+dir+"/src/")
isEOS = re.search( r'eos', diskoutputdir, re.M|re.I)
if isEOS:
    os.system("eos mkdir -p "+diskoutputdir)
    diskoutputdir = "root://eoscms//" + diskoutputdir
else:
    os.system("ssh -o BatchMode=yes -o StrictHostKeyChecking=no pccmsrm25 mkdir -p "+diskoutputmain)
#outputroot = outputmain+"/root/"
#if castordir != "none": 
#    os.system("rfmkdir -p "+outputmain)
#    os.system("rfmkdir -p "+outputroot)
#    os.system("rfchmod 777 "+outputmain)
#    os.system("rfchmod 777 "+outputroot)
#else: os.system("mkdir -p "+outputroot)


#look for the current directory
#######################################
pwd = os.environ['PWD']
#######################################
inputListfile=open(inputlist)
inputfiles = inputListfile.readlines()
ijob=0

#copy the configuration in the actual run directory
#os.system("cp -r config "+dataset_name)

while (len(inputfiles) > 0):
    inputfilename = pwd+"/"+dir+"/input/input_"+str(ijob)+".list"
    inputfile = open(inputfilename,'w')
    for line in range(min(ijobmax,len(inputfiles))):
        ntpfile = inputfiles.pop()
        if ntpfile != '':
            inputfile.write(ntpfile)


    inputfile.close()

    # prepare the script to run
    outputname = dir+"/src/submit_"+str(ijob)+".src"
    outputfile = open(outputname,'w')
    outputfile.write('#!/bin/bash\n')
    outputfile.write('export LANGUAGE=C\n')
    outputfile.write('export LC_ALL=C\n')
    #outputfile.write('cd /afs/cern.ch/user/a/amarini/scratch0/CMSSW_4_2_5/src ; eval `scramv1 runtime -sh` ; cd -\n')
    outputfile.write('cd /afs/cern.ch/work/p/pandolf/CMSSW_5_2_5/src/ ; eval `scramv1 runtime -sh` ; cd -\n')
    outputfile.write('cd $WORKDIR\n')

    if flags=="":
      outputfile.write(pwd+'/'+application+" "+dataset+" "+inputfilename+" "+str(ijob)+"\n")
    else :
      outputfile.write(pwd+'/'+application+" "+dataset+" "+inputfilename+" "+flags+"_"+str(ijob)+"\n")
    if isEOS:
        outputfile.write('ls '+analyzerType+'*.root | xargs -i cmsStage {} '+diskoutputmain+'/{}\n') 
    else:
        outputfile.write('ls '+analyzerType+'*.root | xargs -i scp -o BatchMode=yes -o StrictHostKeyChecking=no {} pccmsrm25:'+diskoutputmain+'/{}\n') 
    outputfile.close
    os.system("echo bsub -q "+queue+" -o "+dir+"/log/"+dataset+"_"+str(ijob)+".log source "+pwd+"/"+outputname)
    os.system("bsub -q "+queue+" -o "+dir+"/log/"+dataset+"_"+str(ijob)+".log source "+pwd+"/"+outputname+" -copyInput="+dataset+"_"+str(ijob))
    ijob = ijob+1
    continue

#!/bin/env python
import subprocess,os,re
nanoaod_path = '/eos/cms/store/group/phys_muon/dmytro/tmp/NanoAOD-skims/mm/'
output_path = 'skim-lists'
version = '505'
for dir in subprocess.check_output('find %s/%s/ -mindepth 1 -maxdepth 1 -type d'%(nanoaod_path,version),shell=True).split("\n"):
    # print os.path.join(output_dir, os.path.basename(input_file).replace(".\ |root","_Skim.root"))
    if not re.search('\S',dir): continue
    name = os.path.basename(dir)
    print "Processing %s" % name
    command = "find %s/%s/%s -mindepth 1 -maxdepth 1 -type f| grep root > %s/%s.txt" % (nanoaod_path,version,name,output_path,name)
    subprocess.call(command,shell=True)



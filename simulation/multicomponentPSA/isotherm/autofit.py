#!/mechthild/software/apps/python-3.6/bin/python3.6

import os
import shutil
import sys
import subprocess
import re
import collections
import pyiast
import pickle
def findloading(srcfile):
    with open(srcfile,'r') as f:
        filelines = f.readlines()
    
        temploading = ""
        frameworkname= ""
        for linefff in filelines:

            if "Average loading absolute [mol/kg framework]" in linefff :
                #print(linefff)
                temploading= linefff[linefff.index("]")+1:linefff.index("+")].strip()
                loadingrelativerro= linefff[linefff.index("/-")+1:linefff.index("[-]")]

            if "Framework name:" in linefff:
                frameworkname= linefff[linefff.index(":")+1:].strip()

            if "External Pressure:" in linefff:
                pressure= linefff[linefff.index(":")+1:linefff.index("[")].strip()

  


    return frameworkname,pressure,temploading,loadingrelativerro

def saveaspickle(file, data):
    with open(file, 'wb') as f:
#    # Pickle the 'data' dictionary using the highest protocol available.
        pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)

def loadpickle(file):
    with open(file, 'rb') as f:
    # The protocol version used is detected automatically, so we do not
    # have to specify it.
        data = pickle.load(f)
    return data


if __name__=="__main__":
###################### post processing simulation data######
    presurelist =[1e3, 1e4, 5e4, 1e5, 5e5, 1e6]
    components =["C2H2","ethene","ethane"]
    componentmodel= ["Fischer", "TraPPE-UA", "TraPPE-UA"]
    mainpath = "/mechthild/home/zhouy/02_IsothermFitting"
    sourcemain="/mechthild/home/zhouy/02_IsothermFitting"
    sourceDatabase="/mechthild/home/zhouy/02_IsothermFitting/DatabaseEqeq"
    sourceTemplates="/mechthild/home/zhouy/02_IsothermFitting/Templates"

    recording_file=os.path.join(mainpath,"isotherm_data")

    dictdata = collections.defaultdict(dict)
    moflist=[]
    with open(recording_file,"w") as fin:
        fin.write(" ")
    for i in components:
        with open(recording_file,"a") as fin:
            fin.write("{}----------\n".format(i))
        for j in presurelist:       
            basedir = os.path.join(mainpath,i+"_"+str(j))
            ScreenDir=os.path.join(basedir,"Screening")

            for sub in os.listdir(ScreenDir):
                finalresultfile_parentfolder= os.path.join(ScreenDir,sub,"cmbc","Output","System_0")
                files=os.listdir(finalresultfile_parentfolder)
                for file in files:
                    filepath = os.path.join(finalresultfile_parentfolder,file)
                    tempmof,temppressure,temploading,temploadingerro = findloading(filepath)
                    dictdata[tempmof][temppressure] = temploading
                    moflist.append(tempmof)
                    #strtoouputfile="{},{},{}".format(tempmof,temppressure,temploading)
        moflist=list(set(moflist))
        list.sort(moflist)
        print(moflist)
        with open(recording_file,"a") as fin:
            for mofitems in moflist:
                fin.write("{} [Pa],{}\n".format(mofitems,dictdata[mofitems])) 
        saveaspickle(os.path.join(mainpath,'isotherms.pickle'),dictdata)



    # mofname,pressure,loading = findloading("/mechthild/home/zhouy/02_IsothermFitting/C2H2_1000.0/Screening/CUVGOQ/cmbc/Output/System_0/output_CUVGOQ_5.3.2_298.000000_1000.data")
    # print(mofname,pressure,loading)
################## isotherm fitting ######################








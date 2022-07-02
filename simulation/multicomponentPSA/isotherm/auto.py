#!/mechthild/software/apps/python-3.6/bin/python3.6

import os
import shutil
import sys
import subprocess



def distribute(src, dst):
    if not os.path.exists(os.path.join(dst,"DatabaseEqeq")):
        shutil.copytree(os.path.join(src,"DatabaseEqeq"), os.path.join(dst,"DatabaseEqeq"))
    if not os.path.exists(os.path.join(dst,"Templates")):
        shutil.copytree(os.path.join(src,"Templates"), os.path.join(dst,"Templates"))
    
    return True


def modifysimulationinput(src, pressure, component,componentmodel):
    # Read in the file
    with open(src, 'r') as file :
        filedata = file.read()

# Replace the target string
    filedata = filedata.replace('targetpressure', str(pressure))
    filedata = filedata.replace('molname', component)
    filedata = filedata.replace('molmodel', componentmodel)

# Write the file out again
    with open(src, 'w') as file:
        file.write(filedata)
    return True

def modifymakefile(src, basedir):
    with open(src, 'r') as file :
        filedata = file.read()

# Replace the target string
    filedata = filedata.replace('NoneBaseDirectory', basedir)

# Write the file out again
    with open(src, 'w') as file:
        file.write(filedata)

def modifybatchsubmitionfile(src,DIR):

    numberofsubfile = len([name for name in os.listdir(DIR) if os.path.isfile(os.path.join(DIR, name))])

    with open(src, 'r') as file :
        filedata = file.read()

# Replace the target string
    filedata = filedata.replace('numberofsubmitionfile', str(numberofsubfile))

# Write the file out again
    with open(src, 'w') as file:
        file.write(filedata)

if __name__=="__main__":
    #set desired pressure, component, and choose molecular models accordingly
    #change the main path varibale to your new project directory
    presurelist =[1e3, 1e4, 5e4, 1e5, 5e5, 1e6]
    component =["C2H2","ethene","ethane"]
    componentmodel= ["Fischer", "TraPPE-UA", "TraPPE-UA"]
    mainpath = "/mechthild/home/zhouy/02_IsothermFitting"
    sourcemain="/mechthild/home/zhouy/02_IsothermFitting"
    sourceDatabase="/mechthild/home/zhouy/02_IsothermFitting/DatabaseEqeq"
    sourceTemplates="/mechthild/home/zhouy/02_IsothermFitting/Templates"
    
    for i,k in zip(component,componentmodel):
        for j in presurelist:       
            basedir = os.path.join(mainpath,i+"_"+str(j))
            
            if not os.path.exists(basedir):  
                os.makedirs(basedir)

            distribute(sourcemain,basedir)

            modifysimulationinput(os.path.join(basedir,"Templates","cmbc","simulation.input"),j,i,k)

            makefilepath=os.path.join(basedir,"Templates","make","makefiles")
            modifymakefile(makefilepath,basedir)

            DatabaseDIR=os.path.join(basedir,"DatabaseEqeq")
            SubmitonFile=os.path.join(basedir,"Templates","make","batch_new")
            modifybatchsubmitionfile(SubmitonFile,DatabaseDIR)


            command=["cd {}".format(os.path.join(basedir,"Templates","make")) + ";"
            + "ls" +";"
            + "chmod +x makefiles" +";"
            + "./makefiles" +";"
            + "cp {} {}".format(os.path.join(basedir,"Templates","make","batch_new"),os.path.join(basedir,"Submit")) +";"
            + "cd {}".format(os.path.join(basedir,"Submit")) +";"
            + "sbatch batch_new"
            
            
            
            
            
            ]
            subprocess.run(command,shell=True)
            #print("cd {}".format(os.path.join(basedir,"Templates","make")))
            #subprocess.run(["ls"],shell=True)
            #os.system("chmod +x makefiles")
            #os.system("./makefiles")
            #os.system("cp {} {}".format(os.path.join(basedir,"Templates","make","batch_new"),os.path.join(basedir,"Submit")))
            #os.system("cd {}".format(os.path.join(basedir,"Submit")))
            #os.system("sbatch batch_new")

    
    templatepath = mainpath + "Templates/"

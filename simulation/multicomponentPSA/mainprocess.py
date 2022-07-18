#!/usr/bin/env python3.7

import util
import onebedmulticomponent
import params
import collections
import os
import numpy as np
import pygaps
from isothermfit import getisothermparam


#Step 1: Specifiy simulation parameters,some parameters are overode
#for simulation make sure you see "params overiding ... " in terminal when running the script,only for qsi,qs0 and bi
isothermdatabase = util.loadpickle(os.path.join(os.getcwd(),'graph','LangmuirIsotherm.pickle'))#load isotherm parameters from pickle(done by isotherfit.py)

moflist=['NEXXEV'] # one specific mof 
#mofid=isothermparameters['C2H2']#simulate for all mofs in the database
adsorbates=["C2H2","ethene","ethane"] #adsorbate molecules

paramoverload =collections.defaultdict() #create a dictionary where parameters that need to be overloaded to create_param function
#other process varibales
paramoverload['feed_pressure'] = 10
paramoverload['PI'] = 1e5
paramoverload['PL'] = 1e4

for mofid in moflist:
    mofdict=getisothermparam(isothermdatabase, mofid, adsorbates)
    paramoverload['qsi']=mofdict['qsi']
    paramoverload['qs0']=mofdict['qs0']
    paramoverload['bi']=mofdict['bi']
    param= params.create_param(paramoverload)

    #Step 2: Simulate the 4-step psa process until cyclic steady state is reached
    #finalcyle, status = onebedmulticomponent.cyclic_steady_state(param,tolerance = 3e-3) # simuation until css is reached

    #data storage as a pickled file
    outputdir = os.path.join(os.getcwd(),'outcome')
    util.safe_mkdir(outputdir)
    filepickle = os.path.join(outputdir, 'statestatus{}.pickle'.format(mofid))
    #util.saveaspickle(filepickle, status) # save for the psa simulation data

    #Step3: postprocessing and plotting
    # compute purity and recovery
    status = util.loadpickle(filepickle)
    param =status.param
    #onebedmulticomponent.plot_data(status.snap[status.cycle],param,os.path.join(outputdir,mof,'{}'.format(mof)))
    param.collect= ['zL_tblw','zL_tads','z0_teva']
    #status.obj = onebedmulticomponent.purity_recovery(status.snap[status.cycle],param)
    if status.obj:
        print(status.obj['purity_ads'], status.obj['recovery'], mofid )
    #util.saveaspickle(filepickle, status) # save for the calculated purity and recovery

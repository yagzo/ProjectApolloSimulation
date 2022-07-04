#!/usr/bin/env python3.7

import util
import onebedmulticomponent
import params
import collections
import os
import numpy as np
#Step 1: Specifiy simulation parameters,some parameters are overode
#for simulation make sure you see "params overiding ... " in terminal when running the script,only for qsi,qs0 and bi
paramoverload =collections.defaultdict() #create a dictionary where parameters that need to be overloaded to create_param function
isothermparameters = util.loadpickle(os.path.join(os.getcwd(),'graph','LangmuirIsotherm.pickle'))#load isotherm parameters from pickle(done by isotherfit.py)

mofid=['NEXXEV'] #
adsorbates=["C2H2","ethene","ethane"]
paramoverload['qsi']=np.zeros((len(adsorbates),)) # intialize with an empty dictionary
paramoverload['qs0']=np.zeros((len(adsorbates),))
paramoverload['bi']=np.zeros((len(adsorbates),))
#optimization varibales
paramoverload['feed_pressure'] = 10
paramoverload['PI'] = 1e5
paramoverload['PL'] = 1e4
for mof in mofid:
    for i in adsorbates:
        paramoverload['qsi'][adsorbates.index(i)]= isothermparameters[i][mof].model.params['n_m']
        paramoverload['bi'][adsorbates.index(i)]=isothermparameters[i][mof].model.params['K']
    paramoverload['qs0']=max(paramoverload['qsi']) #reference loading is set to be the highest saturated loading among all components
    param= params.create_param(paramoverload)

    #Step 2: Simulate the 4-step psa process until cyclic steady state is reached
    #finalcyle, status = onebedmulticomponent.cyclic_steady_state(param,tolerance = 3e-3) # simuation until css is reached

    #data storage as a pickled file
    outputdir = os.path.join(os.getcwd(),'outcome' )
    util.safe_mkdir(outputdir)
    filepickle = os.path.join(outputdir, 'statestatus{}.pickle'.format(mof))
    #util.saveaspickle(filepickle, status)

    #Step3: postprocessing and plotting
    # compute purity and recovery
    status = util.loadpickle(filepickle)
    param =status.param
    onebedmulticomponent.plot_data(status.snap[status.cycle],param,os.path.join(outputdir,mof,'{}'.format(mof)))
    param.collect= ['zL_tblw','zL_tads','z0_teva']
    status.obj = onebedmulticomponent.purity_recovery(status.snap[status.cycle],param)
    print(status.obj)
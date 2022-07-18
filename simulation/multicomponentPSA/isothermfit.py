#!/usr/bin/env python3.7

import collections
import numpy as np
import pygaps
import pickle
import os
import util
import pandas as pd

def unpackisotherm(isoinput,mofid,adsorbates):
    # input:isoinput is a dictionary containing isoinput, ['adsorbate']['mofid']=<class: pygaps.modelisotherm>
    #       mofid = ['mof1','mof2',...]
    #       adsorbate =['C2H2','ethene','ethane']
    # output: iso_ndarray: type: ndarray
    iso_ndarray=collections.defaultdict()
    for mof in mofid:
        iso_ndarray[mof] = collections.defaultdict()
        iso_ndarray[mof]['qsi'] = np.zeros((len(adsorbates),))
        iso_ndarray[mof]['bi'] = np.zeros((len(adsorbates),))
        for i in adsorbates:
            iso_ndarray[mof]['qsi'][adsorbates.index(i)]= isoinput[i][mof].model.params['n_m']
            iso_ndarray[mof]['bi'][adsorbates.index(i)]=isoinput[i][mof].model.params['K']


    return pd.DataFrame.from_dict(iso_ndarray) # isotherm dataframe

def getisothermparam(isothermdb, mofname, adsorbates):
    # This function can be used to get isotherm parameters of a mof named (mofname) from the isothermdatabase (isothermdb)
    # INPUT:  isothermdb: a database for isotherms(fitted) usually stored as a pickle file
    #         mofname: a string defined for the identification of a specific mof e.g."NEXXEV"
    #         adsorbates: adsorbate list e.g. ["C2H2","ethene","ethane"]
    # RETURN: mofdict.qsi=[qsA,qsB,qsC,...]
    #       : mofdict.qs0= max(qsA,qsB,qsC,...)
    #       : mofdict.bi =[bA,bB,bC,...]
    mofdict =collections.defaultdict() #create a dictionary where parameters that need to be overloaded to create_param function
    mofdict['qsi']=np.zeros((len(adsorbates),)) # intialize with an empty dictionary
    mofdict['qs0']=np.zeros(1)
    mofdict['bi']=np.zeros((len(adsorbates),))
    mofdict['name'] = mofname #record the mofname in the dictionary
    for i in adsorbates:
        mofdict['qsi'][adsorbates.index(i)]= isothermdb[i][mofname].model.params['n_m']
        mofdict['bi'][adsorbates.index(i)]=isothermdb[i][mofname].model.params['K']
    mofdict['qs0']=max(mofdict['qsi'])      
    return mofdict

def makedb():
    isothermdict = util.loadpickle('./isotherm/isotherms.pickle')
    #targetmof=['ZADDAJ','FENVOL']
    targetmof=isothermdict["C2H2"]
    #print(targetmof)
    #targetmof=['ORAQUU','FENVOL','ZUQVIQ','CUNXIS','CUNXIS10','GIHBII','NEXXEV','JAVTAC'] #the MOF which will be computed, it is suppose to be a subset of isothermdict[mof]
    targetadsorbate=["C2H2","ethene","ethane"]
    isotherm1=collections.defaultdict()#database for isotherm storage
    isodict = collections.defaultdict()#database for isotherm parameter storage(namely, K and n_m, fit to Langmuir isotherm)
    #print(isothermdict)
    for i in targetadsorbate:
        isotherm1[i] = collections.defaultdict() #for every component create an empty defualtdict
        isodict[i] = collections.defaultdict()
        for mof in targetmof: #unpack mofdict so that it conforms to the input format of the pygaps package
            pressurelist=[] #for each mof iteration empty the previous pressurelist
            loadinglist=[] #same as above except that it is for record ads loading
            for p in isothermdict[i][mof]: #unpack pressure
                pressurelist.append(float(p))
                loadinglist.append(float(isothermdict[i][mof][p]))
            
            materials=mof
            #unpack for one mof finished
            # create an isotherm
            print(pressurelist, loadinglist, materials, i)
            isotherm1[i][mof] = pygaps.PointIsotherm(
                pressure=pressurelist,
                loading=loadinglist,
                material=materials,
                adsorbate= i,
                temperature= 298,
            
                pressure_unit= 'Pa',
                pressure_mode='absolute')
            isotherm1[i][mof].convert_pressure(unit_to='bar')
            #isotherm1[i][mof].plot(save_path='.')
            #model = pygaps.ModelIsotherm.from_pointisotherm(isotherm1, model='Langmuir', verbose=True)
            try:
                isodict[i][mof] =pygaps.modelling.model_iso(isotherm1[i][mof],model='Langmuir',verbose=False)
            except pygaps.utilities.exceptions.CalculationError:
                anotherguess={'n_m': 0.000002, 'K': 1e8}
                isodict[i][mof] =pygaps.modelling.model_iso(isotherm1[i][mof],model='Langmuir',verbose=False,param_guess=anotherguess)

            savepath = os.path.join(os.getcwd(),'graph')
            if not os.path.isdir(savepath):
                os.mkdir(savepath)
            isodict[i][mof].print_info(save_path= os.path.join(savepath,'{}_on_{}'.format(i,mof)))
            

    unpackisotherm(isodict,targetmof,targetadsorbate).to_excel(os.path.join(savepath,'LangmuirParameters.xlsx'))   # store as pd.dataframe to pickle
    util.saveaspickle(os.path.join(savepath,'LangmuirIsotherm.pickle'), isodict)
    return None


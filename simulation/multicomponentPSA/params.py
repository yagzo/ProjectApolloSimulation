import math
import string
import numpy as np
#Copyright (c) 2018 Daniel B. Grunberg
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

'''
Parameters for the finite-volume simulation
of an oxygen concentrator

'''
#constants of the system
class AttrDict(dict):
  def __init__(self, *args, **kwargs):
    super(AttrDict, self).__init__(*args, **kwargs)
    self.__dict__ = self

def compnames(numberofcomponents=2):
    # declare component related state variable names ie.xi,yi given a number of components,
    # input: a number of components (type: an interger between 1 and 26)
    # output: a concatenated list of state variables'names (type: a list of strings)
    # eg. input: numebrofcomponents=2 (default value)
    #     output: ['xA','xB','yA','yB']
    # eg. input: numebrofcomponents=3
    #     output: ['xA','xB','xC','yA','yB','yC']
    assert numberofcomponents<=26 and numberofcomponents >=1
    ylist=[]
    xlist=[]
    alphabet_list = list(string.ascii_uppercase)
    for i in range(numberofcomponents):
        ylist.append("y{}".format(alphabet_list[i]))
        xlist.append("x{}".format(alphabet_list[i]))
    return xlist,ylist

def product_loc(bi,ynames):
  # set the product location parameter, which is used to calculate recovery, from the isotherm parameter Bi
  # Hueristics 1: the larger the Bi, the stronger the interactions
  #            2: Heavy component is collected at the evacuation step ('z0_teva')
  #            3: Light product is collected at the adsorption step ('zL_tads')
  #            4: other are collected at the desorption step ('zL_tblw')
  #
  # test example:
  #               ynames = ['yA', 'yB','yC']
  #               bi = [4, 3, 2]
  #               from params import product_loc
  #               aa= product_loc(bi,ynames) 
  #               print(aa)
  # output: aa= ['yC', 'yB', 'yA']
  #   
  # please note, this is based on purely heuristics, which does necessarily garuntee a correct setting
  # please note, currently, only 3 column ends are coded to be able to collect product streams, they are 'zL_tblw', 'zL_tads','z0_teva'
  orderedynames= [None] * len(bi) #initialization a list with length equal to ynames and bi 
  indices = np.argsort(bi) # sorted in an ascending order
  for k in indices:
    orderedynames[indices.tolist().index(k)] = ynames[k]
  # the lightest component (raffinate)
  # the heaviest component
  return orderedynames

def create_param(mods, print=print,verbose=False):
  # if verbose== True, overloaded parameters are print to screen, default set to False
  #mods is a dict to override the parameter defaults set here
  #override print function if needed for debug statements to a log file
  #any keys with value of None are ignored
  param=AttrDict()





  #Component related parameters
  param.nocomponents = 3 # number of components
  param.feed_yi=[0.05, 0.9, 0.05] #feed gas fraction of component i
  param.ini_yi = [0.1, 0.8, 0.1] # gas composition used for initialization state variables yi before calling the odes solver
  param.isomodel='Langmuir-Freundlich' # types of isotherm model (curruntly available:'Langmuir-Freundlich'; )
  param.bi=[0.0086, 0.25, 0.0086]   # [1/bar] @ Room temp
  param.qsi=[3.09,5.84, 3.09]  #saturation constant for component i  [mol/kg]=[mmol/g] (at infinite pressure)
  param.qs0 = 5.84 #[mol/kg]
  param.ni=[1.00, 1.00, 1.00]     #exponent on bar pressure of component i (note: not 1/ni)
  param.ki=[0.620, 0.620, 0.620]    #mass transfer coefficient
  #not doing Temp and Tw simulation at the moment
  param.Hi=[0, 0, 0]                    #heat of adsorption of component i  [J/mol]

  #Simulation parameters
  #the mode for the difference scheme
  param.mode=1 #differencing scheme 1: forward difference； 2：vanleer； 4：weno 
  param.N= 10  #discretization steps in spatial direction
  param.tN= 51 # number of sampling points over the entire cycle time span (*both ends of the line segement are included; starts from 0; must be an integer greater or equal than 2)
  
  #Physical system parameters
  param.Ta=298    #room temperature in K
  #Bed Size
  param.D=0.30     #Diameter of cylinder [m]
  param.L=1      #Length of cylinder (m)  - 13"
  param.epsilon=0.37  # void fraction (fraction of space filled with gas vs. spheres)
  param.rp=0.7e-3          #particle radius [m]
  
  #Cycle operating parameters
  param.feed_pressure = 1        #pressure in bar (1e5 Pa), absolute #this is also set as the highest pressure
  param.vfeed = 1 # intersitial feed velocity [m/s]
  param.PI= 0.5*1e5   # intermediate pressure in Pa                
  param.PL=0.25*1e5 # low pressure in Pa
 
  #NOTE:Hyper parameters for determine which setup of PSA is adopted
  param.bednum = 1 #number of beds default 1 bed 
  param.bed = "VPSA" #default key denoting operating cycles (4-step VPSA or 6-step DR-PSA)

  if mods['bed'] == "VPSA" or not('bed' in mods.keys()): # 4-step VPSA(default time length)
    param.bed = mods['bed']
    param.lamda_pre=0.5 # [s-1] pressurization profile parameter
    param.lamda_blw=0.5 # [s-1] blowdown pressure profile parameter
    param.lamda_eva=0.5 # [s-1] evacuation pressure profile parameter
    param.tpre = 15 # [s] pressurization duration time
    param.tads = 15 # [s] adsorption duration time
    param.tblw = 30 # [s] blowdown duration time
    param.teva = 40 # [s] evacuation duration time
  
  if mods['bed'] == "DRPSA": # 6-step DR-VPSA(default time length)
    param.bed = mods['bed']
    param.lamda_pre=0.5 # [s-1] pressurization profile parameter
    param.lamda_blw=0.5 # [s-1] blowdown pressure profile parameter
    param.lamda_eva=0.5 # [s-1] evacuation pressure profile parameter
    param.tpre = 15 # [s] pressurization duration time
    param.tads = 15 # [s] adsorption duration time
    param.thr = 30 # [s] heavy reflux duration time
    param.tblw = 30 # [s] blowdown duration time
    param.teva = 40 # [s] evacuation duration time
    param.tlr = 30 # [s] light reflux duration time
    param.R_LR = 0.5 # light reflux ratio = vin_LR/vfeed 
    param.R_HR = 0.5 # heavy reflux ratio = vin_HR/vfeed 

  #Adsorbent and fluid dynmaics parameters 
  param.rho_s=1130       #adsorbent density kg/m3 (was rho_pellet)
  param.rho_w=7800       #wall density kg/m3
  param.Cpg=30.7     #specific heat capacity of gas phase [J/mol/K]
  param.Cpa=30.7     #specific heat capacity of adsorbed phase [J/mol/K]
  #param.Cps=1070     #specific heat capacity of adsorbent [J/kg/K]
  #from Jee paper
  param.Cps=1339     #specific heat capacity of adsorbent [J/kg/K]
  param.Cpw=502      #specific heat capacity of wall [J/kg/K]
  param.mu=1.72e-5   #fluid viscocity
  param.Dm=1.6e-5    #Molecular diffusivity [FOR CO2/N2 paper]
  param.gamma=1.4    #adiabatic constant
  param.Kz=0.09      #Effective gas thermal conductivity [J/m/K/s]
  param.Kw=16        #Effective thermal cond of wall
  param.hin=8.6      #inside heat transfer coeff
  param.hout=2.5     #outside heat transfer coeff
  param.R=8.314      #Gas constant  [m3Pa/mol/K]
  param.eta=.72      #compression/evacuation efficiency
  param.DL=4.9e-4  # m2/sec from another paper


  #####
  #Apply overrides, ignoring None values.  But see cycles/time logic below
  for k in param.keys():
    if mods != None:
      if k in mods and mods[k] is not None:
        param[k]=mods[k]
        if verbose:# if verbose==True(default is false) then print the overloaded parameters to screen
          print('params overriding {}:{}'.format(k, param[k]))
    
  #####
  #Things we have to calculate, need to do this after the overriding
  #check the input parameters to see if they are consitant, e.g. for a ternary mix
  assert len(param.feed_yi) == param.nocomponents
  assert len(param.ini_yi) == param.nocomponents
  assert len(param.bi) == param.nocomponents
  assert len(param.qsi) == param.nocomponents
  assert len(param.ni) == param.nocomponents
  assert len(param.ki) == param.nocomponents
  assert len(param.Hi) == param.nocomponents

  #setup component names for xi, yi
  param.xnames, param.ynames = compnames(param.nocomponents)
  param.state_names=['P','T','Tw'] + param.xnames + param.ynames
  param.state_sizes=[] # will be filled in when we do init()
  for s in param.state_names:
    param.state_sizes.append(param.N)
  param.area=(param.D/2)**2*math.pi     #[m2]
  param.volume=param.area*param.L

  param.rin=param.D/2
  param.rout=param.rin + 0.25
  param.epst=param.epsilon/(1-param.epsilon)
  param.epsb=(1-param.epsilon)/param.epsilon
  #param.cellsize=param.L/param.N
  #param.deltaZ=param.cellsize/param.L   # nondimensional deltaZ
  param.deltaZ=1/param.N  # nondimensional deltaZ
  param.container_vol=param.area*param.L #volume of reactor (m3)


  #reference value of the statevaribales for normalization
  param.PH=param.feed_pressure*1e5   #Pa, will be the normalization pressure
  param.norm_v0=param.vfeed
  param.norm_t0=param.L/param.norm_v0   # so t=time/norm.t0
  param.norm_P0=param.PH
  param.norm_T0=298 # [K] temperature
  assert sum(param.feed_yi) == 1 and len(param.feed_yi) == param.nocomponents

  #parameters used in setting boundary conditions
  if param.bed == "VPSA": #4-step VPSA 
    param.norm_tpre = param.tpre / param.norm_t0
    param.norm_tads = (param.tpre +param.tads) / param.norm_t0
    param.norm_tblw = (param.tpre +param.tads + param.tblw) / param.norm_t0
    param.norm_teva = (param.tpre+ param.tads + param.tblw + param.teva) / param.norm_t0
    param.norm_tend = param.norm_teva # alias for teva used in plots.py
    param.tstep = param.norm_teva/(param.tN-1)      #time step for ODE (dimensionless), -1 because the number of intervals = number of points(of a line segement including both ends) - 1,NOTE: Careful, avoid using arange,because when time step is not an integer,the rightmost time point is not included!

    # the index of descritized time varaible at the end of each step  
    param.index_tpre = math.floor(param.norm_tpre/param.tstep)
    param.index_tads = math.floor(param.norm_tads/param.tstep)
    param.index_tblw = math.floor(param.norm_tblw/param.tstep)
    param.index_teva = math.floor(param.norm_teva/param.tstep)

    #define product stream
    #param.collect = ['zL_tads','z0_teva','zL_tblw'] # where each component is collected, following the format "z+[0/L](z=0/z=L)+_t[pre/ads/blw/eva](three letter code indicating each step)"
    param.collect = product_loc(param.bi,param.ynames)

  if param.bed == "DRPSA": #6-step DR-PSA
      #normalized time durations for a DR-PSA 
    param.norm_tpre =  param.tpre / param.norm_t0 
    param.norm_tads = (param.tpre + param.tads) / param.norm_t0
    param.norm_thr  = (param.tpre + param.tads + param.thr) / param.norm_t0 
    param.norm_tblw = (param.tpre + param.tads + param.thr + param.tblw) / param.norm_t0
    param.norm_teva = (param.tpre + param.tads + param.thr + param.tblw + param.teva) / param.norm_t0
    param.norm_tlr  = (param.tpre + param.tads + param.thr + param.tblw + param.teva + param.tlr ) / param.norm_t0
    param.norm_tend = param.norm_tlr # which step is the last of a cycle
    param.tstep = param.norm_tlr/(param.tN-1)      #time step for ODE (dimensionless), -1 because the number of intervals = number of points(of a line segement including both ends) - 1,NOTE: Careful, avoid using arange,because when time step is not an integer,the rightmost time point is not included!

    # the index of descritized time varaible at the end of each step  
    param.index_tpre = math.floor(param.norm_tpre/param.tstep)
    param.index_tads = math.floor(param.norm_tads/param.tstep)
    param.index_thr =  math.floor(param.norm_thr/param.tstep)
    param.index_tblw = math.floor(param.norm_tblw/param.tstep)
    param.index_teva = math.floor(param.norm_teva/param.tstep)
    param.index_tlr =  math.floor(param.norm_tlr/param.tstep)
    param.lp_yi= param.feed_yi
    param.hp_yi= param.feed_yi
    param.collect = product_loc(param.bi, param.ynames)
  if verbose:
    print("parameter initialization completed")
  return param

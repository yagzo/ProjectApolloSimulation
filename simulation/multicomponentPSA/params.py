import math
import string
import orifice
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


def create_param(mods, print=print,verbose=False):
  # if verbose== True, overloaded parameters are print to screen, default set to False
  #mods is a dict to override the parameter defaults set here
  #override print function if needed for debug statements to a log file
  #any keys with value of None are ignored
  #real_cycle_time and real_vent_time will override cycle_time and vent_time
  #cycles and time are used to set the simulate end time
  param=AttrDict()
  #Physical constants
  param.R=8.314     #J/mol-K
  #Component related parameters
  param.nocomponents = 3 # number of components
  param.feed_yi=[0.05, 0.9, 0.05] #feed gas fraction of component i
  param.ini_yi = [0.1, 0.8, 0.1] # gas composition used for initialization state variables yi before calling the odes solver
  param.collect = ['zL_tads','z0_teva','zL_tblw'] # where each component is collected, following the format "z+[0/L](z=0/z=L)+_t[pre/ads/blw/eva](three letter code indicating each step)"
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
  param.bed=1 #number of operating beds, its value must be either 1 or 2
  param.N= 10  #discretization steps in spatial direction
  param.tN= 51 # number of sampling points over the entire cycle time span (*both ends of the line segement are included; starts from 0; must be an integer greater or equal than 2)
  #param.adsorb=True    #set to True to have adsorption happen
  #cycle_time, vent_time
  #param.cycle_time=36       #time for a half-cycle in dimensionless time units
  #param.vent_time=23.04        #from beginning of cycle in dimensionless time units
  #Physical system parameters
  param.Ta=298    #room temperature in K
  #Approximate Inviracare XL cylinders
  param.D=0.30     #Diameter of cylinder [m]
  param.L=1      #Length of cylinder (m)  - 13"
  #Soda Bottle 2L
  #param.D=0.10     #Diameter of cylinder [m]
  #param.L=0.33      #Length of cylinder (m)  - 13"
  #param.product_tank_volume=2.0/1000     # 2L in m3
  #we set an orifice so that we can 5LPM at a certain pressure (3bar abs) and
  #try to get an output pressure sits around that 3 bar.
  #param.epsilon=0.36  # void fraction (fraction of space filled with gas vs. spheres)
  param.epsilon=0.37
  #param.epsilonp=0.35    #particle voidage
  param.rp=0.7e-3          #particle radius [m]
  #param.tauprime=3.0     # tortuosity
  #flow params
  #param.input_orifice=2.2     # dia mm - pressure feed
  #param.output_orifice=1.40    # dia mm - to output tank
  #we use an orifice that will provide 5LPM at 3 bar abs product tank
  #param.product_orifice=0.475   # dia mm - output from product tank to patient
  #param.blowdown_orifice=2.80     # dia mm
  #param.vent_orifice=1.0        # dia mm
  param.feed_pressure = 1        #pressure in bar (1e5 Pa), absolute #this is also set as the highest pressure
  param.vfeed = 1 # intersitial feed velocity [m/s]

  param.PI= 0.5*1e5   # intermediate pressure in Pa                
  param.PL=0.25*1e5 # low pressure in Pa
 
  
  #parameters used in setting boundary conditions for 4 operating steps 
  param.lamda_pre=0.5 # [s-1] pressurization profile parameter
  param.lamda_blw=0.5 # [s-1] blowdown pressure profile parameter
  param.lamda_eva=0.5 # [s-1] evacuation pressure profile parameter
  param.tpre = 15 # [s] pressurization duration time
  param.tads = 15 # [s] adsorption duration time
  param.tblw = 30 # [s] blowdown duration time
  param.teva = 40 # [s] evacuation duration time

  
  #param.product_pressure=2     # bar in output product tank (not used)
  #Note: random sphere packing is 64%, densest is 74%, void fractions=.36, .26
  #param.DL=1.0        #Axial dispersion coefficient (m2/s)
  #Mol wt O=16, N=14
  #air is about 14e-3 kg/mol
  #These are for my calculation, note different from paper (paper was doing
  #CO2 and N2 separation, and used a more complex Langmuir dual-site model
  #param.qAs=5.26e-3  #saturation constant for O2  mol/cm3 (at infinite pressure)
  #param.qBs=5.26e-3  #saturation constant for N2  mol/cm3 (at infinite pressure)
  #Langmuir constants.  inverse of pressure mol/m3 where we reach 50% of saturation
  #Estimated these from chart in Yang book at 1 ATM.
  
  #param.bA=573.5   # cm3/mol, inverse of 42.5 bar
  #param.bB=2030    # cm3/mol, inverse of 12 bar
  #From Mofahari 2013 paper, for Zeolite 13X
  #param.bA=0.042    # 1/bar @ Room temp
  #param.bB=0.616    # 1/bar @ Room temp
  #From Jee paper, same as Mofahari but *0.10 for k3!!
  #adding correct exponent on Pressure 

  #From Sircar ch22 paper @ 25C
  #param.bA=0.0295    # 1/bar @ Room temp
  #param.bB=0.107    # 1/bar @ Room temp
  #param.qAs=0.00312  #saturation constant for O2  mol/g (at infinite pressure)
  #param.qBs=0.00312  #saturation constant for N2  mol/g (at infinite pressure)

  #For 5A:
  #param.bA=0.052    # 1/bar @ Room temp
  #param.bB=0.165    # 1/bar @ Room temp
  #param.qAs=0.0019  #saturation constant for O2  mol/g (at infinite pressure)
  #param.qBs=0.0025  #saturation constant for N2  mol/g (at infinite pressure)
  #param.kA=.15
  #param.kB=.05

  #From the new paper
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
  #print('canister volume is {} L'.format(param.volume*1000))
  param.rin=param.D/2
  param.rout=param.rin + 0.25
  param.epst=param.epsilon/(1-param.epsilon)
  param.epsb=(1-param.epsilon)/param.epsilon
  param.cellsize=param.L/param.N
  param.deltaZ=param.cellsize/param.L   # nondimensional deltaZ
  param.container_vol=param.area*param.L #volume of reactor (m3)
  #We compute the initial flow rate and use that for velocity normalization
  #in_flow=orifice.orifice_flow2(param.feed_pressure*1e5/1000,100,param.input_orifice,T=param.Ta)/60
  #in_vel=in_flow/param.area    # [m/s]
  #param.norm_v0=in_vel  # [m/s]  feed velocity, which varies

  #reference value of the statevaribales for normalization
  param.PH=param.feed_pressure*1e5   #Pa, will be the normalization pressure
  param.norm_v0=param.vfeed
  param.norm_t0=param.L/param.norm_v0   # so t=time/norm.t0
  param.norm_P0=param.PH
  param.norm_T0=298 # [K] temperature
  assert sum(param.feed_yi) == 1 and len(param.feed_yi) == param.nocomponents
  #molecular weight of air is 29 g/mol or 29e-03 kg/mol
  #param.rho_g=100000/param.R/param.Ta*29e-3   # Needs to be adjusted for pressure/Temp
  
  #param.Dp=param.Dm/param.tauprime


  param.norm_tpre = param.tpre / param.norm_t0
  param.norm_tads = (param.tpre +param.tads) / param.norm_t0
  param.norm_tblw = (param.tpre +param.tads + param.tblw) / param.norm_t0
  param.norm_tend = (param.tpre+ param.tads + param.tblw + param.teva) / param.norm_t0
  param.tstep = param.norm_tend/(param.tN-1)      #time step for ODE (dimensionless), -1 because the number of intervals = number of points(of a line segement including both ends) - 1,NOTE: Careful, avoid using arange,because when time step is not an integer,the rightmost time point is not included!

  # the index of descritized time varaible at the end of each step  
  param.index_tpre = math.floor(param.norm_tpre/param.tstep)
  param.index_tads = math.floor(param.norm_tads/param.tstep)
  param.index_tblw = math.floor(param.norm_tblw/param.tstep)
  param.index_tend = math.floor(param.norm_tend/param.tstep)

  #NOTE: norm.v0 is computed after this, will be set to param.v0
  #Axial dispersion coefficient m2/sec
  #param.DL=0.7*param.Dm+0.5*param.norm_v0*param.rp*2.0

  #print('param.DL={}'.format(param.DL))
  #compute param.end_time for the simulation
  #want exactly 1 to be not None and 1 None
  #assert (mod.cycles is None) + (mod.time is None) == 1
  #total time to simulate (dimensionless).  This is handled differently because
  #we allow 2 ways to override (cycles takes precedence over time)
  #First we check for real_cycle_time and real_vent_time
  if verbose:
    print("parameter initialization completed")
  return param

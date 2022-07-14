#!/usr/bin/env python
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
from cgi import parse_multipart
from cmath import inf
import sys
import math
import collections
import pymoo


import numpy as np
from numpy import zeros_like
from numpy import column_stack
from numpy import hstack
from pymoo.core.problem import ElementwiseProblem
import numpy
import scipy
import scipy.integrate
import matplotlib
import time
import os
import warnings
import importlib
import pickle
##local modules
import plots
import params

#Please use python 3.6 or newer
vinfo=sys.version_info
assert vinfo.major >= 3 and vinfo.minor>=6
#NOTE: matplotlib 3.2 needs python 3.6, at least for easy pip install
'''
This is a simulation of the partial differential equations (PDE) covering gas flow through 
one single cylinder pressure swing adsorption unit for the C2(C2H2,C2H4,C2H6) separation.

NOTE: this version tries to follow: Haghpanah,...Farooq 2013 "Multiobjective
Optimization of a Four-Step Adsorption Process for Postcombustion CO2 Capture 
Via Finite Volume Simulation"  Ind. Eng. Chem. Research 2013, 52, pp 4249-4265.

See README for reference listings

We discretize the PDE in the spatial dimension, then use an ODE solver from scipy
to solve it

'''
#Local modules
#import orifice
#import params_farooq as params
import difference
import mylog
import util

#from params import create_param
#These global variables in this module will control operation
#(since the function that ode sovler' tries to solve doesnt take parameters other than x,t, therefore we must declare global varibale param and dg out side of the function to enable the acess of the param from the oademodel finition function internally)
global param
global dg
param=None
dg=None

class AttrDict(dict):
  def __init__(self, *args, **kwargs):
    super(AttrDict, self).__init__(*args, **kwargs)
    self.__dict__ = self
    
#These are attributes of state variable (keys also)


#All variables are nondimensionalized, normalized by proper values, as given
#by the norm dict:
def create_norm():
  #The normalization values, from Table 1 in Supplementary Info for the Ref paper
  # I don't understand the purpose of this function, It seems it does nothing but "renaming" some of the parameters from the params file,why bother to define another function?  
  norm=AttrDict()
  norm.P0=param.PH           # 300kPa      [Pa]
  norm.T0=293.15       #Room temp = 20C, [K]
  norm.qs0=param.qs0   #The normalization amount mol/kg
  norm.xi=norm.qs0     # qs,0  [mol/kg]
  norm.v0=param.norm_v0            #  [m/s]
  norm.z0=param.L      # Length of cylinder [m]
  norm.t0=param.norm_t0      # so t=time/norm.t0

  return norm

def create_dgroups(param):
  #these are the dimensionless coefficients needed, from SI of paper page 1
  #NOTE that some (like omega's) must be calculated at each state position and
  #are arrays j=1..N (COMMENT: would be better to define another func separately, which would be called many times by the ode solver in solving pdes,
  #while for partial and total mass balance, Pe and psi are simply constant that doesn't depend on state variables, 
  #therefore, I modfiy this function to take one parameter that is param and calculate Pe and psi only for once.
  #or just simply move these two groups calculation to the create_param func in the params.py file!))
  
  g=AttrDict()
  #for partial and total mass balances
  g.Pe=param.norm_v0*param.L/param.DL
  g.psi=param.R*param.norm_T0*param.qs0/param.PH*param.epsb*param.rho_s # qs0[mol/kg] 
  
  # for column energy balance
  #NOTE: probably a typo in Peh definition, units are wrong in paper (qs0 is missing, fixed)
  #NOTE: Need to calculate rho_g from the pressure by gas Law - TODO
  #we just just calculation for 1 ATM pressure, should be redone more accurately
  #num=param.Kz/(norm.v0/param.epsilon/param.L)
  #g.Peh=(param.epsilon*norm.v0*param.L)/(param.Kz/param.rho_g/param.Cpg)*norm.qs0
  #paren1=param.rho_s*param.Cps+norm.qs0*param.Cpa*(sum(state[xi] for xi in param.xnames)) # wrong parenthis position, fixed
  #g.omega1=param.Kz/(norm.v0*param.epsilon*param.L)/(param.epsb*paren1)
  #g.omega2=param.Cpg*norm.P0/param.R/norm.T0/(param.epsb*paren1)
  #g.omega3=param.Cpa*norm.qs0/paren1
  #g.omega4=2*param.hin*param.L/(param.rin*norm.v0)/(1-param.epsilon)/paren1
  #for xi in param.xnames:
  #  g['sigma{}'.format(xi)]=-norm.qs0/norm.T0*param.Hi[param.xanems.index(xi)]/(paren1)
  #g.sigmaB=-norm.qs0/norm.T0*param.HB/(paren1)

  # for wall energy balance
  #g.pi1=param.Kw/(param.rho_w*param.Cpw*norm.v0*param.L)
  #g.pi2=(2*param.rin*param.hin/(param.rout**2-param.rin**2)*
  #       param.L/(param.rho_w*param.Cpw*norm.v0))
  #g.pi3=(2*param.rout*param.hout/(param.rout**2-param.rin**2)*
  #       param.L/(param.rho_w*param.Cpw*norm.v0))

  #dimensionless constant for the product tank evolution
  #g.kprod=1e5/norm.P0*param.L**3/param.product_tank_volume
  return g
  
def compute_alpha(param, verbose=False):
  #compute the alpha for a given c [mol/m3] and qs [saturation mol/m3] at the
  #given pressures
  #k is the constant with dimensions 1/sec
  #we will be using data from a reference paper, not calculating with this
  #formula (where ci is a function of time and space, causing itertation when solving odes) from Haghpanah
  #k=(c/qs)*(15*param.epsilonp*param.Dp/param.rp**2)
  #k is the inverse time constant for the adsorption dynamics
  alphai=np.array(param.ki)*param.L/param.norm_v0 

  if verbose:
    print('dynamics: ki={:5.2f} 1/sec  alphai={:5.2f}   1/alphai={:6.2f}'.format(param.ki,alphai,1/alphai))
  #alpha is the dimensionless constant
  return alphai
  
def adsorp_equilibrium(yi,P,isothermtype='Langmuir-Freundlich'):
  #given yi (0<= yi <= 1) returns adsorption equilibrium in normalized amount
  #e.g. input: yi=[0.2,0.8], P=1, isothermtype='Langmuir-Freundlich' 
  #     (NOTE：1. input P is normalized pressure with a reference pressure of Param.PH=3.5e5 [Pa]
  #      2. currently only Langmuir-Freundlich is implemented)               
  #     output: array([0.00087539, 0.12614982]) normalized with Param.qs0=7.3 mol/kg

    if isothermtype == 'Langmuir-Freundlich':      
        #for those concentrations
        #convert to atm partial pressures
        # sometimes the odes solver would search to a very small negative value of yi such as -0.00145, in order to avoid negative ci when calling the power function (when exponent is not an integer, ci must be a positive value)
        # solution: we first check if the negative value is close to 0, if so we manually change that negative value to zero

        ci= np.array(yi) * P * param.PH/1e5          #in atm

        #these are dimensionless fractions
        #NOTE: can't have a negative pressure because the exponents are not integers
        with warnings.catch_warnings():
            warnings.filterwarnings('error')
        try:
            denom=1+np.sum(np.array(param.bi) * np.power(ci , np.array(param.ni)))
            
            fi=np.array(param.bi)* np.power(ci , np.array(param.ni))/denom
        except Exception:
            sys.stdout.write('WARNING: {} {}\n'.format(ci))
            raise 
        #equilibrium values at these pressures
        qeqi=fi*np.array(param.qsi)   # mol/kg
        return qeqi/param.qs0  #convert to mol/kg and then normalize

    if isothermtype == 'Loading_ratio_correlation': #with paramters (qsi,bi,ni,miu_i)  
        raise NameError("not yet implemented!")

    else:
        raise NameError("unkwon adsorption isotherm type!")
     
def adsorp_equilibrium_array(yi,P,isothermtype='Langmuir-Freundlich'):
  # input:  yi,dimension:[param.N * param.nocomponents]
  #         p, dimension: [param.N *1]
  # output: temp_dict, dimension:[param.N *param.nocomponents]
  temp_dict=np.ndarray(np.shape(yi)) # initial an empty ndarray of the same size as yi
  (yirow,yicol) = np.shape(yi) # get the dimension of input yi and P 
  (Prow,)=np.shape(P)

  assert Prow == yirow #row lens of yi and P must match
  assert yicol == param.nocomponents 

  # the following three loops call func "adsorp_equilibrium" to calculate state.xi under the given pressure state.P[n]
  # P.S. 
  # I don't know if there is a more elegant way of doing this, but this( is obiviously not very efficient due packing and unpacking numpy.arrays n*len(xnames)*len(ynames) times)
  # ,though is best solution that I come up with. Perhaps improvements could be done by adjusting the func "adsorp_equilibrium", so that its output is directly of dimension: [param.N * param.nocomponents]
  for n in range(Prow): 
    input=[]
    xout=[]
    for yname in param.ynames:
        input.append(yi[n][param.ynames.index(yname)])  
    xout = adsorp_equilibrium(input,P[n],isothermtype)
    for xname in param.xnames:
        temp_dict[n][param.xnames.index(xname)]=xout[param.xnames.index(xname)]   #solid phase, temp_dict[xname] == state.xi
  return temp_dict

#Set the initial values for state variables for 1 container
def init(param):
  #return the initialized state object
  #assume 1ATM, equilibrium concentrations of adsorbed components
  #print('initializing container state variables')
  state=AttrDict()

  for yi in param.ynames:
    state[yi]=np.ones((param.N,))*np.array(param.ini_yi[param.ynames.index(yi)])
    #Declare state variables state.yi [value type: numpy array; dimension: n*1] and initialize them with corresponding feed gas fractions of component i: param.ini_yi
    #e.g.  state{'yA': array([0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2]), 
    #            'yB': array([0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8])}
  
  for xi in param.xnames: # similarily,the "for" loop initialize state.xi for the i_th component
    state[xi]=np.ones((param.N,)) 
  
  state.P=np.ones((param.N,))*param.PL/param.norm_P0     #Pressure (total gas phase) 1.0=1 ATM
  state.T=np.ones((param.N,))*param.Ta/param.norm_T0     #Temp of column   # 1.0 = Room Temp (293K)
  #C=(state.P[0]*param.norm_P0)/R/(state.T[0]*param.norm_T0)*1e-6  # total conc [mol/cm3]
  #cA=C*state.yA[0]
  #cB=C*state.yB[0]

  state_yi_concatenated= np.column_stack([state[i] for i in param.ynames]) #pack up param.nocompoent state.yi into one bigger ndarray of dimension
  state_xi_concatenated = adsorp_equilibrium_array(state_yi_concatenated,state.P,'Langmuir-Freundlich')
  for xi in param.xnames:
      state[xi]=state_xi_concatenated[:,param.xnames.index(xi)]   #upacking the calculated solid phase adsorption amount state.xi 
  
  state.Tw=np.ones((param.N,))    #Wall temperature
  #states for the product tank
  return state

def compute_vel(data,param):
  #compute the velocities from the pressure.  The last one will be wrong
  #because we do not have the Pz0 and PzL pressure here.  But should be close
  #enough for plots - last one is the one that is wrong
  r,c=data.P.shape
  V=np.zeros_like(data.P)
  #make the last velocity match the next-to-last one
  for i in range(c):
    PzL=(3*data.P[-1,i]-data.P[-2,i])/2
    V[:,i]=-1/param.deltaZ*4/150*(param.epst)**2*param.rp**2*param.norm_P0/param.mu/param.norm_v0/param.L*difference.diffplus(data.P[:,i],PzL)
  return V

def oadesmodel(t,x):
  # oadesmodel = ordinary algebraic differential equations model
  # input: x is a column vector of dimension [ len(param.state_names)* param.N , 1]
  # output: dxdt time derivatives of state variables (of the same dimension as the input x)

  if param.bed == 1:#equations for single bed simulation
    # defined over 4 operation steps in sequence (i.e. pressurization,adsorption,blowdown,evacuation)
    # you could add in new operation steps by giving self defined boundary conditions and step time duration specified below and in the param.tnewstep in  params.py 
    state1=x_to_state(x, param)
    dstate1= AttrDict() # initialize an empty dictionary that will be storing derivatives of all the state variable
    
    bound=AttrDict() # initialize an empty dictionary that will be storing boundary conditions for different steps

    #set the limit for state variables(would result in oscillative yi profile)
    #for yi in param.ynames:
    #  state1[yi]=util.limit(state1[yi],0.000001, 1)
    #state1.P=util.limit(state1.P,0.2,  4)
    #state1.T=util.limit(state1.T,0.0001, 1.5)
    #state1.Tw=util.limit(state1.Tw,0.0001, 1.5)


    #NOTE: state.any is a arrray with index from 1 to N: [1,2,...N-1,N]
    # other points namely the two boundary points f_0.5 and f_N+0.5 are set by different boundary conditions in each step 
    # other points exceed the boundaries that are f_0 and f_N+1 are calculated using half cell approximation see eq (32)
    # half cell approximation: eq(32)
    # f0 = 2*f0.5-f1
    # fN+1 =2*fN+0.5 -fN
    if 0 <= t <= param.norm_tpre: #pressurization step
      #Boundary condition at Z=0(denoted Xz0),and Z=L (denoted XzL)
      bound.Pz0 = 1/param.norm_P0 *(param.PH - (param.PH-param.PL)*math.exp(-param.lamda_pre* t *param.norm_t0)) # t(=tao in the reference equation (57)) is already normalized
      bound.PzL = state1.P[-1] # lhs:  Pzl=P_N+0.5, rhs: state.P[-1] = P_N
      bound.Tz0 = state1.T[0] # dummy boundary conditions, engergy balances are not actually computed
      bound.TzL = state1.T[-1] #same as above
      bound.Twz0 = state1.Tw[0] #same as above
      bound.TwzL = state1.Tw[-1] #same as above
      bound.vz0 = -2/param.deltaZ*(4/150*param.epst**2)*(param.rp**2)*param.norm_P0/param.mu/param.norm_v0/param.L*(state1.P[0] - bound.Pz0) #inlet veloctiy computed with P1 and P0.5
      for yi in param.ynames:
        bound['{}z0'.format(yi)]=(state1[yi][0]+param.feed_yi[param.ynames.index(yi)]*bound.vz0*dg.Pe*param.deltaZ/2)/(1+bound.vz0*dg.Pe*param.deltaZ/2) # inlet composition is dependent on inlet velocity vz0, which further depends on P1 and P0.5
        bound['{}zL'.format(yi)]= state1[yi][-1]

      forwardflow= True   
    if param.norm_tpre < t <= param.norm_tads: # adsorption step
      #Boundary condition at Z=0(denoted Xz0),and Z=L (denoted XzL)
      bound.vz0 = 1 #inlet veloctiy fixed to vfeed
      bound.Pz0 = state1.P[0] + bound.vz0 * param.deltaZ* param.mu *param.norm_v0 *param.L /(2*4/150 * param.epst**2 * param.rp**2 *param.norm_P0) # t(=tao in the reference equation (57)) is already normalized
      bound.PzL = 1 # = P_H
      bound.Tz0 = state1.T[0] # dummy boundary conditions, engergy balances are not actually computed
      bound.TzL = state1.T[-1] #same as above
      bound.Twz0 = state1.Tw[0] #same as above
      bound.TwzL = state1.Tw[-1] #same as above
      for yi in param.ynames:
        bound['{}z0'.format(yi)]=(state1[yi][0]+param.feed_yi[param.ynames.index(yi)]*bound.vz0*dg.Pe*param.deltaZ/2)/(1+bound.vz0*dg.Pe*param.deltaZ/2) # inlet composition is dependent on inlet velocity vz0, which is a constant now
        bound['{}zL'.format(yi)]= state1[yi][-1]

      forwardflow= True   
    if param.norm_tads < t <= param.norm_tblw: # blowdown step
      #Boundary condition at Z=0(denoted Xz0),and Z=L (denoted XzL)
      bound.Pz0 = state1.P[0]   # t(=tao in the reference equation (57)) is already normalized
      bound.PzL = 1/param.norm_P0 *(param.PI + (param.PH-param.PI)*math.exp(-param.lamda_blw* (t-param.norm_tads) *param.norm_t0))# lhs:  Pzl=P_N+0.5, rhs: state.P[-1] = P_N
      bound.Tz0 = state1.T[0] # dummy boundary conditions, engergy balances are not actually computed
      bound.TzL = state1.T[-1] #same as above
      bound.Twz0 = state1.Tw[0] #same as above
      bound.TwzL = state1.Tw[-1] #same as above
      bound.vzL = -2/param.deltaZ*(4/150*param.epst**2)*(param.rp**2)*param.norm_P0/param.mu/param.norm_v0/param.L*(bound.PzL - state1.P[-1]  )
      for yi in param.ynames:
        bound['{}z0'.format(yi)]= state1[yi][0]
        bound['{}zL'.format(yi)]= state1[yi][-1]

      forwardflow= True     
    if param.norm_tblw < t <= param.norm_tend: # evacuation step
      #this step is a backward flow,  if you modify the flow direction, don't forget to modify the function which computes purity and recovery, especially mole_out,i_evac = f( s_0.5 ), s=v,y,P,T, since now the outlet is at z=0
      #Boundary condition at Z=0(denoted Xz0),and Z=L (denoted XzL)
      bound.Pz0 = 1/param.norm_P0 *(param.PL + (param.PI-param.PL)*math.exp(-param.lamda_eva* (t-param.norm_tblw) *param.norm_t0))  # t(=tao in the reference equation (57)) is already normalized
      bound.PzL = state1.P[-1] # lhs:  Pzl=P_N+0.5, rhs: state.P[-1] = P_N
      bound.Tz0 = state1.T[0] # dummy boundary conditions, engergy balances are not actually computed
      bound.TzL = state1.T[-1] #same as above
      bound.Twz0 = state1.Tw[0] #same as above
      bound.TwzL = state1.Tw[-1] #same as above
      bound.vz0 = -2/param.deltaZ*(4/150*param.epst**2)*(param.rp**2)*param.norm_P0/param.mu/param.norm_v0/param.L*(state1.P[0] - bound.Pz0) #inlet veloctiy computed with P1 and P0.5
      for yi in param.ynames:
        bound['{}z0'.format(yi)]= state1[yi][0]
        bound['{}zL'.format(yi)]= state1[yi][-1]

      forwardflow= False
    elif t > param.norm_tend:
      raise NameError('exceeded the time bound')
  # discretized form of equation (20) -(24) NOTE:engergy balance (23)and(24) is not yet implemented!
    # generate functions for evaluating j+0.5 and j-0.5 according to param.mode
    jplus5,jminus5 = difference.gen_j5_functions(param.mode,forwardflow) 
    # see equation (33)
    vplus5=-1/param.deltaZ*4/150*param.epst**2*param.rp**2*param.norm_P0/param.mu/param.norm_v0/param.L*difference.diffplus(state1.P, bound.PzL)
    vminus5=-1/param.deltaZ*4/150*param.epst**2*param.rp**2*param.norm_P0/param.mu/param.norm_v0/param.L*difference.diffminus(state1.P, bound.Pz0)

    # mass transfer rate equation(22), first calculate equilibrium adsorption xistar amount under yi and P ,then compute eq(22)


    xistar = adsorp_equilibrium_array(np.vstack([state1[yi] for yi in param.ynames]).T, state1['P'], param.isomodel ) # transpose is needed

    for xi in param.xnames:
      dstate1['d{}'.format(xi)] = compute_alpha(param)[param.xnames.index(xi)] * (xistar[:, param.xnames.index(xi)] - state1[xi])

    dstate1['dT']= zeros_like(state1.T) #energy balance equiation (23) assumed to be constant
    dstate1['dTw'] =zeros_like(state1.Tw) #wall energy balance equation (24) assumed constant 

    # total mass balance equation (21)
    dstate1['dP']=(-state1.T/param.deltaZ*(jplus5(state1.P/state1.T,bound.Pz0/bound.Tz0,bound.PzL/bound.TzL)*vplus5-
                          jminus5(state1.P/state1.T,bound.Pz0/bound.Tz0,bound.PzL/bound.TzL)*vminus5)
                          -dg.psi*state1.T*(sum((dstate1['d{}'.format(xi)] for xi in param.xnames))) #sum(dxi) == "sum((state['d{}'.format(xname)] for xname in param.xnames))""
                          +state1.P/state1.T*dstate1.dT)
    
    
    # component mass balance equation (20) for each species   
    for yi in param.ynames:
      dstate1['d{}'.format(yi)]=(1/dg.Pe*state1.T/state1.P/param.deltaZ*(
                    jplus5(state1.P/state1.T,bound.Pz0/bound.Tz0,bound.PzL/bound.TzL) * difference.diffplus(state1[yi],bound['{}zL'.format(yi)])/param.deltaZ-
                    jminus5(state1.P/state1.T,bound.Pz0/bound.Tz0,bound.PzL/bound.TzL) * difference.diffminus(state1[yi],bound['{}z0'.format(yi)])/param.deltaZ)
                    -state1.T/state1.P/param.deltaZ*(
                    jplus5(state1[yi]*state1.P/state1.T, bound['{}z0'.format(yi)]*bound.Pz0/bound.Tz0, bound['{}zL'.format(yi)]*bound.PzL/bound.TzL)*vplus5 -
                    jminus5(state1[yi]*state1.P/state1.T, bound['{}z0'.format(yi)]*bound.Pz0/bound.Tz0, bound['{}zL'.format(yi)]*bound.PzL/bound.TzL)*vminus5)
                    -dg.psi*state1.T/state1.P*dstate1['d{}'.format(param.xnames[param.ynames.index(yi)])] #Dxidt (eg. yname == yA (param.yis == [yA,yB], param.xnames=[xA,xB]) -> "dstate['d{}'.format(param.xnames(param.ynames.index(yi)))]" gives  "dstate[dxA]")
                    -state1[yi]/state1.P* dstate1.dP  #DPdt
                    +state1[yi]/state1.T* dstate1.dT) #DTdt
    

    


    
    
    dstate = np.hstack([dstate1['d{}'.format(item)] for item in param.state_names]) # adsorbed amount
    

    return dstate

  if param.bed == 2 :#equations for two beds simulation with connections
    return myfun(t,x)

def x_to_state(x, param):
  if param.bed == 2:
    state1=AttrDict()
    state2=AttrDict()
    prod=AttrDict()
    q=len(param.state_names)
    #print('x_to_state, x size = {}'.format(x.size))
    start=0
    for j,state in enumerate((state1, state2)):
      for i,name,sz in zip(range(q),param.state_names,param.state_sizes):
        #print('{} {} size {}'.format(i,name,sz))
        state[name]=x[start:start+sz]
        start+=sz
    #now product tank
    prod.Pprod=x[start:start+1]
    start+=1
    prod.yprod=x[start:start+1]
    return state1, state2, prod
  
  if param.bed == 1:
    state1=AttrDict()
    start=0
    for i,name,sz in zip(range(len(param.state_names)),param.state_names,param.state_sizes):

      state1[name]=x[start:start+sz]
      start+=sz
    return state1

def data_to_state(x,param):
  # convert output of simulation (big matrix) to 2 data variables with attributes
  # input: x : 
  #                         [P0,P1,...PN,T0,T1,...TN,Tw0,Tw1, ...,TwN,xA0,..xAN,xB0,..xBN,yA0,...yAN,yB0,...,yBN].transpose, 
  #                         where every variables(P0,P1,..., PN, T0....) is a row vector of dimension param.N,
  #                         the total dimension is [3+2*param.nocomponents, param.N]
  # output: state1:
  #                        state1.P = [P0,P1,...,PN].T
  #                        state1.T = [T0,T1,...,TN].T
  #                        ...
  #                        state1.xB =[xB0,xB1,...,xBN].T
  #                        ...
  ## convert bunch.y back to a dictionary, where state variables such as P, T, Tw, yi, and xi can be aceessed by state1.[state_varibale's name]
  if param.bed == 2:
    state1=AttrDict()
    state2=AttrDict()
    prod=AttrDict()
    N=param.N
    q=len(param.state_names)
    start=0
    for j,state in enumerate((state1, state2)):
      for i,name,sz in zip(range(q),param.state_names,param.state_sizes):
        state[name]=x[start:start+sz,:]
        start+=sz
    #now product tank
    prod.Pprod=x[start:start+1,:]
    start+=1
    prod.yprod=x[start:start+1,:]
    return state1, state2, prod

  if param.bed == 1:
    state1=AttrDict()
    start=0
    for i,name,sz in zip(range(len(param.state_names)),param.state_names,param.state_sizes):

      state1[name]=x[start:start+sz,:]
      start+=sz
    return state1

def print_state(state):
  print('P:{}'.format(state.P))
  print('T:{}'.format(state.T))
  print('Tw:{}'.format(state.Tw))
  for xi in param.xnames:
    print('{}:{}'.format(xi,state[xi]))
  for yi in param.xnames:
    print('{}:{}'.format(yi,state[yi]))


####
#This is the main function of this module, should be the only one exposed
#to the outside

  
def simulate_singlebed(localparam,output_dir='./outcome'):

  #paramsmodule=importlib.import_module(params_file)
  #CORE FUNCTION
  #START
  global param
  global dg

  param = localparam  # set parameters
  dg = create_dgroups(param) # compute constant groups in model equations, if energy balances are computed, so that some of the groups no longer remain constant it would be better to throw this process into the model funtion
  x0_all= init(param)
  x0 = np.hstack([x0_all[item]  for item in param.state_names]) # initialization state variables e.g. for a two component mixture it would be ['P','T','Tw', 'xA','xB','yA','yB']
  ev=np.linspace(0,param.norm_tend,param.tN) # set the intergration time span with a given uniform time step
  
  #integrate the odae functions oadesmodel (discretized using finite volume method)
  bunch=scipy.integrate.solve_ivp(oadesmodel, (0, param.norm_tend), x0, vectorized=False, t_eval=ev,
                                  method='LSODA')
  #CORE FUNCTION
  #END


  #data post prossesing, mainly data storage and drawing plots
  plot_data(bunch, param, output_dir)


  return 0

def cyclic_steady_state(localparam,output_dir='./outcome', maxcycle = 50, tolerance = 1e-4):
  # this function simulate mutiple psa process until cyclic steady state is reached, and store the data generated during simulation
  # return: "status" variable has the following data structure, type(status): AttriDict 
  #          
  #          status.param: "params used in simulation" [type: AttriDict ] e.g. status.param.nocomponents = 3 (corresponds to what you've set in the params_file)
  #          status.snap: "snapshot at each cycle" [type: list] e.g. status.snap[]
  # 
  
  #paramsmodule=importlib.import_module(params_file)
  #CORE FUNCTION
  #START
  global param
  global dg
  status =AttrDict()
  status.snap = AttrDict()
  param = localparam  # set parameters
  status.param = param # bind the particular param file to the output status,this would be convenient when we process the data after simuation for ploting by loading pickles
  dg = create_dgroups(param) # compute constant groups in model equations, if energy balances are computed, so that some of the groups no longer remain constant it would be better to throw this process into the model funtion 
  
  
   # number of max iterations
  for cycle in range(maxcycle):

    if cycle == 0:  
      x0_all= init(param)
      x0 = np.hstack([x0_all[item]  for item in param.state_names]) # initialization state variables e.g. for a two component mixture it would be ['P','T','Tw', 'xA','xB','yA','yB']
    else:
      x0 = xend

    ev=np.linspace(0,param.norm_tend,param.tN) # set the intergration time span with a given uniform time step

    #integrate the odae functions oadesmodel (discretized using finite volume method)
    bunch =scipy.integrate.solve_ivp(oadesmodel, (0, param.norm_tend), x0, vectorized=False, t_eval=ev,
                                  method='LSODA')
    
    x0= bunch.y[:,0] 
    xend = bunch.y[:,-1]
    tol= tolerance * np.ones_like(x0) #iteration tolerance default 1e-4

    status.snap[cycle] = bunch #pack all the data produced by the scipy ode solver in nth cycle into the status AttrDict, which will be used later for plotting
    print('current cylce {}'.format(cycle))
    #print(xend - x0)
    if (np.abs(xend - x0) <= tol).all() : # the css condition whether the difference of state varibales between the beginning of a cylce and the end of a cycle is less than the set tolerance
      status.cycle = cycle # status.cycle is the "final" converged cycle number -1, because python's count starts at 0
      status.converged = True
      print('Converged in the {} cylce'.format(cycle))
      
      return bunch, status
    elif (cycle +1) == maxcycle:
      status.cycle = cycle
      status.converged = False
      return bunch, status

def plot_data(bunch, param, output_dir): 
  # record data and plot
  basedir=os.path.dirname(output_dir)
  if basedir is None:
    dbg=mylog.Log(outdir=None,descr='psa')
    OUTDIR='./'
  else:
    OUTDIR=basedir  # for log files and png files
    dbg=mylog.Log(outdir=OUTDIR,descr='psa')
  #change the print function in the module to print to the log file
  fileprint=dbg.print


  data=AttrDict()
  data.data1=data_to_state(bunch.y,param)
  print(np.shape(bunch.y))
  #print(np.shape(bunch.y))
  data.data1.t = bunch.t
  data.data1.param = param
  fileprint('time :{}'.format(bunch.t))
  fileprint('parameters: {}'.format(param))
  for item in param.state_names:
    fileprint('State variable {}: {}'.format(item, data.data1[item]))

  data.data1.v = compute_vel(data.data1, param) #it seems that we need to pass the entire data dictionary instead of just P, I dont understand
  fileprint('State variable velocity: {}'.format(item, data.data1['v']))

  plots.plot_singlebed(data, out_place=output_dir)
  return 0  
    
def purity_recovery(bunch,param):
  # calculate the purity of the product stream for a 4-step psa process with boundary conditions set the same as in the oadesmodel function
  # input : bunch is a standard ode output
  # return: obj type:[AttrDict]
  #             datastructure: obj['purity_ads'][yi], obj['purity_blw'][yi], obj['purity_eva'][yi], obj[recovery][yi]
  mole= AttrDict()
  mole['z0_tpre']=AttrDict()
  mole['z0_tads']=AttrDict()
  mole['zL_tads']=AttrDict()
  mole['zL_tblw']=AttrDict()
  mole['z0_teva']=AttrDict()

  obj=AttrDict()
  obj['purity_ads']=AttrDict()
  obj['purity_blw']=AttrDict()
  obj['purity_eva']=AttrDict()
  obj['recovery']=AttrDict()
  
  


  #bunch.t == np.linspace(0,param.norm_tend,param.tN)#bunch.t
  data = data_to_state(bunch.y,param) # convert bunch.y back to a dictionary, where state variables such as P, T, Tw, yi, and xi can be aceessed by data.[state_varibale's name][z,t]
  gnorm = param.norm_P0 * param.norm_v0 * param.epsilon *param.area /param.R /param.norm_T0 #normalized constant group
  dg = create_dgroups(param)
  for yi in param.ynames:
    # pressurization step 
    # z=0
    t_pre=bunch.t[0:param.index_tpre]
    pre_Pz0= 1/param.norm_P0 *(param.PH - (param.PH-param.PL)*np.exp(-param.lamda_pre* t_pre *param.norm_t0)) 
    pre_vz0= -2/param.deltaZ*(4/150*param.epst**2)*(param.rp**2)*param.norm_P0/param.mu/param.norm_v0/param.L*(data.P[0, 0:param.index_tpre] - pre_Pz0)
    pre_Tz0= data.T[0,0:param.index_tpre]
    #integration using Simpson's Rule
    mole['z0_tpre'][yi] = gnorm * scipy.integrate.simpson(data[yi][0,0:param.index_tpre] *pre_vz0* pre_Pz0 / pre_Tz0, t_pre)
    
    # adsorption step
    # z=0
    t_ads=bunch.t[param.index_tpre+1:param.index_tads]
    ads_vz0 = param.vfeed/param.norm_v0 * np.ones_like(t_ads)
    ads_Pz0 = data.P[0,param.index_tpre+1:param.index_tads] + ads_vz0 * param.deltaZ* param.mu *param.norm_v0 *param.L /(2*4/150 * param.epst**2 * param.rp**2 *param.norm_P0)
    ads_Tz0 = data.T[0,param.index_tpre+1:param.index_tads]

    # z=L
    ads_PzL = param.PH/param.norm_P0 * np.ones_like(t_ads) # outlet pressure in the adsorption is set to operating at PH
    ads_TzL = data.T[-1,param.index_tpre+1:param.index_tads]
    ads_vzL = -2/param.deltaZ*(4/150*param.epst**2)*(param.rp**2)*param.norm_P0/param.mu/param.norm_v0/param.L*(ads_PzL - data.P[-1,param.index_tpre+1:param.index_tads]  )
    #integration using Simpson's Rule
    mole['z0_tads'][yi] = gnorm * scipy.integrate.simpson((data[yi][0,param.index_tpre+1:param.index_tads]+param.feed_yi[param.ynames.index(yi)]*ads_vz0*dg.Pe*param.deltaZ/2)/(1+ads_vz0*dg.Pe*param.deltaZ/2) *ads_vz0* ads_Pz0 / ads_Tz0, t_ads) 
    mole['zL_tads'][yi]=  gnorm * scipy.integrate.simpson(data[yi][-1,param.index_tpre+1:param.index_tads] *ads_vzL* ads_PzL / ads_TzL, t_ads)
    
    #depressurization step
    # z=L
    t_blw=bunch.t[param.index_tads+1:param.index_tblw]
    blw_PzL = 1/param.norm_P0 *(param.PI + (param.PH-param.PI)*np.exp(-param.lamda_blw* (t_blw-param.norm_tads) *param.norm_t0))# lhs:  Pzl=P_N+0.5, rhs: state.P[-1] = P_N
    blw_TzL = data.T[-1, param.index_tads+1:param.index_tblw]
    blw_vzL = -2/param.deltaZ*(4/150*param.epst**2)*(param.rp**2)*param.norm_P0/param.mu/param.norm_v0/param.L*(blw_PzL - data.P[-1, param.index_tads+1:param.index_tblw]  )
    mole['zL_tblw'][yi]= gnorm * scipy.integrate.simpson(data[yi][-1, param.index_tads+1:param.index_tblw] *blw_vzL* blw_PzL / blw_TzL, t_blw)

    #evacuation step
    # z=0
    #assert (param.index_tend == len(bunch.t)-1)
    t_eva=bunch.t[param.index_tblw+1:param.index_tend]
    eva_Pz0 = 1/param.norm_P0 *(param.PL + (param.PI-param.PL)*np.exp(-param.lamda_eva* (t_eva-param.norm_tblw) *param.norm_t0))  # t(=tao in the reference equation (57)) is already normalized

    eva_Tz0 = data.T[0, param.index_tblw+1:param.index_tend] # dummy boundary conditions, engergy balances are not actually computed
    eva_vz0 = -2/param.deltaZ*(4/150*param.epst**2)*(param.rp**2)*param.norm_P0/param.mu/param.norm_v0/param.L*(data.P[0, param.index_tblw+1:param.index_tend] - eva_Pz0) #inlet veloctiy computed with P1 and P0.5
    #print(eva_vz0)
    if not ((eva_vz0 <= 0).all()): #check if all the velocity is negative, if so apply abs to eva_vz0, because we need a postive number for molar amount
      print("error some velocities at the outlet z=0 of the blowdown step is positive! ")
      return None
    else:
      eva_vz0= np.abs(eva_vz0)
    mole['z0_teva'][yi]= gnorm * scipy.integrate.simpson(data[yi][0, param.index_tblw+1:param.index_tend] *eva_vz0* eva_Pz0 / eva_Tz0, t_eva)

  #now we can compute purity and recovery using eq(60) andeq(61)
  for yi in param.ynames:
    obj['purity_eva'][yi] = mole['z0_teva'][yi]/ np.sum([mole['z0_teva'][yi] for yi in param.ynames])
    obj['purity_ads'][yi] = mole['zL_tads'][yi]/ np.sum([mole['zL_tads'][yi] for yi in param.ynames])
    obj['purity_blw'][yi] = mole['zL_tblw'][yi]/ np.sum([mole['zL_tblw'][yi] for yi in param.ynames])
    for where in param.collect:
      obj['recovery'][yi] = mole['{}'.format(param.collect[param.ynames.index(yi)])][yi] /np.sum([mole['z0_tads'][yi] , mole['z0_tpre'][yi]])

   
  return obj

def optimization_psa():
  #this function trys to optimze the entire psa process on a given adsorbent by adjusting follwoing parameters:

  pass  
  
class PSA_optimization(ElementwiseProblem):

  def __init__(self, mofdict, **kwargs):
    #optimization varibales
    # x[0]=tpre [t] : pressurization time duration
    # x[1]=tads [t] : adsorption time duration 
    # x[2]=tblw [t] : desorption time duration 
    # x[3]=teva [t] : evacuation time duration
    # x[4]=feed_pressure [bar] : feed pressure/the highest pressure/adsorption pressure
    # x[5]=vfeed [m/s]: feed velocity
    # x[6]=PI [Pa] : intermediate pressure
    # x[7]=PL [Pa] : evacuation pressure/ the lowest pressure
    super().__init__(n_var=8,
                     n_obj=2,
                    n_constr=2,
                    xl=np.array([10,10, 10, 10, 1, 0.1, 1e4, 1e4]),
                    xu=np.array([20,100,100,100,10,2, 1e6, 1e6]),
                    **kwargs)
    # this is the parameters for TJT-100
    self.param = collections.defaultdict()
    
    if mofdict == None: #default TJT-100
      self.param['bi']=[3.499, 2.701, 7.185]   # [1/bar] @ Room temp
      self.param['qsi']=[6.29,4.59, 4.94]  #saturation constant for component i  [mol/kg]=[mmol/g] (at infinite pressure)
      self.param['qs0'] = 6.29 #[mol/kg]
      self.mofname="TJT-100"
    else:
      self.param['bi']=mofdict['bi']
      self.param['qsi']=mofdict['qsi'] 
      self.param['qs0'] = mofdict['qs0']
      self.mofname = mofdict['name']


  def _evaluate(self, x, out, *args, **kwargs):
    # given a set of input parameter x (element order see the comments in funcction __ini__)
    # f1:purity and f2:recovery
    # create localparameters with overloaded x input

    paramoverload =collections.defaultdict() #create a dictionary where parameters that need to be overloaded to create_param function
    paramoverload['tpre'] = x[0]
    paramoverload['tads'] = x[1]
    paramoverload['tblw'] = x[2]
    paramoverload['teva'] = x[3]
    paramoverload['feed_pressure'] = x[4]
    paramoverload['vfeed'] =x[5]
    paramoverload['PI'] = x[6]
    paramoverload['PL'] = x[7]
    # overload mof parameters(this is a test )
    paramoverload['bi'] = self.param['bi']
    paramoverload['qsi'] = self.param['qsi']
    paramoverload['qs0'] = self.param['qs0']
    # create a dict of local parameters
    if paramoverload['PL'] < paramoverload['PI'] and (paramoverload['feed_pressure']*1e5 > paramoverload['PI']): # check if PL<PI<PH
      localparam= params.create_param(paramoverload) 

      finalcyle, status = cyclic_steady_state(localparam, tolerance = 3e-3) # simuate until css is reached

      if status.converged == True: # if ccs is reached then compute purity and recovery
        obj = purity_recovery(finalcyle, localparam) 
        if obj: #if purity and recovery are correctly computed then set f1 and f2
          f1 = -1*obj['purity_ads']['yB'] # yB = see params.py; here it denote ethene
          f2 = -1*obj['recovery']['yB']
        else:
          f1= 0
          f2= 0
      else:
        f1= 0
        f2= 0
    else:
      f1= 0
      f2= 0   


    g1 = x[7] - x[6]  # PL <= PI
    g2 = x[6] - x[4]*1e5 #PI <= PH

    out["F"] = [f1, f2]
    out["G"] = [g1, g2]


    




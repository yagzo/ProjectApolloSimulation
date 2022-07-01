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
import sys
import math
import collections
from telnetlib import STATUS
from traceback import print_tb
from xml.dom.minidom import Attr


import matplotlib.pyplot as plt
import numpy as np
from numpy import concatenate
from numpy import zeros_like
from numpy import column_stack
from numpy import hstack
import numpy
import scipy
import scipy.integrate
import matplotlib
import termios
import tty
import pprint
import argparse
import time
import datetime
import os
import warnings
import importlib
import pickle
##local modules
import plots

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
param=None
#param=create_param(None)
norm=None
dg=None
#for start time for debug output
tstart=0
count=0   # number of times myfun is called


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
  print('initializing container state variables')
  state=AttrDict()
  R=param.R
  N=param.N

  for yi in param.ynames:
    state[yi]=np.ones((param.N,))*np.array(param.ini_yi[param.ynames.index(yi)])
    #Declare state variables state.yi [value type: numpy array; dimension: n*1] and initialize them with corresponding feed gas fractions of component i: param.ini_yi
    #e.g.  state{'yA': array([0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2]), 
    #            'yB': array([0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8])}
  
  for xi in param.xnames: # similarily,the "for" loop initialize state.xi for the i_th component
    state[xi]=np.ones((param.N,)) 
  
  state.P=np.ones((N,))*param.PL/param.norm_P0     #Pressure (total gas phase) 1.0=1 ATM
  state.T=np.ones((N,))*param.Ta/param.norm_T0     #Temp of column   # 1.0 = Room Temp (293K)
  C=(state.P[0]*param.norm_P0)/R/(state.T[0]*param.norm_T0)*1e-6  # total conc [mol/cm3]
  cA=C*state.yA[0]
  cB=C*state.yB[0]

  state_yi_concatenated= np.column_stack([state[i] for i in param.ynames]) #pack up param.nocompoent state.yi into one bigger ndarray of dimension
  state_xi_concatenated = adsorp_equilibrium_array(state_yi_concatenated,state.P,'Langmuir-Freundlich')
  for xi in param.xnames:
      state[xi]=state_xi_concatenated[:,param.xnames.index(xi)]   #upacking the calculated solid phase adsorption amount state.xi 
  
  state.Tw=np.ones((N,))    #Wall temperature
  #states for the product tank
  return state

def rate_of_change_vel(state,vin=None,vout=None,
                       in_yA=None,in_yB=None,
                       out_yA=None,out_yB=None,verbose=False,
                       jplus5=None,jminus5=None):
  #compute the rates of change at current state
  #v0=input velocity of feed gas (dimensionless)
  #vL=output velocity of gas  (dimensionless)
  #Note, because paper has j=1...N, we need to be careful of indices
  #Also, paper has i=0 <==> we use A; i=1 <==> we use B for gas components
  #NOTE: out_yA,out_yB are relevant if vout<0, as those are the boundary conditions
  #for the gas coming into the right size (normally output port), this happens
  #during cross-venting
  #Compute dxA, dxB, dP, dT (meaning dT/dt, etc.)
  #jplus5 and jminus5 are the functions to use to compute these quantities, they
  #change because the direction flow changes them (forward vs. backward differences)
  #and we might want to use different difference schemes (van Leer, WENO, etc.)
  assert jplus5
  assert jminus5
  diffplus=difference.diffplus
  diffminus=difference.diffminus
  dstate=AttrDict()
  if verbose:
    print('***rate of change vel')
  deltaZ=param.deltaZ
  #calculate velocity profile.  There will be 2 versions:
  #j+.5 j=[1.5,2.5,...,N+.5] and
  #j-.5 j=[.5, 1.5, ..., N-.5]
  #we these instead of evaluating vj at the .5 points with the finite volume
  #calculations
  
  #Boundary conditions
  P1=state.P[0]
  PN=state.P[-1]
  T1=state.T[0]
  TN=state.T[-1]

  #Boundary condition at Z=0
  Pz0=P1+(vin*deltaZ/2)/(4/150*param.epst**2*param.rp**2*param.norm_P0/param.mu/param.norm_v0/param.L)
  if vin>0:
    #weighted average of input y and the existing y
    for yi in param.ynames:
      dstate['{}z0'.format(yi)]=(state[yi][0]+in_yA*vin*dg.Pe*deltaZ/2)/(1+vin*dg.Pe*deltaZ/2)
      #yBz0=(state.yB[0]+in_yB*vin*dg.Pe*deltaZ/2)/(1+vin*dg.Pe*deltaZ/2)
  else:
    #out flow from input side, so different boundary conditions
    for yi in param.ynames:
      dstate['{}z0'.format(yi)]==state[yi][0]
      #yBz0=state.yB[0]
  Tz0=T1
  Twz0=param.Ta/norm.T0
  #print('yAz0={} yA[0]={} yAfeed={}'.format(yAz0,state.yA[0],yAfeed))
  if verbose:
    print()
    print('Pz0 {} P1 {}'.format(Pz0,P1))
  
  #Boundary condition at Z=1 (L m)
  PzL=PN-(vout*deltaZ/2)/(4/150*param.epst**2*param.rp**2*param.norm_P0/param.mu/param.norm_v0/param.L)
  if vout>0:
    for yi in param.ynames:
      dstate['{}zL'.format(yi)]=state[yi][-1]
    #yBzL=state.yB[-1]
  else:
    #in flow from the output side:
    for yi in param.ynames:
      dstate['{}zL'.format(yi)]=(state[yi][-1]+out_yA*abs(vout)*dg.Pe*deltaZ/2)/(1+abs(vout)*dg.Pe*deltaZ/2)
    #yBzL=(state.yB[-1]+out_yB*abs(vout)*dg.Pe*deltaZ/2)/(1+abs(vout)*dg.Pe*deltaZ/2)
  if (verbose):
    print('rate of change yA[-1] {:6.2f} out_yA {:6.2f} yAzL = {:6.2f}'.format(state.yA[-1],out_yA,yAzL))
  TzL=TN
  TwzL=param.Ta/norm.T0

  #compute the velocities - and double check that pressure difference gives same answer
  vplus5=-1/deltaZ*4/150*param.epst**2*param.rp**2*norm.P0/param.mu/norm.v0/param.L*diffplus(state.P,PzL)
  vminus5=-1/deltaZ*4/150*param.epst**2*param.rp**2*norm.P0/param.mu/norm.v0/param.L*diffminus(state.P,Pz0)
  if verbose:
    print('P = {}'.format(state.P[0:6]))
    print('vplus5 {}'.format(vplus5[0:6]))
    print('vminus5 {}'.format(vminus5[0:6]))
    print(vplus5[1], vout)
    print('comparing vplus5[-1] {:8.4f} to vout {:8.4f}'.format(vplus5[-1],vout))
    print('comparing vminus5[0] {:8.4f} to vin {:8.4f}'.format(vminus5[0],vin))
  
  #now calculate the derivatives with respect to time
  #total mass balance
  g=create_dgroups(norm, param, state)  #some things depend on current state
  #solid phase dynamics
  if param.adsorb:
    alphai=compute_alpha()
    xis=adsorp_equilibrium_array(np.column_stack(state[yname] for yname in param.ynames), state.P) #xis: equilibrium adsorption amount
    for xi in param.xnames:
      dstate['d{}'.format(xi)]=alphai[param.xnames.index(xi)]*(xis[:,param.xnames.index(xi)]-state[xi]) #dxi:calculated by the linear driving force equation
  else:
    #no adsorption occurs
    for xi in param.xnames:
      dstate['d{}'.format(xi)]=zeros_like(state[xi])
  
  dT=np.zeros_like(state.T)
  dTw=np.zeros_like(state.T)
  
  dP=(-state.T/deltaZ*(jplus5(state.P/state.T,Pz0/Tz0,PzL/TzL)*vplus5-
                       jminus5(state.P/state.T,Pz0/Tz0,PzL/TzL)*vminus5)
      -dg.psi*state.T*(sum((dstate['d{}'.format(xi)] for xi in param.xnames))) #sum(dxi) == "sum((state['d{}'.format(xname)] for xname in param.xnames))""
      +state.P/state.T*dT)
  #Column energy balance
  #dT=np.zeros_like(state.T)
  #dT=(dg.omega1/deltaZ*(diffplus(state.T,TzL)/deltaZ-diffminus(state.T,Tz0)/deltaZ)-
  #    dg.omega2/deltaZ*(vplus5*jplus5(state.P,Pz0,PzL)-vminus5*jminus5(state.P,Pz0,PzL))-
  #    dg.omega3*state.T*(dxA+dxB)+
  #    (dg.sigmaA*dxA+dg.sigmaB*dxB)-
  #    dg.omega4*(state.T-state.Tw)-
  #    dg.omega2*dP)
      
  #Wall energy balance
  #dTw=(dg.pi1/deltaZ*(diffplus(state.Tw,TwzL)/deltaZ-diffminus(state.Tw,Twz0)/deltaZ)+
  #     dg.pi2*(state.T-state.Tw)-
  #     dg.pi3*(state.Tw-param.Ta/norm.T0))
  for yi in param.ynames:
    dstate['d{}'.format(yi)]=(1/dg.Pe*state.T/state.P/deltaZ*(
      jplus5(state.P/state.T,Pz0/Tz0,PzL/TzL)*diffplus(state[yi],yAzL)/deltaZ-
      jminus5(state.P/state.T,Pz0/Tz0,PzL/TzL)*diffminus(state[yi],dstate['{}z0'.format(yi)])/deltaZ)
        -state.T/state.P/deltaZ*(
         jplus5(state[yi]*state.P/state.T,dstate['{}z0'.format(yi)]*Pz0/Tz0,yAzL*PzL/TzL)*vplus5 -
         jminus5(state[yi]*state.P/state.T,dstate['{}z0'.format(yi)]*Pz0/Tz0,yAzL*PzL/TzL)*vminus5)
        -dg.psi*state.T/state.P*dstate['d{}'.format(param.xnames(param.ynames.index(yi)))] #eg. yname == yA (param.yis == [yA,yB], param.xnames=[xA,xB]) -> "dstate['d{}'.format(param.xnames(param.ynames.index(yi)))]" gives  "dstate[dxA]"
        -state[yi]/state.P*dP
        +state[yi]/state.T*dT) #yAzl need to be changed  

  #dyB=(1/dg.Pe*state.T/state.P/deltaZ*(
  #  jplus5(state.P/state.T,Pz0/Tz0,PzL/TzL)*diffplus(state.yB,yBzL)/deltaZ-
  #  jminus5(state.P/state.T,Pz0/Tz0,PzL/TzL)*diffminus(state.yB,yBz0)/deltaZ)
  #     -state.T/state.P/deltaZ*(
  #       jplus5(state.yB*state.P/state.T,yBz0*Pz0/Tz0,yBzL*PzL/TzL)*vplus5 -
  #       jminus5(state.yB*state.P/state.T,yBz0*Pz0/Tz0,yBzL*PzL/TzL)*vminus5)
  #     -dg.psi*state.T/state.P*dxB
  #     -state.yB/state.P*dP
  #     +state.yB/state.T*dT)

  dxi=column_stack(dstate[xi] for xi in param.xnames)
  dyi=column_stack(dstate[yi] for yi in param.ynames)
  
  return dP,dT,dTw,dxi,dyi

def rate_of_change(state, in_P=None, out_P=None,
                   in_orifice=0, out_orifice=0,
                   in_yA=None,
                   in_yB=None,
                   out_yA=None,
                   out_yB=None,verbose=False,
                   prod_orifice=0, modehalf=None):
  #calculate in and out velocities and call rate_of_change_vel
  #orifice_flow2(p1, p2, d, C=0.70, T=T_room):
  P_feed=in_P*norm.P0       #convert to Pa abs
  P_1=state.P[0]*norm.P0    #denormalize, convert to Pa
  if verbose:
    sys.stdout.write('state.P: {}\n'.format(state.P))
  P_N=state.P[-1]*norm.P0
  P_exit=out_P*norm.P0
  #compute flows in sN/m3
  if verbose:
    print('P_feed {} P_1 {} P_N {} P_exit {}'.format(P_feed,P_1,P_N,P_exit))
  in_flow=orifice.orifice_flow2(P_feed/1000,P_1/1000,in_orifice,T=param.Ta)/60
  out_flow=orifice.orifice_flow2(P_N/1000,P_exit/1000,out_orifice,T=param.Ta)/60
  #NOTE:
  #HERE WE override the flow rates in m/3
  #if modehalf is not None and modehalf[0]==1:
  #  #we are doing pressurization, override the inflow rate with a constant
  #  in_flow=9.5/60/1000    # LPM to in Nm^3/sec
  #these are flow rates in Nm^3/second
  #now convert to velocity in m/s
  with warnings.catch_warnings():
    warnings.filterwarnings('error')
    try:
      in_vel=in_flow*(1e5/P_1)/param.area    # [m/s]
      out_vel=out_flow*(1e5/P_N)/param.area  # [m/s]
    except Warning:
      sys.stdout.write('WARNING: {} {}\n'.format(in_flow, P_1))
      sys.stdout.write('WARNING: {} {}\n'.format(out_flow, P_N))
      raise 
  #convert prod_flow to dimensionless flow rates
  if verbose:
    print('in_vel={}'.format(in_vel))
    print('out_vel={}'.format(out_vel))
    print('Pressure: {}'.format(state.P))
  #convert to dimensionless
  in_vel/=norm.v0
  out_vel/=norm.v0
  #constant_vel=True
  #if constant_vel and in_vel>0:
  #  in_vel=1.0
  #  out_vel=0
  jplus5,jminus5=difference.gen_j5_functions(param.mode,in_vel>=0 and out_vel>=0)
  #print('flow rates, m/s {} {}'.format(in_vel, out_vel))
    
  return rate_of_change_vel(state, vin=in_vel, vout=out_vel,
                     in_yA=in_yA, in_yB=in_yB,
                            out_yA=out_yA,out_yB=out_yB,verbose=verbose,
                            jplus5=jplus5,jminus5=jminus5)

def rate_of_change_product(prod,state,prod_orifice=0,verbose=False):
  #state is container supplying the product
  P_N=state.P[-1]*norm.P0
  P_tank=prod.Pprod*norm.P0
  #compute flows in sN/m3
  in_prod_flow=orifice.orifice_flow2(P_N/1000,P_tank/1000,prod_orifice,T=param.Ta)/60
  #output from product tank
  out_prod_flow=orifice.orifice_flow2(P_tank/1000,1e5/1000,param.product_orifice,T=param.Ta)/60
  #these are flow rates in Nm^3/second
  #convert prod_flow to dimensionless flow rates
  in_prod_flow/=norm.z0**2*norm.v0
  out_prod_flow/=norm.z0**2*norm.v0
  net_flow=in_prod_flow-out_prod_flow
  #Product tank
  dPprod=dg.kprod*net_flow
  #if dg.kprod > 1e6 or prod.Pprod > 1e6 or state.yA[-1] - prod.yprod) > 1e6 or in_prod_flow>1e6:
  with warnings.catch_warnings():
    warnings.filterwarnings('error')
    try:
      dyprod=dg.kprod/prod.Pprod*(state.yA[-1] - prod.yprod)*in_prod_flow
    except Warning:
      sys.stdout.write('WARNING:{} {} {} {}\n'.format(prod.Pprod, state.yA, prod.yprod, in_prod_flow))
      raise 
  #print('in={} out={}'.format(in_prod_flow,out_prod_flow))
  #print(prod.Pprod,prod.yprod,net_flow,dg.kprod)
  #print('d',dPprod,dyprod)
  return dPprod, dyprod

def find_mode(t,state1,state2,prod):
  #find out what mode we are in
  #returns mode (1,2,3), half (0,1)
  half=0
  while t > param.cycle_time*2:
    t-=2*param.cycle_time
  if t > param.cycle_time:
    t-=param.cycle_time
    half=1
  if t < param.vent_time:
    if half<1:
      if state1.P[-1] > prod.Pprod:
        mode=2
      else:
        mode=1
    else:
      if state2.P[-1] > prod.Pprod:
        mode=2
      else:
        mode=1
  else:
    mode=3    # blowdown
  return mode, half

def connections(stateA,stateB,prod,mode,half):
  #determine the connections for containers
  #This will be called twice, with states flipped and half reversed
  #This determines that connections for stateA, with stateB being the other one
  #mode=1,2,3
  #half=0,1
  #prod is 'Pprod': array(1,), yprod: array(1,)
  assert mode in (1,2,3)
  assert half in (0,1)
  P_feed=param.feed_pressure*1e5/norm.P0
  P_atm=1e5/norm.P0
  if not half and mode==1:
    #pressurize
    in_P=P_feed
    out_P=P_atm
    in_orifice=param.input_orifice
    out_orifice=0  # no output
    in_yA=param.feed_O2
    in_yB=param.feed_N2
    #these don't matter because there is no backwards flow
    out_yA=.5
    out_yB=.5
  elif not half and mode==2:
    #production
    in_P=P_feed
    #output pressure is the product tank pressure
    out_P=prod.Pprod.item()
    in_orifice=param.input_orifice
    out_orifice=param.output_orifice  #feeding into product tank
    in_yA=param.feed_O2
    in_yB=param.feed_N2
    #these don't matter because there is no backwards flow
    out_yA=.5
    out_yB=.5
  elif not half and mode==3:
    #vent
    #produce
    in_P=P_feed
    #the other state:
    out_P=stateB.P[-1]
    in_orifice=param.input_orifice
    out_orifice=param.vent_orifice  #Note: do we have both orifices open?
    #NOTE: might want prod orifice > 0 here?
    in_yA=param.feed_O2
    in_yB=param.feed_N2
    #these don't matter because there is no backwards flow (hopefully)
    out_yA=1.0  # to see if we have a bug
    out_yB=0
  elif half and (mode==1 or mode==2):
    #blowdown
    in_P=P_atm
    out_P=P_atm
    in_orifice=param.blowdown_orifice  # reverse flow
    out_orifice=0  # no output
    #these don't matter because there is reverse flow
    in_yA=param.feed_O2
    in_yB=param.feed_N2
    #don't matter because end flow is off
    out_yA=.5
    out_yB=.5
  elif half and mode==3:
    #blowdown with venting
    in_P=P_atm
    out_P=stateB.P[-1]
    in_orifice=param.blowdown_orifice  # reverse flow
    out_orifice=param.vent_orifice  # no output
    #these don't matter because there is reverse flow
    #NOTE: need to check equations for reverse flow (fix boundary conditions)
    in_yA=param.feed_O2
    in_yB=param.feed_N2
    #output components are the other tank output
    out_yA=stateB.yA[-1]
    out_yB=stateB.yB[-1]
  else:
    raise 'error in mode'
  ret=AttrDict()
  ret.in_P=in_P
  ret.out_P=out_P
  ret.in_orifice=in_orifice
  ret.out_orifice=out_orifice
  ret.in_yA=in_yA
  ret.in_yB=in_yB
  ret.out_yA=out_yA
  ret.out_yB=out_yB
  return ret

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

def myfun(t, x):
  #function used by the ODE solver
  #must return dx/dt
  global count
  count+=1
  N=param.N
  state1,state2,prod=x_to_state(x)
  #fixup some negative values that cause us a problem - ODE solver tries to
  #explore disallowed regions
  fix=True
  if gverbose:
    sys.stdout.write('1 yA {}\n'.format(state1.yA))
    sys.stdout.write('1 yB {}\n'.format(state1.yB))
    sys.stdout.write('1 P {}\n'.format(state1.P))
    sys.stdout.write('2 yA {}\n'.format(state2.yA))
    sys.stdout.write('2 yB {}\n'.format(state2.yB))
    sys.stdout.write('2 P {}\n'.format(state2.P))
  #NOTE: since we have 1/P in our equations, we can't let P be zero, which
  #The ODE solver seems to want to explore
  if fix:
    state1.yA=util.limit(state1.yA,0,1)
    state1.yB=util.limit(state1.yB,0,1)
    state1.P=util.limit(state1.P,0.01,2)
    state1.T=util.limit(state1.T,0.01,1.5)
    state1.Tw=util.limit(state1.Tw,0.01,1.5)
    state2.yA=util.limit(state2.yA,0,1)
    state2.yB=util.limit(state2.yB,0,1)
    state2.P=util.limit(state2.P,0.01,2)
    state2.T=util.limit(state2.T,0.01,1.5)
    state2.Tw=util.limit(state2.Tw,0.01,1.5)    
  #Figure out what state we are in
  #Now calculate the input and output pressures and orifices based on the mode
  mode, half=find_mode(t,state1,state2,prod)
  if count % 5000 == 0:
    sys.stdout.write('***{:10.2f} myfun {:7.4f} half {} mode {} min {:6.2f} max {:6.2f}\n'.format(
      time.time()-tstart,t,half,mode,min(x),max(x)))
  #two pressures we will use
  P_feed=param.feed_pressure*1e5/norm.P0
  P_atm=1e5/norm.P0
  P_product=param.product_pressure*1e5/norm.P0
  conn1=connections(state1,state2,prod,mode,half)
  conn2=connections(state2,state1,prod,mode,1-half)
  prod_orifice=0
  supply_tank=state1
  if mode==2:
    prod_orifice=param.output_orifice
    if half==0:
      supply_tank=state1
    else:
      supply_tank=state2
  rate1=rate_of_change(state1,**conn1, verbose=gverbose, modehalf=(mode,half))
  rate2=rate_of_change(state2,**conn2, modehalf=(mode,1-half))
  prod_rate=rate_of_change_product(prod,supply_tank,prod_orifice,verbose=gverbose)
  ret=np.hstack(rate1+rate2+prod_rate)
  return ret

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

def product_flow_rate(P):
  #given an array of product pressures, compute the corresponding output flow rates
  #in LPM
  p1=P*norm.P0
  p2=1e5   # 1 atm in [Pa]
  #Nm3/min
  f=orifice.orifice_flow2(p1/1000,p2/1000,param.product_orifice)
  return f*1000   # LPM
  
def print_state_ranges(data):
  for n in param.state_names:
    print('{:6}: {:8.2f} - {:8.2f}'.format(n,np.min(data[n]),np.max(data[n])))

####
#This is the main function of this module, should be the only one exposed
#to the outside
gverbose=False
def simulate_and_plot(override_param, do_plots=False, pause=False,
             outdir=None,params_file='params',verbose=False, roi=None):
  data,ret,param,pickle_name,out_place= simulate(override_param, outdir=outdir,
                                                 params_file=params_file, verbose=verbose)
  if do_plots:
    plots.plot(data, roi, pause=pause, out_place=out_place)
  return ret,param

def plot_data(pickle_name,pause=False,basename=None,roi=None, cycles=None):
  #read a pickle file created by simulate and create plots.  Basename is directory
  #plus optional characters to put in front of plot names
  #cycles is (from,to) tuple of the cycle numbers to plot, inclusive, or a single int
  with open(pickle_name,'rb') as fp:
    data=pickle.load(fp)
    plots.plot(data,roi,pause=pause,out_place=basename, cycles=None)
  
def simulate(override_param, outdir=None,params_file='params',verbose=False):
  #params to override the defaults from params.py
  #plots: create output plots
  #pause: if creating plots, pause when they are displayed on screen for <Enter>
  #returns metrics of simulation as a dict, and the full params used
  #writes a log file (and optional plots) into outdir.  use outdir=None for
  #no log file and instead direct all outout to stdout.  plots will go into
  #current directory in that case.
  #params is the file to import (default is params)
  #roi is (t1,t2) time range for the spatial plots
  #returns out_place: directory+logcode prefix that it used for log file
  global gverbose
  gverbose=verbose
  params=importlib.import_module(params_file)

  np.seterr(all='warn')
  global param
  global norm
  global dg
  global count
  global tstart
  #we will alter print statement for this module to go into the log file
  global print
  if outdir is None:
    dbg=mylog.Log(outdir=None,descr='psa')
    OUTDIR='./'
  else:
    OUTDIR=outdir  # for log files and png files
    dbg=mylog.Log(outdir=OUTDIR,descr='psa')
  #change the print function in the module to print to the log file
  old_print=print
  print=dbg.print

  param=params.create_param(override_param, dbg.print)
  norm=create_norm()
  
  print('simulation params: time 0 to {}, step {}'.format(param.end_time,param.tstep))
  print('step time might be increased if long cycle time')
  print('normalization velocity {:8.3f} m/s'.format(param.norm_v0))
  print('normalization of time norm.t0 {:8.3f} sec'.format(norm.t0))
  print('starting with params import file: {}'.format(params))
  state1=init(param,norm)
  state2=init(param,norm)
  prod=AttrDict()
  prod.Pprod=1e5/norm.P0          # 1 atm to start
  prod.yprod=param.feed_yi        # fract of O2 in product tank : eg.[0.2, 0.8]

  #create dimensionless groups for the constants needed in the simulation
  dg=create_dgroups(norm,param,state1)

  np.set_printoptions(precision=5, edgeitems=6, linewidth=95, suppress=False)
  pp=pprint.PrettyPrinter(indent=2)

  print('git revision:{}'.format(util.get_git_commit()))
  print('log code: {}'.format(dbg.log_code))
  sys.stdout.write('log file prefix code: {}\n'.format(dbg.log_code))

  print('version:{}'.format(vinfo))
  print('plt {}'.format(matplotlib.__version__))
  #print out the parameters and numbers for the record
  print('dimensionless groups:\n',pp.pformat(dg.__dict__))
  print('parameters:\n',pp.pformat(param.__dict__))

  #create directory and log_code to prepend to out png files
  util.safe_mkdir(OUTDIR)
  out_place=os.path.join(OUTDIR,dbg.log_code)
  #write to screen, not log file
  sys.stdout.write('out_place={}\n'.format(out_place))
  #output is to 1 ATM pressure (abs)
  x01=np.hstack((state1.P,state1.T,state1.Tw,state1.xA,state1.xB,state1.yA,state1.yB))
  x02=np.hstack((state2.P,state2.T,state2.Tw,state2.xA,state2.xB,state2.yA,state2.yB))
  x03=np.array([prod.Pprod, prod.yprod])
  print(x01.shape,x02.shape,x03.shape)
  x0=np.hstack((x01,x02,x03))
  tstart=time.time()
  #time normalized by norm.t0
  #Vectorized means that myfun can take a matrix of values to evaluate dx/dt at,
  #not that x0 is a matrix
  compute_alpha(verbose=True)  # just print out the values for the log file
  #we add a small number because we want to get full cycle at end
  #also, increase tstep is a long cycle time:
  param.tstep=max(param.tstep, param.cycle_time/100)
  end_simulation=param.end_time+1e-6
  t_span=(0,end_simulation)
  #we want to include the end_simulation time
  ev=np.arange(0,end_simulation,param.tstep)
  print('t_span {}'.format(t_span))
  #print('t_eval {}'.format(ev))
  sys.stdout.write('Calling ODE solver\n')
  bunch=scipy.integrate.solve_ivp(myfun, t_span, x0, vectorized=False, t_eval=ev,
                                  method='BDF')
###                                  method='Radau')
  if not bunch.success:
    print(bunch.message)
    sys.stdout.write(bunch.message)
    return None
  print('evaluated at {} places'.format(bunch.t.size))
  print('num eval of RHS {}'.format(bunch.nfev))
  print('num eval of Jacobian {}'.format(bunch.njev))
  print('total number of elements in data is {}'.format(bunch.y.size))
  print('number of elements > 2 = {}'.format(np.sum(bunch.y>2)))
  print('number of elements > 1 = {}'.format(np.sum(bunch.y>1)))
  print('number of elements < 0 = {}'.format(np.sum(bunch.y<0)))
  
  #t will be array of time points (t_points)
  #y will be array of (n x t_points) values
  N=param.N
  #take the output data from ODE solver and place in data.<name> for plotting
  data=AttrDict()
  data.data1,data.data2,data.datap=data_to_state(bunch.y)
  data.t=bunch.t
  print('bunch.t',data.t)
  data.param=param
  print_state_ranges(data.data1)
  print_state_ranges(data.data2)
  #compute xAs at each point
  data.data1.xAs,data.data1.xBs=adsorp_equilibrium(data.data1.yA,data.data1.yB,data.data1.P)
  data.data2.xAs,data.data2.xBs=adsorp_equilibrium(data.data2.yA,data.data2.yB,data.data2.P)
  #compute velocity at each point
  data.data1.V=compute_vel(data.data1)
  data.data2.V=compute_vel(data.data2)
  #Output key statistics in log file with keywords
  #get indexes for full cycle
  print('CYCLE {:8.2f}'.format(param.cycle_time))
  print('VENT {:8.2f}'.format(param.vent_time))
  #indexes to select over - the last 2 cycles for container 1
  r=np.argmax(bunch.t>(max(bunch.t)-2*param.cycle_time))
  ret=AttrDict()
  ret.time=datetime.datetime.now().replace(microsecond=0).isoformat()
  ret.end_time=end_simulation
  ret.cycle_time=param.cycle_time
  ret.vent_time=param.vent_time
  ret.norm_t0=norm.t0
  ret.real_cycle_time=param.cycle_time*norm.t0
  ret.real_vent_time=param.vent_time*norm.t0
  #these are all min,max,mean for the last cycle of container 1
  #for container 1, we calculate parameters for all the odd cycles and store in
  # a list.  To get last one, use index [-1]
  idx=0
  ret.container_y=[]
  ret.prod_pressure=[]
  ret.prod_y=[]
  ret.prod_flow_rate=[]
  while True:
    start=param.cycle_time*idx
    stop=param.cycle_time*(idx+1)
    if stop > max(bunch.t):
      break
    #create a mask of the times for this odd cycle (starting with 1st one)
    mask=np.logical_and(start <= bunch.t, bunch.t <= stop)
    ret.container_y.append([np.min(data.data1.yA[-1,mask]),
                            np.max(data.data1.yA[-1,mask]), 
                            np.mean(data.data1.yA[-1,mask])])
    ret.prod_pressure.append([np.min(data.datap.Pprod[0,mask]),
                              np.max(data.datap.Pprod[0,mask]),
                              np.mean(data.datap.Pprod[0,mask])])
    ret.prod_y.append([np.min(data.datap.yprod[0,mask]),
                       np.max(data.datap.yprod[0,mask]),
                       np.mean(data.datap.yprod[0,mask])])
    #compute product flow in LPM
    ret.prod_flow_rate.append([product_flow_rate(np.min(data.datap.Pprod[0,mask])),
                               product_flow_rate(np.max(data.datap.Pprod[0,mask])),
                               product_flow_rate(np.mean(data.datap.Pprod[0,mask]))])
    idx+=2
    
  ret.logcode=dbg.log_code
  print('results:\n',pp.pformat(ret))
  print('These are repr() outputs for programmatic reading')
  print('PARAM {}'.format(repr(param)))
  print('RESULT {}'.format(repr(ret)))
  data.ret=ret
  #save to pickle file
  pickle_name='{}-data.pickle'.format(out_place)
  with open(pickle_name, 'wb') as fp:
    pickle.dump(data,fp)
  #close the log file because might be called again before exiting and will
  #want a new log file for that simulation
  dbg.close()
  #restore print function
  print=old_print
  return data, ret, param, pickle_name, out_place

def simulate_singlebed(params_file='params',output_dir='./outcome'):

  paramsmodule=importlib.import_module(params_file)
  #CORE FUNCTION
  #START
  global param
  global dg

  param = paramsmodule.create_param(None)  # set parameters
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

def cyclic_steady_state(params_file='params',output_dir='./outcome', maxcycle = 100, tolerance = 1e-4):
  # this function simulate mutiple psa process until cyclic steady state is reached, and store the data generated during simulation
  # return: "status" variable has the following data structure, type(status): AttriDict 
  #          
  #          status.param: "params used in simulation" [type: AttriDict ] e.g. status.param.nocomponents = 3 (corresponds to what you've set in the params_file)
  #          status.snap: "snapshot at each cycle" [type: list] e.g. status.snap[]
  # 
  
  paramsmodule=importlib.import_module(params_file)
  #CORE FUNCTION
  #START
  global param
  global dg
  status =AttrDict()
  status.snap = AttrDict()
  param = paramsmodule.create_param(None)  # set parameters
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
  if output_dir is None:
    dbg=mylog.Log(outdir=None,descr='psa')
    OUTDIR='./'
  else:
    OUTDIR=output_dir  # for log files and png files
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

  plots.plot_singlebed(data)
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
    assert (param.index_tend == len(bunch.t)-1)
    t_eva=bunch.t[param.index_tblw+1:param.index_tend]
    eva_Pz0 = 1/param.norm_P0 *(param.PL + (param.PI-param.PL)*np.exp(-param.lamda_eva* (t_eva-param.norm_tblw) *param.norm_t0))  # t(=tao in the reference equation (57)) is already normalized

    eva_Tz0 = data.T[0, param.index_tblw+1:param.index_tend] # dummy boundary conditions, engergy balances are not actually computed
    eva_vz0 = -2/param.deltaZ*(4/150*param.epst**2)*(param.rp**2)*param.norm_P0/param.mu/param.norm_v0/param.L*(data.P[0, param.index_tblw+1:param.index_tend] - eva_Pz0) #inlet veloctiy computed with P1 and P0.5
    assert((eva_vz0 <= 0).all()) #check if all the velocity is negative, if so apply abs to eva_vz0, because we need a postive number for molar amount
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

def save_to_pickle():
  pass
  
  


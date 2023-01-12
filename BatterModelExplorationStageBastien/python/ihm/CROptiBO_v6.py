# Auteur : Bastien VENOT
# Date   : 06/09/2022
# Projet : CROpti

# Script appelé par l'IHM pour lancer le calcul du contrôle de charge de la batterie
# et pour la récupération des résultats

from casadi import *
import numpy as np
import math
import time

from CROptiCellConfig import *

##Classe des résultats d'optimisation
class Results():
    def __init__(self,tf,ib,soc,vb,Tb,Nsamp):
        self.tf = tf
        self.ib = ib
        self.soc= soc
        self.vb = vb
        self.Tb = Tb
        self.t  = np.linspace(0,tf,Nsamp+1)


## Définition contenant le modèle d'optim
# Appelé par l'IHM pour réaliser la simulation
def run_opti(cell,sim):
    
    # Intialisation
    Imin = cell.Imin
    Imax = cell.Imax
    Vbmin = cell.Vmin
    Vbmax = cell.Vmax
    Qb = cell.Qb
    L = cell.L
    R = cell.R
    m = cell.m
    Tbmax = cell.Tbmax
    Tbmin = cell.Tbmin
    Cb = cell.Cp
    h = cell.h
    
    S  = 2*np.pi*L*R
        
    soc0 = sim.soc0
    Tb0 = sim.Tb0
    Text = sim.Text
    N = sim.Nsamp
    
    R0 = 0.012
    Rs = 0.016
    Cs = 786.2714
    Rl = 0.0353
    Cl = 6299.5676
    
    a1 = -0.5863
    a2 = 21.9
    a3 = 3.414
    a4 = 0.1102
    a5 = -0.1718
    a6 = 0.008
    
    socf = 0.999

    opti = Opti() # Optimization problem
        
    # ---- decision variables ---------
    X   = opti.variable(4,N+1)  # state trajectory
    soc = X[0,:]                # state of charge
    vs  = X[1,:]                # short-time transient voltage response
    vl  = X[2,:]                # long-time transient voltage response
    Tb  = X[3,:]                # battery temperature
    U = opti.variable(1,N)      # control trajectory (charging current)
    tf = opti.variable()        # final time
    
    
    # ---- objective          ---------
    lag = 0
    may = tf
    J = may + lag
    opti.minimize(J)
    
    # ---- dynamic constraints --------
    # length of a control interval
    dt = tf/N
    
    # state functions : dx/dt = f(x,u)
    f = lambda x,u: vertcat(-u/Qb,\
                            u/Cs-x[1]/(Rs*Cs),\
                            u/Cl-x[2]/(Rl*Cl),\
                            (R0*u**2 + (x[1]**2)/Rs + (x[2]**2)/Rl - h*S*(x[3]-Text) ) / (m*Cb)) 
    
    # loop over control intervals
    for k in range(N):
        #R0  = b1*(X[0,k]**4) + b2*(X[0,k]**3) + b3*(X[0,k]**2) +b4*X[0,k] + b5
        #Rs  = k1*np.exp(-k2*X[0,k]) + k3
        #Cs  = d1*(X[0,k]**4) + d2*(X[0,k]**3) + d3*(X[0,k]**2) + d4*X[0,k] + d5
        
        # Runge-Kutta 4 integration
        kr1 = f(X[:,k],         U[:,k])
        kr2 = f(X[:,k]+dt/2*kr1, U[:,k])
        kr3 = f(X[:,k]+dt/2*kr2, U[:,k])
        kr4 = f(X[:,k]+dt*kr3,   U[:,k])
        x_next = X[:,k] + dt/6*(kr1+2*kr2+2*kr3+kr4)
        opti.subject_to(X[:,k+1]==x_next) # close the gaps
        
        
    voc = a1*np.exp(-a2*soc) + a3 + a4*soc + a5*np.exp(-a6/(1-soc))
    vb = voc[0:N] - vs[0:N] - vl[0:N] - R0*U
    
    
    # ---- states constraints -----------
    opti.subject_to(opti.bounded(0,soc,1))
    opti.subject_to(opti.bounded(Tbmin,Tb,Tbmax))
    opti.subject_to(opti.bounded(-Imax,U,-Imin))
    opti.subject_to(opti.bounded(Vbmin,vb,Vbmax))
    
    # ---- boundary conditions --------
    opti.subject_to(soc[0]==soc0)
    opti.subject_to(Tb[0]==Tb0)
    opti.subject_to(vs[0]==0)
    opti.subject_to(vl[0]==0)
    opti.subject_to(soc[-1]==socf)
    opti.subject_to(U[-1]==0)
    
    
    # ---- misc. constraints  ----------
    opti.subject_to(tf>=0) # Time must be positive
    
    # ---- initial values for solver ---
    opti.set_initial(tf, 20*60)
    
    # ---- solve NLP              ------
    opti.solver("ipopt") # set numerical backend
    
    tic = time.process_time()
    sol = opti.solve()   # actual solve
    toc = time.process_time()
    
    solving_time = toc-tic
    print("solving_time : "+str(solving_time)+" s")
    print("Temperature : "+str(sol.value(Tb[-1]-273.15))+"°C")
    print("Imax : "+str(cell.Imax)+"A")

    print("socf = "+str(sol.value(soc[-1])))
    
    # ---- post-processing        ------
    from pylab import step, figure, show, spy
    import matplotlib.pyplot as plt
    
    t = np.linspace(0,sol.value(tf),N+1)
    
    ib_res = np.zeros(N+1)    
    ib_res[0:N] = sol.value(U)
    ib_res[N] = None
    
    vb_res = np.zeros(N+1)
    vb_res[0:N] = sol.value(vb)
    vb_res[N] = None
    
    soc_res= sol.value(soc)
    Tb_res = sol.value(Tb)
    
    # put optimization results in a class
    res = Results(sol.value(tf),
                  sol.value(U),
                  sol.value(soc),
                  sol.value(vb),
                  sol.value(Tb),
                  N)
    
    # return results class to IHM script
    return res
# Auteur : Bastien VENOT
# Date   : 06/09/2022
# Projet : CROpti

# Programme d'analyse du contrôle optimal du courant de charge d'une batterie
# selon la configuration de la Fonction Objectif

import numpy  as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import datetime
import time

from casadi import *
from pylab  import step, figure, show, spy

from CROptiCellConfig import *

#### Initialisation ####

opti = Opti() # Optimization problem
N = 100 # number of control intervals
nb_tests = 4

Text = 25 + 273.15      # External temperature [K]
soc_min = 0             # Minimal state of charge
soc_max = 0.999         # Maximal state of charge

time_in_solver  = []
total_comput_time = []

soc_table     = np.zeros((N+1,nb_tests), dtype=object)
vb_table      = np.zeros((N+1,nb_tests), dtype=object)
Tb_table      = np.zeros((N+1,nb_tests), dtype=object)
current_table = np.zeros((N+1,nb_tests), dtype=object)
C_rate_table = np.zeros((N+1,nb_tests), dtype=object)

tf_vec      = []  
tf_minute_vec = []  
Tbf_vec     = []
Tb0_vec     = []
DeltaTb_vec = []
socf_vec    = []
Imax_vec    = []
Imean_vec   = []

current_dic = {}
I_dic       = {}
soc_dic     = {}
vb_dic      = {}
Tb_dic      = {}

#### Loop over the test ####

for j in range(nb_tests):
    tstart_comput = time.process_time()

    X   = opti.variable(4,N+1)  # state trajectory
    soc = X[0,:]                # state of charge
    vs  = X[1,:]                # short-term transient voltage response [V]
    vl  = X[2,:]                # long-term transient voltage response [V]
    Tb  = X[3,:]                # Battery temperature [K]
    U = opti.variable(1,N)      # control trajectory (charging current [A])
    tf = opti.variable()        # final time

    soc0 = soc_min
    socf = soc_max
    Tb0  = Text 

    ## ---- objective functions  ---------    
    # lag = Lagrgange part
    # may = Mayer part
    
    lag = 0
    
    if j == 0:
        opti_name = "FO0 : MAY{Tb} + LAG{err2_soc}"
        may = X[3,-1] 
        for l in range(N+1): 
            lag += (X[0,l]-1)**2
    
    elif j == 1:
        opti_name = "FO1 : LAG{Tb + err2_soc}"
        may = 0 
        for l in range(N+1): 
            lag += (X[0,l]-1)**2 + X[3,l] 
    
    elif j == 2:
        opti_name = "FO2 : MAY{Tb + tf}"
        may = X[3,-1] + tf
        lag = 0
    
    elif j == 3:
        opti_name = "FO3 : MAY{tf} + LAG{Tb}"
        may = tf
        for l in range(N+1): 
            lag += X[3,l]
    
    J = may + lag
    opti.minimize(J)
    
    ##### dynamic constraints #### 
    # length of a control interval
    dt = tf/N
    
    # state functions : dx/dt = f(x,u)
    f = lambda x,u: vertcat(-u/Qb,\
                            u/Cs-x[1]/(Rs*Cs),\
                            u/Cl-x[2]/(Rl*Cl),\
                            (R0*u**2 + (x[1]**2)/Rs + (x[2]**2)/Rl - h*S*(x[3]-Text) ) / (m*Cb)) 
    
    # loop over control intervals
    for k in range(N): 
       # Runge-Kutta 4 integration
       k1 = f(X[:,k],         U[:,k])
       k2 = f(X[:,k]+dt/2*k1, U[:,k])
       k3 = f(X[:,k]+dt/2*k2, U[:,k])
       k4 = f(X[:,k]+dt*k3,   U[:,k])
       x_next = X[:,k] + dt/6*(k1+2*k2+2*k3+k4)
       opti.subject_to(X[:,k+1]==x_next) # close the gaps
          
    voc = a1*np.exp(-a2*soc) + a3 + a4*soc + a5*np.exp(-a6/(1-soc))
    vb = voc[0:N] - vs[0:N] - vl[0:N] - R0*U

    #### Constraints ####
    # ---- states constraints -----------
    opti.subject_to(opti.bounded(soc_min,soc,soc_max))
    opti.subject_to(opti.bounded(Tbmin,Tb,Tbmax))
    opti.subject_to(opti.bounded(-Ifast,U,-Imin))
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
    opti.set_initial(U, -Ifast)
    opti.set_initial(tf, 30*60)

    # ---- solve NLP              ------
    opti.solver("ipopt")                 # set numerical backend
    tstart_solver = time.process_time()
    sol = opti.solve()                   # actual solve
    tstop_solver = time.process_time()
            
    
    #### post-processing ####
    
    t = np.linspace(0,sol.value(tf),N+1)
    
    soc_table[:,j]       = sol.value(soc)
    Tb_table[:,j]        = sol.value(Tb-273.15)
    current_table[0:N,j] = abs(sol.value(U))
    C_rate_table[0:N,j]  = abs(sol.value(U))/Qb * 3600
    vb_table[0:N,j]      = sol.value(vb)

    tf_vec.append(sol.value(tf))    
    tf_minute_vec.append(sol.value(tf)/60)   
    Tb0_vec.append(Tb0)
    Tbf_vec.append(sol.value(Tb[-1]-273.15))
    DeltaTb_vec.append(sol.value(Tb[-1]) - sol.value(Tb[0]))
    socf_vec.append(sol.value(soc[-1]))
    Imax_vec.append(max(abs(sol.value(U))))
    Imean_vec.append(np.mean(abs(sol.value(U))))
    
    current_table[N,j]  = None
    C_rate_table[0:N,j] = None
    vb_table[N,j]       = None
    
    time_in_solver.append(tstop_solver - tstart_solver)    
    total_comput_time.append(tstop_solver - tstart_comput)

    # profiles = {"time (s)":t,\
    #             "time (min)":t/60,\
    #             "current (A)":current_table[:,j],\
    #             "current (C-rate)":C_rate_table[:,j],\
    #             "Vbat (V)":vb_table[:,j],\
    #             "SOC":soc_table[:,j],\
    #             "Temperature (°C)":Tb_table[:,j]}
        
    # globals()["dfprof_"+str(j)] = pd.DataFrame(profiles)
    
    res_key = opti_name
    t_key = "t(FO"+str(j)+")"
    
    I_dic[t_key] = vb_dic[t_key] = soc_dic[t_key] = Tb_dic[t_key] = t/60

    I_dic[res_key]   = current_table[:,j]
    soc_dic[res_key] = soc_table[:,j]
    vb_dic[res_key]  = vb_table[:,j]
    Tb_dic[res_key]  = Tb_table[:,j]
    
    # ---- plot results        ------
    
    # figure()
    # plt.plot(t/60,sol.value(voc))
    # plt.xlabel('Temps [min]')
    # plt.ylabel('Voc')
    # plt.grid()
    
    # figure()
    # plt.plot(t/60,sol.value(vs))
    # plt.xlabel('Temps [min]')
    # plt.ylabel('Vs')
    # plt.title('Vs state [-]')
    # plt.grid()
    
    # figure()
    # plt.plot(t/60,sol.value(vl))
    # plt.xlabel('Temps [min]')
    # plt.ylabel('Vl')
    # plt.title('Vl state [-]')
    # plt.grid()
    
    figure()
    plt.plot(t/60,sol.value(Tb-273.15))
    plt.xlabel('Temps [min]')
    plt.ylabel('Tb')
    plt.title('Battery temperature vs time')
    plt.grid()
    
    # figure()
    # plt.plot(t/60,sol.value(soc))
    # plt.xlabel('Temps [min]')
    # plt.ylabel('soc')
    # plt.title('SoC state [-]')
    # plt.grid()
    
    # figure()
    # plt.plot(t/60,abs(sol.value(U)),'k')
    # plt.xlabel('Temps [min]')
    # plt.ylabel('Current')
    # plt.title('Input Signal [A]')
    # plt.grid()
    
    
    fig, ax1 = plt.subplots() 
    plt.grid()
    
    ax1.set_xlabel('Time [min]') 
    ax1.set_ylabel('Current & SoC [C-rate]')
    #ax1.set_ylim([0,1])
    ax1.plot(t[0:N]/60,abs(sol.value(U))/Qb*3600, color = 'red') 
    ax1.plot(t/60,sol.value(soc), color = 'green') 
    ax1.tick_params(axis ='y', labelcolor = 'red') 
    
    ax2 = ax1.twinx()   
    ax2.set_ylabel('Tension [V]', color = 'blue') 
    #ax2.set_ylim([0,4.4])
    ax2.plot(t[0:N]/60,sol.value(vb), color = 'blue') 
    ax2.tick_params(axis ='y', labelcolor = 'blue')    
    
    # figure()
    # spy(sol.value(jacobian(opti.g,opti.x)))
    # figure()
    # spy(sol.value(hessian(opti.f+dot(opti.lam_g,opti.g),opti.x)[0]))
    


## ----- Export to Excel ------

now = datetime.datetime.now()
datetime_str = now.strftime("%Y%m%d_%H%M%S")  
            
filename='CROpti_BO_ObjFunc_Tests_'+datetime_str+'.xlsx'

data = {"Run #"          : np.arange(nb_tests),\
        "Tb0 (°C)"       : Tb0_vec,\
        "Tbf (°C)"       : Tbf_vec,\
        "Delta Tbat (°C)": DeltaTb_vec,\
        "tf (sec)"       : tf_vec,\
        "tf (min)"       : tf_minute_vec,\
        "soc"            : socf_vec,\
        "Imax (A)"       : Imax_vec,\
        "Imean (A)"      : Imean_vec,\
        "solver (sec)"   : time_in_solver,\
        "loop (sec)"     : total_comput_time}
dfval = pd.DataFrame(data)

dfI   = pd.DataFrame(I_dic)
dfsoc = pd.DataFrame(soc_dic)
dfvb  = pd.DataFrame(vb_dic)
dfTb  = pd.DataFrame(Tb_dic)

with pd.ExcelWriter(filename) as writer:
    # for i in range(nb_tests):  
    #     sheetname='Run #'+str(i)         
    #     globals()["dfprof_"+str(i)].to_excel(writer,sheet_name=sheetname,index=False)
     
    dfval.to_excel(writer,sheet_name='analysis',index=False)
    dfI.to_excel(writer,sheet_name='current (A)',index=False)
    dfsoc.to_excel(writer,sheet_name='SOC',index=False)
    dfvb.to_excel(writer,sheet_name='Vbat (V)',index=False)
    dfTb.to_excel(writer,sheet_name='Tbat (°C)',index=False)


plt.show()

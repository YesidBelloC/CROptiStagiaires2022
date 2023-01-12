# Auteur : Bastien VENOT
# Date   : 06/09/2022
# Projet : CROpti

# Programme d'analyse du contrôle optimal du courant de charge d'une batterie
# selon la température intiale de la batterie (analyse paramétrique avec variation
# de la température initiale)

from casadi import *
import numpy as np
import pandas as pd
import math
import time
import datetime

from pylab import step, figure, show, spy
import matplotlib.pyplot as plt

from CROptiCellConfig import *

##### Initialisation ####

opti = Opti()   # Optimization problem
N = 100         # number of control intervals

soc_min = 0     # Minimal state of charge
soc_max = 0.999 # Maximal state of charge

Text = 25 + 273.15      # External temperature [K]
Tb0_min  = 57 + 273.15  # Minimal battery temperature during the test
Tb0_max  = 59 + 273.15  # Maximal battery temperature during the test
Tb0_step = 1            # temperature step for the loop

# Temperature range for the test
rangeTb0 = np.arange(Tb0_min, Tb0_max + Tb0_step, Tb0_step)
nb_tests  = len(rangeTb0)

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

j = 0
for Tb0 in rangeTb0:
    tstart_comput = time.process_time()
    
    ## ---- decision variables ---------
    X   = opti.variable(4,N+1)  # state trajectory
    soc = X[0,:]                # state of charge
    vs  = X[1,:]                # short-term transient voltage response [V]
    vl  = X[2,:]                # long-term transient voltage response [V]
    Tb  = X[3,:]                # Battery Temperature [K]
    U = opti.variable(1,N)      # control trajectory (charging current [A])
    tf = opti.variable()        # final time
    
    soc0 = soc_min
    socf = soc_max
        
    ## ---- objective function ---------
    # lag = Lagrgange part
    # may = Mayer part
    
    lag = 0
    may = tf
    J = may + lag
    opti.minimize(J)
    
    ## ---- dynamic constraints --------
    
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
        kr1 = f(X[:,k],          U[:,k])
        kr2 = f(X[:,k]+dt/2*kr1, U[:,k])
        kr3 = f(X[:,k]+dt/2*kr2, U[:,k])
        kr4 = f(X[:,k]+dt*kr3,   U[:,k])
        x_next = X[:,k] + dt/6*(kr1+2*kr2+2*kr3+kr4)
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
    opti.subject_to(tf>=0)              # Time must be positive
    
    # ---- initial values for solver ---
    opti.set_initial(tf, 40*60)

    # ---- solve NLP              ------
    opti.solver("ipopt")                 # set numerical backend
    tstart_solver = time.process_time()
    sol = opti.solve()                   # actual solve
    tstop_solver = time.process_time()
    
    print("********** Test conclu : "+str(Tb0-273.15)+" °C ****************")
    
    ## ---- post-processing        ------
    t = np.linspace(0,sol.value(tf),N+1)
    
    soc_table[:,j]       = sol.value(soc)
    Tb_table[:,j]        = sol.value(Tb-273.15)
    current_table[0:N,j] = abs(sol.value(U))
    C_rate_table[0:N,j]  = abs(sol.value(U))/Qb * 3600
    vb_table[0:N,j]      = sol.value(vb)

    tf_vec.append(sol.value(tf))    
    tf_minute_vec.append(sol.value(tf)/60)    
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
    
    res_key = "Tb0="+str(Tb0-273.15)+"°C"
    t_key = "t(Tb0="+str(Tb0-273.15)+") [min]"
    
    I_dic[t_key] = vb_dic[t_key] = soc_dic[t_key] = Tb_dic[t_key] = t/60

    I_dic[res_key]   = current_table[:,j]
    soc_dic[res_key] = soc_table[:,j]
    vb_dic[res_key]  = vb_table[:,j]
    Tb_dic[res_key]  = Tb_table[:,j]

    j += 1
    
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
    
    #figure()
    #plt.plot(t/60,sol.value(Tb-273.15))
    #plt.xlabel('Temps [min]')
    #plt.ylabel('Tb')
    #plt.title('Battery temperature vs time')
    #plt.grid()
    
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
    #ax1.plot(t[0:N]/60,abs(sol.value(U))/Qb*3600, color = 'red', label='Current') 
    #ax1.plot(t/60,sol.value(soc), color = 'green', label='SoC') 
    ax1.plot(t[0:N]/60,abs(sol.value(U))/Qb*3600, color = 'red') 
    ax1.plot(t/60,sol.value(soc), color = 'green')
    ax1.tick_params(axis ='y', labelcolor = 'red') 
    
    ax2 = ax1.twinx()   
    ax2.set_ylabel('Tension [V]', color = 'blue') 
    #ax2.set_ylim([0,4.4])
    ax2.plot(t[0:N]/60,sol.value(vb), color = 'blue') 
    ax2.tick_params(axis ='y', labelcolor = 'blue')
    #plt.legend() 
    
    
    
    # figure()
    # spy(sol.value(jacobian(opti.g,opti.x)))
    # figure()
    # spy(sol.value(hessian(opti.f+dot(opti.lam_g,opti.g),opti.x)[0]))
    
    #show()


## ---- Export to Excel        ------

now = datetime.datetime.now()
datetime_str = now.strftime("%Y%m%d_%H%M%S")              
filename='CROpti_BO_Temperature_Analysis'+datetime_str+'.xlsx'

data = {"Run #"          : np.arange(nb_tests),\
        "Tb0 (°C)"       : rangeTb0-273.15,\
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
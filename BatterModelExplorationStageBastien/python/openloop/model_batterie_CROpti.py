# Auteur : Bastien VENOT
# Date   : 06/09/2022
# Projet : CROpti

# Programme de calcul du contrôle optimal du courant de charge d'une batterie
# avec le framework CASADI, en considérant les variables d'état de charge, de 
# tension et de température.

from casadi import *
import numpy as np
import math
import time

from CROptiCellConfig import *

##### Initialisation ####

opti = Opti()

N = 100                 #number of control intervals

Text  = 25 + 273.15     # External temperature [K]
Tb0  = Text             # Initial battery temperature [K]
soc0 = 0                # Initial state of charge
socf = 0.999            # Final state of charge


##### decision variables ####

X   = opti.variable(4,N+1)  # State vector
soc = X[0,:]                # State of charge
vs  = X[1,:]                # short-term transient voltage response [V]
vl  = X[2,:]                # long-term transient voltage response [V]
Tb  = X[3,:]                # battery temperature [K]
U = opti.variable(1,N)      # control variable (charging current [A])
tf = opti.variable()        # final time (charge time)
   
#### Objective Function #####

#Lagrange Integration
lag = 0
#for l in range(N): 
    #lag += a*(X[0,l]-1)**2 + (1-a)*X[3,l]/Tbmax
    #lag += (X[0,l]-1)**2
    #lag += (1-a)*X[3,l]/Tbmax

# Mayer part
may = tf
#may = (soc[-1]-1)**2
#may = (soc[-1]-1)**2 + 1/60 *Tb[-1]
#may = (Tb[-1]-25)**2
#may= 1/60 * Tb[-1]

J = may + lag
opti.minimize(J)

#### dynamic constraints ####

# length of a control interval
dt = tf/N

# State functions : dx/dt = f(x,u)
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


#### Constraints ####
# ---- states constraints -----------
opti.subject_to(opti.bounded(0,soc,1))
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
opti.subject_to(tf>=0)

# ---- initial values for solver ---
opti.set_initial(tf, 30*60)

##### solve NLP ####
opti.solver("ipopt")

tic = time.process_time()
sol = opti.solve()
toc = time.process_time()

solving_time = toc-tic
print("solving_time : "+str(solving_time)+" s")
print("Temperature : "+str(sol.value(Tb[-1]-273.15))+"°C")

print("socf = "+str(sol.value(soc[-1])))


##### post-processing ####

from pylab import step, figure, show, spy
import matplotlib.pyplot as plt

t = np.linspace(0,sol.value(tf),N+1)

# figure()
# plt.plot(t/60,sol.value(voc))
# plt.xlabel('Temps [min]')
# plt.ylabel('Voc')
# plt.grid()

figure()
plt.plot(t/60,sol.value(vs))
plt.xlabel('Temps [min]')
plt.ylabel('Vs')
plt.title('Vs state [-]')
plt.grid()

figure()
plt.plot(t/60,sol.value(vl))
plt.xlabel('Temps [min]')
plt.ylabel('Vl')
plt.title('Vl state [-]')
plt.grid()


figure()
plt.plot(t/60,sol.value(Tb-273.15))
plt.xlabel('Time [min]')
plt.ylabel('Tb [°C]')
plt.title('Battery temperature vs time')
plt.grid()

# figure()
# plt.plot(t/60,sol.value(soc))
# plt.xlabel('Temps [min]')
# plt.ylabel('soc')
# plt.title('SoC state [-]')
# plt.grid()

# figure()
# plt.plot(np.linspace(0,sol.value(T),N)/60,abs(sol.value(U)),'k')
# plt.xlabel('Temps [min]')
# plt.ylabel('Current')
# plt.title('Input Signal [A]')
# plt.grid()


fig, ax1 = plt.subplots() 
plt.grid()

ax1.set_xlabel('Time [min]') 
ax1.set_ylabel('Current (red) & SoC (green) [C-rate]')
ax1.set_ylim([0,4])
ax1.plot(t[0:N]/60,abs(sol.value(U))/Qb*3600, color = 'red') 
ax1.plot(t/60,sol.value(soc), color = 'green') 
ax1.tick_params(axis ='y', labelcolor = 'black') 

ax2 = ax1.twinx()   
ax2.set_ylabel('Tension [V]', color = 'blue') 
ax2.set_ylim([2,4])
ax2.plot(t[0:N]/60,sol.value(vb), color = 'blue') 
ax2.tick_params(axis ='y', labelcolor = 'blue')
#plt.legend() 
plt.show()


# figure()
# spy(sol.value(jacobian(opti.g,opti.x)))
# figure()
# spy(sol.value(hessian(opti.f+dot(opti.lam_g,opti.g),opti.x)[0]))

show()

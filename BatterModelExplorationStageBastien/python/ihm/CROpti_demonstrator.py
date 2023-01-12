# Auteur : Bastien VENOT
# Date   : 06/09/2022
# Projet : CROpti

# Script de l'interface graphique permettant une utilisation plus intuitive du
# simulateur de charge de batterie

from tkinter import *
import runpy
import numpy as np

import matplotlib.pyplot as plt
from   matplotlib.figure import Figure
from   matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

import CROptiCellConfig as config
import CROptiBO_v6 as openloop


## Functions

def layout_adv_param(active_const):
    # Fonction permettant une disposition dynamique de la fenetre de configuration
    # avancée, selon le choix de paramètres constant ou variable du modèle électrique

    global frame_const, frame_var
    ElecM = config.ElecModel()

    in_elec_parameters = []
    in_poly_parameters = []

    constants = ["R0","Rs","Cs","Rl","Cl"]

    factors   = [["b1","b2","b3","b4","b5"],
                 ["k1","k2","k3"],
                 ["d1","d2","d3","d4","d5"],
                 ["g1","g2","g3","g4","g5"],
                 ["h1","h2","h3","h4","h5","h6"]]

    if active_const==True:
        state_const,state_poly = 'normal', 'disabled'
    else:
        state_const,state_poly = 'disabled', 'normal'

    i=0
    for param in constants:
        Label(frame_const,text=param).grid(row=0,column=i, padx=5)

        in_elec_param_i = Entry(frame_const,state=state_const)
        in_elec_param_i.grid(row=1,column=i,padx=10,pady=5)
        in_elec_parameters.append(in_elec_param_i)

        in_elec_param_i.insert(0,str(ElecM.list_param[i]))

        i+=1

    k=0
    for list in factors:
        j=0
        in_poly_parameters.append([])

        column_title = constants[k]+"(soc)"
        Label(frame_var,text=column_title).grid(row=0,column=2*k+1, pady=5)

        for factor in list:
            Label(frame_var,text=factor).grid(row=j+1,column=2*k, padx=5)
            in_poly_param_j = Entry(frame_var,state=state_poly)
            in_poly_param_j.grid(row=j+1,column=2*k+1, padx=10,pady=5)
            in_poly_parameters[k].append(in_poly_param_j)
            j+=1
        k+=1

    frame_const.grid(row=1,column=0, columnspan=2)
    frame_var.grid(row=2,column=0, columnspan=2)


def init_advanced_param(win_param,active_const):
    # Fonction d'initialisation de la fenêtre de paramètres avancés
    global frame_const, frame_var
    frame_const = Frame(win_param)
    frame_var = Frame(win_param)
    layout_adv_param(active_const)

def modify_win_param(active_const):
    # Lance la modification de la dispotion des éléments de al fenêtre de paramètres avancés
    global frame_const, frame_var
    frame_const.grid_forget()
    frame_var.grid_forget()
    layout_adv_param(active_const)

def open_advanced_param():
    # Ouverture de la fenêtre de paramètres avancés et disposition des éléments fixes
    win_param = Toplevel(app)
    win_param.resizable(False,False)
    win_param.title('Advanced Setting')
    win_param.grab_set()

    active_const = BooleanVar(None,True)

    rb_const = Radiobutton(win_param,
                           text="Constant",
                           variable=active_const,
                           value=True,
                           command = lambda:modify_win_param(True),
                           indicatoron=False,
                           width=30)

    rb_var = Radiobutton(win_param,
                         text="Variable",
                         variable=active_const,
                         value=False,
                         command = lambda:modify_win_param(False),
                         indicatoron=False,
                         width=30)

    rb_const.grid(row=0,column=0,sticky='E')
    rb_var.grid(row=0,column=1,sticky='W')

    init_advanced_param(win_param,True)

    frame_equations = Frame(win_param)
    frame_equations.grid(row=3,column=0,columnspan=2,pady=10)

    R0_equation = "R0(soc) = b1*soc^4 + b2*soc^3 + b3*soc^2 +b4*soc + b5"
    Rs_equation = "Rs(soc) = k1*exp(-k2*soc) + k3"
    Cs_equation = "Cs(soc) = d1*soc^4 + d2*soc^3 + d3*soc^2 + d4*soc + d5"
    Rl_equation = "Rl(soc) = g1*exp(-g2*soc) + g3 + g4*soc"
    Cl_equation = "Cl(soc) = h1*soc^5 + h2*soc^4 + h3*soc^3 + h4*soc^2 + h5*soc + h6"

    Label(frame_equations,text=R0_equation).pack(anchor=W)
    Label(frame_equations,text=Rs_equation).pack(anchor=W)
    Label(frame_equations,text=Cs_equation).pack(anchor=W)
    Label(frame_equations,text=Rl_equation).pack(anchor=W)
    Label(frame_equations,text=Cl_equation).pack(anchor=W)

    win_param.mainloop()

def init_setting():
    # Initialisation des valeurs affichées dans les entrées de l'IHM
    cell = config.Cell()
    sim  = config.SimulationParam()

    in_Imin.delete(0,last=END)
    in_Imax.delete(0,last=END)
    in_Vmin.delete(0,last=END)
    in_Vmax.delete(0,last=END)
    in_Qb.delete(0,last=END)
    in_L.delete(0,last=END)
    in_R.delete(0,last=END)
    in_m.delete(0,last=END)

    in_Tbmin.delete(0,last=END)
    in_Tbmax.delete(0,last=END)
    in_Cp.delete(0,last=END)
    in_h.delete(0,last=END)

    in_soc0.delete(0,last=END)
    in_Tb0.delete(0,last=END)
    in_Text.delete(0,last=END)
    in_Nsamp.delete(0,last=END)

    in_Imin.insert(0,str(cell.Imin))
    in_Imax.insert(0,str(cell.Imax))
    in_Vmin.insert(0,str(cell.Vbmin))
    in_Vmax.insert(0,str(cell.Vbmax))
    in_Qb.insert(0,str(round(cell.Qb/3600,2)))
    in_L.insert(0,str(cell.L))
    in_R.insert(0,str(cell.R))
    in_m.insert(0,str(cell.m))

    in_Tbmin.insert(0,str(cell.Tbmin-273.15))
    in_Tbmax.insert(0,str(cell.Tbmax-273.15))
    in_Cp.insert(0,str(cell.Cb))
    in_h.insert(0,str(cell.h))

    in_soc0.insert(0,str(sim.soc0))
    in_Tb0.insert(0,str(sim.Tb0))
    in_Text.insert(0,str(sim.Text))
    in_Nsamp.insert(0,str(sim.Nsamp))


def restore_default_setting():
    # Fonction permettant de revenir aux paramètres par défaut
    init_setting()

def get_entries():
    # Récupère les valeurs indiquées dans les entrées et retourne les valeurs
    cell = config.Cell()
    sim = config.SimulationParam()

    cell.Imin = eval(in_Imin.get())
    cell.Imax = eval(in_Imax.get())
    cell.Vmin = eval(in_Vmin.get())
    cell.Vmax = eval(in_Vmax.get())
    cell.Qb = eval(in_Qb.get())*3600
    cell.L = eval(in_L.get())
    cell.R = eval(in_R.get())
    cell.m = eval(in_m.get())

    cell.Tbmin = eval(in_Tbmin.get())+273.15
    cell.Tbmax = eval(in_Tbmax.get())+273.15
    cell.Cp = eval(in_Cp.get())
    cell.h = eval(in_h.get())

    sim.soc0 = eval(in_soc0.get())/100
    sim.Tb0 = eval(in_Tb0.get())+273.15
    sim.Text = eval(in_Text.get())+273.15
    sim.Nsamp = eval(in_Nsamp.get())

    return cell,sim


def plot_optres(selgraph):
    # Affichage des graphiques selon l'onglet sélectionné dans la zone de graphiques
    global cangraph

    fig = plt.figure(figsize=(5,3.5))
    cangraph = FigureCanvasTkAgg(fig, master=frame_graph)
    cangraph.draw()
    cangraph.get_tk_widget().grid(row=0, column=0, sticky='WESN')


    if selgraph == 0:
        fig.add_subplot(111).plot(optres.t/60, optres.soc*100,'black')
        plt.xlabel('Time [min]')
        plt.ylabel('SOC [%]')
        plt.title('Sate of Charge (SOC) vs time')
    elif selgraph == 1:
        fig.add_subplot(111).plot(optres.t[0:sim.Nsamp]/60, abs(optres.ib),'red')
        plt.xlabel('Time [min]')
        plt.ylabel('Ib [A]')
        plt.title('Charge current (Ib) vs time')
    elif selgraph == 2:
        fig.add_subplot(111).plot(optres.t[0:sim.Nsamp]/60, optres.vb,'blue')
        plt.xlabel('Time [min]')
        plt.ylabel('Vb [V]')
        plt.title('Battery voltage (Vb) vs time')
    elif selgraph == 3:
        fig.add_subplot(111).plot(optres.t/60, optres.Tb-273.15,'orange')
        plt.xlabel('Time [min]')
        plt.ylabel('Tb [A]')
        plt.title('External surface temperature of the battery (Tb) vs time')

    plt.grid()

    toolbarFrame = Frame(frame_graph)
    toolbarFrame.grid(row=1,column=0,sticky='W')
    toolbar = NavigationToolbar2Tk(cangraph, toolbarFrame)
    toolbar.update()


def run_sim():
    # Lance le programme de calcul d'optimisation
    global optres,cell,sim
    cell,sim = get_entries()
    optres   = openloop.run_opti(cell,sim)
    plot_optres(0)




############# MAIN ###############

app = Tk()
app.resizable(False,False)
app.title("CROpti Demonstrator - © EXPLEO Group")

optres = openloop.Results(0,0,0,0,0,0)

## Création et disposition des sections de la fenêtre

frame1 = LabelFrame(app, text= 'CELL CONFIGURATION')
frame2 = LabelFrame(app,text='THERMAL PARAMETERS')
frame3 = Frame(app)
frame_plot = Frame(app,highlightbackground="grey", highlightthickness=1)
frame_output = Frame(app, highlightbackground="grey", highlightthickness=1)

frame3.grid(row=0,column=0,padx=20,pady=10,sticky='E')
frame1.grid(row=1,column=0,padx=20,pady=10,sticky='WESN')
frame2.grid(row=2,column=0,padx=20,pady=10,sticky='WESN')
frame_plot.grid(row=0,column=1,rowspan=2,padx=20,pady=10,sticky='WESN')
frame_output.grid(row=2,column=1,padx=20,pady=10,sticky='WESN')


frame1_L = Frame(frame1)
frame1_R = Frame(frame1)
frame1_L.grid(row=0,column=0,padx=20,pady=10,sticky='WESN')
frame1_R.grid(row=0,column=1,padx=(20,30),pady=10,sticky='WESN')

frame3_L = Frame(frame3)
frame3_R = LabelFrame(frame3,text='SIMULATION CONTROL')
frame3_L.grid(row=0,column=0,padx=20)
frame3_R.grid(row=0,column=1,pady=10,ipadx=20)

frame_selgraph  = Frame(frame_plot)
frame_graph  = Frame(frame_plot)
frame_selgraph.grid(row=0,column=0,sticky='WESN')
frame_graph.grid(row=1,column=0,sticky='WESN')



## Entrées des paramètres principaux

pady_entry = 10
padx_label = 5

Label(frame1_L, text="Imin [A]").grid(row=0,column=0,padx=padx_label,sticky='W')
Label(frame1_L, text="Imax [A]").grid(row=1,column=0,padx=padx_label,sticky='W')
Label(frame1_L, text="Vmin [V]").grid(row=2,column=0,padx=padx_label,sticky='W')
Label(frame1_L, text="Vmax [V]").grid(row=3,column=0,padx=padx_label,sticky='W')

Label(frame1_R, text="Qb [A.h]").grid(row=0,column=0,padx=padx_label,sticky='W')
Label(frame1_R, text="L [m]").grid(row=1,column=0,padx=padx_label,sticky='W')
Label(frame1_R, text="R [m]").grid(row=2,column=0,padx=padx_label,sticky='W')
Label(frame1_R, text="m [kg]").grid(row=3,column=0,padx=padx_label,sticky='W')

in_Imin = Entry(frame1_L)
in_Imax = Entry(frame1_L)
in_Vmin = Entry(frame1_L)
in_Vmax = Entry(frame1_L)

in_Qb = Entry(frame1_R)
in_L  = Entry(frame1_R)
in_R  = Entry(frame1_R)
in_m  = Entry(frame1_R)

in_Imin.grid(row=0,column=1,pady=pady_entry)
in_Imax.grid(row=1,column=1,pady=pady_entry)
in_Vmin.grid(row=2,column=1,pady=pady_entry)
in_Vmax.grid(row=3,column=1,pady=pady_entry)

in_Qb.grid(row=0,column=1,pady=pady_entry)
in_L.grid(row=1,column=1,pady=pady_entry)
in_R.grid(row=2,column=1,pady=pady_entry)
in_m.grid(row=3,column=1,pady=pady_entry)


## Entrées des paramètres thermiques
Label(frame2, text="Minimum Operating Temperature [°C]")\
     .grid(row=0,column=0,padx=padx_label,sticky='W')

Label(frame2, text="Maximum Operating Temperature [°C]")\
     .grid(row=1,column=0,padx=padx_label,sticky='W')

Label(frame2, text="Specific Heat Capacity [J/(kg.K)]")\
     .grid(row=2,column=0,padx=padx_label,sticky='W')

Label(frame2, text="Convective Heat Transfert Coeff. [W/(m².K)]")\
     .grid(row=3,column=0,padx=padx_label,sticky='W')

in_Tbmin = Entry(frame2)
in_Tbmax = Entry(frame2)
in_Cp = Entry(frame2)
in_h  = Entry(frame2)

in_Tbmin.grid(row=0,column=1,pady=pady_entry)
in_Tbmax.grid(row=1,column=1,pady=pady_entry)
in_Cp.grid(row=2,column=1,pady=pady_entry)
in_h.grid(row=3,column=1,pady=pady_entry)



## Boutons d'actions

Button(frame3_L,text='RUN',command=run_sim)\
 .grid(row=0,column=0,pady=10,sticky='WE')

Button(frame3_L,text='Advanced Param.',command=open_advanced_param)\
 .grid(row=1,column=0,pady=10,sticky='WE')

Button(frame3_L,text='Default Setting',command=restore_default_setting)\
 .grid(row=2,column=0,pady=10,sticky='WE')


## Entrées des paramètres de simulation

Label(frame3_R, text="Initial State of Charge [%]")\
     .grid(row=0,column=0,padx=(0,padx_label),sticky='W')

Label(frame3_R, text="Initial Temperature [°C]")\
     .grid(row=1,column=0,padx=(0,padx_label),sticky='W')

Label(frame3_R, text="Ambient Temperature [°C]")\
     .grid(row=2,column=0,padx=(0,padx_label),sticky='W')

Label(frame3_R, text="Number of samples")\
     .grid(row=3,column=0,padx=(0,padx_label),sticky='W')

in_soc0  = Entry(frame3_R)
in_Tb0   = Entry(frame3_R)
in_Text  = Entry(frame3_R)
in_Nsamp = Entry(frame3_R)

in_soc0.grid(row=0,column=1,pady=pady_entry)
in_Tb0.grid(row=1,column=1,pady=pady_entry)
in_Text.grid(row=2,column=1,pady=pady_entry)
in_Nsamp.grid(row=3,column=1,pady=pady_entry)



## Disposition du panel des graphs
rb_graph = IntVar(None,0)

rb_soc = Radiobutton(frame_selgraph,
                     text='SOC',
                     variable=rb_graph,
                     value=0,
                     indicatoron=False,
                     command=lambda:plot_optres(0))

rb_I   = Radiobutton(frame_selgraph,
                     text='Current',
                     variable=rb_graph,
                     value=1,
                     indicatoron=False,
                     command=lambda:plot_optres(1))

rb_V   = Radiobutton(frame_selgraph,
                     text='Voltage',
                     variable=rb_graph,
                     value=2,
                     indicatoron=False,
                     command=lambda:plot_optres(2))

rb_Tb  = Radiobutton(frame_selgraph,
                     text='Temperature',
                     variable=rb_graph,
                     value=3,
                     indicatoron=False,
                     command=lambda:plot_optres(3))

rb_soc.grid(row=0,column=0,sticky='WE')
rb_I.grid(row=0,column=1,sticky='WE')
rb_V.grid(row=0,column=2,sticky='WE')
rb_Tb.grid(row=0,column=3,sticky='WE')

fig = plt.figure(figsize=(5,3.5))
cangraph = FigureCanvasTkAgg(fig, master=frame_graph)
cangraph.draw()
cangraph.get_tk_widget().grid(row=0, column=0, sticky='WESN')


#can_graph = Canvas(frame_graph, width=550, height=380, bg='white')\
#           .grid(row=1,column=0,rowspan=2,columnspan=4)

## Sorties (Résultats)
Label(frame_output, text="Test = 0").grid(row=0)


init_setting()

app.mainloop()


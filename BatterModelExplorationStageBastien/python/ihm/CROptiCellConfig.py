# Auteur : Bastien VENOT
# Date   : 06/09/2022
# Projet : CROpti

# Script d'initialisation des paramètres de la batterie pour l'affichage intiale
# de l'interface graphique

# A123 Systems' APRI8650mI LiFeP04 cells

# Coefficients des polynomes d'estimation des paramètres variables du modèle électrique
class PolynomeVoc:
    a1 = -0.5863
    a2 = 21.9
    a3 = 3.414
    a4 = 0.1102
    a5 = -0.1718
    a6 = 0.008

class PolynomeR0:
    b2 = -0.2518
    b1 = 0.1369
    b3 = 0.1609
    b4 = -0.041
    b5 = 0.0821

class PolynomeRs:
    k1 = 5.896*(10**-10)
    k2 = -18.75
    k3 = 0.01388

class PolynomeCs:
    d1 = -10260
    d2 = 17230
    d3 = -10130
    d4 = 2340
    d5 = 684.9

class PolynomeRl:
    g1 = 8.913*(10**-15)
    g2 = -32.23
    g3 = 0.031
    g4 = 0.007473

class PolynomeCl:
    h1 = -154100
    h2 = 204200
    h3 = -4009
    h4 = -81240
    h5 = 22830
    h6 = 7144

# Paramètres de la batterie
class Cell:
    Qb    = 1.1 * 3600  # Capacité de la cellule [A.s]
    Vbmin = 2           # Tension min
    Vbmax = 3.6         # Tension max
    Vbnom = 3.3         # Tensino nominale
    Imax = 4            # Courant maximal de charge rapide
    Imin  = 0           # Courant minimal de charge
    Tbmax = 60 +273.15  # Température maximale de la cellule
    Tbmin = -30+273.15  # Température minimale de la cellule

    m  = 0.039          # Masse de la cellule (kg)
    L  = 0.065          # Hauteur de la cellule (m)
    R  = 0.009          # Rayon de la cellule (m)

    Cb = 1720           # Specific Heat for LFPO4 cell
    h  = 0.2            # convective heat transfert coefficient

# Paramètres constants du modèle électrique
class ElecModel:
    #R0 = 0.079
    R0 = 0.012
    Rs = 0.016
    Cs = 786.2714
    Rl = 0.0353
    Cl = 6299.5676

    list_param = [R0,Rs,Cs,Rl,Cl]

# Paramètres de simulation
class SimulationParam:
    Tb0 = 25
    soc0 = 0
    Text = 25
    Nsamp = 100





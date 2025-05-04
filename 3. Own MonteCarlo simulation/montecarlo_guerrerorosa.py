# -*- coding: utf-8 -*-
"""
Created on Sun May  4 13:11:04 2025

@author: andre
"""

import numpy as np
import matplotlib.pyplot as plt

# VARIABLES TERMODINAMIQUES FIXADES-----------------------------------------------

N=1000 #nombre de mol.lecules
V=1 #volum [m^3]
T=293 #temperatura [K]
kb=1.4E-23 #constant de Boltzmann [J/K]
k=1 #constant del potencial harmonic E = 1/2*k*r^2
L=1 #llargada de la recta / costat del pla / aresta de la caixa on son els atoms

n_p = 100 #nombre de passos de la simulacio

def E(r): #potencial harmonic
    return 0.5 * k * np.sum(r**2)

def Cv(d):
    return d*N*k/2 #calculem la capacitat teorica
    
#DEMANAR LA DIMENSIO--------------------------------------------------------------
print("Indiqueu la dimensio en que voleu el sistema (1, 2 o 3)")
while True:
    d = int(input())
    if d==1 or d==1 or d==3:
        print("La dimensio escollida es d =", d)
        break
    else:
        print("Si us plau introdueix unicament 1, 2 o 3")




#BUCLE QUE FA LA SIMULACIO---------------------------------------------------
r = np.random.uniform(-L, L, size=(N, d)) #llista on s'emmagatzemmen les posicions de les N particules en d dimensions contingudes en la longitud L

for pas in range(n_p):
    



    

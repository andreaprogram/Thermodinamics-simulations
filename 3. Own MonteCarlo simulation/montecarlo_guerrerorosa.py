# -*- coding: utf-8 -*-
"""
Created on Sun May  4 13:11:04 2025

@author: andrea i david
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
delta=L/2 #desplacament

n_p = 100 #nombre de passos de la simulacio

def m(r):  # modul al quadrat d’un vector posicio r
    return np.sum(r**2)

def E(R):  # energia total del sistema
    E_total = 0
    for i in range(len(R)): #R es el vector de vectors posicions r de les particules
        E_total += 0.5 * k * m(R[i]) #R[i] es el vector posicio de la particula i
    return E_total

def Cv(d):
    return d*N*kb/2 #calculem la capacitat teorica
    
#DEMANAR LA DIMENSIO--------------------------------------------------------------
print("Indiqueu la dimensio en que voleu el sistema (1, 2 o 3)")
while True:
    d = int(input())
    if d==1 or d==2 or d==3:
        print("La dimensio escollida es d =", d)
        break
    else:
        print("Si us plau introdueix unicament 1, 2 o 3")




#BUCLE QUE FA LA SIMULACIO---------------------------------------------------
r = np.random.uniform(-L, L, size=(N, d)) #llista on s'emmagatzemmen les posicions de les N particules en d dimensions contingudes en la longitud L
energies=[] #llista on s'emmagatzemmen les energies del sistema

for pas in range(n_p):
    i = np.random.randint(N) #numero de particula aleatori
    r_nova = r[i] + np.random.uniform(-delta, delta, size = d) #canviem la posicio de manera aleatoria
    delta_E=delta_E= 0.5 * k *(m(r_nova)-m(r[i])) #trobem la diferencia d'energia fent tan sols el canvi per la particula seleccionada (la resta s'eliminen en fer la diferencia)
     
    if delta_E<0 or np.random.rand()<np.exp(-delta_E/(kb*T)):  #REGLA DE METROPOLIS
        r[i]=r_nova
    energies.append(E(r)) #es guarda la nova energia

energies = np.array(energies)
Cv_sim=(np.mean(energies**2)-np.mean(energies)**2)/(kb*T**2) #calcul Cv simulacio

print("Capacitat calorífica simulació",Cv_sim)
print("Capacitat calorífica teorica",Cv(d))
    

    
    



    

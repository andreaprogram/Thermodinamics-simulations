# -*- coding: utf-8 -*-
"""
Created on Sun May  4 13:11:04 2025

@author: andrea i david
"""

import numpy as np
import matplotlib.pyplot as plt

# VARIABLES TERMODINAMIQUES FIXADES-----------------------------------------------

N=1000 #nombre de mol.lecules
T=293 #temperatura [K]
kb=1.4E-23 #constant de Boltzmann [J/K]
L=1 #llargada de la recta / costat del pla / aresta de la caixa on son els atoms
m=4E-3/6E23 #massa de l'Heli

n_p = 100 #nombre de passos de la simulacio

def mod(v):  # modul al quadrat d’un vector velocitat v
    return np.sum(v**2)

def E(V):  # energia total del sistema Ec=0.5*m*v**2
    E_total = 0
    for i in range(len(V)): #V es el vector de vectors velocitats v de les particules
        E_total += 0.5 * m * m(modV[i]) #V[i] es el vector velocitat de la particula i
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
v_rms=np.sqrt(d*kb*T/m) #agafem el maxim que sabemm de l'amplada de la distribucio, no podem agafar des de -infinit a +infinit a la practica
delta= v_rms/2 #variació
v = np.random.uniform(-v_rms, v_rms, size=(N, d)) #llista on s'emmagatzemmen les velocitats de les N particules en d dimensions, que segueix distribucio MB
energies=[] #llista on s'emmagatzemmen les energies del sistema

for pas in range(n_p):
    i = np.random.randint(N) #numero de particula aleatori
    v_nova = v[i] + np.random.uniform(-delta, delta, size = d) #canviem la velocitat de manera aleatoria
    delta_E=delta_E= 0.5 * m *(mod(v_nova)-mod(v[i])) #trobem la diferencia d'energia fent tan sols el canvi per la particula seleccionada (la resta s'eliminen en fer la diferencia)
     
    if delta_E<0 or np.random.rand()<np.exp(-delta_E/(kb*T)):  #REGLA DE METROPOLIS
        v[i]=v_nova
    energies.append(E(v)) #es guarda la nova energia

energies = np.array(energies)
Cv_sim=(np.mean(energies**2)-np.mean(energies)**2)/(kb*T**2) #calcul Cv simulacio

print("Capacitat calorífica simulació",Cv_sim)
print("Capacitat calorífica teorica",Cv(d))
print("Error relatiu", np.abs(Cv_sim-Cv(d))/Cv(d))
    

    
    



    

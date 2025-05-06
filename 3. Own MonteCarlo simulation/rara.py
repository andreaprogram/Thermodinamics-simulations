# -*- coding: utf-8 -*-
"""
Created on Sun May  4 13:11:04 2025

@author: andrea i david
"""

import numpy as np
import matplotlib.pyplot as plt

# VARIABLES TERMODINAMIQUES FIXADES I FUNCIONS NECESSARIES-----------------------------------------------

N=3000 #nombre de mol.lecules
T=293 #temperatura [K]
kb=1.4E-23 #constant de Boltzmann [J/K]
m=4E-3/6E23 #massa de l'Heli



n_p = 10000*N #nombre de passos de la simulacio






#DEMANAR LA DIMENSIO--------------------------------------------------------------
print("Indiqueu la dimensio en que voleu el sistema (1, 2 o 3)")
while True:
    d = int(input())
    if d==1 or d==2 or d==3:
        print("La dimensio escollida es d =", d)
        break
    else:
        print("Si us plau introdueix unicament 1, 2 o 3")



v_rms=np.sqrt(d*kb*T/m) #agafem el maxim que sabemm de l'amplada de la distribucio, no podem agafar des de -infinit a +infinit a la practica
delta= v_rms #variacio maxima en la velocitat

#BUCLE QUE FA LA SIMULACIO---------------------------------------------------
v = np.random.uniform(-v_rms, v_rms, size=(N, d)) #llista on s'emmagatzemmen les velocitats de les N particules en d dimensions, que segueix distribucio MB
energies=[] #llista on s'emmagatzemmen les energies del sistema
energies2=[]
pas=0
while pas<n_p:
    i = np.random.randint(N) #numero de particula aleatori
    v_nova = v[i] + np.random.uniform(-delta, delta, size = (1,d)) #canviem la velocitat de manera aleatoria
    delta_E= 0.5 * m *(np.sum(v_nova**2)-np.sum(v[i]**2)) #trobem la diferencia d'energia fent tan sols el canvi per la particula seleccionada (la resta s'eliminen en fer la diferencia)
     
    if delta_E<0 or np.random.uniform(0,1)<np.exp(-delta_E/(kb*T)):  #REGLA DE METROPOLIS
        v[i]=v_nova

    if pas % 10*N==0 and pas>0:
        E=m*np.sum(v**2)/2
        energies.append(E)
        energies2.append(E**2)#es guarda la nova energia
    pas+=1
    
energies = np.array(energies)
energies2 = np.array(energies2)
cv_sim=(np.mean(energies2)-np.mean(energies)**2)/(kb*T**2) #calcul Cv simulacio

print("Capacitat calorífica simulació",cv_sim)
print("Capacitat calorífica teorica",d*N*kb/2)
print("Error relatiu", np.abs(cv_sim-d*N*kb/2)/(d*N*kb/2)*100, "%")
    

    
    



    

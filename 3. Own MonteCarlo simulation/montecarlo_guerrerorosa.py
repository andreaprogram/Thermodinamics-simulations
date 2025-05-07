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
kb=1.380649E-23 #constant de Boltzmann [J/K]
m=6.646E-27 #massa de l'Heli [kg]
Na=6.022E23
n=N/Na

def E(c): #funcio que calcula l'energia (cinetica), on c es un ARRAY
    return 0.5 * m *np.sum(c**2)

def C(D): #funcio capacitat calorifica teorica
    return 0.5 * D*N*kb 



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



#Alterarem les velocitats de les particules, doncs definim els parametres:
v_rms=np.sqrt(d*kb*T/m) #agafem el maxim que sabemm de l'amplada de la distribucio, no podem agafar des de -infinit a +infinit a la practica
delta= v_rms #variacio maxima en la velocitat


#BUCLE QUE FA LA SIMULACIO---------------------------------------------------
v = np.random.uniform(-v_rms, v_rms, size=(N, d)) #llista on s'emmagatzemmen les velocitats de les N particules en d dimensions, que segueix distribucio MB
energies=[] #llista on s'emmagatzemmen les energies del sistema 
pas=0 #variable pel recompte de passos

while pas<n_p:
    i = np.random.randint(N) #numero de particula aleatori
    v_nova = v[i] + np.random.uniform(-delta, delta, size = d) #canviem la velocitat de manera aleatoria
    delta_E= E(v_nova) - E(v[i]) #trobem la diferencia d'energia fent tan sols el canvi per la particula seleccionada (la resta s'eliminen en fer la diferencia)
     
    if delta_E<=0 or np.random.uniform(0,1) < np.exp(-delta_E/(kb*T)):  #REGLA DE METROPOLIS
        v[i]=v_nova

    #emmagatzem els valors d'energia cada 1000 execucions del bucle ja que no tots els valors son d'interes, per exemple, amb aixo descartem els 1000 primers
    # i obtenim uns resultats mes acurats perque estem interessats en fer estadistica i per aixo necessitem un alt nombre de punts amb que treballar
    
    if pas % 10*N==0 and pas>0: 
        energia_total=E(v)
        energies.append(energia_total)
    pas+=1
    
energies = np.array(energies) 
Cv_sim=(np.mean(energies**2)-(np.mean(energies))**2)/(kb*T**2) #calcul Cv simulacio
cv_teo = C(d)/n #capacitat calorífica molar

print("Capacitat calorífica simulació", f"{Cv_sim/n:.2e}", "J/K*mol")
print("Capacitat calorífica teorica",f"{cv_teo:.2e}", "J/K*mol")
print("Error relatiu", f"{np.abs(Cv_sim/n-cv_teo)/(cv_teo)*100:.2f} %")
    
    




    

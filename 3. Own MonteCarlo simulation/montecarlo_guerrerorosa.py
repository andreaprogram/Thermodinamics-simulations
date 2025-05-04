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
k=1 #constant del potencial harmonic


#DEMANAR LA DIMENSIO--------------------------------------------------------------
print("Indiqueu la dimensio en que voleu el sistema (1, 2 o 3)")
while True:
    d = int(input())
    if d==1 or d==1 or d==3:
        print("La dimensio escollida es d =", d)
        break
    else:
        print("Si us plau introdueix unicament 1, 2 o 3")


C_v_teo=d*N*k/2 #calculem la capacitat teorica
def E(r):
    return 0.5 * k * np.sum(r**2)

    

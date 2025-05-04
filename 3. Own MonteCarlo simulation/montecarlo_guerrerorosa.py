# -*- coding: utf-8 -*-
"""
Created on Sun May  4 13:11:04 2025

@author: andre
"""

import numpy as np
import matplotlib.pyplot as plt

# VARIABLES TERMODINÀMIQUES FIXADES

N=1000 #nombre de mol.lecules
V=1 #volum [m^3]
T=293 #temperatura [K]

print("Indicate the dimension in which you want the system (1, 2 or 3)")
while True:
    d = int(input())
    if d==1 or d==1 or d==3:
        print("The dimension chosen is d =", d)
        break
    else:
        print("Please enter only 1, 2 or 3")

    
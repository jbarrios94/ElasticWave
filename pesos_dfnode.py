#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 14:20:09 2022

@author: cbarrios
"""

import numpy as np
from scipy.special import factorial

o = 2
A = np.zeros((o//2, o//2))

for i in range(o//2):
    for j in range(o//2):
        A[i, j] = (2*(j+1) - 1)**(2*(i+1) - 1)/factorial(2*(i+1) - 1)
        

x = np.zeros(o//2)
x[0] = 1
c = 1/2*np.linalg.solve(A, x)

print(c)
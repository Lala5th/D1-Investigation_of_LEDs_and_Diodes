#!/usr/bin/ipython3

import sys
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import scipy.constants as const
from scipy.optimize import curve_fit

plt.ion()
args = sys.argv
dtype = [('V',np.float64,()),('I',np.float64,())]
data = np.loadtxt(args[1],dtype=dtype,skiprows=0,delimiter=',')
#data = np.sort(data,order = 'V')
R = 15.08
true_V = data['V'] - R*data['I']
plt.figure()
plt.plot(true_V,data['I'],'x-')
plt.yscale('log')
plt.show()

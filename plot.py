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
data = np.sort(data,order = 'V')
T = float(args[2])
R_0 = (data['V'][-1]-data['V'][-2])/(data['I'][-1]-data['I'][-2])
R= 15.08
#true_V = data['V'] - data['I']*R
def shockley(x,C,n):
    true_V = x[0] - x[1]*R
    return true_shockley(true_V,C,n)
def true_shockley(x,C,n):
    return C*(np.exp(const.e*x/(n*T*const.k)) - 1)
params = ['C','n_id']
comb_V = np.array(list(zip(data['V'],data['I']))).transpose()
fit, cov = curve_fit(shockley,comb_V,data['I'],maxfev=10000)

xs = np.arange(data['V'][0],data['V'][-1],0.0001)
plt.plot(data['V'],data['I'],'x-')
plt.plot(data['V'],shockley(comb_V,*fit))
plt.show()
for i in zip(params,fit,np.sqrt(np.diag(cov))):
    print(i[0],':',i[1],'+-',i[2])

#!/usr/bin/ipython3

import sys
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import scipy.constants as const
from scipy.optimize import curve_fit

plt.ion()
args = sys.argv
dtype = [('id',np.int,()),('I',np.float64,()),('misc','<U20',(12,)),('V',np.float64,()),('misc2','<U20',(9,))]
data = np.loadtxt(args[1],dtype=dtype,skiprows=8,delimiter=',')
data = np.sort(data,order = 'V')
T = float(args[2])
R = float(args[3])

def shockley(x,c,Eg):
    true_V = x[0] - x[1]*R
    return (c*T**3)*np.exp(-Eg*1.6e-19/(T*const.k))*(np.exp(const.e*true_V/(T*const.k)) - 1)
def true_shockley(x,c,Eg):
    true_V = x
    return (c*T**3)*np.exp(-Eg*1.6e-19/(T*const.k))*(np.exp(const.e*true_V/(T*const.k)) - 1)
params = ['c','Eg']

fit, cov = curve_fit(shockley,np.array(list(zip(data['V'],data['I']))).transpose(),data['I'],maxfev=1000)
xs = np.arange(data['V'][0],data['V'][-1],0.0001)
plt.plot(data['V'],data['I'],'x-')
fit_shockley = true_shockley(xs,*fit)
lim = np.where(np.logical_and(xs-R*fit_shockley > data['V'].min(), data['V'].max() > xs-R*fit_shockley))
plt.plot(xs[lim]+fit_shockley[lim]*R,fit_shockley[lim])
plt.show()
for i in zip(params,fit,np.sqrt(np.diag(cov))):
    print(i[0],':',i[1],'+-',i[2])

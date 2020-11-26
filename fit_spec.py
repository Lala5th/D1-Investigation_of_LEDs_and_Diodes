#!/usr/bin/ipython3
import sys
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import scipy.constants as const
from scipy.optimize import curve_fit
from scipy.signal import find_peaks,convolve


T=0
def peak(x,A,Eg,sigma):
	global T
	E = x
	sigma = np.abs(sigma)
	pure = np.sqrt((E-Eg)*(1-(E<Eg)))*np.exp(-E*const.e/(const.k*T))
	if pure.max() != 0:
		pure = pure/pure.max()
	pure = pure*A
	dx = (E.max() - E.min())/E.size
	gx = np.arange(-5*sigma, 5*sigma, dx)
	gaussian = np.exp(-(gx/sigma)**2/2)
	conv = convolve(pure, gaussian, mode="same")
	conv = conv/conv.max()
	return A*conv
params = ['A','Eg']
paramunits = ['a.u.','eV']
dtype = [('l',np.float64,()),('I',np.float64,())]
plt.ion()
args = sys.argv
data1 = np.loadtxt(args[1],dtype=dtype,skiprows=33,comments='[',delimiter=';')
data2 = np.loadtxt(args[2],dtype=dtype,skiprows=33,comments='[',delimiter=';')
data1['l'] = const.h*const.c/(data1['l']*1e-9*const.e)
data2['l'] = const.h*const.c/(data2['l']*1e-9*const.e)
peaks1, props1 = find_peaks(data1['I'],np.max(data1['I'])*0.75,distance = 500)
peaks2, props2 = find_peaks(data2['I'],np.max(data2['I'])*0.75,distance = 500)
peak1 = data1['l'][peaks1[0]]
peak2 = data2['l'][peaks2[0]]
Eg1 = peak1 - const.k*295/(2*const.e)
Eg2 = peak2 - const.k*77/(2*const.e)
def fitLN(xmin,xmax,E,sigma=0.025):
	global T
	loc = np.where(np.logical_and(data1['l'] > xmin, data1['l'] < xmax))
	T = 77
	fit2,cov2 = curve_fit(peak,data2['l'][loc],data2['I'][loc],p0=[data2['I'].max(),E,sigma],maxfev=10000)
	plt.plot(data2['l'][loc],peak(data2['l'][loc],*fit2))
	return (fit2,cov2)

plt.figure()
plt.plot(data1['l'],data1['I'],label='RT')
plt.plot(data2['l'],data2['I'],label='LN')
plt.plot(peak1,data1['I'][peaks1[0]],'x',label='RT Fit')
plt.plot(peak2,data2['I'][peaks2[0]],'x',label='LN Fit')
plt.legend()
plt.grid()
plt.xlabel('E [eV]')
plt.ylabel('Intensity [a.u.]')
plt.show()
print(Eg1)
print(Eg2)

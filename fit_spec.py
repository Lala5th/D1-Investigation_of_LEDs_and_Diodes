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
	if sigma >1:
		raise Exception('Too large sigma')
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


def get_chi_squared(func,params,x,y):
    assert(x.shape == y.shape)
    global removed_bg
    #diff = np.abs(x[0]-x)
    #diff_nonzero = diff[np.nonzero(diff)]
    #error = 0.25*diff_nonzero[diff_nonzero.argmin()]/2
    error = np.std(data1['I'][0:200])
    chi_squared = np.sum((y - func(x,*params))**2)/error**2
    ndof = x.size - params.size
    return chi_squared/ndof

def fitLN(xmin,xmax,E,sigma=0.025):
	global T
	loc = np.where(np.logical_and(data2['l'] > xmin, data2['l'] < xmax))
	params = ['A','E','sigma']
	paramunits = ['a.u.','eV','eV']
	T = 77
	fit2,cov2 = curve_fit(peak,data2['l'][loc],data2['I'][loc],p0=[data2['I'].max(),E,sigma],maxfev=10000,sigma=1-data2['I'][loc])
	plt.plot(data2['l'][loc],peak(data2['l'][loc],*fit2),label='LN fit')
	chi_sq = get_chi_squared(peak,fit2,data1['l'],data1['I'])
	for i in zip(params,fit2,np.sqrt(np.diag(cov2)),paramunits):
		print(i[0],':',i[1],'+-',i[2],i[3])
	print('Reduced Chi-Squared:',chi_sq)
	return (fit2,cov2)

def fitLN2(xmin,xmax,E1,E2,sigma=0.025):
	global T
	fitfunc = lambda x,A1,A2,E1,E2,sigma : peak(x,A1,E1,sigma) + peak(x,A2,E2,sigma)
	params = ['A1','A2','E1','E2','sigma']
	paramunits = ['a.u.','a.u.','eV','eV','eV']
	loc = np.where(np.logical_and(data2['l'] > xmin, data2['l'] < xmax))
	T = 77
	fit2,cov2 = curve_fit(fitfunc,data2['l'][loc],data2['I'][loc],p0=[data2['I'].max(),data2['I'].max(),E1,E2,sigma],maxfev=10000,sigma=1-data2['I'][loc])
	plt.plot(data2['l'][loc],fitfunc(data2['l'][loc],*fit2),label='LN Fit')
	chi_sq = get_chi_squared(fitfunc,fit2,data2['l'],data2['I'])
	for i in zip(params,fit2,np.sqrt(np.diag(cov2)),paramunits):
		print(i[0],':',i[1],'+-',i[2],i[3])
	print('E1-E2:',np.abs(fit2[2]-fit2[3]),'+-',np.sqrt(cov2[2,2]+cov2[3,3]),'eV')
	print('Reduced Chi-Squared:',chi_sq)
	return (fit2,cov2)

def fitRT(xmin,xmax,E,sigma=0.050):
	global T
	fitfunc = lambda x,A,E,sigma : peak(x,A,E,sigma)
	params = ['A','E','sigma']
	paramunits = ['a.u.','eV','eV']
	loc = np.where(np.logical_and(data1['l'] > xmin, data1['l'] < xmax))
	T = 295
	fit2,cov2 = curve_fit(fitfunc,data1['l'][loc],data1['I'][loc],p0=[data1['I'].max(),E,sigma],maxfev=10000,sigma=1-data1['I'][loc])
	plt.plot(data1['l'][loc],fitfunc(data1['l'][loc],*fit2),label='RT Fit')
	chi_sq = get_chi_squared(fitfunc,fit2,data1['l'],data1['I'])
	fit2[1] -= const.k*295/(2*const.e)
	for i in zip(params,fit2,np.sqrt(np.diag(cov2)),paramunits):
		print(i[0],':',i[1],'+-',i[2],i[3])
	print('Reduced Chi-Squared:',chi_sq)
	return (fit2,cov2)

def fitRT2(xmin,xmax,E1,E2,sigma=0.050):
	global T
	fitfunc = lambda x,A1,A2,E1,E2,sigma : peak(x,A1,E1,sigma) + peak(x,A2,E2,sigma)
	params = ['A1','A2','E1','E2','sigma']
	paramunits = ['a.u.','a.u.','eV','eV','eV']
	loc = np.where(np.logical_and(data1['l'] > xmin, data1['l'] < xmax))
	T = 295
	fit2,cov2 = curve_fit(fitfunc,data1['l'][loc],data1['I'][loc],p0=[data1['I'].max(),data1['I'].max(),E1,E2,sigma],maxfev=10000,sigma=1-data1['I'][loc])
	plt.plot(data1['l'][loc],fitfunc(data1['l'][loc],*fit2),label='RT Fit')
	chi_sq = get_chi_squared(fitfunc,fit2,data1['l'],data1['I'])
	for i in zip(params,fit2,np.sqrt(np.diag(cov2)),paramunits):
		print(i[0],':',i[1],'+-',i[2],i[3])
	print('Reduced Chi-Squared:',chi_sq)
	return (fit2,cov2)

plt.figure()
plt.plot(data1['l'],data1['I'],label='RT')
plt.plot(data2['l'],data2['I'],label='LN')
#plt.plot(peak1,data1['I'][peaks1[0]],'x',label='RT Fit')
#plt.plot(peak2,data2['I'][peaks2[0]],'x')
plt.legend()
plt.grid()
plt.xlabel('E [eV]')
plt.ylabel('Intensity [a.u.]')
plt.show()
print(Eg1)
print(Eg2)

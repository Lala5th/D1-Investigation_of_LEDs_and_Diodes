#!/usr/bin/ipython3
import sys
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import scipy.constants as const
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

def moving_average(a, n=3) :
	ret = np.cumsum(a, dtype=float)
	ret[n:] = ret[n:] - ret[:-n]
	return ret[n - 1:] / n

def derive(x,y,n=3):
	x_avg = moving_average(x,n)
	y_avg = moving_average(np.log(y),n)
	ret = []
	for i in range(1,x_avg.size):
		dx = x_avg[i] - x_avg[i-1]
		ret.append([(x_avg[i]+x_avg[i-1])/2,(y_avg[i] - y_avg[i-1])/dx])
	ret = np.array(ret)
	ret_avg = interp1d(moving_average(ret[:,0],1),moving_average(ret[:,1],1),fill_value='extrapolate')
	mask = np.ones(ret[:,1].size, dtype=bool)
	for i in range(1,ret[:,1].size):
		if(np.abs(ret[i,1]-ret_avg(ret[i,0]-0.01))> 0.25):
			mask[i] = False
	return ret[mask]

R = 15.08
cum_data = []
I0s = []
ns = []
plt.ion()
dtype = [('V',np.float64,()),('I',np.float64,())]
IStep = 0.00001
V_err = 0.01

def monotonity_check(x,y):
	mask = np.ones(x.size, dtype=bool)
	for i in range(x.size):
		if i == 0:
			continue
		for d in range(1,i+1):
			if not mask[i-d]:
				continue
			if x[i-d] > x[i] or y[i-d] > y[i]:
				mask[i] = False
				break
		if x.size - 2 <= i:
			continue
		if y[i] > y[i+1] and y[i+1] - y[i-1] > IStep:
			mask[i] = False
			continue
	return mask

def load_data(fname):
	data = np.loadtxt(fname,dtype=dtype,skiprows=0,delimiter=',')
	mask = monotonity_check(data['V']-R*data['I'],data['I'])
	return data[mask]

def shockley(x,I0,n):
	global T
	if(n == 0):
		return 0
	ret = I0*np.exp(const.e*x/(n*T*const.k))
	return ret

def fit_shockley(V,LNI):
	fitfunc = lambda x, lnI, n: const.e*x/(n*T*const.k) + lnI
	fit,cov = curve_fit(fitfunc,V,LNI,maxfev=10000)
	fit[0] = np.exp(fit[0])
	cov[0,0] = (fit[0]*np.sqrt(cov[0,0]))**2
	return (fit,cov)

params=['I0','n']
paramunits = ['A','']

prefix = input("Add a prefix [Leave empty for no prefix]: ")
suffix = input("Add a suffix [Leave empty for no suffix]: ")

while(True):
	dataloc = input("Dataset you want to add [Leave blank if finished]: ")
	if(dataloc == ''):
		break
	try:
		data = load_data(prefix+dataloc+suffix)
	except Exception as e:
		print("Could not open file")
		print(e)
		continue
	T = float(input("Temperature of dataset [K]: "))
	true_V = data['V'] - data['I']*R
	cum_data.append([data['V'],data['I'],true_V])
	der = derive(true_V,data['I'],n=1)
	index = der[:,1].argmax()
	n = const.e/(der[index,1]*const.k*T)
	n_err = V_err*const.e/(const.k*T*der[index,1]**2)
	I0 = np.exp(np.log(data['I'][index])-der[index,1]*der[index,0])
	I0_err = V_err*np.sqrt(2)*I0*der[index,1]
	der2 = derive(der[:,0],der[:,1],n=1)
	fig, (ax1,ax2) = plt.subplots(2,1,sharex=True)
	ax1.plot(true_V,data['I'],'-')
	ax1.set_yscale('log')
	ax1.set_ylabel('I [A]')
	ax2.set_xlabel('V [V]')
	ax2.plot(der[:,0],der[:,1])
	ax2.plot(der[index,0],der[index,1],'rx')
	ax1.grid()
	ax2.grid()
	fig.canvas.show()
	fig.canvas.draw()
	fig.canvas.flush_events()
	print('n :',n,'+-',n_err)
	print('I0 :',I0,'+-',I0_err)
	I0s.append([T,I0,I0_err])
	ns.append([T,n,n_err])

def Isat(x,c,Eg):
	return c*x**3 * np.exp(-const.e*Eg/(const.k*x))
pIsat = ['c','Eg']
Iparamunits = ['A/K^3','eV']

def refit_data(dataset_i):
	dataloc = input("Dataset: ")
	if(dataloc == ''):
		return
	try:
		data = load_data(prefix+dataloc+suffix)
	except Exception as e:
		print("Could not open file")
		print(e)
		return
	T = float(input("Temperature of dataset [K]: "))
	true_V = data['V'] - data['I']*R
	cum_data[dataset_i] = [data['V'],data['I'],true_V]
	der = derive(true_V,data['I'],n=1)
	index = der[:,1].argmax()
	n = const.e/(der[index,1]*const.k*T)
	n_err = V_err*const.e/(const.k*T*der[index,1]**2)
	I0 = np.exp(np.log(data['I'][index])-der[index,1]*der[index,0])
	I0_err = V_err*np.sqrt(2)*I0*der[index,1]
	der2 = derive(der[:,0],der[:,1],n=1)
	fig, (ax1,ax2) = plt.subplots(2,1,sharex=True)
	ax1.plot(true_V,data['I'],'-')
	ax1.set_yscale('log')
	ax1.set_ylabel('I [A]')
	ax2.set_xlabel('V [V]')
	ax2.plot(der[:,0],der[:,1])
	ax2.plot(der[index,0],der[index,1],'rx')
	ax1.grid()
	ax2.grid()
	fig.canvas.show()
	fig.canvas.draw()
	fig.canvas.flush_events()
	print('n :',n,'+-',n_err)
	print('I0 :',I0,'+-',I0_err)
	I0s[dataset_i] = [T,I0,I0_err]
	ns[dataset_i] = [T,n,n_err]

def final_plots():
	global ns,I0s
	I0s = np.array(I0s)
	ns = np.array(ns)
	fit,cov = curve_fit(Isat,I0s[:,0],I0s[:,1],sigma=I0s[:,2],bounds=np.array([[0,np.inf],[0,np.inf]]).transpose(),absolute_sigma=True,maxfev = 10000)
	min = I0s[:,0].min()
	max = I0s[:,0].max()
	Is = np.linspace(min,max,1000)
	plt.figure()
	plt.errorbar(I0s[:,0],I0s[:,1],yerr = I0s[:,2],fmt='x',capsize=4)
	plt.plot(Is,Isat(Is,*fit))
	plt.show()
	plt.figure()
	plt.errorbar(ns[:,0],ns[:,1],yerr = ns[:,2],fmt='x',capsize=4)
	plt.ylabel('$n_{id}$ [1]')
	plt.xlabel('T [K]')
	ax1.grid()
	ax2.grid()
	plt.show()
	for i in zip(pIsat,fit,np.sqrt(np.diag(cov)),Iparamunits):
		print(i[0],':',i[1],'+-',i[2],i[3])

final_plots()

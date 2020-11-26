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
		if(np.abs(ret[i,1]-ret_avg(ret[i,0]-0.01))> 1 or np.abs(ret[i,1]-ret_avg(ret[i,0]+0.01))> 1):
			#pass
			mask[i] = False
	return ret[mask]

R = 15.08
cum_data = []
I0s = []
ns = []
plt.ion()
dtype = [('V',np.float64,()),('I',np.float64,())]
IStep = 0.00001
V_err = 0.1

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

prefix = sys.argv[1]
suffix = '.csv'
names = ['LN','0','10','20','30','40','50','60','70','80','90','100']
temps = [77,273,283,293,303,313,323,333,343,353,363,373]

for i in zip(names,temps):
	try:
		try:
			data = load_data(prefix+i[0]+suffix)
		except Exception as e:
			print("Could not open file")
			print(e)
			continue
		T = i[1]
		print(T,'K:')
		true_V = data['V'] - data['I']*R
		cum_data.append([data['V'],data['I'],true_V])
		der = derive(true_V,data['I'],n=1)
		index = der[:,1].argmax()
		V0 = der[index,0]
		fig, (ax1,ax2) = plt.subplots(2,1,sharex=True)
		ax1.plot(true_V,data['I'],'-')
		ax1.set_ylabel('I [A]')
		ax2.set_xlabel('V [V]')
		ax2.plot(der[:,0],der[:,1])
		ax2.plot(der[index,0],der[index,1],'rx')
		ax1.grid()
		ax2.grid()
		ax1.set_yscale('log')
		in_V0 = input("V0: ")
		if(in_V0!=''):
			V0 = float(in_V0)
		loc = np.where(np.abs(V0 - true_V) < V_err)
		fit,cov = fit_shockley(true_V[loc],np.log(data['I'][loc]))
		n = const.e/(der[index,1]*const.k*T)
		n_err = V_err*const.e/(const.k*T*der[index,1]**2)
		I0 = np.exp(np.log(data['I'][index])-der[index,1]*der[index,0])
		I0_err = V_err*np.sqrt(2)*I0*der[index,1]
		der2 = derive(der[:,0],der[:,1],n=1)
		ax1.plot(true_V[loc],shockley(true_V[loc],*fit))
		fig.canvas.show()
		fig.canvas.draw()
		fig.canvas.flush_events()
		for i in zip(params,fit,np.sqrt(np.diag(cov)),paramunits):
			print(i[0],':',i[1],'+-',i[2],i[3])
		#if(input('Should we append this [y] :') != 'y'):
		#	continue
		#if fit[1] > 2.1 or fit[1] < 1:
		#	continue
		I0s.append([T,fit[0],np.sqrt(cov[0,0])])
		ns.append([T,fit[1],np.sqrt(cov[1,1])])
		#print('n :',n,'+-',n_err)
		#print('I0 :',I0,'+-',I0_err)
		#I0s.append([T,I0,I0_err])
		#ns.append([T,n,n_err])
		print('')
	except Exception as e:
		print(e)

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
	der = derive(true_V,data['I'],n=2)
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

nf = None

def final_plots():
	global ns,I0s,nf
	I0s = np.array(I0s)
	ns = np.array(ns)
	nfit, ncov = curve_fit(lambda x, a, b : a*x + b, ns[:,0],ns[:,1],sigma=ns[:,2], absolute_sigma=True)
	nf = interp1d(moving_average(ns[:,0],n=1),moving_average(ns[:,1],n=1),fill_value='extrapolate',kind=1)
	def Isat(x,c,Eg0):
		global nf
		#Eg = lambda T : Eg0 - a*1e-6*T**2/(T+b)
		#return c*x**3 *np.exp(-const.e*Eg0/(const.k*x))
		return np.log(c*x**2) - const.e*Eg0/(const.k*nf(x)*x)
	pIsat = ['c','Eg','a','b']
	Iparamunits = ['A/K^3','eV','ueV/K','K']
	#true_I = I0s[:,1]/ns[:,1]
	#true_I_err = np.sqrt((I0s[:,2]/true_I)**2 + (ns[:,2]/ns[:,1])**2)*true_I
	fit,cov = curve_fit(Isat,I0s[:,0],np.log(I0s[:,1]),sigma=I0s[:,2]/I0s[:,1],bounds=np.array([[0,np.inf],[0,np.inf]]).transpose(),absolute_sigma=True,maxfev = 100000)
	#fit,cov = curve_fit(Isat,I0s[:,0],true_I,sigma=true_I_err,bounds=np.array([[0,np.inf],[0,np.inf],[-np.inf,np.inf],[-np.inf,np.inf]]).transpose(),p0=[1,1,0.5,200],absolute_sigma=True,maxfev = 10000)
	min = I0s[:,0].min()
	max = I0s[:,0].max()
	Is = np.linspace(min,max,1000)
	fig, (ax1,ax2) = plt.subplots(2,1,sharex=True)
	ax1.errorbar(I0s[:,0],I0s[:,1],yerr = I0s[:,2],fmt='x',capsize=4)
	#ax1.errorbar(I0s[:,0],true_I,yerr = true_I_err,fmt='x',capsize=4)
	ax1.plot(Is,np.exp(Isat(Is,*fit)))
	ax1.grid()
	ax1.set_ylabel('$I_0$ [A]')
	ax1.set_yscale('log')
	ax2.errorbar(ns[:,0],ns[:,1],yerr = ns[:,2],fmt='x',capsize=4)
	ax2.plot(Is,nf(Is))
	ax2.set_ylabel('$n_{id}$ [1]')
	ax2.set_xlabel('T [K]')
	ax2.grid()
	plt.show()
	cov[1,1] = cov[1,1] + ns[:,2].mean()*fit[1]/ns[:,1].mean()**2
	for i in zip(pIsat,fit,np.sqrt(np.diag(cov)),Iparamunits):
		print(i[0],':',i[1],'+-',i[2],i[3])

final_plots()

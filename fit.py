#!/usr/bin/ipython3
import sys
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import scipy.constants as const
from scipy.optimize import curve_fit
from scipy.misc import derivative
from scipy.interpolate import interp1d

R = 15.08
cum_data = []
I0s = []
ns = []
plt.ion()
dtype = [('V',np.float64,()),('I',np.float64,())]

def load_data(fname):
    data = np.loadtxt(fname,dtype=dtype,skiprows=0,delimiter=',')
    mask = np.ones(data['V'].size, dtype=bool)
    for i in range(data.size):
        if i == 0:
            continue
        for d in range(1,i+1):
            if not mask[i-d]:
                continue
            # Monotonity check
            if data['V'][i-d] > data['V'][i] or data['I'][i-d] > data['I'][i]:
                mask[i] = False
                break
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
    except:
        print("Could not open file")
        continue
    T = float(input("Temperature of dataset [K]: "))
    true_V = data['V'] - data['I']*R
    cum_data.append([data['V'],data['I'],true_V])
    interp = lambda x: np.exp(interp1d(true_V[9:],np.log(data['I'][9:]),kind=9,fill_value='extrapolate')(x))
    Vs = np.arange(true_V[10],true_V[-1],0.01)
    der = derivative(interp,Vs,n=2,dx=1e-6)
    fig, (ax1,ax2) = plt.subplots(2,1,sharex=True)
    ax1.plot(true_V,data['I'],'-')
    ax1.set_yscale('log')
    ax1.plot(Vs,interp(Vs))
    ax1.set_ylabel('I [A]')
    ax2.set_xlabel('V [V]')
    ax2.plot(Vs,der)
    fig.canvas.show()
    while(True):
        xmin = float(input("Lower limit of fitting: "))
        xmax = float(input("Upper limit of fitting: "))
        try:
            loc = np.where(np.logical_and(xmin < true_V, xmax > true_V))
            #fit,cov = curve_fit(shockley,true_V[loc],data['I'][loc],p0=[1e-50,1],maxfev=10000)
            fit,cov = fit_shockley(true_V[loc],np.log(data['I'][loc]))
        except Exception as e:
            print(e)
            continue
        break
    ax1.plot(true_V[loc],shockley(true_V[loc],*fit))
    fig.canvas.draw()
    fig.canvas.flush_events()
    for i in zip(params,fit,np.sqrt(np.diag(cov)),paramunits):
        print(i[0],':',i[1],'+-',i[2],i[3])
    I0s.append([T,fit[0],np.sqrt(cov[0,0])])
    ns.append([T,fit[1],np.sqrt(cov[1,1])])

def Isat(x,c,Eg):
    return c*x**3 * np.exp(-const.e*Eg/(const.k*x))
pIsat = ['c','Eg']
Iparamunits = ['A/K^3','eV']

def refit_data(index):
    dataloc = input("Dataset: ")
    if(dataloc == ''):
        return
    try:
        data = load_data(prefix+dataloc+suffix)
    except:
        print("Could not open file")
        return
    T = float(input("Temperature of dataset [K]: "))
    true_V = data['V'] - data['I']*R
    cum_data[index] = ([data['V'],data['I'],true_V])
    interp = lambda x: np.exp(interp1d(true_V[9:],np.log(data['I'][9:]),kind=9,fill_value='extrapolate')(x))
    Vs = np.arange(true_V[10],true_V[-1],0.01)
    der = derivative(interp,Vs,n=2,dx=1e-6)
    fig, (ax1,ax2) = plt.subplots(2,1,sharex=True)
    ax1.plot(true_V,data['I'],'-')
    ax1.set_yscale('log')
    ax1.plot(Vs,interp(Vs))
    ax1.set_ylabel('I [A]')
    ax2.set_xlabel('V [V]')
    ax2.plot(Vs,der)
    fig.canvas.show()
    while(True):
        xmin = float(input("Lower limit of fitting: "))
        xmax = float(input("Upper limit of fitting: "))
        try:
            loc = np.where(np.logical_and(xmin < true_V, xmax > true_V))
            #fit,cov = curve_fit(shockley,true_V[loc],data['I'][loc],p0=[1e-50,1],maxfev=10000)
            fit,cov = fit_shockley(true_V[loc],np.log(data['I'][loc]))
        except Exception as e:
            print(e)
            continue
        break
    ax[-1].plot(true_V[loc],shockley(true_V[loc],*fit))
    fig.canvas.draw()
    fig.canvas.flush_events()
    for i in zip(params,fit,np.sqrt(np.diag(cov)),paramunits):
        print(i[0],':',i[1],'+-',i[2],i[3])
    I0s[index] = [T,fit[0],np.sqrt(cov[0,0])]
    ns[index] = [T,fit[1],np.sqrt(cov[1,1])]

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
    plt.show()
    for i in zip(pIsat,fit,np.sqrt(np.diag(cov)),Iparamunits):
        print(i[0],':',i[1],'+-',i[2],i[3])

final_plots()

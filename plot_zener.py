#!/usr/bin/ipython3
import sys
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.interpolate import interp1d

dtype = [('V',np.float64,()),('I',np.float64,())]
plt.ion()
R = 15.08
IStep = 0.00001
args = sys.argv
cmap = cm.get_cmap('gnuplot')

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
			pass
			#mask[i] = False
	return ret[mask]

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

plt.figure()
plt.grid()
plt.xlabel('Voltage [V]')
plt.ylabel('Current [Np]')

Vgs = []

for i in range(11):
	d = i*10
	try:
		data1 = load_data(args[1] + str(d) + '.csv')
	except:
		continue
	true_V = data1['V']-R*data1['I']
	der = derive(true_V,data1['I'],n=2)
	loc = np.where(der[:,1] > 10)[0]
	loc = loc[np.where(loc > 10)[0][0]]
	refloc = np.where(true_V < der[loc,0])[0][-1]

	plt.plot(true_V,np.log(data1['I']/data1['I'][refloc]),label=str(d+273) + ' K',color=cmap(d/150))
	plt.axvline(der[loc,0],linestyle='--',color=cmap(d/150))
	print('V:',der[loc,0],'+-',0.005*der[loc,0],'V')
	Vgs.append([d+273,der[loc,0],0.005*der[loc,0]])

Vgs = np.array(Vgs)
plt.legend()
plt.show()
plt.figure()
plt.errorbar(*Vgs.transpose(),capsize = 4, fmt = 'x')
plt.grid()
plt.ylabel('$V_{breakdown}$ [V]')
plt.xlabel('T [K]')
plt.show()

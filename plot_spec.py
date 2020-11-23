#!/usr/bin/ipython3
import sys
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

dtype = [('l',np.float64,()),('I',np.float64,())]
plt.ion()
args = sys.argv
data1 = np.loadtxt(args[1],dtype=dtype,skiprows=33,comments='[',delimiter=';')
data2 = np.loadtxt(args[2],dtype=dtype,skiprows=33,comments='[',delimiter=';')

plt.plot(data1['l'],data1['I'],label='RT')
plt.plot(data2['l'],data2['I'],label='LN')
plt.legend()
plt.grid()
plt.xlabel('$\lambda$ [nm]')
plt.ylabel('Intensity [a.u.]')
plt.show()

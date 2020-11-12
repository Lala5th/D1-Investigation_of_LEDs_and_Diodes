import pyvisa
import sys
import numpy as np
import matplotlib.pyplot as plt
from time import sleep

R = 15.08
plt.ion()
args = sys.argv
rm = pyvisa.ResourceManager()
inst = rm.open_resource('USB0::0x05E6::0x2450::04392009::INSTR')

if inst.query("*LANG?") != "TSP\n":
	inst.write("*LANG TSP")

inst.write('smu.measure.func = smu.FUNC_DC_CURRENT')
inst.write('smu.source.func  = smu.FUNC_DC_VOLTAGE')
inst.write('smu.source.ilimit.level = 0.2')

V0 = float(args[1])
Vm = float(args[2])
Vstep = float(args[3])

if Vm > 5 or V0 > 5:
	sys.exit()

Vs = np.arange(V0,Vm,Vstep)
data = []
inst.write("defbuffer1.clear()")
for V in Vs:
	inst.write("smu.source.level = " + str(V))
	sleep(0.01)
	inst.write("smu.source.output = smu.ON")
	d = inst.query("print(smu.measure.read(defbuffer1),defbuffer1.sourcevalues[1])")
	inst.write("smu.source.output = smu.OFF")
	ds = d.split('\t')
	inst.write("defbuffer1.clear()")
	if(float(ds[1]) > 1000):
		continue
	print(d)
	data.append([float(ds[1]),float(ds[0])])


save_data = np.array(data)
plt.figure()
plt.plot(save_data[:,0],save_data[:,1])
plt.show()

plt.figure()
plt.plot(save_data[:,0]-R*save_data[:,1],save_data[:,1])
plt.yscale('log')
plt.show()

fname = args[4]
np.savetxt(fname,save_data,delimiter=',')

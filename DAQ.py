import pyvisa
import sys
import numpy as np
import matplotlib.pyplot as plt
from time import sleep

plt.ion()
args = sys.argv
rm = pyvisa.ResourceManager()
inst = rm.open_resource('USB0::0x05E6::0x2450::04392009::INSTR')

if inst.query("*LANG?") != "TSP\n":
	inst.write("*LANG TSP")

inst.write('smu.measure.func = smu.FUNC_DC_VOLTAGE')
inst.write('smu.source.func  = smu.FUNC_DC_CURRENT')
inst.write('smu.source.vlimit.level = 5.5')
inst.write("smu.source.output = smu.ON")

I0 = float(args[1])
Im = float(args[2])
Istep = float(args[3])

if Im >= 0.3 or I0 >= 0.3:
	sys.exit()

Is = np.arange(I0,Im,Istep)
data = []
inst.write("defbuffer1.clear()")
for I in Is:
	inst.write("smu.source.level = " + str(I))
	#sleep(0.001)
	d = inst.query("print(smu.measure.read(defbuffer1),defbuffer1.sourcevalues[1])")
	print(d)
	inst.write("defbuffer1.clear()")
	ds = d.split('\t')
	data.append([float(ds[0]),float(ds[1])])


inst.write("smu.source.output = smu.OFF")
save_data = np.array(data)
plt.plot(save_data[:,0],save_data[:,1])

fname = args[4]
np.savetxt(fname,save_data,delimiter=',')
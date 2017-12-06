import matplotlib
matplotlib.use('Agg')
from time_death import *
from pysb.simulator.base import SimulationResult
import matplotlib.pyplot as plt
import numpy as np

sim_result = SimulationResult.load('EARM_Simulations1k.h5')
tspan = np.linspace(0, 20000, 100)

list_results = []

for sim_array in sim_result.observables:
    observable = sim_array['cPARP']
    xdata = tspan
    time_death = curve_fit_ftn(sig_apop,'cPARP',xdata,sim_array,**{'p0': [100, 100, 100]})
    list_results.append(time_death)

#creating list of tuples for Td:parameter idx
param_idx = list(range(1000))
tuple_list = zip(list_results,param_idx)

#time_death_filt = [x for x in list_results if x>0 and x<20000]
#  #907
tuple_list_filt = [x for x in tuple_list if x>0 and x<20000]
#tuple_list_filt.sort(key=lambda x: float(x[0]))
tuple_list_filt.sort(key=float)
slowest = tuple_list_filt[0:100]
fastest = tuple_list_filt[807:907]

#print(time_death_filt)
#plt.hist(time_death_filt,bins=25)
#y = mlab.normpdf(bins, mu, sigma)
#plt.plot(bins, y, 'r--')
#plt.xlabel('time(s)')
#plt.ylabel('cell count')
#plt.title("Time of Death")
#plt.savefig('Histogram_20ksim12_5_td')


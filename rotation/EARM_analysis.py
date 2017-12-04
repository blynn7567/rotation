import matplotlib
matplotlib.use('Agg')
from time_death import *
from pysb.simulator.base import SimulationResult
import matplotlib.pyplot as plt
import numpy as np

sim_result = SimulationResult.load('EARM_Simulations20k.h5')
tspan = np.linspace(0, 20000, 100)

list_results = []

for sim_array in sim_result.observables:
    observable = sim_array['cPARP']
    xdata = tspan
    time_death = curve_fit_ftn(sig_apop,'cPARP',xdata,sim_array,**{'p0': [100, 100, 100]})
    list_results.append(time_death)

time_death_filt = [x for x in list_results if x>0 and x<20000]
print(time_death_filt)
plt.hist(time_death_filt,bins=25)
plt.title("Histogram Time of Death")
plt.savefig('Histogram_20ksim_td')


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
    non_norm_observable = sim_array['cPARP']
    max_value = max(non_norm_observable)
    observable = non_norm_observable/max_value
    xdata = tspan
    time_death = curve_fit_ftn(sig_apop,observable,xdata,sim_array,**{'p0': [100, 100, 100]})
    list_results.append(time_death)

plt.hist(list_results,bins='auto')
plt.title("Histogram Time of Death")
plt.savefig('Histogram_1ksim_td')


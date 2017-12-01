from time_death import *
from pysb.simulator.base import SimulationResult
import matplotlib.pyplot as plt
import nump as np

sim_result = SimulationResult.load('EARM_Simulations1k.h5')
tspan = np.linspace(0, 20000, 100)

list_results = [0] * 1000

for sim_array in sim_result.observables:
    observable = sim_array['cPARP']
    xdata = tspan
    functions = curve_fit(sig_apop(),'cPARP',xdata,sim_array)
    list.append(list_results)

plt.hist(list_results,bin='auto')
plt.title("Histogram Time of Death")
plt.save('Histogram_1ksim_td')


import matplotlib
matplotlib.use('Agg')
from time_death import *
from pysb.simulator.base import SimulationResult
import matplotlib.pyplot as plt
import numpy as np

sim_result = SimulationResult.load('EARM_Simulations12_6.h5')
tspan = np.linspace(0, 20000, 100)

list_results = []

for sim_array in sim_result.observables:
    observable = sim_array['cPARP']
    xdata = tspan
    time_death = curve_fit_ftn(sig_apop,'cPARP',xdata,sim_array,**{'p0': [100, 100, 100]})
    list_results.append(time_death)

#creating list of tuples for Td:parameter idx
param_idx = list(range(1000))
#tuple_list = np.array(zip(list_results,param_idx))
tuple_list = zip(list_results,param_idx)
#time_death_filt = [x for x in list_results if x>0 and x<20000]
# 907 afer filter
#filter_set = list(range(20000))
tuple_list_filt = [tup for tup in tuple_list if 0 < tup[0] < 20000]
tuple_list_filt.sort(key=lambda x: float(x[0]), reverse=True)
#tuple_list_filt = np.where(np.logical_and(tuple_list>=0,tuple-list<=20000))
#tuple_list_filt.sort(key=float)
slowest = tuple_list_filt[0:100]
slow_td,slow_idx = zip(*slowest)
fastest = tuple_list_filt[807:907]
fast_td,fast_idx = zip(*fastest)

#Now use param index to plot IC of species in tails of histogram
sim_par = np.load('par_dist.npy')
fetched_params_slow = [sim_par[i] for i in slow_idx]
fetched_params_fast = [sim_par[i] for i in fast_idx]

bid_ic_fast = [item[55] for item in fetched_params_fast]
bid_ic_slow = [item[55] for item in fetched_params_slow]
mean_bid_fast = sum(bid_ic_fast)/ float(len(bid_ic_fast))
mean_bid_slow = sum(bid_ic_slow)/ float(len(bid_ic_slow))
min_val_bid_fast = min(bid_ic_fast)
max_val_bid_fast = max(bid_ic_fast)
min_val_bid_slow = min(bid_ic_slow)
max_val_bid_slow = max(bid_ic_slow)

#graph ic bid:
plt.hist(bid_ic_slow,bins=10)
plt.xlabel('[bid]')
plt.ylabel('density')
plt.title("Slow Td: Initial Bid Concetrations")
plt.savefig('Hist_bid_ic2')


#For plotting td histogram
#print(time_death_filt)
#plt.hist(time_death_filt,bins=25)
#y = mlab.normpdf(bins, mu, sigma)
#plt.plot(bins, y, 'r--')
#plt.xlabel('time(s)')
#plt.ylabel('cell count')
#plt.title("Time of Death")
#plt.savefig('Histogram_20ksim12_5_td')


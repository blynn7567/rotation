import matplotlib
matplotlib.use('Agg')
from time_death import *
from pysb.simulator.base import SimulationResult
import matplotlib.pyplot as plt
import numpy as np

sim_result = SimulationResult.load('EARM_Simulations12_6_scaled.h5')
tspan = np.linspace(0, 20000, 100)

list_results = []

for sim_array in sim_result.observables:
    observable = sim_array['cPARP']
    xdata = tspan
    time_death = curve_fit_ftn(sig_apop,'cPARP',xdata,sim_array,**{'p0': [100, 100, 100]})
    list_results.append(time_death)

#creating list of tuples for Td:parameter idx
param_idx = list(range(20000))
#tuple_list = np.array(zip(list_results,param_idx))
tuple_list = zip(list_results,param_idx)
#time_death_filt = [x for x in list_results if x>0 and x<20000]
# 907 afer filter
#filter_set = list(range(20000))
tuple_list_filt = [tup for tup in tuple_list if 0 < tup[0] < 20000]
tuple_list_filt.sort(key=lambda x: float(x[0]), reverse=True)
#tuple_list_filt = np.where(np.logical_and(tuple_list>=0,tuple-list<=20000))
#tuple_list_filt.sort(key=float)
slowest = tuple_list_filt[0:2000]
slow_td,slow_idx = zip(*slowest)
fastest = tuple_list_filt[17993:19993]
fast_td,fast_idx = zip(*fastest)

#Now use param index to plot IC of species in tails of histogram
sim_par = np.load('par_dist_scaled.npy')
fetched_params_slow = [sim_par[i] for i in slow_idx]
fetched_params_fast = [sim_par[i] for i in fast_idx]

#Graphing Bid fast/slow first:
bid_ic_fast = [item[55] for item in fetched_params_fast]
bid_ic_slow = [item[55] for item in fetched_params_slow]
mean_bid_fast = sum(bid_ic_fast)/ float(len(bid_ic_fast))
mean_bid_slow = sum(bid_ic_slow)/ float(len(bid_ic_slow))
min_val_bid_fast = min(bid_ic_fast)
max_val_bid_fast = max(bid_ic_fast)
min_val_bid_slow = min(bid_ic_slow)
max_val_bid_slow = max(bid_ic_slow)

xmin = min(min_val_bid_fast,min_val_bid_slow) # choosing automatic xaxis ranges based off of parsing
xmax = max(max_val_bid_fast,max_val_bid_slow)

plt.figure(1,figsize=(12, 9))
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("molecules Bid", fontsize=16)
plt.ylabel("count", fontsize=16)
plt.title("Fast & Slow Td: Initial Bid Concentrations")
plt.hist(bid_ic_fast, range=[xmin,xmax], color="#3F5D7D", bins=100,alpha=0.5,label='Bid Fast(Td)') #dark blue
plt.hist(bid_ic_slow, range=[xmin,xmax], color="#dbdb8d", bins=100,alpha=0.5,label='Bid Slow(Td)') #light pea?
plt.legend(loc='upper right')
plt.savefig('scaled_bid_c')

#Graphing Bak fast/slow
bak_ic_fast = [item[64] for item in fetched_params_fast]
bak_ic_slow = [item[64] for item in fetched_params_slow]
mean_bak_fast = sum(bak_ic_fast)/ float(len(bak_ic_fast))
mean_bak_slow = sum(bak_ic_slow)/ float(len(bak_ic_slow))
min_val_bak_fast = min(bak_ic_fast)
max_val_bak_fast = max(bak_ic_fast)
min_val_bak_slow = min(bak_ic_slow)
max_val_bak_slow = max(bak_ic_slow)

xmin = min(min_val_bak_fast,min_val_bak_slow) # choosing automatic xaxis ranges based off of parsing
xmax = max(max_val_bak_fast,max_val_bak_slow)

plt.figure(2,figsize=(12, 9))
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("molecules Bak", fontsize=16)
plt.ylabel("count", fontsize=16)
plt.title("Fast & Slow Td: Initial Bak Concentrations")
plt.hist(bak_ic_fast, range=[xmin,xmax], color="#3F5D7D", bins=100,alpha=0.5,label='Bak Fast(Td)') #dark blue
plt.hist(bak_ic_slow, range=[xmin,xmax], color="#dbdb8d", bins=100,alpha=0.5,label='Bak Slow(Td)') #light pea?
plt.legend(loc='upper right')
plt.savefig('scaled_bak_c')

#Graphing Mcl1 fast/slow
mcl1_ic_fast = [item[57] for item in fetched_params_fast]
mcl1_ic_slow = [item[57] for item in fetched_params_slow]
mean_mcl1_fast = sum(mcl1_ic_fast)/ float(len(mcl1_ic_fast))
mean_mcl1_slow = sum(mcl1_ic_slow)/ float(len(mcl1_ic_slow))
min_val_mcl1_fast = min(mcl1_ic_fast)
max_val_mcl1_fast = max(mcl1_ic_fast)
min_val_mcl1_slow = min(mcl1_ic_slow)
max_val_mcl1_slow = max(mcl1_ic_slow)

xmin = min(min_val_mcl1_fast,min_val_mcl1_slow) # choosing automatic xaxis ranges based off of parsing
xmax = max(max_val_mcl1_fast,max_val_mcl1_slow)

plt.figure(3,figsize=(12, 9))
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("molecules Mcl1", fontsize=16)
plt.ylabel("count", fontsize=16)
plt.title("Fast & Slow Td: Initial Mcl1 Concentrations")
plt.hist(mcl1_ic_fast, range=[xmin,xmax], color="#3F5D7D", bins=100,alpha=0.5,label='Mcl1 Fast(Td)') #dark blue
plt.hist(mcl1_ic_slow, range=[xmin,xmax], color="#dbdb8d", bins=100,alpha=0.5,label='Mcl1 Slow(Td)') #light pea?
plt.legend(loc='upper right')
plt.savefig('scaled_mcl1_c')

#Graphing C3 fast/slow
c3_ic_fast = [item[19] for item in fetched_params_fast]
c3_ic_slow = [item[19] for item in fetched_params_slow]
mean_c3_fast = sum(c3_ic_fast)/ float(len(c3_ic_fast))
mean_c3_slow = sum(c3_ic_slow)/ float(len(c3_ic_slow))
min_val_c3_fast = min(c3_ic_fast)
max_val_c3_fast = max(c3_ic_fast)
min_val_c3_slow = min(c3_ic_slow)
max_val_c3_slow = max(c3_ic_slow)

xmin = min(min_val_c3_fast,min_val_c3_slow) # choosing automatic xaxis ranges based off of parsing
xmax = max(max_val_c3_fast,max_val_c3_slow)

plt.figure(4,figsize=(12, 9))
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("molecules Caspase-3", fontsize=16)
plt.ylabel("count", fontsize=16)
plt.title("Fast & Slow Td: Initial Caspase-3 Concentrations")
plt.hist(c3_ic_fast, range=[xmin,xmax], color="#3F5D7D", bins=100,alpha=0.5,label='Caspase-3 Fast(Td)') #dark blue
plt.hist(c3_ic_slow, range=[xmin,xmax], color="#dbdb8d", bins=100,alpha=0.5,label='Caspase-3 Slow(Td)') #light pea?
plt.legend(loc='upper right')
plt.savefig('scaled_c3_c')

#Graphing C6 fast/slow
c6_ic_fast = [item[20] for item in fetched_params_fast]
c6_ic_slow = [item[20] for item in fetched_params_slow]
mean_c6_fast = sum(c6_ic_fast)/ float(len(c6_ic_fast))
mean_c6_slow = sum(c6_ic_slow)/ float(len(c6_ic_slow))
min_val_c6_fast = min(c6_ic_fast)
max_val_c6_fast = max(c6_ic_fast)
min_val_c6_slow = min(c6_ic_slow)
max_val_c6_slow = max(c6_ic_slow)

xmin = min(min_val_c6_fast,min_val_c6_slow) # choosing automatic xaxis ranges based off of parsing
xmax = max(max_val_c6_fast,max_val_c6_slow)

plt.figure(5,figsize=(12, 9))
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("molecules Caspase-6", fontsize=16)
plt.ylabel("count", fontsize=16)
plt.title("Fast & Slow Td: Initial Caspase-6 Concentrations")
plt.hist(c6_ic_fast, range=[xmin,xmax], color="#3F5D7D", bins=100,alpha=0.5,label='Caspase-6 Fast(Td)') #dark blue
plt.hist(c6_ic_slow, range=[xmin,xmax], color="#dbdb8d", bins=100,alpha=0.5,label='Caspase-6 Slow(Td)') #light pea?
plt.legend(loc='upper right')
plt.savefig('scaled_c6_c')
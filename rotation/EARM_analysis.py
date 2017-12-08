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

#Graphing Bid fast first:
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

#graph ic bid:
#plt.hist(bid_ic_slow,bins=10)
#plt.xlabel('[bid]')
#plt.ylabel('density')
#plt.title("Slow Td: Initial Bid Concentrations")
#plt.savefig('Hist_bid_ic2')

# You typically want your plot to be ~1.33x wider than tall.
# Common sizes: (10, 7.5) and (12, 9)
plt.figure(1,figsize=(12, 9))

# Remove the plot frame lines. They are unnecessary chartjunk.
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Ensure that the axis ticks only show up on the bottom and left of the plot.
# Ticks on the right and top of the plot are generally unnecessary chartjunk.
#ax.get_xaxis().tick_bottom()
#ax.get_yaxis().tick_left()

# Make sure your axis ticks are large enough to be easily read.
# You don't want your viewers squinting to read your plot.
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

# Along the same vein, make sure your axis labels are large
# enough to be easily read as well. Make them slightly larger
# than your axis tick labels so they stand out.
plt.xlabel("molecules Bid", fontsize=16)
plt.ylabel("count", fontsize=16)
plt.title("Fast Td: Initial Bid Concentrations")
plt.hist(bid_ic_fast, range=[xmin,xmax], color="#3F5D7D", bins=100) #dark blue
plt.savefig('scaled_bid_f')
#bid slow:
plt.figure(2,figsize=(12, 9))
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("molecules Bid", fontsize=16)
plt.ylabel("count", fontsize=16)
plt.title("Slow Td: Initial Bid Concentrations")
plt.hist(bid_ic_slow, range=[xmin,xmax], color="#3F5D7D", bins=100)
plt.savefig('scaled_bid_s')
#bcl2 fast
bcl2_ic_fast = [item[58] for item in fetched_params_fast]
bcl2_ic_slow = [item[58] for item in fetched_params_slow]
mean_bcl2_fast = sum(bcl2_ic_fast)/ float(len(bcl2_ic_fast))
mean_bcl2_slow = sum(bcl2_ic_slow)/ float(len(bcl2_ic_slow))
min_val_bcl2_fast = min(bcl2_ic_fast)
max_val_bcl2_fast = max(bcl2_ic_fast)
min_val_bcl2_slow = min(bcl2_ic_slow)
max_val_bcl2_slow = max(bcl2_ic_slow)

xmin = min(min_val_bcl2_fast,min_val_bcl2_slow) # choosing automatic xaxis ranges based off of parsing
xmax = max(max_val_bcl2_fast,max_val_bcl2_slow)

plt.figure(3,figsize=(12, 9))
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("molecules Bcl2", fontsize=16)
plt.ylabel("count", fontsize=16)
plt.title("Fast Td: Initial Bcl2 Concentrations")
plt.hist(bcl2_ic_fast, range=[xmin,xmax], color="#3F5D7D", bins=100)
plt.savefig('scaled_bcl2_f')
#bcl2 slow
plt.figure(4,figsize=(12, 9))
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("molecules Bcl2", fontsize=16)
plt.ylabel("count", fontsize=16)
plt.title("Slow Td: Initial Bcl2 Concentrations")
plt.hist(bcl2_ic_slow, range=[xmin,xmax], color="#3F5D7D", bins=100)
plt.savefig('scaled_bcl2_s')

#Bax fast
bax_ic_fast = [item[63] for item in fetched_params_fast]
bax_ic_slow = [item[63] for item in fetched_params_slow]
mean_bax_fast = sum(bax_ic_fast)/ float(len(bax_ic_fast))
mean_bax_slow = sum(bax_ic_slow)/ float(len(bax_ic_slow))
min_val_bax_fast = min(bax_ic_fast)
max_val_bax_fast = max(bax_ic_fast)
min_val_bax_slow = min(bax_ic_slow)
max_val_bax_slow = max(bax_ic_slow)

xmin = min(min_val_bax_fast,min_val_bax_slow) # choosing automatic xaxis ranges based off of parsing
xmax = max(max_val_bax_fast,max_val_bax_slow)

plt.figure(5,figsize=(12, 9))
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("molecules Bax", fontsize=16)
plt.ylabel("count", fontsize=16)
plt.title("Fast Td: Initial Bax Concentrations")
plt.hist(bax_ic_fast, range=[xmin,xmax], color="#3F5D7D", bins=100)
plt.savefig('scaled_bax_f')
#bax slow
plt.figure(6,figsize=(12, 9))
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("molecules Bax", fontsize=16)
plt.ylabel("count", fontsize=16)
plt.title("Slow Td: Initial Bax Concentrations")
plt.hist(bax_ic_slow, range=[xmin,xmax], color="#3F5D7D", bins=100)
plt.savefig('scaled_bax_s')

#Mcl1 fast
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

plt.figure(7,figsize=(12, 9))
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("molecules Mcl1", fontsize=16)
plt.ylabel("count", fontsize=16)
plt.title("Fast Td: Initial Mcl1 Concentrations")
plt.hist(mcl1_ic_fast, range=[xmin,xmax], color="#3F5D7D", bins=100)
plt.savefig('scaled_mcl1_f')
#mcl1 slow
plt.figure(8,figsize=(12, 9))
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("molecules Mcl1", fontsize=16)
plt.ylabel("count", fontsize=16)
plt.title("Slow Td: Initial Mcl1 Concentrations")
plt.hist(mcl1_ic_slow, range=[xmin,xmax], color="#3F5D7D", bins=100)
plt.savefig('scaled_mcl1_s')

#Bad fast
bad_ic_fast = [item[59] for item in fetched_params_fast]
bad_ic_slow = [item[59] for item in fetched_params_slow]
mean_bad_fast = sum(bad_ic_fast)/ float(len(bad_ic_fast))
mean_bad_slow = sum(bad_ic_slow)/ float(len(bad_ic_slow))
min_val_bad_fast = min(bad_ic_fast)
max_val_bad_fast = max(bad_ic_fast)
min_val_bad_slow = min(bad_ic_slow)
max_val_bad_slow = max(bad_ic_slow)

xmin = min(min_val_bad_fast,min_val_bad_slow) # choosing automatic xaxis ranges based off of parsing
xmax = max(max_val_bad_fast,max_val_bad_slow)

plt.figure(9,figsize=(12, 9))
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("molecules Bad", fontsize=16)
plt.ylabel("count", fontsize=16)
plt.title("Fast Td: Initial Bad Concentrations")
plt.hist(bad_ic_fast, range=[xmin,xmax], color="#3F5D7D", bins=100)
plt.savefig('scaled_bad_f')
#bad slow
plt.figure(10,figsize=(12, 9))
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("molecules Bad", fontsize=16)
plt.ylabel("count", fontsize=16)
plt.title("Slow Td: Initial Bad Concentrations")
plt.hist(bad_ic_slow, range=[xmin,xmax], color="#3F5D7D", bins=100)
plt.savefig('scaled_bad_s')

#noxa fast
noxa_ic_fast = [item[60] for item in fetched_params_fast]
noxa_ic_slow = [item[60] for item in fetched_params_slow]
mean_noxa_fast = sum(noxa_ic_fast)/ float(len(noxa_ic_fast))
mean_noxa_slow = sum(noxa_ic_slow)/ float(len(noxa_ic_slow))
min_val_noxa_fast = min(noxa_ic_fast)
max_val_noxa_fast = max(noxa_ic_fast)
min_val_noxa_slow = min(noxa_ic_slow)
max_val_noxa_slow = max(noxa_ic_slow)

xmin = min(min_val_noxa_fast,min_val_noxa_slow) # choosing automatic xaxis ranges based off of parsing
xmax = max(max_val_noxa_fast,max_val_noxa_slow)

plt.figure(11,figsize=(12, 9))
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("molecules Noxa", fontsize=16)
plt.ylabel("count", fontsize=16)
plt.title("Fast Td: Initial Noxa Concentrations")
plt.hist(noxa_ic_fast, range=[xmin,xmax], color="#3F5D7D", bins=100)
plt.savefig('scaled_noxa_f')
#noxa slow
plt.figure(12,figsize=(12, 9))
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("molecules Noxa", fontsize=16)
plt.ylabel("count", fontsize=16)
plt.title("Slow Td: Initial Noxa Concentrations")
plt.hist(noxa_ic_slow, range=[xmin,xmax], color="#3F5D7D", bins=100)
plt.savefig('scaled_noxa_s')

#CytoC fast
cytoc_ic_fast = [item[61] for item in fetched_params_fast]
cytoc_ic_slow = [item[61] for item in fetched_params_slow]
mean_cytoc_fast = sum(cytoc_ic_fast)/ float(len(cytoc_ic_fast))
mean_cytoc_slow = sum(cytoc_ic_slow)/ float(len(cytoc_ic_slow))
min_val_cytoc_fast = min(cytoc_ic_fast)
max_val_cytoc_fast = max(cytoc_ic_fast)
min_val_cytoc_slow = min(cytoc_ic_slow)
max_val_cytoc_slow = max(cytoc_ic_slow)

xmin = min(min_val_cytoc_fast,min_val_cytoc_slow) # choosing automatic xaxis ranges based off of parsing
xmax = max(max_val_cytoc_fast,max_val_cytoc_slow)

plt.figure(13,figsize=(12, 9))
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("molecules Cytoc", fontsize=16)
plt.ylabel("count", fontsize=16)
plt.title("Fast Td: Initial CytoC Concentrations")
plt.hist(cytoc_ic_fast, range=[xmin,xmax], color="#3F5D7D", bins=100)
plt.savefig('scaled_cytoc_f')
#CytoC slow
plt.figure(14,figsize=(12, 9))
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("molecules Cytoc", fontsize=16)
plt.ylabel("count", fontsize=16)
plt.title("Slow Td: Initial CytoC Concentrations")
plt.hist(cytoc_ic_slow, range=[xmin,xmax], color="#3F5D7D", bins=100)
plt.savefig('scaled_cytoc_s')

#Smac fast
smac_ic_fast = [item[62] for item in fetched_params_fast]
smac_ic_slow = [item[62] for item in fetched_params_slow]
mean_smac_fast = sum(smac_ic_fast)/ float(len(smac_ic_fast))
mean_smac_slow = sum(smac_ic_slow)/ float(len(smac_ic_slow))
min_val_smac_fast = min(smac_ic_fast)
max_val_smac_fast = max(smac_ic_fast)
min_val_smac_slow = min(smac_ic_slow)
max_val_smac_slow = max(smac_ic_slow)

xmin = min(min_val_smac_fast,min_val_smac_slow) # choosing automatic xaxis ranges based off of parsing
xmax = max(max_val_smac_fast,max_val_smac_slow)

plt.figure(15,figsize=(12, 9))
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("molecules Smac", fontsize=16)
plt.ylabel("count", fontsize=16)
plt.title("Fast Td: Initial Smac Concentrations")
plt.hist(smac_ic_fast, range=[xmin,xmax], color="#3F5D7D", bins=100)
plt.savefig('scaled_smac_f')
#smacslow
plt.figure(16,figsize=(12, 9))
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("molecules Smac", fontsize=16)
plt.ylabel("count", fontsize=16)
plt.title("Slow Td: Initial Smac Concentrations")
plt.hist(smac_ic_slow, range=[xmin,xmax], color="#3F5D7D", bins=100)
plt.savefig('scaled_smac_s')

#Bak fast
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

plt.figure(17,figsize=(12, 9))
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("molecules Bak", fontsize=16)
plt.ylabel("count", fontsize=16)
plt.title("Fast Td: Initial Bak Concentrations")
plt.hist(bak_ic_fast, range=[xmin,xmax], color="#3F5D7D", bins=100)
plt.savefig('scaled_bak_f')
#smacslow
plt.figure(18,figsize=(12, 9))
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("molecules Bak", fontsize=16)
plt.ylabel("count", fontsize=16)
plt.title("Slow Td: Initial Bak Concentrations")
plt.hist(bak_ic_slow, range=[xmin,xmax], color="#3F5D7D", bins=100)
plt.savefig('scaled_bak_s')

#C3 fast
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

plt.figure(19,figsize=(12, 9))
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("molecules Caspase-3", fontsize=16)
plt.ylabel("count", fontsize=16)
plt.title("Fast Td: Initial Caspase-3 Concentrations")
plt.hist(c3_ic_fast, range=[xmin,xmax], color="#3F5D7D", bins=100)
plt.savefig('scaled_c3_f')
#c3 slow
plt.figure(20,figsize=(12, 9))
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("molecules Caspase-3", fontsize=16)
plt.ylabel("count", fontsize=16)
plt.title("Slow Td: Initial Caspase-3 Concentrations")
plt.hist(c3_ic_slow, range=[xmin,xmax], color="#3F5D7D", bins=100)
plt.savefig('scaled_c3_s')

#C6 fast
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

plt.figure(21,figsize=(12, 9))
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("molecules Caspase-6", fontsize=16)
plt.ylabel("count", fontsize=16)
plt.title("Fast Td: Initial Caspase-6 Concentrations")
plt.hist(c6_ic_fast, range=[xmin,xmax], color="#3F5D7D", bins=100)
plt.savefig('scaled_c6_f')
#c6 slow
plt.figure(22,figsize=(12, 9))
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("molecules Caspase-6", fontsize=16)
plt.ylabel("count", fontsize=16)
plt.title("Slow Td: Initial Caspase-6 Concentrations")
plt.hist(c6_ic_slow, range=[xmin,xmax], color="#3F5D7D", bins=100)
plt.savefig('scaled_c6_s')

#C9 fast
c9_ic_fast = [item[21] for item in fetched_params_fast]
c9_ic_slow = [item[21] for item in fetched_params_slow]
mean_c9_fast = sum(c9_ic_fast)/ float(len(c9_ic_fast))
mean_c9_slow = sum(c9_ic_slow)/ float(len(c9_ic_slow))
min_val_c9_fast = min(c9_ic_fast)
max_val_c9_fast = max(c9_ic_fast)
min_val_c9_slow = min(c9_ic_slow)
max_val_c9_slow = max(c9_ic_slow)

xmin = min(min_val_c9_fast,min_val_c9_slow) # choosing automatic xaxis ranges based off of parsing
xmax = max(max_val_c9_fast,max_val_c9_slow)

plt.figure(23,figsize=(12, 9))
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("molecules Caspase-9", fontsize=16)
plt.ylabel("count", fontsize=16)
plt.title("Fast Td: Initial Caspase-9 Concentrations")
plt.hist(c9_ic_fast, range=[xmin,xmax], color="#3F5D7D", bins=100)
plt.savefig('scaled_c9_f')
#c9 slow
plt.figure(24,figsize=(12, 9))
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("molecules Caspase-9", fontsize=16)
plt.ylabel("count", fontsize=16)
plt.title("Slow Td: Initial Caspase-9 Concentrations")
plt.hist(c9_ic_slow, range=[xmin,xmax], color="#3F5D7D", bins=100)
plt.savefig('scaled_c9_s')

#For plotting td histogram
#print(time_death_filt)
#plt.hist(time_death_filt,bins=25)
#y = mlab.normpdf(bins, mu, sigma)
#plt.plot(bins, y, 'r--')
#plt.xlabel('time(s)')
#plt.ylabel('cell count')
#plt.title("Time of Death")
#plt.savefig('Histogram_20ksim12_5_td')


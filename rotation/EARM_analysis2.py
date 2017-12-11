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

#Graphing Caspase-8 fast/slow first:
bid_ic_fast = [item[55] for item in fetched_params_fast]
bid_ic_slow = [item[55] for item in fetched_params_slow]
mean_bid_fast = sum(bid_ic_fast)/ float(len(bid_ic_fast))
mean_bid_slow = sum(bid_ic_slow)/ float(len(bid_ic_slow))
min_val_bid_fast = min(bid_ic_fast)
max_val_bid_fast = max(bid_ic_fast)
min_val_bid_slow = min(bid_ic_slow)
max_val_bid_slow = max(bid_ic_slow)
sd_bid_fast = stdev(bid_ic_fast)
sd_bid_slow = stdev(bid_ic_slow)

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
plt.title('Bid Fast(Td)'+r'$\sigma$ ='+str(round(sd_bid_fast, 2))+', ' r'$\mu$ ='+str(round(mean_bid_fast, 2)) + '\n'
          +'Bid Slow(Td)'+r'$\sigma$ ='+str(round(sd_bid_slow, 2))+', ' r'$\mu$ ='+str(round(mean_bid_slow, 2)))
plt.hist(bid_ic_fast, range=[xmin,xmax], color="#3F5D7D", bins=100,alpha=0.5,label='Bid Fast(Td)') #dark blue
plt.hist(bid_ic_slow, range=[xmin,xmax], color="#dbdb8d", bins=100,alpha=0.5,label='Bid Slow(Td)') #light pea?
plt.legend(loc='upper right')
plt.savefig('scaled_bid_X')

#Graphing Bak fast/slow
bak_ic_fast = [item[64] for item in fetched_params_fast]
bak_ic_slow = [item[64] for item in fetched_params_slow]
mean_bak_fast = sum(bak_ic_fast)/ float(len(bak_ic_fast))
mean_bak_slow = sum(bak_ic_slow)/ float(len(bak_ic_slow))
min_val_bak_fast = min(bak_ic_fast)
max_val_bak_fast = max(bak_ic_fast)
min_val_bak_slow = min(bak_ic_slow)
max_val_bak_slow = max(bak_ic_slow)
sd_bak_fast = stdev(bak_ic_fast)
sd_bak_slow = stdev(bak_ic_slow)


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
plt.title('Bak Fast(Td)'+r'$\sigma$ ='+str(round(sd_bak_fast, 2))+', ' r'$\mu$ ='+str(round(mean_bak_fast, 2)) + '\n'
          +'Bak Slow(Td)'+r'$\sigma$ ='+str(round(sd_bak_slow, 2))+', ' r'$\mu$ ='+str(round(mean_bak_slow, 2)))
plt.hist(bak_ic_fast, range=[xmin,xmax], color="#3F5D7D", bins=100,alpha=0.5,label='Bak Fast(Td)') #dark blue
plt.hist(bak_ic_slow, range=[xmin,xmax], color="#dbdb8d", bins=100,alpha=0.5,label='Bak Slow(Td)') #light pea?
plt.legend(loc='upper right')
plt.savefig('scaled_bak_X')

#Graphing Mcl1 fast/slow
mcl1_ic_fast = [item[57] for item in fetched_params_fast]
mcl1_ic_slow = [item[57] for item in fetched_params_slow]
mean_mcl1_fast = sum(mcl1_ic_fast)/ float(len(mcl1_ic_fast))
mean_mcl1_slow = sum(mcl1_ic_slow)/ float(len(mcl1_ic_slow))
min_val_mcl1_fast = min(mcl1_ic_fast)
max_val_mcl1_fast = max(mcl1_ic_fast)
min_val_mcl1_slow = min(mcl1_ic_slow)
max_val_mcl1_slow = max(mcl1_ic_slow)
sd_mcl1_fast = stdev(mcl1_ic_fast)
sd_mcl1_slow = stdev(mcl1_ic_slow)

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
plt.title('Mcl1 Fast(Td)'+r'$\sigma$ ='+str(round(sd_mcl1_fast, 2))+', ' r'$\mu$ ='+str(round(mean_mcl1_fast, 2)) + '\n'
          +'Mcl1 Slow(Td)'+r'$\sigma$ ='+str(round(sd_mcl1_slow, 2))+', ' r'$\mu$ ='+str(round(mean_mcl1_slow, 2)))
plt.hist(mcl1_ic_fast, range=[xmin,xmax], color="#3F5D7D", bins=100,alpha=0.5,label='Mcl1 Fast(Td)') #dark blue
plt.hist(mcl1_ic_slow, range=[xmin,xmax], color="#dbdb8d", bins=100,alpha=0.5,label='Mcl1 Slow(Td)') #light pea?
plt.legend(loc='upper right')
plt.savefig('scaled_mcl1_X')

#Graphing C3 fast/slow
c3_ic_fast = [item[19] for item in fetched_params_fast]
c3_ic_slow = [item[19] for item in fetched_params_slow]
mean_c3_fast = sum(c3_ic_fast)/ float(len(c3_ic_fast))
mean_c3_slow = sum(c3_ic_slow)/ float(len(c3_ic_slow))
min_val_c3_fast = min(c3_ic_fast)
max_val_c3_fast = max(c3_ic_fast)
min_val_c3_slow = min(c3_ic_slow)
max_val_c3_slow = max(c3_ic_slow)
sd_c3_fast = stdev(c3_ic_fast)
sd_c3_slow = stdev(c3_ic_slow)

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
plt.title('Caspase-3 Fast(Td)'+r'$\sigma$ ='+str(round(sd_c3_fast, 2))+', ' r'$\mu$ ='+str(round(mean_c3_fast, 2)) + '\n'
          +'Caspase-3 Slow(Td)'+r'$\sigma$ ='+str(round(sd_c3_slow, 2))+', ' r'$\mu$ ='+str(round(mean_c3_slow, 2)))
plt.hist(c3_ic_fast, range=[xmin,xmax], color="#3F5D7D", bins=100,alpha=0.5,label='Caspase-3 Fast(Td)') #dark blue
plt.hist(c3_ic_slow, range=[xmin,xmax], color="#dbdb8d", bins=100,alpha=0.5,label='Caspase-3 Slow(Td)') #light pea?
plt.legend(loc='upper right')
plt.savefig('scaled_c3_X')

#Graphing C6 fast/slow
c6_ic_fast = [item[20] for item in fetched_params_fast]
c6_ic_slow = [item[20] for item in fetched_params_slow]
mean_c6_fast = sum(c6_ic_fast)/ float(len(c6_ic_fast))
mean_c6_slow = sum(c6_ic_slow)/ float(len(c6_ic_slow))
min_val_c6_fast = min(c6_ic_fast)
max_val_c6_fast = max(c6_ic_fast)
min_val_c6_slow = min(c6_ic_slow)
max_val_c6_slow = max(c6_ic_slow)
sd_c6_fast = stdev(c6_ic_fast)
sd_c6_slow = stdev(c6_ic_slow)

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
plt.title('Caspase-6 Fast(Td)'+r'$\sigma$ ='+str(round(sd_c6_fast, 2))+', ' r'$\mu$ ='+str(round(mean_c6_fast, 2)) + '\n'
          +'Caspase-6 Slow(Td)'+r'$\sigma$ ='+str(round(sd_c6_slow, 2))+', ' r'$\mu$ ='+str(round(mean_c6_slow, 2)))
plt.hist(c6_ic_fast, range=[xmin,xmax], color="#3F5D7D", bins=100,alpha=0.5,label='Caspase-6 Fast(Td)') #dark blue
plt.hist(c6_ic_slow, range=[xmin,xmax], color="#dbdb8d", bins=100,alpha=0.5,label='Caspase-6 Slow(Td)') #light pea?
plt.legend(loc='upper right')
plt.savefig('scaled_c6_X')

#Graphing flip fast/slow
flip_ic_fast = [item[2] for item in fetched_params_fast]
flip_ic_slow = [item[2] for item in fetched_params_slow]
mean_flip_fast = sum(flip_ic_fast)/ float(len(flip_ic_fast))
mean_flip_slow = sum(flip_ic_slow)/ float(len(flip_ic_slow))
min_val_flip_fast = min(flip_ic_fast)
max_val_flip_fast = max(flip_ic_fast)
min_val_flip_slow = min(flip_ic_slow)
max_val_flip_slow = max(flip_ic_slow)
sd_flip_fast = stdev(flip_ic_fast)
sd_flip_slow = stdev(flip_ic_slow)

xmin = min(min_val_flip_fast,min_val_flip_slow) # choosing automatic xaxis ranges based off of parsing
xmax = max(max_val_flip_fast,max_val_flip_slow)

plt.figure(6,figsize=(12, 9))
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("molecules Flip", fontsize=16)
plt.ylabel("count", fontsize=16)
plt.title('Flip Fast(Td)'+r'$\sigma$ ='+str(round(sd_flip_fast, 2))+', ' r'$\mu$ ='+str(round(mean_flip_fast, 2)) + '\n'
          +'Flip Slow(Td)'+r'$\sigma$ ='+str(round(sd_flip_slow, 2))+', ' r'$\mu$ ='+str(round(mean_flip_slow, 2)))
plt.hist(flip_ic_fast, range=[xmin,xmax], color="#3F5D7D", bins=100,alpha=0.5,label='Flip Fast(Td)') #dark blue
plt.hist(flip_ic_slow, range=[xmin,xmax], color="#dbdb8d", bins=100,alpha=0.5,label='Flip Slow(Td)') #light pea?
plt.legend(loc='upper right')
plt.savefig('scaled_flip_X')

#Graphing Bclxl fast/slow
bclxl_ic_fast = [item[56] for item in fetched_params_fast]
bclxl_ic_slow = [item[56] for item in fetched_params_slow]
mean_bclxl_fast = sum(bclxl_ic_fast)/ float(len(bclxl_ic_fast))
mean_bclxl_slow = sum(bclxl_ic_slow)/ float(len(bclxl_ic_slow))
min_val_bclxl_fast = min(bclxl_ic_fast)
max_val_bclxl_fast = max(bclxl_ic_fast)
min_val_bclxl_slow = min(bclxl_ic_slow)
max_val_bclxl_slow = max(bclxl_ic_slow)
sd_bclxl_fast = stdev(bclxl_ic_fast)
sd_bclxl_slow = stdev(bclxl_ic_slow)

xmin = min(min_val_bclxl_fast,min_val_bclxl_slow) # choosing automatic xaxis ranges based off of parsing
xmax = max(max_val_bclxl_fast,max_val_bclxl_slow)

plt.figure(7,figsize=(12, 9))
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("molecules Bclxl", fontsize=16)
plt.ylabel("count", fontsize=16)
plt.title('Bclxl Fast(Td)'+r'$\sigma$ ='+str(round(sd_bclxl_fast, 2))+', ' r'$\mu$ ='+str(round(mean_bclxl_fast, 2)) + '\n'
          +'Bclxl Slow(Td)'+r'$\sigma$ ='+str(round(sd_bclxl_slow, 2))+', ' r'$\mu$ ='+str(round(mean_bclxl_slow, 2)))
plt.hist(bclxl_ic_fast, range=[xmin,xmax], color="#3F5D7D", bins=100,alpha=0.5,label='Bclxl Fast(Td)') #dark blue
plt.hist(bclxl_ic_slow, range=[xmin,xmax], color="#dbdb8d", bins=100,alpha=0.5,label='Bclxl Slow(Td)') #light pea?
plt.legend(loc='upper right')
plt.savefig('scaled_bclxl_X')

#Graphing L_0 fast/slow
ligand_ic_fast = [item[56] for item in fetched_params_fast]
ligand_ic_slow = [item[56] for item in fetched_params_slow]
mean_ligand_fast = sum(ligand_ic_fast)/ float(len(ligand_ic_fast))
mean_ligand_slow = sum(ligand_ic_slow)/ float(len(ligand_ic_slow))
min_val_ligand_fast = min(ligand_ic_fast)
max_val_ligand_fast = max(ligand_ic_fast)
min_val_ligand_slow = min(ligand_ic_slow)
max_val_ligand_slow = max(ligand_ic_slow)
sd_ligand_fast = stdev(ligand_ic_fast)
sd_ligand_slow = stdev(ligand_ic_slow)

xmin = min(min_val_ligand_fast,min_val_ligand_slow) # choosing automatic xaxis ranges based off of parsing
xmax = max(max_val_ligand_fast,max_val_ligand_slow)

plt.figure(8,figsize=(12, 9))
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("molecules Ligand", fontsize=16)
plt.ylabel("count", fontsize=16)
plt.title('Ligand Fast(Td)'+r'$\sigma$ ='+str(round(sd_ligand_fast, 2))+', ' r'$\mu$ ='+str(round(mean_ligand_fast, 2)) + '\n'
          +'Ligand Slow(Td)'+r'$\sigma$ ='+str(round(sd_ligand_slow, 2))+', ' r'$\mu$ ='+str(round(mean_ligand_slow, 2)))
plt.hist(ligand_ic_fast, range=[xmin,xmax], color="#3F5D7D", bins=100,alpha=0.5,label='Ligand Fast(Td)') #dark blue
plt.hist(ligand_ic_slow, range=[xmin,xmax], color="#dbdb8d", bins=100,alpha=0.5,label='Ligand Slow(Td)') #light pea?
plt.legend(loc='upper right')
plt.savefig('scaled_ligand_X')

#Graphing R_0 fast/slow
receptor_ic_fast = [item[1] for item in fetched_params_fast]
receptor_ic_slow = [item[1] for item in fetched_params_slow]
mean_receptor_fast = sum(receptor_ic_fast)/ float(len(receptor_ic_fast))
mean_receptor_slow = sum(receptor_ic_slow)/ float(len(receptor_ic_slow))
min_val_receptor_fast = min(receptor_ic_fast)
max_val_receptor_fast = max(receptor_ic_fast)
min_val_receptor_slow = min(receptor_ic_slow)
max_val_receptor_slow = max(receptor_ic_slow)
sd_receptor_fast = stdev(receptor_ic_fast)
sd_receptor_slow = stdev(receptor_ic_slow)

xmin = min(min_val_receptor_fast,min_val_receptor_slow) # choosing automatic xaxis ranges based off of parsing
xmax = max(max_val_receptor_fast,max_val_receptor_slow)

plt.figure(9,figsize=(12, 9))
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("molecules Receptor", fontsize=16)
plt.ylabel("count", fontsize=16)
plt.title('Receptor Fast(Td)'+r'$\sigma$ ='+str(round(sd_receptor_fast, 2))+', ' r'$\mu$ ='+str(round(mean_receptor_fast, 2)) + '\n'
          +'Receptor Slow(Td)'+r'$\sigma$ ='+str(round(sd_receptor_slow, 2))+', ' r'$\mu$ ='+str(round(mean_receptor_slow, 2)))
plt.hist(receptor_ic_fast, range=[xmin,xmax], color="#3F5D7D", bins=100,alpha=0.5,label='receptor Fast(Td)') #dark blue
plt.hist(receptor_ic_slow, range=[xmin,xmax], color="#dbdb8d", bins=100,alpha=0.5,label='receptor Slow(Td)') #light pea?
plt.legend(loc='upper right')
plt.savefig('scaled_receptor_X')

#Graphing BAR fast/slow
bar_ic_fast = [item[4] for item in fetched_params_fast]
bar_ic_slow = [item[4] for item in fetched_params_slow]
mean_bar_fast = sum(bar_ic_fast)/ float(len(bar_ic_fast))
mean_bar_slow = sum(bar_ic_slow)/ float(len(bar_ic_slow))
min_val_bar_fast = min(bar_ic_fast)
max_val_bar_fast = max(bar_ic_fast)
min_val_bar_slow = min(bar_ic_slow)
max_val_bar_slow = max(bar_ic_slow)
sd_bar_fast = stdev(bar_ic_fast)
sd_bar_slow = stdev(bar_ic_slow)

xmin = min(min_val_bar_fast,min_val_bar_slow) # choosing automatic xaxis ranges based off of parsing
xmax = max(max_val_bar_fast,max_val_bar_slow)

plt.figure(10,figsize=(12, 9))
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("molecules BAR", fontsize=16)
plt.ylabel("count", fontsize=16)
plt.title('BAR Fast(Td)'+r'$\sigma$ ='+str(round(sd_bar_fast, 2))+', ' r'$\mu$ ='+str(round(mean_bar_fast, 2)) + '\n'
          +'BAR Slow(Td)'+r'$\sigma$ ='+str(round(sd_bar_slow, 2))+', ' r'$\mu$ ='+str(round(mean_bar_slow, 2)))
plt.hist(bar_ic_fast, range=[xmin,xmax], color="#3F5D7D", bins=100,alpha=0.5,label='BAR Fast(Td)') #dark blue
plt.hist(bar_ic_slow, range=[xmin,xmax], color="#dbdb8d", bins=100,alpha=0.5,label='BAR Slow(Td)') #light pea?
plt.legend(loc='upper right')
plt.savefig('scaled_bar_X')

#Graphing Apaf fast/slow
apaf_ic_fast = [item[18] for item in fetched_params_fast]
apaf_ic_slow = [item[18] for item in fetched_params_slow]
mean_apaf_fast = sum(apaf_ic_fast)/ float(len(apaf_ic_fast))
mean_apaf_slow = sum(apaf_ic_slow)/ float(len(apaf_ic_slow))
min_val_apaf_fast = min(apaf_ic_fast)
max_val_apaf_fast = max(apaf_ic_fast)
min_val_apaf_slow = min(apaf_ic_slow)
max_val_apaf_slow = max(apaf_ic_slow)
sd_apaf_fast = stdev(apaf_ic_fast)
sd_apaf_slow = stdev(apaf_ic_slow)

xmin = min(min_val_apaf_fast,min_val_apaf_slow) # choosing automatic xaxis ranges based off of parsing
xmax = max(max_val_apaf_fast,max_val_apaf_slow)

plt.figure(11,figsize=(12, 9))
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("molecules Apaf", fontsize=16)
plt.ylabel("count", fontsize=16)
plt.title('Apaf Fast(Td)'+r'$\sigma$ ='+str(round(sd_apaf_fast, 2))+', ' r'$\mu$ ='+str(round(mean_apaf_fast, 2)) + '\n'
          +'Apaf Slow(Td)'+r'$\sigma$ ='+str(round(sd_apaf_slow, 2))+', ' r'$\mu$ ='+str(round(mean_apaf_slow, 2)))
plt.hist(apaf_ic_fast, range=[xmin,xmax], color="#3F5D7D", bins=100,alpha=0.5,label='Apaf Fast(Td)') #dark blue
plt.hist(apaf_ic_slow, range=[xmin,xmax], color="#dbdb8d", bins=100,alpha=0.5,label='Apaf Slow(Td)') #light pea?
plt.legend(loc='upper right')
plt.savefig('scaled_apaf_X')

#Graphing XIAP fast/slow
xiap_ic_fast = [item[22] for item in fetched_params_fast]
xiap_ic_slow = [item[22] for item in fetched_params_slow]
mean_xiap_fast = sum(xiap_ic_fast)/ float(len(xiap_ic_fast))
mean_xiap_slow = sum(xiap_ic_slow)/ float(len(xiap_ic_slow))
min_val_xiap_fast = min(xiap_ic_fast)
max_val_xiap_fast = max(xiap_ic_fast)
min_val_xiap_slow = min(xiap_ic_slow)
max_val_xiap_slow = max(xiap_ic_slow)
sd_xiap_fast = stdev(xiap_ic_fast)
sd_xiap_slow = stdev(xiap_ic_slow)

xmin = min(min_val_xiap_fast,min_val_xiap_slow) # choosing automatic xaxis ranges based off of parsing
xmax = max(max_val_xiap_fast,max_val_xiap_slow)

plt.figure(12,figsize=(12, 9))
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("molecules XIAP", fontsize=16)
plt.ylabel("count", fontsize=16)
plt.title('XIAP Fast(Td)'+r'$\sigma$ ='+str(round(sd_xiap_fast, 2))+', ' r'$\mu$ ='+str(round(mean_xiap_fast, 2)) + '\n'
          +'XIAP Slow(Td)'+r'$\sigma$ ='+str(round(sd_xiap_slow, 2))+', ' r'$\mu$ ='+str(round(mean_xiap_slow, 2)))
plt.hist(xiap_ic_fast, range=[xmin,xmax], color="#3F5D7D", bins=100,alpha=0.5,label='XIAP Fast(Td)') #dark blue
plt.hist(xiap_ic_slow, range=[xmin,xmax], color="#dbdb8d", bins=100,alpha=0.5,label='XIAP Slow(Td)') #light pea?
plt.legend(loc='upper right')
plt.savefig('scaled_xiap_X')

#Graphing PARP fast/slow
parp_ic_fast = [item[23] for item in fetched_params_fast]
parp_ic_slow = [item[23] for item in fetched_params_slow]
mean_parp_fast = sum(parp_ic_fast)/ float(len(parp_ic_fast))
mean_parp_slow = sum(parp_ic_slow)/ float(len(parp_ic_slow))
min_val_parp_fast = min(parp_ic_fast)
max_val_parp_fast = max(parp_ic_fast)
min_val_parp_slow = min(parp_ic_slow)
max_val_parp_slow = max(parp_ic_slow)
sd_parp_fast = stdev(parp_ic_fast)
sd_parp_slow = stdev(parp_ic_slow)

xmin = min(min_val_parp_fast,min_val_parp_slow) # choosing automatic xaxis ranges based off of parsing
xmax = max(max_val_parp_fast,max_val_parp_slow)

plt.figure(13,figsize=(12, 9))
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("molecules PARP", fontsize=16)
plt.ylabel("count", fontsize=16)
plt.title('PARP Fast(Td)'+r'$\sigma$ ='+str(round(sd_parp_fast, 2))+', ' r'$\mu$ ='+str(round(mean_parp_fast, 2)) + '\n'
          +'PARP Slow(Td)'+r'$\sigma$ ='+str(round(sd_parp_slow, 2))+', ' r'$\mu$ ='+str(round(mean_parp_slow, 2)))
plt.hist(parp_ic_fast, range=[xmin,xmax], color="#3F5D7D", bins=100,alpha=0.5,label='PARP Fast(Td)') #dark blue
plt.hist(parp_ic_slow, range=[xmin,xmax], color="#dbdb8d", bins=100,alpha=0.5,label='PARP Slow(Td)') #light pea?
plt.legend(loc='upper right')
plt.savefig('scaled_parp_X')

#Graphing caspase-8 fast/slow
c8_ic_fast = [item[3] for item in fetched_params_fast]
c8_ic_slow = [item[3] for item in fetched_params_slow]
mean_c8_fast = sum(c8_ic_fast)/ float(len(c8_ic_fast))
mean_c8_slow = sum(c8_ic_slow)/ float(len(c8_ic_slow))
min_val_c8_fast = min(c8_ic_fast)
max_val_c8_fast = max(c8_ic_fast)
min_val_c8_slow = min(c8_ic_slow)
max_val_c8_slow = max(c8_ic_slow)
sd_c8_fast = stdev(c8_ic_fast)
sd_c8_slow = stdev(c8_ic_slow)

xmin = min(min_val_c8_fast,min_val_c8_slow) # choosing automatic xaxis ranges based off of parsing
xmax = max(max_val_c8_fast,max_val_c8_slow)

plt.figure(14,figsize=(12, 9))
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("molecules Caspase-8", fontsize=16)
plt.ylabel("count", fontsize=16)
plt.title('Caspase-8 Fast(Td)'+r'$\sigma$ ='+str(round(sd_c8_fast, 2))+', ' r'$\mu$ ='+str(round(mean_c8_fast, 2)) + '\n'
          +'Caspase-8 Slow(Td)'+r'$\sigma$ ='+str(round(sd_c8_slow, 2))+', ' r'$\mu$ ='+str(round(mean_c8_slow, 2)))
plt.hist(c8_ic_fast, range=[xmin,xmax], color="#3F5D7D", bins=100,alpha=0.5,label='Caspase-8 Fast(Td)') #dark blue
plt.hist(c8_ic_slow, range=[xmin,xmax], color="#dbdb8d", bins=100,alpha=0.5,label='Caspase-8 Slow(Td)') #light pea?
plt.legend(loc='upper right')
plt.savefig('scaled_c8_X')



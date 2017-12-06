from earm.lopez_embedded import model
import numpy as np
from pysb.simulator.scipyode import ScipyOdeSimulator
from pysb.simulator.cupsoda import CupSodaSimulator
# import matplotlib.pylab as pl
from pysb import Observable
from pysb.bng import generate_equations
from IC_distributions_earm import sample_lognormal

#CytoC = model.monomers['CytoC']
#Observable('CytoC_O', CytoC(bf=None,state='C'))

model.reset_equations()
generate_equations(model)
# This did not work: Observable('CytoC', model.species[65])

# To change Parameter values (two ways)
# change model file and change parameter values
# Or pass new values to simulator object
#Parameter matrix needs to be indexed by set in rows N= 6572 and paramters are columns and pass this set after naming it to the simulator
parameters = np.load('calibrated_6572pars.npy')
pars1 = parameters[0]
tspan = np.linspace(0, 20000, 100)

#Time to change/sample from initial conditions

parameters_ic = {idx: p for idx, p in enumerate(model.parameters) if p in model.parameters_initial_conditions()[1:]}
# this is the dictionary from model.parameters with index number and parameter for IC of each species
samples = 1000 # gives list of arrays for each simulation so you need to index them in the visualizations

repeated_parameter_values = np.tile(pars1, (samples, 1))
# creating an array each row: set of initial conditions and columns: Parameters
# ----This is where Oscar kept parameters by copy paste and changed initial conditions
for idx, par in parameters_ic.items():
    repeated_parameter_values[:, idx] = sample_lognormal(par, size=samples)
np.save('par_dist.npy',repeated_parameter_values)
# can sample from many different initial conditions defined by mean and coefficient of variance.
# We know these proteins have a log normal distribution
# If you want to change the initial conditions you can change it in two different ways:
# 1) Create a dictioary where the key has to be the name of the parameter in pysb (model.parameters)
#    amd the value is the new value that you want to assign to the parameter and then pass this
#    dictionary to the run function (the argumennt has to look like param_values={})
# We have so many parameter sets because they are generated by the runs that fit the cut-off for model calibration eg: 0.01
#  Remember the error bars from the experimental data that we try to fit to the curve-this can be calibrated by both DREAM
# and Particle Swarm Optimization
# Both depend on the cost_function: square difference between experimental data and simulations
# F(x) = Sumsymb (i,j) (Exp data J (species) over i (time)-Sim(j over i)) squared all over variance of experimental data
# (if you want to normalize)
# rememer the landscape topology Oscar drew, there are local and global minimas we must find within the 105 dimentions
#  of the parameter values.
# Sloppy analysis says that even if you fit the curve, 1 parameter set, and you ave defined experimental data
#  for each species and dimension-a linear range will still be better suited
#PSO places particles in locations (random and smart ways to do this) and they converge towards minimas
#  because they have access to the other's cost_function values) may not need to be a global min
# DREAM sounds more efficient-uses particles with their own velocities that spend more time in minimas and
# therefore the minima may be directly related
# So if Oscar ran 20K simulations and 4 past the cut-off, this is the parameter set.
#  Remember our different initial conditions create the log distribution of species concentrations
#Morkof Chain Monte Carlo refers to DREAM algorithm
# sim = ScipyOdeSimulator(model, tspan=tspan,param_values=repeated_parameter_values).run()

integrator_opt = {'rtol': 1e-6,'atol': 1e-6,'mxsteps': 20000}
vol = 1e-19
sim = CupSodaSimulator(model, tspan=tspan,param_values=repeated_parameter_values,obs_species_only=False).run()
#obs species= F lets us look at concentrations of everything after run
sim.save("EARM_Simulations12_6.h5") #.h5 is HDF format for saving
# After getting the simulation visualize using matplotlib and the observables.
# Excersices: Find conditions under which the cell dies faster or slower relative to the 'wild type simulation"
# show some visualizations
# Take a look at how distributions work, and how to sample a distribution in numpy. Take a look aat log normal
# distributions.

#simres = ScipyOdeSimulator(model, tspan=tspan,param_values=pars1).run()
# yout = sim.all  #species and observables in numpy array or simres.dataframe (pandas dataframe) or simres.species (only species)
# pl.plot(tspan, yout[0]['cPARP'], label="cPARP")
# pl.plot(tspan, yout[1]['cPARP'], label="cPARP")
#pl.plot(tspan, yout[2]['cPARP'], label="cPARP")
#pl.plot(tspan, yout[3]['cPARP'], label="cPARP")
#pl.plot(tspan,yout['mBid'], label="mBid")
##pl.plot(tspan,yout[0]['aSmac'],label="aSmac")
#pl.plot(tspan,yout[1]['aSmac'],label="aSmac")
#pl.plot(tspan,yout[2]['aSmac'],label="aSmac")
#pl.plot(tspan,yout[3]['aSmac'],label="aSmac")
#pl.plot(tspan,yout[0]['CytoC_O'],label='CytoC')
#pl.plot(tspan,yout[1]['CytoC_O'],label='CytoC')
#pl.plot(tspan,yout[2]['CytoC_O'],label='CytoC')
#pl.plot(tspan,yout[3]['CytoC_O'],label='CytoC')# must use refernce from npy array to add new observable
# pl.legend()
# pl.show()

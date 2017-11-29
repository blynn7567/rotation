#import
from firstmodel import model
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pylab as pl
from pysb.simulator import ScipyOdeSimulator

#Simulation Results run() and ODEs and time
t = pl.linspace(0,200)
simres = ScipyOdeSimulator(model, tspan=t).run()
yout = simres.all  #species and observables in numpy array or simres.dataframe (pandas dataframe) or simres.species (only species)

pl.ion()  # CTRL forward slash comments everything out
pl.figure()
pl.plot(t, yout['bound_enzyme'], label="bound_enzyme")
pl.plot(t, yout['obsS'], label="Unbound Substrate")
pl.plot(t, yout['obsP'], label="Product")
pl.legend()
pl.xlabel("Time (s)")
pl.ylabel("Molecules/cell")
pl.show()

# import the pysb module and all its methods and functions
from pysb import *

# instantiate a model
Model()

# declare monomers
Monomer('cyclin', ['b', 'state1'],{'state1': ['u','p']})
Monomer('cdc2', ['b', 'state2'],{'state2': ['u','p']})

#Rules
Rule('cyclin-pcdc2-p_cyclin-pcdc2', 'cyclin'(b=1, state1 = 'p') % 'cdc2'(b=1, state2 = 'p') <> 'cyclin'(b=1, state1 = 'p') % 'cdc2'(b=1, state2 = 'u'), kf, kr)

# input the parameter values
Parameter('kf', 0.018) #p-cdc2%p-cyclin loses p(i) from cdc2 phosphatase Y-15
Parameter('kr', 0.0)   #p-cyclin%cdc2 kinase T-167 --> p-cycline%p-cdc2

#initial conditions
Parameter('cyclin_0', 100)
Parameter('cdc2_0', 100)
Initial(cyclin(b=1, state1 = 'p'), cyclin_0)
Initial(cdc2(b=1, state2='u'), cdc2_0)

#observables
Observable('cyclin_T', cyclin(b=WILD))
Observable('cyclin%cdc2', cyclin(b = 1, state1 = 'p') % cdc2(b=1, state2 = 'u'))
Observable('cyclin%cdc2-p', cyclin(b = 1, state1 = 'p') % cdc2(b=1, state2 = 'p'))
Observable('cdc2', cdc2(b = 1, state2 = 'p'))
Observable('cdc2', cdc2(b = 1, state2 = 'u'))

#visualize in seperate .py



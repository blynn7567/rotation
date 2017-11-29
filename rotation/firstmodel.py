# import the pysb module and all its methods and functions
from pysb import *
from pysb.simulator import ScipyOdeSimulator

#using macro for reversible catalysis
def catalyze(enz, sub, site, state1, state2, kf, kr, kc):           # (0) function call
    """Oscar's catalytic process"""                                 # (1) reaction name
    association_rule = '%s_assoc_%s' % (enz.name, sub.name)         # (2) name of association reaction for rule
    dissociation_rule = '%s_diss_%s' % (enz.name, sub.name)         # (3) name of dissociation reaction for rule
    E = enz(b=None)                                                 # (4) define enzyme state in function
    S = sub({'b': None, site: state1})                              # (5) define substrate state in function
    ES = enz(b=1) % sub({'b': 1, site: state1})
                                                                    # (6) define state of enzyme:substrate complex
    P = sub({'b': None, site: state2})                              # (7) define state of product
    Rule(association_rule, E + S <> ES, kf, kr)                     # (8) rule for enzyme + substrate association (bidirectional)
    Rule(dissociation_rule, ES >> E + P, kc)                        # (9) rule for enzyme:substrate dissociation  (unidirectional)

# instantiate a model
Model()

# declare monomers
Monomer('S', ['b', 'state'],{'state': ['U','P']})
Monomer('E', ['b'])

# input the parameter values
Parameter('kf', 1.1e-06)
Parameter('kr', 1.0e-04)
Parameter('kc', 1.0)

#initial conditions
Parameter('E_0', 1111)
Parameter('S_0', 11111)
Initial(E(b=None), E_0)
Initial(S(b=None, state='U'), S_0)

#observables
#Observable('obsE', enz(b=None))   #Oscar says "doing extra work and defining the species, ex: total enzyme vs. free enzyme species"
#Observable('totalE', E()) # enz()  = enz(b = None) + enz(b=1)% sub(b = 1) or  enz(b = None) + enz(b=1) or enz(b=WILD)<-bound or unbound
#Observable('bound_enzyme', E(b=1)% S(b=1))
#Observable('unbound_E',E(b=None))
#Observable('E1', E(b=1)) #this is the same
Observable('E_Wild', E(b=WILD))
Observable('obsS', S(b = None, state = 'U'))
Observable('obsP', S(b = None, state = 'P'))

# Catalysis
catalyze(E, S, 'state', 'U', 'P', kf, kr, kc)

#visualize




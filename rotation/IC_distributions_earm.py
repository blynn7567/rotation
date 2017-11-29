# two parameters that np needs to process the log distributions that Spencer gives us (the mean and coefficient of variance )
# Now we need the standard deviation and the log mean
# We want to pass different initial conditions to solver:  ScipyOdeSolver(model,tspan,param_values) param_values holds the different conditions!

from __future__ import print_function

import math
import os

import numpy as np
from numpy.random import lognormal



def normal_mu_sigma(log_mean, cv):
    """

    :param log_mean:
    :param cv:
    :return:
    """
    sigma_normal = math.sqrt(math.log((cv ** 2)+1)) # relates coefficient of variance to the standard deviation
    mu_normal = math.log(log_mean)  # - 0.5*(sigma_normal ** 2) To get to the log mean. The log normal disrtibution needs the log mean and the standard deviation
    return mu_normal, sigma_normal


def sample_lognormal(parameter_ic, size, cv=0.25): # For sampling each initial condition of the model; size is how many samples you want to do; cv-default cv for all proteins "normal value for protein variance"
    """

    :param parameter_ic: PySB model parameter
    :param size:
    :param cv:
    :return:
    """
    # mean = np.log(parameter_ic.value)
    cv = cv
    if parameter_ic.name == 'C3_0':
        cv = 0.282
    elif parameter_ic.name == 'XIAP_0' or parameter_ic == 'Bid_0':
        cv = 0.288
    elif parameter_ic.name == 'Bax_0':
        cv = 0.271
    elif parameter_ic.name == 'Bcl2_0':
        cv = 0.294

    mean_normal, mu_normal = normal_mu_sigma(parameter_ic.value, cv)

    return lognormal(mean_normal, mu_normal, size)



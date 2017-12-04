from scipy.optimize import curve_fit


def curve_fit_ftn(functions, observable, xdata, ydata, **kwargs):
    """
    Fit simulation data to specific function

    Parameters
    ----------
    functions: list, must be same length as species
        functions that would be used for fitting the data
    species: list-like, must be same length as functions
        species whose trajectories will be fitted to a function
    xdata: list-like,
        x-axis data points (usually time span of the simulation)
    ydata: list-like,
        y-axis data points (usually concentration of species in time)
    kwargs: dict,
        Key arguments to use in curve-fit

    Returns
    -------
    Parameter values of the functions used to fit the data

    """
    if callable(functions):
        functions = [functions]

    results[sim_array] = curve_fit(functions[sim_array], xdata, ydata[observable], p0=kwargs['p0'])[0]
    return results[0]


def sig_apop(t, f, td, ts):
    """Return the amount of substrate cleaved at time t.

    Keyword arguments:
    t -- time
    f -- is the fraction cleaved at the end of the reaction
    td -- is the delay period between TRAIL addition and half-maximal substrate cleavage
    ts -- is the switching time between initial and complete effector substrate  cleavage
    """
    return f - f / (1 + np.exp((t - td) / (4 * ts)))

import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import minimize

from model import *

linearlabels = ['bk', 'beta', 'gamm', 'rho2', 'tau']

# Linearized model definition.
def linear(x, xlist, params, dt=1e-2):
    # Unpack parameters.
    a, b = x
    bk, beta, gamm, rho, tau = [params[coeff] for coeff in linearlabels]
    astr, bstr, _ = roots( params )

    # Unpack data and delay term.
    itau = round( tau/dt )
    aT = 0 if len( xlist ) < itau else xlist[-itau][0]

    # Compute linear rate of change.
    dx = np.array( [
        bk*beta*(bstr*a + astr*b) - (rho + 2*bk*gamm*astr)*a,
        (rho + 2*bk*gamm*astr)*aT - bk*beta*(bstr*a + astr*b),
    ] )

    # Return next state.
    return x + dx*dt

# Coefficients to characteristic equation of linearized model.
def charactercoeff(params):
    # Unpack parameters.
    bk, beta, gamm, rho, tau = [params[coeff] for coeff in linearlabels]
    astr, bstr, _ = roots( params )

    # Coefficient calculations and return.
    c1 = rho + 2*bk*gamm*astr + bk*beta*(astr - bstr)
    c2 = bk*beta*astr*(rho + 2*bk*gamm*astr)
    return c1, c2

# Constraint on coefficients (negative values are invalid).
def characterconstraint(c1, c2, params=None):
    if params is not None:
        c1, c2 = charactercoeff( params )
    return 2*c2 - c1**2

# Characteristic frequency variable.
def characterfreq(params):
    c1, c2 = charactercoeff( params )
    w2 = characterconstraint( c1, c2 )
    if w2 < 0:
        return np.nan, np.nan, np.nan
    return np.sqrt( w2 ), c1, c2

# Phase equation values.
def phaseconstraint(tau, w, c1, c2, offset=2*np.pi):
    p1 = -w*tau + offset
    p2 = np.arctan2( w*c1, (c2 - w**2) )
    return p1, p2

# Objective function associated with delay length.
def bifurcationcost(theta, label, params, offset=2*np.pi, verbose=False):
    th = theta[0]
    cparams = params.copy()
    cparams[label] = th
    if label == 'N':
        cparams['bk'] = params['k0']*th**(params['alpha'] - 1)

    # Compute characteristic parameters.
    w, c1, c2 = characterfreq( cparams )
    if w is None:
        return np.inf

    # Compute phase parameters.
    p1, p2 = phaseconstraint( cparams['tau'], w, c1, c2, offset=offset )
    return (p1 - p2)**2

# Solve for the critical delay implicitly.
def computebifurcation(thetarange, params, label='tau', kmin=-1, kmax=3, nslc=10):
    thmin, thmax = thetarange
    offsetlist = 2*np.pi*np.arange( kmin, kmax+1 )
    theta0list = np.linspace( thmin, thmax, nslc )

    resultsdata = np.array( [[
        minimize( bifurcationcost, theta, args=(label,params,offset), bounds=[(0,None)], method='Nelder-Mead' )
            for offset in offsetlist]
            for theta in theta0list] )

    errdata = np.array( [[results.success for results in resultslist] for resultslist in resultsdata] )
    fundata = np.array( [[results.fun     for results in resultslist] for resultslist in resultsdata] )
    taudata = np.array( [[results.x[0]    for results in resultslist] for resultslist in resultsdata] )


    succdata = errdata & (fundata < 1e-9)
    icrit = np.argmin( taudata[succdata] ) if np.sum( succdata ) > 0 else None

    return (label, resultsdata[0][0] if icrit is None else resultsdata[succdata][icrit])

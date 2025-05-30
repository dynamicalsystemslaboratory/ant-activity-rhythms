
import sys
from os.path import expanduser

homefolder = expanduser('~')
sys.path.insert( 0, homefolder+'/prog/four' )

import numpy as np
from FOUR.RealFourier import *

from scipy.signal import find_peaks

from model import *

# Define model metrics.
def averageactivity(alist, cut):
    return np.mean( alist[-cut:] )

def valleyminimum(alist, cut):
    return np.min( alist[-cut:] )

def peakmaximum(alist, cut):
    return np.max( alist[-cut:] )

def peakvalleydiff(alist, cut):
    amax = np.max( alist[-cut:] )
    amin = np.min( alist[-cut:] )
    return np.abs( amax - amin )

def cycleperiod(alist, cut, dt=1e-2):
    acut = alist[-cut:]
    peaklist = find_peaks( acut, prominence=1/4*(np.max( acut ) - np.min( acut )) )[0]
    if len( peaklist ) < 1:
        return -1
    return dt*np.mean( np.diff( peaklist ) )

def computemetrics(params, T=None, cut=None, dt=1e-2, computeperd=True):
    # Default simulation time.
    T = 100*params['tau'] if T is None else T
    cut = round( 10*params['tau']/dt ) if cut is None else cut

    # Simulate model with parameter choice and extract activity.
    tlist, xlist = simulate( params, T, [1/2,1/2,0], dt=dt )
    alist = xlist.T[0]

    # Compute and return metrics.
    aperd = cycleperiod( alist, cut, dt=dt ) if computeperd else None
    amean = averageactivity( alist, cut )
    apeak = peakmaximum( alist, cut )
    avall = valleyminimum( alist, cut )
    return aperd, amean, apeak, avall

import sys
from os.path import expanduser

from joblib import Parallel, delayed

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
from cycler import cycler

from scipy.signal import find_peaks
from scipy.optimize import minimize
import statsmodels.api as sm

from model import *
from linear import *

from dataproc import *
from labels import *

from stats import *

from FOUR.RealFourier import *

# Plot font.
default_cycler = cycler(
    color=[
        'indianred',        #CD5C5C
        'cornflowerblue',   #6699ff
        'olivedrab',        #6B8E23
        'mediumpurple',     #9370DB
        'sandybrown',       #F4A460
        'steelblue',        #4682B4
    ] )
plt.rcParams.update( {
    'axes.prop_cycle': default_cycler,
    'text.usetex': True,
    'font.family': 'mathptmx',
    'font.size': 14,
    'text.latex.preamble': "\\usepackage{amsmath}",
} )

# Set global number print setting.
np.set_printoptions(precision=3, suppress=False, threshold=np.inf, linewidth=np.inf)

# Fit scaling laws to data.
def fitscalinglaw(x, y):
    return sm.OLS( np.log( y ), sm.add_constant( np.log( x ) ) ).fit()

# Helper function for adding significance annotation.
def addsignificance(axs, t, p, base, offset, tick):
    if p < 0.05:
        label = ('*' if results.pvalue < 0.001 else '%.3f' % results.pvalue) \
            + (' (+)' if t > 0 else (' (-)'))
        axs.plot( [0, i], [base + offset, base + offset], color='k' )
        axs.plot( [0, 0], [base + offset, base + offset - tick], color='k' )
        axs.plot( [i, i], [base + offset, base + offset - tick], color='k' )
        axs.text( i/2, base + offset + tick, label, horizontalalignment='center' )
        return True
    return False

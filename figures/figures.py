
import sys
from os.path import expanduser

sys.path.insert( 0, './..' )

from args import *

# Define data metrics.
def activity(adata):
    return np.mean( adata, axis=1 )

# Average peak height metric.
def getpeakpoints(adata):
    return [find_peaks( alist, prominence=1/4*(np.max( alist ) - np.min( alist )) )[0] for alist in adata]
def peakheight(adata):
    pdata = getpeakpoints( adata )
    return np.array( [np.mean( alist[plist] ) for alist, plist in zip( adata, pdata )] )

# Average valley depth metric.
def getvallpoints(adata):
    return [find_peaks( -alist, prominence=1/4*(np.max( alist ) - np.min( alist )) )[0] for alist in adata]
def valleydepth(adata):
    vdata = getvallpoints( adata )
    return np.array( [np.mean( alist[vlist] ) for alist, vlist in zip( adata, vdata )] )

# Average peak-valley difference metric.
def peakvalleydifference(adata):
    pdata = getpeakpoints( adata )
    vdata = getvallpoints( adata )
    return np.array( [np.mean( [alist[p] - alist[v] for p, v in zip( plist, vlist )] )
        for alist, plist, vlist in zip( adata, pdata, vdata )] )

# Dominant period metric.
def getdomperiod(alist, verbose=False):
    tlist = 1/2*np.arange( len( alist ) )
    fvar = RealFourier( tlist[None], alist[None] - np.mean( alist ) ).dft()
    period = realmaximum( fvar ).period[0][0]
    if verbose:
        return period, fvar
    return period
def dominantperiod(adata):
    return np.array( [getdomperiod( alist ) for alist in adata] )

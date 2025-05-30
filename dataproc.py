
import sys
from os.path import expanduser

import numpy as np

homefolder = expanduser('~')
figurefolder = '/Documents/papers/rsc-delay/figures/'

from labels import *

# Helper function for saving figures.
def savefigure(fig, filename):
    fig.savefig( homefolder + figurefolder + filename,
        dpi=600, bbox_inches='tight', transparent=1 )

# Generally saving data.
def savedata(filename, data, ext='.txt'):
    if ext not in filename:
        filename += ext
    np.savetxt( filename, data, delimiter=',' )
    return filename

def loaddata(filename, ext='.txt'):
    if ext not in filename:
        filename += ext
    data = np.loadtxt( filename, delimiter=',' )
    return data

# Function to unpack simulation data.
def unpackrealdata(filelist):
    datalist = []
    adata = []
    for file in filelist:
        if '.csv' not in file:
            file += '.csv'
        data = np.loadtxt( file, delimiter=',' ).T
        adata += [data]
    return adata

# Simple moving average function.
def sma(T, X, W, dT=0.1, giveT=1):
    # Extract bounds of the new time-series.
    MT = T.shape[0];  MX = X.shape[0]
    Tmin = np.min( T );  Tmax = np.max( T )
    NT = round( (Tmax - Tmin)/dT ) + 1

    # Initialize new time-series.
    Tm = np.empty( (MT, NT) )
    Xm = np.empty( (MX, NT) )

    # Calculate SMA over the data set at selected points.
    for k in range( NT ):
        Tm[:,k] = k*dT + Tmin
        Xm[:,k] = np.mean( X[((Tm[:,k] - W) < T) & (T < (Tm[:,k] + W))] )

    # Return the equally spaced, averaged data.
    if giveT:
        return Tm, Xm
    return Xm

# Import colony-level data.
def importcolonydata(colony='EGTi'):
    # Get data file names.
    colonyfilelist = [homefolder + '/prog/ant-active/' + filename for filename in datafiles[colony]]
    colonynum = np.arange( len( colonyfilelist ) )

    # Unpack experimental data.
    dt = 1/2
    adataraw = np.array( unpackrealdata( colonyfilelist ) )

    # List of colony population sizes.
    Nlist = np.array( [Ndict[filename] for filename in datafiles[colony]] )
    Nsort = np.unique( Nlist )

    # Smooth data.
    tlist = dt*np.arange( adataraw.shape[1] )
    adata = np.array( [sma( tlist[None], alist[None], 2, dT=dt, giveT=0 )[0]
        for N, alist in zip( Nlist, adataraw )] )

    return colonyfilelist, (colonynum, Nlist), (dt, tlist, adata)

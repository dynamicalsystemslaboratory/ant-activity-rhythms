# File: Transforms.py
# Created by: Michael Napoli
# Created on: Jul 11, 2023
# Purpose: To develop and study various methods of expanding
#   and solving real/complex Fourier series.

import numpy as np
from copy import deepcopy
from FOUR.Preprocess import *

# Class: CharacteristicWave()
# Purpose: To be used to save the characteristic wave form found through the Transform() class.
class CharacteristicWave:
    def __init__(self, ampl=None, freq=None, phase=None, wave_type='sin'):
        # Initialize class variables.
        self.ampl = ampl
        self.freq = freq
        self.phase = phase
        self.wave_type = wave_type

    @property
    def period(self):
        return 2*np.pi/self.freq

    def __str__(self):
        line1 = 'Characteristic wave: '
        line2 = 'x(t) = %.3f ' % self.ampl \
            + self.wave_type \
            + '( %.3f t + %.3f )' % (self.freq, self.phase)
        return line1 + line2

    def setCharacteristics(self, ampl=None, freq=None, phase=None):
        # Update components if given values are not None.
        self.ampl = self.ampl if ampl is None else ampl
        self.freq = self.freq if freq is None else freq
        self.phase = self.phase if phase is None else phase

        # Return instance of self.
        return self

    def solve(self, T):
        assert None not in (self.freq, self.phase, self.ampl), \
            'ERROR: One (or all) of CharacteristicWave.freq/phase/ampl is not set.'
        # Return solution of wave form.
        return self.ampl*np.sin( self.freq*T + self.phase )

# Class: Transform()
class Transform:
    def __init__(self, T, X, N=None, dt=None):
        # Initialize transform variables and dimensions.
        self.T = T      # Input data.
        self.X = X      # Output data.

        self.F = None   # List of frequencies.
        self.P = None   # Power spectrum list split by coefficient members.
        self.R = None   # Power spectrum.
        self.Rmax = None    # Max value from power spectrum.

        # Dimensions and frequency parameter setup for the transform.
        self.K = T.shape[0]
        self.N = round( T.shape[1]/2 ) if N is None else N
        self.Nt = X.shape[0]
        self.dt = self.T[0,1] - self.T[0,0] if dt is None else dt
        self.tau = 2*self.N*self.dt
        self.err = None
        self.sort = np.array( [np.arange( 0, self.N+1 ) for _ in range( self.K )] )

        # For user meta-data.
        self.meta = None

    def setDataSets(self, T, X):
        self.__init__( T, X )
        # Return instance of self.
        return self

    def setCoefficientNumber(self, N):
        self.N = N
        # Return instance of self.
        return self

    def frequencies(self):
        # Generate frequency list.
        self.F = 2*np.pi/self.tau*np.arange( self.K*(self.N + 1) )[:,None]

        # Return instance of self.
        return self

    def resError(self, T=None, X=None, save=0):
        self.check

        # Initialize data matrix (default: training data).
        T = self.T if T is None else T
        X = self.X if X is None else X

        # Solve for approximation of set.
        Y = self.solve( T )

        # Calculate residual error.
        err = np.linalg.norm( X - Y )**2

        # Save if requested and return.
        if save:
            self.err = err
        return err

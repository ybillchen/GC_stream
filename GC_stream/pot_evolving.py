"""
BSD 3-Clause License
Copyright (c) 2024 Yingtian Chen
All rights reserved.
"""

import numpy as np
from scipy.interpolate import interp1d

import agama
agama.setUnits(mass=1, length=1, velocity=1) # Msun, kpc, km/s

gravG = 4.30091727003628e-6 # in defult unit
ln10 = 2.302585092994046

def dlnMdlnr_from_single(Menclose_interpolator, radius):
    lnr = np.log(radius)
    dlnr = 0.1
    lnM1 = np.log(Menclose_interpolator((lnr - 0.5*dlnr)/ln10))
    lnM2 = np.log(Menclose_interpolator((lnr + 0.5*dlnr)/ln10))
    dlnM = lnM2 - lnM1
    return dlnM/dlnr

def get_Menclose_interpolator(pot, logrmin=-1., logrmax=3., N=100):
    logrs = np.linspace(logrmin, logrmax, N+1)
    Mencloses = []
    for logr in logrs:
        Mencloses.append(pot.enclosedMass(10**logr))

    return interp1d(logrs, Mencloses)

class Evolving(object):
    """Evolving potential"""

    def __init__(self, pots, times):
        """
        Note: times must be ascending
        """
        super(Evolving, self).__init__()
        self.pots = pots
        self.times = times
        self.Menclose_interpolators = []
        for pot in pots:
            self.Menclose_interpolators.append(
                get_Menclose_interpolator(pot))

    def dlnMdlnr(self, radius, time):
        logr = np.log10(radius)
        if time < self.times[0]:
            return dlnMdlnr_from_single(
                self.Menclose_interpolators[0], logr)
        if time >= self.times[-1]:
            return dlnMdlnr_from_single(
                self.Menclose_interpolators[-1], logr)
        
        for i in range(len(self.times)-1):
            if time >= self.times[i+1]:
                continue
            else:
                r = ((time - self.times[i]) / 
                    (self.times[i+1] - self.times[i]))
                out1 = dlnMdlnr_from_single(
                    self.Menclose_interpolators[i], logr)
                out2 = dlnMdlnr_from_single(
                    self.Menclose_interpolators[i+1], logr)
                return out1*(1-r) + out2*r

    def Menclose(self, radius, time):
        logr = np.log10(radius)
        if time < self.times[0]:
            return self.Menclose_interpolators[0](logr)
        if time >= self.times[-1]:
            return self.Menclose_interpolators[-1](logr)
        
        for i in range(len(self.times)-1):
            if time >= self.times[i+1]:
                continue
            else:
                r = ((time - self.times[i]) / 
                    (self.times[i+1] - self.times[i]))
                out1 = self.Menclose_interpolators[i](logr)
                out2 = self.Menclose_interpolators[i+1](logr)
                return out1*(1-r) + out2*r

    def Vcirc(self, radius, time):
        return np.sqrt(gravG * self.Menclose(radius, time) / radius)

class M17Pot(Evolving):
    """McMillan el al. (2017) potential"""

    def __init__(self):

        disk1 = agama.Potential(type='Disk', 
            surfaceDensity=8.95679e+08, scaleRadius=2.49955, 
            scaleHeight=0.3, innerCutoffRadius=0, sersicIndex=1)
        disk2 = agama.Potential(type='Disk', 
            surfaceDensity=1.83444e+08, scaleRadius=3.02134, 
            scaleHeight=0.9, innerCutoffRadius=0, sersicIndex=1)
        disk3 = agama.Potential(type='Disk', 
            surfaceDensity=5.31319e+07, scaleRadius=7, 
            scaleHeight=-0.085, innerCutoffRadius=4, sersicIndex=1)
        disk4 = agama.Potential(type='Disk', 
            surfaceDensity=2.17995e+09, scaleRadius=1.5, 
            scaleHeight=-0.045, innerCutoffRadius=12, sersicIndex=1)
        halo1 = agama.Potential(type='Spheroid', 
            densityNorm=9.8351e+10, axisRatioZ=0.5, 
            gamma=0, beta=1.8, scaleRadius=0.075, 
            outerCutoffRadius=2.1, alpha=1, axisRatioY=1, 
            cutoffStrongth=2)
        halo2 = agama.Potential(type='Spheroid', 
            densityNorm=8.53702e+06, axisRatioZ=1, 
            gamma=1, beta=3, scaleRadius=19.5725, 
            outerCutoffRadius=1e6, alpha=1, axisRatioY=1, 
            cutoffStrongth=1e6)
        pot = agama.Potential(
            disk1, disk2, disk3, disk4, halo1, halo2)

        super(M17Pot, self).__init__([pot], [0])
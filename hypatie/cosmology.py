"""
Module cosmology
===============
This module has a parent classes (CosModel) to represent a Flat Lambda-CDM
model of the universe.
There are also child classes (suck as Plank18) who provide initial data for
the model.
"""

import numpy as np
import scipy.integrate as integrate

class CosModel:
    """
    Cosmological model by density parameters (Omegas)

    Parameters
    ----------
        H0 (float)  : Current Hubble constant
        Oph0 (float): Current density parameter of photon
        One0 (float): Current density parameter of neutrino
        Obm0 (float): Current density parameter of baryonic matter
        Odm0 (float): Current density parameter of dark matter
        Ode0 (float): Current density parameter of dark energy

    Attributes
    ----------
        sf (np.poly1d): scale factor
        
    Methods
    -------
        scale_factor: Scale factor
        current_age: Current age of the universe
        proper_distance: Proper distance between two moments, observed & emitted
        z: The redshift
    """
    def __init__(self, H0, Oph0, One0, Obm0, Odm0, Ode0):
        
        self.H0 = H0
        self.Oph0 = Oph0
        self.One0 = One0
        self.Obm0 = Obm0
        self.Odm0 = Odm0
        self.Ode0 = Ode0
        self.sf = None

    def initial_calculations(self):
        # Convert unit of H0 from ('km / (Mpc s)') to ('1 / Gyr')
        self.H0 = self.H0 * 0.00102271
        # Density parameters
        self.Or0 = self.Oph0 + self.One0
        self.Om0 = self.Obm0 + self.Odm0
        self.Oall = self.Or0 + self.Om0 + self.Ode0

    def scale_factor(self, a_min=1e-10, a_max=5, n_points=3000, fit_deg=9):
        """
        Scale factor of the universe in terms of time.
        
        Arguments
        ----------
            a_min (float): maximum value for scale factor
            a_max (float): minimum value for scale factor
            n_points (int): number of points used to fit
            fit_deg (int): degree of the polynomiat fit
        
        Returns
        -------
            Scale factor as a np.poly1d object
        """
        if self.sf is None:
            self.initial_calculations()
            a = np.linspace(a_min, a_max, n_points)
            tmp = lambda a: 1 / np.sqrt(
                self.Or0/a**2 + self.Om0/a + self.Ode0*a**2 + 1-self.Oall
                )
            t = np.array([integrate.quad(tmp, 0, i)[0]/self.H0 for i in a])
            coefs = np.polyfit(t, a, fit_deg)
            self.sf = np.poly1d(coefs)
        return self.sf

    def current_age(self):
        """
        Current age of the universe

        Returns
        -------
            Current age as Billion years from the Big bang
        """
        if self.sf is None:
            self.sf = self.scale_factor()
        r = np.roots(self.sf-1)
        return r[np.isreal(r)].real[0]

    def proper_distance(self, to, te):
        """
        Proper distance between two moments, observed (to) and emitted (te)
        
        Arguments
        ----------
            to (float): emitted moment (in Billion years from Big bang)
            te (float): observed moment (in Billion years from Big bang)
        
        Returns
        -------
            Proper distance in Billion light-years
        """
        if self.sf is None:
            self.sf = self.scale_factor()
        dp = integrate.quad(lambda t: 1/self.sf(t), te, to)[0]
        return dp

    def z(self, t):
        """
        Redshift (z value) at a given moment
        
        Arguments
        ----------
            t (float): moment (in Billion years from Big bang)
        
        Returns
        -------
            z value as float
        """
        if self.sf is None:
            self.sf = self.scale_factor()
        return (1/self.sf(t))-1


class Planck18(CosModel):
    """
    Predefined CosModel class with density parameters provided by Planck 2018
    """
    def __init__(self):
        self.H0 = 67.66
        self.Oph0 = 5.402015137139351e-05
        self.One0 = 0.0014396743040845238
        self.Obm0 = 0.04897
        self.Odm0 = 0.26069
        self.Ode0 = 0.6888463055445441
        self.sf = np.poly1d([ 3.14416422e-17, -8.77434127e-15, 1.09165168e-12,
            -7.98435601e-11, 3.81160114e-09, -1.24683679e-07, 2.85512434e-06,
            -4.59502352e-05, 5.14673404e-04, -3.92192424e-03, 1.96940133e-02,
            -6.30581943e-02, 1.85807439e-01,  1.54295181e-02])

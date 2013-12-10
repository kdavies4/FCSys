#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Load, analyze, and plot results from Modelica_ simulations.

This module supports the approach to physical units and quantities established
in FCSys.  Please see
`FCSys.Units <http://kdavies4.github.io/FCSys/FCSys_Units.html>`_ for more
information.

.. _Modelica: http://www.modelica.org/
"""
__author__ = "Kevin Davies"
__email__ = "kdavies4@gmail.com"
__copyright__ = "Copyright 2012-2013, Georgia Tech Research Corporation"
__license__ = "BSD-compatible (see LICENSE.txt)"


import re
import modelicares.base as base

from modelicares import SimRes, label_number

# Establish the default units.
import _config
config = _config.config
default_units = config['Default display units']

SI_prefixes = dict(Y=1e24, # yotta
                   Z=1e21, # zetta
                   E=1e18, # exa
                   P=1e15, # peta
                   T=1e12, # tera
                   G=1e9, # giga
                   M=1e6, # mega
                   k=1e3, # kilo
                   h=100, # hecto
                   da=10, # deca
                   d=0.1, # deci
                   c=0.01, # centi
                   m=1e-3, # milli
                   u=1e-6, # micro
                   n=1e-9, # nano
                   p=1e-12, # pico
                   f=1e-15, # femto
                   a=1e-18, # atto
                   z=1e-21, # zepto
                   y=1e-24) # yocto


class SimRes(SimRes):
    """Base class for Modelica_-based simulation results and methods to analyze
    those results

    .. :Note: This version of :class:`SimRes` is different than that in
       modelicares
       `http://kdavies4.github.com/modelicares/
       <http://kdavies4.github.com/modelicares/>`_.
       The Modelica_ *unit* attribute is used to indicate the dimension of a
       variable according to the method by Davies and Paredis 2012 (`ref
       <https://modelica.org/events/modelica2012/Proceedings>`_).
    """

    # TODO: Update sankey() and bar() from modelicares.SimRes.  See sankey.py
    # in this folder.

    def __init__(self, fname="dsres.mat", constants_only=False):
        """On initialization, load Modelica_ simulation results from a
        MATLAB\ :sup:`®` file in Dymola\ :sup:`®` format.

        **Arguments:**

        - *fname*: Name of the file (may include the path)

             The file extension ('.mat') is optional.  *fname* should have an
             record of base units and constants named *environment.baseUnits*.

        - *constants_only*: *True*, if only the variables from the first data
          table should be loaded

             The first data table typically contains all of the constants,
             parameters, and variables that don't vary.  If only that
             information is needed, it will save some time and memory to set
             *constants_only* to *True*.

        **Example:**

           >>> from modelicares import SimRes
           >>> sim = SimRes('examples/SaturationPressure.mat')
        """
        super(SimRes, self).__init__(fname, constants_only)
        self._set_constants('environment.baseUnits')

        # Functions for units that involve offsets or more complex mappings
        # See "Conversions with offsets" in FCSys/Units.mos.
        self._to_unit = {
            'degC': lambda x: x/self._unitvalue('K') - 273.15,
            'degF': lambda x: (9/5.0)*(x/self._unitvalue('K') - 273.15) + 32,
            'kPag': lambda x: x/self._unitvalue('kPa') - 101.325,
            'Pag': lambda x: x/self._unitvalue('Pa') - 101325,
            }
        # TODO: Add special scaling for the der() function

    def get_unit(self, names):
        """Return the unit(s) for trajectory variable(s).

        The unit is the variable's *displayUnit* attribute unless it is empty
        ('').  In that case, the unit is taken from the global *default_unit*
        dictionary based on the variable's dimension.

        **Arguments:**

        - *names*: Name(s) of the variable(s) from which to get the unit(s)

             This may be a single string or (possibly nested) list of strings
             representing the names of the variables.

        If *names* is a string, then the output will be a single display unit.
        If *names* is a (optionally nested) list of strings, then the output
        will be a (nested) list of units.

        .. Note:: This method (:meth:`fcres.SimRes.get_unit`) is different
           than :meth:`modelicares.SimRes.get_unit`.  In FCSys_, the *unit*
           attribute is used to indicate the **dimension** of a quantity.  The
           **unit** is given by the variable's *displayUnit* attribute or the
           default unit based on the variable's dimension, as described above.

        **Example:**

           >>> from fcres import SimRes

           >>> sim = SimRes('examples/SaturationPressure.mat')
           >>> sim.get_unit('environment.p')
           'kPa'
           >>> sim.get_unit([[['environment.p', 'environment.T']]])
           [[['kPa', 'degC']]]
        """
        def _get_unit(name):
            """Return the unit for a single variable.
            """
            if name == 'Time':
                # Time is special since it's predefined in Modelica.
                return 's'
            # TODO: Need to handle scaling properly for time (don't scale)

            unit = self._traj[name].displayUnit
            if unit == '':
                unit = self._traj[name].unit
                if unit<>'1':
                    unit = default_units[unit]
            return unit

        return self._get(names, _get_unit)

    def get_displayUnit(self):
        raise AttributeError("'SimRes' object has no attribute "
                             "'get_displayUnit'.  Use 'get_unit' instead." )

    def _unitvalue(self, unit):
        """Convert a Modelica_ unit string to a numeric value using the values
        of base the constants and units.

        For example, "m/s2" becomes m/s^2, where m and s are the values of keys
        'm' and 's' in the internal *U* dictionary.

        **Arguments:**

        - *unit*: Unit string in Modelica_ notation

             .. Seealso:: Modelica Specification, version 3.3, p. 235--236
                (https://www.modelica.org/documents)

                In summary, '.' indicates multiplication.  The denominator is
                enclosed in parentheses and begins with '/'.  Exponents
                directly follow the significand (e.g., no carat ('^')).
        """
        splitter = re.compile('([^0-9+-]*)(.*)')

        def _process_unit(unit):
            """Convert a simple Modelica_ unit to LaTeX.
            """
            if unit.isdigit():
                return 1
            base, exponent = splitter.match(unit).groups()
            try:
                base = self.U[base]
            except KeyError:
                try:
                    base = SI_prefixes[base[0]]*self.U[base[1:]]
                except KeyError:
                    print("The unit group %s is not recognized.  "
                          "It will be skipped." % unit)
            return base**int(exponent) if exponent else base

        def _process_group(unit):
            """Convert the numerator or denominator of a unit string to a
            number.
            """
            if unit.startswith('('):
                assert unit.endswith(')'), ("The unit group %s starts with '(' but "
                                            "does not end with ')'." % unit)
                unit = unit[1:-1]
            val = 1
            for u in unit.split('.'):
                val *= _process_unit(u)
            return val

        if unit:
            # Split the numerator and the denominator.
            if '/' in unit:
                try:
                    numer, denom = unit.split('/')
                except ValueError:
                    print("Check that the unit string %s has at most one "
                          "division sign." % unit)
                return _process_group(numer)/_process_group(denom)
            else:
                return _process_group(unit)
        return 1

    def to_unit(self, unit):
        """Return a function to display a quantity in a unit.

        **Arguments:**

        - *unit*: Unit as a string in Modelica_ notation

             .. Seealso:: Modelica Specification, version 3.3, p. 235--236
                (https://www.modelica.org/documents)

                In summary, '.' indicates multiplication.  The denominator is
                enclosed in parentheses and begins with '/'.  Exponents
                directly follow the significand (e.g., no carat ('^')).

        **Example:**

           >>> from fcres import SimRes
           >>> sim = SimRes('examples/SaturationPressure.mat')
           >>> to_kPa = sim.to_unit('kPa')

           >>> # Method 1:
           >>> p = sim.get_IV('environment.p')
           >>> print("The pressure is %.1f kPa." % to_kPa(p))
           The pressure is 101.3 kPa.

           >>> # Method 2:
           >>> p_kPa = sim.get_IV('environment.p', to_kPa)
           >>> print("The pressure is %.1f kPa." % p_kPa)
           The pressure is 101.3 kPa.

           >>> # In other units:
           >>> print("The pressure is %.1f kPag." % sim.to_unit('kPag')(p))
           The pressure is 0.0 kPag.
           >>> print("The pressure is %.0f Pa." % sim.to_unit('Pa')(p))
           The pressure is 101325 Pa.
           >>> print("The pressure is %.3f atm." % sim.to_unit('atm')(p))
           The pressure is 1.000 atm.
           >>> print("The pressure is %.3f bar." % sim.to_unit('bar')(p))
           The pressure is 1.013 bar.
        """
        try:
            # Assume the unit involves an offset or a more complex function.
            return self._to_unit[unit]
        except KeyError:
            # The unit is just a scaling factor.
            return lambda x: x/self._unitvalue(unit)

    def plot(self, ynames1=[], ylabel1=None, yunit1=None, legends1=[],
             leg1_kwargs={'loc': 'best'}, ax1=None,
             ynames2=[], ylabel2=None, yunit2=None, legends2=[],
             leg2_kwargs={'loc': 'best'}, ax2=None,
             xname='Time', xlabel=None, xunit=None,
             title=None, label="xy", incl_prefix=False, suffix=None,
             use_paren=True, **kwargs):
        r"""Plot data as points and/or curves in 2D Cartesian coordinates.

        A new figure is created if necessary.

        **Arguments:**

        - *ynames1*: Names of variables for the primary y axis

             If any names are invalid, then they will be skipped.

        - *ylabel1*: Label for the primary y axis

             If *ylabel1* is *None* (default) and all of the variables have the
             same Modelica_ description string, then the common description
             will be used.  Use '' for no label.

        - *yunit1*: String indicating the unit for the primary y-axis (see note
          for *xunit*)

             If *yunit1* is *None*, the Modelica_ *displayUnit* of the first
             entry of *ynames1* or the default unit (from *fcres.ini*) based on
             that variable's dimension will be used (in decreasing priority).

             .. Note:: Dimension checking is not currently performed, so it is
                important to ensure that a proper unit is chosen.

        - *legends1*: List of legend entries for variables assigned to the
          primary y axis

             If *legends1* is an empty list ([]), ynames1 will be used.  If
             *legends1* is *None* and all of the variables on the primary axis
             have the same unit, then no legend will be shown.

        - *leg1_kwargs*: Dictionary of keyword arguments for the primary legend

        - *ax1*: Primary y axes

             If *ax1* is not provided, then axes will be created in a new
             figure.

        - *ynames2*, *ylabel2*, *yunit2*, *legends2*, *leg2_kwargs*, and *ax2*:
          Similar to *ynames1*, *ylabel1*, *yunit1*, *legends1*, *leg1_kwargs*,
          and *ax1* but for the secondary y axis

        - *xname*: Name of the x-axis data

        - *xlabel*: Label for the x axis

             If *xlabel* is *None* (default), the variable's Modelica_
             description string will be applied.  Use '' for no label.

        - *xunit*: String indicating the unit for the x axis  (see note
          for *yunit1*)

             If *xunit* is *None*, the Modelica_ variable's *displayUnit* or
             the default unit (from *fcres.ini*) based on the variable's
             dimension will be used (in decreasing priority).

        - *title*: Title for the figure

             If *title* is *None* (default), then the title will be the base
             filename.  Use '' for no title.

        - *label*: Label for the figure (ignored if ax is provided)
             This will be used as a base filename if the figure is saved.

        - *incl_prefix*: If *True*, prefix the legend strings with the base
          filename of the class.

        - *suffix*: String that will be added at the end of the legend entries

        - *use_paren*: Add parentheses around the suffix

        - *\*\*kwargs*: Additional arguments for  :meth:`base.plot` (and thus to
          :meth:`matplotlib.pyplot.plot`)

             If both y axes are used (primary and secondary), then the *dashes*
             argument is ignored.  The curves on the primary axis will be solid
             and the curves on the secondary axis will be dotted.

        **Returns:**

        1. *ax1*: Primary y axes

        2. *ax2*: Secondary y axes

        **Example:**

        .. testsetup::
           >>> from fcres import closeall
           >>> closeall()

        .. code-block:: python

           >>> from fcres import SimRes, saveall

           >>> sim = SimRes('examples/SaturationPressure')
           >>> sim.plot(xname="subregion.gas.H2O.T",
           ...          ynames1=["subregion.gas.H2O.p", "p_sat"],
           ...          legends1=["FCSys (from Gibbs equilibrium)",
           ...                    "Modelica.Media (correlated function)"],
           ...          ylabel1='Saturation pressure',
           ...          title="Water Saturation Pressure",
           ...          label='examples/SaturationPressure') # doctest: +ELLIPSIS
           (<matplotlib.axes._subplots.AxesSubplot object at 0x...>, None)

           >>> saveall()
           Saved examples/SaturationPressure.pdf
           Saved examples/SaturationPressure.png

        .. only:: html

           .. image:: ../examples/SaturationPressure.png
              :scale: 70 %
              :alt: plot of water saturation pressure

        .. only:: latex

           .. figure:: ../examples/SaturationPressure.pdf
              :scale: 70 %

              Plot of water saturation pressure
        """
        # Note:  ynames1 is the first argument (besides self) so that plot()
        # can be called with simply a variable name.

        def _ystrings(ynames, ylabel, yunit, legends):
            """Generate a y-axis label and set of legend entries.
            """
            if ynames:
                if ylabel is None: # Try to create a suitable axis label.
                    descriptions = self.get_description(ynames)
                    # If the descriptions are the same, label the y axis with
                    # the 1st one.
                    ylabel = descriptions[0]
                    if len(set(descriptions)) <> 1:
                        print("The y-axis variable descriptions are not all "
                              "the same.  The first has been used.  Please "
                              "provide the proper name via ylabel1 or ylabel2.")
                if legends == []:
                    legends = ynames
                if incl_prefix:
                    legends = [self.fbase + ': ' + leg for leg in legends]
                if suffix:
                    legends = ([leg + ' (%s)' % suffix for leg in legends]
                               if use_paren else
                               [leg + suffix for leg in legends])
                assert len(set(self.get_dimension(ynames))) == 1, \
                    "The variables on the y-axis do not have the same physical dimension."
                if yunit is None:
                    # Use the unit of the 1st variable.
                    yunit = self.get_unit(ynames[0])
                if ylabel <> "":
                    ylabel = label_number(ylabel, yunit)

            return ylabel, yunit, legends

        # Process the inputs.
        ynames1 = base.flatten_list(ynames1)
        ynames2 = base.flatten_list(ynames2)
        assert ynames1 or ynames2, "No signals were provided."
        if title is None:
            title = self.fbase

        # Create primary and secondary axes if necessary.
        if not ax1:
            fig = base.figure(label)
            ax1 = fig.add_subplot(111)
        if ynames2 and not ax2:
            ax2 = ax1.twinx()

        # Generate the x-axis label.
        if xlabel is None:
            xlabel = 'Time' if xname == 'Time' else self.get_description(xname)
            # With Dymola 7.4, the description of the time variable will be
            # "Time in", which isn't good.
        if xunit is None:
            xunit = self.get_unit(xname)
        if xlabel<> "":
            xlabel = label_number(xlabel, xunit)

        # Generate the y-axis labels and sets of legend entries.
        ylabel1, yunit1, legends1 = _ystrings(ynames1, ylabel1, yunit1, legends1)
        ylabel2, yunit2, legends2 = _ystrings(ynames2, ylabel2, yunit2, legends2)

        # Read the data.
        if xname == 'Time':
            t_scale = lambda t: t*self._unitvalue('s')/self._unitvalue(xunit)
            y_1 = self.get_values(ynames1, f=self.to_unit(yunit1))
            y_2 = self.get_values(ynames2, f=self.to_unit(yunit2))
        else:
            x = self.get_values(xname, f=self.to_unit(xunit))
            times = self.get_times(xname)
            y_1 = self.get_values_at_times(ynames1, times, f=self.to_unit(yunit1))
            y_2 = self.get_values_at_times(ynames2, times, f=self.to_unit(yunit2))

        # Plot the data.
        if ynames1:
            if ynames2:
                # Use solid lines for primary axis and dotted lines for
                # secondary.
                kwargs['dashes'] = [(None, None)]
                base.plot(y_1, self.get_times(ynames1, f=t_scale) if xname == 'Time'
                          else x, ax1, label=legends1, **kwargs)
                kwargs['dashes'] = [(3, 3)]
                base.plot(y_2, self.get_times(ynames2, f=t_scale) if xname == 'Time'
                          else x, ax2, label=legends2, **kwargs)
            else:
                base.plot(y_1, self.get_times(ynames1, f=t_scale) if xname == 'Time'
                          else x, ax1, label=legends1, **kwargs)
        elif ynames2:
            base.plot(y_2, self.get_times(ynames2, f=t_scale) if xname == 'Time'
                      else x, ax2, label=legends2, **kwargs)

        # Decorate the figure.
        ax1.set_title(title)
        ax1.set_xlabel(xlabel)
        if ylabel1:
            ax1.set_ylabel(ylabel1)
        if ylabel2:
            ax2.set_ylabel(ylabel2)
        if legends1:
            if legends2:
                # Put the primary legend in the upper left and secondary in
                # upper right.
                leg1_kwargs['loc'] = 2
                leg2_kwargs['loc'] = 1
                ax1.legend(**leg1_kwargs)
                ax2.legend(**leg2_kwargs)
            else:
                ax1.legend(**leg1_kwargs)
        elif legends2:
            ax2.legend(**leg2_kwargs)

        return ax1, ax2

    def get_dimension(self, names):
        """Return the dimension(s) of trajectory variable(s).

        **Arguments:**

        - *names*: Name(s) of the variable(s) from which to get the
          dimension(s)

             This may be a single string or (possibly nested) list of strings
             representing the names of the variables.

        If *names* is a string, then the output will be a single unit.  If
        *names* is a (optionally nested) list of strings, then the output will
        be a (nested) list of units.

        **Example:**

           >>> from fcres import SimRes

           >>> sim = SimRes('examples/SaturationPressure.mat')
           >>> sim.get_dimension('environment.p')
           'm/(l.T2)'
           >>> sim.get_dimension([[['environment.p', 'environment.T']]])
           [[['m/(l.T2)', 'l2.m/(N.T2)']]]
        """
        return self._get(names, lambda name: self._traj[name].unit)

    def _set_constants(self, record='environment.baseUnits'):
        """Establish the values of constants and units.

        **Arguments:**

        - *record*: Name of the record which contains the base units

        There are no return values.  All the constants and units are dependent
        on the base constants and units which are included in the
        environment.baseUnits record.  See the Units package in FCSys/Units.mo
        and the Environment model in FCSys/Conditions.mo.
        """
        # Mathematical constants
        from math import acos, exp
        pi = 2*acos(0)
        self.U = dict(pi=pi)

        # Base constants and units
        record += '.'
        for u in ['rad', 'R_inf', 'c', 'k_J', 'R_K', "'cd'", 'k_F', 'R']:
            value = self.get_IV(record + u)
            if value is None:
                print("The base constants and units were not loaded from the "
                      "simulation.  The defaults from fcres.ini will be used.")
                constants = config['Default base constants']
                constants["'cd'"] = constants.pop('cd')
                for key in constants.iterkeys():
                    constants[key] = float(constants[key])
                self.U.update(constants)
                break
            else:
                self.U.update({u: value})
        rad = self.U['rad']
        R_inf = self.U['R_inf']
        c = self.U['c']
        k_J = self.U['k_J']
        R_K = self.U['R_K']
        cd = self.U["'cd'"]
        k_F = self.U['k_F']
        R = self.U['R']

        # Empirical units
        m = 10973731.568539*rad/R_inf
        s = 299792458*m/c
        Wb = 483597.870e9/k_J
        S = 25812.8074434/R_K
        mol = 96485.3365*Wb*S/k_F
        K = 8.3144621*(Wb*rad)**2*S/(s*mol*R)
        self.U.update(m=m, s=s, Wb=Wb, S=S, mol=mol, K=K)

        # SI base units [BIPM2006, Table 1] and intermediate units
        V = Wb*rad/s
        A = V*S
        C = A*s
        J = V*C
        Gy = (m/s)**2
        kg = J/Gy
        self.U.update(V=V, A=A, C=C, J=J, Gy=Gy, kg=kg)

        # Coherent derived units in the SI with special names and symbols
        cyc = 2*pi*rad
        Hz = cyc/s
        sr = rad**2
        N = J/m
        Pa = N/m**2
        W = J/s
        F = C/V
        ohm = 1/S
        H = V*s/A
        T = Wb/m**2
        lm = cd*sr
        lx = lm/m**2
        Bq = Hz
        Sv = Gy
        kat = mol/s
        g = 0.001*kg
        self.U.update(cyc=cyc, Hz=Hz, sr=sr, N=N, Pa=Pa, W=W, F=F, ohm=ohm,
                      H=H,  T=T, lm=lm, lx=lx, Bq=Bq, Sv=Sv, kat=kat, g=g)

        # Non-SI units accepted for use with SI units [BIPM2006, Table 6]
        minute = 60*s
        hr = 60*minute
        day = 24*hr
        degree = 2*pi*rad/360
        L = (0.1*m)**3
        self.U.update(minute=minute, hr=hr, day=day, degree=degree, L=L)

        # Derived physical constants
        # Electromagnetism
        G_0 = 2/R_K
        Phi_0 = 1/k_J
        q = G_0*Phi_0
        h = 2*q*Phi_0
        alpha = pi*1e-7*c*s*G_0/(m*S)
        Z_0 = 2*R_K*alpha
        mu_0 = Z_0/c
        epsilon_0 = 1/(Z_0*c)
        k_A = mu_0/(4*pi)
        k_e = k_A*c**2
        E_h = 2*R_inf*h*c
        eV = q*V
        self.U.update(G_0=G_0, Phi_0=Phi_0, q=q, h=h, alpha=alpha, Z_0=Z_0,
                      mu_0=mu_0, epsilon_0=epsilon_0, k_A=k_A, k_e=k_e,
                      E_h=E_h, eV=eV)
        # Electrochemistry
        N_A = k_F/q
        self.U.update(N_A=N_A)
        # Thermal physics
        k_B = R/N_A
        c_1 = cyc*h*c**2
        c_2 = h*c/k_B
        c_3_lambda = c_2/4.965114231744276
        c_3_f = 2.821439372122079*k_B/h
        sigma = 2*pi*(k_B*pi)**4/(15*(h*rad)**3*c**2)
        self.U.update(k_B=k_B, c_1=c_1, c_2=c_2, c_3_lambda=c_3_lambda,
                      c_3_f=c_3_f, sigma=sigma)

        # Selected other non-SI units from [BIPM2006, Table 8]
        bar = 1e5*Pa
        angstrom = 1e-10*m
        self.U.update(bar=bar, angstrom=angstrom)

        # Additional units that are useful for fuel cells
        atm = 101325*Pa
        M = mol/L
        cc = (0.01*m)**3
        yr = 365.25*day
        self.U.update(atm=atm, cc=cc, yr=yr)
        self.U.update({'%': 0.01})


class Info:
    """Shortcuts to the "get" methods in :class:`SimRes`
    """
    description = SimRes.get_description
    """Alias for :meth:`SimRes.get_description`"""
    dimension = SimRes.get_dimension
    """Alias for :meth:`SimRes.get_dimension`"""
    indices_wi_times = SimRes.get_indices_wi_times
    """Alias for :meth:`SimRes.get_indices_wi_times`"""
    IV = SimRes.get_IV
    """Alias for :meth:`SimRes.get_IV`"""
    FV = SimRes.get_FV
    """Alias for :meth:`SimRes.get_FV`"""
    times = SimRes.get_times
    """Alias for :meth:`SimRes.get_times`"""
    unit = SimRes.get_unit
    """Alias for :meth:`SimRes.get_unit`"""
    values = SimRes.get_values
    """Alias for :meth:`SimRes.get_values`"""
    values_at_times = SimRes.get_values_at_times
    """Alias for :meth:`SimRes.get_description`"""


if __name__ == '__main__':
    """Test the contents of this file."""
    import doctest
    doctest.testmod()
    exit()

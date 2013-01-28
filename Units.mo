within FCSys;
package Units "Constants and units of physical measure"
  extends Modelica.Icons.Package;
  model Evaluate "Evaluate the values assigned to constants and units"
    extends Modelica.Icons.Example;

    // -----------------------------------------------------------------------
    // Mathematical constants

    final constant Q.Number pi=U.pi "pi";
    final constant Q.Number e=U.e "Euler number";

    // -----------------------------------------------------------------------
    // Base physical constants and units

    final constant Q.Angle rad=U.rad "radian";
    final constant Q.Wavenumber R_inf=U.R_inf "Rydberg constant";
    final constant Q.Velocity c=U.c "speed of light in vacuum";
    final constant Q.MagneticFluxReciprocal k_J=U.k_J "Josephson constant";
    final constant Q.ResistanceElectrical R_K=U.R_K "von Klitzing constant";
    final constant Q.PowerRadiant 'cd'=U.'cd' "candela";
    final constant Q.Number k_F=U.k_F "Faraday constant";
    final constant Q.Number R=U.R "gas constant";

    // -----------------------------------------------------------------------
    // Empirical units

    final constant Q.Length m=U.m "meter";
    final constant Q.Time s=U.s "second";
    final constant Q.MagneticFlux Wb=U.Wb "weber";
    final constant Q.ConductanceElectrical S=U.S "siemen";
    final constant Q.Amount mol=U.mol "mole";
    final constant Q.Potential K=U.K "kelvin";

    // -----------------------------------------------------------------------
    // SI prefixes [BIPM2006, Table 5]

    final constant Q.Number yotta=U.yotta "yotta (Y)";
    final constant Q.Number zetta=U.zetta "zetta (Z)";
    final constant Q.Number exa=U.exa "exa (E)";
    final constant Q.Number peta=U.peta "peta (P)";
    final constant Q.Number tera=U.tera "tera (T)";
    final constant Q.Number giga=U.giga "giga (G)";
    final constant Q.Number mega=U.mega "mega (M)";
    final constant Q.Number kilo=U.kilo "kilo (k)";
    final constant Q.Number hecto=U.hecto "hecto (h)";
    final constant Q.Number deca=U.deca "deca (da)";
    final constant Q.Number deci=U.deci "deci (d)";
    final constant Q.Number centi=U.centi "centi (c)";
    final constant Q.Number milli=U.milli "milli (m)";
    final constant Q.Number micro=U.micro "micro (u)";
    final constant Q.Number nano=U.nano "nano (n)";
    final constant Q.Number pico=U.pico "pico (p)";
    final constant Q.Number femto=U.femto "femto (f)";
    final constant Q.Number atto=U.atto "atto (a)";
    final constant Q.Number zepto=U.zepto "zepto (z)";
    final constant Q.Number yocto=U.yocto "yocto (y)";

    // -----------------------------------------------------------------------
    // SI base units [BIPM2006, Table 1] and intermediate units

    final constant Q.Potential V=U.V "volt";
    final constant Q.Current A=U.A "ampere";
    final constant Q.Amount C=U.C "coulomb";
    final constant Q.Energy J=U.J "joule";
    final constant Q.Velocity2 Sv=U.Sv "sievert";
    final constant Q.Mass kg=U.kg "kilogram ";

    // -----------------------------------------------------------------------
    // Coherent derived units in the SI with special names and symbols
    // [BIPM2006, Table 3]

    final constant Q.Angle cyc=U.cyc "cycle";
    final constant Q.Frequency Hz=U.Hz "hertz";
    final constant Q.Angle2 sr=U.sr "steradian";
    final constant Q.Force N=U.N "newton";
    final constant Q.Pressure Pa=U.Pa "pascal";
    final constant Q.Power W=U.W "watt";
    final constant Q.Capacitance F=U.F "farad";
    final constant Q.ResistanceElectrical ohm=U.ohm "ohm (Omega)";
    final constant Q.Inductance H=U.H "henry";
    final constant Q.MagneticFluxAreic T=U.T "tesla";
    final constant Q.Power lm=U.lm "lumen";
    final constant Q.PowerAreic lx=U.lx "lux";
    final constant Q.Frequency Bq=U.Bq "becquerel";
    final constant Q.Velocity2 Gy=U.Gy "gray";
    final constant Q.Current kat=U.kat "katal";
    final constant Q.Mass g=U.g "gram";

    // -----------------------------------------------------------------------
    // Non-SI units accepted for use with SI units [BIPM2006, Table 6].

    final constant Q.Time min=U.min "minute";
    final constant Q.Time hr=U.hr "hour";
    final constant Q.Time day=U.day "day";
    final constant Q.Angle degree=U.degree "degree";
    final constant Q.Volume L=U.L "liter (L or l)";

    // -----------------------------------------------------------------------
    // Physical constants

    // Electromagnetism
    final constant Q.ConductanceElectrical G_0=U.G_0 "conductance quantum";
    final constant Q.MagneticFlux Phi_0=U.Phi_0 "magnetic flux quantum";
    final constant Q.Amount q=U.q "elementary charge";
    final constant Q.MomentumAngular h=U.h "Planck constant";
    final constant Q.Number alpha=U.alpha "fine-structure constant";
    final constant Q.ResistanceElectrical Z_0=U.Z_0
      "characteristic impedance of vacuum";
    final constant Q.Permeability mu_0=U.mu_0 "magnetic constant";
    final constant Q.Permittivity epsilon_0=U.epsilon_0 "electric constant";
    final constant Q.Permeability k_A=U.k_A "magnetic force constant";
    final constant Q.PermittivityReciprocal k_e=U.k_e "Coulomb constant";
    final constant Q.Energy E_h=U.E_h "Hartree energy";
    final constant Q.Energy eV=U.eV "electron volt";

    // Electrochemistry
    final constant Q.AmountReciprocal N_A=U.N_A "Avogadro constant";

    // Thermal physics
    final constant Q.Amount k_B=U.k_B "Boltzmann constant";
    final constant Q.PowerArea c_1=U.c_1 "first radiation constant";
    final constant Q.PotentialPerWavenumber c_2=U.c_2
      "second radiation constant";
    final constant Q.PotentialPerWavenumber c_3_f=U.c_3_f
      "Wien frequency displacement law constant";
    final constant Q.MagneticFluxReciprocal c_3_lambda=U.c_3_lambda
      "Wien wavelength displacement law constant";
    final constant Q.PowerAreicPerPotential4 sigma=U.sigma
      "Stefan-Boltzmann constant";

    // -----------------------------------------------------------------------
    // Selected other non-SI units from [BIPM2006, Table 8]

    final constant Q.Pressure bar=U.bar "bar";
    final constant Q.Length Aring=U.Aring "angstrom";

    // -----------------------------------------------------------------------
    // Additional units that are useful for fuel cells

    final constant Q.Pressure atm=U.atm "atmosphere";
    final constant Q.Pressure kPa=U.kPa "kilopascal";
    final constant Q.Length cm=U.cm "centimeter";
    final constant Q.Length mm=U.mm "millimeter";
    final constant Q.Number '%'=U.'%' "percent";
    final constant Q.AmountVolumic M=U.M "molar";
    final constant Q.Volume cc=U.cc "cubic centimeter";
    annotation (Documentation(info="<html><p>This model is used by the units script (\"FCSys/resources/scripts/units.mos\") to
  establish the values of the units in order to convert values to numbers for display.</p>
  <p>For more information, see the documentation in the
  <a href=\"modelica://FCSys.Units\">FCSys.Units</a> package.</p></html>"),
        Commands(file="resources/scripts/units.mos" "Re-initialize the units."));

  end Evaluate;

  package Bases "Sets of base constants and units"
    extends Modelica.Icons.Package;
    record ScaledFC
      "Base constants and units that are well-scaled for fuel cell simulation and analysis"
      extends U.Bases.Base(
        final R_inf=1e-1*10973731.568539,
        final c=1e-1*299792458,
        final R_K=1e10*25812.8074434,
        final k_J=483597.870e9*sqrt(S*s/1e4)/m,
        final k_F=1,
        final R=1);
      // Note:  k_J = 483597.870e9*sqrt(S*s/x)/m sets kg = x.
      annotation (Documentation(info="<html><p>The values of this record result in the following values for the base SI units
  (besides cd = 1, which is the default):
       <ul>
       <li>A &asymp; 1e-5 (&rArr; 1e3 C &asymp; 1)
       <li>K &asymp; 8.617</li>
       <li>kg = 1e4 (&rArr; 0.1 g = 1)</li>
       <li>m = 10 (&rArr; 10 cm = 1)</li>
       <li>mol &asymp; 96.485 (&rArr; 0.01036 mol &asymp;1)</li>
       <li>s = 100 (&rArr; 10 ms = 1)</li></ul>
   which are well-scaled for the states and efforts of a single-cell
   PEMFC.  Also, with these settings:
       <ul>
       <li>10 m/s &asymp; 1</li>
       <li>9.872e-5 atm &asymp; 1</li></ul></p>  <p>For more information, see the documentation in the
  <a href=\"modelica://FCSys.Units\">FCSys.Units</a> package.</p></html>"),
          Commands(file="resources/scripts/units.mos"
            "Re-initialize the units."));

    end ScaledFC;

    record AK
      "<html>Base constants and units for SI with <i>k</i><sub>F</sub> and <i>R</i> normalized instead of A and K</html>"

      extends U.Bases.Base(
        final R_inf=10973731.568539,
        final c=299792458,
        final R_K=96485.3365^2*25812.8074434,
        final k_J=483597.870e9*sqrt(S*s)/m,
        final 'cd'=1,
        final k_F=1,
        final R=1);
      annotation (Documentation(info="<html><p>The values of the un-normalized SI base units are (see
  \"FCSys/resources/unit-systems.cdf\"):
  <ul>
  <li>A &asymp; 1.03643e-5</li>
  <li>K &asymp; 8.31446</li>
  </ul>
  </p>  <p>For more information, see the documentation in the
  <a href=\"modelica://FCSys.Units\">FCSys.Units</a> package.</p></html>"),
          Commands(file="resources/scripts/units.mos"
            "Re-initialize the units."));

    end AK;

    record Am
      "<html>Base constants and units for SI with <i>k</i><sub>F</sub> and <i>R</i> normalized instead of A and m</html>"

      extends U.Bases.Base(
        final R_inf=sqrt(8.3144621)*10973731.568539,
        final c=299792458/sqrt(8.3144621),
        final R_K=(96485.3365^2*25812.8074434)/8.3144621,
        final k_J=483597.870e9*sqrt(S*s)/m,
        final 'cd'=1,
        final k_F=1,
        final R=1);
      annotation (Documentation(info="<html><p>The values of the un-normalized SI base units are (see
  \"FCSys/resources/unit-systems.cdf\"):
  <ul>
  <li>A &asymp; 1.03643e-5</li>
  <li>m &asymp; 0.346803</li>
  </ul>
  </p>  <p>For more information, see the documentation in the
  <a href=\"modelica://FCSys.Units\">FCSys.Units</a> package.</p></html>"),
          Commands(file="resources/scripts/units.mos"
            "Re-initialize the units."));

    end Am;

    record As
      "<html>Base constants and units for SI with <i>k</i><sub>F</sub> and <i>R</i> normalized instead of A and s</html>"
      extends U.Bases.Base(
        final R_inf=10973731.568539,
        final c=299792458/sqrt(8.3144621),
        final R_K=(96485.3365^2*25812.8074434)/sqrt(8.3144621),
        final k_J=483597.870e9*sqrt(S*s)/m,
        final 'cd'=1,
        final k_F=1,
        final R=1);
      annotation (Documentation(info="<html><p>The values of the un-normalized SI base units are (see
  \"FCSys/resources/unit-systems.cdf\"):
  <ul>
  <li>A &asymp; 3.59436e-6</li>
  <li>s &asymp; 2.88348</li>
  </ul>
  </p>

  <p>For more information, see the documentation in the
  <a href=\"modelica://FCSys.Units\">FCSys.Units</a> package.</p></html>"),
          Commands(file="resources/scripts/units.mos"
            "Re-initialize the units."));

    end As;

    record Kmol
      "<html>Base constants and units for SI with <i>k</i><sub>F</sub> and <i>R</i> normalized instead of K and mol</html>"
      extends U.Bases.Base(
        final R_inf=10973731.568539,
        final c=299792458,
        final R_K=25812.8074434,
        final k_J=483597.870e9*sqrt(S*s)/m,
        final 'cd'=1,
        final k_F=1,
        final R=1);
      annotation (Documentation(info="<html><p>The values of the un-normalized SI base units are (see
  \"FCSys/resources/unit-systems.cdf\"):
  <ul>
  <li>K &asymp; 8.61733e-5</li>
  <li>mol &asymp; 96485.3</li>
  </ul>
  </p>  <p>For more information, see the documentation in the
  <a href=\"modelica://FCSys.Units\">FCSys.Units</a> package.</p></html>"),
          Commands(file="resources/scripts/units.mos"
            "Re-initialize the units."));

    end Kmol;

    record Ks
      "<html>Base constants and units for SI with <i>k</i><sub>F</sub> and <i>R</i> normalized instead of K and s</html>"
      extends U.Bases.Base(
        final R_inf=10973731.568539,
        final c=96485.3365*299792458,
        final R_K=96485.3365^3*25812.8074434,
        final k_J=483597.870e9*sqrt(S*s)/m,
        final 'cd'=1,
        final k_F=1,
        final R=1);
      annotation (Documentation(info="<html><p>The values of the un-normalized SI base units are (see
  \"FCSys/resources/unit-systems.cdf\"):
  <ul>
  <li>K &asymp; 7.74028e10</li>
  <li>s &asymp; 1.03643e-5</li>
  </ul></p>
  <p>For more information, see the documentation in the
  <a href=\"modelica://FCSys.Units\">FCSys.Units</a> package.</p></html>"),
          Commands(file="resources/scripts/units.mos"
            "Re-initialize the units."));

    end Ks;

    record mmol
      "<html>Base constants and units for SI with <i>k</i><sub>F</sub> and <i>R</i> normalized instead of m and mol</html>"
      extends U.Bases.Base(
        final R_inf=sqrt(8.3144621/96485.3365)*10973731.568539,
        final c=sqrt(96485.3365/8.3144621)*299792458,
        final R_K=(96485.3365*25812.8074434)/8.3144621,
        final k_J=483597.870e9*sqrt(S*s)/m,
        final 'cd'=1,
        final k_F=1,
        final R=1);
      annotation (Documentation(info="<html><p>The values of the un-normalized SI base units are (see
  \"FCSys/resources/unit-systems.cdf\"):
  <ul>
  <li>m &asymp; 107.724</li>
  <li>mol &asymp; 96485.3</li>
  </ul>
  </p>  <p>For more information, see the documentation in the
  <a href=\"modelica://FCSys.Units\">FCSys.Units</a> package.</p></html>"),
          Commands(file="resources/scripts/units.mos"
            "Re-initialize the units."));

    end mmol;

    record ms
      "<html>Base constants and units for SI with <i>k</i><sub>F</sub> and <i>R</i> normalized instead of m and s</html>"
      extends U.Bases.Base(
        final R_inf=96485.3365*sqrt(8.3144621)*10973731.568539,
        final c=299792458/sqrt(8.3144621),
        final R_K=(96485.3365*25812.8074434)/8.3144621,
        final k_J=483597.870e9*sqrt(S*s)/m,
        final 'cd'=1,
        final k_F=1,
        final R=1);
      annotation (Documentation(info="<html><p>The values of the un-normalized SI base units are (see
  \"FCSys/resources/unit-systems.cdf\"):
  <ul>
  <li>m &asymp; 3.59436e-6</li>
  <li>s &asymp; 1.03643e-5</li>
  </ul></p>
  <p>For more information, see the documentation in the
  <a href=\"modelica://FCSys.Units\">FCSys.Units</a> package.</p></html>"),
          Commands(file="resources/scripts/units.mos"
            "Re-initialize the units."));

    end ms;

    record mols
      "<html>Base constants and units for SI with <i>k</i><sub>F</sub> and <i>R</i> normalized instead of mol and s</html>"
      extends U.Bases.Base(
        final R_inf=10973731.568539,
        final c=(96485.3365/8.3144621)^(1/3)*299792458,
        final R_K=(96485.3365*25812.8074434)/8.3144621,
        final k_J=483597.870e9*sqrt(S*s)/m,
        final 'cd'=1,
        final k_F=1,
        final R=1);
      annotation (Documentation(info="<html><p>The values of the un-normalized SI base units are (see
  \"FCSys/resources/unit-systems.cdf\"):
  <ul>
  <li>mol &asymp; 4261.73</li>
  <li>s &asymp; 0.0441697</li>
  </ul></p>
  <p>For more information, see the documentation in the
  <a href=\"modelica://FCSys.Units\">FCSys.Units</a> package.</p></html>"),
          Commands(file="resources/scripts/units.mos"
            "Re-initialize the units."));

    end mols;

    record BasisGaussian
      "<html>Base constants and units for Gaussian units (<i>k</i><sub>A</sub> = <i>k</i><sub>e</sub> = 1)</html>"
      extends U.Bases.Base(final c=1, final R_K=2*pi/alpha);
      annotation (Documentation(info="<html><p>Gaussian systems of units impose:
  <ul>
  <li><i>k</i><sub>A</sub> = 1 &rArr; <i>R</i><sub>K</sub>/<i>c</i> = 2&pi;/&alpha;</li>
  <li><i>k</i><sub>e</sub> = 1 &rArr; <i>R</i><sub>K</sub>*<i>c</i> = 2&pi;/&alpha;</li>
  </ul>
  Together, <i>c</i> = 1 and <i>R</i><sub>K</sub> = 2&pi;/&alpha;</p>  <p>For more information, see the documentation in the
  <a href=\"modelica://FCSys.Units\">FCSys.Units</a> package.</p></html>"));

    end BasisGaussian;

    record BasisLH
      "<html>Base constants and units for Lorentz-Heaviside units (&mu;<sub>0</sub> = &epsilon;<sub>0</sub> = 1)</html>"
      extends U.Bases.Base(final c=1, final R_K=1/(2*U.alpha));
      annotation (Documentation(info="<html><p>Lorentz-Heaviside systems of units impose:
  <ul>
  <li>&mu;<sub>0</sub> = 1 &rArr; <i>R</i><sub>K</sub>/<i>c</i> = 1/(2&alpha;)</li>
  <li>&epsilon;<sub>0</sub> = 1 &rArr; <i>R</i><sub>K</sub>*<i>c</i> = 1/(2&alpha;)</li>
  </ul>
  Together, <i>c</i> = 1 and <i>R</i><sub>K</sub> = 1/(2&alpha;)</p>  <p>For more information, see the documentation in the
  <a href=\"modelica://FCSys.Units\">FCSys.Units</a> package.</p></html>"));

    end BasisLH;

    record Base "Base constants and units"

      final constant Q.Angle rad=1 "radian";
      // SI unit of rotation or planar angle
      // This relation is from [BIPM2006, Table 3].  The radian is currently
      // not adjustable because BIPM doesn't explicitly use angle in the
      // definitions of Hz, sr, etc. and NIST doesn't explicitly use angle
      // in the relations for R_inf, c_3_nu, etc [NIST2010].
      constant Q.Wavenumber R_inf=1
        "<html>Rydberg constant (R<sub>&infin;</sub>)</html>";
      // The SI unit length (meter) is inversely proportional to this value,
      // which should be increased for larger characteristic lengths.
      constant Q.Velocity c=1 "<html>speed of light in vacuum (c)</html>";
      // The SI unit time (second) is inversely proportional to this value (and
      // R_inf), which should be increased for larger characteristic times.
      constant Q.MagneticFluxReciprocal k_J=1
        "<html>Josephson constant (k<sub>J</sub>)</html>";
      // The SI unit of magnetic flux (weber) is inversely proportional to this
      // value, which should be increased for larger magnetic flux numbers.
      // Also, the SI unit of charge (coulomb) is inversely proportional to this
      // value.
      constant Q.ResistanceElectrical R_K=1
        "<html>von Klitzing constant (R<sub>K</sub>)</html>";
      // The SI unit of electrical conductance (siemen) is inversely
      // proportional to this value, which should be increased for larger
      // characteristic conductances.  Also, the SI unit of charge (coulomb) is
      // inversely proportional to this value.
      constant Q.PowerRadiant 'cd'=1 "candela";
      // SI unit of luminous intensity
      // From http://en.wikipedia.org/wiki/Luminous_intensity, accessed 11/5/11:
      // "luminous intensity is a measure of the wavelength-weighted power
      // emitted by a light source in a particular direction per unit solid
      // angle, based on the luminosity function, a standardized model of the
      // sensitivity of the human eye."
      constant Q.Number k_F=1 "<html>Faraday constant (k<sub>F</sub>)</html>";
      // The unit of substance (mole) is inversely proportional to this value,
      // which should be increased for larger particle numbers.  If k_F is set
      // to 1, then charge is considered an to be amount of substance.
      constant Q.Number R=1 "gas constant";
      // The unit of temperature (kelvin) is inversely proportional to this
      // value, which should be increased for larger temperature numbers.  If R
      // is set to 1, then temperature is considered to be a potential.
      annotation (Documentation(info="<html><p>For more information, see the notes in the Modelica code and the documentation in the
  <a href=\"modelica://FCSys.Units\">FCSys.Units</a> package.</p></html>"),
          Commands(file="resources/scripts/units.mos"
            "Re-initialize the units."));

    end Base;
    annotation (Documentation(info="<html><p>The International System of Units (SI)-like
  sets in this package are named by listing (in alphabetical order) the two units that are
  <b>not</b> normalized for the sake of setting the Faraday and gas constants equal to one.
  There are eight possible sets of this type.</p>
  <p>For more information, see the documentation in the
  <a href=\"modelica://FCSys.Units\">FCSys.Units</a> package.</p></html>"),
        Commands(file="resources/scripts/units.mos" "Re-initialize the units."));

  end Bases;

  // -----------------------------------------------------------------------
  // Mathematical constants

  function from_degC "From temperature in degree Celsius"
    extends Modelica.SIunits.Conversions.ConversionIcon;
    input Real T_degC "Temperature in degree Celsius";
    output Q.TemperatureAbsolute T "Thermodynamic temperature";

  algorithm
    T := (T_degC + 273.15)*K annotation (Inline=true);
  end from_degC;

  function to_degC "To temperature in degree Celsius"

    extends Modelica.SIunits.Conversions.ConversionIcon;
    input Q.TemperatureAbsolute T "Thermodynamic temperature";
    output Real T_degC "Temperature in degree Celsius";

  algorithm
    T_degC := T/K - 273.15 annotation (Inline=true);
  end to_degC;

  function from_kPag "From gauge pressure in kilopascals"
    extends Modelica.SIunits.Conversions.ConversionIcon;
    input Real p_kPag "Gauge pressure in kilopascals";
    output Q.PressureAbsolute p "Absolute pressure";

  algorithm
    p := p_kPag*kPa + atm annotation (Inline=true);
  end from_kPag;

  function to_kPag "To gauge pressure in kilopascals"
    extends Modelica.SIunits.Conversions.ConversionIcon;
    input Q.PressureAbsolute p "Absolute pressure";
    output Real p_kPag "Gauge pressure in kilopascals";

  algorithm
    p_kPag := (p - atm)/kPa annotation (Inline=true);
  end to_kPag;
  final constant Q.Number pi=2*arccos(0) "<html>pi (<i>&pi;</i>)</html>";
  // Circumference per unit diameter
  final constant Q.Number e=exp(1) "Euler number";
  // Natural base

  // -----------------------------------------------------------------------
  // Base physical constants and units

  replaceable constant Bases.ScaledFC base constrainedby Bases.Base
    "Scalable base constants and units";
  // Note:  The base constants and units may be replaced to suit the scale
  // of the physical system.

  final constant Q.Angle rad=base.rad "radian";
  final constant Q.Wavenumber R_inf=base.R_inf
    "<html>Rydberg constant (<i>R</i><sub>&infin;</sub>)</html>";
  final constant Q.Velocity c=base.c
    "<html>speed of light in vacuum (<i>c</i>)</html>";
  final constant Q.MagneticFluxReciprocal k_J=base.k_J
    "<html>Josephson constant (<i>k</i><sub>J</sub>)</html>";
  final constant Q.ResistanceElectrical R_K=base.R_K
    "<html>von Klitzing constant (<i>R</i><sub>K</sub>)</html>";
  final constant Q.PowerRadiant 'cd'=base.'cd' "candela";
  final constant Q.Number k_F=base.k_F
    "<html>Faraday constant (<i>k</i><sub>F</sub>)</html>";
  final constant Q.Number R=base.R "gas constant";

  // -----------------------------------------------------------------------
  // Empirical units
  // Note:  The values are currently based on the those from [NIST2010].
  // The measured values are used rather than conventional values (where
  // they exist).

  constant Q.Length m=10973731.568539*rad/R_inf "meter";
  // SI unit of length
  // This is the "Rydberg constant" relation [NIST2010].  The unit radian
  // is included to be explicit, although it's currently one by definition
  // [BIPM2006].  The Rydberg constant may be determined by measuring the
  // spectra of hydrogen, deuterium, and antiprotonic helium
  // (http://en.wikipedia.org/wiki/Rydberg_constant).
  constant Q.Time s=299792458*m/c "second";
  // SI unit of time or duration
  // This is the "speed of light in vacuum" relation [NIST2010].  c may
  // be determined (among other ways) by measuring the time for
  // electromagnetic signals to travel to and from spacecraft.
  constant Q.MagneticFlux Wb=483597.870e9/k_J "weber";
  // SI unit of magnetic flux
  // This is the "Josephson constant" relation [NIST2010].  The Josephson
  // constant can be determined by measurements of supercurrent
  // (http://en.wikipedia.org/wiki/Josephson_effect).
  constant Q.ConductanceElectrical S=25812.8074434/R_K "siemen";
  // SI unit of electrical conductance
  // This is the "von Klitzing constant" relation [NIST2010].  The unit
  // radian is included on the denominator for dimensional consistency, but
  // it's one by the current definition [BIPM2006].  The von Klitzing
  // constant may be determined by measuring the quantum hall effect
  // (http://en.wikipedia.org/wiki/Quantum_Hall_effect).
  constant Q.Amount mol=96485.3365*Wb*S/k_F "mole";
  // SI unit of amount of substance
  // This is the "Faraday constant" relation [NIST2010].  The factor Wb*S
  // is the coulomb, which is defined below.  The ratio may be determined by
  // electrochemical experiments relating the amount of charge and the
  // amount of substance involved in a reaction.
  constant Q.Potential K=8.3144621*(Wb*rad)^2*S/(s*mol*R) "kelvin";
  // This is the "molar gas constant" relation [NIST2010].  The factor
  // (Wb*rad)^2*S/s is the joule, which is defined below.  The ratio may be
  // determined by measuring acoustic resonation (see
  // http://nvlpubs.nist.gov/nistpubs/sp958-lide/339-343.pdf).  The gas
  // constant is directly related to the Boltzmann constant through other
  // constants above (see the definition of k_B below).  The Boltzmann
  // constant may be determined (among other ways) by measuring thermal
  // electronic noise
  // (http://en.wikipedia.org/wiki/Johnson-Nyquist_noise).

  // -----------------------------------------------------------------------
  // SI base units [BIPM2006, Table 1] and intermediate units
  // Note:  Only A and kg  are remaining (s, m, S, K, mol, and cd
  // already defined).

  final constant Q.Potential V=Wb*rad/s "volt";
  // SI unit of EMF and rate of magnetic flux
  // Note that currently rad = 1, so this equation is consistent with that
  // of [BIPM2006].  The factor of rad is necessary to dimensional
  // consistency when considering angle.
  final constant Q.Current A=V*S "ampere";
  // SI unit of MMF and rate of charge
  final constant Q.Amount C=A*s "coulomb";
  // SI unit of charge
  final constant Q.Energy J=V*C "joule";
  // SI unit of energy
  final constant Q.Velocity2 Gy=(m/s)^2 "gray";
  // SI unit of specific energy (imparted by radiation into material)
  final constant Q.Mass kg=J/Sv "kilogram ";
  // SI unit of mass

  // -----------------------------------------------------------------------
  // SI prefixes [BIPM2006, Table 5]

  final constant Q.Number yotta=1e24 "yotta (Y)";
  final constant Q.Number zetta=1e21 "zetta (Z)";
  final constant Q.Number exa=1e18 "exa (E)";
  final constant Q.Number peta=1e15 "peta (P)";
  final constant Q.Number tera=1e12 "tera (T)";
  final constant Q.Number giga=1e9 "giga (G)";
  final constant Q.Number mega=1e6 "mega (M)";
  final constant Q.Number kilo=1e3 "kilo (k)";
  final constant Q.Number hecto=1e2 "hecto (h)";
  final constant Q.Number deca=1e1 "deca (da)";
  final constant Q.Number deci=1e-1 "deci (d)";
  final constant Q.Number centi=1e-2 "centi (c)";
  final constant Q.Number milli=1e-3 "milli (m)";
  final constant Q.Number micro=1e-6 "micro (u)";
  final constant Q.Number nano=1e-9 "nano (n)";
  final constant Q.Number pico=1e-12 "pico (p)";
  final constant Q.Number femto=1e-15 "femto (f)";
  final constant Q.Number atto=1e-18 "atto (a)";
  final constant Q.Number zepto=1e-21 "zepto (z)";
  final constant Q.Number yocto=1e-24 "yocto (y)";

  // -----------------------------------------------------------------------
  // Coherent derived units in the SI with special names and symbols
  // [BIPM2006, Table 3]
  // Note:  rad, S, C, Wb, V, J, and Sv have already been defined.  Degree
  // Celsius is defined in FCSys/resources/scripts/units.mos since it includes
  // an offset.

  final constant Q.Angle cyc=2*pi*rad "cycle";
  // This is defined for convenience (not listed by [BIPM2006]).
  final constant Q.Frequency Hz=cyc/s "hertz";
  // SI unit of frequency
  // Note:  Numerically, this doesn't evaluate to Hz = 1/s as stated by
  // BIPM (that relation isn't dimensionally correct when considering
  // angle as a dimension), but allows the conversion of frequency into as
  // cycles per second (Hz) or radians per second (rad/s).  Since BIPM
  // defines rad = 1 (also dimensionally incorrect) and given that
  // cyc = 2*pi*rad, Hz evaluates numerically to 2*pi/s [BIPM2006].
  final constant Q.Angle2 sr=rad^2 "steradian";
  // SI unit of solid angle
  // BIPM currently defines rad = 1; therefore, this is consistent with
  // Table 3 [BIPM2006].
  final constant Q.Force N=J/m "newton";
  // SI unit of force
  final constant Q.Pressure Pa=N/m^2 "pascal";
  // SI unit of pressure
  final constant Q.Power W=J/s "watt";
  // SI unit of power
  final constant Q.Capacitance F=C/V "farad";
  // SI unit of capacitance
  final constant Q.ResistanceElectrical ohm=1/S "<html>ohm (&Omega;)</html>";
  // SI unit of electrical resistance
  final constant Q.Inductance H=V*s/A "henry";
  // SI unit of inductance or permeance
  // Note:  By this definition, H = Wb*rad/A, which is currently equivalent
  // to the definition in Table 3 since BIPM currently defines rad = 1
  // [BIPM2006].
  final constant Q.MagneticFluxAreic T=Wb/m^2 "tesla";
  // SI unit of magnetic flux density
  final constant Q.Power lm='cd'*sr "lumen";
  // Note:  Table 3 gives lm = cd, but this is only so because SI
  // implicitly assigns sr = 1.
  final constant Q.PowerAreic lx=lm/m^2 "lux";
  final constant Q.Frequency Bq=Hz "becquerel";
  // SI unit of frequency
  final constant Q.Velocity2 Sv=Gy "sievert";
  // SI unit of specific energy (imparted by radiation into biological
  // tissue)
  final constant Q.Current kat=mol/s "katal";
  final constant Q.Mass g=kg/kilo "gram";
  // CGS unit of mass
  // The base SI unit of mass includes a prefix.  See section 3.2 of
  // [BIPM2006].

  // -----------------------------------------------------------------------
  // Non-SI units accepted for use with SI units [BIPM2006, Table 6]

  final constant Q.Time min=60*s "minute";
  final constant Q.Time hr=60*min "hour";
  final constant Q.Time day=24*hr "day";
  final constant Q.Angle degree=2*pi*rad/360 "<html>degree (&deg;)</html>";
  final constant Q.Volume L=(deci*m)^3 "liter (L or l)";

  // -----------------------------------------------------------------------
  // Derived physical constants
  // Note:  These are established by definition, but may include
  // transcendental mathematical constants.

  // Electromagnetism
  final constant Q.ConductanceElectrical G_0=2/R_K
    "<html>conductance quantum (<i>G</i><sub>0</sub>)</html>";
  final constant Q.MagneticFlux Phi_0=1/k_J
    "<html>magnetic flux quantum (&Phi;<sub>0</sub>)</html>";
  final constant Q.Amount q=G_0*Phi_0 "elementary charge";
  final constant Q.MomentumAngular h=2*q*Phi_0 "Planck constant";
  // The Planck constant over 2*pi (hbar) isn't included as a unique
  // variable.  The unit of angle (rad or cyc) should be factored into the
  // variable that represents frequency as a quantity.  Then, it's
  // unnecessary to use hbar, e.g.:
  //     hbar = h ~= 1.0545e-34*J/Hz ~= 6.6260e-34*J*s/cyc,
  // where Hz = rad/s.  Currently, rad = 1 (see U.Bases.Base).
  final constant Q.Number alpha=pi*1e-7*c*s*G_0/(m*S)
    "<html>fine-structure constant (&alpha;)</html>";
  // The fine-structure constant includes the product of the speed of light
  // in vacuum, expressed in meters per second and conductance quantum,
  // expressed in siemens.
  final constant Q.ResistanceElectrical Z_0=2*R_K*alpha
    "<html>characteristic impedance of vacuum (<i>Z</i><sub>0</sub>)</html>";
  // See:  http://en.wikipedia.org/wiki/Characteristic_impedance_of_vacuum.
  final constant Q.Permeability mu_0=Z_0/c
    "<html>magnetic constant (&mu;<sub>0</sub>)</html>";
  // This is also called the vacuum permeability or permeability of free
  // space.
  final constant Q.Permittivity epsilon_0=1/(Z_0*c)
    "<html>electric constant (&epsilon;<sub>0</sub>)</html>";
  // This is also called the vacuum permittivity or permittivity of free
  // space.
  final constant Q.Permeability k_A=mu_0/(4*pi)
    "<html>magnetic force constant (<i>k</i><sub>A</sub>)</html>";
  // The factor of 4*pi is the result of the line integral that is used to
  // derive Ampere's force law
  // (http://en.wikipedia.org/wiki/Ampere's_force_law).
  final constant Q.PermittivityReciprocal k_e=k_A*c^2
    "<html>Coulomb constant (<i>k</i><sub>e</sub>)</html>";
  // This is the coefficient in Coulomb's law; see
  // http://en.wikipedia.org/wiki/Coulomb's_law.
  final constant Q.Energy E_h=2*R_inf*h*c
    "<html>Hartree energy (<i>E</i><sub>h</sub>)</html>";
  final constant Q.Energy eV=q*V "electron volt";

  // Electrochemistry
  final constant Q.AmountReciprocal N_A=k_F/q
    "<html>Avogadro constant (<i>N</i><sub>A</sub>)</html>";

  // Thermal physics
  final constant Q.Amount k_B=R/N_A
    "<html>Boltzmann constant (<i>k</i><sub>B</sub>)</html>";
  final constant Q.PowerArea c_1=cyc*h*c^2
    "<html>first radiation constant (<i>c</i><sub>1</sub>)</html>";
  // Note the factor of cyc, which is currently 2*pi.
  final constant Q.PotentialPerWavenumber c_2=h*c/k_B
    "<html>second radiation constant (<i>c</i><sub>2</sub>)</html>";
  final constant Q.PotentialPerWavenumber c_3_lambda=c_2/4.965114231744276
    "<html>Wien wavelength displacement law constant (<i>c</i><sub>3 &lambda;</sub>)</html>";
  // See the notes for c_3_nu.  The derivation is similar to that of c_3_nu,
  // but here, the value is the solution to exp(x)*(5 - x) =  5.  The value
  // is from Mathematica (FCSys/resources/math-constants.cdf).  Note that the
  // frequency displacement constant isn't directly related to the
  // wavelength displacement constant:  "Because the spectrum resulting from
  // Planck's law of black body radiation takes a different shape in the
  // frequency domain from that of the wavelength domain, the frequency
  // of the peak emission doesn't correspond to the peak wavelength using
  // the simple relation between frequency, wavelength, and the speed of
  // light"
  // (http://en.wikipedia.org/wiki/Wien's_displacement_law, accessed
  // 1/19/10).
  final constant Q.MagneticFluxReciprocal c_3_f=2.821439372122079*k_B/h
    "<html>Wien frequency displacement law constant (<i>c</i><sub>3 <i>f</i></sub>)</html>";
  // The Wien displacement constant can be derived by setting the derivative
  // of Planck's law to zero and solving for h*f/(k_B*T) in order to find
  // the frequency at maximum radiant intensity.  That procedure results in
  // solving the following equation: exp(x)*(3 - x) = 3, where x is h*f/k_B.
  // The value is from Mathematica (FCSys/resources/math-constants.cdf).
  final constant Q.PowerAreicPerPotential4 sigma=2*pi*(k_B*pi)^4/(15*(h*rad)^3*
      c^2) "<html>Stefan-Boltzmann constant (&sigma;)</html>";
  // Total blackbody radiant intensity per 4th power of temperature
  // See http://en.wikipedia.org/wiki/Stefan-Boltzmann_constant.  The
  // equation can be derived from Planck's law for spectral radiance:
  // B = 2*(h*f)^3/(h*rad*c)^2/(exp(h*f/(k_B*T)) - 1).  It can be written
  // as:
  // B*(h*rad*c)^2/(2*((k_B*T))^3) = (h*f/(k_B*T))^3/(exp(h*f/(k_B*T)) - 1).
  // The RHS is multiplied by pi due to integration over the half sphere and
  // the LHS is multiplied by h*rad/(k_B*T) due to substitution prior to
  // integration (http://en.wikipedia.org/wiki/Stefan-Boltzmann_law).  Now,
  // the equation is:
  // B*h*rad*(h*rad*c)^2/(2*(k_B*T)^4) = pi*(h*f/T)^3/(exp(h*f/(k_B*T)) - 1).
  // The integral of the RHS as a function of (h*f/(k_B*T)) over the entire
  // frequency domain (0 to infinity) is pi^4/15.  Finally, this results in:
  // B_tot/T^4 = 2*pi^5*k_B^4/(15*(h*rad)^3*c^2), where the RHS is the
  // Stefan-Boltzmann constant.  Here, the unit rad has been included for
  // dimensional consistency.

  // -----------------------------------------------------------------------
  // Selected other non-SI units from [BIPM2006, Table 8]
  // Note:  Logarithmic ratios have been excluded because they can't be
  // represented in Dymola's unit conversion GUI.

  final constant Q.Pressure bar=1e5*Pa "bar";
  final constant Q.Length Aring=0.1*nano*m "<html>angstrom (&#8491;)</html>";

  // -----------------------------------------------------------------------
  // Additional units that are useful for fuel cells

  final constant Q.Pressure atm=101325*Pa "atmosphere";
  // Value from "standard atmosphere" [NIST2010]
  final constant Q.Pressure kPa=kilo*Pa "kilopascal";
  final constant Q.Length cm=centi*m "centimeter";
  final constant Q.Length mm=milli*m "millimeter";
  final constant Q.Number '%'=centi "percent (%)";
  final constant Q.AmountVolumic M=U.mol/U.L "molar";
  final constant Q.Volume cc=U.cm^3 "cubic centimeter";

  annotation (Documentation(info="<html><p>When a physical variable is assigned a quantity, it is the product of a number
    and a unit [<a href=\"modelica://FCSys.UsersGuide.References\">BIPM2006</a>].  In <a href=\"modelica://FCSys\">FCSys</a>, units are also assigned numeric values in a consistent
    manner.  A unit
    may be the product of any of the units defined in <a href=\"modelica://FCSys.Units\">FCSys.Units</a> raised to any power.
    When a quantity is converted to a number for display, it is divided by its unit.
    This conveniently handles unit conversion&mdash;without preference towards any particular set of units.
     It can also be used to scale
    the floating point values associated with quantities.  In order to provide
    well-scaled quantities, the value of the unit should be increased in inverse proportion
    to a typical value of a quantity expressed as a number in that unit.  For example, if a model
    has a length scale on the order of a nanometer, it may be best to chose a unit system that
    results in a number of 10<sup>9</sup> for the unit meter.
    </p>

    <p>Regardless of a quantity's unit, it is mapped through the relations in <a href=\"modelica://FCSys.Units\">FCSys.Units</a> (<code>U</code>)
    to a number times a product of five independent base constants or units each raised to the appropriate powers.
    These five are: the Rydberg constant (<code>U.R_inf</code>), speed of light in vacuum (<code>U.c</code>),
    Josephson constant (<code>U.k_J</code>), von Klitzing constant (<code>U.R_K</code>), and
    the candela (<code>U.'cd'</code>).  The constants are related to units through the
    experimentally-derived definitions (e.g., <code>U.Wb=483597.870e9/U.k_J</code>).
    Here the designation of \"constant\" is in the physical sense; even the units
    are <code>constant</code> in the Modelica sense.  These experimental relations, along with
    direct relations among the units (e.g., <code>U.V=U.Wb*U.rad/U.s</code>), establish the values
    of all the units of the International System of Units (SI) and other convenient units for
    modeling fuel cells.</p>

    <p>This method of defining the base units from measured physical quantities (rather than
    vice versa) is natural and reflects the way that standards organizations define units.
    Here, the experimentally-derived relations for the physical constants are
    used rather the \"exact\" or \"conventional\" ones, where they exist
    [<a href=\"modelica://FCSys.UsersGuide.References\">NIST2010</a>].</p>

    <p>There are five independent constants or units in <a href=\"modelica://FCSys.Units\">FCSys.Units</a>,
    but the International System of Units (SI) has seven independent base units (<code>U.kg</code>,
    <code>U.m</code>, <code>U.s</code>, <code>U.A</code>, <code>U.mol</code>, <code>U.K</code>, and
    <code>U.'cd'</code>).  In <a href=\"modelica://FCSys\">FCSys</a>, two additional constraints are imposed in order
    to simplify the model equations and allow electrons and chemical species to be to represented by the
    same base <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.
    First, the Faraday constant (<code>U.k_F</code> or <code>96485.3399*U.C/U.mol</code>)
    is normalized to one.
    This implies that the mole (<code>U.mol</code>) is proportional to the coulomb (<code>U.C</code> or <code>U.Wb*U.S</code>),
    which is considered a number of reference particles given the charge number.
    Also, the gas constant (<code>U.R</code> or <code>8.314472*U.J/(U.mol*U.K)</code>) is normalized to one.
    Therefore, the kelvin (<code>U.K</code>) is proportional to the volt
    (<code>U.V</code> or <code>U.J/U.C</code>). In addition, the radian <code>U.rad</code> is defined as a base constant.  However, it must
    be set equal to one in the current version International System of Units (SI)
    [<a href=\"modelica://FCSys.UsersGuide.References\">BIPM2006</a>].</p>

    <p>The user should be aware of several aspects of the implementation.
    A script
    (\"FCSys/resources/scripts/units.mos\") runs a model (<a href=\"modelica://Units.Evaluate\">Units.Evaluate</a>)
    to determine the values of the units and defines the conversions to display values as numbers
    in units.  This script is set to run when <a href=\"modelica://FCSys\">FCSys</a> is loaded via
    \"FCSys/load.mos\" or from a link in the \"Commands\" menu of <a href=\"modelica://FCSys.Units\">FCSys.Units</a> (when viewed within
    Dymola).  Unless the values of the base constants and units are changed during an editing session, it should be sufficient
    to simply allow the script to run when the package is loaded.
    When the <code>der()</code> operator is utilized, it must be divided by
    the unit second (<code>U.s</code>), since this scaling is not automatic in Modelica's specification
    and Dymola's implementation.  Currently, it is not possible to define a function that implements
    <code>der()</code> with built-in scaling.</p>

  <p>Note that the common \"natural\" units systems (Planck, Stoney, Hartree, Rydberg,
  or Natural) cannot be implemented given the constraints that the gas and
  Faraday constants equal one.  Since these systems set the Boltzmann constant equal
  to one, the Planck constant must equal the von Klitzing constant.  That is not the
  case for any of these systems of units
  (<a href=\"http://en.wikipedia.org/wiki/Natural_units\">http://en.wikipedia.org/wiki/Natural_units</a>).
  The structure of <a href=\"modelica://FCSys.Units\">FCSys.Units</a> allows the constraints on the Faraday and gas constants
  to be relaxed, but the models in <a href=\"modelica://FCSys\">FCSys</a> generally do not.</p>

  <p>This package also contains functions (e.g., <a href=\"modelica://FCSys.Units.to_degC\">to_degC</a>) that
  convert quantities from the unit system defined in <a href=\"modelica://FCSys\">FCSys</a> to quantities 
  expressed in units.  Functions are
  included for units that involve offsets<!-- or other functions besides simple scaling-->.  
  For conversions that require just a scaling factor, it is best to use the 
  units directly.  For example, to convert from potential in volts use <code>v = v_V*U.V</code>,
  where <code>v</code> is potential and <code>v_V</code> is potential expressed in volts.</p>

  <p>Although it is not necessary in <a href=\"http://www.modelica.org\">Modelica</a>, the declarations
  in this package are presorted so that they can be easily ported to imperative or causal languages (e.g.,
  <a href=\"http://www.python.org\">Python</a>, C).</p>

  <p>For more information, see the related paper [<a href=\"modelica://FCSys.UsersGuide.References\">Davies and Paredis, 2012</a>].</p>
  <p>
  <b>Licensed by the Georgia Tech Research Corporation under the Modelica License 2</b><br>
  Copyright 2007&ndash;2012, Georgia Tech Research Corporation.
  </p>
  <p>
  <i>This Modelica package is <u>free</u> software and the use is completely at <u>your own risk</u>;
  it can be redistributed and/or modified under the terms of the Modelica License 2. For license conditions (including the
  disclaimer of warranty) see <a href=\"modelica://FCSys.UsersGuide.ModelicaLicense2\">
  FCSys.UsersGuide.ModelicaLicense2</a> or visit <a href=\"http://www.modelica.org/licenses/ModelicaLicense2\">
  http://www.modelica.org/licenses/ModelicaLicense2</a>.</i>
  </p></html>"), Commands(file="resources/scripts/units.mos"
        "Re-initialize the units."));
end Units;

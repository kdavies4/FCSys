within FCSys;
package Units "Constants and units of physical measure"
  function setup "Establish conversions to display quantities in units"
    import Modelica.Utilities.Streams.print;

  algorithm
    print("Establishing display units...");

    // ------------------------------------------------------------------------
    // Default display units
    // ------------------------------------------------------------------------
    // If units other than those in the displayUnit attribute of the
    // quantities in FCSys.Quantities should be used by default, then specify
    // them here.  Be sure that the desired unit is included in a
    // defineUnitConversion command below.
    // Generated from FCSys/resources/quantities.xls, 2013-5-23

    defineDefaultDisplayUnit("l/(T.s)", "cm/s2")
      "for derivative of velocity in Dymola";
    defineDefaultDisplayUnit("l/T2", "cm/s2") "Acceleration";
    defineDefaultDisplayUnit("N", "C") "Amount";
    defineDefaultDisplayUnit("1/N", "1/mol") "Reciprocal amount";
    defineDefaultDisplayUnit("A", "rad") "Angle";
    defineDefaultDisplayUnit("A2", "sr") "Solid angle";
    defineDefaultDisplayUnit("l2", "cm2") "Area";
    defineDefaultDisplayUnit("l2/N", "m2/C") "Specific area";
    defineDefaultDisplayUnit("N2.T2/(l2.m)", "uF") "Capacitance";
    defineDefaultDisplayUnit("N2.T/(l2.m)", "S") "Electrical conductance";
    defineDefaultDisplayUnit("N/s", "A") "for derivative of amount in Dymola";
    defineDefaultDisplayUnit("N/T", "A") "Current";
    defineDefaultDisplayUnit("N/(l2.T)", "A/cm2") "Areic current";
    defineDefaultDisplayUnit("N/(T.s)", "A/s")
      "for derivative of current in Dymola";
    defineDefaultDisplayUnit("N/T2", "A/s") "Rate of current";
    defineDefaultDisplayUnit("N/l3", "C/cc") "Density";
    defineDefaultDisplayUnit("N/(l3.T)", "C/(cc.s)") "DensityRate";
    defineDefaultDisplayUnit("l2.m/T2", "J") "Energy";
    defineDefaultDisplayUnit("l.T/m", "1/(Pa.s)") "Fluidity";
    defineDefaultDisplayUnit("l.m/T2", "N") "Force";
    defineDefaultDisplayUnit("A/T", "rad/s") "Frequency";
    defineDefaultDisplayUnit("l2.m/N2", "uH") "Inductance";
    defineDefaultDisplayUnit("l", "cm") "Length";
    defineDefaultDisplayUnit("l/N", "m/mol") "Specific length";
    defineDefaultDisplayUnit("l2.m/(A.N.T)", "Wb") "Magnetic flux";
    defineDefaultDisplayUnit("m/(A.N.T)", "T") "Areic magnetic flux";
    defineDefaultDisplayUnit("A.N.T/(l2.m)", "1/Wb") "Reciprocal magnetic flux";
    defineDefaultDisplayUnit("m", "g") "Mass";
    defineDefaultDisplayUnit("m/s", "g/s") "for derivative of mass in Dymola";
    defineDefaultDisplayUnit("m/N", "g/mol") "Specific mass";
    defineDefaultDisplayUnit("N.T/m", "cm2/(V.s)") "Mobility";
    defineDefaultDisplayUnit("l2.m/(A.T)", "J.s/rad") "Rotational momentum";
    defineDefaultDisplayUnit("l.m/(N.T)", "kg.m/(C.s)")
      "Specific translational momentum";
    defineDefaultDisplayUnit("l.m/N2", "H/m") "Permeability";
    defineDefaultDisplayUnit("N2.T2/(l3.m)", "F/m") "Permittivity";
    defineDefaultDisplayUnit("l3.m/(N2.T2)", "m/H") "Reciprocal permittivity";
    defineDefaultDisplayUnit("l2.m/(N.T2)", "V") "Potential";
    defineDefaultDisplayUnit("l3.m/(A.N.T2)", "V.m/rad")
      "Potential per wavenumber";
    defineDefaultDisplayUnit("l2.m/(N.T2.s)", "V/s")
      "for derivative of potential in Dymola";
    defineDefaultDisplayUnit("l2.m/(N.T3)", "V/s") "Rate of potential";
    defineDefaultDisplayUnit("l2.m/T3", "W") "Power";
    defineDefaultDisplayUnit("l4.m/T3", "W.m2") "Power times area";
    defineDefaultDisplayUnit("m/T3", "W/m2") "Areic power";
    defineDefaultDisplayUnit("m.T5/l8", "W/(m2.K4)")
      "Areic power per 4th power of potential";
    defineDefaultDisplayUnit("l2.m/(A2.T3)", "cd") "Radiant power";
    defineDefaultDisplayUnit("m/(l.T2)", "kPa") "Pressure";
    defineDefaultDisplayUnit("m/(l.T2.s)", "Pa/s")
      "for derivative of pressure in Dymola";
    defineDefaultDisplayUnit("m/(T2.l.s)", "Pa/s")
      "for derivative of pressure in Dymola";
    defineDefaultDisplayUnit("m/(l.T3)", "Pa/s") "Rate of pressure";
    defineDefaultDisplayUnit("l.T2/m", "1/kPa") "Reciprocal pressure";
    defineDefaultDisplayUnit("l2.m/(N2.T)", "ohm") "Electrical resistance";
    defineDefaultDisplayUnit("l.T/N", "cm/A") "Resistivity";
    defineDefaultDisplayUnit("T/l2", "s/mm2") "Material resistivity";
    defineDefaultDisplayUnit("T", "s") "Time";
    defineDefaultDisplayUnit("l/T", "cm/s") "Velocity";
    defineDefaultDisplayUnit("l2/T2", "Sv") "Squared velocity";
    defineDefaultDisplayUnit("T/l", "s/m") "Reciprocal velocity";
    defineDefaultDisplayUnit("l3", "cc") "Volume";
    defineDefaultDisplayUnit("l3/T", "L/min") "Rate of volume";
    defineDefaultDisplayUnit("l3/N", "cc/C") "Specific volume";
    defineDefaultDisplayUnit("l3/(N.T)", "cc/(C.s)") "Rate of specific volume";
    defineDefaultDisplayUnit("l3/(N.s)", "cc/(C.s)")
      "for derivative of specific volume in Dymola";
    defineDefaultDisplayUnit("A/l", "rad/m") "Wavenumber";

    // -----------------------------------------------------------------------
    // Conversions to display quantities in units
    // -----------------------------------------------------------------------

    defineUnitConversion(
        "l/(T.s)",
        "cm/s2",
        s/cm) "for derivative of velocity in Dymola";
    defineUnitConversion(
        "l/T2",
        "cm/s2",
        s^2/cm) "Acceleration";
    defineUnitConversion(
        "l/(T.s)",
        "m/s2",
        s/m) "for derivative of velocity in Dymola";
    defineUnitConversion(
        "l/T2",
        "m/s2",
        s^2/m) "Acceleration";
    defineUnitConversion(
        "N",
        "C",
        1/C) "Amount";
    defineUnitConversion(
        "N",
        "mol",
        1/mol) "Amount";
    defineUnitConversion(
        "N",
        "q",
        1/q) "Amount";
    defineUnitConversion(
        "N",
        "J/K",
        K/J) "Amount";
    defineUnitConversion(
        "1/N",
        "1/C",
        C) "Reciprocal amount";
    defineUnitConversion(
        "1/N",
        "1/mol",
        mol) "Reciprocal amount";
    defineUnitConversion(
        "A",
        "degree",
        1/degree) "Angle";
    defineUnitConversion(
        "A",
        "rad",
        1/rad) "Angle";
    defineUnitConversion(
        "A2",
        "sr",
        1/sr) "Solid angle";
    defineUnitConversion(
        "l2",
        "cm2",
        1/cm^2) "Area";
    defineUnitConversion(
        "l2",
        "m2",
        1/m^2) "Area";
    defineUnitConversion(
        "l2/N",
        "m2/C",
        C/m^2) "Specific area";
    defineUnitConversion(
        "N2.T2/(l2.m)",
        "F",
        1/F) "Capacitance";
    defineUnitConversion(
        "N2.T2/(l2.m)",
        "uF",
        1/(micro*F)) "Capacitance";
    defineUnitConversion(
        "N2.T/(l2.m)",
        "S",
        1/S) "Electrical conductance";
    defineUnitConversion(
        "N/s",
        "A",
        1/C) "for derivative of amount in Dymola";
    defineUnitConversion(
        "N/T",
        "A",
        1/A) "Current";
    defineUnitConversion(
        "N/s",
        "kat",
        1/mol) "for derivative of amount in Dymola";
    defineUnitConversion(
        "N/T",
        "kat",
        1/kat) "Current";
    defineUnitConversion(
        "N/s",
        "W/K",
        K/J) "for derivative of amount in Dymola";
    defineUnitConversion(
        "N/T",
        "W/K",
        K/W) "Current";
    defineUnitConversion(
        "N/(l2.T)",
        "A/cm2",
        cm^2/A) "Areic current";
    defineUnitConversion(
        "N/(l2.T)",
        "kat/cm2",
        cm^2/kat) "Areic current";
    defineUnitConversion(
        "N/(l2.T)",
        "A/m2",
        m^2/A) "Areic current";
    defineUnitConversion(
        "N/(l2.T)",
        "kat/m2",
        m^2/kat) "Areic current";
    defineUnitConversion(
        "N/(T.s)",
        "A/s",
        1/A) "for derivative of current in Dymola";
    defineUnitConversion(
        "N/T2",
        "A/s",
        s/A) "Rate of current";
    defineUnitConversion(
        "N/(T.s)",
        "kat/s",
        1/kat) "for derivative of current in Dymola";
    defineUnitConversion(
        "N/T2",
        "kat/s",
        s/kat) "Rate of current";
    defineUnitConversion(
        "N/l3",
        "M",
        1/M) "Density";
    defineUnitConversion(
        "N/l3",
        "C/cc",
        cc/C) "Density";
    defineUnitConversion(
        "N/l3",
        "C/m3",
        m^3/C) "Density";
    defineUnitConversion(
        "N/(l3.T)",
        "C/(cc.s)",
        cc*s/C) "DensityRate";
    defineUnitConversion(
        "N/(l3.T)",
        "C/(m3.s)",
        m^3*s/C) "DensityRate";
    defineUnitConversion(
        "N/(l3.T)",
        "M/s",
        s/M) "DensityRate";
    defineUnitConversion(
        "l2.m/T2",
        "J",
        1/J) "Energy";
    defineUnitConversion(
        "l.T/m",
        "1/(Pa.s)",
        Pa*s) "Fluidity";
    defineUnitConversion(
        "l.m/T2",
        "N",
        1/N) "Force";
    defineUnitConversion(
        "A/T",
        "Hz",
        1/Hz) "Frequency";
    defineUnitConversion(
        "A/T",
        "rad/s",
        s/rad) "Frequency";
    defineUnitConversion(
        "l2.m/N2",
        "H",
        1/H) "Inductance";
    defineUnitConversion(
        "l2.m/N2",
        "uH",
        1/(micro*H)) "Inductance";
    defineUnitConversion(
        "l",
        "cm",
        1/cm) "Length";
    defineUnitConversion(
        "l",
        "m",
        1/m) "Length";
    defineUnitConversion(
        "l",
        "mm",
        1/mm) "Length";
    defineUnitConversion(
        "l/N",
        "m/C",
        C/m) "Specific length";
    defineUnitConversion(
        "l/N",
        "m/mol",
        mol/m) "Specific length";
    defineUnitConversion(
        "l2.m/(A.N.T)",
        "Wb",
        1/Wb) "Magnetic flux";
    defineUnitConversion(
        "m/(A.N.T)",
        "T",
        1/T) "Areic magnetic flux";
    defineUnitConversion(
        "A.N.T/(l2.m)",
        "1/Wb",
        Wb) "Reciprocal magnetic flux";
    defineUnitConversion(
        "m",
        "g",
        1/g) "Mass";
    defineUnitConversion(
        "m",
        "kg",
        1/kg) "Mass";
    defineUnitConversion(
        "m/s",
        "g/s",
        s/g) "for derivative of mass in Dymola";
    defineUnitConversion(
        "m/s",
        "kg/s",
        s/kg) "for derivative of mass in Dymola";
    defineUnitConversion(
        "m/N",
        "ug/C",
        C/(micro*g)) "Specific mass";
    defineUnitConversion(
        "m/N",
        "g/mol",
        mol/g) "Specific mass";
    defineUnitConversion(
        "m/N",
        "kg/mol",
        mol/kg) "Specific mass";
    defineUnitConversion(
        "N.T/m",
        "C.s/g",
        g/(C*s)) "Mobility";
    defineUnitConversion(
        "N.T/m",
        "cm2/(V.s)",
        V*s/cm^2) "Mobility";
    defineUnitConversion(
        "l2.m/(A.T)",
        "J.s/rad",
        rad/(J*s)) "Rotational momentum";
    defineUnitConversion(
        "l.m/(N.T)",
        "kg.m/(C.s)",
        C*s/(kg*m)) "Specific translational momentum";
    defineUnitConversion(
        "l.m/(N.T)",
        "kg.m/(C.s)",
        C*s/(kg*m)) "Specific translational momentum";
    defineUnitConversion(
        "1",
        "%",
        1/'%') "Number";
    defineUnitConversion(
        "1",
        "J/(mol.K)",
        mol*K/J) "Absolute number";
    defineUnitConversion(
        "l.m/N2",
        "H/m",
        m/H) "Permeability";
    defineUnitConversion(
        "N2.T2/(l3.m)",
        "F/m",
        m/F) "Permittivity";
    defineUnitConversion(
        "l3.m/(N2.T2)",
        "m/H",
        H/m) "Reciprocal permittivity";
    defineUnitConversion(
        "l2.m/(N.T2)",
        "V",
        1/V) "Potential";
    defineUnitConversion(
        "l2.m/(N.T2)",
        "J/mol",
        mol/J) "Potential";
    defineUnitConversion(
        "l2.m/(N.T2)",
        "K",
        1/K) "Absolute potential";
    defineUnitConversion(
        "l3.m/(A.N.T2)",
        "K.m/rad",
        rad/(K*m)) "Potential per wavenumber";
    defineUnitConversion(
        "l3.m/(A.N.T2)",
        "V.m/rad",
        rad/(V*m)) "Potential per wavenumber";
    defineUnitConversion(
        "l2.m/(N.T2.s)",
        "K/s",
        1/K) "for derivative of potential in Dymola";
    defineUnitConversion(
        "l2.m/(N.T3)",
        "K/s",
        s/K) "Rate of potential";
    defineUnitConversion(
        "l2.m/(N.T2.s)",
        "V/s",
        1/V) "for derivative of potential in Dymola";
    defineUnitConversion(
        "l2.m/(N.T3)",
        "V/s",
        s/V) "Rate of potential";
    defineUnitConversion(
        "l2.m/T3",
        "lm",
        1/lm) "Power";
    defineUnitConversion(
        "l2.m/T3",
        "W",
        1/W) "Power";
    defineUnitConversion(
        "l4.m/T3",
        "W.m2",
        1/(W*m^2)) "Power times area";
    defineUnitConversion(
        "m/T3",
        "lm/m2",
        m^2/lm) "Areic power";
    defineUnitConversion(
        "m/T3",
        "W/m2",
        m^2/W) "Areic power";
    defineUnitConversion(
        "m.T5/l8",
        "W/(m2.K4)",
        m^2*K^4/W) "Areic power per 4th power of potential";
    defineUnitConversion(
        "l2.m/(A2.T3)",
        "cd",
        1/'cd') "Radiant power";
    defineUnitConversion(
        "l2.m/(A2.T3)",
        "W/sr",
        sr/W) "Radiant power";
    defineUnitConversion(
        "m/(l.T2)",
        "atm",
        1/atm) "Pressure";
    defineUnitConversion(
        "m/(l.T2)",
        "bar",
        1/bar) "Pressure";
    defineUnitConversion(
        "m/(l.T2)",
        "kPa",
        1/(kilo*Pa)) "Pressure";
    defineUnitConversion(
        "m/(l.T2)",
        "Pa",
        1/Pa) "Pressure";
    defineUnitConversion(
        "m/(l.T2.s)",
        "Pa/s",
        1/Pa) "for derivative of pressure in Dymola";
    defineUnitConversion(
        "m/(T2.l.s)",
        "Pa/s",
        1/Pa) "for derivative of pressure in Dymola";
    defineUnitConversion(
        "m/(l.T3)",
        "Pa/s",
        s/Pa) "Rate of pressure";
    defineUnitConversion(
        "l.T2/m",
        "1/kPa",
        kilo*Pa) "Reciprocal pressure";
    defineUnitConversion(
        "l.T2/m",
        "1/Pa",
        Pa) "Reciprocal pressure";
    defineUnitConversion(
        "l2.m/(N2.T)",
        "ohm",
        1/ohm) "Electrical resistance";
    defineUnitConversion(
        "l.T/N",
        "cm/A",
        A/cm) "Resistivity";
    defineUnitConversion(
        "l.T/N",
        "m/A",
        A/m) "Resistivity";
    defineUnitConversion(
        "l.T/N",
        "m.K/W",
        W/(m*K)) "Resistivity";
    defineUnitConversion(
        "T/l2",
        "s/m2",
        m^2/s) "Material resistivity";
    defineUnitConversion(
        "T/l2",
        "s/mm2",
        mm^2/s) "Material resistivity";
    defineUnitConversion(
        "T",
        "day",
        1/day) "Time";
    defineUnitConversion(
        "T",
        "hr",
        1/hr) "Time";
    defineUnitConversion(
        "T",
        "us",
        1/(micro*s)) "Time";
    defineUnitConversion(
        "T",
        "ms",
        1/(milli*s)) "Time";
    defineUnitConversion(
        "T",
        "min",
        1/min) "Time";
    defineUnitConversion(
        "T",
        "s",
        1/s) "Time";
    defineUnitConversion(
        "l/T",
        "cm/s",
        s/cm) "Velocity";
    defineUnitConversion(
        "l/T",
        "m/s",
        s/m) "Velocity";
    defineUnitConversion(
        "l/T",
        "mm/s",
        s/mm) "Velocity";
    defineUnitConversion(
        "l2/T2",
        "Sv",
        1/Sv) "Squared velocity";
    defineUnitConversion(
        "T/l",
        "s/m",
        m/s) "Reciprocal velocity";
    defineUnitConversion(
        "l3",
        "cc",
        1/cc) "Volume";
    defineUnitConversion(
        "l3",
        "L",
        1/L) "Volume";
    defineUnitConversion(
        "l3",
        "m3",
        1/m^3) "Volume";
    defineUnitConversion(
        "l3/s",
        "L/min",
        min/(L*s)) "for derivative of volume in Dymola";
    defineUnitConversion(
        "l3/T",
        "L/min",
        min/L) "Rate of volume";
    defineUnitConversion(
        "l3/s",
        "cc/s",
        1/(centi*m)^3) "for derivative of volume in Dymola";
    defineUnitConversion(
        "l3/T",
        "cc/s",
        s/cc) "Rate of volume";
    defineUnitConversion(
        "l3/s",
        "m3/s",
        1/m^3) "for derivative of volume in Dymola";
    defineUnitConversion(
        "l3/T",
        "m3/s",
        s/m^3) "Rate of volume";
    defineUnitConversion(
        "l3/N",
        "cc/C",
        C/cc) "Specific volume";
    defineUnitConversion(
        "l3/N",
        "m3/C",
        C/m^3) "Specific volume";
    defineUnitConversion(
        "l3/N",
        "m3/mol",
        mol/m^3) "Specific volume";

    defineUnitConversion(
        "l3/(N.T)",
        "cc/(C.s)",
        C*s/cc) "Rate of specific volume";
    defineUnitConversion(
        "l3/(N.T)",
        "m3/(C.s)",
        C*s/m^3) "Rate of specific volume";
    defineUnitConversion(
        "l3/(N.T)",
        "m3/(mol.s)",
        mol*s/m^3) "Rate of specific volume";
    defineUnitConversion(
        "l3/(N.s)",
        "cc/(C.s)",
        C*s/cc) "for derivative of specific volume in Dymola";
    defineUnitConversion(
        "l3/(N.s)",
        "m3/(C.s)",
        C*s/m^3) "for derivative of specific volume in Dymola";
    defineUnitConversion(
        "l3/(N.s)",
        "m3/(mol.s)",
        mol*s/m^3) "for derivative of specific volume in Dymola";
    defineUnitConversion(
        "A/l",
        "cyc/m",
        m/cyc) "Wavenumber";
    defineUnitConversion(
        "A/l",
        "rad/m",
        m/rad) "Wavenumber";

    print("Done.");
    annotation (Documentation(info="<html><p>This has no inputs or outputs.  
The <code>defineDefaultDisplayUnit</code> and <code>defineUnitConversion</code> functions
used by this function are not defined in the Modelica language (as of version 3.3) but are 
recognized by Dymola.
For more information, see the documentation in
<a href=\"modelica://FCSys.Units\">FCSys.Units</a>.</p></html>"));
  end setup;

  package Examples "Examples"
    extends Modelica.Icons.ExamplesPackage;

    model Evaluate "Evaluate the values assigned to constants and units"
      extends Modelica.Icons.Example;

      // ------------------------------------------------------------------------
      // Mathematical constants

      final constant Q.Number pi=U.pi "pi";
      final constant Q.Number e=U.e "Euler number";

      // ------------------------------------------------------------------------
      // Base physical constants and units

      final constant Q.Angle rad=U.rad "radian";
      final constant Q.Wavenumber R_inf=U.R_inf "Rydberg constant";
      final constant Q.Velocity c=U.c "speed of light in vacuum";
      final constant Q.MagneticFluxReciprocal k_J=U.k_J "Josephson constant";
      final constant Q.ResistanceElectrical R_K=U.R_K "von Klitzing constant";
      final constant Q.PowerRadiant 'cd'=U.'cd' "candela";
      final constant Q.Number k_F=U.k_F "Faraday constant";
      final constant Q.Number R=U.R "gas constant";

      // ------------------------------------------------------------------------
      // Empirical units

      final constant Q.Length m=U.m "meter";
      final constant Q.Time s=U.s "second";
      final constant Q.MagneticFlux Wb=U.Wb "weber";
      final constant Q.ConductanceElectrical S=U.S "siemen";
      final constant Q.Amount mol=U.mol "mole";
      final constant Q.Potential K=U.K "kelvin";

      // ------------------------------------------------------------------------
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

      // ------------------------------------------------------------------------
      // SI base units [BIPM2006, Table 1] and intermediate units

      final constant Q.Potential V=U.V "volt";
      final constant Q.Current A=U.A "ampere";
      final constant Q.Amount C=U.C "coulomb";
      final constant Q.Energy J=U.J "joule";
      final constant Q.Velocity2 Sv=U.Sv "sievert";
      final constant Q.Mass kg=U.kg "kilogram ";

      // ------------------------------------------------------------------------
      // Coherent derived units in SI with special names and symbols [BIPM2006,
      // Table 3]

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

      // ------------------------------------------------------------------------
      // Non-SI units accepted for use with SI units [BIPM2006, Table 6].

      final constant Q.Time min=U.min "minute";
      final constant Q.Time hr=U.hr "hour";
      final constant Q.Time day=U.day "day";
      final constant Q.Angle degree=U.degree "degree";
      final constant Q.Volume L=U.L "liter (L or l)";

      // ------------------------------------------------------------------------
      // Physical constants

      // Electromagnetism
      final constant Q.ConductanceElectrical G_0=U.G_0 "conductance quantum";
      final constant Q.MagneticFlux Phi_0=U.Phi_0 "magnetic flux quantum";
      final constant Q.Amount q=U.q "elementary charge";
      final constant Q.MomentumRotational h=U.h "Planck constant";
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

      // ------------------------------------------------------------------------
      // Selected other non-SI units from [BIPM2006, Table 8]

      final constant Q.Pressure bar=U.bar "bar";
      final constant Q.Length Aring=U.Aring "angstrom";

      // ------------------------------------------------------------------------
      // Additional units that are useful for fuel cells

      final constant Q.Pressure atm=U.atm "atmosphere";
      final constant Q.Pressure kPa=U.kPa "kilopascal";
      final constant Q.Length cm=U.cm "centimeter";
      final constant Q.Length mm=U.mm "millimeter";
      final constant Q.Number '%'=U.'%' "percent";
      final constant Q.Density M=U.M "molar";
      final constant Q.Volume cc=U.cc "cubic centimeter";
      annotation (Documentation(info="<html><p>This model may be used to calculate the values of the
  constants and units.</p>

<p>For more information, see the documentation in the
  <a href=\"modelica://FCSys.Units\">FCSys.Units</a> package.</p></html>"),
          Commands(executeCall=FCSys.Units.setup() "Re-initialize the units."));

    end Evaluate;

  end Examples;
  import arccos = Modelica.Math.acos;
  // Note:  The command line of Dymola 7.4 recognizes arccos() rather than
  // acos().
  extends Modelica.Icons.Package;

  package Bases "Sets of base constants and units"
    extends Modelica.Icons.Package;
    record Gaussian
      "<html>Base constants and units for Gaussian units (<i>k</i><sub>A</sub> = <i>k</i><sub>e</sub> = 1)</html>"
      extends Base(final c=1,final R_K=25812.8074434/(299792458*1e-7));
      annotation (Documentation(info="<html><p>Gaussian systems of units impose:
  <ul>
  <li><i>k</i><sub>A</sub> = 1 &rArr; <i>R</i><sub>K</sub>/<i>c</i> = 2&pi;/&alpha;</li>
  <li><i>k</i><sub>e</sub> = 1 &rArr; <i>R</i><sub>K</sub>*<i>c</i> = 2&pi;/&alpha;</li>
  </ul>
  Together, <i>c</i> = 1 and <i>R</i><sub>K</sub> = 2&pi;/&alpha;</p>

<p>The Gaussian conditions are not sufficient
to fully establish the values of the base constants and units of the
<a href=\"modelica://FCSys.Units\">Units</a> package.  Gaussian units
encompass other systems of units.</p>

<p>For more information, see the documentation in the
  <a href=\"modelica://FCSys.Units\">Units</a> package.</p></html>"));

    end Gaussian;

    record Hartree "Base constants and units for Hartree units"

      extends Base(
        final R_inf=pi*(299792458*1e-7/25812.8074434)^2,
        final c=25812.8074434/(2*pi*299792458*1e-7),
        final k_J=1/pi,
        final R_K=2*pi);

      // R_inf = pi*(299792458*1e-7/25812.8074434)^2

      // **need to update R_inf
      // Use:
      // h = 2*pi
      // q = 1
      // k_B = 1
      // 2*R_inf*h/(q*c*alpha^2) = 1

      // 2*Phi_0 = 2*pi
      // G_0*Phi_0 = 1
      // k_B = 1
      // 4*R_inf/alpha^2 = k_J
      annotation (Documentation(info="<html><p>The candela (<code>'cd'</code>)
  is not final because luminous intensity is not included in Hartree units.</p>

<p>For more information, see the documentation in the
  <a href=\"modelica://FCSys.Units\">Units</a> package.</p></html>"), Commands(
            executeCall=FCSys.Units.setup() "Re-initialize the units."));

    end Hartree;

    record LH
      "<html>Base constants and units for Lorentz-Heaviside units (&mu;<sup>0</sub> = &epsilon;<sup>0</sub> = 1)</html>"
      extends Base(final c=1,final R_K=25812.8074434/(4*pi*299792458*1e-7));
      annotation (Documentation(info="<html><p>Lorentz-Heaviside systems of units impose:
  <ul>
  <li>&mu;<sup>0</sub> = 1 &rArr; <i>R</i><sub>K</sub>/<i>c</i> = 1/(2&alpha;)</li>
  <li>&epsilon;<sup>0</sub> = 1 &rArr; <i>R</i><sub>K</sub>*<i>c</i> = 1/(2&alpha;)</li>
  </ul>
  Together, <i>c</i> = 1 and <i>R</i><sub>K</sub> = 1/(2&alpha;)</p>

<p>The Lorentz-Heaviside conditions are not sufficient
to fully establish the values of the base constants and units of the
<a href=\"modelica://FCSys.Units\">Units</a> package.  Lorentz-Heaviside units
encompass other systems of units.</p>

<p>For more information, see the documentation in the
  <a href=\"modelica://FCSys.Units\">Units</a> package.</p></html>"));

    end LH;

    record Stoney "Base constants and units for Stoney units"

      extends Gaussian(final k_J=2e-7*299792458/25812.8074434);
      annotation (Documentation(info="<html><p>The Rydberg constant (<code>R_inf</code>)
  is not final because the <a href=\"modelica://FCSys.Units\">Units</a> package does not
  include the gravitational constant.  The candela (<code>'cd'</code>)
  is not final because luminous intensity is not included in Stoney units.</p>

<p>For more information, see the documentation in the
  <a href=\"modelica://FCSys.Units\">Units</a> package.</p></html>"), Commands(
            executeCall=FCSys.Units.setup() "Re-initialize the units."));

    end Stoney;

    record ScaledFC
      "Base constants and units that are well-scaled for fuel cell simulation and analysis"
      extends Base(
        final R_inf=1e-1*10973731.568539,
        final c=1e-1*299792458,
        final R_K=1e10*25812.8074434,
        final k_J=483597.870e9*sqrt(S*s/1e4)/m);
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
       <li>9.872e-5 atm &asymp; 1</li></ul></p>

<p>For more information, see the documentation in the
  <a href=\"modelica://FCSys.Units\">Units</a> package.</p></html>"), Commands(
            executeCall=FCSys.Units.setup() "Re-initialize the units."));

    end ScaledFC;

    record SIAK
      "<html>Base constants and units for SI with <i>k</i><sub>F</sub> and <i>R</i> normalized instead of A and K</html>"

      extends Base(
        final R_inf=10973731.568539,
        final c=299792458,
        final R_K=96485.3365^2*25812.8074434,
        final k_J=483597.870e9*sqrt(S*s)/m,
        final 'cd'=1);
      annotation (Documentation(info="<html><p>The values of the un-normalized SI base units are (see
  \"FCSys/resources/unit-systems.cdf\"):
  <ul>
  <li>A &asymp; 1.03643e-5</li>
  <li>K &asymp; 8.31446</li>
  </ul></p>

<p>For more information, see the documentation in the
  <a href=\"modelica://FCSys.Units\">Units</a> package.</p></html>"), Commands(
            executeCall=FCSys.Units.setup() "Re-initialize the units."));

    end SIAK;

    record SIAm
      "<html>Base constants and units for SI with <i>k</i><sub>F</sub> and <i>R</i> normalized instead of A and m</html>"

      extends Base(
        final R_inf=sqrt(8.3144621)*10973731.568539,
        final c=299792458/sqrt(8.3144621),
        final R_K=(96485.3365^2*25812.8074434)/8.3144621,
        final k_J=483597.870e9*sqrt(S*s)/m,
        final 'cd'=1);
      annotation (Documentation(info="<html><p>The values of the un-normalized SI base units are (see
  \"FCSys/resources/unit-systems.cdf\"):
  <ul>
  <li>A &asymp; 1.03643e-5</li>
  <li>m &asymp; 0.346803</li>
  </ul></p>

<p>For more information, see the documentation in the
  <a href=\"modelica://FCSys.Units\">Units</a> package.</p></html>"), Commands(
            executeCall=FCSys.Units.setup() "Re-initialize the units."));

    end SIAm;

    record SIAs
      "<html>Base constants and units for SI with <i>k</i><sub>F</sub> and <i>R</i> normalized instead of A and s</html>"
      extends Base(
        final R_inf=10973731.568539,
        final c=299792458/sqrt(8.3144621),
        final R_K=(96485.3365^2*25812.8074434)/sqrt(8.3144621),
        final k_J=483597.870e9*sqrt(S*s)/m,
        final 'cd'=1);
      annotation (Documentation(info="<html><p>The values of the un-normalized SI base units are (see
  \"FCSys/resources/unit-systems.cdf\"):
  <ul>
  <li>A &asymp; 3.59436e-6</li>
  <li>s &asymp; 2.88348</li>
  </ul></p>

  <p>For more information, see the documentation in the
  <a href=\"modelica://FCSys.Units\">Units</a> package.</p></html>"), Commands(
            executeCall=FCSys.Units.setup() "Re-initialize the units."));

    end SIAs;

    record SIKmol
      "<html>Base constants and units for SI with <i>k</i><sub>F</sub> and <i>R</i> normalized instead of K and mol</html>"
      extends Base(
        final R_inf=10973731.568539,
        final c=299792458,
        final R_K=25812.8074434,
        final k_J=483597.870e9*sqrt(S*s)/m,
        final 'cd'=1);
      annotation (Documentation(info="<html><p>The values of the un-normalized SI base units are (see
  \"FCSys/resources/unit-systems.cdf\"):
  <ul>
  <li>K &asymp; 8.61733e-5</li>
  <li>mol &asymp; 96485.3</li>
  </ul></p>

<p>For more information, see the documentation in the
  <a href=\"modelica://FCSys.Units\">Units</a> package.</p></html>"), Commands(
            executeCall=FCSys.Units.setup() "Re-initialize the units."));

    end SIKmol;

    record SIKs
      "<html>Base constants and units for SI with <i>k</i><sub>F</sub> and <i>R</i> normalized instead of K and s</html>"
      extends Base(
        final R_inf=10973731.568539,
        final c=96485.3365*299792458,
        final R_K=96485.3365^3*25812.8074434,
        final k_J=483597.870e9*sqrt(S*s)/m,
        final 'cd'=1);
      annotation (Documentation(info="<html><p>The values of the un-normalized SI base units are (see
  \"FCSys/resources/unit-systems.cdf\"):
  <ul>
  <li>K &asymp; 7.74028e10</li>
  <li>s &asymp; 1.03643e-5</li>
  </ul></p>

<p>For more information, see the documentation in the
  <a href=\"modelica://FCSys.Units\">Units</a> package.</p></html>"), Commands(
            executeCall=FCSys.Units.setup() "Re-initialize the units."));

    end SIKs;

    record SImmol
      "<html>Base constants and units for SI with <i>k</i><sub>F</sub> and <i>R</i> normalized instead of m and mol</html>"
      extends Base(
        final R_inf=sqrt(8.3144621/96485.3365)*10973731.568539,
        final c=sqrt(96485.3365/8.3144621)*299792458,
        final R_K=(96485.3365*25812.8074434)/8.3144621,
        final k_J=483597.870e9*sqrt(S*s)/m,
        final 'cd'=1);
      annotation (Documentation(info="<html><p>The values of the un-normalized SI base units are (see
  \"FCSys/resources/unit-systems.cdf\"):
  <ul>
  <li>m &asymp; 107.724</li>
  <li>mol &asymp; 96485.3</li>
  </ul></p>

<p>For more information, see the documentation in the
  <a href=\"modelica://FCSys.Units\">Units</a> package.</p></html>"), Commands(
            executeCall=FCSys.Units.setup() "Re-initialize the units."));

    end SImmol;

    record SIms
      "<html>Base constants and units for SI with <i>k</i><sub>F</sub> and <i>R</i> normalized instead of m and s</html>"
      extends Base(
        final R_inf=96485.3365*sqrt(8.3144621)*10973731.568539,
        final c=299792458/sqrt(8.3144621),
        final R_K=(96485.3365*25812.8074434)/8.3144621,
        final k_J=483597.870e9*sqrt(S*s)/m,
        final 'cd'=1);
      annotation (Documentation(info="<html><p>The values of the un-normalized SI base units are (see
  \"FCSys/resources/unit-systems.cdf\"):
  <ul>
  <li>m &asymp; 3.59436e-6</li>
  <li>s &asymp; 1.03643e-5</li>
  </ul></p>

<p>For more information, see the documentation in the
  <a href=\"modelica://FCSys.Units\">Units</a> package.</p></html>"), Commands(
            executeCall=FCSys.Units.setup() "Re-initialize the units."));

    end SIms;

    record SImols
      "<html>Base constants and units for SI with <i>k</i><sub>F</sub> and <i>R</i> normalized instead of mol and s</html>"
      extends Base(
        final R_inf=10973731.568539,
        final c=(96485.3365/8.3144621)^(1/3)*299792458,
        final R_K=(96485.3365*25812.8074434)/8.3144621,
        final k_J=483597.870e9*sqrt(S*s)/m,
        final 'cd'=1);
      annotation (Documentation(info="<html><p>The values of the un-normalized SI base units are (see
  \"FCSys/resources/unit-systems.cdf\"):
  <ul>
  <li>mol &asymp; 4261.73</li>
  <li>s &asymp; 0.0441697</li>
  </ul></p>

<p>For more information, see the documentation in the
  <a href=\"modelica://FCSys.Units\">Units</a> package.</p></html>"), Commands(
            executeCall=FCSys.Units.setup() "Re-initialize the units."));

    end SImols;

    record Base "Base constants and units"

      final constant Q.Angle rad=1 "radian";
      // SI unit of rotation or planar angle
      // This condition is required by BIPM [BIPM2006, Table 3].  It can't
      // be relaxed  because BIPM doesn't explicitly use angle in the
      // definitions of Hz, sr, etc. and NIST doesn't explicitly use angle
      // in the relations for R_inf, c_3_nu, etc. [NIST2010].
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
      final constant Q.Number k_F=1
        "<html>Faraday constant (k<sub>F</sub>)</html>";
      // The unit of substance (mole) is inversely proportional to this value.
      // The Faraday constant isn't adjustable because the equations of FCSys
      // require that it's one, which means that charge is considered to be
      // an amount of substance.
      final constant Q.Number R=1 "gas constant";
      // The unit of temperature (kelvin) is inversely proportional to this value.
      // The gas constant isn't adjustable because the equations of FCSys
      // require that it's one, which means that temperature is considered to
      // be a potential.
      annotation (Documentation(info="<html><p>For more information, see the notes in the Modelica code and the documentation in the
  <a href=\"modelica://FCSys.Units\">Units</a> package.</p></html>"), Commands(
            executeCall=FCSys.Units.setup() "Re-initialize the units."));

    end Base;
    annotation (Documentation(info="<html>
  <p><a href=\"modelica://FCSys\">FCSys</a> requires that the Faraday and gas constants are
  normalized to one.  The structure of the <a href=\"modelica://FCSys.Units\">Units</a> package allows
  those constants to be relaxed, but the models in <a href=\"modelica://FCSys\">FCSys</a>
  generally do not.</p>

  <p>Some natural systems of units
  are not compatible with <a href=\"modelica://FCSys\">FCys</a>.
  Since the Faraday and gas constants
  are both normalized, it follows that <code>k_B = q</code>.  This is not
  the case for the Planck, Rydberg, and Natural systems of units
  [<a href=\"http://en.wikipedia.org/wiki/Natural_units\">http://en.wikipedia.org/wiki/Natural_units</a>].</p>

  <p>The quasi-SI
  sets in this package are named by listing (in alphabetical order) the two units that are
  <i>not</i> normalized for the sake of setting the Faraday and gas constants equal to one.
  There are eight possible sets of this type (<a href=\"modelica://FCSys.Units.Bases.AK\">AK</a>,
  <a href=\"modelica://FCSys.Units.Bases.SIAm\">SIAm</a>,
  <a href=\"modelica://FCSys.Units.Bases.SIAs\">SIAs</a>,
  <a href=\"modelica://FCSys.Units.Bases.SIKmol\">SIKmol</a>,
  <a href=\"modelica://FCSys.Units.Bases.SIKs\">SIKs</a>,
  <a href=\"modelica://FCSys.Units.Bases.SImmol\">SImmol</a>
  <a href=\"modelica://FCSys.Units.Bases.SIms\">SIms</a>,
  <a href=\"modelica://FCSys.Units.Bases.SImols\">SImols</a>).</p>

  <p>For more information, see the documentation in the
  <a href=\"modelica://FCSys.Units\">Units</a> package.</p></html>"), Commands(
          executeCall=FCSys.Units.setup() "Re-initialize the units."));

  end Bases;

  function from_degC "Convert from temperature in degree Celsius"
    extends Modelica.SIunits.Conversions.ConversionIcon;

    input Real T_degC "Temperature in degree Celsius";
    output Q.TemperatureAbsolute T "Thermodynamic temperature";

  algorithm
    T := (T_degC + 273.15)*K
      annotation (Inline=true, inverse(T_degC=to_degC(T)));

  end from_degC;

  function to_degC "Convert to temperature in degree Celsius"
    extends Modelica.SIunits.Conversions.ConversionIcon;

    input Q.TemperatureAbsolute T "Thermodynamic temperature";
    output Real T_degC "Temperature in degree Celsius";

  algorithm
    T_degC := T/K - 273.15
      annotation (Inline=true, inverse(T=from_degC(T_degC)));

  end to_degC;

  function from_kPag "Convert from gauge pressure in kilopascals"
    extends Modelica.SIunits.Conversions.ConversionIcon;

    input Real p_kPag "Gauge pressure in kilopascals";
    output Q.PressureAbsolute p "Absolute pressure";

  algorithm
    p := p_kPag*kPa + atm annotation (Inline=true,inverse(p_kPag=to_kPag(p)));

  end from_kPag;

  function to_kPag "Convert to gauge pressure in kilopascals"
    extends Modelica.SIunits.Conversions.ConversionIcon;

    input Q.PressureAbsolute p "Absolute pressure";
    output Real p_kPag "Gauge pressure in kilopascals";

  algorithm
    p_kPag := (p - atm)/kPa
      annotation (Inline=true, inverse(p=from_kPag(p_kPag)));

  end to_kPag;

  // ------------------------------------------------------------------------
  // Mathematical constants
  // ------------------------------------------------------------------------

  final constant Q.Number pi=2*arccos(0) "<html>pi (<i>&pi;</i>)</html>";
  // Circumference per unit diameter
  final constant Q.Number e=exp(1) "Euler number";
  // Natural base

  // ------------------------------------------------------------------------
  // Base physical constants and units
  // ------------------------------------------------------------------------

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

  // ------------------------------------------------------------------------
  // Empirical units
  // ------------------------------------------------------------------------
  // Note:  The values are currently based on the those from [NIST2010].
  // The measured values are used rather than conventional values (where
  // they exist).

  constant Q.Length m=10973731.568539*rad/R_inf "meter";
  // SI unit of length
  // This is the "Rydberg constant" relation [NIST2010].  The unit radian
  // is included to be explicit, although it's currently one by definition
  // [BIPM2006].  The Rydberg constant may be determined by measuring the
  // spectra of hydrogen, deuterium, and antiprotonic helium
  // [http://en.wikipedia.org/wiki/Rydberg_constant].
  constant Q.Time s=299792458*m/c "second";
  // SI unit of time or duration
  // This is the "speed of light in vacuum" relation [NIST2010].  c may
  // be determined (among other ways) by measuring the time for
  // electromagnetic signals to travel to and from spacecraft.
  constant Q.MagneticFlux Wb=483597.870e9/k_J "weber";
  // SI unit of magnetic flux
  // This is the "Josephson constant" relation [NIST2010].  The Josephson
  // constant can be determined by measurements of supercurrent
  // [http://en.wikipedia.org/wiki/Josephson_effect].
  constant Q.ConductanceElectrical S=25812.8074434/R_K "siemen";
  // SI unit of electrical conductance
  // This is the "von Klitzing constant" relation [NIST2010].  The unit
  // radian is included on the denominator for dimensional consistency, but
  // it's one by the current definition [BIPM2006].  The von Klitzing
  // constant may be determined by measuring the quantum hall effect
  // [http://en.wikipedia.org/wiki/Quantum_Hall_effect].
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
  // [http://en.wikipedia.org/wiki/Johnson-Nyquist_noise].

  // ------------------------------------------------------------------------
  // SI base units [BIPM2006, Table 1] and intermediate units
  // ------------------------------------------------------------------------
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

  // ------------------------------------------------------------------------
  // SI prefixes [BIPM2006, Table 5]
  // ------------------------------------------------------------------------

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

  // ------------------------------------------------------------------------
  // Derived units in SI with special names and symbols [BIPM2006, Table 3]
  // ------------------------------------------------------------------------
  // Note:  rad, S, C, Wb, V, J, and Sv have already been defined.  Degree
  // Celsius is only defined in setup(), degC(), and from_degC() since it
  // includes an offset.

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

  // ------------------------------------------------------------------------
  // Non-SI units accepted for use with SI units [BIPM2006, Table 6]
  // ------------------------------------------------------------------------

  final constant Q.Time min=60*s "minute";
  final constant Q.Time hr=60*min "hour";
  final constant Q.Time day=24*hr "day";
  final constant Q.Angle degree=2*pi*rad/360 "<html>degree (&deg;)</html>";
  final constant Q.Volume L=(deci*m)^3 "liter (L or l)";

  // ------------------------------------------------------------------------
  // Derived physical constants
  // ------------------------------------------------------------------------
  // Note:  These are established by definition, but may include
  // transcendental mathematical constants.

  // Electromagnetism
  final constant Q.ConductanceElectrical G_0=2/R_K
    "<html>conductance quantum (<i>G</i><sub>0</sub>)</html>";
  final constant Q.MagneticFlux Phi_0=1/k_J
    "<html>magnetic flux quantum (&Phi;<sup>0</sub>)</html>";
  final constant Q.Amount q=G_0*Phi_0 "elementary charge";
  final constant Q.MomentumRotational h=2*q*Phi_0 "Planck constant";
  // The Planck constant over 2*pi (hbar) isn't included as a unique
  // variable.  The unit of angle (rad or cyc) should be factored into the
  // variable that represents frequency as a quantity.  Then, it's
  // unnecessary to use hbar, e.g.:
  //     hbar = h = 1.0545e-34*J/Hz = 6.6260e-34*J*s/cyc,
  // where Hz = rad/s.  Currently, rad = 1 (see U.Bases.Base).
  final constant Q.Number alpha=pi*1e-7*c*s*G_0/(m*S)
    "<html>fine-structure constant (&alpha;)</html>";
  // The fine-structure constant includes the product of the speed of light
  // in vacuum, expressed in meters per second and conductance quantum,
  // expressed in siemens.
  final constant Q.ResistanceElectrical Z_0=2*R_K*alpha
    "<html>characteristic impedance of vacuum (<i>Z</i><sub>0</sub>)</html>";
  // See  http://en.wikipedia.org/wiki/Characteristic_impedance_of_vacuum.
  final constant Q.Permeability mu_0=Z_0/c
    "<html>magnetic constant (&mu;<sup>0</sub>)</html>";
  // This is also called the vacuum permeability or permeability of free
  // space.
  final constant Q.Permittivity epsilon_0=1/(Z_0*c)
    "<html>electric constant (&epsilon;<sup>0</sub>)</html>";
  // This is also called the vacuum permittivity or permittivity of free
  // space.
  final constant Q.Permeability k_A=mu_0/(4*pi)
    "<html>magnetic force constant (<i>k</i><sub>A</sub>)</html>";
  // The factor of 4*pi is the result of the line integral that is used to
  // derive Ampere's force law
  // [http://en.wikipedia.org/wiki/Ampere's_force_law].
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
  // [http://en.wikipedia.org/wiki/Wien's_displacement_law, accessed
  // 1/19/10].
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
  // integration [http://en.wikipedia.org/wiki/Stefan-Boltzmann_law].  Now,
  // the equation is:
  // B*h*rad*(h*rad*c)^2/(2*(k_B*T)^4) = pi*(h*f/T)^3/(exp(h*f/(k_B*T)) - 1).
  // The integral of the RHS as a function of (h*f/(k_B*T)) over the entire
  // frequency domain (0 to infinity) is pi^4/15.  Finally, this results in:
  // B_tot/T^4 = 2*pi^5*k_B^4/(15*(h*rad)^3*c^2), where the RHS is the
  // Stefan-Boltzmann constant.  Here, the unit rad has been included for
  // dimensional consistency.

  // ------------------------------------------------------------------------
  // Selected other non-SI units from [BIPM2006, Table 8]
  // ------------------------------------------------------------------------
  // Note:  Logarithmic ratios have been excluded because they can't be
  // represented in Dymola's unit conversion GUI.

  final constant Q.Pressure bar=1e5*Pa "bar";
  final constant Q.Length Aring=0.1*nano*m "<html>angstrom (&#8491;)</html>";

  // ------------------------------------------------------------------------
  // Additional units that are useful for fuel cells
  // ------------------------------------------------------------------------

  final constant Q.Pressure atm=101325*Pa "atmosphere";
  // Value from "standard atmosphere" [NIST2010]
  final constant Q.Pressure kPa=kilo*Pa "kilopascal";
  final constant Q.Length cm=centi*m "centimeter";
  final constant Q.Length mm=milli*m "millimeter";
  final constant Q.Number '%'=centi "percent (%)";
  final constant Q.Density M=U.mol/U.L "molar";
  final constant Q.Volume cc=U.cm^3 "cubic centimeter";
  annotation (Documentation(info="<html>
  <p>The <a href=\"modelica://FCSys.Units\">Units</a> package is abbreviated as <code>U</code> for convenience throughout
  the rest of <a href=\"modelica://FCSys.FCSys\">FCSys</a>.  For example, an initial pressure might be defined as
  <code>p_IC = U.atm</code>.</p>

  <p>The information below has been updated from
  [<a href=\"modelica://FCSys.UsersGuide.References\">Davies and Paredis, 2012</a>].  That paper
  also offers suggestions as to how the approach might be better integrated in
  <a href=\"http://www.modelica.org\">Modelica</a>.</p>

<p><b>Overview:</b></p>

<p>Models of physical systems involve variables that represent physical quantities.
As stated by the Bureau International des Poids et Mesures (BIPM)
[<a href=\"modelica://FCSys.UsersGuide.References\">BIPM2006</a>, p. 103]:
<blockquote>
  \"The value of a quantity is generally expressed as the product of a number and a unit.  The
  unit is simply a particular example of the quantity concerned which is used as a reference, and
  the number is the ratio of the value of the quantity to the unit.\"
</blockquote>
In general, a unit may be the product of powers of other units, whether they are base units or
units derived from the base units in the same manner.</p>

<p>In Modelica, a physical quantity is generally expressed as an instance of the <code>Real</code>
type.  Its <code>value</code> attribute is typically the number associated with the value of the
quantity (not the value of the quantity, as will be seen).  The <code>unit</code> attribute is a
string that describes the unit by which the value of the quantity is divided to arrive at the
number.<sup><a href=\"#fn1\" id=\"ref1\">1</a></sup>  The <code>displayUnit</code> attribute (also
a string) describes the unit by which the value should be divided to arrive at the number as it
is entered by or presented to the user.  The <code>Real</code> type contains other attributes as
well, including <code>quantity</code> string.</p>

<p>The <a href=\"modelica://Modelica.SIunits\">SIunits</a> package of the Modelica Standard Library contains types that
extend the <code>Real</code> type.  The type definitions modify the
<code>unit</code>, <code>displayUnit</code>, and <code>quantity</code> attributes (among others)
to represent various physical quantities.  The <code>unit</code> and <code>displayUnit</code>
attributes are based on SI.   The <code>quantity</code> string is generally used to
describe the name of the physical quantity.  For example, the <a href=\"modelica://Modelica.SIunits.Velocity\">Velocity</a> type has
a <code>unit</code> of \"m/s\" and a <code>quantity</code> of
\"Velocity\".</p>

<p>If an instance of <a href=\"modelica://Modelica.SIunits.Velocity\">Velocity</a> has
a <code>value</code> of one (<i>v</i> = 1),
then it is meant that \"the value of velocity is one meter per second.\"  Again, the
<code>value</code> attribute represents the number, or the value divided by the unit, not the
value itself.  This apparent conflict is solved in <a href=\"modelica://FCSys\">FCSys</a> by
establishing units (including the meter and the second) as mathematical entities and writing
<i>v</i> = 1&middot;m/s.  Here, the variable <i>v</i> directly represents the value.
Its <code>value</code> attribute is truly the value in the context of the statement by BIPM.</p>

One advantage is that unit conversion is handled
naturally.  The essence of unit conversion is that the phrase \"value in unit\" mathematically means
\"value divided by unit.\"  Continuing with the previous example, <i>v</i>
is divided by m/s in order to display <i>v</i> in meters per second (as a
number).  If another unit of length like the foot is established by the
appropriate relation (ft &asymp; 0.3048&middot;m) and <i>v</i> is divided by
ft/s, the result is velocity in feet per second (&sim;3.2894).</p>

<p>As another example, frequency is sometimes represented by a variable
in hertz or cycles per second (e.g., &nu;) and other times by a variable in radians
per second (e.g., &omega;).  If the variable represents the value directly, then there
is no need to specify which units it is in.  The units are included; they have not been factored
out by division.  A common variable (e.g., <i>f</i>) can be used in both cases, which
standardizes the equations. The math is equivalent due to the relationships
among units (e.g., 1&middot;cycle = 2&pi;&middot;rad).</p>

<p><b>Method:</b></p>

<p>Each unit is represented by a constant.  The values of most units is derived from
the values of other units and constants (e.g., 1&middot;cycle = 2&pi;&middot;rad).  However,
there are several units (in SI, 7) that are independent.  These
base units are established by the \"particular example of the quantity
concerned which is used as a reference\" quoted previously
[<a href=\"modelica://FCSys.UsersGuide.References\">BIPM2006</a>].  The designation of \"base\"
or \"derived\" is somewhat arbitrary [<a href=\"modelica://FCSys.UsersGuide.References\">Fritzson2004</a>, p. 375]
but regardless, there are a number of units that must be defined by example.  Considering only
the immediate physical system, these units are linearly independent.</p>

<p>If only SI will be used, then it is easiest to strictly set each of the base units of
SI equal to one&mdash;the meter (m), kilogram (kg), second (s), ampere (A),
kelvin (K), mole (mol), and candela (cd).  This is implicitly the case in
<code>Modelica.SIunits</code>, but again, it hardly captures the idea that a value is the
product of a number and a unit.</p>

<p>Instead, most of the base units are established by universal physical constants.
The \"particular example of the quantity concerned which is used as a reference\"
[<a href=\"modelica://FCSys.UsersGuide.References\">BIPM2006</a>] is an experiment that yields
precise and repeatable results in determining a constant rather than a prototype which is
carefully controlled and distributed via replicas.
This method of defining the base units from measured physical quantities (rather than
vice versa) is natural and reflects the way that standards organizations (e.g., NIST) define units.
A system of units is considered to be natural if
all of its base units are established by physical constants
[<a href=\"http://en.wikipedia.org/wiki/Natural_units\">http://en.wikipedia.org/wiki/Natural_units</a>].
Often, those physical constants are defined to be equal to one, but the values can be chosen to scale
the numerical values of variables during simulation.</p>

<p>There are systems where typical values are many orders of magnitude larger or smaller than the
related product of powers of base SI units (e.g., the domains of astrophysics and atomic
physics).  In modeling and simulating those systems, it may be advantageous to choose
appropriately small or large values (respectively) for the corresponding base units such that the
product of the number (large or small in magnitude) and the unit (small or large, respectively)
is well-scaled.  Products of this type are often involved in initial conditions or parameter
expressions, which are not time-varying.  Therefore, the number and the unit can be multiplied
before the dynamic simulation.  During the simulation, only the value is important.  After the
simulation, the trajectory of the value may be divided by the unit for display.  This scaling is
usually unnecessary due to the wide range and appropriate distribution of the real numbers that
are representable in floating point space.  The Modelica language specification recommends that
floating point numbers be represented in at least IEEE double precision, which covers magnitudes
from &sim;2.225&times;10<sup>-308</sup> to &sim;1.798&times;10<sup>308</sup>
[<a href=\"modelica://FCSys.UsersGuide.References\">Modelica3.2</a>, p. 13].
However, in some cases it may be preferable to carefully scale the units and use single
precision instead for the sake of computational performance.  There are fields of research where,
even today, simulations are sometimes performed in single precision
[<a href=\"modelica://FCSys.UsersGuide.References\">Brown2011</a>,
<a href=\"modelica://FCSys.UsersGuide.References\">Hess2008</a>]
and where scaling is a concern
[<a href=\"modelica://FCSys.UsersGuide.References\">Rapaport2004</a>, p. 29].</p>

<p>The method is neutral
with regards to not only the values of the base units, but also the choice of the base units and
even the number of base units.  This is an advantage because many systems of units are used besides
SI. As mentioned previously, the choice of base units is somewhat
arbitrary, and different systems of units are based on different choices.  Some systems of units
have fewer base units (lower rank) than SI, since additional constraints are added that
exchange base units for derived units.  For example, the Planck, Stoney, Hartree, and Rydberg
systems of units define the Boltzmann constant to be equal to one (<i>k</i><sub>B</sub> = 1)
[<a href=\"http://en.wikipedia.org/wiki/Natural_units\">http://en.wikipedia.org/wiki/Natural_units</a>].
The unit K is eliminated
[<a href=\"modelica://FCSys.UsersGuide.References\">Greiner1995</a>, p. 386]
or, more precisely, considered a derived unit instead of a base unit.  In SI, the
kelvin would be derived from the units kilogram, meter, and second (K
&asymp; 1.381&times;10<sup>-23</sup>&middot;kg&middot;m<sup>2</sup>/s<sup>2</sup>).</p>

    <p>There are six independent constants or units in the <a href=\"modelica://FCSys.Units\">Units</a> package (see
    <a href=\"modelica://FCSys.Units.Bases\">Units.Bases</a>),
    but SI has seven independent base units (m, kg, s, A, K, mol, and cd).
    In <a href=\"modelica://FCSys\">FCSys</a>, two additional constraints are imposed in order
    to simplify the model equations and allow electrons and chemical species to be to represented by the
    same base <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.
    First, the Faraday constant (k<sub>F</sub> or 96485.3399&middot;C/mol)
    is normalized to one. This implies that the mole (mol) is proportional to the coulomb
    (C or Wb&middot;S), which is considered a number of reference particles given the charge number.
    Also, the gas constant (R or 8.314472&middot;J/(mol&middot;K)) is normalized to one.
    Therefore, the kelvin (K) is proportional to the volt
    (V or J/C). In addition, the radian (rad) is defined as a base constant.
    However, it must be set equal to one in the current version of the International System of Units (SI)
    [<a href=\"modelica://FCSys.UsersGuide.References\">BIPM2006</a>].</p>

<p><b>Implementation:</b><p><p>The units and constants are defined as variables in this
<a href=\"modelica://FCSys.Units\">Units</a> package.  Each is a <code>constant</code> of
the appropriate type from the <a href=\"modelica://FCSys.Quantities\">Quantities</a> package. The
first section of the Modelica definition of this package establishes mathematical constants.  The next
 section establishes the base constants and units, which grouped in a replaceable subpackage.  The third section
 establishes the constants and units which may be derived from the base units and constants using
 accepted empirical relations.  The rest of the code establishes the SI prefixes
 and the remaining derived units and constants.  The SI prefixes are included in their
 unabbreviated form in order to avoid name conflicts.  All of the primary units of SI
 are included (Tables 1 and 3 of
 [<a href=\"modelica://FCSys.UsersGuide.References\">BIPM2006</a>]) except for &deg;C, since
 it involves an offset.  Other convenient units are included for the system at hand (e.g.,
 atm).</p>

<p>The <a href=\"modelica://FCSys.Units.setup\">Units.setup</a> function establishes unit conversions
using the values of the units, constants, and prefixes.  These unit conversions may include offsets.
The function also sets the default display units.  It is automatically called when
<a href=\"modelica://FCSys\">FCSys</a> is
loaded from the \"FCSys/load.mos\" script.  It can also be called manually from the
\"Re-initialize the units\" command available in Dymola from the
<a href=\"modelica://FCSys.Units\">Units</a> package or any subpackage.  A spreadsheet
(\"FCSys/resources/quantities.xls\") is available to help
maintain the quantities, default units, and the setup function.</p>

<p>The values of the units, constants, and prefixes can be evaluated by translating the
<a href=\"modelica://FCSys.Units.Examples.Evaluate\">Units.Examples.Evaluate</a> model.  This
defines the values in the Dymola workspace.  For convenience, the \"FCSys/load.mos\" script automatically
does this and saves the result as \"units.mos\" in the working directory.</p>

  <p>This package also contains functions (e.g., <a href=\"modelica://FCSys.Units.to_degC\">to_degC</a>) that
  convert quantities from the unit system defined in <a href=\"modelica://FCSys\">FCSys</a> to quantities
  expressed in units.  Functions are
  included for units that involve offsets<!-- or other functions besides simple scaling-->.
  For conversions that require just a scaling factor, it is simpler to use the
  units directly.  For example, to convert from potential in volts use <code>v = v_V*V</code>,
  where <code>v</code> is potential and <code>v_V</code> is potential expressed in volts.</p>

<p>An instance of the <a href=\"modelica://FCSys.Conditions.Environment\">Environment</a> model is usually included
at the top level of a model.  It records the base units or constants so that it is possible to re-derive
all of the other units and constants.  This is important in order to properly interpret simulation results if the
base units or constants are later re-adjusted.</p>

<p>Where the <code>der</code> operator is used in models, it is explicitly divided by the unit second
(e.g., <code>der(x)/s</code>).  This is necessary because the global variable <code>time</code>
is in seconds.</p>

<p>In theory, standard Modelica unit checking tools may be used to check the dimensions of equations
in <a href=\"modelica://FCSys\">FCSys</a>.</p>

<p>Some units are defined that include prefixes (e.g., kg, mm, and kPa).  However,
most prefixes must be given as explicit factors (e.g., <code>kilo*m</code>).</p>

  <p>Although it is not necessary in <a href=\"http://www.modelica.org\">Modelica</a>, the declarations
  in this package are presorted so that they can be easily ported to imperative or causal languages (e.g.,
  <a href=\"http://www.python.org\">Python</a>, C).</p>

<hr>

    <small>
    <p id=\"fn1\">1. Hereafter, the value of the quantity is referred to as simply the value, but
    it should not be confused with the <code>value</code> attribute (which, as of version 3.3 of the
    Modelica specification, is the number).<a href=\"#ref1\" title=\"Jump back to footnote 1 in the text.\">&#8629;</a></p>

<p><b>Licensed by the Georgia Tech Research Corporation under the Modelica License 2</b><br>
  Copyright 2007&ndash;2012, Georgia Tech Research Corporation.</p>

<p><i>This Modelica package is <u>free</u> software and the use is completely at <u>your own risk</u>;
  it can be redistributed and/or modified under the terms of the Modelica License 2. For license conditions (including the
  disclaimer of warranty) see <a href=\"modelica://FCSys.UsersGuide.ModelicaLicense2\">
  FCSys.UsersGuide.ModelicaLicense2</a> or visit <a href=\"http://www.modelica.org/licenses/ModelicaLicense2\">
  http://www.modelica.org/licenses/ModelicaLicense2</a>.</i></p></html>"),
      Commands(executeCall=FCSys.Units.setup() "Re-initialize the units."));

end Units;

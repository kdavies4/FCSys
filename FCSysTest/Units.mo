within FCSysTest;
package Units
  extends Modelica.Icons.Package;
  function callAll
    "<html>Call all of the test functions for the <a href=\"modelica://FCSys.Units\">Units</a> package</html>"
    import Modelica.Utilities.Streams.print;
    extends Modelica.Icons.Function;

    input String logFile="FCSysTestLog.txt" "Filename where the log is stored";
    output Boolean ok "true, if all tests passed";

  algorithm
    print("--- Test FCSys.Units");
    print("--- Test FCSys.Units", logFile);
    ok := testValues(logFile) and testConversions(logFile);
    annotation (Documentation(info="<html><p>This function call will fail if any of the functions return an
  incorrect result.  It will return <code>true</code> if all of the functions pass.
  There are no inputs.</p></html>"));
  end callAll;

  function testValues
    "<html>Test the values of units and constants in the <a href=\"modelica://FCSys.Units\">Units</a> package</html>"
    import Modelica.Utilities.Streams.print;
    import FCSys.Units.*;
    extends Modelica.Icons.Function;

    input String logFile="FCSysTestLog.txt" "Filename where the log is stored";
    output Boolean ok "true, if all tests passed";

  protected
    function test = Test.assertValue (final expected=1);

  algorithm
    print("... Test of the values of constants and units");
    print("... Test of the values of constants and units", logFile);

    // ----------------------------------------------------------------------
    // Set 1:  Mathematical constants and relations
    test(pi/3.14159265358979323846264338327950288419716939937510, name=
      "1 in set 1");
    // Value from http://en.wikipedia.org/wiki/Pi#Approximate_value
    test(2*pi*rad/(360*degree), name="3 in set 1");
    test('%'/0.01, name="4 in set 1");
    // ----------------------------------------------------------------------
    // Set 2:  Relations from [BIPM2006]
    test(q/(1.602176487e-19*C), name="1 in set 2");
    test(C*V/J, name="2 in set 2");
    // ----------------------------------------------------------------------
    // Set 3: Coherent derived units in the SI with special names and
    // symbols [BIPM2006]
    test(sr, name="1 in set 3");
    test(Hz/(cyc/s), name="2 in set 3");
    // BIPM implicitly assumes that the unit cycle is 1, but the unit rad is 1
    // too---a discrepancy.
    test(N/(kg*m/s^2), name="3 in set 3");
    test(Pa/(N/m^2), name="4 in set 3");
    test(J/(N*m), name="5 in set 3");
    test(W/(J/s), name="6 in set 3");
    test(C/(A*s), name="7 in set 3");
    test(V/(W/A), name="8 in set 3");
    test(F/(C/V), name="9 in set 3");
    test(ohm/(V/A), name="10 in set 3");
    test(S/(A/V), name="11 in set 3");
    test(Wb/(V*s), name="12 in set 3");
    test(T/(Wb/m^2), name="13 in set 3");
    test(H/(Wb/A), name="14 in set 3");
    test(lm/('cd'*sr), name="15 in set 3");
    test(lx/(lm/m^2), name="16 in set 3");
    test(Bq/(cyc/s), name="17 in set 3");
    // See Hz.
    test(Gy/(J/kg), name="18 in set 3");
    test(Sv/(J/kg), name="19 in set 3");
    test(kat/(mol/s), name="20 in set 3");
    // ----------------------------------------------------------------------
    // Set 4:  Relations from [NIST2010]
    // Generated from Resources/NIST.xls, 2013-1-23
    test(1/alpha/137.035999074, name="inverse fine-structure constant");
    test(1/G_0/(12906.4037217*ohm), name="inverse of conductance quantum");
    test(1/m*h*c/(1.239841930e-6*eV), name=
      "inverse meter-electron volt relationship");
    test(1/m*h*c/(4.556335252755e-8*E_h), name=
      "inverse meter-hartree relationship");
    test(1/m*h*c/(1.986445684e-25*J), name="inverse meter-joule relationship");
    test(1/m*h*c/k_B/(1.4387770e-2*K), name="inverse meter-kelvin relationship");
    test(1/m*h/c/(2.210218902e-42*kg), name=
      "inverse meter-kilogram relationship");
    test(100*kPa/(k_B*273.15*K)/(2.6516462e25/m^3), name=
      "Loschmidt constant (273.15 K, 100 kPa)");
    test(101.325*kPa/(k_B*273.15*K)/(2.6867805e25/m^3), name=
      "Loschmidt constant (273.15 K, 101.325 kPa)");
    test(12*g/mol/(12e-3*kg/mol), name="molar mass of carbon-12");
    test(alpha/7.2973525698e-3, name="fine-structure constant");
    test(atm/(101325*Pa), name="standard atmosphere");
    test(c/(299792458*m/s), name="speed of light in vacuum");
    test(c_1/(3.74177153e-16*W*m^2), name="first radiation constant");
    test(c_2/(1.4387770e-2*m*K), name="second radiation constant");
    test(c_3_f*cyc/(5.8789254e10*Hz/K), name=
      "Wien frequency displacement law constant");
    // Factor of cyc due to inconsistencies b/w rad and cyc in [BIPM2006]
    test(c_3_lambda/(2.8977721e-3*m*K), name=
      "Wien wavelength displacement law constant");
    test(cyc/m*c/(299792458*Hz), name="inverse meter-hertz relationship");
    test(E_h/(4.35974434e-18*J), name="Hartree energy");
    test(E_h/(27.21138505*eV), name="Hartree energy in eV");
    test(E_h/(27.21138505*eV), name="hartree-electron volt relationship");
    test(E_h/(4.35974434e-18*J), name="hartree-joule relationship");
    test(E_h/(h*c)/(2.194746313708e7/m), name=
      "hartree-inverse meter relationship");
    test(E_h/c^2/(4.85086979e-35*kg), name="hartree-kilogram relationship");
    test(E_h/h*cyc/(6.579683920729e15*Hz), name="hartree-hertz relationship");
    // Factor of cyc due to inconsistencies b/w rad and cyc in [BIPM2006]
    test(E_h/k_B/(3.1577504e5*K), name="hartree-kelvin relationship");
    test(epsilon_0/(8.854187817e-12*F/m), name="electric constant");
    test(eV/(1.602176565e-19*J), name="electron volt");
    test(eV/(3.674932379e-2*E_h), name="electron volt-hartree relationship");
    test(eV/(1.602176565e-19*J), name="electron volt-joule relationship");
    test(eV/(h*c)/(8.06554429e5/m), name=
      "electron volt-inverse meter relationship");
    test(eV/c^2/(1.782661845e-36*kg), name=
      "electron volt-kilogram relationship");
    test(eV/h*cyc/(2.417989348e14*Hz), name="electron volt-hertz relationship");
    // Factor of cyc due to inconsistencies b/w rad and cyc in [BIPM2006]
    test(eV/k_B/(1.1604519e4*K), name="electron volt-kelvin relationship");
    test(G_0/(7.7480917346e-5*S), name="conductance quantum");
    test(g/mol/(1e-3*kg/mol), name="molar mass constant");
    test(h/(6.62606957e-34*J*s), name="Planck constant");
    test(h/(4.135667516e-15*eV*s), name="Planck constant in eV s");
    test(h*c/(2*pi)/(197.3269718*mega*eV*femto*m), name=
      "Planck constant over 2 pi times c in MeV fm");
    test(h/(2*pi)/(1.054571726e-34*J*s), name="Planck constant over 2 pi");
    test(h/(2*pi)/(6.58211928e-16*eV*s), name=
      "Planck constant over 2 pi in eV s");
    test(Hz*h/cyc/(4.135667516e-15*eV), name="hertz-electron volt relationship");
    // Factor of 1/cyc due to inconsistencies b/w rad and cyc in [BIPM2006]
    test(Hz*h/cyc/(1.5198298460045e-16*E_h), name="hertz-hartree relationship");
    // Factor of 1/cyc due to inconsistencies b/w rad and cyc in [BIPM2006]
    test(Hz*h/cyc/(6.62606957e-34*J), name="hertz-joule relationship");
    // Factor of 1/cyc due to inconsistencies b/w rad and cyc in [BIPM2006]
    test(Hz*h/c^2/cyc/(7.37249668e-51*kg), name="hertz-kilogram relationship");
    // Factor of 1/cyc due to inconsistencies b/w rad and cyc in [BIPM2006]
    test(Hz*h/k_B/cyc/(4.7992434e-11*K), name="hertz-kelvin relationship");
    // Factor of 1/cyc due to inconsistencies b/w rad and cyc in [BIPM2006]
    test(Hz/c/cyc/(3.335640951e-9/m), name="hertz-inverse meter relationship");
    // Factor of 1/cyc due to inconsistencies b/w rad and cyc in [BIPM2006]
    test(J/(6.24150934e18*eV), name="joule-electron volt relationship");
    test(J/(2.29371248e17*E_h), name="joule-hartree relationship");
    test(J/(h*c)/(5.03411701e24/m), name="joule-inverse meter relationship");
    test(J/c^2/(1.112650056e-17*kg), name="joule-kilogram relationship");
    test(J/h*cyc/(1.509190311e33*Hz), name="joule-hertz relationship");
    test(J/k_B/(7.2429716e22*K), name="joule-kelvin relationship");
    test(k_B/(1.3806488e-23*J/K), name="Boltzmann constant");
    test(k_B/(8.6173324e-5*eV/K), name="Boltzmann constant in eV/K");
    test(k_B/(h*c)/(69.503476/(m*K)), name=
      "Boltzmann constant in inverse meters per kelvin");
    test(k_B/h*cyc/(2.0836618e10*Hz/K), name="Boltzmann constant in Hz/K");
    // Factor of cyc due to inconsistencies b/w rad and cyc in [BIPM2006]
    test(k_F/(96485.3365*C/mol), name="Faraday constant");
    test(k_J*cyc/(483597.9e9*Hz/V), name=
      "conventional value of Josephson constant");
    test(k_J*cyc/(483597.870e9*Hz/V), name="Josephson constant");
    // Factor of cyc due to inconsistencies b/w rad and cyc in [BIPM2006]
    test(K*k_B/(8.6173324e-5*eV), name="kelvin-electron volt relationship");
    test(K*k_B/(3.1668114e-6*E_h), name="kelvin-hartree relationship");
    test(K*k_B/(1.3806488e-23*J), name="kelvin-joule relationship");
    test(K*k_B/(h*c)/(69.503476/m), name="kelvin-inverse meter relationship");
    test(K*k_B/c^2/(1.5361790e-40*kg), name="kelvin-kilogram relationship");
    test(K*k_B/h*cyc/(2.0836618e10*Hz), name="kelvin-hertz relationship");
    // Factor of cyc due to inconsistencies b/w rad and cyc in [BIPM2006]
    test(kg*c^2/(5.60958885e35*eV), name="kilogram-electron volt relationship");
    test(kg*c^2/(2.061485968e34*E_h), name="kilogram-hartree relationship");
    test(kg*c^2/(8.987551787e16*J), name="kilogram-joule relationship");
    test(kg*c^2/h*cyc/(1.356392608e50*Hz), name="kilogram-hertz relationship");
    // Factor of cyc due to inconsistencies b/w rad and cyc in [BIPM2006]
    test(kg*c^2/k_B/(6.5096582e39*K), name="kilogram-kelvin relationship");
    test(kg*c/h/(4.52443873e41/m), name="kilogram-inverse meter relationship");
    test(mu_0/(12.566370614e-7*N/A^2), name="mag. constant");
    test(N_A/(6.02214129e23/mol), name="Avogadro constant");
    test(N_A*h/(3.9903127176e-10*J*s/mol), name="molar Planck constant");
    test(N_A*h*c/(0.119626565779*J*m/mol), name="molar Planck constant times c");
    test(Phi_0/(2.067833758e-15*Wb), name="mag. flux quantum");
    test(q/(1.602176565e-19*C), name="elementary charge");
    test(q/h/(2.417989348e14*A/J), name="elementary charge over h");
    test(R/(8.3144621*J/(mol*K)), name="molar gas constant");
    test(R_inf/(10973731.568539/m), name="Rydberg constant");
    test(R_inf*c*cyc/(3.289841960364e15*Hz), name=
      "Rydberg constant times c in Hz");
    // Factor of cyc due to inconsistencies b/w rad and cyc in [BIPM2006]
    test(R_inf*h*c/(13.60569253*eV), name="Rydberg constant times hc in eV");
    test(R_inf*h*c/(2.179872171e-18*J), name="Rydberg constant times hc in J");
    test(R_K/(25812.807*ohm), name=
      "conventional value of von Klitzing constant");
    test(R_K/(25812.8074434*ohm), name="von Klitzing constant");
    test(R*273.15*K/(100*kPa)/(22.710953e-3*m^3/mol), name=
      "molar volume of ideal gas (273.15 K, 100 kPa)");
    test(R*273.15*K/(101.325*kPa)/(22.413968e-3*m^3/mol), name=
      "molar volume of ideal gas (273.15 K, 101.325 kPa)");
    test(sigma/(5.670373e-8*W/(m^2*K^4)), name="Stefan-Boltzmann constant");
    test(Z_0/(376.730313461*ohm), name="characteristic impedance of vacuum");

    ok := true;
    annotation (Documentation(info="<html><p>This function call will fail if any of the functions return an
  incorrect result.  It will return <code>true</code> if all of the functions pass.
  There are no inputs.</p></html>"));
  end testValues;

  function testConversions
    "<html>Test the unit conversions in the <a href=\"modelica://FCSys.Units\">Units</a> package</html>"
    import Modelica.Utilities.Streams.print;
    import FCSysTest.Test.assertValue;
    import FCSys.Units.*;
    extends Modelica.Icons.Function;

    input String logFile="FCSysTestLog.txt" "Filename where the log is stored";
    output Boolean ok "true, if all tests passed";

  algorithm
    print("... Test of conversion functions");
    print("... Test of conversion functions", logFile);

    // "From" functions
    assertValue(
        from_degC(100),
        373.15*FCSys.Units.K,
        name="from_degC");
    assertValue(
        from_kPag(101.325),
        2*FCSys.Units.atm,
        name="from_kPag");
    // Inverses
    assertValue(
        FCSys.Units.to_degC(from_degC(1)),
        1,
        name="degC");
    assertValue(
        FCSys.Units.to_kPag(from_kPag(1)),
        1,
        name="kPag");

    ok := true;

    annotation (Documentation(info="<html><p>This function call will fail if any of the functions return an
  incorrect result.  It will return <code>true</code> if all of the functions pass.
  There are no inputs.</p></html>"));
  end testConversions;

end Units;

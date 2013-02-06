within FCSys;
package Characteristics
  "Packages with data and functions to correlate physical properties"

  import Modelica.Media.IdealGases.Common.FluidData;
  import Modelica.Media.IdealGases.Common.SingleGasesData;

  extends Modelica.Icons.MaterialPropertiesPackage;

  package Examples "Examples"
    extends Modelica.Icons.ExamplesPackage;

    model Correlations
      "Evaluate the material properties over varying temperature and pressure"
      extends Modelica.Icons.Example;

      // Data
      import 'DataC+' = FCSys.Characteristics.'C+'.Graphite;
      import 'DataC19HF37O5S-' = FCSys.Characteristics.'C19HF37O5S-'.Ionomer;
      import 'Datae-' = FCSys.Characteristics.'e-'.Graphite;
      import DataH2 = FCSys.Characteristics.H2.Gas;

    protected
      package DataH2IG = FCSys.Characteristics.H2.Gas (b_v=[1], specVolPow={-1,
              0}) "H2 as ideal gas";
      import DataH2O = FCSys.Characteristics.H2O.Gas;
      import DataH2OLiquid = FCSys.Characteristics.H2O.Liquid;
      import 'DataH+' = FCSys.Characteristics.'H+'.Ionomer;
      import DataN2 = FCSys.Characteristics.N2.Gas;
      import DataO2 = FCSys.Characteristics.O2.Gas;

      // Conditions

    public
      Connectors.RealOutputInternal T(unit="l2.m/(N.T2)",displayUnit="K")
        "Temperature" annotation (Placement(transformation(extent={{50,16},{70,
                36}}), iconTransformation(extent={{-10,16},{10,36}})));
      Connectors.RealOutputInternal p(unit="m/(l.T2)") "Pressure" annotation (
          Placement(transformation(extent={{50,-36},{70,-16}}),
            iconTransformation(extent={{-10,-36},{10,-16}})));

      // Results
      // -------
      // Isobaric specific heat capacity
      output Q.CapacityThermalSpecific 'c_p_C+'='DataC+'.c_p(T, p);
      output Q.CapacityThermalSpecific 'c_p_C19HF37O5S-'='DataC19HF37O5S-'.c_p(
          T, p);
      output Q.CapacityThermalSpecific 'c_p_e-'='Datae-'.c_p(T, p);
      output Q.CapacityThermalSpecific 'c_p_H+'='DataH+'.c_p(T, p);
      output Q.CapacityThermalSpecific c_p_H2=DataH2.c_p(T, p);
      output Q.CapacityThermalSpecific c_p_H2O=DataH2O.c_p(T, p);
      output Q.CapacityThermalSpecific c_p_H2O_liquid=DataH2OLiquid.c_p(T, p);
      output Q.CapacityThermalSpecific c_p_N2=DataN2.c_p(T, p);
      output Q.CapacityThermalSpecific c_p_O2=DataO2.c_p(T, p);
      //
      // Isochoric specific heat capacity
      output Q.CapacityThermalSpecific 'c_v_C+'='DataC+'.c_v(T, p);
      output Q.CapacityThermalSpecific 'c_v_C19HF37O5S-'='DataC19HF37O5S-'.c_v(
          T, p);
      output Q.CapacityThermalSpecific 'c_v_e-'='Datae-'.c_v(T, p);
      output Q.CapacityThermalSpecific 'c_v_H+'='DataH+'.c_v(T, p);
      output Q.CapacityThermalSpecific c_v_H2=DataH2.c_v(T, p);
      output Q.CapacityThermalSpecific c_v_H2O=DataH2O.c_v(T, p);
      output Q.CapacityThermalSpecific c_v_H2O_liquid=DataH2OLiquid.c_v(T, p);
      output Q.CapacityThermalSpecific c_v_N2=DataN2.c_v(T, p);
      output Q.CapacityThermalSpecific c_v_O2=DataO2.c_v(T, p);
      //
      // Gibbs potential
      output Q.Potential 'g_C+'='DataC+'.g(T, p);
      output Q.Potential 'g_C19HF37O5S-'='DataC19HF37O5S-'.g(T, p);
      output Q.Potential 'g_e-'='Datae-'.g(T, p);
      output Q.Potential 'g_H+'='DataH+'.g(T, p);
      output Q.Potential g_H2=DataH2.g(T, p);
      output Q.Potential g_H2O=DataH2O.g(T, p);
      output Q.Potential g_H2O_liquid=DataH2OLiquid.g(T, p);
      output Q.Potential g_N2=DataN2.g(T, p);
      output Q.Potential g_O2=DataO2.g(T, p);
      //
      // Specific enthalpy
      output Q.Potential 'h_C+'='DataC+'.h(T);
      output Q.Potential 'h_C19HF37O5S-'='DataC19HF37O5S-'.h(T);
      output Q.Potential 'h_e-'='Datae-'.h(T);
      output Q.Potential 'h_H+'=FCSys.Characteristics.'H+'.Ionomer.h(T);
      output Q.Potential h_H2=DataH2.h(T);
      output Q.Potential h_H2O=DataH2O.h(T);
      output Q.Potential h_H2O_liquid=DataH2OLiquid.h(T);
      output Q.Potential h_N2=DataN2.h(T);
      output Q.Potential h_O2=DataO2.h(T);
      //
      // Pressure (indirectly via v_Tp())
      output Q.PressureAbsolute 'p_C+'='DataC+'.p_Tv(T, v_C) if 'DataC+'.isCompressible;
      output Q.PressureAbsolute 'p_C19HF37O5S-'='DataC19HF37O5S-'.p_Tv(T,
          'v_C19HF37O5S-') if 'DataC19HF37O5S-'.isCompressible;
      output Q.PressureAbsolute 'p_e-'='Datae-'.p_Tv(T, 'v_e-') if 'Datae-'.isCompressible;
      output Q.PressureAbsolute 'p_H+'='DataH+'.p_Tv(T, 'v_H+') if 'DataH+'.isCompressible;
      output Q.PressureAbsolute p_H2=DataH2.p_Tv(T, v_H2) if DataH2.isCompressible;
      output Q.PressureAbsolute p_H2O=DataH2O.p_Tv(T, v_H2O) if DataH2O.isCompressible;
      // Note that p_H2O diverges from p in Dymola 7.4 due to the large
      // coefficients in the second row of DataH2O.b_v, which cause numerical
      // errors.
      output Q.PressureAbsolute p_H2O_liquid=DataH2OLiquid.p_Tv(T, v_H2O_liquid)
        if DataH2OLiquid.isCompressible;
      output Q.PressureAbsolute p_N2=DataN2.p_Tv(T, v_N2) if DataN2.isCompressible;
      output Q.PressureAbsolute p_O2=DataO2.p_Tv(T, v_O2) if DataO2.isCompressible;
      output Q.PressureAbsolute p_IG=DataH2IG.p_Tv(T, v_IG) if DataH2IG.isCompressible
        "Pressure of ideal gas";
      //
      // Specific entropy
      output Q.Number 's_C+'='DataC+'.s(T, p);
      output Q.Number 's_C19HF37O5S-'='DataC19HF37O5S-'.s(T, p);
      output Q.Number 's_e-'='Datae-'.s(T, p);
      output Q.Number 's_H+'=FCSys.Characteristics.'H+'.Ionomer.s(T, p);
      output Q.Number s_H2=DataH2.s(T, p);
      output Q.Number s_H2O=DataH2O.s(T, p);
      output Q.Number s_H2O_liquid=DataH2OLiquid.s(T, p);
      output Q.Number s_N2=DataN2.s(T, p);
      output Q.Number s_O2=DataO2.s(T, p);
      //
      // Specific volume
      output Q.VolumeSpecificAbsolute 'v_C+'='DataC+'.v_Tp(T, p);
      output Q.VolumeSpecificAbsolute 'v_C19HF37O5S-'='DataC19HF37O5S-'.v_Tp(T,
          p);
      output Q.VolumeSpecificAbsolute 'v_e-'='Datae-'.v_Tp(T, p);
      output Q.VolumeSpecificAbsolute 'v_H+'='DataH+'.v_Tp(T, p);
      output Q.VolumeSpecificAbsolute v_H2=DataH2.v_Tp(T, p);
      output Q.VolumeSpecificAbsolute v_H2O=DataH2O.v_Tp(T, p);
      output Q.VolumeSpecificAbsolute v_H2O_liquid=DataH2OLiquid.v_Tp(T, p);
      output Q.VolumeSpecificAbsolute v_N2=DataN2.v_Tp(T, p);
      output Q.VolumeSpecificAbsolute v_O2=DataO2.v_Tp(T, p);
      output Q.VolumeSpecificAbsolute v_IG=DataH2IG.v_Tp(T, p)
        "Specific volume of ideal gas";

    protected
      Modelica.Blocks.Sources.Clock clock
        annotation (Placement(transformation(extent={{-70,-10},{-50,10}})));
      Modelica.Blocks.Sources.Constant temperatureOffset(k=273.15*U.K)
        annotation (Placement(transformation(extent={{-20,40},{0,60}})));
      Modelica.Blocks.Sources.Constant pressureOffset(k=U.atm)
        annotation (Placement(transformation(extent={{-20,-60},{0,-40}})));
      Modelica.Blocks.Math.Gain temperatureGain(k=200*U.K)
        annotation (Placement(transformation(extent={{-20,10},{0,30}})));
      Modelica.Blocks.Math.Gain pressureGain(k=0*U.bar)
        annotation (Placement(transformation(extent={{-20,-30},{0,-10}})));
      Modelica.Blocks.Math.Add addTemperature
        annotation (Placement(transformation(extent={{20,16},{40,36}})));
      Modelica.Blocks.Math.Add addPressure
        annotation (Placement(transformation(extent={{20,-36},{40,-16}})));

    equation
      connect(clock.y, temperatureGain.u) annotation (Line(
          points={{-49,6.10623e-16},{-30,6.10623e-16},{-30,20},{-22,20}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(clock.y, pressureGain.u) annotation (Line(
          points={{-49,6.10623e-16},{-30,6.10623e-16},{-30,-20},{-22,-20}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(temperatureOffset.y, addTemperature.u1) annotation (Line(
          points={{1,50},{10,50},{10,32},{18,32}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(pressureOffset.y, addPressure.u2) annotation (Line(
          points={{1,-50},{10,-50},{10,-32},{18,-32}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(temperatureGain.y, addTemperature.u2) annotation (Line(
          points={{1,20},{18,20}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(pressureGain.y, addPressure.u1) annotation (Line(
          points={{1,-20},{18,-20}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(addTemperature.y, T) annotation (Line(
          points={{41,26},{60,26}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(addPressure.y, p) annotation (Line(
          points={{41,-26},{60,-26}},
          color={0,0,127},
          smooth=Smooth.None));
      annotation (
        experiment,
        Commands(file=
              "resources/scripts/Dymola/Characteristics.Examples.Correlations.mos"),

        Icon(graphics));
    end Correlations;

  end Examples;

  package 'C+' "<html>C<sup>+</sup></html>"
    extends Modelica.Icons.Package;
    package Graphite "C+ graphite"

      extends BaseClasses.Characteristic(
        final formula="C+",
        final phase="graphite",
        p0=U.atm,
        specVolPow={0,0},
        b_v=[U.cc*m/(2.2210*U.g)],
        p_min=-Modelica.Constants.inf,
        m=12.0107000*U.g/U.mol - 'e-'.Graphite.m,
        Deltah0_f=0*U.J/U.mol,
        Deltah0=1053.500*U.J/U.mol,
        T_lim_c={200.000,600.000,2000.000,6000.000}*U.K,
        b_c=[1.132856760e5, -1.980421677e3, 1.365384188e1, -4.636096440e-2,
            1.021333011e-4, -1.082893179e-7, 4.472258860e-11; 3.356004410e5, -2.596528368e3,
            6.948841910, -3.484836090e-3, 1.844192445e-6, -5.055205960e-10,
            5.750639010e-14; 2.023105106e5, -1.138235908e3, 3.700279500, -1.833807727e-4,
            6.343683250e-8, -7.068589480e-12, 3.335435980e-16] .* fill({U.K^(3
             - i) for i in 1:7}, size(T_lim_c, 1) - 1),
        B_c=[8.943859760e3*U.K, -7.295824740e1; 1.398412456e4*U.K, -4.477183040e1;
            5.848134850e3, -2.350925275e1] - b_c[:, 2:3]*ln(U.K),
        r=170*U.pico*U.m/U.q);
      annotation (Documentation(info="<html>

     <p>Assumptions:
     <ol>
     <li>Constant specific volume (i.e., incompressible and without
          thermal expansion)</li>
     </ol>
     </p>

     <p>Additional notes:
     <ul>
     <li>The data for this species is for C rather than C<sup>+</sup> (with the exception of specific mass).</li>
     <li>The radius is from <a href=\"http://en.wikipedia.org/wiki/Carbon\">http://en.wikipedia.org/wiki/Carbon</a>.  See also
   <a href=\"http://en.wikipedia.org/wiki/Van_der_Waals_radius\">http://en.wikipedia.org/wiki/Van_der_Waals_radius</a>.</li>
     <li>The default specific volume (<code>v=U.cc*m/(2.210*U.g)</code>) is of pyrolytic graphite
  at 300 K according to [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, p. 909].  Other forms
  are (also at 300 K and based on the same reference) are:
  <ul>
       <li>Amorphous carbon:  <code>v=U.cc*m/(1.950*U.g)</code></li>
       <li>Diamond (type IIa):  <code>v=U.cc*m/(3.500*U.g)</code></li>
       </li>
     </ul>
     </p>
<p>For more information, see the
  <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic\">Characteristic</a> record.</p></html>"));

    end Graphite;

  end 'C+';

  package 'C19HF37O5S-'
    "<html>C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S<sup>-</sup></html>"
    extends Modelica.Icons.Package;
    package Ionomer "C9HF17O5S- ionomer"
      // Note:  HTML formatting isn't used in the description because
      // Dymola 7.4 doesn't support it in the GUI for replaceable lists.
      // The same applies to other species.

      extends BaseClasses.Characteristic(
        final formula="C9HF17O5S-",
        final phase="solid",
        p0=U.atm,
        m=1044.214*U.g/U.mol - 'H+'.Ionomer.m,
        specVolPow={0,0},
        b_v=[U.cc*m/(2.00*U.g)],
        p_min=-Modelica.Constants.inf,
        Deltah0_f=0,
        Deltah0=0,
        specHeatCapPow=0,
        T_lim_c={0,Modelica.Constants.inf},
        b_c=[4188*U.J*m/(U.kg*U.K)],
        B_c=[-298.15*U.K*b_c[1, 1] + Deltah0_f, 0],
        r=(147 + 2259.8/2)*U.pico*U.m/U.q);
      annotation (Documentation(info="<html>
       <p>Assumptions:
     <ol>
     <li>Constant specific volume (i.e., incompressible and without
          thermal expansion)</li>
     </ol>
     </p>

     <p>Additional notes:
     <ul>
     <li>Most of the data for this species is for C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S rather than
     C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S<sup>-</sup> (with the exception of specific mass).</li>
     <li>A form of C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S is
     C<sub>7</sub>HF<sub>13</sub>O<sub>5</sub>S.(C<sub>2</sub>F<sub>4</sub>)<sub>6</sub>, which is a typical
   configuration of Nafion sulfonate resin after hydrolysis
   [<a href=\"modelica://FCSys.UsersGuide.References\">Mark1999</a>, p. 234].</li>
     <li>Thermodynamic data for this material is not available from
     [<a href=\"modelica://FCSys.UsersGuide.References\">McBride2002</a>].  The default specific heat capacity
     (<code>b_c=[4188*U.J*m/(U.kg*U.K)]</code>) based on [<a href=\"modelica://FCSys.UsersGuide.References\">Shah2009</a>, p. B472].</li>
     <li>According to [<a href=\"modelica://FCSys.UsersGuide.References\">Avogadro1.03</a>], the furthest distance
   between two atoms of C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S is 2259.8 pm and is between fluorines.
   The radius of F is 147 pm (<a href=\"http://en.wikipedia.org/wiki/Fluorine\">http://en.wikipedia.org/wiki/Fluorine</a>).</li>
     <li>From [<a href=\"modelica://FCSys.UsersGuide.References\">Avogadro1.03</a>], the molecular weight of C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S
   is 1044.214 g/mol.  However, the \"11\" in Nafion 11x indicates a
   molecular weight of 1100 g/mol.  According to
   <a href=\"http://en.wikipedia.org/wiki/Nafion\">http://en.wikipedia.org/wiki/Nafion</a>,
       \"the molecular weight of Nafion is uncertain due to differences in
        processing and solution morphology.\"</li>
     <li>The specific volume (<code>v = U.cc*m/(2.00*U.g)</code>) is based on
   [<a href=\"modelica://FCSys.UsersGuide.References\">Lin2006</a>, p. A1327].
       </li>
     </ul>
     </p>

<p>For more information, see the
  <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic\">Characteristic</a> record.</p></html>"));

    end Ionomer;

  end 'C19HF37O5S-';

  package 'e-' "<html>e<sup>-</sup></html>"
    extends Modelica.Icons.Package;
    package Gas "e- gas"
      import Data = Modelica.Media.IdealGases.Common.SingleGasesData.eminus;

      extends BaseClasses.Characteristic(
        final formula="e-",
        phase="gas",
        m=Data.MM*U.kg/U.mol,
        b_v=[1],
        specVolPow={-1,0},
        Deltah0_f=Data.MM*Data.Hf*U.J/U.mol,
        Deltah0=Data.MM*Data.H0*U.J/U.mol,
        T_lim_c={200.000,20000.000}*U.K,
        b_c={Data.alow} .* fill({U.K^(3 - i) for i in 1:size(Data.alow, 1)},
            size(T_lim_c, 1) - 1),
        B_c={Data.blow} .* fill({U.K,1}, size(T_lim_c, 1) - 1) - b_c[:, 2:3]*ln(
            U.K),
        r=U.k_A/m);
      annotation (Documentation(info="<html>
     <p>Notes:
     <ul>
     <li>The equation for the radius is the classical radius of an electron (see
  <a href=\"http://en.wikipedia.org/wiki/Classical_electron_radius\">http://en.wikipedia.org/wiki/Classical_electron_radius</a>).</li>
  <li>McBride and Gordon [<a href=\"modelica://FCSys.UsersGuide.References\">McBride1996</a>] provide correlations for the transport
  properties of e<sup>-</sup> gas.  However, they are not entered here, since they
  contain only one temperature range (2000 to 5000 K) which is beyond the expected operating range of the model.</li>
  <li>
  The thermodynamic data [<a href=\"modelica://FCSys.UsersGuide.References\">McBride2002</a>] splits the correlations into three
  temperature ranges, but the coefficients are the same for each.
   Therefore, the temperature limits are set here such that the entire
   range is handled without switching.  The lower temperature limit
   in the source data is 298.150 K, but here it is expanded down to 200 K.
   The constants are independent of temperature anyway.
  </li>
</ul>
</p>
<p>For more information, see the
  <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic\">Characteristic</a> record.</p></html>"));

    end Gas;

    package Graphite "e- in graphite"
      extends Gas(
        final phase="graphite",
        specVolPow='C+'.Graphite.specVolPow,
        b_v='C+'.Graphite.b_v);
      annotation (Documentation(info="<html>
  <p>Assumptions:<ol>
    <li>The density of e<sup>-</sup> is equal to that of C<sup>+</sup> as graphite.
    </li>
    </ol></p>

     <p>For more information, see the
     <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic\">Characteristic</a> record.</p></html>"));

    end Graphite;

  end 'e-';

  package 'H+' "<html>H<sup>+</sup></html>"
    extends Modelica.Icons.Package;
    package Gas "H+ gas"
      import Data = Modelica.Media.IdealGases.Common.SingleGasesData.Hplus;

      extends BaseClasses.Characteristic(
        final formula="H+",
        phase="gas",
        final m=Data.MM*U.kg/U.mol,
        b_v=[1],
        specVolPow={-1,0},
        Deltah0_f=Data.MM*Data.Hf*U.J/U.mol,
        Deltah0=Data.MM*Data.H0*U.J/U.mol,
        specHeatCapPow=-2,
        T_lim_c={200.000,20000.000}*U.K,
        b_c={Data.alow} .* fill({U.K^(3 - i) for i in 1:size(Data.alow, 1)},
            size(T_lim_c, 1) - 1),
        B_c={Data.blow} .* fill({U.K,1}, size(T_lim_c, 1) - 1) - b_c[:, 2:3]*ln(
            U.K),
        r=120*U.pico*U.m/U.q);
      annotation (Documentation(info="<html>
         <p>Assumptions:
     <ol>
  <li>The radius (for the purpose of the rigid-sphere assumption of
  kinetic theory) is that of H
  (<a href=\"http://en.wikipedia.org/wiki/Hydrogen\">http://en.wikipedia.org/wiki/Hydrogen</a>).</li>
     </ol>
     </p>

     <p>Additional notes:
     <ul>
     <li>The source data [<a href=\"modelica://FCSys.UsersGuide.References\">McBride2002</a>] breaks the data into three
   temperature ranges, but the constants are the same for each.
   Therefore, the temperature limits are set here such that the entire
   range is handled without switching.
   The lower temperature limit in the source data
   is 298.150 K, but here it is expanded down to 200 K.  The constants are
   independent of temperature anyway.</li>
   <li>The enthalpy to produce H<sup>+</sup> from H<sub>2</sub>O (H<sub>2</sub>O &#8640; OH<sup>-</sup> + H<sup>+</sup>) is
   -96569.804 J/mol.  The enthalpy to produce H<sup>+</sup> from H<sub>3</sub>O<sup>+</sup>
   (H<sub>3</sub>O<sup>+</sup> &#8640; H<sub>2</sub>O + H<sup>+</sup>) is 839826.0 J/mol.  Based on
   [<a href=\"modelica://FCSys.UsersGuide.References\">Tissandier1998</a>], the
   enthalpy of formation of aqueous H<sup>+</sup> is 1150.1e3 J/mol.</li>
     </ul>
     </p>
  <p>For more information, see the
  <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic\">Characteristic</a> record.</p></html>"));

    end Gas;

    package Ionomer "H+ in ionomer"

      extends Gas(
        final phase="solid",
        specVolPow={0,0},
        b_v=[1/(0.95*U.M)]);
      annotation (Documentation(info="<html>
  <p>The specific volume of protons corresponds to the concentration measured
  by Spry and Fayer (0.95 M) in Nafion<sup>&reg;</sup> at
  &lambda; = 12, where &lambda; is the number of
  H<sub>2</sub>O molecules to SO<sub>3</sub>H
  endgroups.  At &lambda; = 22, the concentration was measured at 0.54 M
  [<a href=\"modelica://FCSys.UsersGuide.References\">Spry2009</a>].</p>
</html>"));

    end Ionomer;

  end 'H+';

  package H2 "<html>H<sub>2</sub></html>"
    extends Modelica.Icons.Package;
    package Gas "H2 gas"
      import Data = Modelica.Media.IdealGases.Common.SingleGasesData.H2;

      extends BaseClasses.CharacteristicNASA(
        final formula=Data.name,
        final phase="gas",
        final m=Data.MM*U.kg/U.mol,
        specVolPow={-1,-3},
        b_v={{0,0,0,1},{8.0282e6*U.K^3,-2.6988e5*U.K^2,-129.26*U.K,17.472}*U.cm
            ^3/U.mol},
        p_min=1e-323*p0,
        Deltah0_f=Data.MM*Data.Hf*U.J/U.mol,
        Deltah0=Data.MM*Data.H0*U.J/U.mol,
        T_lim_c={200.000,Data.Tlimit,6000.000,20000.000}*U.K,
        b_c={Data.alow,Data.ahigh,{4.966884120e8,-3.147547149e5,79.84121880,-8.414789210e-3,
            4.753248350e-7,-1.371873492e-11,1.605461756e-16}} .* fill({U.K^(3
             - i) for i in 1:size(Data.alow, 1)}, size(T_lim_c, 1) - 1),
        B_c={Data.blow,Data.bhigh,{2.488433516e6,-669.5728110}} .* fill({U.K,1},
            size(T_lim_c, 1) - 1) - b_c[:, 2:3]*ln(U.K),
        r=(120 + 100.3/2)*U.pico*U.m/U.q,
        T_lim_F_R={200.0,1000.0,5000.0,15000.0}*U.K,
        b_F={fromNASAViscosity({0.74553182,43.555109,-3.2579340e3,0.13556243}),
            fromNASAViscosity({0.96730605,679.31897,-2.1025179e5,-1.8251697}),
            fromNASAViscosity({1.0126129,1.4973739e3,-1.4428484e6,-2.3254928})},

        b_R={fromNASAThermalConductivity({1.0059461,279.51262,-2.9792018e4,
            1.1996252}),fromNASAThermalConductivity({1.0582450,248.75372,
            1.1736907e4,0.82758695}),fromNASAThermalConductivity({-0.22364420,-6.9650442e3,
            -7.7771313e4,13.189369})});

      // Note:  In Dymola 7.4 ln(1e-323) returns a valid result, but ln(1e-324)
      // doesn't.
      annotation (Documentation(info="<html>
            <p>Notes:
     <ul>
  <li>According to [<a href=\"modelica://FCSys.UsersGuide.References\">Avogadro1.03</a>], the (center-to-center)
   bond length of H-H is 100.3 pm.  The radius of H is from
   <a href=\"http://en.wikipedia.org/wiki/Hydrogen\">http://en.wikipedia.org/wiki/Hydrogen</a>.  See also
   <a href=\"http://en.wikipedia.org/wiki/Van_der_Waals_radius\">http://en.wikipedia.org/wiki/Van_der_Waals_radius</a>.</li>
  <li>The virial coefficients are from [<a href=\"modelica://FCSys.UsersGuide.References\">Dymond2002</a>, p. 41].  The
  temperature range of the coefficients is [60, 500] K, but this is not enforced in the functions.</li>
     </ul>
     </p>

<p>For more information, see the
  <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic\">Characteristic</a> record.</p></html>"));

    end Gas;

  end H2;

  package H2O "<html>H<sub>2</sub>O</html>"
    extends Modelica.Icons.Package;
    package Gas "H2O gas"
      import Data = Modelica.Media.IdealGases.Common.SingleGasesData.H2O;

      extends BaseClasses.CharacteristicNASA(
        final formula=Data.name,
        final phase="gas",
        final m=Data.MM*U.kg/U.mol,
        specVolPow={-1,-3},
        b_v={{0,0,0,1},{-5.6932e10*U.K^3,1.8189e8*U.K^2,-3.0107e5*U.K,158.83}*U.cm
            ^3/U.mol},
        p_min=1e-323*p0,
        Deltah0_f=Data.MM*Data.Hf*U.J/U.mol,
        Deltah0=Data.MM*Data.H0*U.J/U.mol,
        T_lim_c={200.000,Data.Tlimit,6000.000}*U.K,
        b_c={Data.alow,Data.ahigh} .* fill({U.K^(3 - i) for i in 1:size(Data.alow,
            1)}, size(T_lim_c, 1) - 1),
        B_c={Data.blow,Data.bhigh} .* fill({U.K,1}, size(T_lim_c, 1) - 1) - b_c[
            :, 2:3]*ln(U.K),
        r=(282/2)*U.pico*U.m/U.q,
        T_lim_F_R={373.2,1073.2,5000.0,15000.0}*U.K,
        b_F={fromNASAViscosity({0.50019557,-697.12796,8.8163892e4,3.0836508}),
            fromNASAViscosity({0.58988538,-537.69814,5.4263513e4,2.3386375}),
            fromNASAViscosity({0.64330087,-95.668913,-3.7742283e5,1.8125190})},

        b_R={fromNASAThermalConductivity({1.0966389,-555.13429,1.0623408e5,-0.24664550}),
            fromNASAThermalConductivity({0.39367933,-2.2524226e3,6.1217458e5,
            5.8011317}),fromNASAThermalConductivity({-0.41858737,-1.4096649e4,
            1.9179190e7,14.345613})});

      // Note:  In Dymola 7.4 ln(1e-323) returns a valid result, but ln(1e-324)
      // doesn't.
      annotation (Documentation(info="<html>
        <p>Notes:
     <ul>
  <li>The radius of H<sub>2</sub>O is 282 pm
   (<a href=\"http://www.lsbu.ac.uk/water/molecule.html\">http://www.lsbu.ac.uk/water/molecule.html</a>).  Using the radius of H
   from <a href=\"http://en.wikipedia.org/wiki/Hydrogen\">http://en.wikipedia.org/wiki/Hydrogen</a> and the center-to-center
   distance of hydrogen atoms in H<sub>2</sub>O from [<a href=\"modelica://FCSys.UsersGuide.References\">Avogadro1.03</a>],
   156.6 pm, the radius of H<sub>2</sub>O would be (120 + 156.6/2) pm = 198.3 pm.</li>
  <li>The virial coefficients are from [<a href=\"modelica://FCSys.UsersGuide.References\">Dymond2002</a>, p. 4]).  The
  temperature range of the coefficients is [350, 770] K, but this is not enforced in the functions.</li>
     </ul>
     </p>

  <p>For more information, see the
  <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic\">Characteristic</a> record.</p></html>"));

    end Gas;

    package Ionomer "H2O in ionomer"
      extends Gas;
      annotation (Documentation(info="<html>
        <p>Assumptions:
     <ol>
  <li>The properties are the same as H<sub>2</sub>O gas.</li>
     </ol>
     </p>

  <p>For more information, see the
  <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic\">Characteristic</a> record.</p></html>"));

    end Ionomer;

    package Liquid "H2O liquid"

      extends BaseClasses.Characteristic(
        final formula="H2O",
        final phase="liquid",
        final m=0.01801528*U.kg/U.mol,
        p0=U.atm,
        specVolPow={0,0},
        b_v=[U.cc*m/(0.99656*U.g)],
        Deltah0_f=-285830.000*U.J/U.mol,
        Deltah0=13278.000*U.J/U.mol,
        T_lim_c={273.150,373.150,600.000}*U.K,
        b_c=[1.326371304e9, -2.448295388e7, 1.879428776e5, -7.678995050e2,
            1.761556813, -2.151167128e-3, 1.092570813e-6; 1.263631001e9, -1.680380249e7,
            9.278234790e4, -2.722373950e2, 4.479243760e-1, -3.919397430e-4,
            1.425743266e-7] .* fill({U.K^(3 - i) for i in 1:7}, size(T_lim_c, 1)
             - 1),
        B_c=[1.101760476e8*U.K, -9.779700970e5; 8.113176880e7*U.K, -5.134418080e5]
             - b_c[:, 2:3]*ln(U.K),
        r=(282/2)*U.pico*U.m/U.q);
      annotation (Documentation(info="<html>     <p>Assumptions:
     <ol>
     <li>Constant specific volume (i.e., incompressible and without
          thermal expansion)</li>
     </ol>
     </p>
            <p>Additional notes:
     <ul>
     <li>See note in <a href=\"modelica://FCSys.Characteristics.H2O.Gas\">Characteristics.H2O.Gas</a> regarding the radius.</li>
     <li>The default specific volume (<code>b_v=[U.cc*m/(0.99656*U.g)]</code>) is at 300 K based on [<a href=\"modelica://FCSys.UsersGuide.References\">Takenaka1990</a>].</li>
     </ul>
     </p>
  <p>For more information, see the
  <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic\">Characteristic</a> record.</p></html>"));

    end Liquid;

  end H2O;

  package N2 "<html>N<sub>2</sub></html>"
    extends Modelica.Icons.Package;
    package Gas "N2 gas"
      import Data = Modelica.Media.IdealGases.Common.SingleGasesData.N2;

      extends BaseClasses.CharacteristicNASA(
        final formula=Data.name,
        final phase="gas",
        final m=Data.MM*U.kg/U.mol,
        specVolPow={-1,-4},
        b_v={{0,0,0,0,1},{-2.7198e9*U.K^4,6.1253e7*U.K^3,-1.4164e6*U.K^2,-9.3378e3
            *U.K,40.286}*U.cc/U.mol},
        p_min=1e-323*p0,
        Deltah0_f=Data.MM*Data.Hf*U.J/U.mol,
        Deltah0=Data.MM*Data.H0*U.J/U.mol,
        T_lim_c={200.000,Data.Tlimit,6000.000,20000.000}*U.K,
        b_c={Data.alow,Data.ahigh,{8.310139160e8,-6.420733540e5,2.020264635e2,-3.065092046e-2,
            2.486903333e-6,-9.705954110e-11,1.437538881e-15}} .* fill({U.K^(3
             - i) for i in 1:size(Data.alow, 1)}, size(T_lim_c, 1) - 1),
        B_c={Data.blow,Data.bhigh,{4.938707040e6,-1.672099740e3}} .* fill({U.K,
            1}, size(T_lim_c, 1) - 1) - b_c[:, 2:3]*ln(U.K),
        r=(155 + 145.2/2)*U.pico*U.m/U.q,
        T_lim_F_R={200.0,1000.0,5000.0,15000.0}*U.K,
        b_F={fromNASAViscosity({0.62526577,-31.779652,-1.6407983e3,1.7454992}),
            fromNASAViscosity({0.87395209,561.52222,-1.7394809e5,-0.39335958}),
            fromNASAViscosity({0.88503551,909.02171,-7.3129061e5,-0.53503838})},

        b_R={fromNASAThermalConductivity({0.85439436,105.73224,-1.2347848e4,
            0.47793128}),fromNASAThermalConductivity({0.88407146,133.57293,-1.1429640e4,
            0.24417019}),fromNASAThermalConductivity({2.4176185,8.0477749e3,
            3.1055802e6,-14.517761})});

      // Note:  In Dymola 7.4 ln(1e-323) returns a valid result, but ln(1e-324)
      // doesn't.
      annotation (Documentation(info="<html>
                  <p>Notes:
     <ul>
  <li>According to [<a href=\"modelica://FCSys.UsersGuide.References\">Avogadro1.03</a>], the (center-to-center)
   bond length of N-N is 145.2 pm.  The radius of N is from
   <a href=\"http://en.wikipedia.org/wiki/Nitrogen\">http://en.wikipedia.org/wiki/Nitrogen</a>.  See also
   <a href=\"http://en.wikipedia.org/wiki/Van_der_Waals_radius\">http://en.wikipedia.org/wiki/Van_der_Waals_radius</a>.</li>
  <li>The virial coefficients are from [<a href=\"modelica://FCSys.UsersGuide.References\">Dymond2002</a>, p. 69].  The
  temperature range of the coefficients is [75, 745] K, but this is not enforced in the functions.  More precise virial coefficients are available from
  <a href=\"http://www.tpub.com/content/nasa1996/NASA-96-cr4755/NASA-96-cr47550059.htm\">http://www.tpub.com/content/nasa1996/NASA-96-cr4755/NASA-96-cr47550059.htm</a>.</li>
     </ul>
     </p>
  <p>For more information, see the
  <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic\">Characteristic</a> record.</p></html>"));

    end Gas;

  end N2;

  package O2 "<html>O<sub>2</sub></html>"
    extends Modelica.Icons.Package;
    package Gas "O2 gas"
      import Data = Modelica.Media.IdealGases.Common.SingleGasesData.O2;

      extends BaseClasses.CharacteristicNASA(
        final formula=Data.name,
        final phase="gas",
        final m=Data.MM*U.kg/U.mol,
        specVolPow={-1,-4},
        b_v={{0,0,0,0,1},{5.0855e9*U.K^4,-1.6393e8*U.K^3,5.2007e5*U.K^2,-1.7696e4
            *U.K,42.859}*U.cc/U.mol},
        p_min=1e-323*p0,
        Deltah0_f=Data.MM*Data.Hf*U.J/U.mol,
        Deltah0=Data.MM*Data.H0*U.J/U.mol,
        T_lim_c={200.000,Data.Tlimit,6000.000,20000.000}*U.K,
        b_c={Data.alow,Data.ahigh,{4.975294300e8,-2.866106874e5,6.690352250e1,-6.169959020e-3,
            3.016396027e-7,-7.421416600e-12,7.278175770e-17}} .* fill({U.K^(3
             - i) for i in 1:size(Data.alow, 1)}, size(T_lim_c, 1) - 1),
        B_c={Data.blow,Data.bhigh,{2.293554027e6,-5.530621610e2}} .* fill({U.K,
            1}, size(T_lim_c, 1) - 1) - b_c[:, 2:3]*ln(U.K),
        r=(152 + 128.2/2)*U.pico*U.m/U.q,
        T_lim_F_R={200.0,1000.0,5000.0,15000.0}*U.K,
        b_F={fromNASAViscosity({0.60916180,-52.244847,-599.74009,2.0410801}),
            fromNASAViscosity({0.72216486,175.50839,-5.7974816e4,1.0901044}),
            fromNASAViscosity({0.73981127,391.94906,-3.7833168e5,0.90931780})},

        b_R={fromNASAThermalConductivity({0.77229167,6.8463210,-5.8933377e3,
            1.2210365}),fromNASAThermalConductivity({0.90917351,291.24182,-7.9650171e4,
            0.064851631}),fromNASAThermalConductivity({-1.1218262,-1.9286378e4,
            2.3295011e7,20.342043})});

      // Note:  In Dymola 7.4 ln(1e-323) returns a valid result, but ln(1e-324)
      // doesn't.
      annotation (Documentation(info="<html><p>Notes:<ul>
  <li>According to [<a href=\"modelica://FCSys.UsersGuide.References\">Avogadro1.03</a>], the (center-to-center)
   bond length of O-O is 128.2 pm.  The radius of O is from
   <a href=\"http://en.wikipedia.org/wiki/Oxygen\">http://en.wikipedia.org/wiki/Oxygen</a>.  See also
   <a href=\"http://en.wikipedia.org/wiki/Van_der_Waals_radius\">http://en.wikipedia.org/wiki/Van_der_Waals_radius</a>.</li>
  <li>The virial coefficients are from [<a href=\"modelica://FCSys.UsersGuide.References\">Dymond2002</a>, p. 69].  The
  temperature range of the coefficients is [70, 495] K, but this is not enforced in the functions.</li>
     </ul>
     </p>
<p>For more information, see the
  <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic\">Characteristic</a> record.</p></html>"));

    end Gas;

  end O2;

  package BaseClasses "Base classes (not generally for direct use)"
    extends Modelica.Icons.BasesPackage;
    package CharacteristicNASA
      "Thermodynamic record with transport properties based on NASA CEA"

      extends Characteristic;

      constant Q.TemperatureAbsolute T_lim_F_R[:]={0,Modelica.Constants.inf}
        "<html>Temperature limits for the rows of <i>b</i><sub><i>F</i></sub> and <i>b</i><sub><i>R</i></sub> (<i>T</i><sub>lim <i>F</i> <i>R</i></sub>)</html>";
      constant Real b_F[size(T_lim_F_R, 1) - 1, 4]
        "<html>Correlation constants for fluidity (<i>b</i><sub><i>F</i></sub>)</html>";
      constant Real b_R[size(T_lim_F_R, 1) - 1, 4]
        "<html>Correlation constants for thermal resistivity (<i>b</i><sub><i>R</i></sub>)</html>";

    protected
      function fromNASAViscosity
        "Return constants for fluidity given NASA CEA constants for viscosity"
        input Real b_eta[4] "NASA CEA constants for viscosity";
        output Real b_F[4] "Constants for fluidity";

      algorithm
        b_F := {-b_eta[1],-b_eta[2]*U.K,-b_eta[3]*U.K^2,-b_eta[4] + b_eta[1]*ln(
          U.K) + ln(1e4*U.m*U.s/U.g)} annotation (Inline=true);

      end fromNASAViscosity;

      function fromNASAThermalConductivity
        "Return constants for thermal resistivity given NASA CEA constants for thermal conductivity"
        input Real b_lambda[4] "NASA CEA constants for thermal conductivity";
        output Real b_R[4] "Constants for thermal resistivity";

      algorithm
        b_R := {-b_lambda[1],-b_lambda[2]*U.K,-b_lambda[3]*U.K^2,-b_lambda[4]
           + b_lambda[1]*ln(U.K) + ln(1e4*U.m*U.K/U.W)}
          annotation (Inline=true);

      end fromNASAThermalConductivity;

    public
      redeclare function F "Fluidity as a function of temperature"
        extends Modelica.Icons.Function;

        input Q.TemperatureAbsolute T=298.15*U.K "Temperature";
        input Q.PressureAbsolute p=U.atm "Pressure";
        // Note:  Pressure isn't used here but is included for generality.
        output Q.FluidityDynamic F "Dynamic fluidity";

      algorithm
        /*
    assert(T_lim_F_R[1] <= T and T <= T_lim_F_R[size(T_lim_F_R, 1)], "Temperature "
     + String(T/(U.K)) + " K is out of range for the resistivities ([" +
    String(T_lim_F_R[1]/U.K) + ", " + String(T_lim_F_R[size(T_lim_F_R, 1)]
    /U.K) + "] K).");
    */
        // Note:  This is commented out so that the function can be inlined.
        // Note:  In Dymola 7.4 T_lim_F_R[end] can't be used instead of
        // T_lim_F_R[size(T_lim_F_R, 1)] due to:
        //     "Error, not all "end" could be expanded."

        F := smooth(0, exp(sum(if (T_lim_F_R[i] <= T or i == 1) and (T <
          T_lim_F_R[i + 1] or i == size(T_lim_F_R, 1) - 1) then b_F[i, 1]*ln(T)
           + (b_F[i, 2] + b_F[i, 3]/T)/T + b_F[i, 4] else 0 for i in 1:size(
          T_lim_F_R, 1) - 1)))
          annotation (
          InlineNoEvent=true,
          Inline=true,
          smoothOrder=0);
        annotation (Documentation(info="<html><p>This function is based on based on NASA CEA
  [<a href=\"modelica://FCSys.UsersGuide.References\">McBride1996</a>, <a href=\"modelica://FCSys.UsersGuide.References\">Svehla1995</a>]</p>

  <p>For more information, see <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.F\">Characteristic.F</a>().</p>
  </html>"));
      end F;

      redeclare function R
        "Thermal resistivity as a function of temperature and pressure"
        extends Modelica.Icons.Function;

        input Q.TemperatureAbsolute T=298.15*U.K "Temperature";
        input Q.PressureAbsolute p=U.atm "Pressure";
        // Note:  Pressure isn't used here but is included for generality.
        output Q.ResistivityThermal R "Thermal resistivity";

      algorithm
        /*
    assert(T_lim_F_R[1] <= T and T <= T_lim_F_R[size(T_lim_F_R, 1)], "Temperature "
     + String(T/(U.K)) + " K is out of range for the resistivities ([" +
    String(T_lim_F_R[1]/U.K) + ", " + String(T_lim_F_R[size(T_lim_F_R, 1)]
    /U.K) + "] K).");
    */
        // Note:  This is commented out so that the function can be inlined.
        // Note:  In Dymola 7.4 T_lim_F_R[end] can't be used instead of
        // T_lim_F_R[size(T_lim_F_R, 1)] due to:
        //     "Error, not all "end" could be expanded."

        R := smooth(0, exp(sum(if (T_lim_F_R[i] <= T or i == 1) and (T <
          T_lim_F_R[i + 1] or i == size(T_lim_F_R, 1) - 1) then b_R[i, 1]*ln(T)
           + (b_R[i, 2] + b_R[i, 3]/T)/T + b_R[i, 4] else 0 for i in 1:size(
          T_lim_F_R, 1) - 1)))
          annotation (
          InlineNoEvent=true,
          Inline=true,
          smoothOrder=0);
        annotation (Documentation(info="<html><p>This function is based on based on NASA CEA
  [<a href=\"modelica://FCSys.UsersGuide.References\">McBride1996</a>, <a href=\"modelica://FCSys.UsersGuide.References\">Svehla1995</a>]</p>

    <p>For more information, see <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.R\">Characteristic.R</a>().</p>
    </html>"));
      end R;
      annotation (defaultComponentPrefixes="replaceable",Documentation(info="<html><p>The correlations for transport properties are available in
  [<a href=\"modelica://FCSys.UsersGuide.References\">McBride1996</a>,
  <a href=\"modelica://FCSys.UsersGuide.References\">McBride2002</a>]. For more information, see the
  <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic\">Characteristic</a> record.</p></html>"));

    end CharacteristicNASA;

    package Characteristic "Record for thermodynamic and resistive properties"
      import Modelica.Math.BooleanVectors.anyTrue;

      extends Modelica.Icons.MaterialPropertiesPackage;

      constant String formula "Chemical formula";
      constant String phase "Material phase";
      constant Q.MassSpecific m(min=Modelica.Constants.small) "Specific mass";
      // Note:  The positive minimum value prevents a structural singularity
      // when checking FCSys.Subregions.Species.SpeciesInertStagnant in Dymola
      // 7.4.
      constant Q.LengthSpecific r "Specific radius" annotation (Dialog);
      final constant Integer z=Chemistry.charge(formula) "Charge number";
      constant Q.PressureAbsolute p0=U.bar
        "<html>Reference pressure (<i>p</i>&deg;)</html>";
      constant Q.Pressure p_min=-Modelica.Constants.inf
        "<html>Minimum pressure for numerical protection (<i>p</i><sub>min</sub>)</html>";
      constant Real b_v[:, :]=[1]
        "<html>Coefficients of specific volume as a polynomial in <i>p</i>/<i>T</i> and <i>T</i> (<i>b</i><sub><i>v</i></sub>)</html>";
      // Note:  p/T is the argument instead of p so that b_p will have the
      // same size as b_v for the typical definitions of the second virial
      // coefficients in [Dymond2002].
      constant Integer specVolPow[2]={-1,0}
        "<html>Powers of <i>p</i>/<i>T</i> and <i>T</i> for 1<sup>st</sup> row and column of <i>b</i><sub><i>v</i></sub>, respectively</html>";
      constant Q.PotentialChemical Deltah0_f
        "<html>Enthalpy of formation at 298.15 K, <i>p</i>&deg; (&Delta;<i>h</i>&deg;<sub>f</sub>)</html>";
      constant Q.PotentialChemical Deltah0
        "<html><i>h</i>&deg;(298.15 K) - <i>h</i>&deg;(0 K) (&Delta;<i>h</i>&deg;)</html>";
      constant Q.PotentialChemical h_offset=0
        "<html>Additional enthalpy offset (<i>h</i><sub>offset</sub>)</html>";
      constant Integer specHeatCapPow=-2
        "<html>Power of <i>T</i> for 1<sup>st</sup> column of <i>b</i><sub><i>c</i></sub></html>";
      constant Q.TemperatureAbsolute T_lim_c[:]={0,Modelica.Constants.inf}
        "<html>Temperature limits for the rows of <i>b</i><sub><i>c</i></sub> and <i>B</i><sub><i>c</i></sub> (<i>T</i><sub>lim <i>c</i></sub>)</html>";
      constant Real b_c[size(T_lim_c, 1) - 1, :]
        "<html>Coefficients of isobaric specific heat capacity at <i>p</i>&deg; as a polynomial in <i>T</i> (<i>b</i><sub><i>c</i></sub>)</html>";
      constant Real B_c[size(T_lim_c, 1) - 1, 2]
        "<html>Integration constants for specific enthalpy and entropy (<i>B</i><sub><i>c</i></sub>)</html>";

      final constant Boolean isCompressible=anyTrue({anyTrue({abs(b_v[i, j]) >
          Modelica.Constants.small and specVolPow[1] + i - 1 <> 0 for i in 1:
          size(b_v, 1)}) for j in 1:size(b_v, 2)})
        "true, if specific volume depends on pressure";
      final constant Boolean hasThermalExpansion=anyTrue({anyTrue({abs(b_v[i, j])
           > Modelica.Constants.small and specVolPow[2] + j - specVolPow[1] - i
           <> 0 for i in 1:size(b_v, 1)}) for j in 1:size(b_v, 2)})
        "true, if specific volume depends on temperature";

    protected
      constant Integer pressPow[2]={specVolPow[1] - size(b_v, 1) + 1,specVolPow[
          2] + 1}
        "Powers of v and T for 1st row and column of b_p, respectively";
      final constant Real b_p[size(b_v, 1), size(b_v, 2)]=if size(b_v, 1) == 1
           then b_v .^ (-pressPow[1]) else {(if specVolPow[1] + i == 0 or
          specVolPow[1] + i == 1 or size(b_v, 1) == 1 then b_v[i, :] else (if
          specVolPow[1] + i == 2 and specVolPow[1] <= 0 then b_v[i, :] + b_v[i
           - 1, :] .^ 2 else (if specVolPow[1] + i == 3 and specVolPow[1] <= 0
           then b_v[i, :] + b_v[i - 2, :] .* (b_v[i - 2, :] .^ 2 + 3*b_v[i - 1,
          :]) else zeros(size(b_v, 2))))) for i in size(b_v, 1):-1:1}
        "Coefficients of p as a polynomial in v and T";
      // Note:  This is from [Dymond2002, p. 2].  If necessary, additional terms
      // can be computed using FCSys/resources/virial-relations.cdf.

      partial function alpha
        "<html>Ideal base resistivity factor as a function of temperature (&alpha;)</html>"
        extends Modelica.Icons.Function;

        input Q.TemperatureAbsolute T=298.15*U.K "Temperature";
        output Q.Resistivity alpha "Resistivity";

      algorithm
        alpha := 6*U.pi*r^2*U.q*sqrt(U.pi*m/T) annotation (Inline=true);
        annotation (Documentation(info="<html>
  <p>This function is based on the kinetic theory of gases with the rigid-sphere (\"billiard-ball\")
  assumption [<a href=\"modelica://FCSys.UsersGuide.References\">Present1958</a>].  It is
  independent of pressure or specific volume.  This independence very accurately matches the measured
  dynamic fluidity of gases.  However, the dynamic fluidity varies by species and
  generally falls more rapidly with temperature than shown by this function
  [<a href=\"modelica://FCSys.UsersGuide.References\">Present1958</a>, p. 41].</p></html>"));
      end alpha;

    public
      function beta_T
        "<html>Isothermal compressibility as a function of temperature and pressure (&beta;<sub><i>T</i></sub>)</html>"
        extends Modelica.Icons.Function;

        input Q.TemperatureAbsolute T=298.15*U.K "Temperature";
        input Q.PressureAbsolute p=U.atm "Pressure";
        output Q.PressureReciprocal beta_T "Isothermal compressibility";

      algorithm
        beta_T := -dv_Tp(
                T=T,
                p=p,
                dT=0,
                dp=1)/v_Tp(T, p) annotation (Inline=true);
        annotation (Documentation(info="<html>
  <p>Note that the compressibility given by this function is static&mdash;unique
  from the dynamic compressibility given by
  <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.Xi\">Xi</a>().</p>

  <p>For an ideal gas, this function is independent of temperature
  (although temperature remains as a valid input).</p>
  </html>"));
      end beta_T;

      function c_p
        "<html>Isobaric specific heat capacity as a function of temperature and pressure (<i>c</i><sub><i>p</i></sub>)</html>"
        extends Modelica.Icons.Function;

        input Q.TemperatureAbsolute T=298.15*U.K "Temperature";
        input Q.PressureAbsolute p=U.atm "Pressure";
        output Q.CapacityThermalSpecific c_p "Isobaric specific heat capacity";

      protected
        function c0_p
          "Isobaric specific heat capacity at reference pressure as a function of temperature"

          input Q.TemperatureAbsolute T "Temperature";
          output Q.CapacityThermalSpecific c0_p
            "Isobaric specific heat capacity at reference pressure";

        algorithm
          /*
    assert(T_lim_c[1] <= T and T <= T_lim_c[size(T_lim_c, 1)], "Temperature " +
      String(T/U.K) + " K is out of bounds for " + name + " ([" + String(T_lim_c[1]
      /U.K) + ", " + String(T_lim_c[size(T_lim_c, 1)]/U.K) + "] K).");
    */
          // Note:  This is commented out so that the function can be inlined.
          // Note:  In Dymola 7.4 T_lim_c[size(T_lim_c, 1)] must be used
          // instead of T_lim_c[end] due to:
          //    "Error, not all 'end' could be expanded."
          c0_p := smooth(0, sum(if (T_lim_c[i] <= T or i == 1) and (T < T_lim_c[
            i + 1] or i == size(T_lim_c, 1) - 1) then Polynomial.f(
                    T,
                    b_c[i, :],
                    specHeatCapPow) else 0 for i in 1:size(T_lim_c, 1) - 1))
            annotation (
            InlineNoEvent=true,
            Inline=true,
            smoothOrder=0);

        end c0_p;

        function c_p_resid
          "Residual isobaric specific heat capacity for pressure adjustment"
          input Q.TemperatureAbsolute T "Temperature";
          input Q.PressureAbsolute p "Pressure";
          input Integer rowLimits[2]={1,size(b_v, 1)}
            "Beginning and ending indices of rows of b_v to be included";
          output Q.CapacityThermalSpecific c_p_resid
            "Temperature times the partial derivative of the integral of (dels/delp)_T*dp up to p w.r.t. T";

        algorithm
          c_p_resid := Polynomial.F(
                    p,
                    {Polynomial.f(
                      T,
                      b_v[i, :] .* {(specVolPow[2] - specVolPow[1] + j - i)*(
                specVolPow[1] - specVolPow[2] + i - j + 1) for j in 1:size(b_v,
                2)},  specVolPow[2] - specVolPow[1] - i) for i in rowLimits[1]:
              rowLimits[2]},
                    specVolPow[1]) annotation (Inline=true);
          // See s_resid() in Characteristic.s for the integral of (dels/delp)_T*dp.
          // This is temperature times the isobaric partial derivative of that function
          // with respect to temperature.  It is zero for an ideal gas.
        end c_p_resid;

      algorithm
        c_p := c0_p(T) + c_p_resid(T, p) - (if phase <> "gas" then c_p_resid(T,
          p0) else c_p_resid(
                T,
                p0,
                {1,-specVolPow[1]})) annotation (Inline=true);
        // See the notes in the algorithm of Characteristic.s.
        // Note:  [Dymond2002, p.17, eqs. 1.45 & 1.46] may be incorrect.

        annotation (Documentation(info="<html>
  <p>For an ideal gas, this function is independent of pressure
  (although pressure remains as a valid input).</p>
  </html>"));
      end c_p;

      function c_v
        "<html>Isochoric specific heat capacity as a function of temperature and pressure (<i>c</i><sub><i>v</i></sub>)</html>"
        extends Modelica.Icons.Function;

        input Q.TemperatureAbsolute T=298.15*U.K "Temperature";
        input Q.PressureAbsolute p=U.atm "Pressure";
        output Q.CapacityThermalSpecific c_v "Isochoric specific heat capacity";

      algorithm
        c_v := c_p(T, p) - T*dp_Tv(
                T,
                v_Tp(T, p),
                dT=1,
                dv=0)*dv_Tp(
                T,
                p,
                dT=1,
                dp=0) "[Moran2004, p. 546, eq. 11.66]" annotation (Inline=true);
        // Note 1:  This reduces to c_v = c_p - 1 for an ideal gas (where in
        // FCSys 1 = U.R).
        // Note 2:  [Dymond2002, p.17, eqs. 1.43 & 1.44] may be incorrect.
        annotation (Documentation(info="<html>
  <p>For an ideal gas, this function is independent of pressure
  (although pressure remains as a valid input).</p>
  </html>"));
      end c_v;

      function dp_Tv
        "<html>Derivative of pressure as defined by <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.p_Tv\">p_Tv</a>()</html>"

        extends Modelica.Icons.Function;

        input Q.TemperatureAbsolute T=298.15*U.K "Temperature";
        input Q.VolumeSpecificAbsolute v=298.15*U.K/U.atm "Specific volume";
        input Q.Temperature dT=0 "Derivative of temperature";
        input Q.VolumeSpecific dv=0 "Derivative of specific volume";
        output Q.Pressure dp "Derivative of pressure";

      algorithm
        dp := if isCompressible then Polynomial.f(
                v,
                {Polynomial.f(
                  T,
                  b_p[i, :] .* {(pressPow[1] + i - 1)*T*dv + (pressPow[2] + j
               - 1)*v*dT for j in 1:size(b_p, 2)},
                  pressPow[2] - 1) for i in 1:size(b_p, 1)},
                pressPow[1] - 1) else 0 annotation (
          Inline=true,
          inverse(dv=dv_Tp(
                      T,
                      p_Tv(T, v),
                      dT,
                      dp)),
          smoothOrder=999);

      end dp_Tv;

      function dv_Tp
        "<html>Derivative of specific volume as defined by <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.v_Tp\">v_Tp</a>()</html>"
        extends Modelica.Icons.Function;

        input Q.TemperatureAbsolute T=298.15*U.K "Temperature";
        input Q.PressureAbsolute p=U.atm "Pressure";
        input Q.Temperature dT=0 "Derivative of temperature";
        input Q.Pressure dp=0 "Derivative of pressure";
        output Q.VolumeSpecific dv "Derivative of specific volume";

      algorithm
        dv := Polynomial.f(
                p,
                {Polynomial.f(
                  T,
                  b_v[i, :] .* {(specVolPow[1] + i - 1)*T*dp + (specVolPow[2]
               - specVolPow[1] + j - i)*p*dT for j in 1:size(b_v, 2)},
                  specVolPow[2] - specVolPow[1] - i) for i in 1:size(b_v, 1)},
                specVolPow[1] - 1) annotation (Inline=true, inverse(dp=dp_Tv(
                      T,
                      v_Tp(T, p),
                      dT,
                      dv)));

      end dv_Tp;

      replaceable function F "Fluidity as a function of temperature"

        extends Modelica.Icons.Function;

        input Q.TemperatureAbsolute T=298.15*U.K "Temperature";
        input Q.PressureAbsolute p=U.atm "Pressure";
        // Note:  Pressure isn't used here but is included for generality.
        output Q.FluidityDynamic F "Dynamic fluidity";

      algorithm
        F := alpha(T)/m annotation (Inline=true);
        annotation (Documentation(info="<html>
<p>Note that dynamic fluidity is defined as the reciprocal of viscosity&mdash;typically dynamic viscosity
(see <a href=\"http://en.wikipedia.org/wiki/Viscosity#Fluidity\">http://en.wikipedia.org/wiki/Viscosity#Fluidity</a>).
</p>
</html>"));
      end F;

      function g "Gibbs potential as a function of temperature and pressure"

        extends Modelica.Icons.Function;

        input Q.TemperatureAbsolute T=298.15*U.K "Temperature";
        input Q.PressureAbsolute p=U.atm "Pressure";
        input ReferenceEnthalpy referenceEnthalpy=ReferenceEnthalpy.EnthalpyOfFormationAt25degC
          "Choice of enthalpy reference";
        output Q.Potential g "Gibbs potential";

      algorithm
        g := h( T,
                p,
                referenceEnthalpy) - T*s(T, p) annotation (Inline=true);

      end g;

      function h "Specific enthalpy as a function of temperature and pressure"

        extends Modelica.Icons.Function;

        input Q.TemperatureAbsolute T=298.15*U.K "Temperature";
        input Q.PressureAbsolute p=U.atm "Pressure";
        input ReferenceEnthalpy referenceEnthalpy=ReferenceEnthalpy.EnthalpyOfFormationAt25degC
          "Choice of enthalpy reference";
        output Q.Potential h "Specific enthalpy";

      protected
        function h0_i
          "Return h0 as a function of T using one of the temperature intervals"
          input Q.TemperatureAbsolute T "Temperature";
          input Integer i "Index of the temperature interval";
          output Q.Potential h0
            "Specific enthalpy at given temperature relative to enthalpy of formation at 25 degC, both at reference pressure";

        algorithm
          h0 := Polynomial.F(
                    T,
                    b_c[i, :],
                    specHeatCapPow) + B_c[i, 1]
            annotation (Inline=true, derivative=dh0_i);
          // This is the integral of c0_p*dT up to T at p0.  The lower bound is the
          // enthalpy of formation (of ideal gas, if the material is gaseous) at
          // 25 degC [McBride2002, p. 2].

        end h0_i;

        function dh0_i "Derivative of h0_i"
          // Note:  This function is necessary for Dymola 7.4 to differentiate h().

          input Q.TemperatureAbsolute T "Temperature";
          input Integer i "Index of the temperature interval";
          input Q.Temperature dT "Derivative of temperature";
          output Q.Potential dh0
            "Derivative of specific enthalpy at reference pressure";

        algorithm
          dh0 := Polynomial.f(
                    T,
                    b_c[i, :],
                    specHeatCapPow)*dT annotation (Inline=true);

        end dh0_i;

        function h_resid
          "Residual specific enthalpy for pressure adjustment for selected rows of b_v"
          input Q.TemperatureAbsolute T "Temperature";
          input Q.PressureAbsolute p "Pressure";
          input Integer rowLimits[2]={1,size(b_v, 1)}
            "Beginning and ending indices of rows of b_v to be included";
          output Q.Potential h_resid
            "Integral of (delh/delp)_T*dp up to p with zero integration constant (for selected rows)";

        algorithm
          h_resid := Polynomial.F(
                    p,
                    {Polynomial.f(
                      T,
                      b_v[i, :] .* {specVolPow[1] - specVolPow[2] + i - j + 1
                for j in 1:size(b_v, 2)},
                      specVolPow[2] - specVolPow[1] - i + 1) for i in rowLimits[
              1]:rowLimits[2]},
                    specVolPow[1] + rowLimits[1] - 1) annotation (Inline=true);
          // Note:  The partial derivative (delh/delp)_T is equal to v +
          // T*(dels/delp)_T by definition of enthalpy change (dh = T*ds + v*dp)
          // and then to v - T*(delv/delT)_p by applying the appropriate Maxwell
          // relation, (dels/delp)_T = -(delv/delT)_p.
          // Note:  This is zero for an ideal gas.

        end h_resid;

      algorithm
        /*
    assert(T_lim_c[1] <= T and T <= T_lim_c[size(T_lim_c, 1)], "Temperature " +
    String(T/U.K) + " K is out of bounds for " + name + " ([" + String(T_lim_c[1]
    /U.K) + ", " + String(T_lim_c[size(T_lim_c, 1)]/U.K) + "] K).");
  */
        // Note:  This is commented out so that the function can be inlined.
        // Note:  In Dymola 7.4 T_lim_c[size(T_lim_c, 1)] must be used
        // instead of T_lim_c[end] due to:
        //    "Error, not all 'end' could be expanded."

        h := smooth(1, sum(if (T_lim_c[i] <= T or i == 1) and (T < T_lim_c[i +
          1] or i == size(b_c, 1)) then h0_i(T, i) else 0 for i in 1:size(b_c,
          1))) + (if referenceEnthalpy == ReferenceEnthalpy.ZeroAt0K then
          Deltah0 else 0) - (if referenceEnthalpy <> ReferenceEnthalpy.EnthalpyOfFormationAt25degC
           then Deltah0_f else 0) + h_offset + h_resid(T, p) - (if phase <>
          "gas" then h_resid(T, p0) else h_resid(
                T,
                p0,
                {1,-specVolPow[1]}))
          annotation (
          InlineNoEvent=true,
          Inline=true,
          smoothOrder=1);
        // The last two terms adjust for the actual pressure relative to the
        // reference.  In general, the lower limit of the integral of
        // (delh/delp)_T*dp is the reference pressure (p0).  However, if the
        // material is gaseous, then the reference is the corresponding ideal gas.
        // In that case, the lower limit of the real gas terms of the integral is
        // p=0, where a real gas behaves as an ideal gas.  See [Rao 1997, p. 271].
        annotation (Documentation(info="<html>
  <p>For an ideal gas, this function is independent of pressure
  (although pressure remains as a valid input).</p>
    </html>"));
      end h;

      replaceable function nu
        "<html>Dynamic **bulk viscosity as a function of temperature and specific volume (&Xi;)</html>"
        extends Modelica.Icons.Function;

        input Q.TemperatureAbsolute T=298.15*U.K "Temperature";
        input Q.PressureAbsolute p=U.atm "Pressure";
        // Note:  Pressure isn't used here but is included for generality.
        output Real nu "**Dynamic bulk viscosity";
        // **Dimension

      algorithm
        nu := m/alpha(T) annotation (Inline=true);

        // **note in doc here and/or in species model: Dynamic compressibility is a modified self diffusivity (2/2/13 notes)

        annotation (Documentation(info="<html>
<p>\"Dynamic compressibility\" is defined here as the reciprocal of the volume,
second, or bulk dynamic viscosity and specific volume (see
<a href=\"http://en.wikipedia.org/wiki/Volume_viscosity\">http://en.wikipedia.org/wiki/Volume_viscosity</a>).
</p>
</html>"));
      end nu;

      function p_Tv "Pressure as a function of temperature and specific volume"
        extends Modelica.Icons.Function;

        input Q.TemperatureAbsolute T=298.15*U.K "Temperature";
        input Q.VolumeSpecificAbsolute v=298.15*U.K/U.atm "Specific volume";
        output Q.PressureAbsolute p "Pressure";

      algorithm
        // assert(isCompressible,
        //  "The pressure is undefined since the material is incompressible.",
        //  AssertionLevel.warning);
        // Note:  In Dymola 7.4 the assertion level can't be set, although it has
        // been defined as an argument to assert() since Modelica 3.0.

        p := if isCompressible then Polynomial.f(
                v,
                {Polynomial.f(
                  T,
                  b_p[i, :],
                  pressPow[2]) for i in 1:size(b_p, 1)},
                pressPow[1]) else 0 annotation (
          Inline=true,
          inverse(v=v_Tp(T, p)),
          derivative=dp_Tv);
        annotation (Documentation(info="<html><p>If the species is incompressible, then <i>p</i>(<i>T</i>, <i>v</i>) is undefined,
  and the function will return a value of zero.</p>
  <p>The derivative of this function is <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.dp\">dp</a>().</p></html>"));
      end p_Tv;

      replaceable function R
        "Thermal resistivity as a function of temperature and pressure"

        extends Modelica.Icons.Function;

        input Q.TemperatureAbsolute T=298.15*U.K "Temperature";
        input Q.PressureAbsolute p=U.atm "Pressure";
        output Q.ResistivityThermal R "Thermal resistivity";

      algorithm
        R := alpha(T)/c_v(T, p) annotation (Inline=true);

      end R;

      function s "Specific entropy as a function of temperature and pressure"
        extends Modelica.Icons.Function;

        input Q.TemperatureAbsolute T=298.15*U.K "Temperature";
        input Q.PressureAbsolute p=U.atm "Pressure";
        output Q.NumberAbsolute s "Specific entropy";

      protected
        function s0_i
          "Return s0 as a function of T using one of the temperature intervals"
          input Q.TemperatureAbsolute T "Temperature";
          input Integer i "Index of the temperature interval";
          output Q.NumberAbsolute s0
            "Specific entropy at given temperature and reference pressure";

        algorithm
          s0 := Polynomial.F(
                    T,
                    b_c[i, :],
                    specHeatCapPow - 1) + B_c[i, 2]
            annotation (Inline=true, derivative=ds0_i);
          // This is the integral of c0_p/T*dT up to T at p0 with the absolute
          // entropy at the lower bound [McBride2002, p. 2].

        end s0_i;

        function ds0_i "Derivative of s0_i"
          // Note:  This function is necessary for Dymola 7.4 to differentiate s().

          input Q.TemperatureAbsolute T "Temperature";
          input Integer i "Index of the temperature interval";
          input Q.Temperature dT "Derivative of temperature";
          output Q.NumberAbsolute ds0
            "Derivative of specific entropy at given temperature and reference pressure";

        algorithm
          ds0 := Polynomial.f(
                    T,
                    b_c[i, :],
                    specHeatCapPow - 1)*dT annotation (Inline=true);

        end ds0_i;

        function s_resid
          "Residual specific entropy for pressure adjustment for selected rows of b_v"
          input Q.TemperatureAbsolute T "Temperature";
          input Q.PressureAbsolute p "Pressure";
          input Integer rowLimits[2]={1,size(b_v, 1)}
            "Beginning and ending indices of rows of b_v to be included";
          output Q.NumberAbsolute s_resid
            "Integral of (dels/delp)_T*dp up to p with zero integration constant (for selected rows)";

        algorithm
          s_resid := Polynomial.F(
                    p,
                    {Polynomial.f(
                      T,
                      b_v[i, :] .* {specVolPow[1] - specVolPow[2] + i - j for j
                 in 1:size(b_v, 2)},
                      specVolPow[2] - specVolPow[1] - i) for i in rowLimits[1]:
              rowLimits[2]},
                    specVolPow[1] + rowLimits[1] - 1) annotation (Inline=true);
          // Note:  According to the Maxwell relations, (dels/delp)_T = -(delv/delT)_p.

        end s_resid;

      algorithm
        /*
  assert(T_lim_c[1] <= T and T <= T_lim_c[size(T_lim_c, 1)], "Temperature " +
    String(T/U.K) + " K is out of bounds for " + name + " ([" + String(T_lim_c[1]
    /U.K) + ", " + String(T_lim_c[size(T_lim_c, 1)]/U.K) + "] K).");
  */
        // Note:  This is commented out so that the function can be inlined.
        // Note:  In Dymola 7.4 T_lim_c[size(T_lim_c, 1)] must be used
        // instead of T_lim_c[end] due to:
        //    "Error, not all 'end' could be expanded.";

        s := smooth(1, sum(if (T_lim_c[i] <= T or i == 1) and (T < T_lim_c[i +
          1] or i == size(b_c, 1)) then s0_i(T, i) else 0 for i in 1:size(b_c,
          1))) + s_resid(T, p) - (if phase <> "gas" then s_resid(T, p0) else
          s_resid(
                T,
                p0,
                {1,-specVolPow[1]}))
          annotation (
          InlineNoEvent=true,
          Inline=true,
          smoothOrder=1);
        // The first term gives the specific entropy at the given temperature and
        // reference pressure (p0).  The following terms adjust for the actual
        // pressure (p) by integrating (dels/delp)_T*dp.  In general, the
        // integration is from p0 to p.  However, for gases the reference state
        // is the ideal gas at the reference pressure.  Therefore, the lower
        // integration limit for the higher-order real gas terms (i.e., terms with
        // power of p greater than -1) is p=0.  This is pressure limit at which a
        // real gas behaves as an ideal gas, and it implies that those terms are
        // zero.  See [Rao1997, p. 272].  Note that the first V_m inside the curly
        // brackets of the related eq. 1.47 in [Dymond2002, p. 17] should be a
        // subscript rather than a multiplicative factor.
        // Note:  If the phase is gas, the virial equation of state (as defined by
        // b_v and specVolPow) must include an ideal gas term (v = ... + f(T)/p +
        // ...).  Otherwise, an indexing error will occur.

      end s;

      function v_Tp "Specific volume as a function of temperature and pressure"
        extends Modelica.Icons.Function;

        input Q.TemperatureAbsolute T=298.15*U.K "Temperature";
        input Q.PressureAbsolute p=U.atm "Pressure";
        output Q.VolumeSpecificAbsolute v "Specific volume";

      algorithm
        v := Polynomial.f(
                p,
                {Polynomial.f(
                  T,
                  b_v[i, :],
                  specVolPow[2] - specVolPow[1] - i + 1) for i in 1:size(b_v, 1)},
                specVolPow[1]) annotation (
          Inline=true,
          inverse(p=p_Tv(T, v)),
          derivative=dv_Tp);
        annotation (Documentation(info="<html>
  <p>The derivative of this function is
  <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.dv\">dv</a>().</p></html>"));
      end v_Tp;

      replaceable function Xi
        "<html>Dynamic compressibility as a function of temperature and specific volume (&Xi;)</html>"
        extends Modelica.Icons.Function;

        input Q.TemperatureAbsolute T=298.15*U.K "Temperature";
        input Q.PressureAbsolute p=U.atm "Pressure";
        // Note:  Pressure isn't used here but is included for generality.
        output Q.CompressibilityDynamic Xi "**Dynamic compressibility";
        // **Dimension: M.L2/(N.T) Weber?

      algorithm
        Xi := p*alpha(T) annotation (Inline=true);

        // **note in doc here and/or in species model: Dynamic compressibility is a modified self diffusivity (2/2/13 notes)

        annotation (Documentation(info="<html>
<p>\"Dynamic compressibility\" is defined here as the reciprocal of the volume,
second, or bulk dynamic viscosity and specific volume (see
<a href=\"http://en.wikipedia.org/wiki/Volume_viscosity\">http://en.wikipedia.org/wiki/Volume_viscosity</a>).
</p>
</html>"));
      end Xi;
      annotation (defaultComponentPrefixes="replaceable",Documentation(info="<html>
    <p>This package is compatible with NASA CEA thermodynamic data
    [<a href=\"modelica://FCSys.UsersGuide.References\">McBride2002</a>] and the virial equation of state
    [<a href=\"modelica://FCSys.UsersGuide.References\">Dymond2002</a>].  It may be used with
    the assumption of ideal gas or of constant specific volume, although it is more general than
    that.</p>

    <p>Assumptions:
    <ol><li>Specific mass is constant.</li></ol>
    </p>
    <p>The following notes apply to the constants:
    <ul>
    <li>Currently, <code>formula</code> may not contain parentheses or brackets.</li>
    <li><code>r</code> is the Van der Waals radius or the radius for the
    rigid-sphere (\"billiard-ball\") approximation of the kinetic theory of gases
    [<a href=\"modelica://FCSys.UsersGuide.References\">Present1958</a>].</li>
    <li><code>b_v</code>: The powers of <i>p</i>/<i>T</i> increase by row.  The powers of
    <i>T</i> increase by column.  If <code>specVolPow[1] == -1</code>, then the rows
    of <code>b_v</code> correspond to 1, <i>B</i><sup>*</sup><i>T</i>,
    <i>C</i><sup>*</sup><i>T</i><sup>2</sup>, <i>D</i><sup>*</sup><i>T</i><sup>3</sup>, &hellip;,
    where
    1, <i>B</i><sup>*</sup>, <i>C</i><sup>*</sup>, and <i>D</i><sup>*</sup> are
    the first, second, third, and fourth coefficients in the volume-explicit
    virial equation of state
    ([<a href=\"modelica://FCSys.UsersGuide.References\">Dymond2002</a>, pp. 1&ndash;2]).
    Currently,
    virial equations of state are supported up to the fourth coefficient (<i>D</i><sup>*</sup>).
    If additional terms are required, review and modify the definition of <code>b_p</code>.</li>
    <li><code>specVolPow</code> is defined as a <code>Real</code> vector.  However,
    special modifications are necessary if non-integer values are used
    (see <a href=\"modelica://FCSys.Characteristics.'e-'.Graphite\">'e-'.Graphite</a>).
    <li>The environment for <code>b_v</code> and <code>specVolPow</code> represent ideal gas.</li>
    <li><code>b_c</code>: The rows give the coefficients for the temperature intervals bounded
    by the values in <code>T_lim_c</code>.
    The powers of <i>T</i> increase
    by column.
    By default,
    the powers of <i>T</i> for the first column are each -2, which corresponds to [<a href=\"modelica://FCSys.UsersGuide.References\">McBride2002</a>].
    In that case, the dimensionalities of the coefficients are {L4.M2/(N2.T4), L2.M/(N.T2), 1, &hellip;}
    for each row, where L is length, M is mass, N is particle number, and T is time. (In <a href=\"modelica://FCSys\">FCSys</a>,
    temperature is a potential with dimension L2.M/(N.T2); see
    the <a href=\"modelica://FCSys.Units\">Units</a> package.)</li>
    <li><code>B_c</code>: As in <code>b_c</code>, the rows correspond to different
    temperature intervals.  The first column is for specific enthalpy and has dimensionality
    L2.M/(N.T2).  The second is for specific entropy and is dimensionless.
    The integration constants for enthalpy are defined such that the enthalpy at
    25 &deg;C is the specific enthalpy of formation at that temperature and reference pressure
    [<a href=\"modelica://FCSys.UsersGuide.References\">McBride2002</a>, p. 2].
    The integration constants for specific entropy are defined such that specific entropy is absolute.</li>
    <li><code>T_lim_c</code>: The first and last entries are the minimum and
    maximum valid temperatures.  The intermediate entries are the thresholds
    between rows of <code>b_c</code> (and <code>B_c</code>).  Therefore, if there are <i>n</i> temperature intervals
    (and rows in <code>b_c</code> and <code>B_c</code>), then <code>T_lim_c</code> must
    have <i>n</i> + 1 entries.</li>
    <li>The reference pressure is <code>p0</code>.   In the
    NASA CEA data [<a href=\"modelica://FCSys.UsersGuide.References\">McBride2002</a>], it is 1 bar for gases and 1 atm for condensed
    species.  For gases, the reference state is the ideal gas at <code>p0</code>.
    For example, the enthalpy of a non-ideal (real) gas at 25 &deg;C and <code>p0</code> with
    <code>ReferenceEnthalpy.ZeroAt25degC</code> selected is not exactly zero.
    </li>
    <li>If the <i>p</i>-<i>v</i>-<i>T</i> equation of state includes an ideal gas term (nonzero first virial coefficient), then the
    correlation for specific entropy (<a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.s\">s</a>) will involve the logarithm of pressure.
    The <code>p_min</code> parameter is used to guard against the logarithm of non-positive pressures.
    To disable this protection (and simplify the translated code), set <code>p_min</code> to zero or a
    negative pressure.</li>
    <li>If the material is gaseous (<code>phase == \"gas\"</code>), then the first virial coefficient
    must be independent of temperature.  Otherwise, the function for specific enthalpy
    (<a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.h\">h</a>) will be ill-posed.
    Typically the first virial coefficient is one (or equivalently <code>U.R</code>), which satisfies
    this requirement.</li>
    </ul>
    </p></html>"));

    end Characteristic;

    type ReferenceEnthalpy = enumeration(
        ZeroAt0K "Enthalpy at 0 K and p0 is 0 (if no additional offset)",
        ZeroAt25degC
          "Enthalpy at 25 degC and p0 is 0 (if no additional offset)",
        EnthalpyOfFormationAt25degC
          "Enthalpy at 25 degC and p0 is enthalpy of formation at 25 degC and p0 (if no additional offset)")
      "Enumeration for the reference enthalpy of a species";

  end BaseClasses;
  annotation (Documentation(info="<html>
  <p>Each species has a subpackage for each material phase in which the species
  is represented.  The thermodynamic properties are generally different for each phase.
  </p>
  <p>Additional materials may be included as needed.  The thermodynamic data for
  materials that are condensed at standard conditions is available in
  [<a href=\"modelica://FCSys.UsersGuide.References\">McBride2002</a>].
  The thermodynamic data for materials
  that are gases at standard conditions is available in
  <a href=\"modelica://Modelica.Media.IdealGases.Common.SingleGasesData\">Modelica.Media.IdealGases.Common.SingleGasesData</a>
  (and [<a href=\"modelica://FCSys.UsersGuide.References\">McBride2002</a>]). Virial coefficients are available in
  [<a href=\"modelica://FCSys.UsersGuide.References\">Dymond2002</a>].  Transport characteristics are available in
  [<a href=\"modelica://FCSys.UsersGuide.References\">McBride1996</a>,
  <a href=\"modelica://FCSys.UsersGuide.References\">McBride2002</a>].
  </p>
  <p>
<b>Licensed by the Georgia Tech Research Corporation under the Modelica License 2</b><br>
Copyright 2007&ndash;2013, Georgia Tech Research Corporation.
</p>
<p>
<i>This Modelica package is <u>free</u> software and the use is completely at <u>your own risk</u>;
it can be redistributed and/or modified under the terms of the Modelica License 2. For license conditions (including the
disclaimer of warranty) see <a href=\"modelica://FCSys.UsersGuide.ModelicaLicense2\">
FCSys.UsersGuide.ModelicaLicense2</a> or visit <a href=\"http://www.modelica.org/licenses/ModelicaLicense2\">
http://www.modelica.org/licenses/ModelicaLicense2</a>.</i>
</p></html>"));

end Characteristics;

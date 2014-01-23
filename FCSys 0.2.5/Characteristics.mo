within FCSys;
package Characteristics "Data and functions to correlate physical properties"
  import Modelica.Media.IdealGases.Common.FluidData;
  import Modelica.Media.IdealGases.Common.SingleGasesData;
  extends Modelica.Icons.Package;
  package Examples "Examples"
    extends Modelica.Icons.ExamplesPackage;
    model Correlations
      "Evaluate the material properties over varying temperature and pressure"
      extends Modelica.Icons.Example;

      parameter Q.TemperatureAbsolute T_start=273.15*U.K "Initial temperature"
        annotation (Dialog(__Dymola_label=
              "<html><i>T</i><sub>start</sub></html>"));
      parameter Q.TemperatureAbsolute T_stop=373.15*U.K "Final temperature"
        annotation (Dialog(__Dymola_label=
              "<html><i>T</i><sub>stop</sub></html>"));
      parameter Q.PressureAbsolute p_start=U.atm "Initial pressure" annotation
        (Dialog(__Dymola_label="<html><i>p</i><sub>start</sub></html>"));
      parameter Q.PressureAbsolute p_stop=U.atm "Final pressure" annotation (
          Dialog(__Dymola_label="<html><i>p</i><sub>stop</sub></html>"));

      // Property models
      PropertiesRT 'C+'(redeclare package Data = Characteristics.'C+'.Graphite)
        annotation (Placement(transformation(extent={{30,38},{50,58}})));
      PropertiesRT 'SO3-'(redeclare package Data =
            Characteristics.'SO3-'.Ionomer)
        annotation (Placement(transformation(extent={{30,26},{50,46}})));
      PropertiesRT 'e-'(redeclare package Data = Characteristics.'e-'.Graphite)
        annotation (Placement(transformation(extent={{30,14},{50,34}})));
      PropertiesRT 'H+'(redeclare package Data = Characteristics.'H+'.Ionomer)
        annotation (Placement(transformation(extent={{30,2},{50,22}})));
      PropertiesRT H2(redeclare package Data = FCSys.Characteristics.H2.Gas)
        annotation (Placement(transformation(extent={{30,-10},{50,10}})));
      PropertiesRT H2IG(redeclare package Data = FCSys.Characteristics.H2.Gas (
              b_v=[1],n_v={-1,0})) "H2 as ideal gas"
        annotation (Placement(transformation(extent={{30,-22},{50,-2}})));
      PropertiesRT H2O(redeclare package Data = FCSys.Characteristics.H2O.Gas)
        annotation (Placement(transformation(extent={{30,-34},{50,-14}})));
      // Note that H2O.p diverges from p in Dymola 2014 due to the large
      // coefficients in the second row of H2O.Data.b_v, which cause numerical
      // errors.
      PropertiesRT H2OLiquid(redeclare package Data =
            Characteristics.H2O.Liquid)
        annotation (Placement(transformation(extent={{30,-46},{50,-26}})));
      PropertiesRT N2(redeclare package Data = FCSys.Characteristics.N2.Gas)
        annotation (Placement(transformation(extent={{30,-58},{50,-38}})));
      PropertiesRT O2(redeclare package Data = FCSys.Characteristics.O2.Gas)
        annotation (Placement(transformation(extent={{30,-70},{50,-50}})));

      // Conditions

    protected
      Connectors.RealOutputInternal T(unit="L2.M/(N.T2)",displayUnit="K")
        "Temperature" annotation (Placement(transformation(extent={{-10,10},{10,
                30}}), iconTransformation(extent={{-10,16},{10,36}})));
      Connectors.RealOutputInternal p(unit="M/(L.T2)") "Pressure" annotation (
          Placement(transformation(extent={{-10,-30},{10,-10}}),
            iconTransformation(extent={{-10,-36},{10,-16}})));

    protected
      Modelica.Blocks.Sources.Ramp rampTemperature(
        final height=T_stop - T_start,
        final offset=T_start,
        duration=1)
        annotation (Placement(transformation(extent={{-40,10},{-20,30}})));
      Modelica.Blocks.Sources.Ramp rampPressure(
        final height=p_stop - p_start,
        final offset=p_start,
        duration=1)
        annotation (Placement(transformation(extent={{-40,-30},{-20,-10}})));

      model PropertiesRT
        "FCSys.Characteristics.Examples.Properties, with round-trip pressure calculation"
        extends FCSys.Characteristics.Examples.Properties;

        Q.PressureAbsolute p_RT=Data.p_Tv(T, v) if Data.isCompressible;

      end PropertiesRT;

    equation
      connect(rampTemperature.y, T) annotation (Line(
          points={{-19,20},{5.55112e-16,20}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(rampPressure.y, p) annotation (Line(
          points={{-19,-20},{5.55112e-16,-20}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(T, 'C+'.T) annotation (Line(
          points={{5.55112e-16,20},{20,20},{20,50},{29,50}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(p, 'C+'.p) annotation (Line(
          points={{5.55112e-16,-20},{10,-20},{10,46},{29,46}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(T, 'SO3-'.T) annotation (Line(
          points={{5.55112e-16,20},{20,20},{20,38},{29,38}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(p, 'SO3-'.p) annotation (Line(
          points={{5.55112e-16,-20},{10,-20},{10,34},{29,34}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(T, 'e-'.T) annotation (Line(
          points={{5.55112e-16,20},{20,20},{20,26},{29,26}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(p, 'e-'.p) annotation (Line(
          points={{5.55112e-16,-20},{10,-20},{10,22},{29,22}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(T, 'H+'.T) annotation (Line(
          points={{5.55112e-16,20},{20,20},{20,14},{29,14}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(p, 'H+'.p) annotation (Line(
          points={{5.55112e-16,-20},{10,-20},{10,10},{29,10}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(T, H2.T) annotation (Line(
          points={{5.55112e-16,20},{20,20},{20,2},{29,2}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(p, H2.p) annotation (Line(
          points={{5.55112e-16,-20},{10,-20},{10,-2},{29,-2}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(T, H2IG.T) annotation (Line(
          points={{5.55112e-16,20},{20,20},{20,-10},{29,-10}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(p, H2IG.p) annotation (Line(
          points={{5.55112e-16,-20},{10,-20},{10,-14},{29,-14}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(T, H2O.T) annotation (Line(
          points={{5.55112e-16,20},{20,20},{20,-22},{29,-22}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(p, H2O.p) annotation (Line(
          points={{5.55112e-16,-20},{10,-20},{10,-26},{29,-26}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(T, H2OLiquid.T) annotation (Line(
          points={{5.55112e-16,20},{20,20},{20,-34},{29,-34}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(p, H2OLiquid.p) annotation (Line(
          points={{5.55112e-16,-20},{10,-20},{10,-38},{29,-38}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(T, N2.T) annotation (Line(
          points={{5.55112e-16,20},{20,20},{20,-46},{29,-46}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(p, N2.p) annotation (Line(
          points={{5.55112e-16,-20},{10,-20},{10,-50},{29,-50}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(T, O2.T) annotation (Line(
          points={{5.55112e-16,20},{20,20},{20,-58},{29,-58}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(p, O2.p) annotation (Line(
          points={{5.55112e-16,-20},{10,-20},{10,-62},{29,-62}},
          color={0,0,127},
          smooth=Smooth.None));
      annotation (experiment, Commands(file=
              "Resources/Scripts/Dymola/Characteristics.Examples.Correlations.mos"
            "Characteristics.Examples.Correlations.mos"));
    end Correlations;

    model Properties
      "<html>Model that implements the functions of the <a href=\"modelica://FCSys.UsersGuide.Characteristic\">Characteristic</a> package</html>"
      extends FCSys.Icons.Blocks.ContinuousShort;

      replaceable package Data = Characteristics.BaseClasses.Characteristic
        constrainedby Characteristics.BaseClasses.Characteristic
        "Characteristic data" annotation (
        choicesAllMatching=true,
        __Dymola_choicesFromPackage=true,
        Placement(transformation(extent={{-60,40},{-40,60}}),
            iconTransformation(extent={{-10,90},{10,110}})));

      Connectors.RealInput T(unit="L2.M/(N.T2)") "Temperature"
        annotation (Placement(transformation(extent={{-120,10},{-100,30}})));
      Connectors.RealInput p(unit="M/(L.T2)") "Pressure"
        annotation (Placement(transformation(extent={{-120,-30},{-100,-10}})));

      Q.MassSpecific m "Specific mass";
      Q.LengthSpecific d "Specific diameter";
      Q.VolumeSpecific v "Specific volume";
      Q.CapacityThermalSpecific c_p "Isobaric specific heat capacity";
      Q.CapacityThermalSpecific c_v "Isobaric specific heat capacity";
      Q.Potential g "Gibbs potential";
      Q.Potential h "Specific enthalpy";
      Q.NumberAbsolute s "Specific entropy";
      Q.Resistivity eta "Fluidity";
      Q.ResistivityThermal theta "Thermal resistivity";
      Q.PressureReciprocal beta "Isothermal compressibility";
      Q.TimeAbsolute tauprime "Phase change interval";
      Q.Mobility mu "Mobility";
      Q.TimeAbsolute nu "Thermal independity";

    equation
      m = Data.m;
      d = Data.m;
      v = Data.v_Tp(T, p);
      c_p = Data.c_p(T, p);
      c_v = Data.c_v(T, p);
      g = Data.g(T, p);
      h = Data.h(T, p);
      s = Data.s(T, p);
      eta = Data.eta(T, v);
      theta = Data.theta(T, v);
      beta = Data.beta(T, p);
      tauprime = Data.tauprime(T, v);
      mu = Data.mu(T, v);
      nu = Data.nu(T, v);
      annotation (Diagram(graphics));
    end Properties;

    model SaturationPressure
      "<html>Evaluate the saturation pressure curve of H<sub>2</sub>O</html>"
      import FCSys.Characteristics.H2O.Liquid;
      import FCSys.Characteristics.H2O.Gas;
      package IdealGas = FCSys.Characteristics.H2O.Gas (b_v=[1], n_v={-1,0});

      extends Modelica.Icons.Example;

      Q.TemperatureAbsolute T "Temperature";
      Q.PressureAbsolute p_sat(start=U.kPa)
        "Saturation pressure via chemical equilibrium";
      Q.PressureAbsolute p_sat_IG(start=U.kPa)
        "Saturation pressure of ideal gas via chemical equilibrium";
      output Q.PressureAbsolute p_sat_MSL=Characteristics.H2O.p_sat(T)
        "Saturation pressure via Modelica.Media";
      output Q.Number T_degC=U.to_degC(T) "Temperature in degree Celsius";

      Modelica.Blocks.Sources.Ramp temperatureSet(
        height=99*U.K,
        duration=10,
        offset=274.15*U.K)
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

    equation
      T = temperatureSet.y;
      Liquid.g(T, p_sat) = Gas.g(T, p_sat);
      Liquid.g(T, p_sat_IG) = IdealGas.g(T, p_sat_IG);

      annotation (
        Documentation(info=
              "<html><p>See also <a href=\"modelica://FCSys.Subregions.Examples.PhaseChange.Evaporation\">Subregions.Examples.PhaseChange.Evaporation</a>.</p></html>"),

        experiment(StopTime=10),
        Commands(file(ensureTranslated=true) =
            "Resources/Scripts/Dymola/Characteristics.Examples.SaturationPressure.mos"
            "Characteristics.Examples.SaturationPressure.mos"));
    end SaturationPressure;

    model HydrationLevel
      "Evaluate the equilibrium hydration of the ionomer as a function of relative humidity"
      extends Modelica.Icons.Example;
      import Absorbed = FCSys.Characteristics.H2O.Ionomer;
      package Gas = FCSys.Characteristics.H2O.Gas (b_v=[1], n_v={-1,0});
      import Solid = FCSys.Characteristics.'SO3-'.Ionomer;
      import FCSys.Characteristics.H2O.p_sat;

      Q.Number RH "Relative humidity";
      output Q.Number lambda=Solid.v_Tp()/Absorbed.v_Tp(environment.T, p_abs)
        "Hydration due to chemical equilibrium";
      output Q.Number lambda_Springer=FCSys.Characteristics.H2O.lambda_eq(RH)
        if environment.analysis
        "Hydration according to Springer et al. correlation";
      Q.PressureAbsolute p_abs "H2O pressure in the ionomer";

      Modelica.Blocks.Sources.Ramp humiditySet(
        height=-0.999,
        offset=1,
        duration=10) "Set the relative humidity"
        annotation (Placement(transformation(extent={{-10,-40},{10,-20}})));

      inner Conditions.Environment environment(T=303.15*U.K)
        annotation (Placement(transformation(extent={{-10,10},{10,30}})));

    equation
      RH = humiditySet.y;
      Absorbed.g(environment.T, p_abs) = Gas.g(environment.T, RH*p_sat(
        environment.T));

      annotation (
        Documentation(info=
              "<html><p>See also <a href=\"modelica://FCSys.Subregions.Examples.PhaseChange.Hydration\">Subregions.Examples.PhaseChange.Hydration</a>.</p></html>"),

        experiment(StopTime=10),
        Commands(file=
              "Resources/Scripts/Dymola/Characteristics.Examples.HydrationLevel.mos"
            "Characteristics.Examples.HydrationLevel.mos"),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics));
    end HydrationLevel;

    model CellPotential
      "<html>Evaluate the potential of an H<sub>2</sub>/O<sub>2</sub> cell as a function of temperature</html>"
      import FCSys.Characteristics.H2O.Liquid;
      import H2O = FCSys.Characteristics.H2O.Gas;
      import H2 = FCSys.Characteristics.H2.Gas;
      import O2 = FCSys.Characteristics.O2.Gas;
      import H2OLiquid = FCSys.Characteristics.H2O.Liquid;
      package H2OIG = FCSys.Characteristics.H2O.Gas (b_v=[1], n_v={-1,0});
      package O2IG = FCSys.Characteristics.O2.Gas (b_v=[1], n_v={-1,0});
      package H2IG = FCSys.Characteristics.H2.Gas (b_v=[1], n_v={-1,0});

      extends Modelica.Icons.Example;

      output Real T_degC=U.to_degC(T) "Temperature in deg C";
      output Q.Potential w_gas=0.5*H2.g(T, environment.p_dry) + 0.25*O2.g(T,
          environment.p_O2) - 0.5*H2O.g(T, environment.p_H2O)
        "Cell potential with H2O as gas";
      output Q.Potential w_IG=0.5*H2IG.g(T, environment.p_dry) + 0.25*O2IG.g(T,
          environment.p_O2) - 0.5*H2OIG.g(T, environment.p_H2O)
        "Cell potential with ideal gases";
      output Q.Potential w_liq=0.5*H2.g(T, environment.p_dry) + 0.25*O2.g(T,
          environment.p_O2) - 0.5*H2OLiquid.g(T, environment.p)
        "Cell potential with H2O as liquid";

      Q.TemperatureAbsolute T "Temperature";

      Modelica.Blocks.Sources.Ramp temperatureSet(
        height=99*U.K,
        duration=10,
        offset=274.15*U.K)
        annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

      inner Conditions.Environment environment
        annotation (Placement(transformation(extent={{-10,0},{10,20}})));

    equation
      T = temperatureSet.y;

      annotation (experiment(StopTime=10), Commands(file=
              "Resources/Scripts/Dymola/Characteristics.Examples.CellPotential.mos"
            "Characteristics.Examples.CellPotential.mos"));
    end CellPotential;

    model Leverett
      "<html>Evaluate the Leverett J function from [<a href=\"modelica://FCSys.UsersGuide.References.Wang2001\">Wang2001</a>]</html>"
      extends Modelica.Icons.Example;

      output Q.NumberAbsolute s=saturationSet.y "Liquid saturation";
      output Q.NumberAbsolute J=FCSys.Characteristics.H2O.J(s)
        "Result of Leverett correlation";

      Modelica.Blocks.Sources.Ramp saturationSet(height=1, duration=1)
        "Set the saturation"
        annotation (Placement(transformation(extent={{-10,-40},{10,-20}})));

      inner Conditions.Environment environment(T=303.15*U.K)
        annotation (Placement(transformation(extent={{-10,10},{10,30}})));

      annotation (
        experiment,
        Documentation(info=
              "<html><p>Please see <a href=\"modelica://FCSys.Characteristics.H2O.J\">H2O.J</a>().</p></html>"),

        Commands(file=
              "Resources/Scripts/Dymola/Characteristics.Examples.Leverett.mos"
            "Characteristics.Examples.Leverett.mos"),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics));

    end Leverett;

    model LatentHeat
      "<html>Evaluate the latent heat of vaporization of H<sub>2</sub>O</html>"
      import FCSys.Characteristics.H2O.Liquid;
      import FCSys.Characteristics.H2O.Gas;

      extends Modelica.Icons.Example;

      parameter Q.PressureAbsolute p(displayUnit="atm") = U.atm
        "Pressure of the liquid (and total pressure of the gas)";
      Q.PressureAbsolute p_sat(displayUnit="atm") "Saturation pressure";
      Q.TemperatureAbsolute T "Temperature";
      output Q.Number T_degC=U.to_degC(T) "Temperature in degree Celsius";
      output Q.Potential h_g(displayUnit="J/mol") = Gas.h(T, p_sat)
        "Specific enthalpy of the saturated vapor";
      output Q.Potential h_l(displayUnit="J/mol") = Liquid.h(T, p)
        "Specific enthalpy of the liquid";
      output Q.Potential h_gl(displayUnit="J/mol") = h_g - h_l
        "Specific enthalpy of vaporization";
      output Q.Velocity2 hbar_gl(displayUnit="J/g") = h_gl/Liquid.m
        "Massic enthalpy of vaporization";

      Modelica.Blocks.Sources.Ramp temperatureSet(
        height=99*U.K,
        duration=10,
        offset=274.15*U.K)
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

    equation
      T = temperatureSet.y;
      Gas.g(T, p_sat) = Liquid.g(T, p);

      annotation (
        Documentation(info=
              "<html><p>See also <a href=\"modelica://FCSys.Subregions.Examples.PhaseChange.Condensation\">Subregions.Examples.PhaseChange.Condensation</a>.</p></html>"),

        experiment(StopTime=10),
        Commands(file(ensureTranslated=true) =
            "Resources/Scripts/Dymola/Characteristics.Examples.LatentHeat.mos"
            "Characteristics.Examples.LatentHeat.mos"));
    end LatentHeat;

    model MobilityFactors "Test the mobility factors"
      extends Modelica.Icons.Example;
      import FCSys.Characteristics.MobilityFactors.*;

      output Q.NumberAbsolute H2_H2O=k_H2_H2O() "Between H2 and H2O";
      output Q.NumberAbsolute H2O_N2=k_H2O_N2() "Between H2O and N2";
      output Q.NumberAbsolute H2O_O2=k_H2O_O2() "Between H2O and O2";
      output Q.NumberAbsolute N2_O2=k_N2_O2() "Between N2 and O2";
      annotation (Commands(file=
              "Resources/Scripts/Dymola/Characteristics.Examples.MobilityFactors.mos"
            "Characteristics.Examples.MobilityFactors.mos"));
    end MobilityFactors;

    model SurfaceTension
      "<html>Evaluate the surface tension of H<sub>2</sub>O using the model of Garai</html>"
      import FCSys.Characteristics.H2O.Liquid;
      import FCSys.Characteristics.H2O.Gas;

      extends Modelica.Icons.Example;

      parameter Q.PressureAbsolute p(displayUnit="atm") = U.atm
        "Pressure of the liquid (and total pressure of the gas)";
      Q.PressureAbsolute p_sat(displayUnit="atm") "Saturation pressure";
      Q.TemperatureAbsolute T "Temperature";
      Q.TemperatureAbsolute T_sat(start=373.15*U.K) "Saturation temperature";
      output Q.Number T_degC=U.to_degC(T) "Temperature in degree Celsius";
      output Q.Potential h_gl(displayUnit="J/mol") = Gas.h(T, p_sat) - Liquid.h(
        T, p) "Specific enthalpy of vaporization";
      output Q.SurfaceTension gamma=(h_gl - T_sat)/(2*U.pi*Liquid.d^2*U.q)
        "Surface tension";

      Modelica.Blocks.Sources.Ramp temperatureSet(
        height=99*U.K,
        duration=10,
        offset=274.15*U.K)
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

    equation
      T = temperatureSet.y;
      Gas.g(T_sat, p) = Liquid.g(T_sat, p);
      Gas.g(T, p_sat) = Liquid.g(T, p);

      annotation (
        Documentation(info=
              "<html><p>See also <a href=\"modelica://FCSys.Subregions.Examples.PhaseChange.Condensation\">Subregions.Examples.PhaseChange.Condensation</a>.</p></html>"),

        experiment(StopTime=10),
        Commands);
    end SurfaceTension;

  end Examples;

  package 'C+' "<html>C<sup>+</sup></html>"
    extends Modelica.Icons.Package;
    package Graphite "C+ graphite"
      // Note:  HTML formatting isn't used in the description because Dymola 2014
      // doesn't support it in the GUI for replaceable lists.  The same applies
      // to other species.

      extends BaseClasses.Characteristic(
        final formula="C+",
        final phase=Phase.solid,
        p0=U.atm,
        n_v={0,0},
        b_v=[U.cc*m/(2.2210*U.g)],
        m=12.0107000*U.g/U.mol - 'e-'.Gas.m,
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
            5.848134850e3, -2.350925275e1] - b_c[:, 2:3]*log(U.K),
        d=340*U.pico*U.m/U.q);
      annotation (Documentation(info="<html>

     <p>Assumptions:</p>
     <ol>
     <li>Constant specific volume (i.e., incompressible and without
          thermal expansion)</li>
     </ol>

     <p>Additional notes:</p>
     <ul>
     <li>The data for this species is for C rather than C<sup>+</sup>, with the exception of specific mass.</li>
     <li>The radius is from <a href=\"http://en.wikipedia.org/wiki/Carbon\">http://en.wikipedia.org/wiki/Carbon</a>.  See also
   <a href=\"http://en.wikipedia.org/wiki/Van_der_Waals_radius\">http://en.wikipedia.org/wiki/Van_der_Waals_radius</a>.</li>
     <li>The default specific volume (<code>v = U.cc*m/(2.210*U.g)</code>) is of pyrolytic graphite
  at 300&nbsp;K according to [<a href=\"modelica://FCSys.UsersGuide.References.Incropera2002\">Incropera2002</a>, p.&nbsp;909].  Other forms
  are (Ibid., also at 300&nbsp;K) are:
  <ul>
       <li>Amorphous carbon:  <code>v = U.cc*m/(1.950*U.g)</code></li>
       <li>Diamond (type IIa):  <code>v = U.cc*m/(3.500*U.g)</code></li>
       </ul>
     </ul>

<p>For more information, please see the
  <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic\">Characteristic</a> package.</p></html>"));

    end Graphite;

  end 'C+';

  package 'SO3-'
    "<html>C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S<sup>-</sup> (abbreviated as SO<sub>3</sub><sup>-</sup>)</html>"
    extends Modelica.Icons.Package;
    package Ionomer "C9HF17O5S- ionomer"

      extends BaseClasses.Characteristic(
        final formula="C9HF17O5S-",
        final phase=Phase.solid,
        p0=U.atm,
        m=1044.214*U.g/U.mol - 'H+'.Gas.m,
        n_v={0,0},
        b_v=[U.cc*m/(2.00*U.g)],
        Deltah0_f=0,
        Deltah0=0,
        n_c=0,
        T_lim_c={0,Modelica.Constants.inf},
        b_c=[4188*U.J*m/(U.kg*U.K)],
        B_c=[Deltah0_f - 298.15*U.K*b_c[1, 1], 0],
        d=(294 + 2259.8)*U.pico*U.m/U.q);

      annotation (Documentation(info="<html>
       <p>Assumptions:</p>
     <ol>
     <li>Constant specific volume (i.e., incompressible and without
          thermal expansion)</li>
     </ol>

    <p>The specific volume (<code>v = U.cc*m/(2.00*U.g)</code>) is based on
   [<a href=\"modelica://FCSys.UsersGuide.References.Lin2006\">Lin2006</a>, p.&nbsp;A1327].  Note that this is

   approximately 1.912 M, which does not match

   the default density of

  <a href=\"modelica://FCSys.Characteristics.'H+'.Ionomer\">'H+'.Ionomer</a> (0.95 M), but it
  simplifies the model by requiring only C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S<sup>-</sup>
  (not C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S) for charge neutrality.</p>

     <p>Additional notes:</p>
     <ul>
     <li>Most of the data for this species is for C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S rather than
     C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S<sup>-</sup> (with the exception of specific mass).</li>
     <li>A form of C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S is
     C<sub>7</sub>HF<sub>13</sub>O<sub>5</sub>S.(C<sub>2</sub>F<sub>4</sub>)<sub>6</sub>, which is a typical
   configuration of Nafion sulfonate (after hydrolysis)
   [<a href=\"modelica://FCSys.UsersGuide.References.Mark1999\">Mark1999</a>, p.&nbsp;234].</li>
     <li>Thermodynamic data for this material is not available from
     [<a href=\"modelica://FCSys.UsersGuide.References.McBride2002\">McBride2002</a>].  The default specific heat capacity
     (<i>b<sub>c</sub></i> = <code>[4188*U.J*m/(U.kg*U.K)]</code>) is based on [<a href=\"modelica://FCSys.UsersGuide.References.Shah2009\">Shah2009</a>, p.&nbsp;B472].</li>
     <li>According to [<a href=\"modelica://FCSys.UsersGuide.References.Avogadro\">Avogadro</a>], the furthest distance
   between two atoms of C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S is 2259.8 pm and is between fluorines.
   The radius of F is 147 pm (<a href=\"http://en.wikipedia.org/wiki/Fluorine\">http://en.wikipedia.org/wiki/Fluorine</a>).</li>
     <li>From [<a href=\"modelica://FCSys.UsersGuide.References.Avogadro\">Avogadro</a>], the molecular weight of C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S
   is 1044.214 g/mol.  However, the \"11\" in Nafion 11x indicates a
   molecular weight of 1100 g/mol.  According to
   <a href=\"http://en.wikipedia.org/wiki/Nafion\">http://en.wikipedia.org/wiki/Nafion</a>,
       \"the molecular weight of Nafion is uncertain due to differences in
        processing and solution morphology.\"</li>
     </ul>

<p>For more information, please see the
  <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic\">Characteristic</a> package.</p></html>"));

    end Ionomer;

  end 'SO3-';

  package 'e-' "<html>e<sup>-</sup></html>"
    extends Modelica.Icons.Package;
    package Gas "e- gas"
      import Data = Modelica.Media.IdealGases.Common.SingleGasesData.eminus;

      extends BaseClasses.Characteristic(
        final formula="e-",
        final z=-1,
        phase=Phase.gas,
        m=Data.MM*U.kg/U.mol,
        Deltah0_f=Data.MM*Data.Hf*U.J/U.mol,
        Deltah0=Data.MM*Data.H0*U.J/U.mol,
        T_lim_c={200.000,20000.000}*U.K,
        b_c={Data.alow} .* fill({U.K^(3 - i) for i in 1:size(Data.alow, 1)},
            size(T_lim_c, 1) - 1),
        B_c={Data.blow} .* fill({U.K,1}, size(T_lim_c, 1) - 1) - b_c[:, 2:3]*
            log(U.K),
        d=U.alpha^3/(2*U.pi*U.R_inf*U.q));

      annotation (Documentation(info="<html>
     <p>Notes:</p>
     <ul>
     <li>The specific mass (<i>m</i>) is also

     2<i>k</i><sub>A</sub>/<i>d</i>, where <i>k</i><sub>A</sub> (in the <a href=\"modelica://FCSys.Units\">Units</a> package) is the magnetic force constant and

     <i>d</i> (in this package) is the specific classical diameter of an electron,

     &alpha;<sup>3</sup>/(2&pi;&nbsp;<i>R</i><sub>&infin;</sub>&nbsp;<i>q</i>).</li>
  <li>McBride and Gordon [<a href=\"modelica://FCSys.UsersGuide.References.McBride1996\">McBride1996</a>] provide correlations for the transport
  properties of e<sup>-</sup> gas.  However, they are not entered here, since they
  contain only one temperature range (2000 to 5000&nbsp;K) which is beyond the expected operating range of the model.</li>
     <li>The equation for the radius is the classical radius of an electron (see
  <a href=\"http://en.wikipedia.org/wiki/Classical_electron_radius\">http://en.wikipedia.org/wiki/Classical_electron_radius</a>).</li>
  <li>McBride and Gordon [<a href=\"modelica://FCSys.UsersGuide.References.McBride1996\">McBride1996</a>] provide correlations for the transport
  properties of e<sup>-</sup> gas.  However, they are not entered here, since they
  contain only one temperature range (2000 to 5000&nbsp;K) which is beyond the expected operating range of the model.</li>
  <li>
  The thermodynamic data [<a href=\"modelica://FCSys.UsersGuide.References.McBride2002\">McBride2002</a>] splits the correlations into three
  temperature ranges, but the coefficients are the same for each.
   Therefore, the temperature limits are set here such that the entire
   range is handled without switching.  The lower temperature limit
   in the source data is 298.150&nbsp;K, but here it is expanded down to 200&nbsp;K.
   The constants are independent of temperature anyway.
  </li>
</ul>

<p>For more information, please see the
  <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic\">Characteristic</a> package.</p></html>"));

    end Gas;

    package Graphite "e- in graphite"
      extends Gas(n_v={0,0}, b_v='C+'.Graphite.b_v);
      annotation (Documentation(info="<html><p>Assumptions:</p><ol>
      <li>There is a 1:1 ratio of free (conductance band) electrons and carbon atoms.  The density of the carbon is set by

      <a href=\"modelica://FCSys.Characteristics.'C+'.Graphite\">'C+'.Graphite</a>.</li></ol></html>"));

    end Graphite;

  end 'e-';

  package 'H+' "<html>H<sup>+</sup></html>"
    extends Modelica.Icons.Package;
    package Gas "H+ gas"
      import Data = Modelica.Media.IdealGases.Common.SingleGasesData.Hplus;

      extends BaseClasses.Characteristic(
        final formula="H+",
        final z=1,
        phase=Phase.gas,
        final m=Data.MM*U.kg/U.mol,
        referenceEnthalpy=ReferenceEnthalpy.zeroAt25degC,
        Deltah0_f=Data.MM*Data.Hf*U.J/U.mol,
        Deltah0=Data.MM*Data.H0*U.J/U.mol,
        n_c=-2,
        T_lim_c={200.000,20000.000}*U.K,
        b_c={Data.alow} .* fill({U.K^(3 - i) for i in 1:size(Data.alow, 1)},
            size(T_lim_c, 1) - 1),
        B_c={Data.blow} .* fill({U.K,1}, size(T_lim_c, 1) - 1) - b_c[:, 2:3]*
            log(U.K),
        d=240*U.pico*U.m/U.q);

      annotation (Documentation(info="<html>
         <p>Assumptions:</p>
     <ol>
  <li>The radius (for the purpose of the rigid-sphere assumption of
  kinetic theory) is that of H
  (<a href=\"http://en.wikipedia.org/wiki/Hydrogen\">http://en.wikipedia.org/wiki/Hydrogen</a>).</li>
     </ol>

     <p>Additional notes:</p>
     <ul>
     <li>The default reference enthalpy is zero at 25&nbsp;&deg;C.</li>
     <li>The source data [<a href=\"modelica://FCSys.UsersGuide.References.McBride2002\">McBride2002</a>] breaks the

     data into three
   temperature ranges, but the constants are the same for each.
   Therefore, the temperature limits are set here such that the entire
   range is handled without switching.
   The lower temperature limit in the source data
   is 298.150&nbsp;K, but here it is expanded down to 200&nbsp;K.  The constants are
   independent of temperature anyway.</li>
   <li>The enthalpy to produce H<sup>+</sup> from H<sub>2</sub>O

   (H<sub>2</sub>O&nbsp;&#8640;&nbsp;OH<sup>-</sup>&nbsp;+&nbsp;H<sup>+</sup>) is
   -96569.804&nbsp;J/mol.  The enthalpy to produce H<sup>+</sup> from H<sub>3</sub>O<sup>+</sup>
   (H<sub>3</sub>O<sup>+</sup>&nbsp;&#8640;&nbsp;H<sub>2</sub>O&nbsp;+&nbsp;H<sup>+</sup>) is

   839826.0&nbsp;J/mol.  Based on
   [<a href=\"modelica://FCSys.UsersGuide.References.Tissandier1998\">Tissandier1998</a>], the
   enthalpy of formation of aqueous H<sup>+</sup> is 1150.1e3&nbsp;J/mol.</li>
     </ul>

<p>For more information, please see the
  <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic\">Characteristic</a> package.</p></html>"));

    end Gas;

    package Ionomer "H+ in ionomer"
      extends Gas(n_v={0,0}, b_v=[1/(0.95*U.M)]);
      annotation (Documentation(info="<html><p>The initial density corresponds to the measurement
  by Spry and Fayer (0.95 M) in Nafion<sup>&reg;</sup> at
  &lambda; = 12, where &lambda; is the number of
  H<sub>2</sub>O molecules to SO<sub>3</sub>H
  endgroups.  At &lambda; = 22, the density was measured at 0.54 M
  [<a href=\"modelica://FCSys.UsersGuide.References.Spry2009\">Spry2009</a>].</p></html>"));

    end Ionomer;

  end 'H+';

  package H2 "<html>H<sub>2</sub></html>"
    extends Modelica.Icons.Package;
    package Gas "H2 gas"
      import Data = Modelica.Media.IdealGases.Common.SingleGasesData.H2;

      extends BaseClasses.CharacteristicNASA(
        final formula=Data.name,
        final phase=Phase.gas,
        final m=Data.MM*U.kg/U.mol,
        n_v={-1,-3},
        b_v={{0,0,0,1},{8.0282e6*U.K^3,-2.6988e5*U.K^2,-129.26*U.K,17.472}*U.cm
            ^3/U.mol},
        Deltah0_f=Data.MM*Data.Hf*U.J/U.mol,
        Deltah0=Data.MM*Data.H0*U.J/U.mol,
        T_lim_c={200.000,Data.Tlimit,6000.000,20000.000}*U.K,
        b_c={Data.alow,Data.ahigh,{4.966884120e8,-3.147547149e5,79.84121880,-8.414789210e-3,
            4.753248350e-7,-1.371873492e-11,1.605461756e-16}} .* fill({U.K^(3
             - i) for i in 1:size(Data.alow, 1)}, size(T_lim_c, 1) - 1),
        B_c={Data.blow,Data.bhigh,{2.488433516e6,-669.5728110}} .* fill({U.K,1},
            size(T_lim_c, 1) - 1) - b_c[:, 2:3]*log(U.K),
        d=(240 + 100.3)*U.pico*U.m/U.q,
        T_lim_eta_theta={200.0,1000.0,5000.0,15000.0}*U.K,
        b_eta={fromNASAViscosity({0.74553182,43.555109,-3.2579340e3,0.13556243}),
            fromNASAViscosity({0.96730605,679.31897,-2.1025179e5,-1.8251697}),
            fromNASAViscosity({1.0126129,1.4973739e3,-1.4428484e6,-2.3254928})},

        b_theta={fromNASAThermalConductivity({1.0059461,279.51262,-2.9792018e4,
            1.1996252}),fromNASAThermalConductivity({1.0582450,248.75372,
            1.1736907e4,0.82758695}),fromNASAThermalConductivity({-0.22364420,-6.9650442e3,
            -7.7771313e4,13.189369})});

      annotation (Documentation(info="<html>
            <p>Notes:</p>
     <ul>
  <li>According to [<a href=\"modelica://FCSys.UsersGuide.References.Avogadro\">Avogadro</a>], the (center-to-center)
   bond length of H-H is 100.3 pm.  The radius of H is from
   <a href=\"http://en.wikipedia.org/wiki/Hydrogen\">http://en.wikipedia.org/wiki/Hydrogen</a>.  See also
   <a href=\"http://en.wikipedia.org/wiki/Van_der_Waals_radius\">http://en.wikipedia.org/wiki/Van_der_Waals_radius</a>.</li>
  <li>The virial coefficients are from [<a href=\"modelica://FCSys.UsersGuide.References.Dymond2002\">Dymond2002</a>, p.&nbsp;41].  The
  temperature range of the coefficients is [60, 500] K, but this is not enforced in the functions.</li>
     </ul>

<p>For more information, please see the
  <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic\">Characteristic</a> package.</p></html>"));

    end Gas;

  end H2;

  package H2O "<html>H<sub>2</sub>O</html>"
    extends Modelica.Icons.Package;
    package Gas "H2O gas"
      import Data = Modelica.Media.IdealGases.Common.SingleGasesData.H2O;

      extends BaseClasses.CharacteristicNASA(
        final formula=Data.name,
        final phase=Phase.gas,
        final m=Data.MM*U.kg/U.mol,
        n_v={-1,-3},
        b_v={{0,0,0,1},{-5.6932e10*U.K^3,1.8189e8*U.K^2,-3.0107e5*U.K,158.83}*U.cm
            ^3/U.mol},
        Deltah0_f=Data.MM*Data.Hf*U.J/U.mol,
        Deltah0=Data.MM*Data.H0*U.J/U.mol,
        T_lim_c={200.000,Data.Tlimit,6000.000}*U.K,
        b_c={Data.alow,Data.ahigh} .* fill({U.K^(3 - i) for i in 1:size(Data.alow,
            1)}, size(T_lim_c, 1) - 1),
        B_c={Data.blow,Data.bhigh} .* fill({U.K,1}, size(T_lim_c, 1) - 1) - b_c[
            :, 2:3]*log(U.K),
        d=282*U.pico*U.m/U.q,
        T_lim_eta_theta={373.2,1073.2,5000.0,15000.0}*U.K,
        b_eta={fromNASAViscosity({0.50019557,-697.12796,8.8163892e4,3.0836508}),
            fromNASAViscosity({0.58988538,-537.69814,5.4263513e4,2.3386375}),
            fromNASAViscosity({0.64330087,-95.668913,-3.7742283e5,1.8125190})},

        b_theta={fromNASAThermalConductivity({1.0966389,-555.13429,1.0623408e5,
            -0.24664550}),fromNASAThermalConductivity({0.39367933,-2.2524226e3,
            6.1217458e5,5.8011317}),fromNASAThermalConductivity({-0.41858737,-1.4096649e4,
            1.9179190e7,14.345613})});

      annotation (Documentation(info="<html>
        <p>Notes:</p>
     <ul>
  <li>The radius of H<sub>2</sub>O is 282 pm
   (<a href=\"http://www.lsbu.ac.uk/water/molecule.html\">http://www.lsbu.ac.uk/water/molecule.html</a>).  Using the radius of H
   from <a href=\"http://en.wikipedia.org/wiki/Hydrogen\">http://en.wikipedia.org/wiki/Hydrogen</a> and the center-to-center
   distance of hydrogen atoms in H<sub>2</sub>O from [<a href=\"modelica://FCSys.UsersGuide.References.Avogadro\">Avogadro</a>],
   156.6 pm, the radius of H<sub>2</sub>O would be (120 + 156.6/2) pm = 198.3 pm.</li>
  <li>The virial coefficients are from [<a href=\"modelica://FCSys.UsersGuide.References\">Dymond2002</a>, p.&nbsp;4].  The
  temperature range of the coefficients is [350, 770] K, but this is not enforced in the functions.</li>
     </ul>

  <p>For more information, please see the
  <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic\">Characteristic</a> package.</p></html>"));

    end Gas;

    package Ionomer "H2O in ionomer"
      extends Gas(
        b_v=[1],
        n_v={-1,0},
        b_c=Gas.b_c + fill({0,0,1,-2*0.2580123533308264/U.K,6*
            4.841882910380711e-4/U.K^2,-12*3.328156493413594e-7/U.K^3,0}, 2),
        B_c=Gas.B_c + fill({0,47.15136731353458 - log('SO3-'.Ionomer.b_v[1, 1]*
            U.atm/14)}, 2));
      // These coefficients are based on the saturation pressure correlation
      // (Eq. 15) from Springer1991.  The factor of 1/14 in the 2nd column of
      // B_c gives lambda = 14 in equilibrium with saturated vapor.

      annotation (Documentation(info="<html><p>The thermodynamic data is set such that the
  density in equilibrium with H<sub>2</sub>O vapor will match the Springer et al.

   hydration curve
  [<a href=\"modelica://FCSys.UsersGuide.References.Springer1991\">Springer1991</a>] (see <a href=\"modelica://FCSys.Characteristics.H2O.lambda_eq\">lambda_eq</a>())
  at &lambda; = 14 and &lambda; = 0.
  Otherwise, the properties are the same as H<sub>2</sub>O as an ideal gas.</p>

  <p>For more information, please see the
  <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic\">Characteristic</a> package.</p></html>"));

    end Ionomer;

    package Liquid "H2O liquid"

      extends BaseClasses.Characteristic(
        final formula="H2O",
        final phase=Phase.liquid,
        final m=0.01801528*U.kg/U.mol,
        p0=U.atm,
        n_v={0,0},
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
             - b_c[:, 2:3]*log(U.K),
        d=282*U.pico*U.m/U.q);
      annotation (Documentation(info="<html>     <p>Assumptions:</p>
     <ol>
     <li>Constant specific volume (i.e., incompressible and without
          thermal expansion)</li>
     </ol>

<p>Additional notes:</p>
     <ul>
     <li>See <a href=\"modelica://FCSys.Characteristics.H2O.Gas\">Characteristics.H2O.Gas</a> regarding the radius.</li>
     <li>The default specific volume (<i>b<sub>v</sub></i> = <code>[U.cc*m/(0.99656*U.g)]</code>) is at 300&nbsp;K based on [<a href=\"modelica://FCSys.UsersGuide.References.Takenaka1990\">Takenaka1990</a>].</li>
     </ul>

<p>For more information, please see the
  <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic\">Characteristic</a> package.</p></html>"));

    end Liquid;

    function p_sat
      "<html>Saturation pressure (<i>p</i><sub>sat</sub>) as a function of temperature</html>"
      extends Modelica.Icons.Function;

      input Q.TemperatureAbsolute T "Temperature"
        annotation (Dialog(__Dymola_label="<html><i>T</i></html>"));
      output Q.PressureAbsolute p_sat "Saturation pressure" annotation (Dialog(
            __Dymola_label="<html><i>p</i><sub>sat</sub></html>"));

    algorithm
      p_sat := Modelica.Media.Air.MoistAir.saturationPressureLiquid(T/U.K)*U.Pa;
      annotation (Inline=true);
    end p_sat;

    function lambda_eq
      "<html>Equilibrium hydration of ionomer in contact with vapor (&lambda;<sub>eq</sub>) as a function of relative humidity</html>"
      extends Modelica.Icons.Function;

      input Q.NumberAbsolute RH "Relative humidity";
      output Real lambda
        "<html>Mole ratio of H<sub>2</sub>O to SO<sub>3</sub><sup>-</sup></html>"
        annotation (Dialog(__Dymola_label="<html>&lambda;<sub>eq</sub></html>"));

    algorithm
      lambda := 0.043 + 17.81*RH - 39.85*RH^2 + 36*RH^3;

      annotation (Inline=true,Documentation(info="<html><p>This implements the correlation by Springer et al. [<a href=\"modelica://FCSys.UsersGuide.References.Springer1991\">Springer1991</a>]
  for the ratio of H<sub>2</sub>O molecules to SO<sub>3</sub><sup>-</sup> units of
  Nafion&reg; EW 1100 series ionomer.</p></html>"));
    end lambda_eq;

  end H2O;

  package N2 "<html>N<sub>2</sub></html>"
    extends Modelica.Icons.Package;
    package Gas "N2 gas"
      import Data = Modelica.Media.IdealGases.Common.SingleGasesData.N2;

      extends BaseClasses.CharacteristicNASA(
        final formula=Data.name,
        final phase=Phase.gas,
        final m=Data.MM*U.kg/U.mol,
        n_v={-1,-4},
        b_v={{0,0,0,0,1},{-2.7198e9*U.K^4,6.1253e7*U.K^3,-1.4164e6*U.K^2,-9.3378e3
            *U.K,40.286}*U.cc/U.mol},
        Deltah0_f=Data.MM*Data.Hf*U.J/U.mol,
        Deltah0=Data.MM*Data.H0*U.J/U.mol,
        T_lim_c={200.000,Data.Tlimit,6000.000,20000.000}*U.K,
        b_c={Data.alow,Data.ahigh,{8.310139160e8,-6.420733540e5,2.020264635e2,-3.065092046e-2,
            2.486903333e-6,-9.705954110e-11,1.437538881e-15}} .* fill({U.K^(3
             - i) for i in 1:size(Data.alow, 1)}, size(T_lim_c, 1) - 1),
        B_c={Data.blow,Data.bhigh,{4.938707040e6,-1.672099740e3}} .* fill({U.K,
            1}, size(T_lim_c, 1) - 1) - b_c[:, 2:3]*log(U.K),
        d=(310 + 145.2)*U.pico*U.m/U.q,
        T_lim_eta_theta={200.0,1000.0,5000.0,15000.0}*U.K,
        b_eta={fromNASAViscosity({0.62526577,-31.779652,-1.6407983e3,1.7454992}),
            fromNASAViscosity({0.87395209,561.52222,-1.7394809e5,-0.39335958}),
            fromNASAViscosity({0.88503551,909.02171,-7.3129061e5,-0.53503838})},

        b_theta={fromNASAThermalConductivity({0.85439436,105.73224,-1.2347848e4,
            0.47793128}),fromNASAThermalConductivity({0.88407146,133.57293,-1.1429640e4,
            0.24417019}),fromNASAThermalConductivity({2.4176185,8.0477749e3,
            3.1055802e6,-14.517761})});

      annotation (Documentation(info="<html>
                  <p>Notes:</p>
     <ul>
  <li>According to [<a href=\"modelica://FCSys.UsersGuide.References.Avogadro\">Avogadro</a>], the (center-to-center)
   bond length of N-N is 145.2 pm.  The radius of N is from
   <a href=\"http://en.wikipedia.org/wiki/Nitrogen\">http://en.wikipedia.org/wiki/Nitrogen</a>.  See also
   <a href=\"http://en.wikipedia.org/wiki/Van_der_Waals_radius\">http://en.wikipedia.org/wiki/Van_der_Waals_radius</a>.</li>
  <li>The virial coefficients are from [<a href=\"modelica://FCSys.UsersGuide.References.Dymond2002\">Dymond2002</a>, p.&nbsp;69].  The
  temperature range of the coefficients is [75, 745] K, but this is not enforced in the functions.  More precise virial coefficients are available from
  <a href=\"http://www.tpub.com/content/nasa1996/NASA-96-cr4755/NASA-96-cr47550059.htm\">http://www.tpub.com/content/nasa1996/NASA-96-cr4755/NASA-96-cr47550059.htm</a>.</li>
     </ul>

<p>For more information, please see the
  <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic\">Characteristic</a> package.</p></html>"));

    end Gas;

  end N2;

  package O2 "<html>O<sub>2</sub></html>"
    extends Modelica.Icons.Package;
    package Gas "O2 gas"
      import Data = Modelica.Media.IdealGases.Common.SingleGasesData.O2;

      extends BaseClasses.CharacteristicNASA(
        final formula=Data.name,
        final phase=Phase.gas,
        final m=Data.MM*U.kg/U.mol,
        n_v={-1,-4},
        b_v={{0,0,0,0,1},{5.0855e9*U.K^4,-1.6393e8*U.K^3,5.2007e5*U.K^2,-1.7696e4
            *U.K,42.859}*U.cc/U.mol},
        Deltah0_f=Data.MM*Data.Hf*U.J/U.mol,
        Deltah0=Data.MM*Data.H0*U.J/U.mol,
        T_lim_c={200.000,Data.Tlimit,6000.000,20000.000}*U.K,
        b_c={Data.alow,Data.ahigh,{4.975294300e8,-2.866106874e5,6.690352250e1,-6.169959020e-3,
            3.016396027e-7,-7.421416600e-12,7.278175770e-17}} .* fill({U.K^(3
             - i) for i in 1:size(Data.alow, 1)}, size(T_lim_c, 1) - 1),
        B_c={Data.blow,Data.bhigh,{2.293554027e6,-5.530621610e2}} .* fill({U.K,
            1}, size(T_lim_c, 1) - 1) - b_c[:, 2:3]*log(U.K),
        d=(304 + 128.2)*U.pico*U.m/U.q,
        T_lim_eta_theta={200.0,1000.0,5000.0,15000.0}*U.K,
        b_eta={fromNASAViscosity({0.60916180,-52.244847,-599.74009,2.0410801}),
            fromNASAViscosity({0.72216486,175.50839,-5.7974816e4,1.0901044}),
            fromNASAViscosity({0.73981127,391.94906,-3.7833168e5,0.90931780})},

        b_theta={fromNASAThermalConductivity({0.77229167,6.8463210,-5.8933377e3,
            1.2210365}),fromNASAThermalConductivity({0.90917351,291.24182,-7.9650171e4,
            0.064851631}),fromNASAThermalConductivity({-1.1218262,-1.9286378e4,
            2.3295011e7,20.342043})});

      annotation (Documentation(info="<html><p>Notes:</p><ul>
  <li>According to [<a href=\"modelica://FCSys.UsersGuide.References.Avogadro\">Avogadro</a>], the (center-to-center)
   bond length of O-O is 128.2 pm.  The radius of O is from
   <a href=\"http://en.wikipedia.org/wiki/Oxygen\">http://en.wikipedia.org/wiki/Oxygen</a>.  See also
   <a href=\"http://en.wikipedia.org/wiki/Van_der_Waals_radius\">http://en.wikipedia.org/wiki/Van_der_Waals_radius</a>.</li>
  <li>The virial coefficients are from [<a href=\"modelica://FCSys.UsersGuide.References.Dymond2002\">Dymond2002</a>, p.&nbsp;69].  The
  temperature range of the coefficients is [70, 495] K, but this is not enforced in the functions.</li>
     </ul>

<p>For more information, please see the
  <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic\">Characteristic</a> package.</p></html>"));

    end Gas;

  end O2;

  package IdealGas "Ideal gas"
    extends BaseClasses.CharacteristicEOS(final n_v={-1,0},final b_v=[1]);
    annotation (defaultComponentPrefixes="replaceable", Documentation(info="<html>
    <p>This package contains the pressure-volume-temperature relationship of an ideal gas.  Thermal
    data (e.g., heat capacity) is not included.</p></html>"));

  end IdealGas;

  package BaseClasses "Base classes (generally not for direct use)"
    extends Modelica.Icons.BasesPackage;
    package CharacteristicNASA
      "Thermodynamic package with diffusive properties based on NASA CEA"
      extends Characteristic;
      constant Q.TemperatureAbsolute T_lim_eta_theta[:]={0,Modelica.Constants.inf}
        "<html>Temperature limits for the rows of <i>b</i><sub>&eta;</sub> and <i>b</i><sub>&theta;</sub> (<i>T</i><sub>lim &eta; &theta;</sub>)</html>";
      constant Real b_eta[size(T_lim_eta_theta, 1) - 1, 4]
        "<html>Correlation constants for fluidity (<i>b</i><sub>&eta;</sub>)</html>";
      constant Real b_theta[size(T_lim_eta_theta, 1) - 1, 4]
        "<html>Correlation constants for thermal resistivity (<i>b</i><sub>&theta;</sub>)</html>";

    protected
      function fromNASAViscosity
        "Return constants for fluidity given NASA CEA constants for viscosity"
        extends Modelica.Icons.Function;

        input Real b[4] "NASA CEA constants for viscosity"
          annotation (Dialog(__Dymola_label="<html><i>b</i></html>"));
        output Real b_eta[4] "Constants for fluidity" annotation (Dialog(
              __Dymola_label="<html><i>b</i><sub>&eta;<sub></html>"));

      algorithm
        b_eta := {-b[1],-b[2]*U.K,-b[3]*U.K^2,-b[4] + b[1]*log(U.K) + log(1e4*U.m
          *U.s/U.g)};
        annotation (Inline=true);
      end fromNASAViscosity;

      function fromNASAThermalConductivity
        "Return constants for thermal resistivity given NASA CEA constants for thermal conductivity"
        extends Modelica.Icons.Function;

        input Real b[4] "NASA CEA constants for thermal conductivity"
          annotation (Dialog(__Dymola_label="<html><i>b</i></html>"));
        output Real b_theta[4] "Constants for thermal resistivity" annotation (
            Dialog(__Dymola_label="<html><i>b</i><sub>&theta;<sub></html>"));

      algorithm
        b_theta := {-b[1],-b[2]*U.K,-b[3]*U.K^2,-b[4] + b[1]*log(U.K) + log(1e4
          *U.m*U.K/U.W)};
        annotation (Inline=true);
      end fromNASAThermalConductivity;

    public
      redeclare function eta
        "<html>Fluidity (&eta;) as a function of temperature</html>"

        extends Modelica.Icons.Function;
        input Q.TemperatureAbsolute T=298.15*U.K "Temperature"
          annotation (Dialog(__Dymola_label="<html><i>T</i></html>"));
        input Q.VolumeSpecific v=298.15*U.K/U.atm "Specific volume";
        // Note:  Specific volume isn't used here but is included for generality.
        output Q.Resistivity eta "Fluidity"
          annotation (Dialog(__Dymola_label="<html>&eta;</html>"));

      algorithm
        /*
    assert(T_lim_eta_theta[1] <= T and T <= T_lim_eta_theta[end], "Temperature "
     + String(T/(U.K)) + " K is out of range for the resistivities ([" +
    String(T_lim_eta_theta[1]/U.K) + ", " + String(T_lim_eta_theta[end]
    /U.K) + "] K).");
  */
        // Note:  This is commented out so that the function can be inlined.
        eta := smooth(0, exp(sum(if (T_lim_eta_theta[i] <= T or i == 1) and (T
           < T_lim_eta_theta[i + 1] or i == size(T_lim_eta_theta, 1) - 1) then
          b_eta[i, 1]*log(T) + (b_eta[i, 2] + b_eta[i, 3]/T)/T + b_eta[i, 4]
           else 0 for i in 1:size(T_lim_eta_theta, 1) - 1)));
        annotation (
          InlineNoEvent=true,
          Inline=true,
          smoothOrder=0,
          Documentation(info="<html><p>This function is based on based on NASA CEA
  [<a href=\"modelica://FCSys.UsersGuide.References.McBride1996\">McBride1996</a>, <a href=\"modelica://FCSys.UsersGuide.References.Svehla1995\">Svehla1995</a>]</p>

  <p>Although specific volume is an input to this function, the result is independent of
  specific volume.</p>
  </html>"));
      end eta;

      redeclare function theta
        "<html>Thermal resistivity (&theta;) as a function of temperature</html>"

        extends Modelica.Icons.Function;
        input Q.TemperatureAbsolute T=298.15*U.K "Temperature"
          annotation (Dialog(__Dymola_label="<html><i>T</i></html>"));
        input Q.VolumeSpecific v=298.15*U.K/U.atm "Specific volume";
        // Note:  Specific volume isn't used here but is included for generality.
        output Q.ResistivityThermal theta "Thermal resistivity"
          annotation (Dialog(__Dymola_label="<html>&theta;</html>"));

      algorithm
        /*
    assert(T_lim_eta_theta[1] <= T and T <= T_lim_eta_theta[end], "Temperature "
     + String(T/(U.K)) + " K is out of range for the resistivities ([" +
    String(T_lim_eta_theta[1]/U.K) + ", " + String(T_lim_eta_theta[end]
    /U.K) + "] K).");
  */
        // Note:  This is commented out so that the function can be inlined.
        theta := smooth(0, exp(sum(if (T_lim_eta_theta[i] <= T or i == 1) and (
          T < T_lim_eta_theta[i + 1] or i == size(T_lim_eta_theta, 1) - 1)
           then b_theta[i, 1]*log(T) + (b_theta[i, 2] + b_theta[i, 3]/T)/T +
          b_theta[i, 4] else 0 for i in 1:size(T_lim_eta_theta, 1) - 1)));
        annotation (
          InlineNoEvent=true,
          Inline=true,
          smoothOrder=0,
          Documentation(info="<html><p>This function is based on based on NASA CEA
  [<a href=\"modelica://FCSys.UsersGuide.References.McBride1996\">McBride1996</a>, <a href=\"modelica://FCSys.UsersGuide.References.Svehla1995\">Svehla1995</a>]</p>

  <p>Although specific volume is an input to this function, the result is independent of
  specific volume.</p>
  </html>"));
      end theta;
      annotation (defaultComponentPrefixes="replaceable",Documentation(info="<html><p>The correlations for transport properties are available in
  [<a href=\"modelica://FCSys.UsersGuide.References.McBride1996\">McBride1996</a>,
  <a href=\"modelica://FCSys.UsersGuide.References.McBride2002\">McBride2002</a>]. For more information, please see the
  <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic\">Characteristic</a> package.</p></html>"));

    end CharacteristicNASA;

    package Characteristic "Package of thermodynamic and diffusive properties"
      import FCSys.Utilities.Chemistry.charge;
      extends CharacteristicEOS;

      constant String formula "Chemical formula";
      constant Phase phase "Material phase";
      constant Q.MassSpecific m "Specific mass";
      constant Q.LengthSpecific d "Specific diameter" annotation (Dialog);
      constant Integer z=charge(formula) "Charge number";
      constant ReferenceEnthalpy referenceEnthalpy=ReferenceEnthalpy.enthalpyOfFormationAt25degC
        "Choice of enthalpy reference";
      constant Q.PotentialChemical Deltah0_f
        "<html>Enthalpy of formation at 298.15 K, <i>p</i><sup>o</sup> (&Delta;<i>h</i><sup>o</sup><sub>f</sub>)</html>";
      constant Q.PotentialChemical Deltah0
        "<html><i>h</i><sup>o</sup>(298.15 K) - <i>h</i><sup>o</sup>(0 K) (&Delta;<i>h</i><sup>o</sup>)</html>";
      constant Q.PotentialChemical h_offset=0
        "<html>Additional enthalpy offset (<i>h</i><sub>offset</sub>)</html>";
      constant Integer n_c=-2
        "<html>Power of <i>T</i> for 1<sup>st</sup> column of <i>b</i><sub><i>c</i></sub> (<i>n</i><sub><i>c</i></sub>)</html>";
      constant Q.TemperatureAbsolute T_lim_c[:]={0,Modelica.Constants.inf}
        "<html>Temperature limits for the rows of <i>b</i><sub><i>c</i></sub> and <i>B</i><sub><i>c</i></sub> (<i>T</i><sub>lim <i>c</i></sub>)</html>";
      constant Real b_c[size(T_lim_c, 1) - 1, :]
        "<html>Coefficients of isobaric specific heat capacity at <i>p</i><sup>o</sup> as a polynomial in <i>T</i> (<i>b</i><sub><i>c</i></sub>)</html>";
      constant Real B_c[size(T_lim_c, 1) - 1, 2]
        "<html>Integration constants for specific enthalpy and entropy (<i>B</i><sub><i>c</i></sub>)</html>";

    protected
      constant Q.AreaSpecific alpha=1.5*U.pi^1.5*d^2*U.q
        "Scaled specific intercept area";
      // The intercept area is U.pi*d^2*U.q, and the additional factor is
      // 3*sqrt(U.pi)/2.

      function omega
        "<html>Root mean square of thermal velocity in one dimension as a function of temperature (&omega; = &radic;<span style=\"text-decoration:overline;\">&nbsp;<i>T</i>/<i>m</i>&nbsp;</span>)</html>"
        extends Modelica.Icons.Function;

        input Q.TemperatureAbsolute T=298.15*U.K "Temperature"
          annotation (Dialog(__Dymola_label="<html><i>T</i></html>"));
        output Q.Velocity omega
          "Root mean square of thermal velocity in one dimension"
          annotation (Dialog(__Dymola_label="<html>&omega;</html>"));

      algorithm
        omega := sqrt(T/m);
        annotation (Inline=true,Documentation(info="<html>
  <p>This function outputs the root mean square of the thermal
  velocity in any one dimension, assuming the speeds of the particles follow the
  Maxwell-Boltzmann distribution.  In combination with &alpha;, this refactorization is
  convenient for calculating
  <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.tauprime\">&tau;&prime;</a>,
  <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.mu\">&mu;</a>,
  <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.nu\">&nu;</a>,
  <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.zeta\">&zeta;</a>,
  <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.eta\">&eta;</a>, and
  <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.theta\">&theta;</a>.</p>
  </html>"));
      end omega;

    public
      function c_p
        "<html>Isobaric specific heat capacity (<i>c</i><sub><i>p</i></sub>) as a function of temperature and pressure</html>"
        import FCSys.Utilities.Polynomial;
        extends Modelica.Icons.Function;
        input Q.TemperatureAbsolute T=298.15*U.K "Temperature"
          annotation (Dialog(__Dymola_label="<html><i>T</i></html>"));
        input Q.PressureAbsolute p=p0 "Pressure"
          annotation (Dialog(__Dymola_label="<html><i>p</i></html>"));
        output Q.CapacityThermalSpecific c_p "Isobaric specific heat capacity"
          annotation (Dialog(__Dymola_label="<html><i>c<sub>p</sub></i></html>"));

      protected
        function c0_p
          "Isobaric specific heat capacity at reference pressure as a function of temperature"
          import FCSys.Utilities.Polynomial;
          input Q.TemperatureAbsolute T "Temperature";
          output Q.CapacityThermalSpecific c0_p
            "Isobaric specific heat capacity at reference pressure";

        algorithm
          /*
  assert(T_lim_c[1] <= T and T <= T_lim_c[end], "Temperature " + String(T/U.K) + " K is out of bounds for "
     + formula + " ([" + String(T_lim_c[1]/U.K) + ", " + String(T_lim_c[size(
    T_lim_c, 1)]/U.K) + "] K).");
  */
          // Note:  This is commented out so that the function can be inlined.
          c0_p := smooth(0, sum(if (T_lim_c[i] <= T or i == 1) and (T < T_lim_c[
            i + 1] or i == size(T_lim_c, 1) - 1) then Polynomial.f(
                    T,
                    b_c[i, :],
                    n_c) else 0 for i in 1:size(T_lim_c, 1) - 1));
          annotation (
            InlineNoEvent=true,
            Inline=true,
            smoothOrder=0);
        end c0_p;

        function c_p_resid
          "Residual isobaric specific heat capacity for pressure adjustment"
          import FCSys.Utilities.Polynomial;
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
                      b_v[i, :] .* {(n_v[2] - n_v[1] + j - i)*(n_v[1] - n_v[2]
                 + i - j + 1) for j in 1:size(b_v, 2)},
                      n_v[2] - n_v[1] - i) for i in rowLimits[1]:rowLimits[2]},
                    n_v[1]);
          // See s_resid() in Characteristic.s for the integral of (dels/delp)_T*dp.
          // This is temperature times the isobaric partial derivative of that
          // function with respect to temperature.  It is zero for an ideal gas.

          annotation (Inline=true);
        end c_p_resid;

      algorithm
        c_p := c0_p(T) + c_p_resid(T, p) - (if phase <> Phase.gas then
          c_p_resid(T, p0) else c_p_resid(
                T,
                p0,
                {1,-n_v[1]}));
        // See the notes in the algorithm of Characteristic.s.
        // Note:  [Dymond2002, p.17, eqs. 1.45 & 1.46] may be incorrect.
        annotation (Inline=true,Documentation(info="<html>
  <p>For an ideal gas, this function is independent of pressure
  (although pressure remains as a valid input).</p>
  </html>"));
      end c_p;

      function c_v
        "<html>Isochoric specific heat capacity (<i>c</i><sub><i>v</i></sub>) as a function of temperature and pressure</html>"
        extends Modelica.Icons.Function;
        input Q.TemperatureAbsolute T=298.15*U.K "Temperature"
          annotation (Dialog(__Dymola_label="<html><i>T</i></html>"));
        input Q.PressureAbsolute p=p0 "Pressure"
          annotation (Dialog(__Dymola_label="<html><i>p</i></html>"));
        output Q.CapacityThermalSpecific c_v "Isochoric specific heat capacity"
          annotation (Dialog(__Dymola_label="<html><i>c<sub>v</sub></i></html>"));

      algorithm
        c_v := c_p(T, p) - T*dp_Tv(
                T,
                v_Tp(T, p),
                dT=1,
                dv=0)*dv_Tp(
                T,
                p,
                dT=1,
                dp=0) "[Moran2004, p.&nbsp;546, eq. 11.66]";
        // Note 1:  This reduces to c_v = c_p - 1 for an ideal gas (where in
        // FCSys 1 = U.R).
        // Note 2:  [Dymond2002, p.17, eqs. 1.43 & 1.44] may be incorrect.
        annotation (Inline=true,Documentation(info="<html>
  <p>For an ideal gas, this function is independent of pressure
  (although pressure remains as a valid input).</p>
  </html>"));
      end c_v;

      replaceable function D
        "Diffusivity as a function of temperature and specific volume"

        extends Modelica.Icons.Function;
        input Q.TemperatureAbsolute T=298.15*U.K "Temperature"
          annotation (Dialog(__Dymola_label="<html><i>T</i></html>"));
        input Q.VolumeSpecific v=298.15*U.K/p0 "Specific volume"
          annotation (Dialog(__Dymola_label="<html><i>v</i></html>"));
        output Q.Diffusivity D "Diffusivity"
          annotation (Dialog(__Dymola_label="<html><i>D</i></html>"));

      algorithm
        D := v*omega(T)/alpha;
        annotation (Inline=true,Documentation(info="<html>
  <p>This function is based on the kinetic theory of gases under the following assumptions
  [<a href=\"modelica://FCSys.UsersGuide.References.Present1958\">Present1958</a>]:</p>
  <ol>
    <li>The particles are smooth and rigid but elastic spheres with identical radii.  This is the
    \"billiard-ball\"
    assumption, and it implies that the collisions are instantaneous and conserve kinetic
    energy.</li>
    <li>Between collisions, particles have no influence on one another.</li>
    <li>The mean free path, or average distance a particle travels between collisions, is much larger than the
    diameter of a particle.</li>
    <li>The properties carried by a particle depend only on those of the last particle with which it collided.</li>
    <li>The speeds of the particles follow the Maxwell-Boltzmann distribution.</li>
  </ol>
</html>"));
      end D;

      function g "Gibbs potential as a function of temperature and pressure"
        extends Modelica.Icons.Function;
        input Q.TemperatureAbsolute T=298.15*U.K "Temperature"
          annotation (Dialog(__Dymola_label="<html><i>T</i></html>"));
        input Q.PressureAbsolute p=p0 "Pressure"
          annotation (Dialog(__Dymola_label="<html><i>p</i></html>"));
        output Q.Potential g "Gibbs potential"
          annotation (Dialog(__Dymola_label="<html><i>g</i></html>"));

      algorithm
        g := h(T, p) - T*s(T, p);
        annotation (Inline=true);
      end g;

      function h "Specific enthalpy as a function of temperature and pressure"
        import FCSys.Utilities.Polynomial;
        extends Modelica.Icons.Function;
        input Q.TemperatureAbsolute T=298.15*U.K "Temperature"
          annotation (Dialog(__Dymola_label="<html><i>T</i></html>"));
        input Q.PressureAbsolute p=p0 "Pressure"
          annotation (Dialog(__Dymola_label="<html><i>p</i></html>"));
        output Q.Potential h "Specific enthalpy"
          annotation (Dialog(__Dymola_label="<html><i>h</i></html>"));

      protected
        function h0_i
          "Return h0 as a function of T using one of the temperature intervals"
          import FCSys.Utilities.Polynomial;
          input Q.TemperatureAbsolute T "Temperature";
          input Integer i "Index of the temperature interval";
          output Q.Potential h0
            "Specific enthalpy at given temperature relative to enthalpy of formation at 25 degC, both at reference pressure";

        algorithm
          h0 := Polynomial.F(
                    T,
                    b_c[i, :],
                    n_c) + B_c[i, 1] annotation (Inline=true, derivative=dh0_i);
          // This is the integral of c0_p*dT up to T at p0.  The lower bound is the
          // enthalpy of formation (of ideal gas, if the material is gaseous) at
          // 25 degC [McBride2002, p. 2].

        end h0_i;

        function dh0_i "Derivative of h0_i"
          // Note:  This function is necessary for Dymola 7.4 to differentiate h().
          import FCSys.Utilities.Polynomial;
          input Q.TemperatureAbsolute T "Temperature";
          input Integer i "Index of the temperature interval";
          input Q.Temperature dT "Derivative of temperature";
          output Q.Potential dh0
            "Derivative of specific enthalpy at reference pressure";

        algorithm
          dh0 := Polynomial.f(
                    T,
                    b_c[i, :],
                    n_c)*dT annotation (Inline=true);

        end dh0_i;

        function h_resid
          "Residual specific enthalpy for pressure adjustment for selected rows of b_v"
          import FCSys.Utilities.Polynomial;
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
                      b_v[i, :] .* {n_v[1] - n_v[2] + i - j + 1 for j in 1:size(
                b_v, 2)},
                      n_v[2] - n_v[1] - i + 1) for i in rowLimits[1]:rowLimits[
              2]},  n_v[1] + rowLimits[1] - 1) annotation (Inline=true);
          // Note:  The partial derivative (delh/delp)_T is equal to v +
          // T*(dels/delp)_T by definition of enthalpy change (dh = T*ds + v*dp)
          // and then to v - T*(delv/delT)_p by applying the appropriate Maxwell
          // relation, (dels/delp)_T = -(delv/delT)_p.
          // Note:  This is zero for an ideal gas.

        end h_resid;

      algorithm
        /*
    assert(T_lim_c[1] <= T and T <= T_lim_c[end], "Temperature " +
    String(T/U.K) + " K is out of bounds for " + name + " ([" + String(T_lim_c[1]
    /U.K) + ", " + String(T_lim_c[size(T_lim_c, 1)]/U.K) + "] K).");
  */
        // Note:  This is commented out so that the function can be inlined.
        h := smooth(1, sum(if (T_lim_c[i] <= T or i == 1) and (T < T_lim_c[i +
          1] or i == size(b_c, 1)) then h0_i(T, i) else 0 for i in 1:size(b_c,
          1))) + (if referenceEnthalpy == ReferenceEnthalpy.zeroAt0K then
          Deltah0 else 0) - (if referenceEnthalpy == ReferenceEnthalpy.enthalpyOfFormationAt25degC
           then 0 else Deltah0_f) + h_offset + h_resid(T, p) - (if phase <>
          Phase.gas then h_resid(T, p0) else h_resid(
                T,
                p0,
                {1,-n_v[1]}));
        // The last two terms adjust for the actual pressure relative to the
        // reference.  In general, the lower limit of the integral of
        // (delh/delp)_T*dp is the reference pressure (p0).  However, if the
        // material is gaseous, then the reference is the corresponding ideal gas.
        // In that case, the lower limit of the real gas terms of the integral is
        // p=0, where a real gas behaves as an ideal gas.  See [Rao1997, p. 271].
        annotation (
          InlineNoEvent=true,
          Inline=true,
          smoothOrder=1,
          Documentation(info="<html>
  <p>For an ideal gas, this function is independent of pressure
  (although pressure remains as a valid input).</p>
    </html>"));
      end h;

      function s "Specific entropy as a function of temperature and pressure"
        import FCSys.Utilities.Polynomial;
        extends Modelica.Icons.Function;
        input Q.TemperatureAbsolute T=298.15*U.K "Temperature"
          annotation (Dialog(__Dymola_label="<html><i>T</i></html>"));
        input Q.PressureAbsolute p=p0 "Pressure"
          annotation (Dialog(__Dymola_label="<html><i>p</i></html>"));
        output Q.NumberAbsolute s "Specific entropy"
          annotation (Dialog(__Dymola_label="<html><i>s</i></html>"));

      protected
        function s0_i
          "Return s0 as a function of T using one of the temperature intervals"
          import FCSys.Utilities.Polynomial.F;
          input Q.TemperatureAbsolute T "Temperature";
          input Integer i "Index of the temperature interval";
          output Q.NumberAbsolute s0
            "Specific entropy at given temperature and reference pressure";

        algorithm
          s0 := F(  T,
                    b_c[i, :],
                    n_c - 1) + B_c[i, 2];
          // This is the integral of c0_p/T*dT up to T at p0 with the absolute
          // entropy at the lower bound [McBride2002, p. 2].

          annotation (Inline=true, derivative=ds0_i);
        end s0_i;

        function ds0_i "Derivative of s0_i"
          // Note:  This function is necessary for Dymola 7.4 to differentiate s().
          import FCSys.Utilities.Polynomial.f;
          input Q.TemperatureAbsolute T "Temperature";
          input Integer i "Index of the temperature interval";
          input Q.Temperature dT "Derivative of temperature";
          output Q.Number ds0
            "Derivative of specific entropy at given temperature and reference pressure";

        algorithm
          ds0 := f( T,
                    b_c[i, :],
                    n_c - 1)*dT;
          annotation (Inline=true);
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
                      b_v[i, :] .* {n_v[1] - n_v[2] + i - j for j in 1:size(b_v,
                2)},  n_v[2] - n_v[1] - i) for i in rowLimits[1]:rowLimits[2]},
                    n_v[1] + rowLimits[1] - 1);
          // Note:  According to the Maxwell relations,
          // (dels/delp)_T = -(delv/delT)_p.

          annotation (Inline=true, derivative=ds_resid);
        end s_resid;

        function ds_resid "Derivative of s_resid"
          // Note:  This function is necessary for Dymola 7.4 to differentiate s().
          input Q.TemperatureAbsolute T "Temperature";
          input Q.PressureAbsolute p "Pressure";
          input Integer rowLimits[2]={1,size(b_v, 1)}
            "Beginning and ending indices of rows of b_v to be included";
          input Q.Temperature dT "Derivative of temperature";
          input Q.Pressure dp "Derivative of pressure";
          output Q.Number ds_resid
            "Derivative of integral of (dels/delp)_T*dp up to p with zero integration constant (for selected rows)";

        algorithm
          ds_resid := Polynomial.dF(
                    p,
                    {Polynomial.df(
                      T,
                      b_v[i, :] .* {n_v[1] - n_v[2] + i - j for j in 1:size(b_v,
                2)},  n_v[2] - n_v[1] - i,
                      dT) for i in rowLimits[1]:rowLimits[2]},
                    n_v[1] + rowLimits[1] - 1,
                    dp);
          annotation (Inline=true);
        end ds_resid;

      algorithm
        /*
  assert(T_lim_c[1] <= T and T <= T_lim_c[end], "Temperature " +
    String(T/U.K) + " K is out of bounds for " + name + " ([" + String(T_lim_c[1]
    /U.K) + ", " + String(T_lim_c[size(T_lim_c, 1)]/U.K) + "] K).");
  */
        // Note:  This is commented out so that the function can be inlined.
        s := smooth(1, sum(if (T_lim_c[i] <= T or i == 1) and (T < T_lim_c[i +
          1] or i == size(b_c, 1)) then s0_i(T, i) else 0 for i in 1:size(b_c,
          1))) + s_resid(T, p) - (if phase <> Phase.gas then s_resid(T, p0)
           else s_resid(
                T,
                p0,
                {1,-n_v[1]}));
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
        // b_v and n_v) must include an ideal gas term (v = ... + f(T)/p +
        // ...).  Otherwise, an indexing error will occur.
        annotation (
          InlineNoEvent=true,
          Inline=true,
          smoothOrder=1);
      end s;

      replaceable function zeta
        "<html>** (&zeta;) as a function of temperature</html>"

        extends Modelica.Icons.Function;
        input Q.TemperatureAbsolute T=298.15*U.K "Temperature"
          annotation (Dialog(__Dymola_label="<html><i>T</i></html>"));
        input Q.VolumeSpecific v=298.15*U.K/p0 "Specific volume"
          annotation (Dialog(__Dymola_label="<html><i>v</i></html>"));
        output Q.DiffusivityMassSpecific zeta "**"
          annotation (Dialog(__Dymola_label="<html>&zeta;</html>"));

      algorithm
        zeta := m*D(T, v);

        annotation (Inline=true,Documentation(info="<html>
  <p>**Continuity is a property is defined in <a href=\"modelica://FCSys\">FCSys</a>
  resistivity to axial compression or material storage during transport.</p>
  </html>"));
      end zeta;

      replaceable function eta
        "<html>Fluidity (&eta;) as a function of temperature</html>"

        extends Modelica.Icons.Function;
        input Q.TemperatureAbsolute T=298.15*U.K "Temperature"
          annotation (Dialog(__Dymola_label="<html><i>T</i></html>"));
        input Q.VolumeSpecific v=298.15*U.K/p0 "Specific volume"
          annotation (Dialog(__Dymola_label="<html><i>v</i></html>"));
        // Note:  Specific volume isn't used here but is included for generality.
        output Q.Fluidity eta "Fluidity"
          annotation (Dialog(__Dymola_label="<html>&eta;</html>"));

      algorithm
        eta := alpha/(m*omega(T));
        annotation (Inline=true,Documentation(info="<html>
  <p>Fluidity is defined as the reciprocal of dynamic viscosity
  (see <a href=\"http://en.wikipedia.org/wiki/Viscosity#Fluidity\">http://en.wikipedia.org/wiki/Viscosity#Fluidity</a>).</p>

  <p>Although specific volume is an input to this function, the result is independent of
  specific volume.  According to Present
  [<a href=\"modelica://FCSys.UsersGuide.References.Present1958\">Present1958</a>], this independence very accurately
  matches the measured fluidity of gases.  However, the fluidity varies by species and
  generally falls more rapidly with temperature than indicated
  [<a href=\"modelica://FCSys.UsersGuide.References.Present1958\">Present1958</a>, p.&nbsp;41].</p>

  <p>This function is based on the kinetic theory of gases under the following assumptions
  [<a href=\"modelica://FCSys.UsersGuide.References.Present1958\">Present1958</a>]:</p>
  <ol>
    <li>The particles are smooth and rigid but elastic spheres with identical radii.  This is the
    \"billiard-ball\"
    assumption, and it implies that the collisions are instantaneous and conserve kinetic
    energy.</li>
    <li>Between collisions, particles have no influence on one another.</li>
    <li>The mean free path, or average distance a particle travels between collisions, is much larger than the
    diameter of a particle.</li>
    <li>The properties carried by a particle depend only on those of the last particle with which it collided.</li>
    <li>The speeds of the particles follow the Maxwell-Boltzmann distribution.</li>
  </ol>
  </html>"));
      end eta;

      replaceable function theta
        "<html>Thermal resistivity (&theta;) as a function of temperature and specific volume</html>"

        extends Modelica.Icons.Function;
        input Q.TemperatureAbsolute T=298.15*U.K "Temperature"
          annotation (Dialog(__Dymola_label="<html><i>T</i></html>"));
        input Q.VolumeSpecific v=298.15*U.K/p0 "Specific volume"
          annotation (Dialog(__Dymola_label="<html><i>v</i></html>"));
        output Q.ResistivityThermal theta "Thermal resistivity"
          annotation (Dialog(__Dymola_label="<html>&theta;</html>"));

      algorithm
        theta := alpha/(c_v(T, p_Tv(T, v))*omega(T));
        annotation (Inline=true,Documentation(info="<html>
  <p>This function is based on the kinetic theory of gases under the following assumptions
  [<a href=\"modelica://FCSys.UsersGuide.References.Present1958\">Present1958</a>]:</p>
  <ol>
    <li>The particles are smooth and rigid but elastic spheres with identical radii.  This is the
    \"billiard-ball\"
    assumption, and it implies that the collisions are instantaneous and conserve kinetic
    energy.</li>
    <li>Between collisions, particles have no influence on one another.</li>
    <li>The mean free path, or average distance a particle travels between collisions, is much larger than the
    diameter of a particle.</li>
    <li>The properties carried by a particle depend only on those of the last particle with which it collided.</li>
    <li>The speeds of the particles follow the Maxwell-Boltzmann distribution.</li>
  </ol>
</html>"));
      end theta;

      replaceable function tauprime
        "<html>Phase change interval (&tau;&prime;) as a function of temperature and specific volume</html>"

        extends Modelica.Icons.Function;
        input Q.TemperatureAbsolute T=298.15*U.K "Temperature"
          annotation (Dialog(__Dymola_label="<html><i>T</i></html>"));
        input Q.VolumeSpecific v=298.15*U.K/p0 "Specific volume"
          annotation (Dialog(__Dymola_label="<html><i>v</i></html>"));
        output Q.TimeAbsolute tauprime "Phase change interval"
          annotation (Dialog(__Dymola_label="<html>&tau;&prime;</html>"));

      algorithm
        tauprime := v/(alpha*omega(T));
        annotation (Inline=true,Documentation(info="<html>
  <p>This function is based on the kinetic theory of gases under the following assumptions
  [<a href=\"modelica://FCSys.UsersGuide.References.Present1958\">Present1958</a>]:</p>
  <ol>
    <li>The particles are smooth and rigid but elastic spheres with identical radii.  This is the
    \"billiard-ball\"
    assumption, and it implies that the collisions are instantaneous and conserve kinetic
    energy.</li>
    <li>Between collisions, particles have no influence on one another.</li>
    <li>The mean free path, or average distance a particle travels between collisions, is much larger than the
    diameter of a particle.</li>
    <li>The properties carried by a particle depend only on those of the last particle with which it collided.</li>
    <li>The speeds of the particles follow the Maxwell-Boltzmann distribution.</li>
  </ol>
  <p>Also, it is assumed that the Einstein relation applies.</p>

  <p>Please see [<a href=\"modelica://FCSys.UsersGuide.References\">Davies2013</a>, Ch. 3] for a derivation
  of the rate of phase change from kinetic theory.</p>

  <p>Although specific volume is an input to this function, the result is independent of
  specific volume.</p>
</html>"));
      end tauprime;

      replaceable function mu
        "<html>Mobility (&mu;) as a function of temperature and specific volume</html>"

        extends Modelica.Icons.Function;
        input Q.TemperatureAbsolute T=298.15*U.K "Temperature"
          annotation (Dialog(__Dymola_label="<html><i>T</i></html>"));
        input Q.VolumeSpecific v=298.15*U.K/p0 "Specific volume"
          annotation (Dialog(__Dymola_label="<html><i>v</i></html>"));
        output Q.Mobility mu "Mobility"
          annotation (Dialog(__Dymola_label="<html>&mu;</html>"));

      algorithm
        mu := v/(m*alpha*omega(T));
        annotation (Inline=true,Documentation(info="<html>
  <p>This function is based on the kinetic theory of gases under the following assumptions
  [<a href=\"modelica://FCSys.UsersGuide.References.Present1958\">Present1958</a>]:</p>
  <ol>
    <li>The particles are smooth and rigid but elastic spheres with identical radii.  This is the
    \"billiard-ball\"
    assumption, and it implies that the collisions are instantaneous and conserve kinetic
    energy.</li>
    <li>Between collisions, particles have no influence on one another.</li>
    <li>The mean free path, or average distance a particle travels between collisions, is much larger than the
    diameter of a particle.</li>
    <li>The properties carried by a particle depend only on those of the last particle with which it collided.</li>
    <li>The speeds of the particles follow the Maxwell-Boltzmann distribution.</li>
  </ol>
  <p>Also, it is assumed that the Einstein relation applies.</p>
</html>"));
      end mu;

      replaceable function nu
        "<html>Thermal independity (&nu;) as a function of temperature and specific volume</html>"

        extends Modelica.Icons.Function;
        input Q.TemperatureAbsolute T=298.15*U.K "Temperature"
          annotation (Dialog(__Dymola_label="<html><i>T</i></html>"));
        input Q.VolumeSpecific v=298.15*U.K/p0 "Specific volume"
          annotation (Dialog(__Dymola_label="<html><i>v</i></html>"));
        output Q.TimeAbsolute nu "Thermal independity"
          annotation (Dialog(__Dymola_label="<html>&nu;</html>"));

      algorithm
        nu := v/(c_p(T, p_Tv(T, v))*alpha*omega(T));
        annotation (Inline=true,Documentation(info="<html>
<p><i>Thermal independity</i> describes the extent to which an exchange of thermal energy between species causes or requires a
temperature difference.</p>

  <p>This function is based on the kinetic theory of gases under the following assumptions
  [<a href=\"modelica://FCSys.UsersGuide.References.Present1958\">Present1958</a>]:</p>
  <ol>
    <li>The particles are smooth and rigid but elastic spheres with identical radii.  This is the
    \"billiard-ball\"
    assumption, and it implies that the collisions are instantaneous and conserve kinetic
    energy.</li>
    <li>Between collisions, particles have no influence on one another.</li>
    <li>The mean free path, or average distance a particle travels between collisions, is much larger than the
    diameter of a particle.</li>
    <li>The properties carried by a particle depend only on those of the last particle with which it collided.</li>
    <li>The speeds of the particles follow the Maxwell-Boltzmann distribution.</li>
  </ol>
  <p>Also, it is assumed that the Einstein relation applies.</p>
</html>"));
      end nu;

      annotation (defaultComponentPrefixes="replaceable",Documentation(info="<html>
    <p>This package is compatible with NASA CEA thermodynamic data
    [<a href=\"modelica://FCSys.UsersGuide.References.McBride2002\">McBride2002</a>] and the virial equation of state
    [<a href=\"modelica://FCSys.UsersGuide.References.Dymond2002\">Dymond2002</a>].</p>

<p>Notes regarding the constants:</p>
    <ul>
    <li>Currently, <code>formula</code> may not contain parentheses or brackets.</li>

    <li><i>d</i> is the Van der Waals diameter or the diameter for the
    rigid-sphere (\"billiard-ball\") approximation of the kinetic theory of gases
    [<a href=\"modelica://FCSys.UsersGuide.References.Present1958\">Present1958</a>].</li>

    <li><i>b<sub>c</sub></i>: The rows give the coefficients for the temperature intervals bounded
    by the values in <i>T</i><sub>lim <i>c</i></sub>.
    The powers of <i>T</i> increase
    by column.
    By default,
    the powers of <i>T</i> for the first column are each -2, which corresponds to [<a href=\"modelica://FCSys.UsersGuide.References.McBride2002\">McBride2002</a>].
    In that case, the dimensionalities of the coefficients are {L4.M2/(N2.T4), L2.M/(N.T2), 1, &hellip;}
    for each row, where L is length, M is mass, N is particle number, and T is time. (In <a href=\"modelica://FCSys\">FCSys</a>,
    temperature is a potential with dimension L2.M/(N.T2); see
    the <a href=\"modelica://FCSys.Units\">Units</a> package.)</li>

    <li><i>B<sub>c</sub></i>: As in <i>b<sub>c</sub></i>, the rows correspond to different
    temperature intervals.  The first column is for specific enthalpy and has dimensionality
    L2.M/(N.T2).  The second is for specific entropy and is dimensionless.
    The integration constants for enthalpy are defined such that the enthalpy at
    25&nbsp;&deg;C is the specific enthalpy of formation at that temperature and reference pressure
    [<a href=\"modelica://FCSys.UsersGuide.References.McBride2002\">McBride2002</a>, p.&nbsp;2].
    The integration constants for specific entropy are defined such that specific entropy is absolute.</li>

    <li><i>T</i><sub>lim <i>c</i></sub>: The first and last entries are the minimum and
    maximum valid temperatures.  The intermediate entries are the thresholds
    between rows of <i>b<sub>c</sub></i> (and <i>B<sub>c</sub></i>).  Therefore, if there are <i>n</i> temperature intervals
    (and rows in <i>b<sub>c</sub></i> and <i>B<sub>c</sub></i>), then <i>T</i><sub>lim <i>c</i></sub> must
    have <i>n</i> + 1 entries.</li>

    <li>The reference pressure is <i>p</i><sup>o</sup>.   In the
    NASA CEA data [<a href=\"modelica://FCSys.UsersGuide.References.McBride2002\">McBride2002</a>], it is 1&nbsp;bar for gases and 1&nbsp;atm for condensed
    species.  For gases, the reference state is the ideal gas at <i>p</i><sup>o</sup>.
    For example, the enthalpy of a non-ideal (real) gas at 25&nbsp;&deg;C and <i>p</i><sup>o</sup> with
    <code>ReferenceEnthalpy.zeroAt25degC</code> is not exactly zero.</li>

    <li>If the material is gaseous (<code>phase == Phase.gas</code>), then the first virial coefficient
    must be independent of temperature.  Otherwise, the function for specific enthalpy
    (<a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.h\"><i>h</i></a>) will be ill-posed.
    Typically, the first virial coefficient is one (or equivalently <code>U.R</code>), which satisfies
    this requirement.</li>
    </ul></html>"));

    end Characteristic;

    package CharacteristicEOS
      "<html>Base thermodynamic package with only the <i>p</i>-<i>v</i>-<i>T</i> relations</html>"
      import Modelica.Math.BooleanVectors.anyTrue;
      extends Modelica.Icons.MaterialPropertiesPackage;

      constant Q.PressureAbsolute p0=U.bar
        "<html>Reference pressure (<i>p</i><sup>o</sup>)</html>";
      constant Integer n_v[2]={-1,0}
        "<html>Powers of <i>p</i>/<i>T</i> and <i>T</i> for 1<sup>st</sup> row and column of <i>b</i><sub><i>v</i></sub> (<i>n</i><sub><i>v</i></sub>)</html>";
      constant Real b_v[:, :]=[1]
        "<html>Coefficients for specific volume as a polynomial in <i>p</i>/<i>T</i> and <i>T</i> (<i>b</i><sub><i>v</i></sub>)</html>";
      // Note:  p/T is the argument instead of p so that b_p will have the same
      // size as b_v for the typical definitions of the second virial
      // coefficients in [Dymond2002].
      final constant Boolean isCompressible=anyTrue({anyTrue({abs(b_v[i, j]) >
          Modelica.Constants.small and n_v[1] + i - 1 <> 0 for i in 1:size(b_v,
          1)}) for j in 1:size(b_v, 2)})
        "<html><code>true</code>, if density depends on pressure</html>";
      final constant Boolean hasThermalExpansion=anyTrue({anyTrue({abs(b_v[i, j])
           > Modelica.Constants.small and n_v[2] + j - n_v[1] - i <> 0 for i
           in 1:size(b_v, 1)}) for j in 1:size(b_v, 2)})
        "<html><code>true</code>, if density depends on temperature</html>";

    protected
      final constant Integer n_p[2]={n_v[1] - size(b_v, 1) + 1,n_v[2] + 1}
        "<html>Powers of <i>v</i> and <i>T</i> for 1<sup>st</sup> row and column of <i>b<sub>p</sub></i></html>";
      final constant Real b_p[size(b_v, 1), size(b_v, 2)]=if size(b_v, 1) == 1
           then b_v .^ (-n_p[1]) else {(if n_v[1] + i == 0 or n_v[1] + i == 1
           or size(b_v, 1) == 1 then b_v[i, :] else (if n_v[1] + i == 2 and n_v[
          1] <= 0 then b_v[i, :] + b_v[i - 1, :] .^ 2 else (if n_v[1] + i == 3
           and n_v[1] <= 0 then b_v[i, :] + b_v[i - 2, :] .* (b_v[i - 2, :] .^
          2 + 3*b_v[i - 1, :]) else zeros(size(b_v, 2))))) for i in size(b_v, 1)
          :-1:1}
        "<html>Coefficients of <i>p</i> as a polynomial in <i>v</i> and <i>T</i></html>";
      // Note:  This is from [Dymond2002, p.&nbsp;2].  If necessary, additional terms
      // can be computed using FCSys/Resources/virial-relations.cdf.

    public
      function dp_Tv
        "<html>Derivative of pressure as defined by <a href=\"modelica://FCSys.Characteristics.BaseClasses.CharacteristicEOS.p_Tv\"><i>p<sub>T v</sub></i></a>()</html>"
        import FCSys.Utilities.Polynomial;
        extends Modelica.Icons.Function;

        input Q.TemperatureAbsolute T=298.15*U.K "Temperature"
          annotation (Dialog(__Dymola_label="<html><i>T</i></html>"));
        input Q.VolumeSpecificAbsolute v=298.15*U.K/p0 "Specific volume"
          annotation (Dialog(__Dymola_label="<html><i>v</i></html>"));
        input Q.Temperature dT=0 "Derivative of temperature"
          annotation (Dialog(__Dymola_label="<html>d<i>T</i></html>"));
        input Q.VolumeSpecific dv=0 "Derivative of specific volume"
          annotation (Dialog(__Dymola_label="<html>d<i>v</i></html>"));
        output Q.Pressure dp "Derivative of pressure"
          annotation (Dialog(__Dymola_label="<html>d<i>p</i></html>"));

      algorithm
        dp := if isCompressible then Polynomial.f(
                v,
                {Polynomial.f(
                  T,
                  b_p[i, :] .* {(n_p[1] + i - 1)*T*dv + (n_p[2] + j - 1)*v*dT
              for j in 1:size(b_p, 2)},
                  n_p[2] - 1) for i in 1:size(b_p, 1)},
                n_p[1] - 1) else 0;

        annotation (
          Inline=true,
          inverse(dv=dv_Tp(
                      T,
                      p_Tv(T, v),
                      dT,
                      dp)),
          smoothOrder=999);
      end dp_Tv;

      function dv_Tp
        "<html>Derivative of specific volume as defined by <a href=\"modelica://FCSys.Characteristics.BaseClasses.CharacteristicEOS.v_Tp\"><i>v<sub>T p</sub></i></a>()</html>"
        import FCSys.Utilities.Polynomial;
        extends Modelica.Icons.Function;
        input Q.TemperatureAbsolute T=298.15*U.K "Temperature"
          annotation (Dialog(__Dymola_label="<html><i>T</i></html>"));
        input Q.PressureAbsolute p=p0 "Pressure"
          annotation (Dialog(__Dymola_label="<html><i>p</i></html>"));
        input Q.Temperature dT=0 "Derivative of temperature"
          annotation (Dialog(__Dymola_label="<html>d<i>T</i></html>"));
        input Q.Pressure dp=0 "Derivative of pressure"
          annotation (Dialog(__Dymola_label="<html>d<i>p</i></html>"));
        output Q.VolumeSpecific dv "Derivative of specific volume"
          annotation (Dialog(__Dymola_label="<html>d<i>v</i></html>"));

      algorithm
        dv := Polynomial.f(
                p,
                {Polynomial.f(
                  T,
                  b_v[i, :] .* {(n_v[1] + i - 1)*T*dp + (n_v[2] - n_v[1] + j -
              i)*p*dT for j in 1:size(b_v, 2)},
                  n_v[2] - n_v[1] - i) for i in 1:size(b_v, 1)},
                n_v[1] - 1);

        annotation (Inline=true, inverse(dp=dp_Tv(
                      T,
                      v_Tp(T, p),
                      dT,
                      dv)));
      end dv_Tp;

      function p_Tv
        "<html>Pressure as a function of temperature and specific volume (<i>p<sub>T v</sub></i>())</html>"
        import FCSys.Utilities.Polynomial;
        extends Modelica.Icons.Function;
        input Q.TemperatureAbsolute T=298.15*U.K "Temperature"
          annotation (Dialog(__Dymola_label="<html><i>T</i></html>"));
        input Q.VolumeSpecificAbsolute v=298.15*U.K/p0 "Specific volume"
          annotation (Dialog(__Dymola_label="<html><i>v</i></html>"));
        output Q.PressureAbsolute p "Pressure"
          annotation (Dialog(__Dymola_label="<html><i>p</i></html>"));

      algorithm
        // assert(isCompressible,
        //  "The pressure is undefined since the material is incompressible.",
        //  AssertionLevel.warning);
        // Note:  This isn't used because it creates an error instead of a warning
        // in Dymola 2014
        p := if isCompressible then Polynomial.f(
                v,
                {Polynomial.f(
                  T,
                  b_p[i, :],
                  n_p[2]) for i in 1:size(b_p, 1)},
                n_p[1]) else p0;
        annotation (
          Inline=true,
          inverse(v=v_Tp(T, p)),
          derivative=dp_Tv,
          smoothOrder=999,
          Documentation(info="<html><p>If the species is incompressible then <i>p</i>(<i>T</i>, <i>v</i>) is undefined,
  and the function will return a value of zero.</p>

<p>The derivative of this function is <a href=\"modelica://FCSys.Characteristics.BaseClasses.CharacteristicEOS.dp_Tv\">dp_Tv</a>().</p></html>"));
      end p_Tv;

      function v_Tp
        "<html>Specific volume as a function of temperature and pressure (<i>v<sub>T p</sub></i>())</html>"
        import FCSys.Utilities.Polynomial;
        extends Modelica.Icons.Function;
        input Q.TemperatureAbsolute T=298.15*U.K "Temperature"
          annotation (Dialog(__Dymola_label="<html><i>T</i></html>"));
        input Q.PressureAbsolute p=p0 "Pressure"
          annotation (Dialog(__Dymola_label="<html><i>p</i></html>"));
        output Q.VolumeSpecificAbsolute v "Specific volume";

      algorithm
        v := Polynomial.f(
                p,
                {Polynomial.f(
                  T,
                  b_v[i, :],
                  n_v[2] - n_v[1] - i + 1) for i in 1:size(b_v, 1)},
                n_v[1]);
        annotation (
          Inline=true,
          inverse(p=p_Tv(T, v)),
          derivative=dv_Tp,
          Documentation(info="<html>
  <p>The derivative of this function is
  <a href=\"modelica://FCSys.Characteristics.BaseClasses.CharacteristicEOS.dv_Tp\">dv_Tp</a>().</p></html>"));
      end v_Tp;

      function beta
        "<html>Isothermal compressibility as a function of temperature and pressure (&beta;)</html>"
        extends Modelica.Icons.Function;

        input Q.TemperatureAbsolute T=298.15*U.K "Temperature"
          annotation (Dialog(__Dymola_label="<html><i>T</i></html>"));
        input Q.PressureAbsolute p=p0 "Pressure"
          annotation (Dialog(__Dymola_label="<html><i>p</i></html>"));

        output Q.PressureReciprocal beta "Isothermal compressibility";

      algorithm
        beta := -dv_Tp(
                T=T,
                p=p,
                dT=0,
                dp=1)/v_Tp(T, p);
        annotation (Inline=true,Documentation(info="<html>
  <p>For an ideal gas, this function is independent of temperature
  (although temperature remains as a valid input).</p>
  </html>"));
      end beta;
      annotation (defaultComponentPrefixes="replaceable",Documentation(info="<html>
    <p>This package may be used with
    the assumption of ideal gas or of constant specific volume, although it is more general than
    that.</p>

<p>Notes regarding the constants:</p>
    <ul>
    <li><i>b<sub>v</sub></i>: The powers of <i>p</i>/<i>T</i> increase by row.  The powers of
    <i>T</i> increase by column.  If <code>n_v[1] == -1</code>, then the rows
    of <i>b<sub>v</sub></i> correspond to 1, <i>B</i><sup>*</sup><i>T</i>,
    <i>C</i><sup>*</sup><i>T</i><sup>2</sup>, <i>D</i><sup>*</sup><i>T</i><sup>3</sup>, &hellip;,
    where
    1, <i>B</i><sup>*</sup>, <i>C</i><sup>*</sup>, and <i>D</i><sup>*</sup> are
    the first, second, third, and fourth coefficients in the volume-explicit
    virial equation of state
    [<a href=\"modelica://FCSys.UsersGuide.References.Dymond2002\">Dymond2002</a>, pp.&nbsp;1&ndash;2].
    Currently,
    virial equations of state are supported up to the fourth coefficient (<i>D</i><sup>*</sup>).
    If additional terms are required, review and modify the definition of <i>b<sub>p</sub></i>.</li>

    <li>The defaults for <i>b<sub>v</sub></i> and <i>n<sub>v</sub></i> represent ideal gas.</li>
    </ul></html>"));

    end CharacteristicEOS;

    type ReferenceEnthalpy = enumeration(
        zeroAt0K "Enthalpy at 0 K and p0 is 0 (if no additional offset)",
        zeroAt25degC
          "Enthalpy at 25 degC and p0 is 0 (if no additional offset)",
        enthalpyOfFormationAt25degC
          "Enthalpy at 25 degC and p0 is enthalpy of formation at 25 degC and p0 (if no additional offset)")
      "Enumeration for the reference enthalpy of a species";
    type Phase = enumeration(
        gas "Gas",
        solid "Solid",
        liquid "Liquid") "Enumeration for material phases";

  end BaseClasses;

  package MobilityFactors "Mobility factors"
    extends Modelica.Icons.Package;
    function k_H2_H2O
      "<html>Binary mobility factor for H<sub>2</sub> and H<sub>2</sub>O</html>"
      import Modelica.Media.IdealGases.Common.FluidData;

      extends BaseClasses.k(
        final T_crit=sqrt(FluidData.H2.criticalTemperature*FluidData.H2O.criticalTemperature)
            *U.K,
        final p_crit=sqrt(FluidData.H2.criticalPressure*FluidData.H2O.criticalPressure)
            *U.Pa,
        redeclare package A = H2.Gas,
        redeclare package B = H2O.Gas,
        redeclare function pDstar = BaseClasses.pDstar_polar);

      annotation (Inline=true);
    end k_H2_H2O;

    function k_H2O_N2
      "<html>Binary mobility factor for H<sub>2</sub>O and N<sub>2</sub></html>"
      import Modelica.Media.IdealGases.Common.FluidData;

      extends BaseClasses.k(
        final T_crit=sqrt(FluidData.H2O.criticalTemperature*FluidData.N2.criticalTemperature)
            *U.K,
        final p_crit=sqrt(FluidData.H2O.criticalPressure*FluidData.N2.criticalPressure)
            *U.Pa,
        redeclare package A = H2O.Gas,
        redeclare package B = N2.Gas,
        redeclare function pDstar = BaseClasses.pDstar_polar);

      annotation (Inline=true);
    end k_H2O_N2;

    function k_H2O_O2
      "<html>Binary mobility factor for H<sub>2</sub>O and O<sub>2</sub></html>"
      import Modelica.Media.IdealGases.Common.FluidData;

      extends BaseClasses.k(
        final T_crit=sqrt(FluidData.H2O.criticalTemperature*FluidData.O2.criticalTemperature)
            *U.K,
        final p_crit=sqrt(FluidData.H2O.criticalPressure*FluidData.O2.criticalPressure)
            *U.Pa,
        redeclare package A = H2O.Gas,
        redeclare package B = O2.Gas,
        redeclare function pDstar = BaseClasses.pDstar_polar);

      annotation (Inline=true);
    end k_H2O_O2;

    function k_N2_O2
      "<html>Binary mobility factor for N<sub>2</sub> and O<sub>2</sub></html>"
      import Modelica.Media.IdealGases.Common.FluidData;

      extends BaseClasses.k(
        final T_crit=sqrt(FluidData.N2.criticalTemperature*FluidData.O2.criticalTemperature)
            *U.K,
        final p_crit=sqrt(FluidData.N2.criticalPressure*FluidData.O2.criticalPressure)
            *U.Pa,
        redeclare package A = N2.Gas,
        redeclare package B = O2.Gas);

      annotation (Inline=true);
    end k_N2_O2;

    package BaseClasses "Base classes (generally not for direct use)"
      extends Modelica.Icons.BasesPackage;
      function k "Binary mobility factor"
        import harmonicMean = FCSys.Utilities.Means.harmonic;

        input Q.TemperatureAbsolute T=298.15*U.K "Temperature"
          annotation (Dialog(__Dymola_label="<html><i>T</i></html>"));
        input Q.PressureAbsolute p_A=U.atm
          "<html>Pressure of the 1<sup>st</sup> species</html>" annotation (
            Dialog(__Dymola_label="<html><i>p</i><sub>A</sub></html>"));
        input Q.PressureAbsolute p_B=U.atm
          "<html>Pressure of the 2<sup>nd</sup> species</html>" annotation (
            Dialog(__Dymola_label="<html><i>p</i><sub>B</sub></html>"));
        input Q.VolumeSpecific v_A=298.15*U.K/A.p0
          "<html>Specific volume of the 1<sup>st</sup> species</html>"
          annotation (Dialog(__Dymola_label="<html><i>v</i><sub>A</sub></html>"));
        input Q.VolumeSpecific v_B=298.15*U.K/B.p0
          "<html>Specific volume of the 2<sup>nd</sup> species</html>"
          annotation (Dialog(__Dymola_label="<html><i>v</i><sub>B</sub></html>"));
        replaceable function pDstar =
            MobilityFactors.BaseClasses.pDstar_nonpolar
          "Reduced pressure-diffusivity product";
        replaceable package A = Characteristics.BaseClasses.Characteristic
          "<html>Characteristic data of the 2<sup>nd</sup> species</html>";
        replaceable package B = Characteristics.BaseClasses.Characteristic
          "<html>Characteristic data of the 2<sup>nd</sup> species</html>";

        input Q.TemperatureAbsolute T_crit
          "Geometric mean of the critical temperatures" annotation (Dialog(
              __Dymola_label="<html><i>T</i><sub>crit</sub></html>"));
        input Q.PressureAbsolute p_crit
          "Geometric mean of the critical pressures" annotation (Dialog(
              __Dymola_label="<html><i>p</i><sub>crit</sub></html>"));
        output Q.NumberAbsolute k_Phi "Binary mobility factor" annotation (
            Dialog(__Dymola_label="<html><i>k</i><sub>&Phi;</sub></html>"));

      algorithm
        k_Phi := pDstar(T/T_crit)*p_crit^(2/3)*(T_crit*U.mol/U.K)^(5/6)*(U.atm/
          U.q)^(1/3)*U.cm^2/U.s/sqrt(harmonicMean({A.m,B.m})*U.g)/(p_A*B.D(T,
          v_A) + p_B*A.D(T, v_B))
          "Based on [Slattery1958, eq. 5] and the exchange equations in FCSys.Species.Species";

        annotation (Inline=true,Documentation(info="<html><p><i>v</i><sub>A</sub> and <i>v</i><sub>B</sub> are given as inputs even though they can be calculated
 from <i>T</i>, <i>p</i><sub>A</sub>, and <i>p</i><sub>B</sub> because it may be desirable to leave <i>v</i><sub>A</sub> and <i>v</i><sub>B</sub> constant while
  varying <i>p</i><sub>A</sub> and <i>p</i><sub>B</sub>.</p>

  <p>This function is based on eq. 5 from [<a href=\"modelica://FCSys.UsersGuide.References.Slattery1958\">Slattery1958</a>].</p></html>"));
      end k;

      function pDstar_nonpolar
        "[Slattery1958] reduced pressure-diffusivity product for nonpolar species"
        input Q.TemperatureAbsolute Tstar=1
          "Ratio of temperature to the geometric mean of the critical temperatures"
          annotation (Dialog(__Dymola_label=
                "<html><i>T</i><sub>c AB</sub></html>"));
        output Q.NumberAbsolute pDstar "Reduced pressure-diffusivity product";

      algorithm
        pDstar := 3.882e-4*Tstar^1.823;
        annotation (Inline=true,Documentation(info=
                "<html><p>This function is based on eq. 4 from [<a href=\"modelica://FCSys.UsersGuide.References.Slattery1958\">Slattery1958</a>].</p></html>"));
      end pDstar_nonpolar;

      function pDstar_polar
        "[Slattery1958] reduced pressure-diffusivity product for polar species"
        input Q.TemperatureAbsolute Tstar=1
          "Ratio of temperature to the geometric mean of the critical temperatures"
          annotation (Dialog(__Dymola_label=
                "<html><i>T</i><sub>c AB</sub></html>"));
        output Q.NumberAbsolute pDstar "Reduced pressure-diffusivity product";

      algorithm
        pDstar := 5.148e-4*Tstar^2.334;
        annotation (Inline=true,Documentation(info=
                "<html><p>This function is based on eq. 9 from [<a href=\"modelica://FCSys.UsersGuide.References.Slattery1958\">Slattery1958</a>].</p></html>"));
      end pDstar_polar;
    end BaseClasses;
  end MobilityFactors;
  annotation (Documentation(info="<html>
  <p>Each species has a subpackage for each material phase in which the species
  is represented.  The thermodynamic properties are generally different for each phase.</p>

<p>Additional materials may be included as needed.  The thermodynamic data for
  materials that are condensed at standard conditions is available in
  [<a href=\"modelica://FCSys.UsersGuide.References.McBride2002\">McBride2002</a>].
  The thermodynamic data for materials
  that are gases at standard conditions is available in
  <a href=\"modelica://Modelica.Media.IdealGases.Common.SingleGasesData\">Modelica.Media.IdealGases.Common.SingleGasesData</a>
  (and [<a href=\"modelica://FCSys.UsersGuide.References.McBride2002\">McBride2002</a>]). Virial coefficients are available in
  [<a href=\"modelica://FCSys.UsersGuide.References.Dymond2002\">Dymond2002</a>].  Transport characteristics are available in
  [<a href=\"modelica://FCSys.UsersGuide.References.McBride1996\">McBride1996</a>,
  <a href=\"modelica://FCSys.UsersGuide.References.McBride2002\">McBride2002</a>].</p>

  <p><b>Licensed by the Hawaii Natural Energy Institute under the Modelica License 2</b><br>
Copyright &copy; 2007&ndash;2014, <a href=\"http://www.hnei.hawaii.edu/\">Hawaii Natural Energy Institute</a> and <a href=\"http://www.gtrc.gatech.edu/\">Georgia Tech Research Corporation</a>.</p>

<p><i>This Modelica package is <u>free</u> software and the use is completely at <u>your own risk</u>;
it can be redistributed and/or modified under the terms of the Modelica License 2. For license conditions (including the
disclaimer of warranty) see <a href=\"modelica://FCSys.UsersGuide.License\">
FCSys.UsersGuide.License</a> or visit <a href=\"http://www.modelica.org/licenses/ModelicaLicense2\">
http://www.modelica.org/licenses/ModelicaLicense2</a>.</i></p>
</html>"), Icon(graphics={
        Line(
          points={{-76,-80},{-62,-30},{-32,40},{4,66},{48,66},{73,45},{62,-8},{
              48,-50},{38,-80}},
          color={64,64,64},
          smooth=Smooth.Bezier),
        Line(
          points={{-40,20},{68,20}},
          color={175,175,175},
          smooth=Smooth.None),
        Line(
          points={{-40,20},{-44,88},{-44,88}},
          color={175,175,175},
          smooth=Smooth.None),
        Line(
          points={{68,20},{86,-58}},
          color={175,175,175},
          smooth=Smooth.None),
        Line(
          points={{-60,-28},{56,-28}},
          color={175,175,175},
          smooth=Smooth.None),
        Line(
          points={{-60,-28},{-74,84},{-74,84}},
          color={175,175,175},
          smooth=Smooth.None),
        Line(
          points={{56,-28},{70,-80}},
          color={175,175,175},
          smooth=Smooth.None),
        Line(
          points={{-76,-80},{38,-80}},
          color={175,175,175},
          smooth=Smooth.None),
        Line(
          points={{-76,-80},{-94,-16},{-94,-16}},
          color={175,175,175},
          smooth=Smooth.None)}));
end Characteristics;

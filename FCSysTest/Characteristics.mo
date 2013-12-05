within FCSysTest;
package Characteristics
  extends Modelica.Icons.ExamplesPackage;
  import FCSys.Characteristics.BaseClasses.ReferenceEnthalpy;
  model TestCellPotentialsGas
    "<html>Test the potentials of the reaction 2H<sub>2</sub> + O<sub>2</sub> &#8652; 2H<sub>2</sub>O<sub>(g)</sub></html>"
    import FCSys.Characteristics.*;
    import FCSysTest.Test.assertValues;
    extends Modelica.Icons.Example;
    parameter FCSys.Quantities.TemperatureAbsolute T[:]={298,373.15,473.15,
        673.15,873.15,1073.15,1273.15}*FCSys.Units.K "Temperatures"
      annotation (Dialog(__Dymola_label="<html><i>T</i></html>"));
    final parameter FCSys.Quantities.Potential w_OC_model[:]=w_OC(T)
      "Correlated open circuit potentials" annotation (Dialog(__Dymola_label=
            "<html><i>w</i><sub>OC model</sub></html>"));
    // Note:  The potentials are scaled in terms of electrons.
    parameter FCSys.Quantities.Potential w_OC_table[size(T, 1)]=0.5*{-228590,-225.2e3,
        -220.4e3,-210.3e3,-199.6e3,-188.6e3,-177.4e3}*FCSys.Units.J/FCSys.Units.mol
      "Tabulated open circuit potentials" annotation (Dialog(__Dymola_label=
            "<html><i>w</i><sub>OC table</sub></html>"));
    // The first entry is based on [Moran2004, p. 803].  The others are
    // from [Larminie2003, p. 28].

  protected
    replaceable function w_OC "Open-circuit voltage"
      input FCSys.Quantities.TemperatureAbsolute T "Temperature";
      input FCSys.Quantities.PressureAbsolute p=FCSys.Units.atm "Pressure";
      output FCSys.Quantities.Potential w_OC "Potential";

    algorithm
      w_OC := 0.5*(H2O.Gas.g(T, p) - H2.Gas.g(T, p) - 0.5*O2.Gas.g(T, p))
        annotation (Inline=true);

    end w_OC;

  initial equation
    assertValues(
        w_OC_model,
        w_OC_table,
        1e-3*FCSys.Units.V,
        name="open circuit potential");
    // Note:  In Dymola 7.4, the w_OC() function call can't be used
    // directly here.  Instead, intermediate variables must be used.
    // Otherwise, the result is different.

  end TestCellPotentialsGas;

  model TestCellPotentialsLiquid
    "<html>Test the potentials of the reaction 2H<sub>2</sub> + O<sub>2</sub> &#8652; 2H<sub>2</sub>O<sub>(l)</sub></html>"
    import FCSys.Characteristics.*;
    import FCSysTest.Test.assertValues;
    extends TestCellPotentialsGas(T={298,298.15,353.15}*FCSys.Units.K,
        w_OC_table=0.5*{-237180,-237.2e3,-228.2e3}*FCSys.Units.J/FCSys.Units.mol);
    final parameter FCSys.Quantities.Potential w_therm_model[:]=w_therm(T)
      "Correlated thermodynamic potentials";
    parameter FCSys.Quantities.Potential w_therm_table[size(T, 1)]=0.5*{-285830,
        -237.2e3/0.83,-228.2e3/0.80}*FCSys.Units.J/FCSys.Units.mol
      "Tabulated thermodynamic potentials" annotation (Dialog(__Dymola_label=
            "<html><i>w</i><sub>therm table</sub></html>"));
    // The first entry is based on [Moran2004, p. 803].  The others are
    // from [Larminie2003, pp. 28 & 33].

  protected
    redeclare function w_OC "Open-circuit voltage"
      input FCSys.Quantities.TemperatureAbsolute T "Temperature";
      input FCSys.Quantities.PressureAbsolute p=FCSys.Units.atm "Pressure";
      output FCSys.Quantities.Potential w_OC "Potential";

    algorithm
      w_OC := 0.5*(H2O.Liquid.g(T, p) - H2.Gas.g(T, p) - 0.5*O2.Gas.g(T, p))
        annotation (Inline=true);

    end w_OC;

    function w_therm "Thermodynamic potential"
      input FCSys.Quantities.TemperatureAbsolute T "Temperature";
      input FCSys.Quantities.PressureAbsolute p=FCSys.Units.atm "Pressure";
      output FCSys.Quantities.Potential w_therm "Potential";

    algorithm
      w_therm := 0.5*(H2O.Liquid.h(T, p) - H2.Gas.h(T, p) - 0.5*O2.Gas.h(T, p))
        annotation (Inline=true);

    end w_therm;

  initial equation
    assertValues(
        w_therm_model,
        w_therm_table,
        1e-2*FCSys.Units.V,
        name="thermodynamic potential");

  end TestCellPotentialsLiquid;

  package H2
    extends Modelica.Icons.ExamplesPackage;
    package Gas
      extends Modelica.Icons.ExamplesPackage;

      model c_p
        "<html>Test the isobaric specific heat capacity of H<sub>2</sub> gas against [<a href=\"modelica://FCSys.UsersGuide.References\">Moran2004</a>, p. 794]</html>"
        import FCSysTest.Test.assertValue;
        extends Modelica.Icons.Example;

        // Note:  To work in Dymola 7.4., this and other tests of the
        // Characteristics package must be models rather than functions.

        replaceable package Data = FCSys.Characteristics.H2.Gas (b_v=[1], n_v={
                -1,0}) "Ideal gas properties";
        parameter FCSys.Quantities.NumberAbsolute eps=2e-3
          "Relative error tolerance";
        parameter FCSys.Quantities.TemperatureAbsolute T[:]={250,300,350,400,
            600,800,1000}*FCSys.Units.K "Temperatures"
          annotation (Dialog(__Dymola_label="<html><i>T</i></html>"));
        final parameter FCSys.Quantities.CapacityThermalSpecific c_p_model[:]=
            Data.c_p(T) "Correlated isobaric specific heat capacity";
        parameter FCSys.Quantities.CapacityThermalSpecific c_p_table[size(T, 1)]
          =Data.m*{14.051,14.307,14.427,14.476,14.546,14.695,14.983}*FCSys.Units.J
            /(FCSys.Units.g*FCSys.Units.K)
          "Tabulated isobaric specific heat capacity" annotation (Dialog(
              __Dymola_label="<html><i>c</i><sub><i>p</i> table</sub></html>"));

      initial equation
        for i in 1:size(T, 1) loop
          assertValue(
                  actual=c_p_model[i],
                  expected=c_p_table[i],
                  eps=eps*c_p_table[i],
                  name="of isobaric specific heat capacity of " + Data.formula
               + " as ideal gas at " + String(T[i]/FCSys.Units.K) + " K");
        end for;

      end c_p;

      model c_v
        "<html>Test the isochoric specific heat capacity of H<sub>2</sub> gas against [<a href=\"modelica://FCSys.UsersGuide.References\">Moran2004</a>, p. 794]</html>"
        import FCSysTest.Test.assertValue;
        extends Modelica.Icons.Example;
        replaceable package Data = FCSys.Characteristics.H2.Gas (b_v=[1], n_v={
                -1,0}) "Ideal gas properties";
        parameter FCSys.Quantities.NumberAbsolute eps=2e-3
          "Relative error tolerance";
        parameter FCSys.Quantities.TemperatureAbsolute T[:]={250,300,350,400,
            600,800,1000}*FCSys.Units.K "Temperatures"
          annotation (Dialog(__Dymola_label="<html><i>T</i></html>"));
        final parameter FCSys.Quantities.CapacityThermalSpecific c_v_model[:]=
            Data.c_v(T) "Correlated isochoric specific heat capacity";
        parameter FCSys.Quantities.CapacityThermalSpecific c_v_table[size(T, 1)]
          =Data.m*{9.927,10.183,10.302,10.352,10.422,10.570,10.859}*FCSys.Units.J
            /(FCSys.Units.g*FCSys.Units.K)
          "Tabulated isochoric specific heat capacity" annotation (Dialog(
              __Dymola_label="<html><i>c</i><sub><i>v</i> table</sub></html>"));

      initial equation
        for i in 1:size(T, 1) loop
          assertValue(
                  actual=c_v_model[i],
                  expected=c_v_table[i],
                  eps=eps*c_v_table[i],
                  name="of isochoric specific heat capacity of " + Data.formula
               + " as ideal gas at " + String(T[i]/FCSys.Units.K) + " K");
        end for;

      end c_v;

      model zeta
        "<html>Test the fluidity of H<sub>2</sub> gas against [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 919&ndash;920]</html>"
        import FCSysTest.Test.assertValue;
        extends Modelica.Icons.Example;

        replaceable package Data = FCSys.Characteristics.H2.Gas
          "Material characteristics";
        parameter FCSys.Quantities.NumberAbsolute eps=0.1
          "Relative error tolerance";
        parameter FCSys.Quantities.TemperatureAbsolute T[:]={200,250,300,350,
            400,600,800,1000,2000}*FCSys.Units.K "Temperatures"
          annotation (Dialog(__Dymola_label="<html><i>T</i></html>"));
        parameter FCSys.Quantities.Fluidity zeta_table[size(T, 1)]={1/68.1e-7,1
            /78.9e-7,1/89.6e-7,1/98.8e-7,1/108.2e-7,1/142.4e-7,1/172.4e-7,1/
            201.3e-7,1/318.2e-7}/(FCSys.Units.Pa*FCSys.Units.s)
          "Tabulated fluidity" annotation (Dialog(__Dymola_label=
                "<html>&zeta;<sub>table</sub></html>"));

      initial equation
        for i in 1:size(T, 1) loop
          assertValue(
                  actual=Data.zeta(T[i]),
                  expected=zeta_table[i],
                  eps=eps*zeta_table[i],
                  name="of fluidity of " + Data.formula + " at " + String(T[i]/
              FCSys.Units.K) + " K");
        end for;

      end zeta;

      model theta
        "<html>Test the thermal resistivity of H<sub>2</sub> gas against [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 919&ndash;920]</html>"
        import FCSysTest.Test.assertValue;
        extends Modelica.Icons.Example;
        replaceable package Data = FCSys.Characteristics.H2.Gas
          "Material characteristics";
        parameter FCSys.Quantities.NumberAbsolute eps=0.2
          "Relative error tolerance";
        parameter FCSys.Quantities.TemperatureAbsolute T[:]={200,250,300,350,
            400,600,800,1000,2000}*FCSys.Units.K "Temperatures"
          annotation (Dialog(__Dymola_label="<html><i>T</i></html>"));
        parameter FCSys.Quantities.ResistivityThermal theta_table[size(T, 1)]={
            1/0.131,1/0.131,1/0.183,1/0.204,1/0.226,1/0.305,1/0.378,1/0.448,1/
            0.878}*FCSys.Units.m*FCSys.Units.K/FCSys.Units.W
          "Tabulated thermal resistivity" annotation (Dialog(__Dymola_label=
                "<html>&theta;<sub>table</sub></html>"));

      initial equation
        for i in 1:size(T, 1) loop
          assertValue(
                  actual=Data.theta(T[i]),
                  expected=theta_table[i],
                  eps=eps*theta_table[i],
                  name="of thermal resistivity of " + Data.formula + " at " +
              String(T[i]/FCSys.Units.K) + " K");
        end for;

      end theta;

    end Gas;

  end H2;

  package H2O
    extends Modelica.Icons.ExamplesPackage;

    model TestSaturationPressure
      "<html>Test the saturation pressure of H<sub>2</sub>O against [<a href=\"modelica://FCSys.UsersGuide.References\">Moran2004</a>, pp. 760-761]</html>"
      import FCSysTest.Test.assertValue;
      import FCSys.Characteristics.H2O;
      extends Modelica.Icons.Example;

      parameter FCSys.Quantities.TemperatureAbsolute T[:]=FCSys.Units.from_degC(
          {0.01,25,50,80,100,150,200}) "Temperatures"
        annotation (Dialog(__Dymola_label="<html><i>T</i></html>"));
      parameter FCSys.Quantities.PressureAbsolute p_sat[size(T, 1)]={0.00611,
          0.03169,0.1235,0.4739,1.014,4.758,15.54}*FCSys.Units.bar
        "Saturation pressures";
      FCSys.Quantities.PressureAbsolute p[size(T, 1)](each start=FCSys.Units.atm)
        "Pressures" annotation (Dialog(__Dymola_label="<html><i>p</i></html>"));

    initial equation
      for i in 1:size(T, 1) loop
        assertValue(
              p[i],
              p_sat[i],
              eps=0.01*p_sat[i],
              name="of saturation pressure at " + String(FCSys.Units.to_degC(T[
            i])) + " deg C");
      end for;

    equation
      H2O.Gas.g(T, p) = H2O.Liquid.g(T, p) "Chemical/phase equilibrium";

    end TestSaturationPressure;

    package Gas
      extends Modelica.Icons.ExamplesPackage;

      model h
        "<html>Test the specific enthalpy of H<sub>2</sub>O gas against [<a href=\"modelica://FCSys.UsersGuide.References\">Moran2004</a>, pp. 799&ndash;801]</html>"
        import FCSysTest.Test.assertValue;
        extends Modelica.Icons.Example;
        replaceable package Data = FCSys.Characteristics.H2O.Gas (
            b_v=[1],
            n_v={-1,0},
            referenceEnthalpy=ReferenceEnthalpy.zeroAt0K)
          "Ideal gas properties w/ 0K enthalpy reference";
        parameter FCSys.Quantities.NumberAbsolute eps=0.02
          "Relative error tolerance";
        parameter FCSys.Quantities.TemperatureAbsolute T[:]={220,300,400,600,
            800,1000,2000,3250}*FCSys.Units.K "Temperatures"
          annotation (Dialog(__Dymola_label="<html><i>T</i></html>"));
        final parameter FCSys.Quantities.Potential h_model[:]=Data.h(T)
          "Correlated specific enthalpy";
        parameter FCSys.Quantities.Potential h_table[size(T, 1)]={7295,9966,
            13356,20402,27896,35882,82593,150272}*FCSys.Units.J/FCSys.Units.mol
          "Tabulated specific enthalpy" annotation (Dialog(__Dymola_label=
                "<html><i>h</i><sub>table</sub></html>"));

      initial equation
        for i in 1:size(T, 1) loop
          assertValue(
                  actual=h_model[i],
                  expected=h_table[i],
                  eps=eps*h_table[i],
                  name="of specific enthalpy of " + Data.formula +
              " as ideal gas at " + String(T[i]/FCSys.Units.K) + " K");
        end for;

      end h;

      model s
        "<html>Test the specific entropy of H<sub>2</sub>O gas against [<a href=\"modelica://FCSys.UsersGuide.References\">Moran2004</a>, pp. 799&ndash;801]</html>"
        import FCSysTest.Test.assertValue;
        extends Modelica.Icons.Example;
        replaceable package Data = FCSys.Characteristics.H2O.Gas (b_v=[1], n_v=
                {-1,0}) "Ideal gas properties w/ 0K enthalpy reference";
        parameter FCSys.Quantities.NumberAbsolute eps=3e-3
          "Relative error tolerance";
        parameter FCSys.Quantities.TemperatureAbsolute T[:]={220,300,400,600,
            800,1000,2000,3250}*FCSys.Units.K "Temperatures"
          annotation (Dialog(__Dymola_label="<html><i>T</i></html>"));
        final parameter FCSys.Quantities.NumberAbsolute s_model[:]=Data.s(T)
          "Correlated specific entropy";
        parameter FCSys.Quantities.NumberAbsolute s_table[size(T, 1)]={178.576,
            188.928,198.673,212.920,223.693,232.597,264.571,290.756}*FCSys.Units.J
            /(FCSys.Units.mol*FCSys.Units.K) "Tabulated specific entropy"
          annotation (Dialog(__Dymola_label=
                "<html><i>s</i><sub>table</sub></html>"));

      initial equation
        for i in 1:size(T, 1) loop
          assertValue(
                  actual=s_model[i],
                  expected=s_table[i],
                  eps=eps*s_table[i],
                  name="of specific entropy of " + Data.formula +
              " as ideal gas at " + String(T[i]/FCSys.Units.K) + " K");
        end for;

      end s;

      model zeta
        "<html>Test the fluidity of H<sub>2</sub>O gas against [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 924&ndash;925]</html>"
        extends H2.Gas.zeta(
          redeclare package Data = FCSys.Characteristics.H2O.Gas,
          eps=0.1,
          T={373.15,400,600}*FCSys.Units.K,
          zeta_table={1/12.02e-6,1/13.05e-6,1/22.7e-6}/(FCSys.Units.Pa*FCSys.Units.s));

      end zeta;

      model theta
        "<html>Test the thermal resistivity of H<sub>2</sub>O gas against [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 924&ndash;925]</html>"
        extends H2.Gas.theta(
          redeclare package Data = FCSys.Characteristics.H2O.Gas,
          eps=1.1,
          T={373.15,400,600}*FCSys.Units.K,
          theta_table={1/24.8e-3,1/27.2e-3,1/92.9e-3}*FCSys.Units.m*FCSys.Units.K
              /FCSys.Units.W);
        // Note:  Tolerance must be very large to pass check (due to value
        // at 600 K).

      end theta;

    end Gas;

  end H2O;

  package N2
    extends Modelica.Icons.ExamplesPackage;
    package Gas
      extends Modelica.Icons.ExamplesPackage;

      model c_p
        "<html>Test the isobaric specific heat capacity of N<sub>2</sub> gas against [<a href=\"modelica://FCSys.UsersGuide.References\">Moran2004</a>, p. 794]</html>"
        extends H2.Gas.c_p(
          redeclare package Data = FCSys.Characteristics.N2.Gas (b_v=[1], n_v={
                  -1,0}),
          eps=1e-3,
          c_p_table=Data.m*{1.039,1.039,1.041,1.044,1.075,1.121,1.167}*FCSys.Units.J
              /(FCSys.Units.g*FCSys.Units.K));

      end c_p;

      model c_v
        "<html>Test the isochoric specific heat capacity of N<sub>2</sub> gas against [<a href=\"modelica://FCSys.UsersGuide.References\">Moran2004</a>, p. 794]</html>"
        extends H2.Gas.c_v(
          redeclare package Data = FCSys.Characteristics.N2.Gas (b_v=[1], n_v={
                  -1,0}),
          eps=0.02,
          c_v_table=Data.m*{0.742,0.743,0.744,0.757,0.778,0.825,0.870}*FCSys.Units.J
              /(FCSys.Units.g*FCSys.Units.K));

      end c_v;

      model h
        "<html>Test the specific enthalpy of N<sub>2</sub> gas against [<a href=\"modelica://FCSys.UsersGuide.References\">Moran2004</a>, pp. 799&ndash;801]</html>"
        extends H2O.Gas.h(
          redeclare package Data = FCSys.Characteristics.N2.Gas (
              b_v=[1],
              n_v={-1,0},
              referenceEnthalpy=ReferenceEnthalpy.zeroAt0K),
          eps=1e-3,
          h_table={6391,8723,11640,17563,23714,30129,64810,110690}*FCSys.Units.J
              /FCSys.Units.mol);

      end h;

      model s
        "<html>Test the specific entropy of N<sub>2</sub> gas against [<a href=\"modelica://FCSys.UsersGuide.References\">Moran2004</a>, pp. 799&ndash;801]</html>"
        extends H2O.Gas.s(
          redeclare package Data = FCSys.Characteristics.N2.Gas (b_v=[1], n_v={
                  -1,0}),
          eps=1e-4,
          s_table={182.638,191.682,200.071,212.066,220.907,228.057,251.969,
              269.763}*FCSys.Units.J/(FCSys.Units.mol*FCSys.Units.K));

      end s;

      model zeta
        "<html>Test the fluidity of H<sub>2</sub>O gas against [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 920&ndash;921]</html>"
        extends FCSysTest.Characteristics.H2.Gas.zeta(
          redeclare package Data = FCSys.Characteristics.N2.Gas,
          eps=0.1,
          T={200,250,300,350,400,600,800,1000,1200}*FCSys.Units.K,
          zeta_table={1/129.2e-7,1/154.9e-7,1/178.2e-7,1/200.0e-7,1/220.4e-7,1/
              290.8e-7,1/349.1e-7,1/399.9e-7,1/445.3e-7}/(FCSys.Units.Pa*FCSys.Units.s));

      end zeta;

      model theta
        "<html>Test the thermal resistivity of H<sub>2</sub>O gas against [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 920&ndash;921]</html>"
        extends H2.Gas.theta(
          redeclare package Data = FCSys.Characteristics.N2.Gas,
          eps=0.02,
          T={200,250,300,350,400,600,800,1000,1200}*FCSys.Units.K,
          theta_table={1/18.3e-3,1/22.2e-3,1/25.9e-3,1/29.3e-3,1/32.7e-3,1/
              44.6e-3,1/54.8e-3,1/64.7e-3,1/75.8e-3}*FCSys.Units.m*FCSys.Units.K
              /FCSys.Units.W);

      end theta;

    end Gas;

  end N2;

  package O2
    extends Modelica.Icons.ExamplesPackage;
    package Gas
      extends Modelica.Icons.ExamplesPackage;

      model c_p
        "<html>Test the isobaric specific heat capacity of O<sub>2</sub> gas against [<a href=\"modelica://FCSys.UsersGuide.References\">Moran2004</a>, p. 794]</html>"
        extends H2.Gas.c_p(
          redeclare package Data = FCSys.Characteristics.O2.Gas (b_v=[1], n_v={
                  -1,0}),
          eps=1e-3,
          c_p_table=Data.m*{0.913,0.918,0.928,0.941,1.003,1.054,1.090}*FCSys.Units.J
              /(FCSys.Units.g*FCSys.Units.K));

      end c_p;

      model c_v
        "<html>Test the isochoric specific heat capacity of O<sub>2</sub> gas against [<a href=\"modelica://FCSys.UsersGuide.References\">Moran2004</a>, p. 794]</html>"
        extends H2.Gas.c_v(
          redeclare package Data = FCSys.Characteristics.O2.Gas (b_v=[1], n_v={
                  -1,0}),
          eps=1e-3,
          c_v_table=Data.m*{0.653,0.658,0.668,0.681,0.743,0.794,0.830}*FCSys.Units.J
              /(FCSys.Units.g*FCSys.Units.K));

      end c_v;

      model h
        "<html>Test the specific enthalpy of O<sub>2</sub> gas against [<a href=\"modelica://FCSys.UsersGuide.References\">Moran2004</a>, pp. 794, 799&ndash;801]</html>"
        extends H2O.Gas.h(
          redeclare package Data = FCSys.Characteristics.O2.Gas (
              b_v=[1],
              n_v={-1,0},
              referenceEnthalpy=ReferenceEnthalpy.zeroAt0K),
          eps=0.01,
          h_table={6404,8736,11711,17929,24523,31389,67881,116827}*FCSys.Units.J
              /FCSys.Units.mol);

      end h;

      model s
        "<html>Test the specific entropy of O<sub>2</sub> gas against [<a href=\"modelica://FCSys.UsersGuide.References\">Moran2004</a>, pp. 794, 799&ndash;801]</html>"
        extends H2O.Gas.s(
          redeclare package Data = FCSys.Characteristics.O2.Gas (b_v=[1], n_v={
                  -1,0}),
          eps=1e-4,
          s_table={196.171,205.213,213.765,226.346,235.810,243.471,268.655,
              287.614}*FCSys.Units.J/(FCSys.Units.mol*FCSys.Units.K));

      end s;

      model zeta
        "<html>Test the fluidity of O<sub>2</sub>O gas against [<a href=\"modelica://FCSys.UsersGuide.References\">Moran2004</a>, p. 924&ndash;925]</html>"
        extends H2.Gas.zeta(
          redeclare package Data = FCSys.Characteristics.O2.Gas,
          eps=0.1,
          T={200,250,300,350,400,600,800,1000,1200}*FCSys.Units.K,
          zeta_table={1/147.5e-7,1/178.6e-7,1/207.2e-7,1/233.5e-7,1/258.2e-7,1/
              343.7e-7,1/415.2e-7,1/477.0e-7,1/532.5e-7}/(FCSys.Units.Pa*FCSys.Units.s));

      end zeta;

      model theta
        "<html>Test the thermal resistivity of O<sub>2</sub>O gas against [<a href=\"modelica://FCSys.UsersGuide.References\">Moran2004</a>, p. 924&ndash;925]</html>"
        extends H2.Gas.theta(
          redeclare package Data = FCSys.Characteristics.O2.Gas,
          eps=0.1,
          T={200,250,300,350,400,600,800,1000,1200}*FCSys.Units.K,
          theta_table={1/18.3e-3,1/22.6e-3,1/26.8e-3,1/29.6e-3,1/33.0e-3,1/
              47.3e-3,1/58.9e-3,1/71.0e-3,1/81.9e-3}*FCSys.Units.m*FCSys.Units.K
              /FCSys.Units.W);

      end theta;

    end Gas;

  end O2;

  package BaseClasses
    extends Modelica.Icons.ExamplesPackage;
    package Characteristic
      extends Modelica.Icons.ExamplesPackage;

      // TODO: Add test for s().

      model TestProperties
        "<html>Test the <code>z</code>, <code>isCompressible</code>, <code>hasThermalExpansion</code> constants of the <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic\">Characteristic</a> package</html>"
        import FCSys.Characteristics.*;
        extends Modelica.Icons.Example;
        // Note:  In Dymola 7.4, this test must be implemented as a model
        // (not as a function) due to the following error when checking
        // isCompressible and hasThermalExpansion:
        //     "Error: Cannot handle unexpanded expression with iterators."

      initial equation
        // z
        assert(H2O.Gas.z == 0, "z failed on test 1.");
        assert('C+'.Graphite.z == 1, "z failed on test 2.");
        assert('SO3-'.Ionomer.z == -1, "z failed on test 3.");
        // isCompressible
        assert(H2O.Gas.isCompressible, "isCompressible failed on test 1.");
        assert(not H2O.Liquid.isCompressible,
          "isCompressible failed on test 2.");
        // hasThermalExpansion
        assert(H2O.Gas.hasThermalExpansion,
          "hasThermalExpansion failed on test 1.");
        assert(not H2O.Liquid.hasThermalExpansion,
          "hasThermalExpansion failed on test 2.");
        annotation (Documentation(info="<html><p>If this model simulates without failure,
    then the test has passed.</p></html>"));
      end TestProperties;

      model c_p
        "<html>Test <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.c_p\">c_p</a>() based on its relation to <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.s\">s</a>()</html>"
        // This approach is based on [Dassault2010, vol. 2, pp. 300-301].
        extends Modelica.Icons.Example;

        package Data = FCSys.Characteristics.H2O.Gas;
        // The data choice is arbitrary but the b_c values must have sufficient
        // richness.
        // Arguments to functions
        FCSys.Quantities.TemperatureAbsolute T=(300 + 100*time)*FCSys.Units.K
          "Temperature";
        FCSys.Quantities.PressureAbsolute p=FCSys.Units.atm
          "Pressure (must be constant)";
        // Note:  The values are arbitrary but must have sufficient richness.
        // Results of functions
        FCSys.Quantities.Potential s "Specific entropy";
        FCSys.Quantities.Potential y "Integral of (c_p/T)*dT";

      initial equation
        y = s;

      equation
        s = Data.s(T, p);
        T*der(y) = Data.c_p(T, p)*der(T) "c_p = T*(dels/delT)_p";
        assert(abs(s - y)/y < 1e-7, "The relationship is incorrect.");
        // Note:  The simulation tolerance is set to 1e-8.
        annotation (Documentation(info="<html><p>If this model simulates without failure,
    then the test has passed.</p></html>"), experiment(Tolerance=1e-8));
      end c_p;

      model c_v
        "<html>Test <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.c_v\">c_v</a>() based on its relation to <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.h\">h</a>() - p&middot;v</html>"
        // This approach is based on [Dassault2010, vol. 2, pp. 300-301].
        extends Modelica.Icons.Example;

        package Data = FCSys.Characteristics.H2O.Gas;
        // The data choice is arbitrary but the b_c values must have sufficient
        // richness.
        // Arguments to functions
        FCSys.Quantities.TemperatureAbsolute T=(300 + 100*time)*FCSys.Units.K
          "Temperature";
        FCSys.Quantities.VolumeSpecific v=Data.v_Tp(300*FCSys.Units.K, FCSys.Units.atm)
          "Specific volume (must be constant)";
        // Note:  The values are arbitrary but must have sufficient richness.
        // Intermediate variables
        FCSys.Quantities.PressureAbsolute p "Pressure";
        // Results of functions
        FCSys.Quantities.Potential u "Internal potential";
        FCSys.Quantities.Potential y "Integral of c_v*dT";

      initial equation
        y = u;

      equation
        p = Data.p_Tv(T, v);
        u = Data.h(T, p) - p*v;
        der(y) = Data.c_v(T, p)*der(T) "c_v = (delu/delT)_v";
        assert(abs(u - y) < 1e-5*FCSys.Units.V,
          "The relationship is incorrect.");
        // Note:  The simulation tolerance is set to 1e-8.
        annotation (Documentation(info="<html><p>If this model simulates without failure,
    then the test has passed.</p></html>"), experiment(Tolerance=1e-8));
      end c_v;

      model dp
        "<html>Test <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.dp_Tv\">dp_Tv</a>() based on its relation to <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.p_Tv\">p_Tv</a>()</html>"
        // This approach is based on [Dassault2010, vol. 2, pp. 300-301].
        extends Modelica.Icons.Example;
        package Data = FCSys.Characteristics.H2O.Gas;
        // The data choice is arbitrary but the b_v values must have sufficient
        // richness.
        // Arguments to functions
        FCSys.Quantities.TemperatureAbsolute T=(300 + 100*time)*FCSys.Units.K
          "Temperature";
        FCSys.Quantities.VolumeSpecific v=(T/FCSys.Units.atm)*(1 + time)
          "Specific volume";
        // Note:  The values are arbitrary but must have sufficient richness.
        // Results of functions
        FCSys.Quantities.Pressure y1 "Direct result of function";
        FCSys.Quantities.Pressure y2 "Integral of derivative of y1";

      initial equation
        y2 = y1;

      equation
        y1 = Data.p_Tv(T, v);
        der(y2) = Data.dp_Tv(
                T,
                v,
                der(T),
                der(v));
        // Note:  This is equivalent to der(y2) = der(y1), but it must be
        // explicit to ensure that the translator uses the defined derivative
        // instead of the automatically derived one.
        assert(abs(y1 - y2) < 1e-7*FCSys.Units.atm,
          "The derivative is incorrect.");
        // Note:  The simulation tolerance is set to 1e-8.
        annotation (Documentation(info="<html><p>If this model simulates without failure,
    then the test has passed.</p></html>"), experiment(Tolerance=1e-8));
      end dp;

      model p_Tv
        "<html>Test <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.p_Tv\">p_Tv</a>() based on its relation to <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.v_Tp\">v_Tp</a>()</html>"
        extends Modelica.Icons.Example;
        package Data = FCSys.Characteristics.H2.Gas;
        // The data choice is arbitrary but the b_v values must have sufficient
        // richness.
        // Note:  This test fails for H2O.Gas due to numerics (not mathematics).
        // See FCSys.Examples.Correlations.
        // Arguments to functions
        FCSys.Quantities.TemperatureAbsolute T=(300 + 100*time)*FCSys.Units.K
          "Temperature";
        FCSys.Quantities.PressureAbsolute p=FCSys.Units.atm*(1 + time^2)
          "Pressure";
        // Note:  The values are arbitrary but must have sufficient richness.
        // Results of functions
        FCSys.Quantities.Pressure y "Indirectly calculated pressure";

      equation
        y = Data.p_Tv(T, Data.v_Tp(T, p))
          "p_Tv and v_Tp are inverses w.r.t. p and v";
        assert(abs(p - y) < 1e-5*FCSys.Units.atm,
          "The relationship is incorrect.");
        // Note:  The simulation tolerance is set to 1e-8.
        annotation (Documentation(info="<html><p>If this model simulates without failure,
    then the test has passed.</p></html>"), experiment(Tolerance=1e-8));
      end p_Tv;

      model dv
        "<html>Test <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.dv_Tp\">dv_Tp</a>() based on its relation to <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.v_Tp\">v_Tp</a>()</html>"
        // This approach is based on [Dassault2010, vol. 2, pp. 300-301].
        extends Modelica.Icons.Example;
        package Data = FCSys.Characteristics.H2O.Gas;

        // The data choice is arbitrary but the b_v values must have sufficient
        // richness.
        // Arguments to functions
        FCSys.Quantities.TemperatureAbsolute T=(300 + 100*time)*FCSys.Units.K
          "Temperature";
        FCSys.Quantities.PressureAbsolute p=FCSys.Units.atm*(1 + time^2)
          "Pressure";
        // Note:  The values are arbitrary but must have sufficient richness.
        // Results of functions
        FCSys.Quantities.VolumeSpecific y1 "Direct result of function";
        FCSys.Quantities.VolumeSpecific y2 "Integral of derivative of y1";

      initial equation
        y2 = y1;

      equation
        y1 = Data.v_Tp(T, p);
        der(y2) = Data.dv_Tp(
                T,
                p,
                der(T),
                der(p));
        // Note:  This is equivalent to der(y2) = der(y1), but it must be
        // explicit to ensure that the translator uses the defined derivative
        // instead of the automatically derived one.
        assert(abs(y1 - y2) < 1e-6*(300*FCSys.Units.K/FCSys.Units.atm),
          "The derivative is incorrect.");
        // Note:  The simulation tolerance is set to 1e-8.
        annotation (Documentation(info="<html><p>If this model simulates without failure,
    then the test has passed.</p></html>"), experiment(Tolerance=1e-8));
      end dv;

      model h
        "<html>Test <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.h\">h</a>(<i>T</i>, <i>p</i>) based on its relation to <i>T</i>, <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.s\">s</a>(<i>T</i>, <i>p</i>), <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.v_Tp\">v_Tp</a>(<i>T</i>, <i>p</i>), and <i>p</i></html>"
        // This approach is based on [Dassault2010, vol. 2, pp. 300-301].
        extends Modelica.Icons.Example;
        package Data = FCSys.Characteristics.H2O.Gas;
        // The data choice is arbitrary but the b_c values must have sufficient
        // richness.
        // Arguments to functions
        FCSys.Quantities.TemperatureAbsolute T=(300 + 100*time)*FCSys.Units.K
          "Temperature";
        FCSys.Quantities.PressureAbsolute p=FCSys.Units.atm*(1 + time^2)
          "Pressure";
        // Note:  The values are arbitrary but must have sufficient richness.
        // Results of functions
        FCSys.Quantities.Potential dh "Direct derivative of h";
        FCSys.Quantities.Potential y "Indirect derivative of h";

      equation
        dh = der(Data.h(T, p));
        y = T*der(Data.s(T, p)) + Data.v_Tp(T, p)*der(p) "dh = T*ds + v*dp";
        assert(abs(dh - y) < 1e-16*FCSys.Units.V,
          "The relationship is incorrect.");
        // Note:  The simulation tolerance is set to 1e-8.
        annotation (Documentation(info="<html><p>If this model simulates without failure,
    then the test has passed.</p></html>"), experiment(Tolerance=1e-8));
      end h;

      model zeta
        "<html>Test the rigid-sphere estimate of <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.zeta\">&zeta;</a>() against the actual fluidity of H<sub>2</sub></html>"
        import FCSysTest.Test.assertLogValue;
        extends Modelica.Icons.Example;
        import DataH2 = FCSys.Characteristics.H2.Gas;

        package Data = FCSys.Characteristics.BaseClasses.Characteristic (m=
                DataH2.m, d=DataH2.d)
          "Properties to estimate fluidity via rigid-sphere assumption";
        constant FCSys.Quantities.Fluidity zeta=Data.zeta(300*FCSys.Units.K);

      initial equation
        assertLogValue(
                actual=zeta,
                expected=1/(89.6e-7*FCSys.Units.Pa*FCSys.Units.s),
                o=0.4);
        // The fluidity is from [Incropera2002, p. 919].
        annotation (Documentation(info="<html><p>If this model simulates without failure,
    then the test has passed.</p></html>"), experiment(Tolerance=1e-8));
      end zeta;

      model theta
        "<html>Test the rigid-sphere estimate of <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.theta\">&theta;</a>() against the actual thermal resistivity of H<sub>2</sub></html>"
        import FCSysTest.Test.assertLogValue;
        extends Modelica.Icons.Example;
        import DataH2 = FCSys.Characteristics.H2.Gas;
        package Data = FCSys.Characteristics.BaseClasses.Characteristic (
            formula=DataH2.formula,
            phase=DataH2.phase,
            m=DataH2.m,
            d=DataH2.d,
            b_c=DataH2.b_c)
          "Properties to estimate thermal resistivity via rigid-sphere assumption";
        constant FCSys.Quantities.ResistivityThermal theta=Data.theta(300*FCSys.Units.K);

      initial equation
        assertLogValue(
                actual=theta,
                expected=(1/0.183)*FCSys.Units.m*FCSys.Units.K/FCSys.Units.W,
                o=0.7);
        // The thermal resistivity is from [Incropera2002, p. 919].
        annotation (Documentation(info="<html><p>If this model simulates without failure,
    then the test has passed.</p></html>"), experiment(Tolerance=1e-8));
      end theta;

    end Characteristic;

  end BaseClasses;

end Characteristics;

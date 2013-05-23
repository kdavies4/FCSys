within FCSys;
package Tests "Models and functions for test and validation"
  extends Modelica.Icons.Package;
  function callAll
    "<html>Call all of the test functions for <a href=\"modelica://FCSys\">FCSys</a></html>"
    extends Modelica.Icons.Function;
    output Boolean ok "true, if all tests passed";
  algorithm
    ok := Units.callAll() and BaseClasses.Utilities.callAll();
    annotation (Documentation(info="<html><p>This function call will fail if any of the functions return an
  incorrect result.  It will return <code>true</code> if all of the functions pass.
  There are no inputs.</p></html>"));
  end callAll;

  model RunAll
    "<html>Run all of the test models for <a href=\"modelica://FCSys\">FCSys</a></html>"
    extends Modelica.Icons.Example;
    Subregions.RunAll testSubregions;
    Characteristics.RunAll testCharacteristics;
    BaseClasses.Utilities.Polynomial.RunAll testBaseClassesUtilities;
    annotation (Documentation(info="<html><p>If this model simulates without failure,
  then the test has passed.</p></html>"), experiment(Tolerance=1e-8));
  end RunAll;

  package Subregions
    extends Modelica.Icons.Package;
    model RunAll
      "<html>Run all of the test models for the <a href=\"modelica://FCSys.Subregions\">Subregions</a> package</html>"
      extends Modelica.Icons.Example;
      Subregion testSubregion;
      Test2Subregions test2Subregions;
      annotation (Documentation(info="<html><p>If this model simulates without failure,
  then the test has passed.</p></html>"), experiment(Tolerance=1e-8));
    end RunAll;

    model Subregion "Test a single subregion"
      extends FCSys.Subregions.Examples.Subregion(
        'inclC+'=true,
        'inclC19HF37O5S-'=true,
        inclH2=true,
        inclN2=true,
        inclO2=true);
      // **fix singularity and include these:
      // **'incle-'=true,
      // **inclH2O=true,
      // **excluded to prevent reactions: 'inclH+'=true (create separate model to test reactions)
      // Currently, there are no assertions.  This model just checks that the
      // simulation runs.
    end Subregion;

    model Test2Subregions
      "Test two subregions with an initial pressure gradient"
      extends FCSys.Subregions.Examples.Subregions(
        n_x=0,
        'inclC+'=true,
        'inclC19HF37O5S-'=true,
        inclH2=true,
        inclN2=true,
        inclO2=true,
        environment(final analysis=true));
      // **fix singularity and include these:
      // 'incle-'=true,
      // **inclH2O=true,
      // **excluded to prevent reactions: 'inclH+'=true,
      output Q.Amount S(stateSelect=StateSelect.never) = subregion1.graphite.
        'C+'.S + subregion2.graphite.'C+'.S + subregion1.ionomer.'C19HF37O5S-'.S
         + subregion2.ionomer.'C19HF37O5S-'.S + subregion1.gas.H2.S +
        subregion2.gas.H2.S + subregion1.gas.N2.S + subregion2.gas.N2.S +
        subregion1.gas.O2.S + subregion2.gas.O2.S "Total entropy";
    equation
      assert(der(S) >= 0, "Entropy may not decrease.");
    end Test2Subregions;
  end Subregions;

  package Characteristics
    extends Modelica.Icons.Package;
    model RunAll
      "<html>Run all of the test models for the <a href=\"modelica://FCSys.Characteristics\">Characteristics</a> package</html>"
      extends Modelica.Icons.Example;
      TestCellPotentialsGas testCellPotentialsGas;
      TestCellPotentialsLiquid testCellPotentialsLiquid;
      H2.Gas.RunAll testH2Gas;
      H2O.RunAll testH2O;
      N2.Gas.RunAll testN2Gas;
      O2.Gas.RunAll testO2Gas;
      BaseClasses.Characteristic.RunAll testBaseClassesCharacteristic;
      annotation (Documentation(info="<html><p>If this model simulates without failure,
  then the test has passed.</p></html>"), experiment(Tolerance=1e-8));
    end RunAll;

    model TestCellPotentialsGas
      "<html>Test the potentials of the reaction 2H<sub>2</sub> + O<sub>2</sub> &#8652; 2H<sub>2</sub>O<sub>(g)</sub></html>"
      import FCSys.Characteristics.*;
      import FCSys.Test.assertValues;
      extends Modelica.Icons.Example;
      parameter Q.TemperatureAbsolute T[:]={298,373.15,473.15,673.15,873.15,
          1073.15,1273.15}*U.K "Temperatures";
      final parameter Q.Potential v_OC_model[:]=v_OC(T)
        "Correlated open circuit potentials";
      // Note:  The potentials are scaled in terms of electrons.
      parameter Q.Potential v_OC_table[size(T, 1)]=0.5*{-228590,-225.2e3,-220.4e3,
          -210.3e3,-199.6e3,-188.6e3,-177.4e3}*U.J/U.mol
        "Tabulated open circuit potentials";
      // The first entry is based on [Moran2004, p. 803].  The others are
      // from [Larminie2003, p. 28].
    protected
      replaceable function v_OC "Open-circuit voltage"
        input Q.TemperatureAbsolute T "Temperature";
        input Q.PressureAbsolute p=1*U.atm "Pressure";
        output Q.Potential v_OC "Potential";
      algorithm
        v_OC := 0.5*(H2O.Gas.g(T, p) - H2.Gas.g(T, p) - 0.5*O2.Gas.g(T, p))
          annotation (Inline=true);
      end v_OC;
    initial equation
      assertValues(
            v_OC_model,
            v_OC_table,
            1e-3*U.V,
            name="open circuit potential");
      // Note:  In Dymola 7.4, the v_OC() function call can't be used
      // directly here.  Instead, intermediate variables must be used.
      // Otherwise, the result is different.
    end TestCellPotentialsGas;

    model TestCellPotentialsLiquid
      "<html>Test the potentials of the reaction 2H<sub>2</sub> + O<sub>2</sub> &#8652; 2H<sub>2</sub>O<sub>(l)</sub></html>"
      import FCSys.Characteristics.*;
      import FCSys.Test.assertValues;
      extends TestCellPotentialsGas(T={298,298.15,353.15}*U.K, v_OC_table=0.5*{
            -237180,-237.2e3,-228.2e3}*U.J/U.mol);
      final parameter Q.Potential v_therm_model[:]=v_therm(T)
        "Correlated thermodynamic potentials";
      parameter Q.Potential v_therm_table[size(T, 1)]=0.5*{-285830,-237.2e3/
          0.83,-228.2e3/0.80}*U.J/U.mol "Tabulated thermodynamic potentials";
      // The first entry is based on [Moran2004, p. 803].  The others are
      // from [Larminie2003, pp. 28 & 33].
    protected
      redeclare function v_OC "Open-circuit voltage"
        input Q.TemperatureAbsolute T "Temperature";
        input Q.PressureAbsolute p=1*U.atm "Pressure";
        output Q.Potential v_OC "Potential";
      algorithm
        v_OC := 0.5*(H2O.Liquid.g(T, p) - H2.Gas.g(T, p) - 0.5*O2.Gas.g(T, p))
          annotation (Inline=true);
      end v_OC;

      function v_therm "Thermodynamic potential"
        input Q.TemperatureAbsolute T "Temperature";
        input Q.PressureAbsolute p=1*U.atm "Pressure";
        output Q.Potential v_therm "Potential";
      algorithm
        v_therm := 0.5*(H2O.Liquid.h(T, p) - H2.Gas.h(T, p) - 0.5*O2.Gas.h(T, p))
          annotation (Inline=true);
      end v_therm;
    initial equation
      assertValues(
            v_therm_model,
            v_therm_table,
            1e-2*U.V,
            name="thermodynamic potential");
    end TestCellPotentialsLiquid;

    package H2
      extends Modelica.Icons.Package;
      package Gas
        extends Modelica.Icons.Package;
        model RunAll
          "<html>Run all of the test models for the <a href=\"modelica://FCSys.Characteristics.N2.Gas\">N2.Gas</a> package</html>"
          extends Modelica.Icons.Example;
          c_p testc_p;
          c_v testc_v;
          F testF;
          R testR;
          annotation (Documentation(info="<html><p>If this model simulates without failure,
  then the test has passed.</p></html>"), experiment(Tolerance=1e-8));
        end RunAll;

        model c_p
          "<html>Test the isobaric specific heat capacity of H<sub>2</sub> gas against [<a href=\"modelica://FCSys.UsersGuide.References\">Moran2004</a>, p. 794]</html>"
          import FCSys.Test.assertValue;
          extends Modelica.Icons.Example;
          replaceable package Data = FCSys.Characteristics.H2.Gas (b_v=[1], n_v
                ={-1,0}) "Ideal gas properties";
          parameter Q.NumberAbsolute eps=2e-3 "Relative error tolerance";
          parameter Q.TemperatureAbsolute T[:]={250,300,350,400,600,800,1000}*U.K
            "Temperatures";
          final parameter Q.CapacityThermalSpecific c_p_model[:]=Data.c_p(T)
            "Correlated isobaric specific heat capacity";
          parameter Q.CapacityThermalSpecific c_p_table[size(T, 1)]=Data.m*{
              14.051,14.307,14.427,14.476,14.546,14.695,14.983}*U.J/(U.g*U.K)
            "Tabulated isobaric specific heat capacity";
        initial equation
          for i in 1:size(T, 1) loop
            assertValue(
                      actual=c_p_model[i],
                      expected=c_p_table[i],
                      eps=eps*c_p_table[i],
                      name="of isobaric specific heat capacity of " + Data.formula
                 + " as ideal gas at " + String(T[i]/U.K) + " K");
          end for;
        end c_p;

        model F
          "<html>Test the fluidity of H<sub>2</sub> gas against [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 919&ndash;920]</html>"
          import FCSys.Test.assertValue;
          extends Modelica.Icons.Example;
          replaceable package Data = FCSys.Characteristics.H2.Gas
            "Material characteristics";
          parameter Q.NumberAbsolute eps=0.1 "Relative error tolerance";
          parameter Q.TemperatureAbsolute T[:]={200,250,300,350,400,600,800,
              1000,2000}*U.K "Temperatures";
          parameter Q.FluidityDynamic F_table[size(T, 1)]={1/68.1e-7,1/78.9e-7,
              1/89.6e-7,1/98.8e-7,1/108.2e-7,1/142.4e-7,1/172.4e-7,1/201.3e-7,1
              /318.2e-7}/(U.Pa*U.s) "Tabulated fluidity";
        initial equation
          for i in 1:size(T, 1) loop
            assertValue(
                      actual=Data.zeta(T[i]),
                      expected=F_table[i],
                      eps=eps*F_table[i],
                      name="of fluidity of " + Data.formula + " at " + String(T[
                i]/U.K) + " K");
          end for;
        end F;

        model c_v
          "<html>Test the isochoric specific heat capacity of H<sub>2</sub> gas against [<a href=\"modelica://FCSys.UsersGuide.References\">Moran2004</a>, p. 794]</html>"
          import FCSys.Test.assertValue;
          extends Modelica.Icons.Example;
          replaceable package Data = FCSys.Characteristics.H2.Gas (b_v=[1], n_v
                ={-1,0}) "Ideal gas properties";
          parameter Q.NumberAbsolute eps=2e-3 "Relative error tolerance";
          parameter Q.TemperatureAbsolute T[:]={250,300,350,400,600,800,1000}*U.K
            "Temperatures";
          final parameter Q.CapacityThermalSpecific c_v_model[:]=Data.c_v(T)
            "Correlated isochoric specific heat capacity";
          parameter Q.CapacityThermalSpecific c_v_table[size(T, 1)]=Data.m*{
              9.927,10.183,10.302,10.352,10.422,10.570,10.859}*U.J/(U.g*U.K)
            "Tabulated isochoric specific heat capacity";
        initial equation
          for i in 1:size(T, 1) loop
            assertValue(
                      actual=c_v_model[i],
                      expected=c_v_table[i],
                      eps=eps*c_v_table[i],
                      name="of isochoric specific heat capacity of " + Data.formula
                 + " as ideal gas at " + String(T[i]/U.K) + " K");
          end for;
        end c_v;

        model R
          "<html>Test the thermal resistivity of H<sub>2</sub> gas against [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 919&ndash;920]</html>"
          import FCSys.Test.assertValue;
          extends Modelica.Icons.Example;
          replaceable package Data = FCSys.Characteristics.H2.Gas
            "Material characteristics";
          parameter Q.NumberAbsolute eps=0.2 "Relative error tolerance";
          parameter Q.TemperatureAbsolute T[:]={200,250,300,350,400,600,800,
              1000,2000}*U.K "Temperatures";
          parameter Q.ResistivityThermal theta_table[size(T, 1)]={1/0.131,1/0.131,1
              /0.183,1/0.204,1/0.226,1/0.305,1/0.378,1/0.448,1/0.878}*U.m*U.K/U.W
            "Tabulated thermal resistivity";
        initial equation
          for i in 1:size(T, 1) loop
            assertValue(
                      actual=Data.R(T[i]),
                      expected=R_table[i],
                      eps=eps*R_table[i],
                      name="of thermal resistivity of " + Data.formula + " at "
                 + String(T[i]/U.K) + " K");
          end for;
        end R;
      end Gas;
    end H2;

    package H2O
      extends Modelica.Icons.Package;
      model RunAll
        "<html>Run all of the test models for the <a href=\"modelica://FCSys.Characteristics.H2O\">H2O</a> package</html>"
        extends Modelica.Icons.Example;
        Gas.RunAll testGas;
        TestSaturationPressure testSaturationPressure;
        annotation (Documentation(info="<html><p>If this model simulates without failure,
  then the test has passed.</p></html>"), experiment(Tolerance=1e-8));
      end RunAll;

      model TestSaturationPressure
        "<html>Test the saturation pressure of H<sub>2</sub>O against [<a href=\"modelica://FCSys.UsersGuide.References\">Moran2004</a>, pp. 760-761]</html>"
        import FCSys.Test.assertValue;
        import FCSys.Characteristics.H2O;
        extends Modelica.Icons.Example;
        parameter Q.TemperatureAbsolute T[:]=U.from_degC({0.01,25,50,80,100,150,
            200}) "Temperatures";
        parameter Q.PressureAbsolute p_sat[size(T, 1)]={0.00611,0.03169,0.1235,
            0.4739,1.014,4.758,15.54}*U.bar "Saturation pressures";
        Q.PressureAbsolute p[size(T, 1)](each start=U.atm) "Pressures";
      initial equation
        for i in 1:size(T, 1) loop
          assertValue(
                  p[i],
                  p_sat[i],
                  eps=0.01*p_sat[i],
                  name="of saturation pressure at " + String(U.to_degC(T[i]))
               + " deg C");
        end for;
      equation
        H2O.Gas.g(T, p) = H2O.Liquid.g(T, p) "Chemical/phase equilibrium";
      end TestSaturationPressure;

      package Gas
        extends Modelica.Icons.Package;
        model RunAll
          "<html>Run all of the test models for the <a href=\"modelica://FCSys.Characteristics.H2O.Gas\">H2O.Gas</a> package</html>"
          extends Modelica.Icons.Example;
          Gas.h testh;
          Gas.s tests;
          F testF;
          R testR;
          annotation (Documentation(info="<html><p>If this model simulates without failure,
  then the test has passed.</p></html>"), experiment(Tolerance=1e-8));
        end RunAll;

        model h
          "<html>Test the specific enthalpy of H<sub>2</sub>O gas against [<a href=\"modelica://FCSys.UsersGuide.References\">Moran2004</a>, pp. 799&ndash;801]</html>"
          import FCSys.Test.assertValue;
          extends Modelica.Icons.Example;
          replaceable package Data = FCSys.Characteristics.H2O.Gas (
              b_v=[1],
              n_v={-1,0},
              h(referenceEnthalpy=ReferenceEnthalpy.ZeroAt0K))
            "Ideal gas properties w/ 0K enthalpy reference";
          parameter Q.NumberAbsolute eps=0.02 "Relative error tolerance";
          parameter Q.TemperatureAbsolute T[:]={220,300,400,600,800,1000,2000,
              3250}*U.K "Temperatures";
          final parameter Q.Potential h_model[:]=Data.h(T)
            "Correlated specific enthalpy";
          parameter Q.Potential h_table[size(T, 1)]={7295,9966,13356,20402,
              27896,35882,82593,150272}*U.J/U.mol "Tabulated specific enthalpy";
        initial equation
          for i in 1:size(T, 1) loop
            assertValue(
                      actual=h_model[i],
                      expected=h_table[i],
                      eps=eps*h_table[i],
                      name="of specific enthalpy of " + Data.formula +
                " as ideal gas at " + String(T[i]/U.K) + " K");
          end for;
        end h;

        model s
          "<html>Test the specific entropy of H<sub>2</sub>O gas against [<a href=\"modelica://FCSys.UsersGuide.References\">Moran2004</a>, pp. 799&ndash;801]</html>"
          import FCSys.Test.assertValue;
          extends Modelica.Icons.Example;
          replaceable package Data = FCSys.Characteristics.H2O.Gas (b_v=[1],
                n_v={-1,0}) "Ideal gas properties w/ 0K enthalpy reference";
          parameter Q.NumberAbsolute eps=3e-3 "Relative error tolerance";
          parameter Q.TemperatureAbsolute T[:]={220,300,400,600,800,1000,2000,
              3250}*U.K "Temperatures";
          final parameter Q.NumberAbsolute s_model[:]=Data.s(T)
            "Correlated specific entropy";
          parameter Q.NumberAbsolute s_table[size(T, 1)]={178.576,188.928,
              198.673,212.920,223.693,232.597,264.571,290.756}*U.J/(U.mol*U.K)
            "Tabulated specific entropy";
        initial equation
          for i in 1:size(T, 1) loop
            assertValue(
                      actual=s_model[i],
                      expected=s_table[i],
                      eps=eps*s_table[i],
                      name="of specific entropy of " + Data.formula +
                " as ideal gas at " + String(T[i]/U.K) + " K");
          end for;
        end s;

        model F
          "<html>Test the fluidity of H<sub>2</sub>O gas against [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 924&ndash;925]</html>"
          extends H2.Gas.F(
            redeclare package Data = FCSys.Characteristics.H2O.Gas,
            eps=0.1,
            T={373.15,400,600}*U.K,
            F_table={1/12.02e-6,1/13.05e-6,1/22.7e-6}/(U.Pa*U.s));
        end F;

        model R
          "<html>Test the thermal resistivity of H<sub>2</sub>O gas against [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 924&ndash;925]</html>"
          extends H2.Gas.R(
            redeclare package Data = FCSys.Characteristics.H2O.Gas,
            eps=1.1,
            T={373.15,400,600}*U.K,
            R_table={1/24.8e-3,1/27.2e-3,1/92.9e-3}*U.m*U.K/U.W);
          // Note:  Tolerance must be very large to pass check (due to value
          // at 600 K).
        end R;
      end Gas;
    end H2O;

    package N2
      extends Modelica.Icons.Package;
      package Gas
        extends Modelica.Icons.Package;
        model RunAll
          "<html>Run all of the test models for the <a href=\"modelica://FCSys.Characteristics.N2.Gas\">N2.Gas</a> package</html>"
          extends Modelica.Icons.Example;
          c_p testc_p;
          c_v testc_v;
          eta testeta;
          h testh;
          s tests;
          F testF;
          R testR;
          annotation (Documentation(info="<html><p>If this model simulates without failure,
  then the test has passed.</p></html>"), experiment(Tolerance=1e-8));
        end RunAll;

        model c_p
          "<html>Test the isobaric specific heat capacity of N<sub>2</sub> gas against [<a href=\"modelica://FCSys.UsersGuide.References\">Moran2004</a>, p. 794]</html>"
          extends H2.Gas.c_p(
            redeclare package Data = FCSys.Characteristics.N2.Gas (b_v=[1], n_v
                  ={-1,0}),
            eps=1e-3,
            c_p_table=Data.m*{1.039,1.039,1.041,1.044,1.075,1.121,1.167}*U.J/(U.g
                *U.K));
        end c_p;

        model c_v
          "<html>Test the isochoric specific heat capacity of N<sub>2</sub> gas against [<a href=\"modelica://FCSys.UsersGuide.References\">Moran2004</a>, p. 794]</html>"
          extends H2.Gas.c_v(
            redeclare package Data = FCSys.Characteristics.N2.Gas (b_v=[1], n_v
                  ={-1,0}),
            eps=0.02,
            c_v_table=Data.m*{0.742,0.743,0.744,0.757,0.778,0.825,0.870}*U.J/(U.g
                *U.K));
        end c_v;

        model eta
          "<html>Test the material resistivity of O<sub>2</sub>O gas against [<a href=\"modelica://FCSys.UsersGuide.References\">Present1958</a>, p. 263]</html>"
          import FCSys.BaseClasses.Utilities.average;
          import FCSys.Test.assertValue;
          extends Modelica.Icons.Example;
          replaceable package Data = FCSys.Characteristics.O2.Gas
            "Material characteristics";
          parameter Q.NumberAbsolute eps=0.7 "Relative error tolerance";
          parameter Q.TemperatureAbsolute T[:]={77.7,194.7,273.2,353.2}*U.K
            "Temperatures";
          parameter Real D_table[size(T, 1)]={0.0168,0.104,average({0.185,0.172}),
              0.287}*U.cm^2/U.s "Tabulated self diffusivity";
          // **Dimension: L2/T
        initial equation
          for i in 1:size(T, 1) loop
            assertValue(
                      actual=U.atm*Data.v_Tp(T[i], U.atm)/Data.eta(T[i], U.atm),
                      expected=D_table[i],
                      eps=eps*D_table[i],
                      name="of material resistivity of " + Data.formula +
                " at " + String(T[i]/U.K) + " K");
          end for;
        end eta;

        model h
          "<html>Test the specific enthalpy of N<sub>2</sub> gas against [<a href=\"modelica://FCSys.UsersGuide.References\">Moran2004</a>, pp. 799&ndash;801]</html>"
          extends H2O.Gas.h(
            redeclare package Data = FCSys.Characteristics.N2.Gas (
                b_v=[1],
                n_v={-1,0},
                h(referenceEnthalpy=ReferenceEnthalpy.ZeroAt0K)),
            eps=1e-3,
            h_table={6391,8723,11640,17563,23714,30129,64810,110690}*U.J/U.mol);
        end h;

        model s
          "<html>Test the specific entropy of N<sub>2</sub> gas against [<a href=\"modelica://FCSys.UsersGuide.References\">Moran2004</a>, pp. 799&ndash;801]</html>"
          extends H2O.Gas.s(
            redeclare package Data = FCSys.Characteristics.N2.Gas (b_v=[1], n_v
                  ={-1,0}),
            eps=1e-4,
            s_table={182.638,191.682,200.071,212.066,220.907,228.057,251.969,
                269.763}*U.J/(U.mol*U.K));
        end s;

        model F
          "<html>Test the fluidity of H<sub>2</sub>O gas against [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 920&ndash;921]</html>"
          extends H2.Gas.F(
            redeclare package Data = FCSys.Characteristics.N2.Gas,
            eps=0.1,
            T={200,250,300,350,400,600,800,1000,1200}*U.K,
            F_table={1/129.2e-7,1/154.9e-7,1/178.2e-7,1/200.0e-7,1/220.4e-7,1/
                290.8e-7,1/349.1e-7,1/399.9e-7,1/445.3e-7}/(U.Pa*U.s));
        end F;

        model R
          "<html>Test the thermal resistivity of H<sub>2</sub>O gas against [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 920&ndash;921]</html>"
          extends H2.Gas.R(
            redeclare package Data = FCSys.Characteristics.N2.Gas,
            eps=0.02,
            T={200,250,300,350,400,600,800,1000,1200}*U.K,
            R_table={1/18.3e-3,1/22.2e-3,1/25.9e-3,1/29.3e-3,1/32.7e-3,1/
                44.6e-3,1/54.8e-3,1/64.7e-3,1/75.8e-3}*U.m*U.K/U.W);
        end R;
      end Gas;
    end N2;

    package O2
      extends Modelica.Icons.Package;
      package Gas
        extends Modelica.Icons.Package;
        model RunAll
          "<html>Run all of the test models for the <a href=\"modelica://FCSys.Characteristics.O2.Gas\">O2.Gas</a> package</html>"
          extends Modelica.Icons.Example;
          c_p testc_p;
          c_v testc_v;
          h testh;
          s tests;
          F testF;
          R testR;
          annotation (Documentation(info="<html><p>If this model simulates without failure,
  then the test has passed.</p></html>"), experiment(Tolerance=1e-8));
        end RunAll;

        model c_p
          "<html>Test the isobaric specific heat capacity of O<sub>2</sub> gas against [<a href=\"modelica://FCSys.UsersGuide.References\">Moran2004</a>, p. 794]</html>"
          extends H2.Gas.c_p(
            redeclare package Data = FCSys.Characteristics.O2.Gas (b_v=[1], n_v
                  ={-1,0}),
            eps=1e-3,
            c_p_table=Data.m*{0.913,0.918,0.928,0.941,1.003,1.054,1.090}*U.J/(U.g
                *U.K));
        end c_p;

        model c_v
          "<html>Test the isochoric specific heat capacity of O<sub>2</sub> gas against [<a href=\"modelica://FCSys.UsersGuide.References\">Moran2004</a>, p. 794]</html>"
          extends H2.Gas.c_v(
            redeclare package Data = FCSys.Characteristics.O2.Gas (b_v=[1], n_v
                  ={-1,0}),
            eps=1e-3,
            c_v_table=Data.m*{0.653,0.658,0.668,0.681,0.743,0.794,0.830}*U.J/(U.g
                *U.K));
        end c_v;

        model eta
          "<html>Test the fluidity of O<sub>2</sub>O gas against [<a href=\"modelica://FCSys.UsersGuide.References\">Moran2004</a>, p. 924&ndash;925]</html>"
          import FCSys.BaseClasses.Utilities.average;
          extends N2.Gas.eta(
            redeclare package Data = FCSys.Characteristics.O2.Gas,
            eps=0.7,
            D_table={0.0153,0.104,average({0.187,0.175}),0.301}*U.cm^2/U.s);
        end eta;

        model h
          "<html>Test the specific enthalpy of O<sub>2</sub> gas against [<a href=\"modelica://FCSys.UsersGuide.References\">Moran2004</a>, pp. 794, 799&ndash;801]</html>"
          extends H2O.Gas.h(
            redeclare package Data = FCSys.Characteristics.O2.Gas (
                b_v=[1],
                n_v={-1,0},
                h(referenceEnthalpy=ReferenceEnthalpy.ZeroAt0K)),
            eps=0.01,
            h_table={6404,8736,11711,17929,24523,31389,67881,116827}*U.J/U.mol);
        end h;

        model s
          "<html>Test the specific entropy of O<sub>2</sub> gas against [<a href=\"modelica://FCSys.UsersGuide.References\">Moran2004</a>, pp. 794, 799&ndash;801]</html>"
          extends H2O.Gas.s(
            redeclare package Data = FCSys.Characteristics.O2.Gas (b_v=[1], n_v
                  ={-1,0}),
            eps=1e-4,
            s_table={196.171,205.213,213.765,226.346,235.810,243.471,268.655,
                287.614}*U.J/(U.mol*U.K));
        end s;

        model F
          "<html>Test the fluidity of O<sub>2</sub>O gas against [<a href=\"modelica://FCSys.UsersGuide.References\">Moran2004</a>, p. 924&ndash;925]</html>"
          extends H2.Gas.F(
            redeclare package Data = FCSys.Characteristics.O2.Gas,
            eps=0.1,
            T={200,250,300,350,400,600,800,1000,1200}*U.K,
            F_table={1/147.5e-7,1/178.6e-7,1/207.2e-7,1/233.5e-7,1/258.2e-7,1/
                343.7e-7,1/415.2e-7,1/477.0e-7,1/532.5e-7}/(U.Pa*U.s));
        end F;

        model R
          "<html>Test the thermal resistivity of O<sub>2</sub>O gas against [<a href=\"modelica://FCSys.UsersGuide.References\">Moran2004</a>, p. 924&ndash;925]</html>"
          extends H2.Gas.R(
            redeclare package Data = FCSys.Characteristics.O2.Gas,
            eps=0.1,
            T={200,250,300,350,400,600,800,1000,1200}*U.K,
            R_table={1/18.3e-3,1/22.6e-3,1/26.8e-3,1/29.6e-3,1/33.0e-3,1/
                47.3e-3,1/58.9e-3,1/71.0e-3,1/81.9e-3}*U.m*U.K/U.W);
        end R;
      end Gas;
    end O2;

    package BaseClasses
      extends Modelica.Icons.Package;
      package Characteristic
        extends Modelica.Icons.Package;
        model RunAll
          "<html>Run all of the test models for the <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic\">Characteristic</a> package</html>"
          extends Modelica.Icons.Example;
          c_p testc_p;
          c_v testc_v;
          dp testdp;
          dv testdv;
          F testF;
          h testh;
          p_Tv testp;
          R testR;
          TestProperties testProperties;
          annotation (Documentation(info="<html><p>If this model simulates without failure,
  then the test has passed.</p></html>"), experiment(Tolerance=1e-8));
        end RunAll;

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
          assert('C19HF37O5S-'.Ionomer.z == -1, "z failed on test 3.");
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
          Q.TemperatureAbsolute T=(300 + 100*time)*U.K "Temperature";
          Q.PressureAbsolute p=U.atm "Pressure (must be constant)";
          // Note:  The values are arbitrary but must have sufficient richness.
          // Results of functions
          Q.Potential s "Specific entropy";
          Q.Potential y "Integral of (c_p/T)*dT";
        initial equation
          y = s;
        equation
          s = Data.s(T, p);
          T*der(y) = Data.c_p(T, p)*der(T) "c_p = T*(dels/delT)_p";
          assert(abs(s - y) < 1e-7*U.V, "The relationship is incorrect.");
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
          Q.TemperatureAbsolute T=(300 + 100*time)*U.K "Temperature";
          Q.VolumeSpecific v=Data.v_Tp(300*U.K, U.atm)
            "Specific volume (must be constant)";
          // Note:  The values are arbitrary but must have sufficient richness.
          // Intermediate variables
          Q.PressureAbsolute p "Pressure";
          // Results of functions
          Q.Potential u "Internal potential";
          Q.Potential y "Integral of c_v*dT";
        initial equation
          y = u;
        equation
          p = Data.p_Tv(T, v);
          u = Data.h(T, p) - p*v;
          der(y) = Data.c_v(T, p)*der(T) "c_v = (delu/delT)_v";
          assert(abs(u - y) < 1e-5*U.V, "The relationship is incorrect.");
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
          Q.TemperatureAbsolute T=(300 + 100*time)*U.K "Temperature";
          Q.VolumeSpecific v=(T/U.atm)*(1 + time) "Specific volume";
          // Note:  The values are arbitrary but must have sufficient richness.
          // Results of functions
          Q.Pressure y1 "Direct result of function";
          Q.Pressure y2 "Integral of derivative of y1";
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
          assert(abs(y1 - y2) < 1e-7*U.atm, "The derivative is incorrect.");
          // Note:  The simulation tolerance is set to 1e-8.
          annotation (Documentation(info="<html><p>If this model simulates without failure,
    then the test has passed.</p></html>"), experiment(Tolerance=1e-8));
        end dp;

        model dv
          "<html>Test <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.dv_Tp\">dv_Tp</a>() based on its relation to <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.v_Tp\">v_Tp</a>()</html>"
          // This approach is based on [Dassault2010, vol. 2, pp. 300-301].
          extends Modelica.Icons.Example;
          package Data = FCSys.Characteristics.H2O.Gas;
          // The data choice is arbitrary but the b_v values must have sufficient
          // richness.
          // Arguments to functions
          Q.TemperatureAbsolute T=(300 + 100*time)*U.K "Temperature";
          Q.PressureAbsolute p=U.atm*(1 + time^2) "Pressure";
          // Note:  The values are arbitrary but must have sufficient richness.
          // Results of functions
          Q.VolumeSpecific y1 "Direct result of function";
          Q.VolumeSpecific y2 "Integral of derivative of y1";
        initial equation
          y2 = y1;
        equation
          y1 = Data.v_Tp(T, p);
          der(y2) = FCSys.Characteristics.BaseClasses.Characteristic.dv_Tp(
                    T,
                    p,
                    der(T),
                    der(p));
          // Note:  This is equivalent to der(y2) = der(y1), but it must be
          // explicit to ensure that the translator uses the defined derivative
          // instead of the automatically derived one.
          assert(abs(y1 - y2) < 1e-6*(300*U.K/U.atm),
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
          Q.TemperatureAbsolute T=(300 + 100*time)*U.K "Temperature";
          Q.PressureAbsolute p=U.atm*(1 + time^2) "Pressure";
          // Note:  The values are arbitrary but must have sufficient richness.
          // Results of functions
          Q.Potential dh "Direct derivative of h";
          Q.Potential y "Indirect derivative of h";
        equation
          dh = der(Data.h(T, p));
          y = T*der(Data.s(T, p)) + Data.v_Tp(T, p)*der(p) "dh = T*ds + v*dp";
          assert(abs(dh - y) < 1e-16*U.V, "The relationship is incorrect.");
          // Note:  The simulation tolerance is set to 1e-8.
          annotation (Documentation(info="<html><p>If this model simulates without failure,
    then the test has passed.</p></html>"), experiment(Tolerance=1e-8));
        end h;

        model F
          "<html>Test the rigid-sphere estimate of <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.F\">F</a>() against the actual fluidity of H<sub>2</sub></html>"
          import FCSys.Test.assertLogValue;
          extends Modelica.Icons.Example;
          import DataH2 = FCSys.Characteristics.H2.Gas;
          package Data = FCSys.Characteristics.BaseClasses.Characteristic (m=
                  DataH2.m, r=DataH2.r)
            "Properties to estimate fluidity via rigid-sphere assumption";
          constant Q.Fluidity zeta=Data.zeta(300*U.K);
        initial equation
          assertLogValue(
                    actual=F,
                    expected=1/(89.6e-7*U.Pa*U.s),
                    o=0.4);
          // The fluidity is from [Incropera2002, p. 919].
          annotation (Documentation(info="<html><p>If this model simulates without failure,
    then the test has passed.</p></html>"), experiment(Tolerance=1e-8));
        end F;

        model p_Tv
          "<html>Test <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.p_Tv\">p_Tv</a>() based on its relation to <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.v_Tp\">v_Tp</a>()</html>"
          extends Modelica.Icons.Example;
          package Data = FCSys.Characteristics.H2.Gas;
          // The data choice is arbitrary but the b_v values must have sufficient
          // richness.
          // Note:  This test fails for H2O.Gas due to numerics (not mathematics).
          // See FCSys.Examples.Correlations.
          // Arguments to functions
          Q.TemperatureAbsolute T=(300 + 100*time)*U.K "Temperature";
          Q.PressureAbsolute p=U.atm*(1 + time^2) "Pressure";
          // Note:  The values are arbitrary but must have sufficient richness.
          // Results of functions
          Q.Pressure y "Indirectly calculated pressure";
        equation
          y = Data.p_Tv(T, Data.v_Tp(T, p))
            "p_Tv and v_Tp are inverses w.r.t. p and v";
          assert(abs(p - y) < 1e-5*U.atm, "The relationship is incorrect.");
          // Note:  The simulation tolerance is set to 1e-8.
          annotation (Documentation(info="<html><p>If this model simulates without failure,
    then the test has passed.</p></html>"), experiment(Tolerance=1e-8));
        end p_Tv;

        model R
          "<html>Test the rigid-sphere estimate of <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.R\">R</a>() against the actual thermal resistivity of H<sub>2</sub></html>"
          import FCSys.Test.assertLogValue;
          extends Modelica.Icons.Example;
          import DataH2 = FCSys.Characteristics.H2.Gas;
          package Data = FCSys.Characteristics.BaseClasses.Characteristic (
              formula=DataH2.formula,
              phase=DataH2.phase,
              m=DataH2.m,
              r=DataH2.r,
              b_c=DataH2.b_c)
            "Properties to estimate thermal resistivity via rigid-sphere assumption";
          constant Q.CapacityThermalSpecific R=Data.R(300*U.K);
        initial equation
          assertLogValue(
                    actual=R,
                    expected=(1/0.183)*U.m*U.K/U.W,
                    o=0.7);
          // The thermal resistivity is from [Incropera2002, p. 919].
          annotation (Documentation(info="<html><p>If this model simulates without failure,
    then the test has passed.</p></html>"), experiment(Tolerance=1e-8));
        end R;
      end Characteristic;
    end BaseClasses;
  end Characteristics;

  package Units
    extends Modelica.Icons.Package;
    function callAll
      "<html>Call all of the test functions for the <a href=\"modelica://FCSys.Units\">Units</a> package</html>"
      extends Modelica.Icons.Function;
      output Boolean ok "true, if all tests passed";
    algorithm
      ok := testValues() and testConversions();
      annotation (Documentation(info="<html><p>This function call will fail if any of the functions return an
  incorrect result.  It will return <code>true</code> if all of the functions pass.
  There are no inputs.</p></html>"));
    end callAll;

    function testValues
      "<html>Test the values of units and constants in the <a href=\"modelica://FCSys.Units\">Units</a> package</html>"
      import FCSys.Units.*;
      extends Modelica.Icons.Function;
      output Boolean ok "true, if all tests passed";
    protected
      function test = FCSys.Test.assertValue (final expected=1);
    algorithm
      ok := false;
      // ----------------------------------------------------------------------
      // Set 1:  Mathematical constants and relations
      test(pi/3.14159265358979323846264338327950288419716939937510, name=
        "1 in set 1");
      // Value from http://en.wikipedia.org/wiki/Pi#Approximate_value
      test(e/2.71828182845904523536028747135266249775724709369995, name=
        "2 in set 1");
      // Value from http://en.wikipedia.org/wiki/E_(mathematical_constant)
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
      // Generated from FCSys/resources/NIST.xls, 2013-1-23
      test(1/alpha/137.035999074, name="inverse fine-structure constant");
      test(1/G_0/(12906.4037217*ohm), name="inverse of conductance quantum");
      test(1/m*h*c/(1.239841930e-6*eV), name=
        "inverse meter-electron volt relationship");
      test(1/m*h*c/(4.556335252755e-8*E_h), name=
        "inverse meter-hartree relationship");
      test(1/m*h*c/(1.986445684e-25*J), name="inverse meter-joule relationship");
      test(1/m*h*c/k_B/(1.4387770e-2*K), name=
        "inverse meter-kelvin relationship");
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
      test(eV/h*cyc/(2.417989348e14*Hz), name=
        "electron volt-hertz relationship");
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
      test(Hz*h/cyc/(4.135667516e-15*eV), name=
        "hertz-electron volt relationship");
      // Factor of 1/cyc due to inconsistencies b/w rad and cyc in [BIPM2006]
      test(Hz*h/cyc/(1.5198298460045e-16*E_h), name=
        "hertz-hartree relationship");
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
      test(kg*c^2/(5.60958885e35*eV), name=
        "kilogram-electron volt relationship");
      test(kg*c^2/(2.061485968e34*E_h), name="kilogram-hartree relationship");
      test(kg*c^2/(8.987551787e16*J), name="kilogram-joule relationship");
      test(kg*c^2/h*cyc/(1.356392608e50*Hz), name="kilogram-hertz relationship");
      // Factor of cyc due to inconsistencies b/w rad and cyc in [BIPM2006]
      test(kg*c^2/k_B/(6.5096582e39*K), name="kilogram-kelvin relationship");
      test(kg*c/h/(4.52443873e41/m), name="kilogram-inverse meter relationship");
      test(mu_0/(12.566370614e-7*N/A^2), name="mag. constant");
      test(N_A/(6.02214129e23/mol), name="Avogadro constant");
      test(N_A*h/(3.9903127176e-10*J*s/mol), name="molar Planck constant");
      test(N_A*h*c/(0.119626565779*J*m/mol), name=
        "molar Planck constant times c");
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
      import FCSys.Test.assertValue;
      import FCSys.Units.*;
      extends Modelica.Icons.Function;
      output Boolean ok "true, if all tests passed";
    algorithm
      ok := false;
      // "From" functions
      assertValue(
            from_degC(100),
            373.15*U.K,
            name="from_degC");
      assertValue(
            from_kPag(101.325),
            2*U.atm,
            name="from_kPag");
      // Inverses
      assertValue(
            U.to_degC(from_degC(1)),
            1,
            name="degC");
      assertValue(
            U.to_kPag(from_kPag(1)),
            1,
            name="kPag");
      ok := true;
      annotation (Documentation(info="<html><p>This function call will fail if any of the functions return an
  incorrect result.  It will return <code>true</code> if all of the functions pass.
  There are no inputs.</p></html>"));
    end testConversions;
  end Units;

  package BaseClasses
    extends Modelica.Icons.Package;
    package Utilities
      extends Modelica.Icons.Package;
      function callAll
        "<html>Call all of the test functions for the <a href=\"modelica://FCSys.BaseClasses.Utilities\">Utilities</a> package (recursive)</html>"
        extends Modelica.Icons.Function;
        output Boolean ok "true, if all tests passed";
      algorithm
        ok := Chemistry() and Polynomial.f() and testFunctions();
        annotation (Documentation(info="<html><p>This function call will fail if any of the functions return an
  incorrect result.  It will return <code>true</code> if all of the functions pass.
  There are no inputs.</p></html>"));
      end callAll;

      function Chemistry
        "<html>Test the <a href=\"modelica://FCSys.BaseClasses.Utilities.Chemistry\">Chemistry</a> package</html>"
        import FCSys.BaseClasses.Utilities.Chemistry.*;
        extends Modelica.Icons.Function;
        output Boolean ok "true, if all tests passed";

      protected
        String strings[6];
        Integer integers[6];

      algorithm
        ok := false;

        // charge()
        for i in 1:4 loop
          assert((charge({"H2O","e-","Hg2+2",""}))[i] == {0,-1,2,0}[i],
            "The charge function failed on entry " + String(i) + ".");
        end for;
        // Note:  As of Modelica 3.2 and Dymola 7.4, assert() doesn't accept
        // vectorized equations.

        // countElements()
        for i in 1:4 loop
          assert((countElements({"H2O","H+","C19HF37O5S-",""}))[i] == {2,2,6,0}
            [i], "The countElements function failed on entry " + String(i) +
            ".");
        end for;

        // readSpecies()
        (strings[1:6],integers[1:6]) := readSpecies("C19HF37O5S-");
        for i in 1:6 loop
          assert(strings[i] == {"C","H","F","O","S","e-"}[i],
            "The elements function failed on element name of entry " + String(i)
             + ".");
          assert(integers[i] == {19,1,37,5,1,1}[i],
            "The elements function failed on element stoichiometric coefficient of entry "
             + String(i) + ".");
        end for;

        // readElement()
        (strings[1],integers[1],integers[2],strings[2]) := readElement("H2");
        assert(strings[1] == "H",
          "The readElement function failed on the element output.");
        assert(integers[1] == 2,
          "The readElement function failed on the n output.");
        assert(integers[2] == 0,
          "The readElement function failed on the z output.");
        assert(strings[2] == "",
          "The readElement function failed on the remainder output.");
        (strings[1],integers[1],integers[2],strings[2]) := readElement("Hg2+2");
        assert(strings[1] == "Hg",
          "The readElement function failed on the element output.");
        assert(integers[1] == 2,
          "The readElement function failed on the n output.");
        assert(integers[2] == 2,
          "The readElement function failed on the z output.");
        assert(strings[2] == "",
          "The readElement function failed on the remainder output.");

        // stoich()
        for i in 1:3 loop
          assert((stoich({"e-","H+","H2"}))[i] == {-2,-2,1}[i],
            "The stoich function failed on entry " + String(i) + ".");
        end for;
        for i in 1:4 loop
          assert((stoich({"e-","H+","O2","H2O"}))[i] == {-4,-4,-1,2}[i],
            "The stoich function failed on entry " + String(i) + ".");
        end for;
        ok := true;
        annotation (Documentation(info="<html><p>This function call will fail if any of the functions return an
  incorrect result.  It will return <code>true</code> if all of the functions pass.
  There are no inputs.</p></html>"));
      end Chemistry;

      package Polynomial
        extends Modelica.Icons.Package;
        import FCSys.BaseClasses.Utilities.Polynomial.*;
        model RunAll
          "<html>Run all of the test models for the <a href=\"modelica://FCSys.BaseClasses.Utilities.Polynomial\">Polynomial</a> package</html>"
          extends Modelica.Icons.Example;
          constant Boolean x[:]={f()} "Function tests";
          F testF;
          dF testdF;
          df testdf;
          d2f testd2f;
          Translatef translatef;
          annotation (Documentation(info="<html><p>If this model simulates without failure,
  then the test has passed.</p></html>"), experiment(Tolerance=1e-8));
        end RunAll;

        model Translatef
          "<html>Evaluate the translated version of <a href=\"modelica://FCSys.BaseClasses.Utilities.Polynomial.f\">f</a>()</html>"
          import FCSys.BaseClasses.Utilities.Polynomial.f;
          extends Modelica.Icons.Example;
          output Real x1=f(
                      time,
                      {1,1,1,1},
                      0);
          // Manually check the translated model to be sure that the polynomial is
          // written in nested form (for efficiency).  In Dymola 7.4 turn on the
          // option "Generate listing of translated Modelica code in dsmodel.mof".
          // dsmodel.mof should contain:
          //     x1 := 1 + time*(1 + time*(1 + time)).
          output Real x2=f(time, ones(31));
          output Real x3=f(time, ones(32));
          // The function is only unrolled to a limited depth (currently 10th order
          // polynomial).  In Dymola 7.4 f(time, ones(31)) is implemented fully
          // recursively, but f(time, ones(32)) isn't.
          output Real x4=f(
                      time + 1,
                      {1,1,1,1},
                      -3);
          // Note:  The an offset must be applied to time to prevent division
          // by zero.
        end Translatef;

        model F
          "<html>Test <a href=\"modelica://FCSys.BaseClasses.Utilities.Polynomial.F\">F</a>() based on its relation to <a href=\"modelica://FCSys.BaseClasses.Utilities.Polynomial.f\">f</a>()</html>"
          // This approach is based on [Dassault2010, vol. 2, pp. 300-301].
          import FCSys.BaseClasses.Utilities.Polynomial.*;
          extends Modelica.Icons.Example;
          parameter Integer n=-1 "Power of the first polynomial term";
          Real u1=1 + time
            "Real arguments to function (must have sufficient richness)";
          parameter Real u2[:]=1:3
            "Real arguments to function (must have sufficient richness)";
          // u2 must not be time-varying.  Otherwise, there's no requirement
          // that y1 == y2.
          Real y1 "Direct result of function";
          Real y2 "Integral of derivative of y1";
        initial equation
          y2 = y1;
        equation
          y1 = F(   u1,
                    u2,
                    n);
          f(        u1,
                    u2,
                    n) = der(y2);
          assert(abs(y1 - y2) < 1e-6, "The derivative is incorrect.");
          // Note:  The simulation tolerance is set to 1e-8.
          annotation (Documentation(info="<html><p>If this model simulates without failure,
    then the test has passed.</p></html>"), experiment(Tolerance=1e-8));
        end F;

        model dF
          "<html>Test <a href=\"modelica://FCSys.BaseClasses.Utilities.Polynomial.dF\">dF</a>() based on its relation to <a href=\"modelica://FCSys.BaseClasses.Utilities.Polynomial.F\">F</a>()</html>"
          // This approach is based on [Dassault2010, vol. 2, pp. 300-301].
          import FCSys.BaseClasses.Utilities.Polynomial.*;
          extends Modelica.Icons.Example;
          parameter Integer n=-1 "Power of the first polynomial term";
          Real u1=1 + time
            "Real arguments to function (must have sufficient richness)";
          Real u2[:]=(1 + time^2)*(1:3)
            "Real arguments to function (must have sufficient richness)";
          Real y1 "Direct result of function";
          Real y2 "Integral of derivative of y1";
        initial equation
          y2 = y1;
        equation
          y1 = F(   u1,
                    u2,
                    n);
          dF(       u1,
                    u2,
                    n,
                    der(u1),
                    der(u2)) = der(y2);
          // Note:  This is equivalent to der(y1) = der(y2), but it must be
          // explicit to ensure that the translator uses the defined derivative
          // instead of the automatically derived one.
          assert(abs(y1 - y2) < 1e-6, "The derivative is incorrect.");
          // Note:  The simulation tolerance is set to 1e-8.
          annotation (Documentation(info="<html><p>If this model simulates without failure,
  then the test has passed.</p></html>"), experiment(Tolerance=1e-8));
        end dF;

        function f
          "<html>Test <a href=\"modelica://FCSys.BaseClasses.Utilities.Polynomial.f\">f</a>()</html>"
          import FCSys.BaseClasses.Utilities.Polynomial.f;
          extends Modelica.Icons.Function;
          output Boolean ok "true, if all tests passed";
        algorithm
          ok := false;
          assert(f( 2,
                    {1,2,1},
                    0) == 1 + 2*2 + 1*2^2, "The f function failed.");
          assert(f(2, zeros(0)) == 0, "The f function failed.");
          assert(f(2, {1}) == 1, "The f function failed.");
          assert(f(2, {0,0,1}) == 4, "The f function failed.");
          assert(f(2, ones(8)) == 2^8 - 1, "The f function failed.");
          assert(f( 2,
                    {1,0,0},
                    -3) == 1/8, "The f function failed.");
          // Note:  F(), dF(), df(), and d2f() are not tested here.  They can be
          // tested by simulating TestF, TestdF, Testdf, and Testd2f.
          ok := true;
          annotation (Documentation(info="<html><p>This function call will fail if any of the functions return an
  incorrect result.  It will return <code>true</code> if all of the functions pass.
  There are no inputs.</p></html>"));
        end f;

        model df
          "<html>Test <a href=\"modelica://FCSys.BaseClasses.Utilities.Polynomial.df\">df</a>() based on its relation to <a href=\"modelica://FCSys.BaseClasses.Utilities.Polynomial.f\">f</a>()</html>"
          // This approach is based on [Dassault2010, vol. 2, pp. 300-301].
          import FCSys.BaseClasses.Utilities.Polynomial.*;
          extends Modelica.Icons.Example;
          parameter Integer n=-1 "Power of the first polynomial term";
          Real u1=1 + time
            "Real arguments to function (must have sufficient richness)";
          Real u2[:]=(1 + time^2)*(1:3)
            "Real arguments to function (must have sufficient richness)";
          Real y1 "Direct result of function";
          Real y2 "Integral of derivative of y1";
        initial equation
          y2 = y1;
        equation
          y1 = f(   u1,
                    u2,
                    n);
          df(       u1,
                    u2,
                    n,
                    der(u1),
                    der(u2)) = der(y2);
          // Note:  This is equivalent to der(y1) = der(y2), but it must be
          // explicit to ensure that the translator uses the defined derivative
          // instead of the automatically derived one.
          assert(abs(y1 - y2) < 1e-6, "The derivative is incorrect.");
          // Note:  The simulation tolerance is set to 1e-8.
          annotation (Documentation(info="<html><p>If this model simulates without failure,
  then the test has passed.</p></html>"), experiment(Tolerance=1e-8));
        end df;

        model d2f
          "<html>Test <a href=\"modelica://FCSys.BaseClasses.Utilities.Polynomial.d2f\">d2f</a>() based on its relation to <a href=\"modelica://FCSys.BaseClasses.Utilities.Polynomial.df\">df</a>()</html>"
          // This approach is based on [Dassault2010, vol. 2, pp. 300-301].
          import FCSys.BaseClasses.Utilities.Polynomial.*;
          extends Modelica.Icons.Example;
          parameter Integer n=-1 "Power of the first polynomial term";
          Real u1=1 + time
            "Real arguments to function (must have sufficient richness)";
          Real u2[:]=(1 + time^2)*(1:3)
            "Real arguments to function (must have sufficient richness)";
          Real y1 "Direct result of function";
          Real y2 "Integral of derivative of y1";
        protected
          final Real du1=der(u1) "Derivative of u1";
          final Real du2[:]=der(u2) "Derivative of u2";
          // In Dymola 7.4, it's necessary to explicitly define these intermediate
          // variables (since there are second-order derivatives).
        initial equation
          y2 = y1;
        equation
          y1 = df(  u1,
                    u2,
                    n,
                    du1,
                    du2);
          d2f(      u1,
                    u2,
                    n,
                    du1,
                    du2,
                    der(du1),
                    der(du2)) = der(y2);
          // Note:  This is equivalent to der(y1) = der(y2), but it must be
          // explicit to ensure that the translator uses the defined derivative
          // instead of the automatically derived one.
          assert(abs(y1 - y2) < 1e-6, "The derivative is incorrect.");
          // Note:  The simulation tolerance is set to 1e-8.
          annotation (Documentation(info="<html><p>If this model simulates without failure,
  then the test has passed.</p></html>"), experiment(Tolerance=1e-8));
        end d2f;
      end Polynomial;

      function testFunctions
        "<html>Test the functions in the <a href=\"modelica://FCSys.BaseClasses.Utilities\">Utilities</a> package (non-recursive)</html>"
        import FCSys.BaseClasses.Utilities.*;
        extends Modelica.Icons.Function;
        output Boolean ok "true, if all tests passed";
      protected
        String strings[6];
        Integer integers[6];
      algorithm
        ok := false;

        // average()
        assert(average({1,2,3}) == 2, "The average function failed.");

        // cartWrap()
        assert(cartWrap(0) == 3, "The cartWrap function failed.");
        assert(cartWrap(4) == 1, "The cartWrap function failed.");

        // countTrue()
        assert(countTrue({true,false,true}) == 2,
          "The countTrue function failed.");

        // Delta()
        assert(Delta({1,2}) == 1, "The Delta function failed.");
        for i in 1:2 loop
          assert((Delta([1, 2; 3, 4]))[i] == 1,
            "The Delta function failed on entry " + String(i) + ".");
        end for;

        // enumerate()
        for i in 1:3 loop
          assert((enumerate({true,false,true}))[i] == {1,0,2}[i],
            "The enumerate function failed on entry " + String(i) + ".");
        end for;

        // index()
        for i in 1:2 loop
          assert((index({true,false,true}))[i] == {1,3}[i],
            "The index function failed on entry " + String(i) + ".");
        end for;

        // inSign()
        assert(inSign(FCSys.BaseClasses.Side.n) == 1,
          "The inSign function failed.");
        assert(inSign(FCSys.BaseClasses.Side.p) == -1,
          "The inSign function failed.");

        // mod1()
        assert(mod1(4, 3) == 1, "The mod1 function failed.");
        assert(mod1(3, 3) == 3, "The mod1 function failed.");
        // Compare mod1() to mod():
        assert(mod(3, 3) == 0, "The mod function failed.");

        // round()
        for i in 1:5 loop
          assert((round({-1.6,-0.4,1.4,1.6,5}))[i] == {-2,0,1,2,5}[i],
            "The round function failed on entry " + String(i) + ".");
        end for;

        // Sigma()
        assert(Sigma({1,2}) == 3, "The Sigma function failed.");
        for i in 1:2 loop
          assert((Sigma([1, 2; 3, 4]))[i] == {3,7}[i],
            "The Sigma function failed on entry " + String(i) + ".");
        end for;

        // Compare Sigma() to sum():
        assert(sum([1, 2; 3, 4]) == 10, "The sum function failed.");
        ok := true;
        annotation (Documentation(info="<html><p>This function call will fail if any of the functions return an
  incorrect result.  It will return <code>true</code> if all of the functions pass.
  There are no inputs.</p></html>"));
      end testFunctions;
    end Utilities;
  end BaseClasses;
  annotation (Documentation(info="<html>
<p>This package may be safely removed from the
<a href=\"modelica://FCSys\">FCSys</a> distribution (along with resources/NIST.xls), but it may be helpful for debugging.
The structure  of the subpackages matches that of <a href=\"modelica://FCSys\">FCSys</a>,
although not all packages are represented.</p>

<p>The <a href=\"modelica://FCSys.Tests.CallAll\">CallAll</a>() function calls of the test
functions and the <a href=\"modelica://FCSys.Tests.RunAll\">RunAll</a> model includes all of the
test models.  Both should be tested to verify the <a href=\"modelica://FCSys\">FCSys</a>
package.</p>
</html>"));
end Tests;

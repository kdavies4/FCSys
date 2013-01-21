within FCSys;
package Test
  "<html>Models and functions to validate <a href=\"modelica://FCSys\">FCSys</a></html>"
  model TestAll "Test all the functions and models in this package (recursive)"
    extends Modelica.Icons.Example;

    Characteristics.TestAll CharacteristicsTestAll;
    BaseClasses.Utilities.TestAll BaseClassesUtilitiesTestAll;
    annotation (Documentation(info="<html><p>If this model simulates without failure,
  then the test has passed.</p></html>"), experiment(Tolerance=1e-8));
    // Note:  The simulation tolerance is set to 1e-8.
  end TestAll;
  extends Modelica.Icons.Package;
  package Characteristics
    extends Modelica.Icons.MaterialPropertiesPackage;

    model TestAll "Run all of the tests on this package"
      extends Modelica.Icons.Example;

      Testc_p testc_p;
      Testc_V testc_V;
      Testdp testdp;
      Testdv testdv;
      Testh testh;
      Testp testp;
      TestProperties testProperties;

      // Note:  The simulation tolerance is set to 1e-8.
      annotation (Documentation(info="<html><p>If this model simulates without failure,
  then the test has passed.</p></html>"), experiment(Tolerance=1e-8));
    end TestAll;

    model Testc_p
      "<html>Verify that <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.c_p\">c_p</a>() is the isobaric partial derivative of <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.h\">h</a>() with respect to temperature</html>"
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
      Q.Potential h "Specific enthalpy";
      Q.Potential y "Integral of c_p*dT";

    initial equation
      y = h;

    equation
      h = Data.h(T, p);
      der(y) = Data.c_p(T, p)*der(T);

      assert(abs(h - y) < 1e-7*U.V, "The relationship is incorrect.");
      // Note:  The simulation tolerance is set to 1e-8.
      annotation (Documentation(info="<html><p>If this model simulates without failure,
    then the test has passed.</p></html>"), experiment(Tolerance=1e-8));
    end Testc_p;

    model Testc_V
      "<html>Verify that <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.c_V\">c_V</a>() is the isobaric partial derivative of <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.h\">h</a>() - p&middot;v with respect to temperature</html>"
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
      Q.Pressure p "Pressure";

      // Results of functions
      Q.Potential u "Internal potential";
      Q.Potential y "Integral of c_V*dT";

    initial equation
      y = u;

    equation
      p = Data.p_Tv(T, v);
      u = Data.h(T, p) - p*v;
      der(y) = Data.c_V(T, p)*der(T);

      assert(abs(u - y) < 1e-5*U.V, "The relationship is incorrect.");
      // Note:  The simulation tolerance is set to 1e-8.
      annotation (Documentation(info="<html><p>If this model simulates without failure,
    then the test has passed.</p></html>"), experiment(Tolerance=1e-8));
    end Testc_V;

    model Testdp
      "<html>Verify that <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.dp\">dp</a>() is the derivative of <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.p_Tv\">p_Tv</a>()</html>"
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
      der(y2) = Data.dp(
            T,
            v,
            der(T),
            der(v));
      // Note:  This is equivalent to der(y1) = der(y2), but it must be
      // explicit to ensure that the translator uses the defined derivative
      // instead of the automatically derived one.

      assert(abs(y1 - y2) < 1e-7*U.atm, "The derivative is incorrect.");
      // Note:  The simulation tolerance is set to 1e-8.
      annotation (Documentation(info="<html><p>If this model simulates without failure,
    then the test has passed.</p></html>"), experiment(Tolerance=1e-8));
    end Testdp;

    model Testdv
      "<html>Verify that <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.dv\">dv</a>() is the derivative of <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.v_Tp\">v_Tp</a>()</html>"
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
      der(y2) = Data.dv(
            T,
            p,
            der(T),
            der(p));
      // Note:  This is equivalent to der(y1) = der(y2), but it must be
      // explicit to ensure that the translator uses the defined derivative
      // instead of the automatically derived one.

      assert(abs(y1 - y2) < 1e-6*(300*U.K/U.atm),
        "The derivative is incorrect.");
      // Note:  The simulation tolerance is set to 1e-8.

      annotation (Documentation(info="<html><p>If this model simulates without failure,
    then the test has passed.</p></html>"), experiment(Tolerance=1e-8));
    end Testdv;

    model Testh
      "<html>Verify that d<a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.h\">h</a>() = T&middot;d<a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.s\">s</a>() + v&middot;d<a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.p\">p</a>()</html>"
      // This approach is based on [Dassault2010, vol. 2, pp. 300-301].

      extends Modelica.Icons.Example;

      package Data = FCSys.Characteristics.H2O.Gas;
      // The data choice is arbitrary but the b_c values must have sufficient
      // richness.

      Q.TemperatureAbsolute T=(300 + 100*time)*U.K "Temperature";
      Q.PressureAbsolute p=U.atm "Pressure";
      Q.Potential y1 "Direct result of function";
      Q.Potential y2 "Integral of derivative of y1";

    initial equation
      y2 = y1;

    equation
      y1 = Data.h(T, p);
      // **update:
      der(y2) = Data.c_V(T, p)*der(T);
      // Note:  This is equivalent to der(y1) = der(y2), but it must be
      // explicit to ensure that the translator uses the defined derivative
      // instead of the automatically derived one.

      // **assert(abs(y1 - y2) < 1e-4*U.V, "The relationship is incorrect.");
      // Note:  The simulation tolerance is set to 1e-8.

      annotation (Documentation(info="<html><p>If this model simulates without failure,
    then the test has passed.</p></html>"), experiment(Tolerance=1e-8));
    end Testh;

    model Testp
      "<html>Verify that <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.p_Tv\">p_Tv</a>() is the inverse of <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.v_Tp\">v_Tp</a>() with respect to specific volume and pressure</html>"

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
      y = Data.p_Tv(T, Data.v_Tp(T, p));

      assert(abs(p - y) < 1e-5*U.atm, "The relationship is incorrect.");
      // Note:  The simulation tolerance is set to 1e-8.
      annotation (Documentation(info="<html><p>If this model simulates without failure,
    then the test has passed.</p></html>"), experiment(Tolerance=1e-8));
    end Testp;

    model TestProperties
      "<html>Test the <code>z</code>, <code>isCompressible</code>, <code>hasThermalExpansion</code> properties of the <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic\">Characteristic</a> package</html>"
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
      assert(not H2O.Liquid.isCompressible, "isCompressible failed on test 2.");

      // hasThermalExpansion
      assert(H2O.Gas.hasThermalExpansion,
        "hasThermalExpansion failed on test 1.");
      assert(not H2O.Liquid.hasThermalExpansion,
        "hasThermalExpansion failed on test 2.");
      annotation (Documentation(info="<html><p>If this model simulates without failure,
    then the test has passed.</p></html>"));
    end TestProperties;
  end Characteristics;

  package BaseClasses
    extends Modelica.Icons.BasesPackage;
    package Utilities
      extends Modelica.Icons.Package;

      model TestAll "Run all of the tests on this package (recursive)"
        extends Modelica.Icons.Example;

        Polynomial.TestAll PolynomialTestAll;

      initial equation
        assert(Chemistry.TestFunctions(), "The Chemistry subpackage failed.");
        assert(TestFunctions(), "TestFunctions failed.");

        // Note:  The simulation tolerance is set to 1e-8.

        annotation (Documentation(info="<html><p>If this model simulates without failure,
  then the test has passed.</p></html>"), experiment(Tolerance=1e-8));
      end TestAll;

      function TestFunctions
        "Test the functions in this package (non-recursive)"
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
        assert(Delta({1,2}) == -1, "The Delta function failed.");
        for i in 1:2 loop
          assert((Delta([1, 2; 3, 4]))[i] == -1,
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
        annotation (Documentation(info="<html><p>
  This function call will fail if any of the functions return an
  incorrect result.  It will return <code>true</code> if all of the functions pass.
  There are no inputs.</p></html>"));
      end TestFunctions;

      package Chemistry
        extends Modelica.Icons.Package;
        function TestFunctions "Test the functions in this package"

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
            assert((countElements({"H2O","H+","C19HF37O5S-",""}))[i] == {2,2,6,
              0}[i], "The countElements function failed on entry " + String(i)
               + ".");
          end for;

          // elements()
          (strings[1:6],integers[1:6]) :=
            FCSys.BaseClasses.Utilities.Chemistry.elements("C19HF37O5S-");
          for i in 1:6 loop
            assert(strings[i] == {"C","H","F","O","S","e-"}[i],
              "The elements function failed on entry " + String(i) + ".");
            assert(integers[i] == {19,1,37,5,1,1}[i],
              "The elements function failed on entry " + String(i) + ".");
          end for;

          // readElement()
          (strings[1],integers[1],integers[2],integers[3]) :=
            FCSys.BaseClasses.Utilities.Chemistry.readElement("H2");
          assert(strings[1] == "H",
            "The readElement function failed on the element output.");
          assert(integers[1] == 2,
            "The readElement function failed on the coeff output.");
          assert(integers[2] == 0,
            "The readElement function failed on the z output.");
          assert(integers[3] == 3,
            "The readElement function failed on the nextindex output.");
          (strings[1],integers[1],integers[2],integers[3]) :=
            FCSys.BaseClasses.Utilities.Chemistry.readElement("Hg2+2");
          assert(strings[1] == "Hg",
            "The readElement function failed on the element output.");
          assert(integers[1] == 2,
            "The readElement function failed on the coeff output.");
          assert(integers[2] == 2,
            "The readElement function failed on the z output.");
          assert(integers[3] == 6,
            "The readElement function failed on the nextindex output.");

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
          annotation (Documentation(info="<html><p>
  This function call will fail if any of the functions return an
  incorrect result.  It will return <code>true</code> if all of the functions pass.
  There are no inputs.</p></html>"));
        end TestFunctions;
      end Chemistry;

      package Polynomial
        extends Modelica.Icons.Package;

        import FCSys.BaseClasses.Utilities.Polynomial.*;

        model TestAll "Run all of the tests on this package"
          extends Modelica.Icons.Example;

          TestF testF;
          TestdF testdF;
          Testdf testdf;
          Testd2f testd2f;
          Translatef translatef;

        initial equation
          assert(Testf(), "Testf failed.");

          // Note:  The simulation tolerance is set to 1e-8.

          annotation (Documentation(info="<html><p>If this model simulates without failure,
  then the test has passed.</p></html>"), experiment(Tolerance=1e-8));
        end TestAll;

        model TestF
          "<html>Verify that <a href=\"modelica://FCSys.BaseClasses.Utilities.Polynomial.F\">F</a>() is correctly related to <a href=\"modelica://FCSys.BaseClasses.Utilities.Polynomial.f\">f</a>()</html>"
          // This approach is based on [Dassault2010, vol. 2, pp. 300-301].

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
        end TestF;

        model TestdF
          "<html>Verify that <a href=\"modelica://FCSys.BaseClasses.Utilities.Polynomial.dF\">dF</a>() is the correct derivative of <a href=\"modelica://FCSys.BaseClasses.Utilities.Polynomial.F\">F</a>()</html>"
          // This approach is based on [Dassault2010, vol. 2, pp. 300-301].

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
        end TestdF;

        model Testdf
          "<html>Verify that <a href=\"modelica://FCSys.BaseClasses.Utilities.Polynomial.df\">df</a>() is the correct derivative of <a href=\"modelica://FCSys.BaseClasses.Utilities.Polynomial.f\">f</a>()</html>"
          // This approach is based on [Dassault2010, vol. 2, pp. 300-301].

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
        end Testdf;

        model Testd2f
          "<html>Verify that <a href=\"modelica://FCSys.BaseClasses.Utilities.Polynomial.d2f\">d2f</a>() is the correct derivative of <a href=\"modelica://FCSys.BaseClasses.Utilities.Polynomial.df\">df</a>()</html>"
          // This approach is based on [Dassault2010, vol. 2, pp. 300-301].

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
        end Testd2f;

        model Translatef
          "<html>Evaluate the translated version of <a href=\"modelica://FCSys.BaseClasses.Utilities.Polynomial.f\">f</a>()</html>"
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
          // Note:  The an offest must be applied to time to prevent division
          // by zero.
        end Translatef;

        function Testf
          "<html>Test <a href=\"modelica://FCSys.BaseClasses.Utilities.Polynomial.f\">f</a>()</html>"

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
          annotation (Documentation(info="<html><p>
  This function call will fail if any of the functions return an
  incorrect result.  It will return <code>true</code> if all of the functions pass.
  There are no inputs.</p></html>"));
        end Testf;
      end Polynomial;
    end Utilities;
  end BaseClasses;
  annotation (Documentation(info="<html>
<p>This package contains tests that are not useful for most users of
<a href=\"modelica://FCSys\">FCSys</a>.  It can be safely removed from the
distribution.  The structure of the subpackages
matches that of <a href=\"modelica://FCSys\">FCSys</a>, although not all
packages are represented.
</p>
</html>"));
  package Units
    extends Modelica.Icons.Package;
    model TestUnits "Test the unit relations"
      import FCSys.Units.*;
      extends Modelica.Icons.Example;

    protected
      function Test "Test a single value or ratio"

        input Real x "Actual value";
        input Real tolerance=1e-7 "Error tolerance";
        input String name "Name of test";

      algorithm
        assert(max(abs(x - 1)) < tolerance, "Test " + name +
          " is incorrect.\nIt must be approximately 1.0 (tolerance " + String(
          tolerance) + ") but is " + String(x) + ".");
      end Test;

    initial equation
      // Set 1:  Mathematical constants and relations
      Test(pi/3.14159265358979323846264338327950288419716939937510, name=
        "1 in set 1");
      // Value from http://en.wikipedia.org/wiki/Pi#Approximate_value
      Test(e/2.71828182845904523536028747135266249775724709369995, name=
        "2 in set 1");
      // Value from http://en.wikipedia.org/wiki/E_(mathematical_constant)
      Test(2*pi*rad/(360*degree), name="3 in set 1");
      Test('%'/0.01, name="4 in set 1");

      // Set 2:  Relations from BIPM (2006)
      Test(q/(1.602176487E-19*C), name="1 in set 2");
      Test(C*V/J, name="2 in set 2");

      // Set 3:  Coherent derived units in the SI with special names and
      // symbols (BIPM, 2006)
      // These must be correct to the computer's numerical accuracy.
      Test(sr, name="1 in set 3");
      Test(Hz/(cyc/s), name="2 in set 3");
      // BIPM implicitly assumes that the unit cycle is 1, but the unit rad is 1
      // too---a discrepency.
      Test(N/(kg*m/s^2), name="3 in set 3");
      Test(Pa/(N/m^2), name="4 in set 3");
      Test(J/(N*m), name="5 in set 3");
      Test(W/(J/s), name="6 in set 3");
      Test(C/(A*s), name="7 in set 3");
      Test(V/(W/A), name="8 in set 3");
      Test(F/(C/V), name="9 in set 3");
      Test(ohm/(V/A), name="10 in set 3");
      Test(S/(A/V), name="11 in set 3");
      Test(Wb/(V*s), name="12 in set 3");
      Test(T/(Wb/m^2), name="13 in set 3");
      Test(H/(Wb/A), name="14 in set 3");
      Test(lm/('cd'*sr), name="15 in set 3");
      Test(lx/(lm/m^2), name="16 in set 3");
      Test(Bq/(cyc/s), name="17 in set 3");
      // See Hz.
      Test(Gy/(J/kg), name="18 in set 3");
      Test(Sv/(J/kg), name="19 in set 3");
      Test(kat/(mol/s), name="20 in set 3");

      // Set 4:  Relations from NIST (2010)
      // This has been generated from FCSys/misc/NIST.xls.
      Test(1/alpha/137.035999074, name="inverse fine-structure constant");
      Test(1/G_0/(12906.4037217*ohm), name="inverse of conductance quantum");
      Test(1/m*h*c/(1.239841930E-6*eV), name=
        "inverse meter-electron volt relationship");
      Test(1/m*h*c/(4.556335252755E-8*E_h), name=
        "inverse meter-hartree relationship");
      Test(1/m*h*c/(1.986445684E-25*J), name="inverse meter-joule relationship");
      Test(1/m*h*c/k_B/(1.4387770E-2*K), name=
        "inverse meter-kelvin relationship");
      Test(1/m*h/c/(2.210218902E-42*kg), name=
        "inverse meter-kilogram relationship");
      Test(100*kPa/(k_B*273.15*K)/(2.6516462E25/m^3), name=
        "Loschmidt constant (273.15 K, 100 kPa)");
      Test(101.325*kPa/(k_B*273.15*K)/(2.6867805E25/m^3), name=
        "Loschmidt constant (273.15 K, 101.325 kPa)");
      Test(12*g/mol/(12E-3*kg/mol), name="molar mass of carbon-12");
      Test(alpha/7.2973525698E-3, name="fine-structure constant");
      Test(atm/(101325*Pa), name="standard atmosphere");
      Test(c/(299792458*m/s), name="speed of light in vacuum");
      Test(c_1/(3.74177153E-16*W*m^2), name="first radiation constant");
      Test(c_2/(1.4387770E-2*m*K), name="second radiation constant");
      Test(c_3_f*cyc/(5.8789254E10*Hz/K), name=
        "Wien frequency displacement law constant");
      // Factor of cyc due to inconsistencies b/w rad and cyc in BIPM (2006)
      Test(c_3_lambda/(2.8977721E-3*m*K), name=
        "Wien wavelength displacement law constant");
      Test(cyc/m*c/(299792458*Hz), name="inverse meter-hertz relationship");
      Test(E_h/(4.35974434E-18*J), name="Hartree energy");
      Test(E_h/(27.21138505*eV), name="Hartree energy in eV");
      Test(E_h/(27.21138505*eV), name="hartree-electron volt relationship");
      Test(E_h/(4.35974434E-18*J), name="hartree-joule relationship");
      Test(E_h/(h*c)/(2.194746313708E7/m), name=
        "hartree-inverse meter relationship");
      Test(E_h/c^2/(4.85086979E-35*kg), name="hartree-kilogram relationship");
      Test(E_h/h*cyc/(6.579683920729E15*Hz), name="hartree-hertz relationship");
      // Factor of cyc due to inconsistencies b/w rad and cyc in BIPM (2006)
      Test(E_h/k_B/(3.1577504E5*K), name="hartree-kelvin relationship");
      Test(epsilon_0/(8.854187817E-12*F/m), name="electric constant");
      Test(eV/(1.602176565E-19*J), name="electron volt");
      Test(eV/(3.674932379E-2*E_h), name="electron volt-hartree relationship");
      Test(eV/(1.602176565E-19*J), name="electron volt-joule relationship");
      Test(eV/(h*c)/(8.06554429E5/m), name=
        "electron volt-inverse meter relationship");
      Test(eV/c^2/(1.782661845E-36*kg), name=
        "electron volt-kilogram relationship");
      Test(eV/h*cyc/(2.417989348E14*Hz), name=
        "electron volt-hertz relationship");
      // Factor of cyc due to inconsistencies b/w rad and cyc in BIPM (2006)
      Test(eV/k_B/(1.1604519E4*K), name="electron volt-kelvin relationship");
      Test(G_0/(7.7480917346E-5*S), name="conductance quantum");
      Test(g/mol/(1E-3*kg/mol), name="molar mass constant");
      Test(h/(6.62606957E-34*J*s), name="Planck constant");
      Test(h/(4.135667516E-15*eV*s), name="Planck constant in eV s");
      Test(h*c/(2*pi)/(197.3269718*mega*eV*femto*m), name=
        "Planck constant over 2 pi times c in MeV fm");
      Test(h/(2*pi)/(1.054571726E-34*J*s), name="Planck constant over 2 pi");
      Test(h/(2*pi)/(6.58211928E-16*eV*s), name=
        "Planck constant over 2 pi in eV s");
      Test(Hz*h/cyc/(4.135667516E-15*eV), name=
        "hertz-electron volt relationship");
      // Factor of 1/cyc due to inconsistencies b/w rad and cyc in BIPM (2006)
      Test(Hz*h/cyc/(1.5198298460045E-16*E_h), name=
        "hertz-hartree relationship");
      // Factor of 1/cyc due to inconsistencies b/w rad and cyc in BIPM (2006)
      Test(Hz*h/cyc/(6.62606957E-34*J), name="hertz-joule relationship");
      // Factor of 1/cyc due to inconsistencies b/w rad and cyc in BIPM (2006)
      Test(Hz*h/c^2/cyc/(7.37249668E-51*kg), name="hertz-kilogram relationship");
      // Factor of 1/cyc due to inconsistencies b/w rad and cyc in BIPM (2006)
      Test(Hz*h/k_B/cyc/(4.7992434E-11*K), name="hertz-kelvin relationship");
      // Factor of 1/cyc due to inconsistencies b/w rad and cyc in BIPM (2006)
      Test(Hz/c/cyc/(3.335640951E-9/m), name="hertz-inverse meter relationship");
      // Factor of 1/cyc due to inconsistencies b/w rad and cyc in BIPM (2006)
      Test(J/(6.24150934E18*eV), name="joule-electron volt relationship");
      Test(J/(2.29371248E17*E_h), name="joule-hartree relationship");
      Test(J/(h*c)/(5.03411701E24/m), name="joule-inverse meter relationship");
      Test(J/c^2/(1.112650056E-17*kg), name="joule-kilogram relationship");
      Test(J/h*cyc/(1.509190311E33*Hz), name="joule-hertz relationship");
      Test(J/k_B/(7.2429716E22*K), name="joule-kelvin relationship");
      Test(k_B/(1.3806488E-23*J/K), name="Boltzmann constant");
      Test(k_B/(8.6173324E-5*eV/K), name="Boltzmann constant in eV/K");
      Test(k_B/(h*c)/(69.503476/(m*K)), name=
        "Boltzmann constant in inverse meters per kelvin");
      Test(k_B/h*cyc/(2.0836618E10*Hz/K), name="Boltzmann constant in Hz/K");
      // Factor of cyc due to inconsistencies b/w rad and cyc in BIPM (2006)
      Test(k_F/(96485.3365*C/mol), name="Faraday constant");
      Test(k_J*cyc/(483597.9E9*Hz/V), name=
        "conventional value of Josephson constant");
      Test(k_J*cyc/(483597.870E9*Hz/V), name="Josephson constant");
      // Factor of cyc due to inconsistencies b/w rad and cyc in BIPM (2006)
      Test(K*k_B/(8.6173324E-5*eV), name="kelvin-electron volt relationship");
      Test(K*k_B/(3.1668114E-6*E_h), name="kelvin-hartree relationship");
      Test(K*k_B/(1.3806488E-23*J), name="kelvin-joule relationship");
      Test(K*k_B/(h*c)/(69.503476/m), name="kelvin-inverse meter relationship");
      Test(K*k_B/c^2/(1.5361790E-40*kg), name="kelvin-kilogram relationship");
      Test(K*k_B/h*cyc/(2.0836618E10*Hz), name="kelvin-hertz relationship");
      // Factor of cyc due to inconsistencies b/w rad and cyc in BIPM (2006)
      Test(kg*c^2/(5.60958885E35*eV), name=
        "kilogram-electron volt relationship");
      Test(kg*c^2/(2.061485968E34*E_h), name="kilogram-hartree relationship");
      Test(kg*c^2/(8.987551787E16*J), name="kilogram-joule relationship");
      Test(kg*c^2/h*cyc/(1.356392608E50*Hz), name="kilogram-hertz relationship");
      // Factor of cyc due to inconsistencies b/w rad and cyc in BIPM (2006)
      Test(kg*c^2/k_B/(6.5096582E39*K), name="kilogram-kelvin relationship");
      Test(kg*c/h/(4.52443873E41/m), name="kilogram-inverse meter relationship");
      Test(mu_0/(12.566370614E-7*N/A^2), name="mag. constant");
      Test(N_A/(6.02214129E23/mol), name="Avogadro constant");
      Test(N_A*h/(3.9903127176E-10*J*s/mol), name="molar Planck constant");
      Test(N_A*h*c/(0.119626565779*J*m/mol), name=
        "molar Planck constant times c");
      Test(Phi_0/(2.067833758E-15*Wb), name="mag. flux quantum");
      Test(q/(1.602176565E-19*C), name="elementary charge");
      Test(q/h/(2.417989348E14*A/J), name="elementary charge over h");
      Test(R/(8.3144621*J/(mol*K)), name="molar gas constant");
      Test(R_inf/(10973731.568539/m), name="Rydberg constant");
      Test(R_inf*c*cyc/(3.289841960364E15*Hz), name=
        "Rydberg constant times c in Hz");
      // Factor of cyc due to inconsistencies b/w rad and cyc in BIPM (2006)
      Test(R_inf*h*c/(13.60569253*eV), name="Rydberg constant times hc in eV");
      Test(R_inf*h*c/(2.179872171E-18*J), name="Rydberg constant times hc in J");
      Test(R_K/(25812.807*ohm), name=
        "conventional value of von Klitzing constant");
      Test(R_K/(25812.8074434*ohm), name="von Klitzing constant");
      Test(R*273.15*K/(100*kPa)/(22.710953E-3*m^3/mol), name=
        "molar volume of ideal gas (273.15 K, 100 kPa)");
      Test(R*273.15*K/(101.325*kPa)/(22.413968E-3*m^3/mol), name=
        "molar volume of ideal gas (273.15 K, 101.325 kPa)");
      Test(sigma/(5.670373E-8*W/(m^2*K^4)), name="Stefan-Boltzmann constant");
      Test(Z_0/(376.730313461*ohm), name="characteristic impedance of vacuum");

      annotation (Documentation(info="<html><p>If this model simulates without failure,
  then the test has passed.</p></html>"));
    end TestUnits;
  end Units;
end Test;

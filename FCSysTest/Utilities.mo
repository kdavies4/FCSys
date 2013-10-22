within FCSysTest;
package Utilities
  extends Modelica.Icons.Package;
  function callAll
    "<html>Call all of the test functions for the <a href=\"modelica://FCSys.Utilities\">Utilities</a> package (recursive)</html>"
    import Modelica.Utilities.Streams.print;
    extends Modelica.Icons.Function;

    input String logFile="FCSysTestLog.txt" "Filename where the log is stored";
    output Boolean ok "true, if all tests passed";

  algorithm
    print("--- Test FCSys.Utilities");
    print("--- Test FCSys.Utilities", logFile);

    ok := Chemistry() and Polynomial.f() and testFunctions();
    annotation (Documentation(info="<html><p>This function call will fail if any of the functions return an
  incorrect result.  It will return <code>true</code> if all of the functions pass.
  There are no inputs.</p></html>"));
  end callAll;

  function Chemistry
    "<html>Test the <a href=\"modelica://FCSys.Utilities.Chemistry\">Chemistry</a> package</html>"
    import Modelica.Utilities.Streams.print;
    import FCSys.Utilities.Chemistry.*;
    import FCSys.Utilities.arrayIntegerEqual;
    import FCSys.Utilities.arrayStringEqual;
    extends Modelica.Icons.Function;

    input String logFile="FCSysTestLog.txt" "Filename where the log is stored";
    output Boolean ok "true, if all tests passed";

  protected
    String strings[6];
    Integer integers[6];

  algorithm
    print("... Test of Utilities.Chemistry");
    print("... Test of Utilities.Chemistry", logFile);

    // charge()
    assert(arrayIntegerEqual(charge({"H2O","e-","Hg2+2",""}), {0,-1,2,0}),
      "The charge function failed.");

    // countElements()
    assert(arrayIntegerEqual(countElements({"H2O","H+","C19HF37O5S-",""}), {2,2,
      6,0}), "The countElements function failed on entry.");

    // readSpecies()
    (strings[1:6],integers[1:6]) := readSpecies("C19HF37O5S-");
    assert(arrayStringEqual(strings[1:6], {"C","H","F","O","S","e-"}),
      "The readSpecies function failed on the element names.");
    assert(arrayIntegerEqual(integers[1:6], {19,1,37,5,1,1}),
      "The readSpecies function failed on the element stoichiometric coefficients.");

    // readElement()
    (strings[1],integers[1],integers[2],strings[2]) := readElement("H2");
    assert(strings[1] == "H",
      "The readElement function failed on the element output.");
    assert(integers[1] == 2, "The readElement function failed on the n output.");
    assert(integers[2] == 0, "The readElement function failed on the z output.");
    assert(strings[2] == "",
      "The readElement function failed on the remainder output.");
    (strings[1],integers[1],integers[2],strings[2]) := readElement("Hg2+2");
    assert(strings[1] == "Hg",
      "The readElement function failed on the element output.");
    assert(integers[1] == 2, "The readElement function failed on the n output.");
    assert(integers[2] == 2, "The readElement function failed on the z output.");
    assert(strings[2] == "",
      "The readElement function failed on the remainder output.");

    // stoich()
    assert(arrayIntegerEqual(stoich({"e-","H+","H2"}), {-2,-2,1}),
      "The stoich function failed on test 1.");
    assert(arrayIntegerEqual(stoich({"e-","H+","O2","H2O"}), {-4,-4,-1,2}),
      "The stoich function failed on test 2.");

    ok := true;
    annotation (Documentation(info="<html><p>This function call will fail if any of the functions return an
  incorrect result.  It will return <code>true</code> if all of the functions pass.
  There are no inputs.</p></html>"));
  end Chemistry;

  package Polynomial
    extends Modelica.Icons.Package;
    import FCSys.Utilities.Polynomial.*;

    model Translatef
      "<html>Evaluate the translated version of <a href=\"modelica://FCSys.Utilities.Polynomial.f\">f</a>()</html>"
      import FCSys.Utilities.Polynomial.f;
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
      "<html>Test <a href=\"modelica://FCSys.Utilities.Polynomial.F\">F</a>() based on its relation to <a href=\"modelica://FCSys.Utilities.Polynomial.f\">f</a>()</html>"
      // This approach is based on [Dassault2010, vol. 2, pp. 300-301].
      import FCSys.Utilities.Polynomial.*;
      extends Modelica.Icons.Example;

      parameter Integer n=-1 "Power of the first polynomial term"
        annotation (Dialog(__Dymola_label="<html><i>n</i></html>"));
      Real u_1=1 + time
        "Real arguments to function (must have sufficient richness)";
      parameter Real u_2[:]=1:3
        "Real arguments to function (must have sufficient richness)"
        annotation (Dialog(__Dymola_label="<html><i>u</i><sub>2</sub></html>"));
      // u_2 must not be time-varying.  Otherwise, there's no requirement
      // that y_1 == y_2.
      Real y_1 "Direct result of function";
      Real y_2 "Integral of derivative of y_1";

    initial equation
      y_2 = y_1;

    equation
      y_1 = F(
            u_1,
            u_2,
            n);
      f(    u_1,
            u_2,
            n) = der(y_2);
      assert(abs(y_1 - y_2) < 1e-6, "The derivative is incorrect.");
      // Note:  The simulation tolerance is set to 1e-8.
      annotation (Documentation(info="<html><p>If this model simulates without failure,
    then the test has passed.</p></html>"), experiment(Tolerance=1e-8));
    end F;

    model dF
      "<html>Test <a href=\"modelica://FCSys.Utilities.Polynomial.dF\">dF</a>() based on its relation to <a href=\"modelica://FCSys.Utilities.Polynomial.F\">F</a>()</html>"
      // This approach is based on [Dassault2010, vol. 2, pp. 300-301].
      import FCSys.Utilities.Polynomial.*;
      extends Modelica.Icons.Example;
      parameter Integer n=-1 "Power of the first polynomial term"
        annotation (Dialog(__Dymola_label="<html><i>n</i></html>"));
      Real u_1=1 + time
        "Real arguments to function (must have sufficient richness)";
      Real u_2[:]=(1 + time^2)*(1:3)
        "Real arguments to function (must have sufficient richness)";
      Real y_1 "Direct result of function";
      Real y_2 "Integral of derivative of y_1";

    initial equation
      y_2 = y_1;

    equation
      y_1 = F(
            u_1,
            u_2,
            n);
      dF(   u_1,
            u_2,
            n,
            der(u_1),
            der(u_2)) = der(y_2);
      // Note:  This is equivalent to der(y_1) = der(y_2), but it must be
      // explicit to ensure that the translator uses the defined derivative
      // instead of the automatically derived one.
      assert(abs(y_1 - y_2) < 1e-6, "The derivative is incorrect.");
      // Note:  The simulation tolerance is set to 1e-8.
      annotation (Documentation(info="<html><p>If this model simulates without failure,
  then the test has passed.</p></html>"), experiment(Tolerance=1e-8));
    end dF;

    function f
      "<html>Test <a href=\"modelica://FCSys.Utilities.Polynomial.f\">f</a>()</html>"
      import Modelica.Utilities.Streams.print;
      import FCSys.Utilities.Polynomial.f;
      extends Modelica.Icons.Function;

      input String logFile="FCSysTestLog.txt"
        "Filename where the log is stored";
      output Boolean ok "true, if all tests passed";

    algorithm
      print("... Test of Utilities.Polynomial.f()");
      print("... Test of Utilities.Polynomial.f()", logFile);

      assert(f(
            2,
            {1,2,1},
            0) == 1 + 2*2 + 1*2^2, "The f function failed.");
      assert(f(2, zeros(0)) == 0, "The f function failed.");
      assert(f(2, {1}) == 1, "The f function failed.");
      assert(f(2, {0,0,1}) == 4, "The f function failed.");
      assert(f(2, ones(8)) == 2^8 - 1, "The f function failed.");
      assert(f(
            2,
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
      "<html>Test <a href=\"modelica://FCSys.Utilities.Polynomial.df\">df</a>() based on its relation to <a href=\"modelica://FCSys.Utilities.Polynomial.f\">f</a>()</html>"
      // This approach is based on [Dassault2010, vol. 2, pp. 300-301].
      import FCSys.Utilities.Polynomial.*;
      extends Modelica.Icons.Example;
      parameter Integer n=-1 "Power of the first polynomial term"
        annotation (Dialog(__Dymola_label="<html><i>n</i></html>"));
      Real u_1=1 + time
        "Real arguments to function (must have sufficient richness)";
      Real u_2[:]=(1 + time^2)*(1:3)
        "Real arguments to function (must have sufficient richness)";
      Real y_1 "Direct result of function";
      Real y_2 "Integral of derivative of y_1";

    initial equation
      y_2 = y_1;

    equation
      y_1 = f(
            u_1,
            u_2,
            n);
      df(   u_1,
            u_2,
            n,
            der(u_1),
            der(u_2)) = der(y_2);
      // Note:  This is equivalent to der(y_1) = der(y_2), but it must be
      // explicit to ensure that the translator uses the defined derivative
      // instead of the automatically derived one.
      assert(abs(y_1 - y_2) < 1e-6, "The derivative is incorrect.");
      // Note:  The simulation tolerance is set to 1e-8.
      annotation (Documentation(info="<html><p>If this model simulates without failure,
  then the test has passed.</p></html>"), experiment(Tolerance=1e-8));
    end df;

    model d2f
      "<html>Test <a href=\"modelica://FCSys.Utilities.Polynomial.d2f\">d2f</a>() based on its relation to <a href=\"modelica://FCSys.Utilities.Polynomial.df\">df</a>()</html>"
      // This approach is based on [Dassault2010, vol. 2, pp. 300-301].
      import FCSys.Utilities.Polynomial.*;
      extends Modelica.Icons.Example;
      parameter Integer n=-1 "Power of the first polynomial term"
        annotation (Dialog(__Dymola_label="<html><i>n</i></html>"));
      Real u_1=1 + time
        "Real arguments to function (must have sufficient richness)";
      Real u_2[:]=(1 + time^2)*(1:3)
        "Real arguments to function (must have sufficient richness)";
      Real y_1 "Direct result of function";
      Real y_2 "Integral of derivative of y_1";

    protected
      final Real du_1=der(u_1) "Derivative of u_1";
      final Real du_2[:]=der(u_2) "Derivative of u_2";
      // In Dymola 7.4, it's necessary to explicitly define these intermediate
      // variables (since there are second-order derivatives).

    initial equation
      y_2 = y_1;

    equation
      y_1 = df(
            u_1,
            u_2,
            n,
            du_1,
            du_2);
      d2f(  u_1,
            u_2,
            n,
            du_1,
            du_2,
            der(du_1),
            der(du_2)) = der(y_2);
      // Note:  This is equivalent to der(y_1) = der(y_2), but it must be
      // explicit to ensure that the translator uses the defined derivative
      // instead of the automatically derived one.
      assert(abs(y_1 - y_2) < 1e-6, "The derivative is incorrect.");
      // Note:  The simulation tolerance is set to 1e-8.
      annotation (Documentation(info="<html><p>If this model simulates without failure,
  then the test has passed.</p></html>"), experiment(Tolerance=1e-8));
    end d2f;

  end Polynomial;

  function testFunctions
    "<html>Test the functions in the <a href=\"modelica://FCSys.Utilities\">Utilities</a> package (non-recursive)</html>"
    import FCSys.Species.Enumerations.Side;
    import FCSys.Utilities.*;
    extends Modelica.Icons.Function;

    output Boolean ok "true, if all tests passed";

  protected
    String strings[6];
    Integer integers[6];

  algorithm
    // arrayIntegerEqual()
    assert(arrayIntegerEqual({1,2}, {1,2}),
      "The arrayIntegerEqual function failed on test 1.");
    assert(not arrayIntegerEqual({1,2}, {1,3}),
      "The arrayIntegerEqual function failed on test 2.");
    assert(not arrayIntegerEqual({1,2}, {1,2,3}),
      "The arrayIntegerEqual function failed on test 3.");

    // arrayRealEqual()
    assert(arrayRealEqual({1,2}, {1,2}),
      "The arrayRealEqual function failed on test 1.");
    assert(not arrayRealEqual({1,2}, {1,2.001}),
      "The arrayRealEqual function failed on test 2.");
    assert(not arrayRealEqual({1,2}, {1,2,3}),
      "The arrayRealEqual function failed on test 3.");

    // arrayStringEqual()
    assert(arrayStringEqual({"a","bc"}, {"a","bc"}),
      "The arrayStringEqual function failed on test 1.");
    assert(not arrayStringEqual({"a","bc"}, {"a","b"}),
      "The arrayStringEqual function failed on test 2.");
    assert(not arrayStringEqual({"a","b"}, {"a","b","c"}),
      "The arrayStringEqual function failed on test 3.");

    // average()
    assert(average({1,2,3}) == 2, "The average function failed.");

    // cartWrap()
    assert(cartWrap(0) == 3, "The cartWrap function failed on test 1.");
    assert(cartWrap(4) == 1, "The cartWrap function failed on test 2.");

    // Delta()
    assert(Delta({1,2}) == 1, "The Delta function failed on test 1.");
    assert(arrayRealEqual(Delta([1, 2; 3, 4]), {1,1}),
      "The Delta function failed on test 2.");

    // inSign()
    assert(inSign(Side.n) == 1, "The inSign function failed on test 1.");
    assert(inSign(Side.p) == -1, "The inSign function failed on test 2.");

    // mod1()
    assert(mod1(4, 3) == 1, "The mod1 function failed on test 1.");
    assert(mod1(3, 3) == 3, "The mod1 function failed on test 2.");
    // Compare mod1() to mod():
    assert(mod(3, 3) == 0, "The mod function failed.");

    // round()
    assert(arrayRealEqual(round({-1.6,-0.4,1.4,1.6,5}), {-2,0,1,2,5}),
      "The round function failed on entry.");

    // Sigma()
    assert(Sigma({1,2}) == 3, "The Sigma function failed on test 1.");
    assert(arrayRealEqual(Sigma([1, 2; 3, 4]), {3,7}),
      "The Sigma function failed on test 2.");
    // Compare Sigma() to sum():
    assert(sum([1, 2; 3, 4]) == 10, "The sum function failed.");

    ok := true;
    annotation (Documentation(info="<html><p>This function call will fail if any of the functions return an
  incorrect result.  It will return <code>true</code> if all of the functions pass.
  There are no inputs.</p></html>"));
  end testFunctions;

end Utilities;

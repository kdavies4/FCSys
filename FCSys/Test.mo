within FCSys;
package Test "Library to facilitate assertion-based testing of Modelica code"
  extends Modelica.Icons.Package;
  function assertValue "Assert that a value is within specification"
    extends Modelica.Icons.Function;
    input Real actual "Actual value";
    input Real expected "Expected value";
    input Real eps=1e-7 "Error tolerance";
    input String name="" "Name of test";

  algorithm
    assert(abs(actual - expected) <= eps, (if name <> "" then "Test " + name +
      " failed.\n" else "") + "The actual value (" + String(actual) +
      ") was not within " + String(eps) + " of the expected value (" + String(
      expected) + ").");

  end assertValue;

  function assertValues "Assert that values are within specification"
    extends Modelica.Icons.Function;
    input Real actual[:] "Actual values";
    input Real expected[size(actual, 1)] "Expected values";
    input Real eps=1e-7 "Error tolerance";
    input String name="" "Name of test";

  algorithm
    for i in 1:size(actual, 1) loop
      assert(abs(actual[i] - expected[i]) <= eps, "Test " + String(i) + (if
        name <> "" then " of " + name else "") + " failed.\n" +
        "The actual value (" + String(actual[i]) + ") was not within " + String(
        eps) + " of the expected value (" + String(expected[i]) + ").")
        annotation (Inline=true);
    end for;

  end assertValues;

  function assertLogValue
    "Assert that a value is within orders of magnitude of specification"
    import Modelica.Math.log10;
    import Modelica.Constants.small;
    extends Modelica.Icons.Function;
    input Real actual "Actual value";
    input Real expected "Expected value";
    input Real o=1 "Error tolerance in orders of magnitude";
    input String name="" "Name of test";

  algorithm
    assert((abs(actual) < small and abs(expected) < small) or abs(log10(actual/
      expected)) <= o, (if name <> "" then "Test " + name + " failed.\n" else
      "") + "The actual value (" + String(actual) + ") was not within " +
      String(o) + " orders of magnitude of the expected value (" + String(
      expected) + ").") annotation (Inline=true);

  end assertLogValue;

  function assertLogValues
    "Assert that values are within orders of magnitude of specification"
    import Modelica.Math.log10;
    import Modelica.Constants.small;
    extends Modelica.Icons.Function;
    input Real actual[:] "Actual values";
    input Real expected[size(actual, 1)] "Expected values";
    input Real o=1 "Error tolerance in orders of magnitude";
    input String name="" "Name of test";

  algorithm
    for i in 1:size(actual, 1) loop
      assert((abs(actual[i]) < small and abs(expected[i]) < small) or abs(log10
        (actual[i]/expected[i])) <= o, "Test " + String(i) + (if name <> ""
         then " of " + name else "") + " failed.\n" + "The actual value (" +
        String(actual[i]) + ") was not within " + String(o) +
        " orders of magnitude of the expected value (" + String(expected[i]) +
        ").") annotation (Inline=true);
    end for;

  end assertLogValues;

  model AssertTrajectory
    parameter Real expected[:, 2];
    parameter Real eps=1e-6;
    input Real actual;

  protected
    Integer cur(start=1,fixed=true);

  algorithm
    when initial() then
      assert(size(expected, 1) > 0,
        "The expected trajectory contains no points.");
      assert(expected[1, 1] >= time,
        "Some trajectory points precede the simulation.");
      cur := 1;
      if (expected[1, 1] >= time and expected[1, 1] <= time) then
        assertValue(
            actual,
            expected[cur, 2],
            eps);
        cur := 2;
      end if;
    end when;
    when cur <= size(expected, 1) and time >= expected[cur, 1] then
      assertValue(
          actual,
          expected[cur, 2],
          eps);
      cur := pre(cur) + 1;
    end when;
    when terminal() then
      assert(cur > size(expected, 1),
        "The simulation ended before all trajectory points could be checked.");
      // Note:  In Dymola 7.4, the simulation log may state "Integration
      // terminated successfully" and then the assertion statement below it.
    end when;

  end AssertTrajectory;

  model AssertBecomesTrueAt
    parameter Modelica.SIunits.Time at;
    parameter Modelica.SIunits.Time eps=1e-6;
    input Boolean event;

  algorithm
    when initial() then
      assert(at > time + eps,
        "The expected crossing time is before the start of the simulation.");
    end when;
    when time >= at - eps then
      assert(not event, "The signal became true before expected crossing time.");
    end when;
    when time >= at + eps then
      assert(event, "The signal was not true by the expected crossing time.");
    end when;
    when terminal() then
      assert(time >= at + eps,
        "The expected crossing time is after the end of the simulation.");
      // Note:  In Dymola 7.4, the simulation log may state "Integration
      // terminated successfully" and then the assertion statement below it.
    end when;

  end AssertBecomesTrueAt;

  model AssertInitial "Assert the initial value of a signal"
    parameter Real expected;
    parameter Real eps=1e-6;
    input Real actual;

  algorithm
    when initial() then
      assertValue(
          actual,
          expected,
          eps);
    end when;

  end AssertInitial;

  model AssertFinal "Assert the final value of a signal"
    parameter Real expected;
    parameter Real eps=1e-6;
    input Real actual;

  algorithm
    when terminal() then
      assertValue(
          actual,
          expected,
          eps);
      // Note:  In Dymola 7.4, the simulation log may state "Integration
      // terminated successfully" and then the assertion statement below it.
    end when;

  end AssertFinal;

  model AssertAverageBetween "Assert an average value between two times."
    parameter Real average;
    parameter Real start;
    parameter Real finish;
    parameter Real eps=1e-6;
    input Real signal;

  protected
    Real integral;

  initial equation
    integral = 0;
    assert(finish > start,
      "The end of interval must be after start of interval.");
    assert(time <= start, "The simulation started after the interval.");

  equation
    der(integral) = if (time < start) then 0 else signal;
    when time >= finish then
      assert(abs((integral/(finish - start)) - average) < eps,
        "The average value between times " + String(start) + " and " + String(
        finish) + " was " + String(average) + " but should have been within "
         + String(eps) + " of " + String(integral/(finish - start)) + ".");
    end when;

  algorithm
    when terminal() then
      assert(time >= finish,
        "The simulation terminated before the interval was completed.");
      // Note:  In Dymola 7.4, the simulation log may state "Integration
      // terminated successfully" and then the assertion statement below it.
    end when;

  end AssertAverageBetween;

  model AssertValueAt "Assert the initial value of a signal at a specific time"
    parameter Real expected;
    parameter Modelica.SIunits.Time at;
    parameter Real eps=1e-6;
    input Real actual;

  algorithm
    when initial() then
      assert(at >= time,
        "The specified time is before the start of the simulation.");
    end when;
    when time >= at then
      assertValue(
          actual,
          expected,
          eps);
    end when;
    when terminal() then
      assert(time >= at,
        "The specified time is after the end of the simulation.");
      // Note:  In Dymola 7.4, the simulation log may state "Integration
      // terminated successfully" and then the assertion statement below it.
    end when;

  end AssertValueAt;

  package Tests "A library to test the assertion primitives in this library"
    extends Modelica.Icons.Package;
    package Trajectory "Tests on the AssertTrajectory model"
      extends Modelica.Icons.Package;
      model CheckSuccess
        Real x=time^2;
        AssertTrajectory check_x(actual=x, expected=[0, 0; 1, 1; 2, 4; 3, 9]);
        annotation (TestCase(action="simulate",result="success"), experiment(
              StopTime=4));

      end CheckSuccess;

      model CheckFailure1 "Check for failure when first point is before start"
        extends CheckSuccess(check_x(expected=[-1, 1; 0, 0; 1, 1]));
        annotation (TestCase(action="simulate",result="failure"), experiment(
              StopTime=4));

      end CheckFailure1;

      model CheckFailure2
        "Check for failure when values don't match during simulation"
        extends CheckSuccess(x=time);
        annotation (TestCase(action="simulate",result="failure"), experiment(
              StopTime=4));

      end CheckFailure2;

      model CheckFailure3 "Check for failure when all points aren't checked"
        extends CheckSuccess(check_x(expected=[0, 0; 1, 1; 5, 25]));
        annotation (TestCase(action="simulate",result="failure"), experiment(
              StopTime=4));

      end CheckFailure3;

    end Trajectory;

    package BecomesTrueAt "Tests on the AssertBecomesTrueAt model"
      extends Modelica.Icons.Package;
      model CheckSuccess
        Real x=time;
        AssertBecomesTrueAt check_event(event=(x > 2), at=2);
        annotation (TestCase(action="simulate",result="success"), experiment(
              StopTime=4));

      end CheckSuccess;

      model CheckFailure1
        "Check for failure when expected transition is before simulation start"
        extends CheckSuccess(check_event(at=-1));
        annotation (TestCase(action="simulate",result="failure"), experiment(
              StopTime=4));

      end CheckFailure1;

      model CheckFailure2 "Check for failure when transition is early"
        extends CheckSuccess(check_event(event=(x > 1)));
        annotation (TestCase(action="simulate",result="failure"), experiment(
              StopTime=4));

      end CheckFailure2;

      model CheckFailure3 "Check for failure when transition is late"
        extends CheckSuccess(check_event(event=(x > 3)));
        annotation (TestCase(action="simulate",result="success"), experiment(
              StopTime=4));

      end CheckFailure3;

      model CheckFailure4
        "Check for failure when expected transition is after simulation end"
        extends CheckSuccess(check_event(at=5)) annotation (TestCase(action=
                "simulate", result="failure"), experiment(StopTime=4));

      end CheckFailure4;

    end BecomesTrueAt;

    package Initial "Tests associated with AssertInitial model"
      extends Modelica.Icons.Package;
      model CheckSuccess
        Real x=2*time + 1;
        AssertInitial check_x(actual=x, expected=1);
        annotation (TestCase(action="simulate",result="success"));

      end CheckSuccess;

      model CheckFailure1 "Check for failure when initial value is incorrect"
        extends CheckSuccess(x=time);
        annotation (TestCase(action="simulate",result="failure"));

      end CheckFailure1;

    end Initial;

    package Final "Tests associated with AssertFinal model"
      extends Modelica.Icons.Package;
      model CheckSuccess
        Real x=2*time + 1;
        AssertFinal check_x(actual=x, expected=9);
        annotation (TestCase(action="simulate",result="success"), experiment(
              StopTime=4));

      end CheckSuccess;

      model CheckFailure1 "Check for failure when final value is incorrect"
        extends CheckSuccess(x=time);
        annotation (TestCase(action="simulate",result="failure"), experiment(
              StopTime=4));

      end CheckFailure1;

    end Final;

    package Average "Tests associated with AssertAverage model"
      extends Modelica.Icons.Package;
      constant Real pi=3.141592653589793238462643383279502884197169399;
      model CheckSuccess "Check for successful value of average"
        Real x=sin(2*time);
        AssertAverageBetween check_x(
          average=0,
          start=0,
          finish=2*pi,
          signal=x,
          eps=1e-4);
        annotation (TestCase(action="simulate",result="success"), experiment(
              StopTime=8));

      end CheckSuccess;

      model CheckFailure1
        "Check for failure when starting in the middle of the interval"
        extends CheckSuccess(check_x(start=-pi,finish=pi));
        annotation (TestCase(action="simulate",result="failure"), experiment(
              StopTime=8));

      end CheckFailure1;

      model CheckFailure2
        "Check for failure when simulation ends before interval"
        extends CheckSuccess(check_x(start=pi,finish=3*pi));
        annotation (TestCase(action="simulate",result="failure"), experiment(
              StopTime=8));

      end CheckFailure2;

      model CheckFailure3 "Check for failure when values don't agree"
        extends CheckSuccess(check_x(finish=7*pi/8));
        annotation (TestCase(action="simulate",result="failure"), experiment(
              StopTime=8));

      end CheckFailure3;

    end Average;

    package ValueAt "Tests associated with AssertValueAt model"
      extends Modelica.Icons.Package;
      model CheckSuccess
        Real x=2*time + 1;
        AssertValueAt check_x(
          actual=x,
          expected=5,
          at=2);
        annotation (TestCase(action="simulate",result="success"), experiment(
              StopTime=4));

      end CheckSuccess;

      model CheckFailure1
        "Check for failure when value is specified before simulation start"
        extends CheckSuccess(check_x(at=-1));
        annotation (TestCase(action="simulate",result="failure"), experiment(
              StopTime=4));

      end CheckFailure1;

      model CheckFailure2 "Check for failure when value is incorrect"
        extends CheckSuccess(x=time);
        annotation (TestCase(action="simulate",result="failure"), experiment(
              StopTime=4));

      end CheckFailure2;

      model CheckFailure3
        "Check for failure when value is specified after simulation end"
        extends CheckSuccess(check_x(at=5));
        annotation (TestCase(action="simulate",result="failure"), experiment(
              StopTime=4));

      end CheckFailure3;

    end ValueAt;

    package Values "Tests associated with assertValues function"
      extends Modelica.Icons.Package;
      function CheckSuccess

      algorithm
        assertValues(
                actual={1,1 + 1e-4},
                expected={1,1},
                eps=1e-4);
        annotation (TestCase(action="call",result="success"));
      end CheckSuccess;

      function CheckFailure1 "Check for failure when a value is incorrect"

      algorithm
        assertValues(actual={1,2}, expected={1,1});
        annotation (TestCase(action="call",result="failure"));
      end CheckFailure1;

    end Values;

    package LogValue "Tests associated with assertLogValue function"
      extends Modelica.Icons.Package;
      function CheckSuccess

      algorithm
        assertLogValue(actual=10, expected=1);
        annotation (TestCase(action="call",result="success"));
      end CheckSuccess;

      function CheckFailure1 "Check for failure when value is incorrect"

      algorithm
        assertLogValue(actual=11, expected=1);
        annotation (TestCase(action="call",result="failure"));
      end CheckFailure1;

    end LogValue;

    package LogValues "Tests associated with assertLogValues function"
      extends Modelica.Icons.Package;
      function CheckSuccess

      algorithm
        assertLogValues(actual={0,1,100}, expected={0,10,10});
        annotation (TestCase(action="call",result="success"));
      end CheckSuccess;

      function CheckFailure1 "Check for failure when a value is incorrect"

      algorithm
        assertLogValues(actual={1,101}, expected={10,10});
        annotation (TestCase(action="call",result="failure"));
      end CheckFailure1;

    end LogValues;

  end Tests;
  annotation (Documentation(info="<html><p>This package is modified from
  XogenyTest version 1.1 by Michael Tiller of Xogeny, Inc.  XogenyTest is available at
  <a href=\"https://github.com/xogeny/XogenyTest\">https://github.com/xogeny/XogenyTest</a> under
  a <a href=\"https://creativecommons.org/licenses/by/3.0/deed.en_US\">Creative Commons
  Attribution 3.0 Unported License</a>.  If the <a href=\"modelica://FCSys.Tests\">Tests</a>
  package is removed from the <a href=\"modelica://FCSys\">FCSys</a> distribution, then this
  package can be safely removed as well.</p></html>"));

end Test;

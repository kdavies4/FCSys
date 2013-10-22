within ;
package FCSysTest "Library to test components of FCSys"
extends Modelica.Icons.Package;


function callAll
  "<html>Call all of the test functions for <a href=\"modelica://FCSys\">FCSys</a></html>"
  import Modelica.Utilities.Files;
  import Modelica.Utilities.Streams.print;
  extends Modelica.Icons.Function;

  input String logFile="FCSysTestLog.txt"
    "Filename where the log of all functions is stored";
  output Boolean ok "true, if all tests passed";

protected
  String file;

  algorithm
  file := Files.fullPathName(logFile);
  print("... callAll(...) is logged in " + file);

  if file <> "" then
    Files.removeFile(file);
  end if;

  ok := Units.callAll(logFile) and Utilities.callAll(logFile);
  annotation (Documentation(info="<html><p>This function call will fail if any of the functions return an
  incorrect result.  It will return <code>true</code> if all of the functions pass.
  There are no inputs.</p></html>"));
  end callAll;


annotation (Icon(graphics={Polygon(
        points={{-70,0},{-44,0},{-24,-34},{50,56},{78,56},{-24,-74},{-70,0}},
        lineColor={75,138,73},
        smooth=Smooth.None,
        fillColor={75,138,73},
        fillPattern=FillPattern.Solid)}), uses(Modelica(version="3.2.1")));
end FCSysTest;

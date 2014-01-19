within ;
package FCSys "Modelica fuel cell library"

// Maximum line width before a new word is wrapped in the code listing in
// the LaTeX document (76 characters, including leading spaces and // )
// -------------------------------------------------------------------------    


extends Modelica.Icons.Package;
// extends FCSys.Icons.Cell;

// Units and quantities
import U = FCSys.Units;
import Q = FCSys.Quantities;

// Enumerations
import FCSys.Characteristics.BaseClasses.Phase;
import FCSys.Characteristics.BaseClasses.ReferenceEnthalpy;
import FCSys.Species.Enumerations.Axis;
import FCSys.Species.Enumerations.Orient;
import FCSys.Species.Enumerations.Side;
import FCSys.Species.Enumerations.ConsThermo;
import FCSys.Species.Enumerations.ConsTrans;
import FCSys.Species.Enumerations.Init;


















annotation (
  preferredView="info",
  uses(Modelica(version="3.2.1")),
  Commands(executeCall=FCSys.Units.setup() "Re-initialize the units."),
  Documentation(info="<html>
    <p><a href=\"modelica://FCSys\">FCSys</a> is a free, open-source library of
    equation-based, object-oriented (EOO) models of proton exchange membrane
    fuel cells (PEMFCs) in the <a href = \"http://www.modelica.org/\">Modelica</a>
    language.  The models are:</p>

    <ul>
    <li><b>Dynamic</b></li>
    
    <li><b>Multi-domain</b>: 
    Chemical, electrical, fluid, and thermal phenomena are included.</li>
    
    <li><b>Multi-phase</b>: 
    Water is included and transported independently as vapor, liquid, and absorbed in the ionomer.  
    Phase change is represented as a dynamic process.</li>

    <li><b>Multi-dimensional</b></li>
    
    <li><b>Highly reconfigurable</b>: 
    There are options to adjust the assumptions, dimensionality (1D, 2D, or 3D), and spatial discretization 
    (i.e., resolution).  Species may be independently enabled at instantiation, unlike 
    the <a href=\"modelica://Modelica.Media\">Modelica media library</a>.  The framework is generic and can be extended 
    to other fluidic or electrochemical devices like batteries.</li>

    <li><b>Highly modular</b>: 
    Each layer of the cell is a separate model which is hierarchically constructed from graphical models of 
    subregions, phases, and species.  At each level, EOO (i.e., effort/flow) connectors are used to combine 
    the various components.</li>
                        
    <li><b>Fully declarative</b>:
    There are no causal connectors besides those used to apply boundary conditions.  Functions are only 
    used to simplify subexpressions of equations.</li>
    
    <li><b>Physics-based</b>: 
    The equations are based on first principles, with explicit conservation of material, momentum, and energy
    in every control volume and across every interface.  A unique and physically appropriate method of 
    upstream discretization is used to describe coupled advective and diffusive transfer.  All physical
    quantities are mapped to universal physical constants using a novel, flexible implementation of natural units.</li>
    
    <li><b>Computationally efficient</b>:  
    There are minimal switching events and no nonlinear systems of equations after appropriate translation.  
    A typical polarization curve can be simulated in less than two seconds.</li>

    </ul>
    
    <p><a href=\"#Fig1\">Figure 1</a> shows the seven primary layers of a typical PEMFC, which are also the components of the 
    fuel cell model shown in <a href=\"#Fig2\">Figure 2</a>.
    Fluid enters and exits the cell through channels in the flow plates (FPs).  It spreads through
    the gas diffusion diffusion layers (GDLs) and reacts in the catalyst layers (CLs) according to the following chemical equations:</p>
      
          <table border=0 cellspacing=0 cellpadding=2 align=center style=\"margin-left: auto;
margin-right: auto;\" class=noBorder>
      <tr>
        <td align=right style=\"white-space:nowrap; text-align:right;\" class=noBorder>
          2(H<sub>2</sub>
        </td>
        <td align=center style=\"white-space:nowrap; text-align:center;\" class=noBorder>
          &rarr;
        </td>
        <td align=left style=\"white-space:nowrap;\" class=noBorder>
          2e<sup>-</sup> + 2H<sup>+</sup>)
        </td>
        <td class=noBorder>
          (anode)
        </td>
      </tr>
      <tr>
        <td align=right style=\"white-space:nowrap; text-align=right;\" class=noBorder>
          4e<sup>-</sup> + 4H<sup>+</sup> + O<sub>2</sub>
        </td>
        <td align=center style=\"white-space:nowrap; text-align:center;\" class=noBorder>
          &rarr;
        </td>
        <td align=left style=\"white-space:nowrap;\" class=noBorder>
          2H<sub>2</sub>O
        </td>
        <td class=noBorder>
          (cathode)
        </td>
      </tr>
      <tr>
        <td colspan=4 class=noBorder>
          <hr>
        </td>
      </tr>
      <tr>
        <td align=right style=\"white-space:nowrap; text-align=right;\" class=noBorder>
          2H<sub>2</sub> + O<sub>2</sub>
        </td>
        <td align=center style=\"white-space:nowrap; text-align:center;\" class=noBorder>
          &rarr;
        </td>
        <td align=left style=\"white-space:nowrap;\" class=noBorder>
          2H<sub>2</sub>O
        </td>
        <td class=noBorder>
          (net)
        </td>
      </tr>
    </table>
      <p>The
    proton exchange membrane (PEM) prevents electronic transport; therefore, electrons must
    pass through an external load to sustain the net reaction.</p>

    <p align=center id=\"Fig1\"><img src=\"modelica://FCSys/Resources/Documentation/CellFlows.png\">
<br>Figure 1: Layers and primary flows of a PEMFC.</p>

    <!--<p align=center id=\"Fig2\"><img src=\"modelica://FCSys/help/FCSys.Assemblies.Cells.CellD.png\" width=600>-->
    <p align=center id=\"Fig2\"><a href=\"modelica://FCSys.Assemblies.Cells.Cell\"><img src=\"modelica://FCSys/Resources/Documentation/FCSys.Assemblies.Cells.CellD.png\"></a>
<br>Figure 2: Diagram of the <a href=\"modelica://FCSys.Assemblies.Cells.Cell\">PEMFC model</a>.</p>

    <p>The fuel cell model can be exercised using the test stand shown in <a href=\"#Fig3\">Figure 3</a> or connected to the <a href=\"modelica://Modelica.Fluid\">Modelica fluid library</a>
    using <a href=\"modelica://FCSys.Conditions.Adapters.MSL\">available adapters</a>.
    Please see the <a href=\"modelica://FCSys.UsersGuide.SampleResults.Cell\">sample cell results</a> for examples and the
    <a href=\"modelica://FCSys.UsersGuide.GettingStarted\">getting started page</a> to information about using the library.</p>

    <!--<p align=center id=\"Fig3\"><img src=\"modelica://FCSys/help/FCSys.Assemblies.Cells.Examples.TestStandD.png\" width=500>-->
    <p align=center id=\"Fig3\"><a href=\"modelica://FCSys.Assemblies.Cells.Examples.TestStand\"><img src=\"modelica://FCSys/Resources/Documentation/FCSys.Assemblies.Cells.Examples.TestStandD.png\"></a>
<br>Figure 3: Diagram of the <a href=\"modelica://FCSys.Assemblies.Cells.Examples.TestStand\">test stand model</a>.</p>

    <p><b>Licensed by the Georgia Tech Research Corporation under the Modelica License 2</b>
<br>Copyright 2007&ndash;2013, <a href=\"http://www.gtrc.gatech.edu/\">Georgia Tech Research Corporation</a>.</p>

    <p><i>This Modelica package is <u>free</u> software and the use is completely
    at <u>your own risk</u>; it can be redistributed and/or modified under the
    terms of the Modelica License 2. For license conditions (including the
    disclaimer of warranty) see
    <a href=\"modelica://FCSys.UsersGuide.License\">
    FCSys.UsersGuide.License</a>
    or visit <a href=\"http://www.modelica.org/licenses/ModelicaLicense2\">
    http://www.modelica.org/licenses/ModelicaLicense2</a>.</i></p></html>",
      revisions="<html>
    <ul>
    <li><a href=\"mailto:kdavies4@gmail.com\">Kevin Davies</a>, x/x/2013:<br>Version x.x.x (initial release)</li>
    </ul>
    </html>"),
  Icon(graphics={
      Polygon(
        points={{-4,52},{-14,42},{6,42},{16,52},{-4,52}},
        lineColor={0,0,0},
        fillColor={255,255,255},
        fillPattern=FillPattern.Solid),
      Polygon(
        points={{-30,52},{-40,42},{-20,42},{-10,52},{-30,52}},
        lineColor={0,0,0},
        fillPattern=FillPattern.HorizontalCylinder,
        fillColor={0,192,0}),
      Polygon(
        points={{-10,52},{-20,42},{-14,42},{-4,52},{-10,52}},
        lineColor={0,0,0},
        fillPattern=FillPattern.HorizontalCylinder,
        fillColor={0,0,0}),
      Rectangle(
        extent={{6,42},{12,-52}},
        fillPattern=FillPattern.Solid,
        fillColor={0,0,0},
        pattern=LinePattern.None),
      Polygon(
        points={{16,52},{6,42},{12,42},{22,52},{16,52}},
        lineColor={0,0,0},
        fillPattern=FillPattern.HorizontalCylinder,
        fillColor={0,0,0}),
      Line(
        points={{-40,42},{-40,-52}},
        pattern=LinePattern.None,
        smooth=Smooth.None),
      Polygon(
        points={{-46,64},{-66,44},{-46,44},{-26,64},{-46,64}},
        lineColor={0,0,0},
        fillColor={135,135,135},
        fillPattern=FillPattern.Solid),
      Rectangle(
        extent={{-39.6277,31.7996},{-67.912,17.6573}},
        lineColor={0,0,0},
        fillPattern=FillPattern.HorizontalCylinder,
        rotation=45,
        fillColor={255,255,255},
        origin={56.5067,67.5353}),
      Rectangle(
        extent={{-14,42},{6,-52}},
        lineColor={0,0,0},
        fillPattern=FillPattern.VerticalCylinder,
        fillColor={255,255,255}),
      Line(points={{-30,52},{32,52}}, color={0,0,0}),
      Rectangle(
        extent={{-5.21738,-5.21961},{-33.5017,-33.5041}},
        lineColor={0,0,170},
        fillPattern=FillPattern.VerticalCylinder,
        rotation=45,
        fillColor={0,0,240},
        origin={31.9983,69.3803}),
      Rectangle(
        extent={{12,42},{52,-52}},
        lineColor={0,0,170},
        fillPattern=FillPattern.VerticalCylinder,
        fillColor={0,0,240}),
      Polygon(
        points={{-26,64},{-46,44},{-46,-64},{-26,-44},{-26,64}},
        lineColor={0,0,0},
        fillColor={95,95,95},
        fillPattern=FillPattern.Solid),
      Rectangle(
        extent={{-5.21774,-5.2196},{-33.502,-33.5042}},
        lineColor={196,11,40},
        fillPattern=FillPattern.HorizontalCylinder,
        rotation=45,
        fillColor={253,52,56},
        origin={-30.001,79.3803}),
      Rectangle(
        extent={{-60,42},{-20,-52}},
        lineColor={196,11,40},
        fillPattern=FillPattern.VerticalCylinder,
        fillColor={253,52,56}),
      Rectangle(
        extent={{-60,42},{-40,-54}},
        fillPattern=FillPattern.Solid,
        fillColor={95,95,95},
        pattern=LinePattern.None,
        lineColor={0,0,0}),
      Rectangle(
        extent={{-76.648,66.211},{-119.073,52.0689}},
        lineColor={95,95,95},
        fillPattern=FillPattern.HorizontalCylinder,
        rotation=45,
        fillColor={135,135,135},
        origin={65.0166,81.3801}),
      Rectangle(
        extent={{-66,44},{-46,-64}},
        lineColor={95,95,95},
        fillPattern=FillPattern.VerticalCylinder,
        fillColor={135,135,135}),
      Polygon(
        points={{46,64},{34,52},{-26,52},{-26,64},{46,64}},
        smooth=Smooth.None,
        fillPattern=FillPattern.Solid,
        fillColor={230,230,230},
        pattern=LinePattern.None,
        lineColor={0,0,0}),
      Rectangle(
        extent={{-76.648,66.211},{-119.073,52.0689}},
        lineColor={95,95,95},
        fillPattern=FillPattern.HorizontalCylinder,
        rotation=45,
        fillColor={135,135,135},
        origin={157.017,81.3801}),
      Rectangle(
        extent={{26,44},{46,-64}},
        lineColor={95,95,95},
        fillPattern=FillPattern.VerticalCylinder,
        fillColor={135,135,135}),
      Polygon(
        points={{-26,64},{-26,52},{-30,52},{-30,60},{-26,64}},
        smooth=Smooth.None,
        fillColor={95,95,95},
        fillPattern=FillPattern.Solid,
        pattern=LinePattern.None,
        lineColor={0,0,0}),
      Ellipse(
        extent={{-44,62},{-36,58}},
        lineColor={135,135,135},
        fillColor={253,52,56},
        fillPattern=FillPattern.Sphere),
      Ellipse(
        extent={{36,50},{44,46}},
        lineColor={135,135,135},
        fillColor={0,0,240},
        fillPattern=FillPattern.Sphere),
      Polygon(
        points={{-26,64},{-26,52},{-30,52},{-40,42},{-46,44},{-26,64}},
        smooth=Smooth.None,
        fillColor={95,95,95},
        fillPattern=FillPattern.Solid,
        pattern=LinePattern.None,
        lineColor={0,0,0}),
      Line(
        points={{-30,52},{-40,42}},
        color={0,0,0},
        smooth=Smooth.None),
      Polygon(
        points={{66,64},{46,44},{46,-64},{66,-44},{66,64}},
        lineColor={0,0,0},
        fillColor={95,95,95},
        fillPattern=FillPattern.Solid),
      Rectangle(extent={{26,44},{46,-64}}, lineColor={0,0,0}),
      Rectangle(extent={{-66,44},{-46,-64}}, lineColor={0,0,0}),
      Line(
        points={{-26,64},{-26,52}},
        color={0,0,0},
        smooth=Smooth.None),
      Line(points={{-30,52},{34,52}}, color={0,0,0}),
      Rectangle(
        extent={{-46,74},{66,64}},
        pattern=LinePattern.None,
        fillColor={230,230,230},
        fillPattern=FillPattern.Solid,
        lineColor={0,0,0}),
      Polygon(
        points={{-46,64},{-26,64},{-46,44},{-66,44},{-46,64}},
        lineColor={0,0,0},
        smooth=Smooth.None),
      Polygon(
        points={{46,64},{66,64},{46,44},{26,44},{46,64}},
        lineColor={0,0,0},
        smooth=Smooth.None),
      Rectangle(extent={{-40,42},{26,-52}}, lineColor={0,0,0}),
      Rectangle(
        extent={{-20,42},{-14,-52}},
        fillPattern=FillPattern.Solid,
        fillColor={0,0,0},
        pattern=LinePattern.None)}),
  version="0.2.3",
  dateModified="2014-01-19 21:25:35Z",
  revisionID="SHA: b36f15d");
end FCSys;

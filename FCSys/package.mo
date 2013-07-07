within ;
package FCSys "Modelica fuel cell library"

// Maximum line width before a new word is wrapped in the code listing in
// the LaTeX document (76 characters, including leading spaces and // )
// ------------------------------------------------------------------------


extends Modelica.Icons.Package;
// extends FCSys.BaseClasses.Icons.Cell;

// Units and quantities
import U = FCSys.Units;
import Q = FCSys.Quantities;

// Enumerations
import FCSys.BaseClasses.Axis;
import FCSys.BaseClasses.Orientation;
import FCSys.BaseClasses.Side;
import FCSys.Characteristics.BaseClasses.Phase;
import FCSys.Characteristics.BaseClasses.ReferenceEnthalpy;
import FCSys.Subregions.BaseClasses.InitCapacitor;
import FCSys.Subregions.Species.BaseClasses.Conservation;
import FCSys.Subregions.Species.BaseClasses.InitScalar;
import FCSys.Subregions.Species.BaseClasses.InitTranslational;


package UsersGuide "User's Guide"
  extends Modelica.Icons.Information;

  model GettingStarted "Getting started"
    extends Modelica.Icons.Information;
    annotation (
      preferedView="info",
      DocumentationClass=true,
      Documentation(info="<html>
    <table border=1 cellspacing=0 cellpadding=4 bordercolor=#86989b >
      <tr bgcolor=#afc1c4>
        <td style=\"border: 1px solid #000000;\">
        <b><font color=#ffffff>Note</font></b>
        </td>
      </tr>
      <tr bgcolor=#f7f7f7>
        <td>
        FCSys should be compatible with any
        modeling environment that supports Modelica 3.1.  The following tools have been tested:
        <ul>
          <li>Dymola: Supported by versions 7.4, 2012 FD01, and 2013 FD01.
          Dymola's annotations for parameter dialogs and replaceable choices are included.</li>
          <li>MapleSim: Not supported as of version 4.5</li>
          <li>MWorks: Not supported as of version 2.6.10</li>
          <li>OpenModelica: Not supported as of version 1.8.1</li>
          <li>SystemModeler: Not supported as of version 3.0</li>
          <!---<li>JModelica: 9/17/12: Not able to get version 1.8 running on Windows or Linux </li>--->
        </ul>
        </td>
      </tr>
    </table>

    <p>These are the suggested steps to begin using FCSys:
    <ol>
        <li>Read the overview in the <a href=\"modelica://FCSys\">top-level documentation of FCSys</a>.</li>
        <li>Browse the subpackages of FCSys.  In general, the subpackages are
        ordered by the level of the model and the physical hierarchy (high-level at the top).
        <li>Call <a href=\"modelica://FCSys.Units.setup\">FCSys.Units.setup</a>() to
        establish the display units.  This is automatic if FCSys
        is loaded via the <a href=\"modelica://FCSys/../load.mos\">load.mos</a> script.
        <li>Simulate the <a href=\"modelica://FCSys.Assemblies.Cells.Examples.CellProfile\">FCSys.Assemblies.Cells.Examples.CellProfile</a>
        model.
        There are scripts in <a href=\"modelica://FCSys/Resources/Scripts/Dymola/README.md\">Resources/Scripts/Dymola/</a> to
        create useful plots of that model and others.
        The scripts should be accessible from the \"Command\" menu of the Modelica environment.
        For more detailed
        analysis, including spatial property distributions and vector plots,
        a Python module called FCRes is available in
        <a href=\"modelica://FCSys/Resources/Source/Python/README.md\">Resources/Source/Python/</a>
        (HTML and PDF documentation <a href=\"modelica://FCSys/Resources/Source/Python/doc/index.html\">here</a> and
        <a href=\"modelica://FCSys/Resources/Source/Python/doc/FCRes.pdf\">here</a>).</li>
        <li>Read the documentation of the classes. In particular, these may be of interest:
        <ul>
            <li><a href=\"modelica://FCSys.Units\">FCSys.Units</a> package:
            Information about the system of units, which is different
            than <a href=\"modelica://Modelica.SIunits\">Modelica.SIunits</a></li>
            <li><a href=\"modelica://FCSys.Connectors\">FCSys.Connectors</a> package:
            Overview of the connectors</li>
            <li><a href=\"modelica://FCSys.Subregions.Volume\">FCSys.Subregions.Volume</a> model:
            Discussion of how species are mixed within phases and how phases are combined within subregions</li>
            <li><a href=\"modelica://FCSys.Subregions.Species.Species\">FCSys.Subregions.Species.Species</a> model:
            Details about the exchange, transport, and storage of material, translational momentum, and
            energy</li>
            <li>A region or layer model such as <a href=\"modelica://FCSys.Regions.AnFPs.AnFP\">FCSys.Regions.AnFPs.AnFP</a>:
            Information about the geometric orientation of the cell</li>
            <li><a href=\"modelica://FCSys.Subregions.DepletionLayer\">FCSys.Subregions.DepletionLayer</a> model:
            Discussion of the electrochemical reaction</li>
        </ul>
        In general, overviews are given in the documentation of containing packages and
        detailed information is given at the appropriate level of inheritance.  If a model does not
        have sufficient documentation, please look at its base model(s) and the package(s) that
        contain it.  Assumptions are only listed at the lowest level of inheritance at which they apply.  Therefore, the
        list of assumptions in a model should be considered in conjunction with the assumptions in all
        the models it inherits from.
        </li>
        <li>Create and simulate examples of other usage scenarios.  Many of the
        nested models (regions, subregions, species) are replaceable.
        Their parameters are often not propagated to the cell level, but may
        be accessed through the parameter dialog by editing the replaceable model.  Note that
        many models have auxiliary output variables for analysis and diagnostics.
        These may be included by setting <code>analysis=true</code> in the outer environment model (instance
        of <a href=\"modelica://FCSys.Conditions.Environment\">Environment</a>).</li>
        <li>Develop your own classes.  It should be possible to model other electrochemical
        devices (solid oxide fuel cells, lithium ion batteries, flow batteries/regenerative fuel cells, etc.) by
        extending the existing classes and
        following the existing framework.  It will be necessary to add species models
        (Li<sup>+</sup>, O<sup>2-</sup>, etc.).</li>
        <li>Please share your additions or modifications to the source code so that the library
        can be improved and others may benefit.  The best way is to create a fork from the
        development page at <a href=\"https://github.com/kdavies4/FCSys\">https://github.com/kdavies4/FCSys</a>.
        Please also feel free to use the <a href=\"modelica://FCSys.UsersGuide.Contact\">contact information</a>.</li>
    </ol></p>
    </html>"));

    end GettingStarted;

  model SampleResults "Sample results"
    extends Modelica.Icons.Information;
    annotation (
      preferedView="info",
      DocumentationClass=true,
      Documentation(info="<html><p>The figures below show results from some of the examples in
    FCSys.  The figures are ordered from low-level tests to evaluations
    of the complete cell.  For more information about any of the results, please
    follow the links to the associated models.  The models have been simulated using Dymola 7.4 and the plots have been
    generated using <a href=\"http://kdavies4.github.io/ModelicaRes/\">ModelicaRes</a> and the commands in the
    <a href=\"modelica://FCSys/Resources/Documentation/UsersGuide/SampleResults/plot-samples.py\">plot-samples.py</a> script.</p>

    <p align=center id=\"Fig1\"><a href=\"modelica://FCSys.Subregions.Examples.SubregionsSound\"><img src=\"modelica://FCSys/Resources/Documentation/UsersGuide/SampleResults/SubregionsSound.png\"></a>
    <br>Figure 1: Pressure waves reflecting across two 1 cm<sup>3</sup> regions with an initial pressure difference.  The nonlinearity is due to upstream discretization (<a href=\"modelica://FCSys.Subregions.Examples.SubregionsSound\">FCSys.Subregions.Examples.SubregionsSound</a>).</p>

    <p align=center id=\"Fig2\"><a href=\"modelica://FCSys.Subregions.Examples.ThermalConduction\"><img src=\"modelica://FCSys/Resources/Documentation/UsersGuide/SampleResults/ThermalConduction.png\"></a>
    <br>Figure 2: Thermal conduction through a graphite bar divided into 1 cm<sup>3</sup> subregions (<a href=\"modelica://FCSys.Subregions.Examples.ThermalConduction\">FCSys.Subregions.Examples.ThermalConduction</a>).</p>

    <p align=center id=\"Fig3\"><a href=\"modelica://FCSys.Subregions.Examples.ThermalConductionConvection\"><img src=\"modelica://FCSys/Resources/Documentation/UsersGuide/SampleResults/ThermalConductionConvection.png\"></a>
    <br>Figure 3: Velocity induced in gas in contact with graphite undergoing thermal conduction (<a href=\"modelica://FCSys.Subregions.Examples.ThermalConductionConvection\">FCSys.Subregions.Examples.ThermalConductionConvection</a>).</p>

    <p align=center id=\"Fig4\"><a href=\"modelica://FCSys.Subregions.Examples.SaturationPressure\"><img src=\"modelica://FCSys/Resources/Documentation/UsersGuide/SampleResults/SaturationPressure.png\"></a>
    <br>Figure 4: Validation of the H<sub>2</sub>O saturation pressure curve derived from Gibbs equilibrium in <a href=\"modelica://FCSys\">FCSys</a>, as compared to
    <a href=\"modelica://Modelica.Media.Air.MoistAir\">Modelica.Media.Air.MoistAir</a> (<a href=\"modelica://FCSys.Subregions.Examples.SaturationPressure\">FCSys.Subregions.Examples.SaturationPressure</a>).</p>

    <p align=center id=\"Fig5\"><a href=\"modelica://FCSys.Subregions.Examples.SubregionEvaporation\"><img src=\"modelica://FCSys/Resources/Documentation/UsersGuide/SampleResults/SubregionEvaporation.png\"></a>
    <br>Figure 5: Dynamic H<sub>2</sub>O evaporation and condensation.  From 1 to 2 s, additional vapor is injected (<a href=\"modelica://FCSys.Subregions.Examples.SubregionEvaporation\">FCSys.Subregions.Examples.SubregionEvaporation</a>).</p>

    <p>**Add various polarization curves.</p>

    </html>"));

    end SampleResults;

  model Glossary "Glossary"
    extends Modelica.Icons.Information;
    annotation (
      preferedView="info",
      DocumentationClass=true,
      Documentation(info="<html><p>The following terms appear in this library with new or somewhat unique meanings:
    <ul>
      <li><b>configuration</b>&mdash;a species in a certain phase within a subregion<br>

      <li><b>dynamic compressibility</b>&mdash;reciprocal of volume, second, or bulk viscosity
      (see <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.beta\">&beta;</a>())<br>

      <li><b>exchange</b>&mdash;transfer of a conserved quantity among configurations within a region<br>

      <li><b>specific</b>&mdash;adjective that indicates the quotient of the following quantity and its associated particle number<br>

      <li><b>thermal independity</b>&mdash;extent to which an exchange of thermal energy between species causes or requires a temperature difference
      (see <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.nu\">&nu;</a>())<br>

      <li><b>translational Nusselt number</b>&mdash;correction to Newton's law of viscous shear for the shape of the flow profile<br>

      <li><b>transport</b>&mdash;transfer of a conserved quantity between adjacent subregions
    </ul></p>
    </table>
    </html>"));

    end Glossary;

  class References "References"
    extends Modelica.Icons.References;
    annotation (
      preferedView="info",
      DocumentationClass=true,
      Documentation(info="<html>
    <p>These papers describe work leading to the current FCSys package (most recent at the top):
    <ol>
    <li>K. L. Davies, \"<a href=\"modelica://FCSys/Resources/Documentation/UsersGuide/References/Declarative Modeling of PEMFCS (dissertation).pdf\">Declarative Modeling of Proton Exchange Membrane Fuel Cells for System Design</a>,\" Ph.D. dissertation, Georgia Institute of Technology, Aug. 2013.</li>
    <li>K. L. Davies, C. L. Haynes, and C. J. Paredis, \"<a href=\"http://www.ep.liu.se/ecp_article/index.en.aspx?issue=076;article=010\">Library for First-Principle Models of Proton Exchange Membrane Fuel Cells in Modelica</a>,\" in <i>Modelica Conference</i> (Munich, Germany), Modelica Assoc., Sep. 2012.</li>
    <li>K. L. Davies, \"<a href=\"http://www.ep.liu.se/ecp_article/index.en.aspx?issue=076;article=082\">Natural Unit Representation in Modelica</a>,\" in <i>Modelica Conference</i> (Munich, Germany), Modelica Assoc., Sep. 2012 (<a href=\"modelica://FCSys/Resources/Documentation/UsersGuide/References/Natural Unit Representation in Modelica (poster).pdf\">poster</a>).</li>
    <li>K. L. Davies, C. L. Haynes, and C. J. Paredis, \"<a href=\"http://www.modelica.org/events/modelica2009/Proceedings/memorystick/pages/papers/0106/0106.pdf\">Modeling Reaction and Diffusion Processes of Fuel Cells within Modelica</a>,\" in <i>Modelica Conference</i> (Como, Italy), Modelica Assoc., Sep. 2009.</li>
    <li>K. L. Davies, R. M. Moore, and G. Bender, \"<a href=\"http://www.modelica.org/events/modelica2009/Proceedings/memorystick/pages/papers/0107/0107.pdf\">Model Library of Polymer Electrolyte Membrane Fuel Cells for System Hardware and Control Design</a>,\" in <i>Modelica Conference</i> (Como, Italy), Modelica Assoc., Sep. 2009.</li>
    <li>K. L. Davies and R. M. Moore, \"<a href=\"http://link.aip.org/link/abstract/ECSTF8/v11/i1/p797/s1\">Object-Oriented Fuel Cell Model Library</a>,\" <i>Electrochem. Soc. T.</i>, vol. 11, no. 1, pp. 797&ndash;808, 2007.</li>
    </ol></p>

    <p>These references are cited throughout the library:</p>
    <table border=0 cellspacing=0 cellpadding=2 class=noBorder>
    <tr>
      <td valign=top class=noBorder>[<a name=\"Aronsson2009\" class=bookmark>Aronsson2009</a>]</td>
      <td valign=top class=noBorder>P. Aronsson and D. Broman, \"<a href=\"http://www.ep.liu.se/ecp_article/index.en.aspx?issue=043;article=105\">Extendable Physical Unit Checking with Understandable Error Reporting</a>,\"  in <i>Modelica Conference</i> (Como, Italy), Modelica Assoc., Sep. 2009.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"Avogadro1.03\" class=bookmark>Avogadro1.03</a>]</td>
      <td valign=top class=noBorder>Avogadro: An Open-Source Molecular Builder and Visualization Tool, ver. 1.03. <a href=\"http://avogadro.openmolecules.net\">http://avogadro.openmolecules.net</a>.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"Bejan2006\" class=bookmark>Bejan2006</a>]</td>
      <td valign=top class=noBorder>A. Bejan, <i>Advanced Engineering Thermodynamics</i>, John Wiley &amp; Sons, 3rd ed., 2006.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"Broman2008\" class=bookmark>Broman2008</a>]</td>
      <td valign=top class=noBorder>D. Broman and P. Aronsson and P. Fritzson, \"<a href=\"http://dx.doi.org/10.1149/1.2221251\">Design Considerations for Dimensional Inference and Unit Consistency Checking in Modelica</a>,\"  in <i>Modelica Conference</i> (Bielefeld, Germany), Modelica Assoc., Mar. 2008.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"Brown2011\" class=bookmark>Brown2011</a>]</td>
      <td valign=top class=noBorder>W. M. Brown, P. Wang, S. J. Plimpton, and A. N. Tharrington, \"<a href=\"http://dx.doi.org/10.1016/j.cpc.2010.12.021\">Implementing Molecular Dynamics on Hybrid High Performance Computers&mdash;Short Range Forces</a>,\" <i>Comput. Phys. Commun.</i>, vol. 182, no. 4, pp. 898&ndash;911, 2011.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"BIPM2006\" class=bookmark>BIPM2006</a>]</td>
      <td valign=top class=noBorder>International Bureau of Weights and Measures (BIPM), \"<a href=\"http://www.bipm.org/utils/common/pdf/si_brochure_8_en.pdf\">The International System of Units (SI)</a>,\" 8th ed., 2006.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"Bernardi1992\" class=bookmark>Bernardi1992</a>]</td>
      <td valign=top class=noBorder>D. M. Bernardi and M. W. Verbrugge, \"<a href=\"http://dx.doi.org/10.1149/1.2221251\">A Mathematical Model of the Solid-Polymer-Electrolyte Fuel Cell</a>,\" <i>J. Electrochem. Soc.</i>, vol. 139, no. 9, pp. 2477&ndash;2491, Sep. 1992.</td>
    </tr>
    <tr>
      <td valign=top class=\"noBorder\">[Dassault2010]</td>
      <td valign=top class=\"noBorder\">Dassault Syst&eacute;mes, <i>Dymola (Dynamic Modeling Laboratory) User Manual</i>, Lund, Sweden, Jan. 2010.  Distributed with Dymola version 7.4.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"DuPont2004N\" class=bookmark>DuPont2004N</a>]</td>
      <td valign=top class=noBorder>DuPont, \"Nafion&reg; PFSA Membranes N-112, NE-1135, N-115, N-117, NE-1110</a>,\" <a href = http://www.fuelcells.dupont.com>http://www.fuelcells.dupont.com</a>, Feb. 2004.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"DuPont2004NRE\" class=bookmark>DuPont2004NRE</a>]</td>
      <td valign=top class=noBorder>DuPont, \"Nafion&reg; PFSA Membranes NRE-211 and NRE-212</a>,\" <a href = http://www.fuelcells.dupont.com>http://www.fuelcells.dupont.com</a>, Feb. 2004.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"DuPont2005\" class=bookmark>DuPont2005</a>]</td>
      <td valign=top class=noBorder>DuPont, \"Nafion&reg; PFSA Membranes NE-1135, N-115, N-117, NE-1110</a>,\" <a href = http://www.fuelcells.dupont.com>http://www.fuelcells.dupont.com</a>, Feb. 2005.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"Dymond2002\" class=bookmark>Dymond2002</a>]</td>
      <td valign=top class=noBorder>J. H. Dymond, K. N. Marsh, R. C. Wilhoit, and K. C. Wong, <i>Virial Coefficients of Pure Gases</i>, Springer-Verlag, 2002.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"Entegris2012\" class=bookmark>Entegris2012</a>]</td>
      <td valign=top class=noBorder>Entegris, \"Industrial Graphite,\" <a href=\"http://www.entegris.com/Resources/assets/6204-7085-0312.pdf\">http://www.entegris.com/Resources/assets/6204-7085-0312.pdf</a>, Apr. 2012.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"Fritzson2004\" class=bookmark>Fritzson2004</a>]</td>
      <td valign=top class=noBorder>P. Fritzson, <i>Principles of Object-Oriented Modeling and Simulation with Modelica 2.1</i>, IEEE Press (Piscataway, NJ), 2004.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"Greiner1995\" class=bookmark>Greiner1995</a>]</td>
      <td valign=top class=noBorder>W. Greiner, L. Neisem and H. St&ouml;cker, \"<a href=\"http://books.google.com/books?id=12DKsFtFTgYC\">A Mathematical Model of the Solid-Polymer-Electrolyte Fuel Cell</a>,\" <i>J. Electrochem. Soc.</i>, vol. 139, no. 9, pp. 2477&ndash;2491, Sep. 1992.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"Gurau1998\" class=bookmark>Gurau1998</a>]</td>
      <td valign=top class=noBorder>V. Gurau, H. Liu, and S. Kaka, \"<a href=\"http://dx.doi.org/10.1002/aic.690441109\">Two-Dimensional Model for Proton Exchange Membrane Fuel Cells</a>,\" <i>AIChE J.</i>, vol. 44, no. 11, pp. 2410&ndash;2422, Nov. 1998.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"Hess2008\" class=bookmark>Hess2008</a>]</td>
      <td valign=top class=noBorder>B. Hess, C. Kutzner, D. van der Spoel, and E. Lindahl, \"GROMACS 4: Algorithms for Highly Efficient, Load-Balanced, and Scalable Molecular Simulation,\" <i>J. Chem. Theory Comput.</i>, vol. 4, no. 3, pp. 435&ndash;447, 2008.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"Incropera2002\" class=bookmark>Incropera2002</a>]</td>
      <td valign=top class=noBorder>F. P. Incropera and D. P. DeWitt, <i>Fundamentals of Heat and Mass Transport</i>, 5th ed., John Wiley &amp; Sons, 2002.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"Kandlikar2009\" class=bookmark>Kandlikar2009</a>]</td>
      <td valign=top class=noBorder>S. G. Kandlikar and Z. Lu, \"<a href=\"http://dx.doi.org/10.1016/j.applthermaleng.2008.05.009\">Thermal Management Issues in a PEMFC Stack&mdash;A Brief Review of Current Status</a>,\" <i>Appl. Therm. Eng.</i>, vol. 29, no. 7, pp. 1276&ndash;1280, 2009.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"Larminie2003\" class=bookmark>Larminie2003</a>]</td>
      <td valign=top class=noBorder>J. Larminie and A. Dicks, <i>Fuel Cell Systems Explained</i>, John Wiley &amp; Sons, 2003.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"Lin2006\" class=bookmark>Lin2006</a>]</td>
      <td valign=top class=noBorder>J. Lin, J. K. Lee, M. Kellner, R. Wycisk, and P. N. Pintauroa, \"<a href=\"http://dx.doi.org/10.1149/1.2196687\">Nafion-Flourinated Ethylene-Propylene Resin Membrane Blends for Direct Methanol Fuel Cells</a>,\" <i>J. Electrochem. Soc.</i>, vol. 153, no. 7, pp. A1325&ndash;A1331, 2006.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"Mark1999\" class=bookmark>Mark1999</a>]</td>
      <td valign=top class=noBorder>J. E. Mark, <i>Polymer Data Handbook</i>, Oxford University Press, 1999.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"Mattsson2008\" class=bookmark>Mattsson2008</a>]</td>
      <td valign=top class=noBorder>S. E. Mattsson and H. Elmqvist, \"Unit Checking and Quantity Conservation,\" in <i>Modelica Conf.</i> (Bielefeld, Germany), Modelica Assoc., Mar. 2008.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"Mattsson1993B\" class=bookmark>Mattsson1993B</a>]</td>
      <td valign=top class=noBorder>S. E. Mattsson and G. Soderlind, \"Index Reduction in Differential-Algebraic Equations Using Dummy Derivatives,\" <i>SIAM J. Sci. Comput.</i>, vol. 14, no. 3, pp. 677&ndash;692, May 1993.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"McBride1996\" class=bookmark>McBride1996</a>]</td>
      <td valign=top class=noBorder>B. J. McBride and S. Gordon, \"<a href=\"http://www.grc.nasa.gov/WWW/CEAWeb/RP-1311P2.htm\">Computer Program for Calculating Complex Chemical Equilibrium Compositions and Applications II. Users Manual and Program Description</a>,\" NASA Reference Publication 1311, Jun. 1996.  Recent data available at <a href=\"http://www.grc.nasa.gov/WWW/CEAWeb/ceaHome.htm\">http://www.grc.nasa.gov/WWW/CEAWeb/ceaHome.htm</a>.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"McBride2002\" class=bookmark>McBride2002</a>]</td>
      <td valign=top class=noBorder>B. J. McBride, M. J. Zehe, and S. Gordon, \"<a href=\"http://gltrs.grc.nasa.gov/cgi-bin/GLTRS/browse.pl?2002/TP-2002-211556.html\">NASA Glenn Coefficients for Calculating Thermodynamic Properties of Individual Species</a>,\" NASA report TP-2002-211556, Sep. 2002.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"Modelica3.2\" class=bookmark>Modelica3.2</a>]</td>
      <td valign=top class=\"noBorder\">Modelica Association, <i><a href=\"https://www.modelica.org/documents/ModelicaSpec32.pdf\">Modelica: A Unified Object-Oriented Language for Physical Systems Modeling: Language Specification</i></a>, Ver. 3.2, Mar. 2010.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"Moran2004\" class=bookmark>Moran2004</a>]</td>
      <td valign=top class=noBorder>M. J. Moran and H. N. Shapiro, <i>Fundamentals of Engineering Thermodynamics</i>, 5th ed., John Wiley &amp; Sons, 2004.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"NIST2010\" class=bookmark>NIST2010</a>]</td>
      <td valign=top class=noBorder>National Institute of Standards and Technology (NIST), \"Fundamental Physical Constants: Complete Listing,\" <a href=\"http://physics.nist.gov/cuu/Constants/Table/allascii.txt\">http://physics.nist.gov/cuu/Constants/Table/allascii.txt</a>, 2010.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"Nitta2008\" class=bookmark>Nitta2008</a>]</td>
      <td valign=top class=noBorder>I. Nitta, O. Himanen, and M. Mikkola, \"<a href=\"http://dx.doi.org/10.1002/fuce.200700054\">Thermal Conductivity and Contact Resistance of Compressed Gas Diffusion Layer of PEM Fuel Cell</a>,\" <i>Fuel Cells</i>, vol. 8, pp. 111&ndash;119, 2008.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"Present1958\" class=bookmark>Present1958</a>]</td>
      <td valign=top class=noBorder>R. D. Present, <i>Kinetic Theory of Gases</i>, McGraw Hill, 1958.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"Rao1997\" class=bookmark>Rao1997</a>]</td>
      <td valign=top class=noBorder>Y. V. C. Rao, <i>Chemical Engineering Thermodynamics</i>, Hyderabad, India:  Universities Press, 1997.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"Rapaport2004\" class=bookmark>Rapaport2004</a>]</td>
      <td valign=top class=noBorder>D. C. Rapaport, <i>The Art of Molecular Dynamics Simulation</i>, Cambridge University Press, 2nd ed., 2004.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"Reichert2010\" class=bookmark>Reichert2010</a>]</td>
      <td valign=top class=noBorder>P. Reichert and N. Schuwirth, \"<a href=\"http://dx.doi.org/10.1016/j.envsoft.2010.03.002\">A Generic Framework for Deriving Process Stoichiometry in Environmental Models</a>,\" <i>Environmental Modelling &amp; Software</i>, vol. 25, pp. 1241&ndash;1251, 2010.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"Salzman2004\" class=bookmark>Salzman2004</a>]</td>
      <td valign=top class=noBorder>W. R. Salzman, \"The Virial Expansion,\" <a href=\"http://www.chem.arizona.edu/~salzmanr/480a/480ants/VIRIAL/virial.html\">http://www.chem.arizona.edu/~salzmanr/480a/480ants/VIRIAL/virial.html</a>, Course notes for Physical Chemistry (Chemistry 480A), University of Arizona, Jul. 2004.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"Schetz1996\" class=bookmark>Schetz1996]</td>
      <td valign=top class=noBorder>J. A. Schetz and A. E. Fuhs (ed.), <i>Handbook of Fluid Dynamics and Fluid Machinery</i>, John Wiley &amp; Sons, vol. 1: Fundamentals of Fluid Dynamics, 1996.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"SGL2004\" class=bookmark>SGL2004</a>]</td>
      <td valign=top class=noBorder>SGL Carbon Group, \"Sigracet&reg; 24 &amp; 25 Series Gas Diffusion Layer,\" <a href=\"http://www.sglcarbon.com\">http://www.sglcarbon.com</a>, Sep. 2004.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"SGL2007\" class=bookmark>SGL2007</a>]</td>
      <td valign=top class=noBorder>SGL Carbon Group, \"Sigracet&reg; 10 Series Gas Diffusion Layer,\" <a href=\"http://www.sglcarbon.com\">http://www.sglcarbon.com</a>, Apr. 2007.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"Shah2009\" class=bookmark>Shah2009</a>]</td>
      <td valign=top class=noBorder>A. A. Shah, T. R. Ralph, and F. C. Walsh, \"Modeling and Simulation of the Degradation of Perfluorinated Ion-Exchange Membranes in PEM Fuel Cells,\" <i>J. Electrochem. Soc.</i>, vol. 156, no. 4, pp. B465&ndash;B484, 2009.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"Springer1991\" class=bookmark>Springer1991</a>]</td>
      <td valign=top class=noBorder>T. E. Springer, T. A. Zawodzinski, and S. Gottesfeld, \"<a href=\"http://dx.doi.org/10.1149/1.2085971\">Polymer Electrolyte Fuel Cell Model</a>,\" <i>J. Electrochem. Soc.</i>, vol. 138, pp. 2334&ndash;2342, Aug. 1991.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"Spry2009\" class=bookmark>Spry2009</a>]</td>
      <td valign=top class=noBorder>D. B. Spry and M. D. Fayer, \"<a href=\"http://www.stanford.edu/group/fayer/articles/384-392/385.pdf\">Proton Transfer and Proton Concentrations in Protonated Nafion Fuel Cell Membranes</a>,\" <i>J. Phys. Chem. B</i>, vol. 113, no. 30, pp. 10210&ndash;10221, 2009.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"\" class=bookmark>Svehla1995</a>]</td>
      <td valign=top class=noBorder>R. A. Svehla, \"<a href=\"http://www.grc.nasa.gov/WWW/CEAWeb/TM-4647.htm\">Transport Coefficients for the NASA Lewis Chemical Equilibrium Program</a>,\" NASA Technical Memorandum 4647, Apr. 1995.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"Takenaka1990\" class=bookmark>Takenaka1990</a>]</td>
      <td valign=top class=noBorder>M. Takenaka and R. Masui, \"<a href=\"http://iopscience.iop.org/0026-1394/27/4/001\">Measurement of the Thermal Expansion of Pure Water in the Temperature Range 0 &deg;C&ndash;85&deg;C</a>,\" <i>Metrologia</i>, vol. 27, pp. 165&ndash;171, 1990.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"Tissandier1998\" class=bookmark>Tissandier1998</a>]</td>
      <td valign=top class=noBorder>M. D. Tissandier, K. A. Cowen, W. Y. Feng, E. Gundlach, M. H. Cohen, A. D. Earhart, J. V. Coe, and T. R. Tuttle, Jr., \"<a href=\"http://www.uh.edu/~chembi/single_ion_hydration.PDF\">The Proton's Absolute Aqueous Enthalpy and Gibbs Free Energy of Solvation from Cluster-Ion Solvation Data</a>,\" <i>J. Phys. Chem. A</i>, vol. 102, pp. 7787&ndash;7794, 1998.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"Toray2010\" class=bookmark>Toray2010</a>]</td>
      <td valign=top class=noBorder>Toray Industries, Inc., \"Carbon Paper,\" <a href = http://www.torayca.com/index2.html>http://www.torayca.com/index2.html</a>, accessed 2010.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"Weber2004\" class=bookmark>Weber2004</a>]</td>
      <td valign=top class=noBorder>A. Z. Weber and J. Newman, \"<a href=\"http://dx.doi.org/10.1021/cr020729l\">Modeling Transport in Polymer-Electrolyte Fuel Cells</a>,\" <i>Chem. Rev.</i>, vol. 104, pp. 4679&ndash;4726, 2004.</td>
    </tr>
    <tr>
      <td valign=top class=noBorder>[<a name=\"Woo1995\" class=bookmark>Woo1995</a>]</td>
      <td valign=top class=noBorder>K. W. Woo and S. I. Yeo, \"Dalton's Law vs. Amagat's Law for the Mixture of Real Gases,\" <i>SNU J. Educ. Res.</i>, vol. 5, pp. 127&ndash;134, 1995.</td>
    </tr>
    </table>
    </html>"));

    end References;

  class Contact "Contact"
    extends Modelica.Icons.Contact;
    annotation (
      preferedView="info",
      DocumentationClass=true,
      Documentation(info="<html>
    <p>Updates to this package may be available online at the
    <a href=\"http://kdavies4.github.io/FCSys/\">main project site</a> or the
    <a href=\"https://modelica.org/libraries\">Modelica libraries page</a>.
    The development page is
    <a href=\"https://github.com/kdavies4/FCSys\">https://github.com/kdavies4/FCSys</a>.
    Please report any problems using the <a href=\"https://github.com/kdavies4/FCSys/issues\">issues</a>
    link on that page.</p>

    <p><b>Author:</b></p>
    <dd>Kevin Davies</dd>
    <dd><a href=\"http://www.me.gatech.edu\">George W. Woodruff School of Mechanical Engineering</a></dd>
    <dd><a href=\"http://www.gatech.edu\">Georgia Institute of Technology</a></dd>
    <dd>813 Ferst Drive</dd>
    <dd>Atlanta, GA 30332-0405</dd>
    <dd>USA</dd>
    <dd><a href=\"mailto:kdavies4@gmail.com\">kdavies4@gmail.com</a></dd>

    <p><b>Acknowledgments:</b><ul>
    <li>Guidance from Robert Moore, Comas Haynes, and Chris Paredis
    <li>Technical discussions and insight from Mike Angelo, Guido Bender, Severine Busquet,
    Chris Ford, Sebastian Herzig, Ben Lee, George Nelson, Mike Tiller, Hubertus Tummescheit, and Mebs Virji</li>
    <li>Source-code contributions and bug fixes from Mohammad Ali, Kevin Bandy, Martin Sj&ouml;lund, Francois Steinmetz, and Joerg Weiss-Ungeth&uuml;m</li>

    <li>Financial support from:
    <ul>
    <li>Presidential Fellowship from the <a href=\"http://www.me.gatech.edu\">George W. Woodruff
    School of Mechanical Engineering</a> and the <a href=\"http://www.gatech.edu\">Georgia Institute
    of Technology</a></li>
    <li>Robert G. Shackelford Fellowship from the <a href=\"http://www.gtri.gatech.edu\">Georgia
    Tech Research Institute</a></li>
    <li>Grant #N00014-04-0682 from the <a href=\"http://www.onr.navy.mil\">Office of Naval Research</a>
    to the <a href=\"http://www.hnei.hawaii.edu\">Hawaii Natural Energy Institute</a></li>
    </ul>
    </ul></p>

</html>"));

    end Contact;

  class ModelicaLicense2 "Modelica License 2"
    extends Modelica.Icons.Information;
    annotation (
      preferedView="info",
      DocumentationClass=true,
      Documentation(info="<html>
<style type=\"text/css\">
*       { font-size: 10pt; font-family: Arial,sans-serif; }
code    { font-size:  9pt; font-family: Courier,monospace;}
h6      { font-size: 10pt; font-weight: bold; color: green; }
h5      { font-size: 11pt; font-weight: bold; color: green; }
h4      { font-size: 13pt; font-weight: bold; color: green; }
address {                  font-weight: normal}
td      { solid #000; vertical-align:top; }
th      { solid #000; vertical-align:top; font-weight: bold; }
table   { solid #000; border-collapse: collapse;}
</style>

<p>All files in this directory (\"FCSys\") and all subdirectories are licensed by
<b><u>Georgia Tech Research Corporation</u></b> under the <b><u>Modelica
License 2</u></b> (below) with the additional condition:<ul>
  <li>This software is controlled under the jurisdiction of the United States
      Department of Commerce and subject to Export Administration Regulations.
      By downloading or using the Software, you are agreeing to comply with
      U. S. export controls.  Diversion contrary to law is prohibited.  The
      software cannot be exported or reexported to sanctioned countries that
      are controlled for Anti-Terrorism (15 CFR Part 738 Supplement 1) or to
      denied parties, <a
      href=\"http://www.bis.doc.gov/index.php/policy-guidance/lists-of-parties-of-concern\">
      http://www.bis.doc.gov/index.php/policy-guidance/lists-of-parties-of-concern</a>.
      EAR99 items cannot be exported or reexported to Iraq for a military
      purpose or to a military end-user (15 CFR Part 746.3).  Export and
      reexport include any release of technology to a foreign national within
      the United States.  Technology is released for export when it is
      available to foreign nationals for visual inspection, when technology is
      exchanged orally or when technology is made available by practice or
      application under the guidance of persons with knowledge of the
      technology.</li></ul></p>

<p><a href=\"#The_Modelica_License_2-outline\">The Modelica License 2</a><br>
<a href=\"#Frequently_Asked_Questions-outline\">Frequently Asked Questions</a><br></p>

<hr>

<h4><a name=\"The_Modelica_License_2-outline\"></a>The Modelica License 2</h4>

<p><b>Preamble.</b> The goal of this license is that Modelica related
model libraries, software, images, documents, data files etc. can be
used freely in the original or a modified form, in open source and in
commercial environments (as long as the license conditions below are
fulfilled, in particular sections 2c) and 2d). The Original Work is
provided free of charge and the use is completely at your own risk.
Developers of free Modelica packages are encouraged to utilize this
license for their work.</p>

<p>The Modelica License applies to any Original Work that contains the
following licensing notice adjacent to the copyright notice(s) for
this Original Work:</p>

<p><b>Licensed by the Georgia Tech Research Corporation under the Modelica License 2</b></p>

<p><b>1. Definitions.</b></p>
<ol type=\"a\">
        <li>&ldquo;License&rdquo; is this Modelica License.</li>

        <li>&ldquo;Original Work&rdquo; is any work of authorship, including
        software, images, documents, data files, that contains the above
        licensing notice or that is packed together with a licensing notice
        referencing it.</li>

        <li>&ldquo;Licensor&rdquo; is the provider of the Original Work who has
        placed this licensing notice adjacent to the copyright notice(s) for
        the Original Work. The Original Work is either directly provided by
        the owner of the Original Work, or by a licensee of the owner.</li>

        <li>&ldquo;Derivative Work&rdquo; is any modification of the Original
        Work which represents, as a whole, an original work of authorship.
        For the matter of clarity and as examples:

        <ol  type=\"A\">
                <li>Derivative Work shall not include work that remains separable from
                the Original Work, as well as merely extracting a part of the
                Original Work without modifying it.</li>

                <li>Derivative Work shall not include (a) fixing of errors and/or (b)
                adding vendor-specific Modelica annotations and/or (c) using a
                subset of the classes of a Modelica package, and/or (d) using a
                different representation, e.g., a binary representation.</li>

                <li>Derivative Work shall include classes that are copied from the
                Original Work where declarations, equations or the documentation
                are modified.</li>

                <li>Derivative Work shall include executables to simulate the models
                that are generated by a Modelica translator based on the Original
                Work (of a Modelica package).</li>
        </ol>

        <li>&ldquo;Modified Work&rdquo; is any modification of the Original Work
        with the following exceptions: (a) fixing of errors and/or (b)
        adding vendor-specific Modelica annotations and/or (c) using a
        subset of the classes of a Modelica package, and/or (d) using a
        different representation, e.g., a binary representation.</li>

        <li>&quot;Source Code&quot; means the preferred form of the Original
        Work for making modifications to it and all available documentation
        describing how to modify the Original Work.</li>

        <li>&ldquo;You&rdquo; means an individual or a legal entity exercising
        rights under, and complying with all of the terms of, this License.</li>

        <li>&ldquo;Modelica package&rdquo; means any Modelica library that is
        defined with the &ldquo;<code><b>package</b>&nbsp;&lt;Name&gt;&nbsp;&hellip;&nbsp;<b>end</b>&nbsp;&lt;Name&gt;;</code>&rdquo; Modelica language element.</li>

</ol>

<p><b>2. Grant of Copyright License.</b> Licensor grants You a
worldwide, royalty-free, non-exclusive, sublicensable license, for
the duration of the copyright, to do the following:</p>

<ol type=\"a\">
        <li><p>To reproduce the Original Work in copies, either alone or as part of
        a collection.</p></li>
        <li><p>To create Derivative Works according to Section 1d) of this License.</p></li>
        <li><p>To distribute or communicate to the public copies of the <u>Original
        Work</u> or a <u>Derivative Work</u> under <u>this License</u>. No
        fee, neither as a copyright-license fee, nor as a selling fee for
        the copy as such may be charged under this License. Furthermore, a
        verbatim copy of this License must be included in any copy of the
        Original Work or a Derivative Work under this License.<br>
        For the matter of clarity, it is permitted A) to distribute or
        communicate such copies as part of a (possible commercial)
        collection where other parts are provided under different licenses
        and a license fee is charged for the other parts only and B) to
        charge for mere printing and shipping costs.</p></li>
        <li><p>To distribute or communicate to the public copies of a <u>Derivative
        Work</u>, alternatively to Section 2c), under <u>any other license</u>
        of your choice, especially also under a license for
        commercial/proprietary software, as long as You comply with Sections
        3, 4 and 8 below.<br>      For the matter of clarity, no
        restrictions regarding fees, either as to a copyright-license fee or
        as to a selling fee for the copy as such apply.</p></li>
        <li><p>To perform the Original Work publicly.</p></li>
        <li><p>To display the Original Work publicly.</p></li>
</ol>

<p><b>3. Acceptance.</b> Any use of the Original Work or a
Derivative Work, or any action according to either Section 2a) to 2f)
above constitutes Your acceptance of this License.</p>

<p><b>4. Designation of Derivative Works and of Modified Works.
</b>The identifying designation of Derivative Work and of Modified
Work must be different to the corresponding identifying designation
of the Original Work. This means especially that the (root-level)
name of a Modelica package under this license must be changed if the
package is modified (besides fixing of errors, adding vendor-specific
Modelica annotations, using a subset of the classes of a Modelica
package, or using another representation, e.g. a binary
representation).</p>

<p><b>5. Grant of Patent License.</b>
Licensor grants You a worldwide, royalty-free, non-exclusive, sublicensable license,
under patent claims owned by the Licensor or licensed to the Licensor by
the owners of the Original Work that are embodied in the Original Work
as furnished by the Licensor, for the duration of the patents,
to make, use, sell, offer for sale, have made, and import the Original Work
and Derivative Works under the conditions as given in Section 2.
For the matter of clarity, the license regarding Derivative Works covers
patent claims to the extent as they are embodied in the Original Work only.</p>

<p><b>6. Provision of Source Code.</b> Licensor agrees to provide
You with a copy of the Source Code of the Original Work but reserves
the right to decide freely on the manner of how the Original Work is
provided.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;For the matter of clarity, Licensor might provide only a binary
representation of the Original Work. In that case, You may (a) either
reproduce the Source Code from the binary representation if this is
possible (e.g., by performing a copy of an encrypted Modelica
package, if encryption allows the copy operation) or (b) request the
Source Code from the Licensor who will provide it to You.</p>

<p><b>7. Exclusions from License Grant.</b> Neither the names of
Licensor, nor the names of any contributors to the Original Work, nor
any of their trademarks or service marks, may be used to endorse or
promote products derived from this Original Work without express
prior permission of the Licensor. Except as otherwise expressly
stated in this License and in particular in Sections 2 and 5, nothing
in this License grants any license to Licensor&rsquo;s trademarks,
copyrights, patents, trade secrets or any other intellectual
property, and no patent license is granted to make, use, sell, offer
for sale, have made, or import embodiments of any patent claims.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;No license is granted to the trademarks of
Licensor even if such trademarks are included in the Original Work,
except as expressly stated in this License. Nothing in this License
shall be interpreted to prohibit Licensor from licensing under terms
different from this License any Original Work that Licensor otherwise
would have a right to license.</p>

<p><b>8. Attribution Rights.</b> You must retain in the Source
Code of the Original Work and of any Derivative Works that You
create, all author, copyright, patent, or trademark notices, as well
as any descriptive text identified therein as an &quot;Attribution
Notice&quot;. The same applies to the licensing notice of this
License in the Original Work. For the matter of clarity, &ldquo;author
notice&rdquo; means the notice that identifies the original
author(s).<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;You must cause the Source Code for any Derivative
Works that You create to carry a prominent Attribution Notice
reasonably calculated to inform recipients that You have modified the
Original Work.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;In case the Original Work or Derivative Work is not provided in
Source Code, the Attribution Notices shall be appropriately
displayed, e.g., in the documentation of the Derivative Work.</p>

<p><b>9. Disclaimer
of Warranty.<br><u>The Original Work is provided under this
License on an &quot;as is&quot; basis and without warranty, either
express or implied, including, without limitation, the warranties of
non-infringement, merchantability or fitness for a particular
purpose. The entire risk as to the quality of the Original Work is
with You.</b></u> This disclaimer of warranty constitutes an
essential part of this License. No license to the Original Work is
granted by this License except under this disclaimer.</p>

<p><b>10. Limitation of Liability.</b> Under no circumstances and
under no legal theory, whether in tort (including negligence),
contract, or otherwise, shall the Licensor, the owner or a licensee
of the Original Work be liable to anyone for any direct, indirect,
general, special, incidental, or consequential damages of any
character arising as a result of this License or the use of the
Original Work including, without limitation, damages for loss of
goodwill, work stoppage, computer failure or malfunction, or any and
all other commercial damages or losses. This limitation of liability
shall not apply to the extent applicable law prohibits such
limitation.</p>

<p><b>11. Termination.</b> This License conditions your rights to
undertake the activities listed in Section 2 and 5, including your
right to create Derivative Works based upon the Original Work, and
doing so without observing these terms and conditions is prohibited
by copyright law and international treaty. Nothing in this License is
intended to affect copyright exceptions and limitations. This License
shall terminate immediately and You may no longer exercise any of the
rights granted to You by this License upon your failure to observe
the conditions of this license.</p>

<p><b>12. Termination for Patent Action.</b> This License shall
terminate automatically and You may no longer exercise any of the
rights granted to You by this License as of the date You commence an
action, including a cross-claim or counterclaim, against Licensor,
any owners of the Original Work or any licensee alleging that the
Original Work infringes a patent. This termination provision shall
not apply for an action alleging patent infringement through
combinations of the Original Work under combination with other
software or hardware.</p>

<p><b>13. Jurisdiction.</b> Any action or suit relating to this
License may be brought only in the courts of a jurisdiction wherein
the Licensor resides and under the laws of that jurisdiction
excluding its conflict-of-law provisions. The application of the
United Nations Convention on Contracts for the International Sale of
Goods is expressly excluded. Any use of the Original Work outside the
scope of this License or after its termination shall be subject to
the requirements and penalties of copyright or patent law in the
appropriate jurisdiction. This section shall survive the termination
of this License.</p>

<p><b>14. Attorneys&rsquo; Fees.</b> In any action to enforce the
terms of this License or seeking damages relating thereto, the
prevailing party shall be entitled to recover its costs and expenses,
including, without limitation, reasonable attorneys' fees and costs
incurred in connection with such action, including any appeal of such
action. This section shall survive the termination of this License.</p>

<p><b>15. Miscellaneous.</b></p>
<ol type=\"a\">
        <li>If any
        provision of this License is held to be unenforceable, such
        provision shall be reformed only to the extent necessary to make it
        enforceable.</li>

        <li>No verbal
        ancillary agreements have been made. Changes and additions to this
        License must appear in writing to be valid. This also applies to
        changing the clause pertaining to written form.</li>

        <li>You may use the
        Original Work in all ways not otherwise restricted or conditioned by
        this License or by law, and Licensor promises not to interfere with
        or be responsible for such uses by You.</li>
</ol>

<hr>

<h4><a name=\"Frequently_Asked_Questions-outline\"></a>
Frequently Asked Questions</h4>
<p>This
section contains questions/answer to users and/or distributors of
Modelica packages and/or documents under Modelica License 2. Note,
the answers to the questions below are not a legal interpretation of
the Modelica License 2. In case of a conflict, the language of the
license shall prevail.</p>

<h6>Using or Distributing a Modelica <u>Package</u> under the Modelica License 2</h6>

<p><b>What are the main
differences to the previous version of the Modelica License?</b></p>
<ol>
        <li><p>Modelica License 1 is unclear whether the licensed Modelica package
        can be distributed under a different license. Version 2 explicitly
        allows that &ldquo;Derivative Work&rdquo; can be distributed under
        any license of Your choice, see examples in Section 1d) as to what
        qualifies as Derivative Work (so, version 2 is clearer).</p>
        <li><p>If You modify a Modelica package under Modelica License 2 (besides
        fixing of errors, adding vendor-specific Modelica annotations, using
        a subset of the classes of a Modelica package, or using another
        representation, e.g., a binary representation), you must rename the
        root-level name of the package for your distribution. In version 1
        you could keep the name (so, version 2 is more restrictive). The
        reason of this restriction is to reduce the risk that Modelica
        packages are available that have identical names, but different
        functionality.</p>
        <li><p>Modelica License 1 states that &ldquo;It is not allowed to charge a
        fee for the original version or a modified version of the software,
        besides a reasonable fee for distribution and support&rdquo;.
        Version 2 has a similar intention for all Original Work under
        <u>Modelica License 2</u> (to remain free of charge and open source)
        but states this more clearly as &ldquo;No fee, neither as a
        copyright-license fee, nor as a selling fee for the copy as such may
        be charged&rdquo;. Contrary to version 1, Modelica License 2 has no
        restrictions on fees for Derivative Work that is provided under a
        different license (so, version 2 is clearer and has fewer
        restrictions).</p>
        <li><p>Modelica License 2 introduces several useful provisions for the
        licensee (articles 5, 6, 12), and for the licensor (articles 7, 12,
        13, 14) that have no counter part in version 1.</p>
        <li><p>Modelica License 2 can be applied to all type of work, including
        documents, images and data files, contrary to version 1 that was
        dedicated for software only (so, version 2 is more general).</p>
</ol>

<p><b>Can I distribute a
Modelica package (under Modelica License 2) as part of my commercial
Modelica modeling and simulation environment?</b></p>

<p>Yes,
according to Section 2c). However, you are not allowed to charge a
fee for this part of your environment. Of course, you can charge for
your part of the environment.</p>

<p><b>Can I distribute a
Modelica package (under Modelica License 2) under a different
license?</b></p>

<p>No.
The license of an unmodified Modelica package cannot be changed
according to Sections 2c) and 2d). This means that you cannot <u>sell</u>
copies of it, any distribution has to be free of charge.</p>

<p><b>Can I distribute a
Modelica package (under Modelica License 2) under a different license
when I first encrypt the package?</b></p>

<p>No.
Merely encrypting a package does not qualify for Derivative Work and
therefore the encrypted package has to stay under Modelica License 2.</p>

<p><b>Can I distribute a
Modelica package (under Modelica License 2) under a different license
when I first add classes to the package?</b></p>

<p>No.
The package itself remains unmodified, i.e., it is Original Work, and
therefore the license for this part must remain under Modelica
License 2. The newly added classes can be, however, under a different
license.</p>

<p><b>Can
I copy a class out of a Modelica package (under Modelica License 2)
and include it <u>unmodified</u> in a Modelica package
under a <u>commercial/proprietary license</u>?</b></p>

<p>No,
according to article 2c). However, you can include model, block,
function, package, record and connector classes in your Modelica
package under <u>Modelica License 2</u>. This means that your
Modelica package could be under a commercial/proprietary license, but
one or more classes of it are under Modelica License 2.<br>Note, a
&ldquo;type&rdquo; class (e.g., type Angle = Real(unit=&rdquo;rad&rdquo;))
can be copied and included unmodified under a commercial/proprietary
license (for details, see the next question).</p>

<p><b>Can
I copy a type class or <u>part</u> of a model, block,
function, record, connector class, out of a Modelica package (under
Modelica License 2) and include it modified or unmodified in a
Modelica package under a <u>commercial/proprietary</u>
license</b></p>

<p>Yes,
according to article 2d), since this will in the end usually qualify
as Derivative Work. The reasoning is the following: A type class or
part of another class (e.g., an equation, a declaration, part of a
class description) cannot be utilized &ldquo;by its own&rdquo;. In
order to make this &ldquo;usable&rdquo;, you have to add additional
code in order that the class can be utilized. This is therefore
usually Derivative Work and Derivative Work can be provided under a
different license. Note, this only holds, if the additional code
introduced is sufficient to qualify for Derivative Work. Merely, just
copying a class and changing, say, one character in the documentation
of this class would be no Derivative Work and therefore the copied
code would have to stay under Modelica License 2.</p>

<p><b>Can
I copy a class out of a Modelica package (under Modelica License 2)
and include it in <u>modified</u> form in a
<u>commercial/proprietary</u> Modelica package?</b></p>

<p>Yes.
If the modification can be seen as a &ldquo;Derivative Work&rdquo;,
you can place it under your commercial/proprietary license. If the
modification does not qualify as &ldquo;Derivative Work&rdquo; (e.g.,
bug fixes, vendor-specific annotations), it must remain under
Modelica License 2. This means that your Modelica package could be
under a commercial/proprietary license, but one or more parts of it
are under Modelica License 2.</p>

<p><b>Can I distribute a
&ldquo;save total model&rdquo; under my commercial/proprietary
license, even if classes under Modelica License 2 are included?</b></p>

<p>Your
classes of the &ldquo;save total model&rdquo; can be distributed
under your commercial/proprietary license, but the classes under
Modelica License 2 must remain under Modelica License 2. This means
you can distribute a &ldquo;save total model&rdquo;, but some parts
might be under Modelica License 2.</p>

<p><b>Can I distribute a
Modelica package (under Modelica License 2) in encrypted form?</b></p>

<p>Yes.
Note, if the encryption does not allow &ldquo;copying&rdquo; of
classes (in to unencrypted Modelica source code), you have to send
the Modelica source code of this package to your customer, if he/she
wishes it, according to article&nbsp;6.</p>

<p><b>Can I distribute an
executable under my commercial/proprietary license, if the model from
which the executable is generated uses models from a Modelica package
under Modelica License 2?</b></p>

<p>Yes,
according to article 2d), since this is seen as Derivative Work. The
reasoning is the following: An executable allows the simulation of a
concrete model, whereas models from a Modelica package (without
pre-processing, translation, tool run-time library) are not able to
be simulated without tool support. By the processing of the tool and
by its run-time libraries, significant new functionality is added (a
model can be simulated whereas previously it could not be simulated)
and functionality available in the package is removed (e.g., to build
up a new model by dragging components of the package is no longer
possible with the executable).</p>

<p><b>Is my modification to
a Modelica package (under Modelica License 2) a Derivative Work?</b></p>

<p>It
is not possible to give a general answer to it. To be regarded as &quot;an
original work of authorship&quot;, a derivative work must be
different enough from the original or must contain a substantial
amount of new material. Making minor changes or additions of little
substance to a preexisting work will not qualify the work as a new
version for such purposes.</p>

<h6>Using or Distributing a Modelica <u>Document</u> under the Modelica License 2</h6>

<p>This
section is devoted especially for the following applications:</p>
<ol type=\"a\">
        <li><p>A Modelica tool extracts information out of a Modelica package and
        presents the result in form of a &ldquo;manual&rdquo; for this
        package in, e.g., html, doc, or pdf format.</p>
        <li><p>The Modelica language specification is a document defining the
        Modelica language. It will be licensed under Modelica License 2.</p>
        <li><p>Someone writes a book about the Modelica language and/or Modelica
        packages and uses information which is available in the Modelica
        language specification and/or the corresponding Modelica package.</p>
</ol>

<p><b>Can I sell a manual
that was basically derived by extracting information automatically
from a Modelica package under Modelica License 2 (e.g., a &ldquo;reference
guide&rdquo; of the Modelica Standard Library):</b></p>

<p>Yes.
Extracting information from a Modelica package, and providing it in a
human readable, suitable format, like html, doc or pdf format, where
the content is significantly modified (e.g. tables with interface
information are constructed from the declarations of the public
variables) qualifies as Derivative Work and there are no restrictions
to charge a fee for Derivative Work under alternative 2d).</p>

<p><b>Can
I copy a text passage out of a Modelica document (under Modelica
License 2) and use it <u>unmodified</u> in my document
(e.g. the Modelica syntax description in the Modelica Specification)?</b></p>

<p>Yes.
In case you distribute your document, the copied parts are still
under Modelica License 2 and you are not allowed to charge a license
fee for this part. You can, of course, charge a fee for the rest of
your document.</p>

<p><b>Can
I copy a text passage out of a Modelica document (under Modelica
License 2) and use it in <u>modified</u> form in my
document?</b></p>

<p>Yes,
the creation of Derivative Works is allowed. In case the content is
significantly modified this qualifies as Derivative Work and there
are no restrictions to charge a fee for Derivative Work under
alternative 2d).</p>

<p><b>Can I sell a printed
version of a Modelica document (under Modelica License 2), e.g., the
Modelica Language Specification?</b></p>

<p>No,
if you are not the copyright-holder, since article 2c) does not allow
a selling fee for a (in this case physical) copy. However, mere
printing and shipping costs may be recovered.</p>
</body>
</html>"));

    end ModelicaLicense2;
  annotation (preferedView="info", DocumentationClass=true);

  end UsersGuide;


annotation (
  preferedView="info",
  uses(Modelica(version="3.2")),
  Commands(executeCall=FCSys.Units.setup() "Re-initialize the units."),
  Documentation(info="<html>
    <p><a href=\"modelica://FCSys\">FCSys</a> is an free, open-source library of
    declarative, dynamic, and flexible models of proton exchange membrane
    fuel cells (PEMFCs) in the <a href = \"http://www.modelica.org/\">Modelica</a>
    language.  Chemical, electrical, fluid, and thermal
    phenomena are included.
    There are options to adjust the assumptions, spatial discretization
    and dimensionality (1D, 2D, or 3D), and the present chemical species and material
    phases.  The framework is generic and could be extended to other electrochemical
    devices like batteries.</p>

    <p>A fuel cell is similar to a battery except that the reactants
    are externally stored or drawn
    from the environment.  The electrochemical reactions of a PEMFC are:
    <table border=0 cellspacing=0 cellpadding=2 align=center style=\"margin-left: auto;
margin-right: auto;\" class=noBorder>
      <tr>
        <td align=right style=\"white-space:nowrap; text-align:right;\" class=noBorder>
          2(H<sub>2</sub>
        </td>
        <td align=center style=\"white-space:nowrap; text-align:center;\" class=noBorder>
          &#8652;
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
          &#8652;
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
          &#8652;
        </td>
        <td align=left style=\"white-space:nowrap;\" class=noBorder>
          2H<sub>2</sub>O
        </td>
        <td class=noBorder>
          (net)
        </td>
      </tr>
    </table></p>

    <p><a href=\"#Fig1\">Figure 1</a> shows the seven primary layers of a typical PEMFC.
    Fluid enters and exits the cell through channels in the flow plates (FPs).  It spreads through
    the gas diffusion diffusion layers (GDLs) and reacts in the catalyst layers (CLs).  The
    proton exchange membrane (PEM) prevents electronic transport; therefore, electrons must
    pass through an external load to sustain the net reaction.  As
    shown in <a href=\"#Fig2\">Figure 2</a>, a PEMFC model can be constructed from models
    of the same layers in <a href=\"modelica://FCSys\">FCSys</a>.
    The model is modular; the gas diffusion and catalyst layers could be combined,
    or microporous layers could be inserted.</p>

    <p align=center id=\"Fig1\"><img src=\"modelica://FCSys/Resources/Documentation/CellFlows.png\">
<br>Figure 1: Layers and primary flows of a PEMFC.</p>

    <!--<p align=center id=\"Fig2\"><img src=\"modelica://FCSys/help/FCSys.Assemblies.Cells.CellD.png\" width=600>-->
    <p align=center id=\"Fig2\"><a href=\"modelica://FCSys.Assemblies.Cells.Cell\"><img src=\"modelica://FCSys/Resources/Documentation/FCSys.Assemblies.Cells.CellD.png\"></a>
<br>Figure 2: Diagram of the <a href=\"modelica://FCSys.Assemblies.Cells.Cell\">PEMFC model</a>.</p>

    <p>The models are primarily based on first principles&mdash;the advection, diffusion, and storage of
    material, momentum, and energy.  Upstream
    discretization is applied in a manner that reduces to pure
    diffusion (e.g., Fick's law, Newton's law of viscosity, or Fourier's law) when the bulk velocity is zero.
    The transport equations do
    not use the <a href = \"http://www.modelica.org/\">Modelica</a>
    <code>stream</code> operator since both advection and diffusion are
    important in fuel cells.</p>

    <p>Each layer may be subdivided in any direction for increased fidelity.  Regions may
    be directly connected without producing nonlinear systems of
    equations, and species may be independently included in each region.  This is
    different than
    <a href=\"modelica://Modelica.Media\">Modelica.Media</a>,
    where each media model contains a predefined set of species.</p>

    <p>A cell may be simulated under specified boundary conditions or connected to
    <a href=\"modelica://Modelica.Fluid\">Modelica.Fluid</a> components using
    <a href=\"modelica://FCSys.Conditions.Adapters\">adapters</a> that are included.
    <a href=\"#Fig3\">Figure 3</a> shows a series of polarization curves generated
    from the <a href=\"modelica://FCSys.Assemblies.Cells.Examples.CellProfile\">CellProfile</a> model **Update with new plot (currently old result).
    Please see the <a href=\"modelica://FCSys.UsersGuide.SampleResults\">sample results</a> for more plots.</p>

    <p align=center id=\"Fig3\"><a href=\"modelica://FCSys.Assemblies.Cells.Examples.CellProfile\"><img src=\"modelica://FCSys/Resources/Documentation/Polarization.png\"></a>
<br>Figure 3: Polarization curves under various conditions.</p>

    <p><b>Licensed by the Georgia Tech Research Corporation under the Modelica License 2</b>
<br>Copyright 2007&ndash;2013, Georgia Tech Research Corporation.</p>

    <p><i>This Modelica package is <u>free</u> software and the use is completely
    at <u>your own risk</u>; it can be redistributed and/or modified under the
    terms of the Modelica License 2. For license conditions (including the
    disclaimer of warranty) see
    <a href=\"modelica://FCSys.UsersGuide.ModelicaLicense2\">
    FCSys.UsersGuide.ModelicaLicense2</a>
    or visit <a href=\"http://www.modelica.org/licenses/ModelicaLicense2\">
    http://www.modelica.org/licenses/ModelicaLicense2</a>.</i></p></html>",
      revisions="<html>
    <ul>
    <li><a href=\"mailto:kdavies4@gmail.com\">Kevin Davies</a>, 8/1/2007:<br>Version 1.0 (as <code>PEMFCSim</code>)</li>
    <li><a href=\"mailto:kdavies4@gmail.com\">Kevin Davies</a>, 8/20/2008:<br>Version 1.1 (to model CO contamination for <a href=\"http://www.hnei.hawaii.edu\">HNEI</a>)</li>
    <li><a href=\"mailto:kdavies4@gmail.com\">Kevin Davies</a>, TBD/2013:<br>Version 2.0.0-alpha (released in the public domain)</li>
    </ul>
    </html>"),
  Icon(graphics={
      Polygon(
        points={{-14,42},{-24,32},{-4,32},{6,42},{-14,42}},
        lineColor={0,0,0},
        fillColor={255,255,255},
        fillPattern=FillPattern.Solid),
      Polygon(
        points={{-40,42},{-50,32},{-30,32},{-20,42},{-40,42}},
        lineColor={0,0,0},
        fillPattern=FillPattern.HorizontalCylinder,
        fillColor={0,192,0}),
      Polygon(
        points={{-20,42},{-30,32},{-24,32},{-14,42},{-20,42}},
        lineColor={0,0,0},
        fillPattern=FillPattern.HorizontalCylinder,
        fillColor={0,0,0}),
      Rectangle(
        extent={{-4,32},{2,-62}},
        fillPattern=FillPattern.Solid,
        fillColor={0,0,0},
        pattern=LinePattern.None),
      Polygon(
        points={{6,42},{-4,32},{2,32},{12,42},{6,42}},
        lineColor={0,0,0},
        fillPattern=FillPattern.HorizontalCylinder,
        fillColor={0,0,0}),
      Line(
        points={{-50,32},{-50,-62}},
        pattern=LinePattern.None,
        smooth=Smooth.None),
      Polygon(
        points={{-56,54},{-76,34},{-56,34},{-36,54},{-56,54}},
        lineColor={0,0,0},
        fillColor={135,135,135},
        fillPattern=FillPattern.Solid),
      Rectangle(
        extent={{-39.6277,31.7996},{-67.912,17.6573}},
        lineColor={0,0,0},
        fillPattern=FillPattern.HorizontalCylinder,
        rotation=45,
        fillColor={255,255,255},
        origin={46.5067,57.5353}),
      Rectangle(
        extent={{-24,32},{-4,-62}},
        lineColor={0,0,0},
        fillPattern=FillPattern.VerticalCylinder,
        fillColor={255,255,255}),
      Line(points={{-40,42},{22,42}}, color={0,0,0}),
      Rectangle(
        extent={{-5.21738,-5.21961},{-33.5017,-33.5041}},
        lineColor={0,0,170},
        fillPattern=FillPattern.VerticalCylinder,
        rotation=45,
        fillColor={0,0,240},
        origin={21.9983,59.3803}),
      Rectangle(
        extent={{2,32},{42,-62}},
        lineColor={0,0,170},
        fillPattern=FillPattern.VerticalCylinder,
        fillColor={0,0,240}),
      Polygon(
        points={{-36,54},{-56,34},{-56,-74},{-36,-54},{-36,54}},
        lineColor={0,0,0},
        fillColor={95,95,95},
        fillPattern=FillPattern.Solid),
      Rectangle(
        extent={{-5.21774,-5.2196},{-33.502,-33.5042}},
        lineColor={196,11,40},
        fillPattern=FillPattern.HorizontalCylinder,
        rotation=45,
        fillColor={253,52,56},
        origin={-40.001,69.3803}),
      Rectangle(
        extent={{-70,32},{-30,-62}},
        lineColor={196,11,40},
        fillPattern=FillPattern.VerticalCylinder,
        fillColor={253,52,56}),
      Rectangle(
        extent={{-70,32},{-50,-64}},
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
        origin={55.0166,71.3801}),
      Rectangle(
        extent={{-76,34},{-56,-74}},
        lineColor={95,95,95},
        fillPattern=FillPattern.VerticalCylinder,
        fillColor={135,135,135}),
      Polygon(
        points={{36,54},{24,42},{-36,42},{-36,54},{36,54}},
        smooth=Smooth.None,
        fillPattern=FillPattern.Solid,
        fillColor={240,240,240},
        pattern=LinePattern.None,
        lineColor={0,0,0}),
      Rectangle(
        extent={{-76.648,66.211},{-119.073,52.0689}},
        lineColor={95,95,95},
        fillPattern=FillPattern.HorizontalCylinder,
        rotation=45,
        fillColor={135,135,135},
        origin={147.017,71.3801}),
      Rectangle(
        extent={{16,34},{36,-74}},
        lineColor={95,95,95},
        fillPattern=FillPattern.VerticalCylinder,
        fillColor={135,135,135}),
      Polygon(
        points={{-36,54},{-36,42},{-40,42},{-40,50},{-36,54}},
        smooth=Smooth.None,
        fillColor={95,95,95},
        fillPattern=FillPattern.Solid,
        pattern=LinePattern.None,
        lineColor={0,0,0}),
      Ellipse(
        extent={{-54,52},{-46,48}},
        lineColor={135,135,135},
        fillColor={253,52,56},
        fillPattern=FillPattern.Sphere),
      Ellipse(
        extent={{26,40},{34,36}},
        lineColor={135,135,135},
        fillColor={0,0,240},
        fillPattern=FillPattern.Sphere),
      Polygon(
        points={{-36,54},{-36,42},{-40,42},{-50,32},{-56,34},{-36,54}},
        smooth=Smooth.None,
        fillColor={95,95,95},
        fillPattern=FillPattern.Solid,
        pattern=LinePattern.None,
        lineColor={0,0,0}),
      Line(
        points={{-40,42},{-50,32}},
        color={0,0,0},
        smooth=Smooth.None),
      Polygon(
        points={{56,54},{36,34},{36,-74},{56,-54},{56,54}},
        lineColor={0,0,0},
        fillColor={95,95,95},
        fillPattern=FillPattern.Solid),
      Rectangle(extent={{16,34},{36,-74}}, lineColor={0,0,0}),
      Rectangle(extent={{-76,34},{-56,-74}}, lineColor={0,0,0}),
      Line(
        points={{-36,54},{-36,42}},
        color={0,0,0},
        smooth=Smooth.None),
      Line(points={{-40,42},{24,42}}, color={0,0,0}),
      Rectangle(
        extent={{-56,64},{56,54}},
        pattern=LinePattern.None,
        fillColor={240,240,240},
        fillPattern=FillPattern.Solid,
        lineColor={0,0,0}),
      Polygon(
        points={{-56,54},{-36,54},{-56,34},{-76,34},{-56,54}},
        lineColor={0,0,0},
        smooth=Smooth.None),
      Polygon(
        points={{36,54},{56,54},{36,34},{16,34},{36,54}},
        lineColor={0,0,0},
        smooth=Smooth.None),
      Rectangle(extent={{-50,32},{16,-62}}, lineColor={0,0,0}),
      Rectangle(
        extent={{-30,32},{-24,-62}},
        fillPattern=FillPattern.Solid,
        fillColor={0,0,0},
        pattern=LinePattern.None)}),
  version=",
  dateModified="",
  revisionID="");
end FCSys;

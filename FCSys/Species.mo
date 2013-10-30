within FCSys;
package Species "Dynamic models of chemical species"
  extends Modelica.Icons.Package;
  package 'C+' "C"
    extends Modelica.Icons.Package;
    package Graphite "<html>C<sup>+</sup> graphite</html>"
      extends Modelica.Icons.Package;

      model Correlated "Correlated properties"
        extends Solid(redeclare replaceable package Data =
              Characteristics.'C+'.Graphite);

        // TODO: Update this to pull properties and settings from Fixed.
        // Do the same for other species.
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="'C+'",
          Documentation(info=
                "<html><p>Please see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"),

          Icon(graphics));

      end Correlated;

      model Fixed "Fixed properties"

        extends Solid(
          redeclare replaceable package Data = Characteristics.'C+'.Graphite (
              n_c=0,
              T_lim_c={0,Modelica.Constants.inf},
              b_c=[935*U.J*Data.m/(U.kg*U.K)],
              B_c=[Data.Deltah0_f - (935*U.J*Data.m/U.kg)*298.15, 154.663*U.J/(
                  U.mol*U.K) - Data.b_c[1, 1]*log(298.15*U.K)]),
          final k_intra,
          final mu=0,
          redeclare parameter Q.TimeAbsolute nu=Data.nu(),
          redeclare parameter Q.ResistivityThermal theta=U.m*U.K/(11.1*U.W));

        // Note:  Parameter expressions (e.g., nu=Data.nu(environment.T)) are not
        // used here since they would render the parameters unadjustable in Dymola
        // 7.4.  This also applies to the other species.

        // See the documentation layer for a table of values for the specific heat
        // capacity and thermal resistivity.
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="'C+'",
          Documentation(info="<html><p>Assumptions:<ol>
    <li>The specific heat capacity is fixed (independent of temperature).</li>
    <li>The thermal independity and thermal resistivity are fixed (e.g., independent of temperature).</li>
    <li>Mobility is zero.</li>
    </ol></p>

   <p>The default isobaric specific heat capacity (<i>b<sub>c</sub></i> = <code>[935*U.J*Data.m/(U.kg*U.K)]</code>)
   and thermal
   resistivity (&theta; = <code>U.m*U.K/(11.1*U.W)</code>) are for graphite fiber epoxy (25% vol)
   composite (with heat flow parallel to the fibers) at 300&nbsp;K
   [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, p. 909].
   The integration offset for specific entropy is set such that
   the specific entropy is 154.663&nbsp;J/(mol&middot;K) at 25&nbsp;&deg;C and <i>p</i><sup>o</sup> (1&nbsp;atm).
   This is the value from Table B in [<a href=\"modelica://FCSys.UsersGuide.References\">McBride2002</a>].
   Additional thermal data is listed in <a href=\"#Tab1\">Table 1</a>.</p>

  <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
    <caption align=\"top\" id=\"Tab1\">Table 1: Properties of forms of C [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, p. 909].</caption>
    <tr>
      <th rowspan=3 valign=\"middle\"><i>T</i><br><code>/U.K</code></th>
      <th rowspan=1 colspan=2 width=1 valign=\"middle\">Diamond (type IIa)</th>
      <th rowspan=1 colspan=1 width=1 valign=\"middle\">Amorphous<br>carbon</th>
      <th rowspan=1 colspan=3 width=1 valign=\"middle\">Graphite (pyrolytic)</th>
      <th rowspan=1 colspan=3 width=1>Graphite fiber epoxy (25% vol)<br>composite</th>
    </tr>
    <tr>
      <th rowspan=2 valign=\"middle\"><i>c<sub>p</sub></i><code>*U.kg<br>*U.K<br>/(U.J<br>*Data.m)</code></th>
      <th rowspan=2 valign=\"middle\">&theta;<code><br>*U.W<br>/(U.m<br>*U.K)</code></th>
      <th rowspan=2 valign=\"middle\">&theta;<code><br>*U.W<br>/(U.m<br>*U.K)</code></th>
      <th rowspan=2 valign=\"middle\"><i>c<sub>p</sub></i><code>*U.kg<br>*U.W<br>/(U.J<br>*Data.m)</code></th>
      <th rowspan=1 colspan=2>&theta;<code>*U.W/(U.m*U.K)</code></th>
      <th rowspan=2 valign=\"middle\"><i>c<sub>p</sub></i><code>*U.kg<br>*U.K<br>/(U.J<br>*Data.m)</code></th>
      <th rowspan=1 colspan=2>&theta;<code>*U.W/(U.m*U.K)</code></th>
    </tr>
    <tr>
      <th valign=\"middle\">Parallel<br>to layers</th>
      <th valign=\"middle\">Perpendicular<br>to layers</th>
      <th valign=\"middle\">Parallel<br>to layers</th>
      <th valign=\"middle\">Perpendicular<br>to layers</th>
    </tr>
<tr><td>100</td><td>21</td><td>1/10000</td><td>1/0.67</td><td>136</td><td>1/4970</td><td>1/16.8</td><td>337</td><td>1/5.7</td><td>1/0.46</td></tr>
<tr><td>200</td><td>194</td><td>1/4000</td><td>1/1.18</td><td>411</td><td>1/3230</td><td>1/9.23</td><td>642</td><td>1/8.7</td><td>1/0.68</td></tr>
<tr><td>300</td><td>509</td><td>1/2300</td><td>1/1.89</td><td>709</td><td>1/1950</td><td>1/5.70</td><td>935</td><td>1/11.1</td><td>1/0.87</td></tr>
<tr><td>400</td><td>853</td><td>1/1540</td><td>1/2.19</td><td>992</td><td>1/1390</td><td>1/4.09</td><td>1216</td><td>1/13.0</td><td>1/1.1</td></tr>
<tr><td>600</td><td>-</td><td>-</td><td>1/2.37</td><td>1406</td><td>1/892</td><td>1/2.68</td><td>-</td><td>-</td><td>-</td></tr>
<tr><td>800</td><td>-</td><td>-</td><td>1/2.53</td><td>1650</td><td>1/667</td><td>1/2.01</td><td>-</td><td>-</td><td>-</td></tr>
<tr><td>1000</td><td>-</td><td>-</td><td>1/2.84</td><td>1793</td><td>1/534</td><td>1/1.60</td><td>-</td><td>-</td><td>-</td></tr>
<tr><td>1200</td><td>-</td><td>-</td><td>1/3.48</td><td>1890</td><td>1/448</td><td>1/1.34</td><td>-</td><td>-</td><td>-</td></tr>
<tr><td>1500</td><td>-</td><td>-</td><td>-</td><td>1974</td><td>1/357</td><td>1/1.08</td><td>-</td><td>-</td><td>-</td></tr>
<tr><td>2000</td><td>-</td><td>-</td><td>-</td><td>2043</td><td>1/262</td><td>1/0.81</td><td>-</td><td>-</td><td>-</td></tr>
  </table>

  <p>For more information, please see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"));

      end Fixed;

    end Graphite;

  end 'C+';

  package 'SO3-'
    "<html>C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S<sup>-</sup> (abbreviated as SO<sub>3</sub><sup>-</sup>)</html>"
    extends Modelica.Icons.Package;
    package Ionomer
      "<html>C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S<sup>-</sup> ionomer</html>"
      extends Modelica.Icons.Package;

      model Correlated "Correlated properties"
        extends Solid(redeclare replaceable package Data =
              Characteristics.'SO3-'.Ionomer);
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="'SO3-'",
          Documentation(info=
                "<html><p>Please see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"));

      end Correlated;

      model Fixed "Fixed properties"
        extends Solid(
          redeclare replaceable package Data = Characteristics.'SO3-'.Ionomer,
          redeclare parameter Q.Mobility mu=Data.mu(),
          redeclare parameter Q.TimeAbsolute nu=Data.nu(),
          redeclare parameter Q.ResistivityThermal theta=U.m*U.K/(0.16*U.W));

        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="'SO3-'",
          Documentation(info="<html><p>Assumptions:
    <ol>
    <li>The thermal independity and thermal resistivity are fixed (e.g., independent of temperature)</li>
    </ol></p>

    <p>The default thermal resistivity (&theta; = <code>U.m*U.K/(0.16*U.W)</code>) is of dry
  Nafion 115 [<a href=\"modelica://FCSys.UsersGuide.References\">Kandlikar2009</a>, p. 1277].</p>

<p>For more information, please see the
    <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"));

      end Fixed;

    end Ionomer;

  end 'SO3-';

  package 'e-' "<html>e<sup>-</sup></html>"
    extends Modelica.Icons.Package;
    package Graphite "<html>e<sup>-</sup> in graphite</html>"
      extends Modelica.Icons.Package;

      model Fixed "Fixed properties"
        extends Ion(
          redeclare final package Data = Characteristics.'e-'.Graphite,
          final Nu_Phi,
          final Nu_Q,
          final consRot,
          final upstreamX=false,
          final upstreamY=false,
          final upstreamZ=false,
          final k_intra,
          final consMaterial=ConsThermo.steady,
          final consTransX=ConsMom.steady,
          final consTransY=ConsMom.steady,
          final consTransZ=ConsMom.steady,
          final consEnergy=ConsThermo.steady,
          final initMaterial=Init.none,
          final initEnergy=Init.none,
          final N_IC,
          final p_IC,
          final h_IC,
          final V_IC,
          final rho_IC,
          final g_IC,
          final T_IC,
          final nu=1,
          final tauprime=V/(v*Io),
          final theta=Modelica.Constants.inf);

        parameter Q.Current Io(min=0) = 10*U.A "Exchange current" annotation (
            Dialog(__Dymola_label="<html><i>I</i><sup>o</sup></html>"));

        output Q.Potential wprime(stateSelect=StateSelect.never) =
          electrochemical.w - g if environment.analysis
          "Reaction overpotential";

        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="'e-'",
          Documentation(info="<html>

    <p>Assumptions:<ol>
          <li>The thermal resistivity is infinite.  All of the thermal conductance is attributed to 
          the substrate
          (e.g., <a href=\"modelica://FCSys.Species.'C+'.Graphite\">C+</a>).<li>
          <li>The exchange current density (<i>I</i><sup>o</sup>) is mapped to the reaction interval (&tau;&prime;) of electrons assuming that 
          the reaction interval is zero for other species involved in the reaction.</li>
          <li>The conductivity is mapped to the mobility of the electrons by assuming that
          the mobility of the substrate (e.g., 
          <a href=\"modelica://FCSys.Species.'C+'.Graphite\">C+</a>) is zero.</li>
    </ol></p>

    <p>For more information, please see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"),

          Diagram(graphics),
          Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                  {100,100}}), graphics));

      end Fixed;

    end Graphite;

  end 'e-';

  package 'H+' "<html>H<sup>+</sup></html>"
    extends Modelica.Icons.Package;
    package Ionomer "<html>H<sup>+</sup> in ionomer</html>"
      extends Modelica.Icons.Package;

      model Fixed "Fixed properties"
        extends Ion(
          redeclare replaceable package Data = Characteristics.'H+'.Ionomer,
          final initEnergy=Init.none,
          redeclare parameter Q.ResistivityThermal theta=U.m*U.K/(0.1661*U.W),
          sigma=0.083*U.S/U.cm,
          final Nu_Phi,
          final Nu_Q,
          final consRot,
          final upstreamX=false,
          final upstreamY=false,
          final upstreamZ=false,
          final alpha,
          final consMaterial=ConsThermo.steady,
          final initMaterial=Init.none,
          final N_IC,
          final p_IC,
          final h_IC,
          final V_IC,
          final rho_IC,
          final g_IC,
          final T_IC,
          final nu=1,
          final tauprime=0,
          final consEnergy=ConsThermo.steady);

        // See the documentation for a table of values.
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="'H+'",
          Documentation(info="<html>  
<p>Assumptions:<ol>
    <li>The generalized resistivities (&eta;, &theta;) are fixed (e.g., independent of temperature).</li>
    <li>The electrochemical reaction rate is governed by the electrons.  Therefore, the reaction interval (&tau;&prime;) is zero for protons.</li>
              <li>The conductivity is mapped to the mobility of the protons by assuming that
          the mobility of the substrate (e.g., 
          <a href=\"modelica://FCSys.Species.'SO3-'.Ionomer\">C19HF37O5S-</a>) is zero.</li>
    </ol></p>

<p>The default conductivity (&sigma; = <code>0.083*U.S/U.cm</code>)
  is for DuPont<sup>TM</sup> Nafion&reg; N-112 [<a href=\"modelica://FCSys.Regions.PEMs.DuPontN112\">DuPontN112</a>].</p>

  <p>The default thermal resistivity (&theta; = <code>U.m*U.K/(0.1661*U.W)</code>) is of H gas
  (rather than H<sup>+</sup>) at 300&nbsp;K from [<a href=\"modelica://FCSys.UsersGuide.References\">Schetz1996</a>, p. 139].
  <a href=\"#Tab1\">Table 1</a> lists the properties at other temperatures.</p>

    <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
  <caption align=\"top\" id=\"Tab1\">Table 1: Properties of H gas (not H<sup>+</sup>) [<a href=\"modelica://FCSys.UsersGuide.References\">Schetz1996</a>, p. 139]</caption>
<tr>
      <th valign=\"middle\"><i>T</i><br><code>/U.K</code></th>
      <th width=1>&eta;<code><br>*U.Pa*U.s</code></th>
      <th width=1>&theta;<code>*U.W<br>/(U.m*U.K)</code></th>
    </tr>
<tr><td>200</td><td>1/3.8e-6</td><td>1/0.1197</td></tr>
<tr><td>300</td><td>1/5.3e-6</td><td>1/0.1661</td></tr>
<tr><td>400</td><td>1/6.7e-6</td><td>1/0.2094</td></tr>
<tr><td>500</td><td>1/8.1e-6</td><td>1/0.2507</td></tr>
<tr><td>600</td><td>1/9.3e-6</td><td>1/0.2906</td></tr>
<tr><td>700</td><td>1/10.6e-6</td><td>1/0.3292</td></tr>
<tr><td>800</td><td>1/11.8e-6</td><td>1/0.3670</td></tr>
<tr><td>900</td><td>1/13.0e-6</td><td>1/0.4040</td></tr>
<tr><td>1000</td><td>1/14.2e-6</td><td>1/0.4403</td></tr>
<tr><td>1100</td><td>1/15.3e-6</td><td>1/0.4761</td></tr>
<tr><td>1200</td><td>1/16.5e-6</td><td>1/0.5114</td></tr>
<tr><td>1300</td><td>1/17.6e-6</td><td>1/0.5462</td></tr>
<tr><td>1400</td><td>1/18.7e-6</td><td>1/0.5807</td></tr>
<tr><td>1500</td><td>1/19.8e-6</td><td>1/0.6149</td></tr>
<tr><td>1600</td><td>1/20.9e-6</td><td>1/0.6487</td></tr>
<tr><td>1700</td><td>1/22.0e-6</td><td>1/0.6823</td></tr>
<tr><td>1800</td><td>1/23.1e-6</td><td>1/0.7156</td></tr>
<tr><td>1900</td><td>1/24.2e-6</td><td>1/0.7488</td></tr>
<tr><td>2000</td><td>1/25.2e-6</td><td>1/0.7817</td></tr>
<tr><td>2100</td><td>1/26.3e-6</td><td>1/0.8144</td></tr>
<tr><td>2200</td><td>1/27.3e-6</td><td>1/0.8470</td></tr>
<tr><td>2300</td><td>1/28.4e-6</td><td>1/0.8794</td></tr>
<tr><td>2400</td><td>1/29.4e-6</td><td>1/0.9117</td></tr>
<tr><td>2500</td><td>1/30.5e-6</td><td>1/0.9438</td></tr>
<tr><td>2600</td><td>1/31.5e-6</td><td>1/0.9758</td></tr>
<tr><td>2700</td><td>1/32.5e-6</td><td>1/1.0077</td></tr>
<tr><td>2800</td><td>1/33.6e-6</td><td>1/1.0395</td></tr>
<tr><td>2900</td><td>1/34.6e-6</td><td>1/1.0711</td></tr>
<tr><td>3000</td><td>1/35.6e-6</td><td>1/1.1027</td></tr>
<tr><td>3100</td><td>1/36.6e-6</td><td>1/1.1347</td></tr>
<tr><td>3200</td><td>1/37.7e-6</td><td>1/1.1664</td></tr>
<tr><td>3300</td><td>1/38.7e-6</td><td>1/1.1978</td></tr>
<tr><td>3400</td><td>1/39.7e-6</td><td>1/1.2288</td></tr>
<tr><td>3500</td><td>1/40.7e-6</td><td>1/1.2592</td></tr>
<tr><td>3600</td><td>1/41.6e-6</td><td>1/1.2884</td></tr>
<tr><td>3700</td><td>1/42.5e-6</td><td>1/1.3171</td></tr>
<tr><td>3800</td><td>1/43.4e-6</td><td>1/1.3455</td></tr>
<tr><td>3900</td><td>1/44.4e-6</td><td>1/1.3735</td></tr>
<tr><td>4000</td><td>1/45.2e-6</td><td>1/1.4012</td></tr>
<tr><td>4100</td><td>1/46.1e-6</td><td>1/1.4290</td></tr>
<tr><td>4200</td><td>1/47.0e-6</td><td>1/1.4566</td></tr>
<tr><td>4300</td><td>1/47.9e-6</td><td>1/1.4842</td></tr>
<tr><td>4400</td><td>1/48.8e-6</td><td>1/1.5116</td></tr>
<tr><td>4500</td><td>1/49.7e-6</td><td>1/1.5389</td></tr>
<tr><td>4600</td><td>1/50.6e-6</td><td>1/1.5661</td></tr>
<tr><td>4700</td><td>1/51.5e-6</td><td>1/1.5933</td></tr>
<tr><td>4800</td><td>1/52.3e-6</td><td>1/1.6204</td></tr>
<tr><td>4900</td><td>1/53.2e-6</td><td>1/1.6477</td></tr>
<tr><td>5000</td><td>1/54.1e-6</td><td>1/1.6750</td></tr>
  </table></p>

<p>For more information, please see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"));

      end Fixed;

    end Ionomer;

  end 'H+';

  package H2 "<html>H<sub>2</sub></html>"
    extends Modelica.Icons.Package;
    package Gas "<html>H<sub>2</sub> gas</html>"
      extends Modelica.Icons.Package;

      model Correlated "Correlated properties"
        extends Fluid(redeclare replaceable package Data =
              Characteristics.H2.Gas (b_v=[1], n_v={-1,0}));
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="H2",
          Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol></p>

<p>For more information, please see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"));

      end Correlated;

      model Fixed "Fixed properties"
        extends Fluid(
          redeclare replaceable package Data = FCSys.Characteristics.H2.Gas (
                b_v=[1], n_v={-1,0}),
          final tauprime,
          redeclare parameter Q.Mobility mu=Data.mu(),
          redeclare parameter Q.TimeAbsolute nu=Data.nu(),
          redeclare parameter Q.Fluidity eta=1/(8.96e-6*U.Pa*U.s),
          redeclare parameter Q.ResistivityThermal theta=U.m*U.K/(0.183*U.W),
          final k_intra,
          final alpha);

        // See the documentation for a table of values.
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="H2",
          Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
        <li>The generalized resistivities (&eta;, &theta;) are fixed (e.g., independent of temperature).</li>

    </ol></p>

<p>Additional notes:<ul>
<li>
  The specific heat capacity is not fixed because it would
  affect the chemical potential and result in an incorrect cell
  potential.</li></ul></p>

<p>The default resistivities (&eta; = <code>1/(89.6e-7*U.Pa*U.s)</code>
and &theta; = <code>U.m*U.K/(183e-3*U.W)</code>) are based on data of H<sub>2</sub> gas at 1&nbsp;atm and
  300&nbsp;K from Incropera and DeWitt [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 919&ndash;920].
  <a href=\"#Tab1\">Table 1</a> lists the properties at other temperatures.</p>

  <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
    <caption align=\"top\" id=\"Tab1\">Table 1: Properties of H<sub>2</sub> gas at 1&nbsp;atm [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 919&ndash;920].</caption>
  <tr>
      <th valign=\"middle\"><i>T</i><br><code>/U.K</code></th>
      <th width=1><i>c<sub>p</sub></i><code>*U.kg*U.K<br>/(U.J*Data.m)</code></th>
      <th width=1>&eta;<code><br>*U.Pa*U.s</code></th>
      <th width=1>&theta;<code>*U.W<br>/(U.m*U.K)</code></th>
    </tr>
<tr><td>100</td><td>11.23e3</td><td>1/42.1e-7</td><td>1/67.0e-3</td></tr>
<tr><td>150</td><td>12.60e3</td><td>1/56.0e-7</td><td>1/101e-3</td></tr>
<tr><td>200</td><td>13.54e3</td><td>1/68.1e-7</td><td>1/131e-3</td></tr>
<tr><td>250</td><td>14.06e3</td><td>1/78.9e-7</td><td>1/157e-3</td></tr>
<tr><td>300</td><td>14.31e3</td><td>1/89.6e-7</td><td>1/183e-3</td></tr>
<tr><td>350</td><td>14.43e3</td><td>1/98.8e-7</td><td>1/204e-3</td></tr>
<tr><td>400</td><td>14.48e3</td><td>1/108.2e-7</td><td>1/226e-3</td></tr>
<tr><td>450</td><td>14.50e3</td><td>1/117.2e-7</td><td>1/247e-3</td></tr>
<tr><td>500</td><td>14.52e3</td><td>1/126.4e-7</td><td>1/266e-3</td></tr>
<tr><td>550</td><td>14.53e3</td><td>1/134.3e-7</td><td>1/285e-3</td></tr>
<tr><td>600</td><td>14.55e3</td><td>1/142.4e-7</td><td>1/305e-3</td></tr>
<tr><td>700</td><td>14.61e3</td><td>1/157.8e-7</td><td>1/342e-3</td></tr>
<tr><td>800</td><td>14.70e3</td><td>1/172.4e-7</td><td>1/378e-3</td></tr>
<tr><td>900</td><td>14.83e3</td><td>1/186.5e-7</td><td>1/412e-3</td></tr>
<tr><td>1000</td><td>14.99e3</td><td>1/201.3e-7</td><td>1/448e-3</td></tr>
<tr><td>1100</td><td>15.17e3</td><td>1/213.0e-7</td><td>1/488e-3</td></tr>
<tr><td>1200</td><td>15.37e3</td><td>1/226.2e-7</td><td>1/528e-3</td></tr>
<tr><td>1300</td><td>15.59e3</td><td>1/238.5e-7</td><td>1/568e-3</td></tr>
<tr><td>1400</td><td>15.81e3</td><td>1/250.7e-7</td><td>1/610e-3</td></tr>
<tr><td>1500</td><td>16.02e3</td><td>1/262.7e-7</td><td>1/655e-3</td></tr>
<tr><td>1600</td><td>16.28e3</td><td>1/273.7e-7</td><td>1/697e-3</td></tr>
<tr><td>1700</td><td>16.58e3</td><td>1/284.9e-7</td><td>1/742e-3</td></tr>
<tr><td>1800</td><td>16.96e3</td><td>1/296.1e-7</td><td>1/786e-3</td></tr>
<tr><td>1900</td><td>17.49e3</td><td>1/307.2e-7</td><td>1/835e-3</td></tr>
<tr><td>2000</td><td>18.25e3</td><td>1/318.2e-7</td><td>1/878e-3</td></tr>
    </tr>
  </table>
<p>For more information, please see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"),

          Diagram(graphics),
          Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                  {100,100}}),graphics));

      end Fixed;

    end Gas;

  end H2;

  package H2O "<html>H<sub>2</sub>O</html>"
    extends Modelica.Icons.Package;
    package Gas "<html>H<sub>2</sub>O gas</html>"
      extends Modelica.Icons.Package;

      model Correlated "Correlated properties"

        extends Fluid(redeclare replaceable package Data =
              Characteristics.H2O.Gas (b_v=[1], n_v={-1,0}));
        output Q.NumberAbsolute RH(
          stateSelect=StateSelect.never,
          displayUnit="%") = p/(
          Modelica.Media.Air.MoistAir.saturationPressureLiquid(T/U.K)*U.Pa) if
          environment.analysis "Relative humidity (approximate)";
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="H2O",
          Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol></p>

<p>For more information, please see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"));

      end Correlated;

      model Fixed "Fixed properties"

        extends Fluid(
          redeclare replaceable package Data = FCSys.Characteristics.H2O.Gas (
                b_v=[1],n_v={-1,0}),
          final tauprime,
          redeclare parameter Q.Mobility mu=Data.mu(),
          redeclare parameter Q.TimeAbsolute nu=Data.nu(),
          redeclare parameter Q.Fluidity eta=1/(9.09e-6*U.Pa*U.s),
          redeclare parameter Q.ResistivityThermal theta=U.m*U.K/(19.6e-3*U.W),

          final k_intra,
          final alpha);

        // See the documentation for tables of values.

        output Q.NumberAbsolute RH(
          stateSelect=StateSelect.never,
          displayUnit="%") = p/(
          Modelica.Media.Air.MoistAir.saturationPressureLiquid(T/U.K)*U.Pa) if
          environment.analysis "Relative humidity (approximate)";

        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="H2O",
          Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
        <li>The generalized resistivities (&eta;, &theta;) are fixed (e.g., independent of temperature).</li>
                <li>The reaction interval (&tau;&prime;) is zero.  The rate of phase change is governed by the
        other configurations (e.g., liquid).</li>
    </ol></p>

  <p>Notes:<ul>
  <li>The specific heat capacity is not fixed because it would
  affect the chemical potential and result in an incorrect cell
  potential.</li></ul></p>

<p>The default resistivities (&eta; = <code>1/(9.09e-6*U.Pa*U.s)</code>
and &theta; = <code>U.m*U.K/(19.6e-3*U.W)</code>) are of H<sub>2</sub>O gas at saturation pressure and
  300&nbsp;K from Incropera and DeWitt [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, p. 921].  <a href=\"#Tab1\">Table 1</a> lists the properties at
  saturation pressure and other temperatures.  <a href=\"#Tab2\">Table 2</a> lists the properties of H<sub>2</sub>O gas at 1&nbsp;atm.
  See also
  <a href=\"http://www.engineeringtoolbox.com/water-dynamic-kinematic-viscosity-d_596.html\">http://www.engineeringtoolbox.com/water-dynamic-kinematic-viscosity-d_596.html</a>.</p>

  <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
    <caption align=\"top\" id=\"Tab1\">Table 1: Properties of H<sub>2</sub>O gas at saturation pressure [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 924&ndash;925].</caption>
  <tr>
      <th valign=\"middle\"><i>T</i><br><code>/U.K</code></th>
      <th width=1><i>c<sub>p</sub></i><code>*U.kg*U.K<br>/(U.J*Data.m)</code></th>
      <th width=1>&eta;<code><br>*U.Pa*U.s</code></th>
      <th width=1>&theta;<code>*U.W<br>/(U.m*U.K)</code></th>
    </tr>
<tr><td>273.15</td><td>1854</td><td>1/8.02e-6</td><td>1/18.2e-3</td></tr>
<tr><td>275</td><td>1855</td><td>1/8.09e-6</td><td>1/18.3e-3</td></tr>
<tr><td>280</td><td>1858</td><td>1/8.29e-6</td><td>1/18.6e-3</td></tr>
<tr><td>285</td><td>1861</td><td>1/8.49e-6</td><td>1/18.9e-3</td></tr>
<tr><td>290</td><td>1864</td><td>1/8.69e-6</td><td>1/19.3e-3</td></tr>
<tr><td>295</td><td>1868</td><td>1/8.89e-6</td><td>1/19.5e-3</td></tr>
<tr><td>300</td><td>1872</td><td>1/9.09e-6</td><td>1/19.6e-3</td></tr>
<tr><td>305</td><td>1877</td><td>1/9.29e-6</td><td>1/20.1e-3</td></tr>
<tr><td>310</td><td>1882</td><td>1/9.49e-6</td><td>1/20.4e-3</td></tr>
<tr><td>315</td><td>1888</td><td>1/9.69e-6</td><td>1/20.7e-3</td></tr>
<tr><td>320</td><td>1895</td><td>1/9.89e-6</td><td>1/21.0e-3</td></tr>
<tr><td>325</td><td>1903</td><td>1/10.09e-6</td><td>1/21.3e-3</td></tr>
<tr><td>330</td><td>1911</td><td>1/10.29e-6</td><td>1/21.7e-3</td></tr>
<tr><td>335</td><td>1920</td><td>1/10.49e-6</td><td>1/22.0e-3</td></tr>
<tr><td>340</td><td>1930</td><td>1/10.69e-6</td><td>1/22.3e-3</td></tr>
<tr><td>345</td><td>1941</td><td>1/10.89e-6</td><td>1/22.6e-3</td></tr>
<tr><td>350</td><td>1954</td><td>1/11.09e-6</td><td>1/23.0e-3</td></tr>
<tr><td>355</td><td>1968</td><td>1/11.29e-6</td><td>1/23.3e-3</td></tr>
<tr><td>360</td><td>1983</td><td>1/11.49e-6</td><td>1/23.7e-3</td></tr>
<tr><td>365</td><td>1999</td><td>1/11.69e-6</td><td>1/24.1e-3</td></tr>
<tr><td>370</td><td>2017</td><td>1/11.89e-6</td><td>1/24.5e-3</td></tr>
<tr><td>373.15</td><td>2029</td><td>1/12.02e-6</td><td>1/24.8e-3</td></tr>
<tr><td>375</td><td>2036</td><td>1/12.09e-6</td><td>1/24.9e-3</td></tr>
<tr><td>380</td><td>2057</td><td>1/12.29e-6</td><td>1/25.4e-3</td></tr>
<tr><td>385</td><td>2080</td><td>1/12.49e-6</td><td>1/25.8e-3</td></tr>
<tr><td>390</td><td>2104</td><td>1/12.69e-6</td><td>1/26.3e-3</td></tr>
<tr><td>400</td><td>2158</td><td>1/13.05e-6</td><td>1/27.2e-3</td></tr>
<tr><td>410</td><td>2221</td><td>1/13.42e-6</td><td>1/28.2e-3</td></tr>
<tr><td>420</td><td>2291</td><td>1/13.79e-6</td><td>1/29.8e-3</td></tr>
<tr><td>430</td><td>2369</td><td>1/14.14e-6</td><td>1/30.4e-3</td></tr>
<tr><td>440</td><td>2460</td><td>1/14.50e-6</td><td>1/31.7e-3</td></tr>
<tr><td>450</td><td>2560</td><td>1/14.85e-6</td><td>1/33.1e-3</td></tr>
<tr><td>460</td><td>2680</td><td>1/15.19e-6</td><td>1/34.6e-3</td></tr>
<tr><td>470</td><td>2790</td><td>1/15.54e-6</td><td>1/36.3e-3</td></tr>
<tr><td>480</td><td>2940</td><td>1/15.88e-6</td><td>1/38.1e-3</td></tr>
<tr><td>490</td><td>3100</td><td>1/16.23e-6</td><td>1/40.1e-3</td></tr>
<tr><td>500</td><td>3270</td><td>1/16.59e-6</td><td>1/42.3e-3</td></tr>
<tr><td>510</td><td>3470</td><td>1/16.95e-6</td><td>1/44.7e-3</td></tr>
<tr><td>520</td><td>3700</td><td>1/17.33e-6</td><td>1/47.5e-3</td></tr>
<tr><td>530</td><td>3960</td><td>1/17.72e-6</td><td>1/50.6e-3</td></tr>
<tr><td>540</td><td>4270</td><td>1/18.1e-6</td><td>1/54.0e-3</td></tr>
<tr><td>550</td><td>4640</td><td>1/18.6e-6</td><td>1/58.3e-3</td></tr>
<tr><td>560</td><td>5090</td><td>1/19.1e-6</td><td>1/63.7e-3</td></tr>
<tr><td>570</td><td>5670</td><td>1/19.7e-6</td><td>1/76.7e-3</td></tr>
<tr><td>580</td><td>6400</td><td>1/20.4e-6</td><td>1/76.7e-3</td></tr>
<tr><td>590</td><td>7350</td><td>1/21.5e-6</td><td>1/84.1e-3</td></tr>
<tr><td>600</td><td>8750</td><td>1/22.7e-6</td><td>1/92.9e-3</td></tr>
<tr><td>610</td><td>11100</td><td>1/24.1e-6</td><td>1/103e-3</td></tr>
<tr><td>620</td><td>15400</td><td>1/25.9e-6</td><td>1/114e-3</td></tr>
<tr><td>635</td><td>18300</td><td>1/27.0e-6</td><td>1/121e-3</td></tr>
<tr><td>630</td><td>22100</td><td>1/28.0e-6</td><td>1/130e-3</td></tr>
<tr><td>635</td><td>27600</td><td>1/30.0e-6</td><td>1/141e-3</td></tr>
<tr><td>640</td><td>42000</td><td>1/32.0e-6</td><td>1/155e-3</td></tr>
  </table>

<br>

    <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
    <caption align=\"top\" id=\"Tab2\">Table 2: Properties of H<sub>2</sub>O gas at 1&nbsp;atm [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, p. 921].</caption>
  <tr>
      <th valign=\"middle\"><i>T</i><br><code>/U.K</code></th>
      <th width=1><i>c<sub>p</sub></i><code>*U.kg*U.K<br>/(U.J*Data.m)</code></th>
      <th width=1>&eta;<code><br>*U.Pa*U.s</code></th>
      <th width=1>&theta;<code>*U.W<br>/(U.m*U.K)</code></th>
    </tr>
<tr><td>380</td><td>2.060e3</td><td>1/127.1e-7</td><td>1/24.6e-3</td></tr>
<tr><td>400</td><td>2.014e3</td><td>1/134.4e-7</td><td>1/26.1e-3</td></tr>
<tr><td>450</td><td>1.980e3</td><td>1/152.5e-7</td><td>1/29.9e-3</td></tr>
<tr><td>500</td><td>1.985e3</td><td>1/170.4e-7</td><td>1/33.9e-3</td></tr>
<tr><td>550</td><td>1.997e3</td><td>1/188.4e-7</td><td>1/37.9e-3</td></tr>
<tr><td>600</td><td>2.206e3</td><td>1/206.7e-7</td><td>1/42.2e-3</td></tr>
<tr><td>650</td><td>2.056e3</td><td>1/224.7e-7</td><td>1/46.4e-3</td></tr>
<tr><td>700</td><td>2.085e3</td><td>1/242.6e-7</td><td>1/50.5e-3</td></tr>
<tr><td>750</td><td>2.119e3</td><td>1/260.4e-7</td><td>1/54.9e-3</td></tr>
<tr><td>800</td><td>2.152e3</td><td>1/278.6e-7</td><td>1/59.2e-3</td></tr>
<tr><td>850</td><td>2.186e3</td><td>1/296.9e-7</td><td>1/63.7e-3</td></tr>
  </table></ul></p>

<p>For more information, please see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"));

      end Fixed;

    end Gas;

    package Ionomer "<html>H<sub>2</sub>O in ionomer</html>"
      extends Modelica.Icons.Package;

      model Correlated "Correlated properties"

        extends Fluid(redeclare replaceable package Data =
              Characteristics.H2O.Ionomer, final tauprime);

        // Auxiliary variables (for analysis)
        output Q.NumberAbsolute lambda(stateSelect=StateSelect.never) = rho*
          Characteristics.'SO3-'.Ionomer.b_v[1, 1] if environment.analysis
          "Ratio of H2O molecules to SO3- end-groups";
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="H2O",
          Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
    <li>The reaction interval (&tau;&prime;) is zero.  The rate of phase change is governed by the
   other configurations (e.g., gas).</li>
          </ol></p>

<p>For more information, please see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"));

      end Correlated;

      model Fixed "Fixed properties"

        extends Fluid(
          redeclare replaceable package Data = Characteristics.H2O.Ionomer,
          redeclare parameter Q.TimeAbsolute tauprime=2e10*Data.tauprime(),
          redeclare parameter Q.Mobility mu=Data.mu(),
          redeclare parameter Q.TimeAbsolute nu=Data.nu(),
          redeclare parameter Q.Fluidity eta=Data.eta(),
          redeclare parameter Q.ResistivityThermal theta=Data.theta(),
          final N_IC,
          final h_IC,
          final g_IC,
          final rho_IC,
          final p_IC,
          final V_IC,
          final T_IC,
          final alpha,
          final consMaterial,
          final initMaterial=Init.none,
          final initEnergy=Init.none);

        parameter Q.NumberAbsolute lambda_IC=14
          "<html>Initial ratio of H<sub>2</sub>O molecules to SO<sub>3</sub><sup>-</sup> end-groups</html>"
          annotation (Dialog(tab="Initialization", __Dymola_label=
                "<html>&lambda;<sub>IC</sub></html>"));

        // Auxiliary variables (for analysis)
        Q.NumberAbsolute lambda(
          start=lambda_IC,
          fixed=true,
          stateSelect=StateSelect.never)
          "Ratio of H2O molecules to SO3- end-groups";

      equation
        lambda = rho*Characteristics.'SO3-'.Ionomer.b_v[1, 1];

        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="H2O",
          Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
        <li>The generalized resistivities (&eta;, &theta;) are fixed (e.g., independent of temperature).</li>
    </ol></p></p>

<p>For more information, please see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"),

          Diagram(graphics),
          Icon(graphics));

      end Fixed;

    end Ionomer;

    package Liquid "<html>H<sub>2</sub>O liquid</html>"
      extends Modelica.Icons.Package;

      model Correlated "Correlated properties"

        extends Fluid(
          redeclare replaceable package Data = Characteristics.H2O.Liquid,
          final tauprime,
          initMaterial=Init.volume);
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="H2O",
          Documentation(info="<html><p>Assumptions:<ol>
        <li>The reaction interval (&tau;&prime;) is zero.  The rate of phase change is governed by the
        other configurations (e.g., gas).</li>
    </ol></p>

<p>For more information, please see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"));

      end Correlated;

      model Fixed "Fixed properties"

        // Initialization
        parameter Q.Number epsilon_IC=0.01 "Initial volumetric fill fraction"
          annotation (Dialog(tab="Initialization", __Dymola_label=
                "<html>&epsilon;<sub>IC</sub></html>"));

        extends Fluid(
          redeclare replaceable package Data = Characteristics.H2O.Liquid,
          redeclare parameter Q.TimeAbsolute tauprime=5e10*Data.tauprime(),
          redeclare parameter Q.Mobility mu=Data.mu(),
          redeclare parameter Q.TimeAbsolute nu=Data.nu(),
          redeclare parameter Q.Fluidity eta=1/(855e-6*U.Pa*U.s),
          redeclare parameter Q.ResistivityThermal theta=U.m*U.K/(0.613*U.W),
          final k_intra,
          final initEnergy=Init.temperature,
          final N_IC,
          final h_IC,
          final g_IC,
          final rho_IC,
          final p_IC,
          final V_IC=epsilon_IC*product(L),
          final initMaterial=Init.volume,
          final alpha);

        // See the documentation for tables of values.

        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="H2O",
          Documentation(info="<html><p>Assumptions:<ol>
        <li>The generalized resistivities (&eta;, &theta;) are fixed (e.g., independent of temperature).</li>
    </ol></p>

          <p>Notes:<ul>
  <li>
  The specific heat capacity is not fixed because it would
  affect the chemical potential and result in an incorrect cell
  potential.</li></ul></p>

<p>The default resistivities (&eta; = <code>1/(855e-6*U.Pa*U.s)</code>
and &theta; = <code>U.m*U.K/(613e-3*U.W)</code>) are of H<sub>2</sub>O liquid at saturation pressure and
  300&nbsp;K from Incropera and DeWitt [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, p. 921].  <a href=\"#Tab1\">Table 1</a> lists the properties at
  saturation pressure and other temperatures.
  See also
  <a href=\"http://www.engineeringtoolbox.com/water-dynamic-kinematic-viscosity-d_596.html\">http://www.engineeringtoolbox.com/water-dynamic-kinematic-viscosity-d_596.html</a>.</p>

  <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
    <caption align=\"top\" id=\"Tab1\">Table 1: Properties of H<sub>2</sub>O liquid at saturation pressure [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 924&ndash;925].</caption>
  <tr>
      <th valign=\"middle\"><i>T</i><br><code>/U.K</code></th>
      <th width=1><i>c<sub>p</sub></i><code>*U.kg*U.K<br>/(U.J*Data.m)</code></th>
      <th width=1>&eta;<code><br>*U.Pa*U.s</code></th>
      <th width=1>&theta;<code>*U.W<br>/(U.m*U.K)</code></th>
    </tr>
<tr><td>273.15</td><td>4217</td><td>1/1750e-6</td><td>1/569e-3</td></tr>
<tr><td>275</td><td>4211</td><td>1/1652e-6</td><td>1/574e-3</td></tr>
<tr><td>280</td><td>4198</td><td>1/1422e-6</td><td>1/582e-3</td></tr>
<tr><td>285</td><td>4189</td><td>1/1225e-6</td><td>1/590e-3</td></tr>
<tr><td>290</td><td>4184</td><td>1/1080e-6</td><td>1/598e-3</td></tr>
<tr><td>295</td><td>4181</td><td>1/959e-6</td><td>1/606e-3</td></tr>
<tr><td>300</td><td>4179</td><td>1/855e-6</td><td>1/613e-3</td></tr>
<tr><td>305</td><td>4178</td><td>1/769e-6</td><td>1/620e-3</td></tr>
<tr><td>310</td><td>4178</td><td>1/695e-6</td><td>1/628e-3</td></tr>
<tr><td>315</td><td>4179</td><td>1/631e-6</td><td>1/634e-3</td></tr>
<tr><td>320</td><td>4180</td><td>1/577e-6</td><td>1/640e-3</td></tr>
<tr><td>325</td><td>4182</td><td>1/528e-6</td><td>1/645e-3</td></tr>
<tr><td>330</td><td>4184</td><td>1/489e-6</td><td>1/650e-3</td></tr>
<tr><td>335</td><td>4186</td><td>1/453e-6</td><td>1/656e-3</td></tr>
<tr><td>340</td><td>4188</td><td>1/420e-6</td><td>1/660e-3</td></tr>
<tr><td>345</td><td>4191</td><td>1/389e-6</td><td>1/668e-3</td></tr>
<tr><td>350</td><td>4195</td><td>1/365e-6</td><td>1/668e-3</td></tr>
<tr><td>355</td><td>4199</td><td>1/343e-6</td><td>1/671e-3</td></tr>
<tr><td>360</td><td>4203</td><td>1/324e-6</td><td>1/674e-3</td></tr>
<tr><td>365</td><td>4209</td><td>1/306e-6</td><td>1/677e-3</td></tr>
<tr><td>370</td><td>4214</td><td>1/289e-6</td><td>1/679e-3</td></tr>
<tr><td>373.15</td><td>4217</td><td>1/279e-6</td><td>1/680e-3</td></tr>
<tr><td>375</td><td>4220</td><td>1/274e-6</td><td>1/681e-3</td></tr>
<tr><td>380</td><td>4226</td><td>1/260e-6</td><td>1/683e-3</td></tr>
<tr><td>385</td><td>4232</td><td>1/248e-6</td><td>1/685e-3</td></tr>
<tr><td>390</td><td>4239</td><td>1/237e-6</td><td>1/686e-3</td></tr>
<tr><td>400</td><td>4256</td><td>1/217e-6</td><td>1/688e-3</td></tr>
<tr><td>410</td><td>4278</td><td>1/200e-6</td><td>1/688e-3</td></tr>
<tr><td>420</td><td>4302</td><td>1/185e-6</td><td>1/688e-3</td></tr>
<tr><td>430</td><td>4331</td><td>1/173e-6</td><td>1/685e-3</td></tr>
<tr><td>440</td><td>4360</td><td>1/162e-6</td><td>1/682e-3</td></tr>
<tr><td>450</td><td>4400</td><td>1/152e-6</td><td>1/678e-3</td></tr>
<tr><td>460</td><td>4440</td><td>1/143e-6</td><td>1/673e-3</td></tr>
<tr><td>470</td><td>4480</td><td>1/136e-6</td><td>1/667e-3</td></tr>
<tr><td>480</td><td>4530</td><td>1/129e-6</td><td>1/660e-3</td></tr>
<tr><td>490</td><td>4590</td><td>1/124e-6</td><td>1/651e-3</td></tr>
<tr><td>500</td><td>4660</td><td>1/118e-6</td><td>1/642e-3</td></tr>
<tr><td>510</td><td>4740</td><td>1/113e-6</td><td>1/631e-3</td></tr>
<tr><td>520</td><td>4840</td><td>1/108e-6</td><td>1/621e-3</td></tr>
<tr><td>530</td><td>4950</td><td>1/104e-6</td><td>1/608e-3</td></tr>
<tr><td>540</td><td>5080</td><td>1/101e-6</td><td>1/594e-3</td></tr>
<tr><td>550</td><td>5240</td><td>1/97e-6</td><td>1/580e-3</td></tr>
<tr><td>560</td><td>5430</td><td>1/94e-6</td><td>1/563e-3</td></tr>
<tr><td>570</td><td>5680</td><td>1/91e-6</td><td>1/548e-3</td></tr>
<tr><td>580</td><td>6000</td><td>1/88e-6</td><td>1/528e-3</td></tr>
<tr><td>590</td><td>6410</td><td>1/84e-6</td><td>1/513e-3</td></tr>
<tr><td>600</td><td>7000</td><td>1/81e-6</td><td>1/497e-3</td></tr>
<tr><td>610</td><td>7850</td><td>1/77e-6</td><td>1/467e-3</td></tr>
<tr><td>620</td><td>9350</td><td>1/72e-6</td><td>1/444e-3</td></tr>
<tr><td>635</td><td>10600</td><td>1/70e-6</td><td>1/430e-3</td></tr>
<tr><td>630</td><td>12600</td><td>1/67e-6</td><td>1/412e-3</td></tr>
<tr><td>635</td><td>16400</td><td>1/64e-6</td><td>1/392e-3</td></tr>
<tr><td>640</td><td>26000</td><td>1/59e-6</td><td>1/367e-3</td></tr>
  </table>

  <p>For more information, please see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"),

          Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                  {100,100}}), graphics));

      end Fixed;

    end Liquid;

  end H2O;

  package N2 "<html>N<sub>2</sub></html>"
    extends Modelica.Icons.Package;
    package Gas "<html>N<sub>2</sub> gas</html>"
      extends Modelica.Icons.Package;

      model Correlated "Correlated properties"

        extends Fluid(redeclare replaceable package Data =
              Characteristics.N2.Gas (b_v=[1], n_v={-1,0}));
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="N2",
          Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol></p>

<p>For more information, please see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"));

      end Correlated;

      model Fixed "Fixed properties"
        import FCSys.Utilities.Polynomial;

        extends Fluid(
          redeclare replaceable package Data = FCSys.Characteristics.N2.Gas (
              b_v=[1],
              n_v={-1,0},
              n_c=0,
              T_lim_c={0,Modelica.Constants.inf},
              b_c=[1041*U.J*Data.m/(U.kg*U.K)],
              B_c=[Data.Deltah0_f - (1041*U.J*Data.m/U.kg)*298.15, 191.610*U.J/
                  (U.mol*U.K) - (1041*U.J*Data.m/(U.kg*U.K))*log(298.15*U.K)]),

          final tauprime,
          redeclare parameter Q.Mobility mu=Data.mu(),
          redeclare parameter Q.TimeAbsolute nu=Data.nu(),
          redeclare parameter Q.Fluidity eta=1/(17.82e-6*U.Pa*U.s),
          redeclare parameter Q.ResistivityThermal theta=U.m*U.K/(25.9e-3*U.W),

          final k_intra,
          final alpha);

        // See the documentation for a table of values.
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="N2",
          Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
    <li>Fixed specific heat capacity (independent of temperature)</li>
        <li>The generalized resistivities (&eta;, &theta;) are fixed (e.g., independent of temperature).</li>

    </ol></p>

<p>The default specific heat capacity (<i>b<sub>c</sub></i> = <code>[1.041e3*U.J*Data.m/(U.kg*U.K)]</code>) and resistivities
(&eta; = <code>1/(17.82e-6*U.Pa*U.s)</code> and &theta; = <code>U.m*U.K/(25.9e-3*U.W))</code>) are based on data of gas at 1&nbsp;atm and
  300&nbsp;K from Incropera and DeWitt [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, p. 920].
   The integration offset for specific entropy is set such that
   the specific entropy is 191.610&nbsp;J/(mol&middot;K) at 25&nbsp;&deg;C and <i>p</i><sup>o</sup> (1&nbsp;bar).
   This is the value from Table B in [<a href=\"modelica://FCSys.UsersGuide.References\">McBride2002</a>].
   Additional thermal data is listed in <a href=\"#Tab1\">Table 1</a>.

  <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
  <caption align=\"top\" id=\"Tab1\">Table 1: Properties of N<sub>2</sub> gas at 1&nbsp;atm [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, p. 920]</caption>
  <tr>
      <th valign=\"middle\"><i>T</i><br><code>/U.K</code></th>
      <th width=1><i>c<sub>p</sub></i><code>*U.kg*U.K<br>/(U.J*Data.m)</code></th>
      <th width=1>&eta;<code><br>*U.Pa*U.s</code></th>
      <th width=1>&theta;<code>*U.W<br>/(U.m*U.K)</code></th>
    </tr>
<tr><td>100</td><td>1.070e3</td><td>1/68.8e-7</td><td>1/9.58e-3</td></tr>
<tr><td>150</td><td>1.050e3</td><td>1/100.6e-7</td><td>1/13.9e-3</td></tr>
<tr><td>200</td><td>1.043e3</td><td>1/129.2e-7</td><td>1/18.3e-3</td></tr>
<tr><td>250</td><td>1.042e3</td><td>1/154.9e-7</td><td>1/22.2e-3</td></tr>
<tr><td>300</td><td>1.041e3</td><td>1/178.2e-7</td><td>1/25.9e-3</td></tr>
<tr><td>350</td><td>1.042e3</td><td>1/200.0e-7</td><td>1/29.3e-3</td></tr>
<tr><td>400</td><td>1.045e3</td><td>1/220.4e-7</td><td>1/32.7e-3</td></tr>
<tr><td>450</td><td>1.050e3</td><td>1/239.6e-7</td><td>1/35.8e-3</td></tr>
<tr><td>500</td><td>1.056e3</td><td>1/257.7e-7</td><td>1/38.9e-3</td></tr>
<tr><td>550</td><td>1.065e3</td><td>1/274.7e-7</td><td>1/41.7e-3</td></tr>
<tr><td>600</td><td>1.075e3</td><td>1/290.8e-7</td><td>1/44.6e-3</td></tr>
<tr><td>700</td><td>1.098e3</td><td>1/320.1e-7</td><td>1/49.9e-3</td></tr>
<tr><td>800</td><td>1.220e3 [sic]</td><td>1/349.1e-7</td><td>1/54.8e-3</td></tr>
<tr><td>900</td><td>1.146e3</td><td>1/375.3e-7</td><td>1/59.7e-3</td></tr>
<tr><td>1000</td><td>1.167e3</td><td>1/399.9e-7</td><td>1/64.7e-3</td></tr>
<tr><td>1100</td><td>1.187e3</td><td>1/423.2e-7</td><td>1/70.0e-3</td></tr>
<tr><td>1200</td><td>1.204e3</td><td>1/445.3e-7</td><td>1/75.8e-3</td></tr>
<tr><td>1300</td><td>1.219e3</td><td>1/466.2e-7</td><td>1/81.0e-3</td></tr>
  </table>

  <p>The fluidity of air at 15.0&nbsp;&deg;C and 1&nbsp;atm is given by
       &eta; = <code>1/(17.8e-6*U.Pa*U.s)</code>
   (<a href=\"http://en.wikipedia.org/wiki/Viscosity\">http://en.wikipedia.org/wiki/Viscosity</a>).</p>

<p>For more information, please see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"));

      end Fixed;

    end Gas;

  end N2;

  package O2 "<html>O<sub>2</sub></html>"
    extends Modelica.Icons.Package;
    package Gas "<html>O<sub>2</sub> gas</html>"
      extends Modelica.Icons.Package;

      model Correlated "Correlated properties"

        extends Fluid(redeclare replaceable package Data =
              Characteristics.O2.Gas (b_v=[1], n_v={-1,0}));
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="O2",
          Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol></p>

<p>For more information, please see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"));

      end Correlated;

      model Fixed "Fixed properties"

        extends Fluid(
          redeclare replaceable package Data = FCSys.Characteristics.O2.Gas (
                b_v=[1], n_v={-1,0}),
          final tauprime,
          redeclare parameter Q.Mobility mu=Data.mu(),
          redeclare parameter Q.TimeAbsolute nu=Data.nu(),
          redeclare parameter Q.Fluidity eta=1/(20.72e-6*U.Pa*U.s),
          redeclare parameter Q.ResistivityThermal theta=U.m*U.K/(26.8e-3*U.W),

          final k_intra,
          final alpha);

        // See the documentation for a table of values.

        parameter Q.PressureAbsolute p_stop=5*U.Pa
          "Pressure below which the simulation should terminate" annotation (
            Dialog(tab="Advanced", __Dymola_label=
                "<html><i>p</i><sub>stop</sub></html>"));

      equation
        when p < p_stop then
          terminate("There is no more " + Data.formula + ".");
        end when;

        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="O2",
          Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
        <li>The generalized resistivities (&eta;, &theta;) are fixed (e.g., independent of temperature).</li>

    </ol></p>

<p>Additional notes:
<ul>
          <li>
  The specific heat capacity is not fixed because it would
  affect the chemical potential and result in an incorrect cell
  potential.</li></ul></p>

  <p>The default resistivities (&eta; = <code>1/(207.2e-7*U.Pa*U.s)</code> and &theta; = <code>U.m*U.K/(26.8e-3*U.W)</code>) are based on data of gas at 1&nbsp;atm and
  300&nbsp;K from Incropera and DeWitt [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 920&ndash;921].
  <a href=\"#Tab1\">Table 1</a> lists the properties at other temperatures.</p>

  <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
  <caption align=\"top\" id=\"Tab1\">Table 1: Properties of O<sub>2</sub> gas at 1&nbsp;atm
  [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 920&ndash;921]</caption>
  <tr>
      <th valign=\"middle\"><i>T</i><br><code>/U.K</code></th>
      <th width=1><i>c<sub>p</sub></i><code>*U.kg*U.K<br>/(U.J*Data.m)</code></th>
      <th width=1>&eta;<code><br>*U.Pa*U.s</code></th>
      <th width=1>&theta;<code>*U.W<br>/(U.m*U.K)</code></th>
    </tr>
<tr><td>100</td><td>0.962e3</td><td>1/76.4e-7</td><td>1/9.25e-3</td></tr>
<tr><td>150</td><td>0.921e3</td><td>1/114.8e-7</td><td>1/13.8e-3</td></tr>
<tr><td>200</td><td>0.915e3</td><td>1/147.5e-7</td><td>1/18.3e-3</td></tr>
<tr><td>250</td><td>0.915e3</td><td>1/178.6e-7</td><td>1/22.6e-3</td></tr>
<tr><td>300</td><td>0.920e3</td><td>1/207.2e-7</td><td>1/26.8e-3</td></tr>
<tr><td>350</td><td>0.929e3</td><td>1/233.5e-7</td><td>1/29.6e-3</td></tr>
<tr><td>400</td><td>0.942e3</td><td>1/258.2e-7</td><td>1/33.0e-3</td></tr>
<tr><td>450</td><td>0.956e3</td><td>1/281.4e-7</td><td>1/36.3e-3</td></tr>
<tr><td>500</td><td>0.972e3</td><td>1/303.3e-7</td><td>1/41.2e-3</td></tr>
<tr><td>550</td><td>0.988e3</td><td>1/324.0e-7</td><td>1/44.1e-3</td></tr>
<tr><td>600</td><td>1.003e3</td><td>1/343.7e-7</td><td>1/47.3e-3</td></tr>
<tr><td>700</td><td>1.031e3</td><td>1/380.8e-7</td><td>1/52.8e-3</td></tr>
<tr><td>800</td><td>1.054e3</td><td>1/415.2e-7</td><td>1/58.9e-3</td></tr>
<tr><td>900</td><td>1.074e3</td><td>1/447.2e-7</td><td>1/64.9e-3</td></tr>
<tr><td>1000</td><td>1.090e3</td><td>1/477.0e-7</td><td>1/71.0e-3</td></tr>
<tr><td>1100</td><td>1.103e3</td><td>1/505.5e-7</td><td>1/75.8e-3</td></tr>
<tr><td>1200</td><td>1.115e3</td><td>1/532.5e-7</td><td>1/81.9e-3</td></tr>
<tr><td>1300</td><td>1.125e3</td><td>1/588.4e-7</td><td>1/87.1e-3</td></tr>
  </table></p>

<p>For more information, please see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"),

          Icon(graphics));

      end Fixed;

    end Gas;

  end O2;

  model Ion "Base model for an ion"
    import assert = FCSys.Utilities.assertEval;
    extends Fluid(
      alpha=0.5,
      final mu=sigma*v,
      dalton(final V=(N - electrostatic.Z/Data.z)*v));
    // Note:  The amount of material includes the additional amount in the
    // charge layer, which may be positive or negative.

    parameter Q.ConductivityElectrical sigma=Data.mu()/Data.v_Tp()
      "Electrical conductivity" annotation (Dialog(group="Material properties",
          __Dymola_label="<html>&sigma;</html>"));

    Connectors.Electrostatic electrostatic "Interface with the dielectric"
      annotation (Placement(transformation(extent={{-30,10},{-10,30}}),
          iconTransformation(extent={{-80,60},{-60,80}})));

  initial equation
    assert(Data.z <> 0, "The Ion model can only be used for charged species.");

  equation
    electrochemical.w = g + Data.z*electrostatic.w
      "The electrochemical potential is the sum of chemical and electrostatic contributions.";

    annotation (
      defaultComponentPrefixes="replaceable",
      Documentation(info="<html>
                <p>Please see the
     <a href=\"modelica://FCSys.Species.Fluid\">Fluid</a> model.</p></html>"),
      Diagram(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          initialScale=0.1),graphics),
      Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{100,
              100}}), graphics));
  end Ion;

  model Fluid "Base model for a fluid species"

    import Modelica.Math.asinh;
    import FCSys.Utilities.Coordinates.after;
    import FCSys.Utilities.Coordinates.before;
    import FCSys.Utilities.Coordinates.cartWrap;
    import FCSys.Utilities.inSign;
    import FCSys.Utilities.Delta;
    import FCSys.Utilities.selectBooleans;
    import FCSys.Utilities.selectIntegers;
    import assert = FCSys.Utilities.assertEval;

    // Geometric parameters
    parameter Q.NumberAbsolute alpha=1 "Charge transfer coefficient"
      annotation (Dialog(group="Geometry", __Dymola_label=
            "<html>&alpha;</html>"));

    // Initialization parameters
    parameter Init initMaterial=Init.pressure
      "Method of initializing the material state"
      annotation (Evaluate=true, Dialog(tab="Initialization"));
    parameter Init initEnergy=Init.temperature
      "Method of initializing the thermal state"
      annotation (Evaluate=true, Dialog(tab="Initialization"));
    extends Species(N(stateSelect=if consMaterial == ConsThermo.dynamic then
            StateSelect.always else StateSelect.prefer));
    // Note:  The extension is after these parameters so that they appear first
    // in the parameter dialog.
    // Note:  StateSelect.always is not ideal, but it is necessary to avoid dynamic
    // state selection in Dymola 2014.  In some cases it isn't appropriate (e.g.,
    // an incompressible liquid that fills the entire subregion), and it those
    // cases it can be modified at instantiation.

    // Material properties
    Q.TimeAbsolute tauprime(nominal=1e-6*U.s) = 0
      "Interval for chemical exchange" annotation (Dialog(group=
            "Material properties", __Dymola_label="<html>&tau;&prime;</html>"));
    Q.Resistivity eta(nominal=10*U.cm/U.A) = Data.eta(T, v) "Fluidity"
      annotation (Dialog(group="Material properties", __Dymola_label=
            "<html>&eta;</html>"));
    Q.ResistivityThermal theta(nominal=10*U.cm/U.A) = Data.theta(T, v)
      "Thermal resistivity" annotation (Dialog(group="Material properties",
          __Dymola_label="<html>&theta;</html>"));

    // Assumptions
    // -----------
    // Upstream discretization
    parameter Boolean upstreamX=true "X" annotation (
      Evaluate=true,
      Dialog(
        tab="Assumptions",
        group="Axes with upstream discretization",
        enable=inclTrans[1],
        compact=true),
      choices(__Dymola_checkBox=true));
    parameter Boolean upstreamY=true "Y" annotation (
      Evaluate=true,
      Dialog(
        tab="Assumptions",
        group="Axes with upstream discretization",
        enable=inclTrans[2],
        compact=true),
      choices(__Dymola_checkBox=true));
    parameter Boolean upstreamZ=true "Z" annotation (
      Evaluate=true,
      Dialog(
        tab="Assumptions",
        group="Axes with upstream discretization",
        enable=inclTrans[3],
        compact=true),
      choices(__Dymola_checkBox=true));
    //
    // Dynamics
    parameter ConsThermo consMaterial=ConsThermo.dynamic "Material" annotation
      (Evaluate=true, Dialog(tab="Assumptions", group=
            "Formulation of the conservation equations"));
    parameter Boolean consRot=false "Conserve rotational momentum" annotation (
      Evaluate=true,
      Dialog(tab="Assumptions", group=
            "Formulation of the conservation equations"),
      choices(__Dymola_checkBox=true));

    parameter ConsMom consTransX=ConsMom.dynamic
      "X-axis translational momentum" annotation (Evaluate=true, Dialog(
        tab="Assumptions",
        group="Formulation of the conservation equations",
        enable=inclTrans[1]));
    parameter ConsMom consTransY=ConsMom.dynamic
      "Y-axis translational momentum" annotation (Evaluate=true, Dialog(
        tab="Assumptions",
        group="Formulation of the conservation equations",
        enable=inclTrans[2]));
    parameter ConsMom consTransZ=ConsMom.dynamic
      "Z-axis translational momentum" annotation (Evaluate=true, Dialog(
        tab="Assumptions",
        group="Formulation of the conservation equations",
        enable=inclTrans[3]));
    parameter ConsThermo consEnergy=ConsThermo.dynamic "Energy" annotation (
        Evaluate=true, Dialog(tab="Assumptions", group=
            "Formulation of the conservation equations"));
    //
    // Flow conditions
    parameter Q.NumberAbsolute Nu_Phi[Axis]={4,4,4}
      "Translational Nusselt numbers" annotation (Dialog(
        tab="Assumptions",
        group="Flow conditions",
        __Dymola_label="<html><b><i>Nu</i><sub>&Phi;</sub></b></html>"));
    parameter Q.NumberAbsolute Nu_Q=1 "Thermal Nusselt number" annotation (
        Dialog(
        tab="Assumptions",
        group="Flow conditions",
        __Dymola_label="<html><i>Nu</i><sub><i>Q</i></sub></html>"));

    // Aliases (for common terms)
    Q.Current I[n_trans](each nominal=U.A, each stateSelect=StateSelect.never)
      "Current";
    Q.Velocity phi_faces[n_trans, Side](each nominal=10*U.cm/U.s, each
        stateSelect=StateSelect.never) "Normal velocities at the faces";

    // Auxiliary variables (for analysis)
    // ----------------------------------
    // Misc. conditions
    output Q.Density rho_faces[n_trans, Side](each stateSelect=StateSelect.never)
       = fill(
        1,
        n_trans,
        2) ./ Data.v_Tp(faces.T, faces.p) if environment.analysis
      "Densities at the faces";
    output Q.VolumeRate Vdot_faces[n_trans, Side](each stateSelect=StateSelect.never)
       = faces.Ndot ./ rho_faces if environment.analysis
      "Volume flow rates into the faces";
    output Q.PressureAbsolute q[n_trans](each stateSelect=StateSelect.never) =
      Data.m*phi .* I ./ (2*A[cartTrans]) if environment.analysis
      "Dynamic pre ssure";
    output Q.Velocity phi_chemical[n_trans](each stateSelect=StateSelect.never)
       = actualStream(electrochemical.phi) if environment.analysis
      "Velocity of the chemical stream";
    output Q.PotentialAbsolute sT_chemical(stateSelect=StateSelect.never) =
      actualStream(electrochemical.sT) if environment.analysis
      "Specific entropy-temperature product of the chemical stream";
    //
    // Potentials and current
    output Q.Potential g_faces[n_trans, Side](each stateSelect=StateSelect.never)
       = Data.g(faces.T, faces.p) if environment.analysis and not Data.isCompressible
      "Gibbs potentials at the faces";
    output Q.Potential Deltag[n_trans](each stateSelect=StateSelect.never) =
      Delta(g_faces) if environment.analysis and not Data.isCompressible
      "Differences in Gibbs potentials across the faces";
    // Note:  If a face is left unconnnected, then it is possible that its pressure
    // may become negative.  If the equation of state has ideal-gas term,
    // then the Gibbs energy will involve a logarithm of pressure.  Therefore,
    // to be safe, these variables are included only for incompressible species.
    //
    // Time constants
    output Q.TimeAbsolute tau_NT[n_trans](
      each stateSelect=StateSelect.never,
      each start=U.s) = fill(beta*beta, n_trans) ./ k[cartTrans] if environment.analysis
      "Time constants for material transport";
    output Q.TimeAbsolute tau_PhiT_para[n_trans](
      each stateSelect=StateSelect.never,
      each start=U.s) = (M*eta) ./ (Lprime .* Nu_Phi[cartTrans]) if environment.analysis
      "Time constants for transverse translational transport (through the whole subregion)";
    output Q.TimeAbsolute tau_QT[n_trans](
      each stateSelect=StateSelect.never,
      each start=U.s) = fill(N*c_v*theta/Nu_Q, n_trans) ./ Lprime if
      environment.analysis
      "Time constants for thermal transport (through the whole subregion)";
    //
    // Peclet numbers (only for the transport axes; others are zero)
    output Q.Number Pe_Phi_para[n_trans](each stateSelect=StateSelect.never) =
      (eta*Data.m/2)*I ./ Lprime if environment.analysis
      "Translational Peclet numbers";
    output Q.Number Pe_Q[n_trans](each stateSelect=StateSelect.never) = (theta*
      Data.c_v(T, p)/2)*I ./ Lprime if environment.analysis
      "Thermal Peclet numbers";
    //
    // Bulk flow rates
    output Q.Force mphiI[n_trans, n_trans](each stateSelect=StateSelect.never)
       = outerProduct(I, Data.m*phi) if environment.analysis
      "Bulk rate of translational advection (1st index: transport axis, 2nd index: translational component)";
    output Q.VolumeRate Vdot[n_trans](each stateSelect=StateSelect.never) = v*I
      if environment.analysis "Bulk volumetric flow rate";
    output Q.Power hI[n_trans](each stateSelect=StateSelect.never) = h*I if
      environment.analysis "Bulk enthalpy flow rate";
    //
    // Translational momentum balance
    output Q.Force Ma[n_trans](each stateSelect=StateSelect.never) = M*(der(phi)
      /U.s + environment.a[cartTrans]) + N*Data.z*environment.E[cartTrans] if
      environment.analysis
      "Acceleration force (including acceleration due to body forces)";
    output Q.Force f_thermo[n_trans](each stateSelect=StateSelect.never) = {(
      if inclTrans[cartTrans[i]] then -Delta(faces[transCart[cartTrans[i]], :].p)
      *A[cartTrans[i]] else 0) for i in 1:n_trans} if environment.analysis
      "Thermodynamic force";
    output Q.Force f_AE[n_trans](each stateSelect=StateSelect.never) = Data.m*(
      actualStream(electrochemical.phi) - phi) .* electrochemical.Ndot if
      environment.analysis "Acceleration force due to advective exchange";
    output Q.Force f_DE[n_trans](each stateSelect=StateSelect.never) = direct.translational.mPhidot
       + {sum(inter[:].mPhidot[i]) for i in 1:n_trans} + {sum(intra[:].mPhidot[
      i]) for i in 1:n_trans} if environment.analysis
      "Friction from other configurations (diffusive exchange)";
    // Note:  The [:] is necessary in Dymola 2014.
    output Q.Force f_AT[n_trans](each stateSelect=StateSelect.never) = {sum(((
      if i == j then phi_faces[j, :] else faces[j, :].phi[cartWrap(cartTrans[i]
       - cartTrans[j])]) - {phi[i],phi[i]})*faces[j, :].Ndot*Data.m for j in 1:
      n_trans) for i in 1:n_trans} if environment.analysis
      "Acceleration force due to advective transport";
    output Q.Force f_DT[n_trans](each stateSelect=StateSelect.never) = {sum(sum(
      if i == j then {0,0} else faces[j, :].mPhidot[cartWrap(cartTrans[i] -
      cartTrans[j])]) for j in 1:n_trans) for i in 1:n_trans} if environment.analysis
      "Shear force from other subregions (diffusive transport)";
    //
    // Energy balance
    output Q.Power Ndere(stateSelect=StateSelect.never) = (N*T*der(s) + M*phi*
      der(phi))/U.s if environment.analysis
      "Rate of energy storage (internal and kinetic) and boundary work at constant mass";
    // Note that T*der(s) = der(u) + p*der(v).
    output Q.Power Edot_AE(stateSelect=StateSelect.never) = (electrochemical.w
       + actualStream(electrochemical.sT) - h + (actualStream(electrochemical.phi)
      *actualStream(electrochemical.phi) - phi*phi)*Data.m/2)*electrochemical.Ndot
      if environment.analysis
      "Relative rate of energy (internal, flow, and kinetic) due to phase change and reaction";
    output Q.Power Edot_DE(stateSelect=StateSelect.never) = direct.translational.phi
      *direct.translational.mPhidot + sum(inter[i].phi*inter[i].mPhidot for i
       in 1:n_inter) + sum(intra[i].phi*intra[i].mPhidot for i in 1:n_intra) +
      direct.thermal.Qdot + sum(intra.Qdot) + sum(inter.Qdot) if environment.analysis
      "Rate of diffusion of energy from other configurations";
    output Q.Power Edot_AT(stateSelect=StateSelect.never) = sum((Data.h(faces[i,
      :].T, faces[i, :].p) - {h,h} + (phi_faces[i, :] .^ 2 + sum(faces[i, :].phi[
      orient] .^ 2 for orient in Orient) - fill(phi*phi, 2))*(Data.m/2))*faces[
      i, :].Ndot for i in 1:n_trans) if environment.analysis
      "Relative rate of energy (internal, flow, and kinetic) due to advective transport";
    output Q.Power Edot_DT(stateSelect=StateSelect.never) = sum(sum(faces[i, :].phi[
      orient]*faces[i, :].mPhidot[orient] for orient in Orient) for i in 1:
      n_trans) + sum(faces.Qdot) if environment.analysis
      "Rate of diffusion of energy from other subregions";
    // Note:  The structure of the problem should not change if these
    // auxiliary variables are included (hence, StateSelect.never).

    Connectors.Face faces[n_trans, Side](
      each p(start=p_IC),
      each T(start=T_IC),
      each Ndot(start=0,stateSelect=StateSelect.never))
      "Connectors for transport" annotation (Placement(transformation(extent={{
              -10,-10},{10,10}}), iconTransformation(extent={{-10,-10},{10,10}})));
    Connectors.Electrochemical electrochemical(
      final n_trans=n_trans,
      w(start=g_IC, final fixed=false),
      sT(start=h_IC - g_IC, final fixed=false))
      "Connector for reactions and phase change" annotation (Placement(
          transformation(extent={{-10,30},{10,50}}), iconTransformation(extent=
              {{-30,80},{-50,100}})));

    // Geometric parameters
  protected
    outer Q.Length Lprime[:]
      "Effective area divided by length of the transport axes" annotation (
        missingInnerMessage="This model should be used within a phase model.");
    outer parameter Boolean inclRot[3]
      "true, if each axis of rotation has all its tangential faces included"
      annotation (missingInnerMessage="This model should be used within a subregion model.
");
    outer parameter Boolean inclTrans[3]
      "true, if each transport axis is included" annotation (
        missingInnerMessage="This model should be used within a subregion model.
");
    // Note:  The size of inclRot and inclTrans is also Axis, but it isn't
    // specified here due to an error in Dymola 2014.
    outer parameter Integer cartRot[:]
      "Cartesian-axis indices of the components of rotational momentum"
      annotation (missingInnerMessage="This model should be used within a subregion model.
");
    outer parameter Integer cartTrans[:]
      "Cartesian-axis indices of the transport axes" annotation (
        missingInnerMessage="This model should be used within a subregion model.
");
    // Note:  The size of cartRot and cartTrans is n_trans,
    // but it isn't specified here due to an error in Dymola 2014.
    outer parameter Integer transCart[3]
      "Face-pair indices of the Cartesian axes" annotation (missingInnerMessage
        ="This model should be used within a subregion model.
");
    // Note:  The size is also Axis, but it isn't specified here due
    // to an error in Dymola 2014.
    final parameter Boolean upstream[n_trans]=selectBooleans({upstreamX,
        upstreamY,upstreamZ}, cartTrans)
      "true, if each transport axis uses upstream discretization"
      annotation (HideResult=true);
    final parameter ConsMom consTrans[n_trans]=selectIntegers({consTransX,
        consTransY,consTransZ}, cartTrans)
      "Formulation of the translational conservation equations for the transport axes"
      annotation (HideResult=true);
    outer parameter Q.NumberAbsolute k[Axis]
      "Scaling factor for diffusive transport" annotation (missingInnerMessage=
          "This model should be used within a phase model.");

    // Additional aliases (for common terms)
    Q.Force mPhidot[n_trans](each nominal=U.N, each stateSelect=StateSelect.never)
      "Force, excluding effects of thermodynamic pressure, dynamic pressure, and unsteady current";
    Q.Force faces_mPhidot[n_trans, Side, Orient](each nominal=U.N, each
        stateSelect=StateSelect.never) "Directly-calculated shear forces";

    outer Conditions.Environment environment "Environmental conditions";

  initial equation
    // Check the initial conditions.
    assert(V >= 0, "The volume of " + getInstanceName() + " is negative.
Check that the volumes of the other phases are set properly.");
    assert(initMaterial <> initEnergy or initMaterial == Init.none or
      consMaterial == ConsThermo.steady or consEnergy == ConsThermo.steady,
      "The initialization methods for material and energy must be different (unless None).");

    // Material
    if consMaterial == ConsThermo.IC then
      // Ensure that a condition is selected since the state is prescribed.
      assert(initMaterial <> Init.none, "The material state of " +
        getInstanceName() + " is prescribed, yet its condition is not defined.
Choose any condition besides None.");
    elseif consMaterial == ConsThermo.dynamic then
      // Initialize since there's a time-varying state.
      if initMaterial == Init.amount then
        N = N_IC;
      elseif initMaterial == Init.amountSS then
        der(N) = 0;
      elseif initMaterial == Init.density then
        1/v = rho_IC;
        assert(Data.isCompressible or Data.hasThermalExpansion, getInstanceName()
           +
          " is isochoric, yet its material initial condition is based on density.");
      elseif initMaterial == Init.densitySS then
        der(1/v) = 0;
        assert(Data.isCompressible or Data.hasThermalExpansion, getInstanceName()
           +
          " is isochoric, yet its material initial condition is based on density.");
      elseif initMaterial == Init.volume then
        V = V_IC;
      elseif initMaterial == Init.volumeSS then
        der(V) = 0;
      elseif initMaterial == Init.volumeSS then
        der(V) = 0;
      elseif initMaterial == Init.pressure then
        p = p_IC;
        assert(Data.isCompressible, getInstanceName() +
          " is incompressible, yet its material initial condition is based on pressure.");
      elseif initMaterial == Init.pressureSS then
        der(p) = 0;
        assert(Data.isCompressible, getInstanceName() +
          " is incompressible, yet its material initial condition is based on pressure.");
      elseif initMaterial == Init.temperature then
        T = T_IC;
      elseif initMaterial == Init.temperatureSS then
        der(T) = 0;
      elseif initMaterial == Init.specificEnthalpy then
        h = h_IC;
      elseif initMaterial == Init.specificEnthalpySS then
        der(h) = 0;
      elseif initMaterial == Init.Gibbs then
        g = g_IC;
      elseif initMaterial == Init.GibbsSS then
        der(g) = 0;
        // Else, there's no initial equation since
        // initMaterial == Init.none or
        // consMaterial == ConsThermo.steady.
      end if;
    end if;

    // Energy
    if consEnergy == ConsThermo.IC then
      // Ensure that a condition is selected since the state is prescribed.
      assert(initEnergy <> Init.none, "The energy state of " + getInstanceName()
         + " is prescribed, yet its condition is not defined.
Choose any condition besides None.");
    elseif consEnergy == ConsThermo.dynamic then
      // Initialize since there's a time-varying state.
      if initEnergy == Init.amount then
        N = N_IC;
      elseif initEnergy == Init.amountSS then
        der(N) = 0;
      elseif initEnergy == Init.density then
        1/v = rho_IC;
        assert(Data.isCompressible or Data.hasThermalExpansion, getInstanceName()
           +
          " is isochoric, yet its thermal initial condition is based on density.");
      elseif initEnergy == Init.densitySS then
        der(1/v) = 0;
        assert(Data.isCompressible or Data.hasThermalExpansion, getInstanceName()
           +
          " is isochoric, yet its thermal initial condition is based on density.");
      elseif initEnergy == Init.volume then
        V = V_IC;
      elseif initEnergy == Init.volumeSS then
        der(V) = 0;
      elseif initEnergy == Init.pressure then
        p = p_IC;
        assert(Data.isCompressible, getInstanceName() +
          " is incompressible, yet its thermal initial condition is based on pressure.");
      elseif initEnergy == Init.pressureSS then
        der(p) = 0;
        assert(Data.isCompressible, getInstanceName() +
          " is incompressible, yet its thermal initial condition is based on pressure.");
      elseif initEnergy == Init.temperature then
        T = T_IC;
      elseif initEnergy == Init.temperatureSS then
        der(T) = 0;
      elseif initEnergy == Init.specificEnthalpy then
        h = h_IC;
      elseif initEnergy == Init.specificEnthalpySS then
        der(h) = 0;
      elseif initEnergy == Init.Gibbs then
        g = g_IC;
      elseif initEnergy == Init.GibbsSS then
        der(g) = 0;
        // Else, there's no initial equation since
        // initEnergy == Init.none or
        // consEnergy == ConsThermo.steady.
      end if;
    end if;

  equation
    // Aliases (only to clarify and simplify other equations)
    I .* L[cartTrans] = N*phi;
    phi_faces = faces.Ndot .* Data.v_Tp(faces.T, faces.p) ./ {A[cartTrans[i]]*{
      1,-1} for i in 1:n_trans};
    mPhidot + M*environment.a[cartTrans] + Data.z*N*environment.E[cartTrans] =
      Data.m*actualStream(electrochemical.phi)*electrochemical.Ndot + direct.translational.mPhidot
       + {sum(intra[:].mPhidot[j]) + sum(inter[:].mPhidot[j]) + sum((if i == j
       then 0 else faces[i, :].phi[cartWrap(cartTrans[j] - cartTrans[i])]*faces[
      i, :].Ndot*Data.m + sum(faces[i, :].mPhidot[cartWrap(cartTrans[j] -
      cartTrans[i])])) for i in 1:n_trans) for j in 1:n_trans};

    // Properties upon outflow due to reaction and phase change
    electrochemical.phi = phi;
    electrochemical.sT = h - g;

    // Material exchange
    if abs(alpha - 0.5) > Modelica.Constants.eps then
      v*tauprime*electrochemical.Ndot = V*(exp(alpha*(electrochemical.w - g)/T)
         - exp((alpha - 1)*(electrochemical.w - g)/T));
      // Note:  V/v is different than N for the Ion model, which inherits from
      // this one.
    else
      electrochemical.w = g + 2*T*asinh(electrochemical.Ndot*tauprime*v/(2*V))
        "Explicitly inverted to avoid nonlinear system of equations";
    end if;

    // Transport
    for i in 1:n_trans loop
      for side in Side loop
        // Material/conservation of translational momentum
        (if consTrans[i] == ConsMom.dynamic then Data.m*der(faces[i, side].Ndot)
          /U.s else 0) = Lprime[i]*(faces[i, side].p - p)*2 + Data.m*faces[i,
          side].Ndot^2/N + inSign(side)*mPhidot[i]/L[cartTrans[i]];

        // Transverse translational momentum
        eta*faces_mPhidot[i, side, Orient.after] = Nu_Phi[after(cartTrans[i])]*
          Lprime[i]*(faces[i, side].phi[Orient.after] - (if inclTrans[after(
          cartTrans[i])] then phi[transCart[after(cartTrans[i])]] else 0))*(if
          upstream[i] then 1 + exp(-eta*Lprime[i]*faces[i, side].Ndot/2) else 2)
          "1st transverse";
        eta*faces_mPhidot[i, side, Orient.before] = Nu_Phi[before(cartTrans[i])]
          *Lprime[i]*(faces[i, side].phi[Orient.before] - (if inclTrans[before(
          cartTrans[i])] then phi[transCart[before(cartTrans[i])]] else 0))*(
          if upstream[i] then 1 + exp(-eta*Lprime[i]*faces[i, side].Ndot/2)
           else 2) "2nd transverse";

        // Thermal energy
        theta*faces[i, side].Qdot = Nu_Q*Lprime[i]*(faces[i, side].T - T)*(if
          upstream[i] then 1 + exp(-theta*Lprime[i]*Data.c_v(T, p)*faces[i,
          side].Ndot/2) else 2);
      end for;

      // Direct mapping of shear forces (calculated above)
      if not (consRot and inclRot[before(cartTrans[i])]) then
        faces[i, :].mPhidot[Orient.after] = faces_mPhidot[i, :, Orient.after];
        // Else, the force must be mapped for zero torque (below).
      end if;
      if not (consRot and inclRot[after(cartTrans[i])]) then
        faces[i, :].mPhidot[Orient.before] = faces_mPhidot[i, :, Orient.before];
        // Else, the force must be mapped for zero torque (below).
      end if;
    end for;

    // Zero-torque mapping of shear forces
    if consRot then
      for axis in cartRot loop
        4*cat(
            1,
            faces[transCart[after(axis)], :].mPhidot[Orient.after],
            faces[transCart[before(axis)], :].mPhidot[Orient.before]) = {{3,1,L[
          before(axis)]/L[after(axis)],-L[before(axis)]/L[after(axis)]},{1,3,-L[
          before(axis)]/L[after(axis)],L[before(axis)]/L[after(axis)]},{L[after(
          axis)]/L[before(axis)],-L[after(axis)]/L[before(axis)],3,1},{-L[after(
          axis)]/L[before(axis)],L[after(axis)]/L[before(axis)],1,3}}*cat(
            1,
            faces_mPhidot[transCart[after(axis)], :, Orient.after],
            faces_mPhidot[transCart[before(axis)], :, Orient.before]);
      end for;
    end if;

    // Material dynamics
    if consMaterial == ConsThermo.IC then
      // Apply the IC forever (material not conserved).
      if initMaterial == Init.amount then
        N = N_IC;
      elseif initMaterial == Init.amountSS then
        der(N) = 0;
      elseif initMaterial == Init.density then
        1/v = rho_IC;
      elseif initMaterial == Init.densitySS then
        der(1/v) = 0;
      elseif initMaterial == Init.volume then
        V = V_IC;
      elseif initMaterial == Init.volumeSS then
        der(V) = 0;
      elseif initMaterial == Init.pressure then
        p = p_IC;
      elseif initMaterial == Init.pressureSS then
        der(p) = 0;
      elseif initMaterial == Init.temperature then
        T = T_IC;
      elseif initMaterial == Init.temperatureSS then
        der(T) = 0;
      elseif initMaterial == Init.specificEnthalpy then
        h = h_IC;
      elseif initMaterial == Init.specificEnthalpySS then
        der(h) = 0;
      elseif initMaterial == Init.Gibbs then
        g = g_IC;
      else
        // if initMaterial == Init.GibbsSS then
        der(g) = 0;
        // Note:  initMaterial == Init.none can't occur due to an
        // assertion.
      end if;
    else
      (if consMaterial == ConsThermo.dynamic then der(N)/U.s else 0) =
        electrochemical.Ndot + sum(faces.Ndot) "Material conservation";
    end if;

    // Linear current profile
    for i in 1:n_trans loop
      0 = 2*I[i] + Delta(faces[i, :].Ndot);
    end for;

    // Thermal dynamics
    if consEnergy == ConsThermo.IC then
      // Apply the IC forever (energy not conserved).
      if initEnergy == Init.amount then
        N = N_IC;
      elseif initEnergy == Init.amountSS then
        der(N) = 0;
      elseif initEnergy == Init.density then
        1/v = rho_IC;
      elseif initEnergy == Init.densitySS then
        der(1/v) = 0;
      elseif initEnergy == Init.volume then
        V = V_IC;
      elseif initEnergy == Init.volumeSS then
        der(V) = 0;
      elseif initEnergy == Init.pressure then
        p = p_IC;
      elseif initEnergy == Init.pressureSS then
        der(p) = 0;
      elseif initEnergy == Init.temperature then
        T = T_IC;
      elseif initEnergy == Init.temperatureSS then
        der(T) = 0;
      elseif initEnergy == Init.specificEnthalpy then
        h = h_IC;
      elseif initEnergy == Init.specificEnthalpySS then
        der(h) = 0;
      elseif initEnergy == Init.Gibbs then
        g = g_IC;
      else
        // if initEnergy == Init.GibbsSS then
        der(g) = 0;
        // Note:  initEnergy == Init.none can't occur due to an
        // assertion.
      end if;
    else
      (if consEnergy == ConsThermo.dynamic then (N*T*der(Data.s(T, p)) + der(M*
        phi*phi)/2)/U.s else 0) = (electrochemical.w + actualStream(
        electrochemical.sT) - h + actualStream(electrochemical.phi)*
        actualStream(electrochemical.phi)*Data.m/2)*electrochemical.Ndot +
        direct.translational.phi*direct.translational.mPhidot + sum(intra[i].phi
        *intra[i].mPhidot for i in 1:n_intra) + sum(inter[i].phi*inter[i].mPhidot
        for i in 1:n_inter) + direct.thermal.Qdot + sum(intra.Qdot) + sum(inter.Qdot)
         + sum((Data.h(faces[i, :].T, faces[i, :].p) - {h,h} + (phi_faces[i, :]
         .^ 2 + sum(faces[i, :].phi[orient] .^ 2 for orient in Orient))*(Data.m
        /2))*faces[i, :].Ndot + sum(faces[i, :].phi[orient]*faces[i, :].mPhidot[
        orient] for orient in Orient) for i in 1:n_trans) + sum(faces.Qdot)
        "Conservation of energy";
      // Note:  In Dymola 2014, der(Data.s(T, p)) is better than der(s) because it avoids
      // dynamic state selection.  Dymola sometimes chooses s as a state even though its
      // stateSelect is StateSelect.never.
    end if;
    annotation (
      defaultComponentPrefixes="replaceable",
      Documentation(info="<html>
    <p>Notes regarding the parameters:
    <ol>    
    <li>If consRot is <code>true</code>, then rotational momentum is conserved without storage
    (i.e., steady).  This means that the shear forces are mapped so that there is no net torque around any
    rotational axis that has all its faces included (i.e., all the faces around the perimeter).  Rotational 
    momentum is not exchanged among species or directly transported (i.e., uniform or shaft rotation).</li></ol></p>
    
    <p>For more information, please see the
     <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p>
    </html>"),
      Diagram(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          initialScale=0.1), graphics),
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
              100}}), graphics));
  end Fluid;

  model Solid "Base model for a solid species"

    // Geometry
    parameter Q.Number epsilon=0.25 "Volumetric fill fraction" annotation (
        Dialog(group="Geometry", __Dymola_label="<html>&epsilon;</html>"));
    extends Species(
      final N_IC,
      final V_IC=epsilon*product(L),
      final rho_IC=1/Data.v_Tp(),
      final p_IC=environment.p,
      final h_IC,
      final g_IC);

    // Material properties
    Q.ResistivityThermal theta(nominal=10*U.cm/U.A) = Data.theta(T, v)
      "Thermal resistivity" annotation (Dialog(group="Material properties",
          __Dymola_label="<html>&theta;</html>"));

    // Assumptions
    parameter ConsThermo consEnergy=ConsThermo.dynamic "Energy" annotation (
        Evaluate=true, Dialog(tab="Assumptions", group=
            "Formulation of the conservation equations"));

    Connectors.ThermalDiffusion faces[n_trans, Side](T(each start=T_IC))
      "Connectors for transport" annotation (Placement(transformation(extent={{
              -10,-10},{10,10}}), iconTransformation(extent={{-10,-10},{10,10}})));

  protected
    outer Q.Length Lprime[:]
      "Effective area divided by length of the transport axes" annotation (
        missingInnerMessage="This model should be used within a phase model.");
    outer parameter Integer cartTrans[:]
      "Cartesian-axis indices of the pairs of faces" annotation (
        missingInnerMessage="This model should be used within a subregion model.
");

  initial equation
    if consEnergy == ConsThermo.dynamic then
      T = T_IC;
    end if;

  equation
    // Assumptions
    phi = zeros(n_trans) "Zero velocity";
    V = epsilon*product(L) "Prescribed volume";

    // Thermal transport (conduction)
    for i in 1:n_trans loop
      for side in Side loop
        faces[i, side].Qdot = 2*Lprime[i]*(faces[i, side].T - T)/theta;
      end for;
    end for;

    // Thermal dynamics
    if consEnergy == ConsThermo.IC then
      // Apply the IC forever (energy not conserved).
      T = T_IC;
    else
      (if consEnergy == ConsThermo.dynamic then N*T*der(s)/U.s else 0) = sum(
        intra[i].phi*intra[i].mPhidot for i in 1:n_intra) + sum(inter[i].phi*
        inter[i].mPhidot for i in 1:n_inter) + direct.thermal.Qdot + sum(intra.Qdot)
         + sum(inter.Qdot) + sum(faces.Qdot) "Conservation of energy";
    end if;

    annotation (
      defaultComponentPrefixes="replaceable",
      Documentation(info="<html><p>Assumptions:<ol>
  <li>There is no material transport.</li>
  <li>Velocity is zero.</li>
  <li>There are no chemical reactions or phase change.</li>
  <li>The volume is constant (determined by &epsilon; and the total volume of the subregion).</li>
  </ol></p>

  <p>For more information, please see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"),

      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
              100,100}}), graphics));

  end Solid;

protected
  partial model Species "Base model for one chemical species in one phase"

    import assert = FCSys.Utilities.assertEval;
    //extends FCSys.Icons.Names.Top5;

    // Geometry
    parameter Q.NumberAbsolute k_intra[n_intra]=ones(n_intra)
      "Coupling factors within the phase" annotation (Dialog(group="Geometry",
          __Dymola_label="<html><i><b>k</b></i><sub>intra</sub></html>"));

    // Material properties
    replaceable package Data = Characteristics.BaseClasses.Characteristic
      constrainedby Characteristics.BaseClasses.Characteristic
      "Characteristic data" annotation (
      Evaluate=true,
      Dialog(group="Material properties"),
      choicesAllMatching=true,
      __Dymola_choicesFromPackage=true);
    Q.Mobility mu(nominal=0.1*U.C*U.s/U.kg) = Data.mu(T, v) "Mobility"
      annotation (Dialog(group="Material properties", __Dymola_label=
            "<html>&mu;</html>"));
    Q.TimeAbsolute nu(nominal=1e-9*U.s) = Data.nu(T, v) "Thermal independity"
      annotation (Dialog(group="Material properties", __Dymola_label=
            "<html>&nu;</html>"));

    // Assumptions
    parameter Integer n_trans=1 "Number of transport axes" annotation (
        HideResult=true, Dialog(tab="Assumptions", __Dymola_label=
            "<html><i>n</i><sub>trans</sub></html>"));
    // This can't be an outer parameter in Dymola 2014.
    parameter Integer n_intra=0
      "Number of exchange connections within the phase"
      annotation (HideResult=true,Dialog(connectorSizing=true));
    parameter Integer n_inter=0
      "Number of exchange connections with other phases" annotation (Dialog(
          __Dymola_label="<html><i>n</i><sub>interM/sub></html>"), HideResult=
          true);

    // Initialization parameters
    parameter Q.Amount N_IC(start=V_IC*rho_IC) "Initial amount of material"
      annotation (Dialog(tab="Initialization", __Dymola_label=
            "<html><i>N</i><sub>IC</sub></html>"));
    parameter Q.Density rho_IC(start=1/Data.v_Tp(T_IC, p_IC)) "Initial density"
      annotation (Dialog(tab="Initialization", __Dymola_label=
            "<html>&rho;<sub>IC</sub></html>"));
    parameter Q.Volume V_IC(start=product(L)) "Initial volume" annotation (
        Dialog(tab="Initialization", __Dymola_label=
            "<html><i>V</i><sub>IC</sub></html>"));
    parameter Q.PressureAbsolute p_IC(start=environment.p) "Initial pressure"
      annotation (Dialog(tab="Initialization", __Dymola_label=
            "<html><i>p</i><sub>IC</sub></html>"));
    parameter Q.TemperatureAbsolute T_IC(start=environment.T)
      "Initial temperature" annotation (Dialog(tab="Initialization",
          __Dymola_label="<html><i>T</i><sub>IC</sub></html>"));
    parameter Q.Potential h_IC(start=Data.h(T_IC, p_IC), displayUnit="kJ/mol")
      "Initial specific enthalpy" annotation (Dialog(tab="Initialization",
          __Dymola_label="<html><i>h</i><sub>IC</sub></html>"));
    parameter Q.Potential g_IC(start=Data.g(T_IC, p_IC), displayUnit="kJ/mol")
      "Initial Gibbs potential" annotation (Dialog(tab="Initialization",
          __Dymola_label="<html><i>g</i><sub>IC</sub></html>"));

    // Preferred states
    // Note:  The start values for these variable aren't fixed because
    // the initial equation section will be used instead.
    Q.Amount N(
      final min=Modelica.Constants.small,
      nominal=4*U.C,
      final start=N_IC,
      final fixed=false,
      stateSelect=StateSelect.prefer) "Amount of material";
    Q.TemperatureAbsolute T(
      nominal=300*U.K,
      final start=T_IC,
      final fixed=false,
      stateSelect=StateSelect.prefer) "Temperature";
    Q.Velocity phi[n_trans](
      each nominal=10*U.cm/U.s,
      each start=0,
      each stateSelect=StateSelect.prefer) "Velocity";

    // Aliases (for common terms)
    // Note:  StateSelect.never helps avoid dynamic state selection of these
    // variables in Dymola 2014.
    Q.PressureAbsolute p(
      nominal=U.atm,
      final start=p_IC,
      final fixed=false,
      stateSelect=StateSelect.never) "Pressure";
    Q.Potential g(
      nominal=U.V,
      final start=g_IC,
      final fixed=false,
      stateSelect=StateSelect.never) "Specific Gibbs energy";
    Q.Mass M(
      nominal=1e-3*U.g,
      final start=Data.m*N_IC,
      stateSelect=StateSelect.never) "Mass";
    Q.VolumeSpecific v(
      nominal=U.cc/(4*U.C),
      final start=1/rho_IC,
      final fixed=false,
      stateSelect=StateSelect.never) "Specific volume";
    Q.Potential h(
      nominal=U.V,
      final start=h_IC,
      final fixed=false,
      stateSelect=StateSelect.never) "Specific enthalpy";
    Q.NumberAbsolute s(
      nominal=10,
      final start=(h_IC - g_IC)/T_IC,
      stateSelect=StateSelect.never) "Specific entropy";

    // Auxiliary variables (for analysis)
    // ----------------------------------
    // Thermodynamic properties
    output Q.Density rho(stateSelect=StateSelect.never) = 1/v if environment.analysis
      "Density";
    output Q.MassVolumic mrho(stateSelect=StateSelect.never) = Data.m*rho if
      environment.analysis "Volumic mass";
    output Q.Amount S(stateSelect=StateSelect.never) = N*s if environment.analysis
      "Entropy";
    output Q.CapacityThermalSpecific c_p(stateSelect=StateSelect.never) =
      Data.c_p(T, p) if environment.analysis "Isobaric specific heat capacity";
    output Q.CapacityThermalSpecific c_v(stateSelect=StateSelect.never) =
      Data.c_v(T, p) if environment.analysis "Isochoric specific heat capacity";
    output Q.PressureReciprocal beta(stateSelect=StateSelect.never) = Data.beta(
      T, p) if environment.analysis "Isothermal compressibility";
    //
    // Time constants
    output Q.TimeAbsolute tau_PhiE_intra[n_intra](
      each stateSelect=StateSelect.never,
      each start=U.s) = fill(Data.m*mu, n_intra) ./ k_intra if environment.analysis
      "Time constant for translational intra-phase exchange";
    output Q.TimeAbsolute tau_PhiE_inter[n_inter](
      each stateSelect=StateSelect.never,
      each start=U.s) = fill(Data.m*mu, n_inter) ./ k_inter if environment.analysis
       and n_inter > 0 "Time constant for translational inter-phase exchange";
    output Q.TimeAbsolute tau_QE_intra[n_intra](
      each stateSelect=StateSelect.never,
      each start=U.s) = fill(c_p*nu, n_intra) ./ k_intra if environment.analysis
      "Time constant for thermal intra-phase exchange";
    output Q.TimeAbsolute tau_QE_inter[n_inter](
      each stateSelect=StateSelect.never,
      each start=U.s) = fill(c_p*nu, n_inter) ./ k_inter if environment.analysis
       and n_inter > 0 "Time constant for thermal inter-phase exchange";
    // Note:  The time contants for the direct connector are zero because it has
    // no resistance.
    // Note:  The structure of the problem should not change if these
    // auxiliary variables are included (hence, StateSelect.never).

    Connectors.Direct direct(
      final n_trans=n_trans,
      final inclTrans=true,
      final inclThermal=true,
      thermal(T(final start=T_IC, final fixed=false)))
      "Connector to directly couple velocity or temperature with other species"
      annotation (Placement(transformation(extent={{-10,-50},{10,-30}}),
          iconTransformation(extent={{28,-80},{48,-100}})));
    Connectors.Intra intra[n_intra](each final n_trans=n_trans, each T(final
          start=T_IC, final fixed=false))
      "Connectors to exchange translational momentum and energy within the phase"
      annotation (Placement(transformation(extent={{10,-30},{30,-10}}),
          iconTransformation(extent={{60,-60},{80,-80}})));
    Connectors.Inter inter[n_inter](each final n_trans=n_trans, each T(final
          start=T_IC, final fixed=false))
      "Connectors to exchange translational momentum and energy with all other species"
      annotation (Placement(transformation(extent={{30,-10},{50,10}}),
          iconTransformation(extent={{80,-28},{100,-48}})));

    Connectors.Dalton dalton(V(
        min=0,
        final start=V_IC,
        final fixed=false) = v*N, p(final start=p_IC, final fixed=false))
      "Connector for additivity of pressure" annotation (Placement(
          transformation(extent={{-50,-10},{-30,10}}), iconTransformation(
            extent={{-100,50},{-80,30}})));

    // Geometric parameters
  protected
    outer Q.Volume V "Volume of the phase (not of the subregion)" annotation (
        missingInnerMessage="This model should be used within a phase model.
");
    outer parameter Q.Length A[Axis] "Cross-sectional areas of the subregion"
      annotation (missingInnerMessage="This model should be used within a subregion model.
");
    outer parameter Q.Length L[Axis] "Lengths of the subregion" annotation (
        missingInnerMessage="This model should be used within a subregion model.
");
    outer parameter Q.NumberAbsolute k_inter[:]
      "Coupling factor for diffusive exchange with other phases" annotation (
        missingInnerMessage="This model should be used within a phase model");

    outer Conditions.Environment environment "Environmental conditions";

  initial equation
    // Check the initial conditions.
    assert(V >= 0, "The volume of " + getInstanceName() + " is negative.
Check that the volumes of the other phases are set properly.");

  equation
    // Aliases (only to clarify and simplify other equations)
    p = dalton.p;
    T = direct.thermal.T;
    h = g + T*s;
    M = Data.m*N;
    phi = direct.translational.phi;

    // Thermodynamic correlations
    if Data.isCompressible then
      p = Data.p_Tv(T, v);
    else
      v = Data.v_Tp(T, p);
    end if;
    h = Data.h(T, p);
    s = Data.s(T, p);

    // Exchange
    // --------
    // Translational momentum
    for i in 1:n_trans loop
      v*mu*intra.mPhidot[i] = k_intra*V .* (intra.phi[i] - fill(phi[i], n_intra));
      if n_inter > 0 then
        v*mu*inter.mPhidot[i] = k_inter*V .* (inter.phi[i] - fill(phi[i],
          n_inter));
      end if;
    end for;
    //
    // Thermal energy
    v*nu*intra.Qdot = k_intra*V .* (intra.T - fill(T, n_intra));
    if n_inter > 0 then
      v*nu*inter.Qdot = k_inter*V .* (inter.T - fill(T, n_inter));
    end if;
    // Note:  V/v is different than N for the Ion model, which inherits from
    // this one.
    annotation (
      defaultComponentPrefixes="replaceable",
      Documentation(info="<html>
    <p>All of the details below are pertinent to the <a href=\"modelica://FCSys.Fluid\">Fluid</a>
    model (and the derived <a href=\"modelica://FCSys.Gas\">Gas</a> and 
    <a href=\"modelica://FCSys.Liquid\">Liquid</a> models) which inherits from this model.  
    Only some of the details apply to the 
    <a href=\"modelica://FCSys.Ion\">Ion</a> model because it excludes shear forces between neighboring 
    regions.
    The <a href=\"modelica://FCSys.Solid\">Solid</a> also excludes the transport and exchange of 
    material.</p>
    
    <p>This model is based on the following fixed assumptions:
    <ol>
       <li>All faces are rectangular.
       <li>The material is orthorhombic.  This implies that a gradient which induces diffusion
       along an axis does not induce diffusion along axes orthogonal to it
       [<a href=\"modelica://FCSys.UsersGuide.References\">Bejan2006</a>,
       pp. 691&ndash;692].</li>
       <li>The coordinate system (x, y, z) is aligned with the principle
       axes of transport.  For example, if the material is stratified, then the
       layers must be parallel to one of the planes in the rectilinear
       grid.</li>
       <li>The factors that describe anisotropic effects (<b><i>k</i></b>)
          are common to material, translational, and thermal transport for every species 
          in the phase.</li>
       <li>There is no radiative heat transfer (or else it must be linearized).</li>
       <li>Rotational momentum is not exchanged, transported, or stored.</li>
    </ol>
    Other assumptions are optional via the parameters.  Additional assumptions may be 
    applied in models that inherit from this one.</p>

    <p><a href=\"#Fig1\">Figure 1</a> shows how instances of
    <a href=\"modelica://FCSys.Species\">Species</a> models are
    connected within a <a href=\"modelica://FCSys.Subregions\">Subregion</a>.  A single species in
    a single phase is called a <i>configuration</i>. The
    generalized resistances (<i>R</i>) affect the force and rates of chemical exchange and 
    heat flow
    associated with differences in activity, velocity, and temperature (respectively) between
    each configuration and a common node.  These exchange processes are diffusive.
    The resistors generate heat 
    in the <a href=\"modelica://FCSys.Species\">Species</a> instance within which they are
    included.</p>
    
    <p align=center id=\"Fig1\"><img src=\"modelica://FCSys/Resources/Documentation/Subregions/Species/Species/Exchange.png\">
<br>Figure 1:  Exchange of a quantity (material, translational momentum, or thermal energy) among configurations
    (A, B, and C) within a subregion.</p>

    <p>Translational momentum and thermal energy are advected as material is exchanged
    due to phase change or reactions.  This occurs at the velocity (&phi;) and the specific entropy-temperature
    product (<i>sT</i>) of the reactants (source configurations), where the reactant/product designation
    depends on the current conditions.</p>

    <p>The advective exchange is modeled using <code>stream</code> connectors
    (<a href=\"modelica://FCSys.Connectors.Physical\">Physical</a> and
    <a href=\"modelica://FCSys.Connectors.Electrochemical\">Electrochemical</a>).
  The rate of advection of translational momentum is the
  product of the velocity of the source (&phi;) and the mass flow rate
  (<i>M&#775;</i> or <i>m</i><i>N&#775;</i>).  The rate of thermal advection is the
  specific entropy-temperature product of the source (<i>sT</i>) times the rate of
  material exchange
  (<i>N&#775;</i>).  If there are multiple sources, then
  their contributions are additive.  If there are multiple sinks, then
  translational momentum is split on a mass basis and the thermal stream is split
  on a particle-number basis.</p>
  
    <p><a href=\"#Fig2\">Figure 2</a> shows how
    a configuration
    is connected between neighboring instances of a
    <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a>.
    Material, translational momentum, and thermal energy are transported by both advection and diffusion.
    Upstream discretization is applied if it is enabled via the <code>upstreamX</code>,
    etc. parameters.  Like for exchange, the transport resistances are inside the
    <a href=\"modelica://FCSys.Species\">Species</a> model.</p>

    <p align=center id=\"Fig2\"><img src=\"modelica://FCSys/Resources/Documentation/Subregions/Species/Species/Transport.png\">
<br>Figure 2:  Transport of a quantity associated with the same configuration
    between subregions (1 and 2).</p>

<p>The <a href=\"modelica://FCSys.Species\">Species</a> instances
    within a <a href=\"modelica://FCSys.Phases\">Phase</a> are combined by Dalton's law of
    partial pressures (see the
    <a href=\"modelica://FCSys.Connectors.Dalton\">Dalton</a> connector), as shown
    in Figure 3a.  The pressures are additive, and each species is assumed to exist at the
    total extensive volume of the phase.  Within a 
    <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a>,
    the <a href=\"modelica://FCSys.Phases\">Phases</a> are combined by Amagat's law of partial volumes
    (see the <a href=\"modelica://FCSys.Connectors.Amagat\">Amagat</a> connector), as shown
    in Figure 3b.  The volumes are additive, and each species is assumed to exist at the
    total pressure in the subregion.</p>

    <table border=0 cellspacing=0 cellpadding=2 align=center class=noBorder style=\"margin-left: auto; margin-right: auto;\">
      <tr align=center class=noBorder>
        <td align=center class=noBorder style=\"margin-left: auto; margin-right: auto;\">
          <img src=\"modelica://FCSys/Resources/Documentation/Subregions/Species/Species/SharePressure.png\">
<br>a:  Pressures of species (A, B, and C) are additive within a phase.
        </td>
        <td align=center class=noBorder style=\"margin-left: auto; margin-right: auto;\">
          <img src=\"modelica://FCSys/Resources/Documentation/Subregions/Species/Species/ShareVolume.png\">
<br>b:  Volumes of phases (I, II, and III) are additive within a subregion.
        </td>
      </tr>
      <tr align=center class=noBorder style=\"margin-left: auto; margin-right: auto;\">
        <td colspan=2 align=center class=noBorder>Figure 3: Methods of attributing pressure and volume.</td>
      </tr>
    </table>

    <p>Notes regarding the parameters:
    <ol>
    <li>In general, if the interval for chemical exchange (&tau;&prime;), mobility (&mu;),
     thermal independity (&nu;), fluidity (&eta;), or thermal resistivity (&theta;)
      is zero, then
    it should be set as <code>final</code> so that index reduction may be performed.
    If two configurations
    are connected through their <code>intra</code>, <code>inter</code>, or <code>faces</code>
    and both have zero generalized resistivities for a
    quantity, then index reduction [<a href=\"modelica://FCSys.UsersGuide.References\">Mattsson1993B</a>] is necessary.</li>
    
    <li>Even if an initialization parameter is not selected for explicit use,
    it may be used a guess value.</li>
    
    <li>If <code>ConsThermo.IC</code> is used for a state (via
    <code>consMaterial</code> or <code>consEnergy</code>),
    then the associated initial condition (IC) will be applied forever instead of the
    corresponding conservation equation.</li>
    
    <li>If <code>consTransX</code>, <code>consTransY</code>, or <code>consTransZ</code> is
    <code>ConsMom.steady</code>, then the derivatives of the corresponding boundary currents 
    are treated as zero and removed from the translational momentum balances at the associated boundaries.  
    If <code>consEnergy</code> is
    <code>ConsThermo.steady</code>, then <i>T</i>&part;<i>s</i>/&part;<i>t</i> + <i>M</i>&phi;&part;&phi;/&part;<i>t</i> is treated as
    zero and removed from the energy balance.</li>
    
    <li>If a transport axis is not included (via the outer <code>inclTrans[:]</code> parameter
    which maps to <code>{inclTransX, inclTransY, inclTransZ}</code> in the
    <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a> model), then the velocity 
    along that axis is considered to be zero.</li>
    
    <li>The <code>start</code> values of the initial conditions for pressure and temperature
    (<i>p</i><sub>IC</sub> and <i>T</i><sub>IC</sub>) are the global default pressure and
    temperature (via the <code>outer</code> instance of the <a href=\"modelica://FCSys.Conditions.Environment\">Environment</a> model).
    The <code>start</code> values of the initial conditions for
    other intensive properties (&rho;<sub>IC</sub>, <i>h</i><sub>IC</sub>, and
    <i>g</i><sub>IC</sub>) are related to the initial pressure and temperature
    by the characteristics of the species.  The <code>start</code> value of the
    initial condition for the extensive volume (<i>V</i><sub>IC</sub>) is the volume of the
    subregion.  The <code>start</code> value for particle number (<i>N</i><sub>IC</sub>)
    is related to it via the material characteristics and the initial pressure and temperature.
    In order to apply other values for any of these initial conditions,
    it may be necessary to do so before translating the model.</li>
    
    <li>Upstream discretization may be applied to translational and thermal transport
    using (<code>upstreamX=true</code>, etc.).  Otherwise, the central difference
    scheme is used.  The central difference scheme
    is always used for material diffusion.</li>
        
    <li>The default thermal Nusselt number is one, which represents pure conduction through the gas.  Use 
    3.66 for internal flow where the boundaries are uniform in temperature or 48/11 or approximately 4.36 
    if the heat flux is uniform [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>].</li>
    
    <li>The indices of the translational Nusselt number (<i>Nu</i><sub>&Phi;</sub>)
    correspond to the axes of material advection, not the axes of
    transport of linear momentum.</li>
   
    <p>In the <code>faces</code> connector array, the transverse translational flow (<i>m</i>&Phi;dot) is only the
    force due to diffusion.  Translational advection is calculated from the velocity and the current.
    The thermal flow (<i>Q&#775;</i>) is only the rate of heat transfer due to diffusion.  The advection of
    thermal energy is determined from the thermodynamic state at the boundary and the current.</p>

    <p>In evaluating the dynamics of a phase, it is typically assumed that all of the species
    exist at the same velocity and temperature.  The translational and thermal time constants
    are usually much shorter than the time span of interest due to the very small coupling
    resistances.  If this is the case, connect the <code>direct</code>
    connectors of the species.  This will reduce the index of the problem.</p>

    <p>For the variables that relate to transport,
    the first index is the axis and the second index is the side.  The sides
    are ordered from negative to positive, according to the
    <a href=\"modelica://FCSys.Species.Enumerations.Side\">Side</a> enumeration.
    Velocity and force are additionally indexed by
    the orientation of the momentum with respect to the face.
    The orientations are ordered following the normal axis in Cartesian space,
    according to the
    <a href=\"modelica://FCSys.Species.Enumerations.Orient\">Orient</a> enumeration.</p>
    </html>"),
      Diagram(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          initialScale=0.1), graphics),
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
              100}}), graphics={Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={127,127,127},
            pattern=LinePattern.Dash,
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid), Text(
            extent={{-100,-20},{100,20}},
            textString="%name",
            lineColor={0,0,0},
            origin={-40,40},
            rotation=45)}));
  end Species;

public
  package Enumerations "Choices of options"
    extends Modelica.Icons.TypesPackage;

    type Axis = enumeration(
        x "X",
        y "Y",
        z "Z") "Enumeration for Cartesian axes";
    type Orient = enumeration(
        after "Axis following the normal axis in Cartesian coordinates",
        before "Axis preceding the normal axis in Cartesian coordinates")
      "Enumeration for orientations relative to a face" annotation (
        Documentation(info="
    <html><p><code>Orient.after</code> indicates the axis following the axis normal to the face
    in Cartesian coordinates (x, y, z).
    <code>Orient.before</code> indicates the axis preceding the normal axis
    in Cartesian coordinates (or following it twice).</p></html>"));
    type Side = enumeration(
        n "Negative",
        p "Positive (greater position along the Cartesian axis)")
      "Enumeration for sides of a region or subregion";
    type ConsThermo = enumeration(
        IC "Initial condition imposed forever (no conservation)",
        steady "Steady (conservation with steady state)",
        dynamic "Dynamic (conservation with storage)")
      "Options for the conservation of material or energy";
    type ConsMom = enumeration(
        steady "Steady (conservation with steady state)",
        dynamic "Dynamic (conservation with storage)")
      "Options for the conservation of momentum";
    type Init = enumeration(
        none "No initialization",
        amount "Prescribed amount",
        amountSS "Steady-state amount",
        density "Prescribed density",
        densitySS "Steady-state density",
        volume "Prescribed volume",
        volumeSS "Steady-state volume",
        pressure "Prescribed pressure",
        pressureSS "Steady-state pressure",
        temperature "Prescribed temperature",
        temperatureSS "Steady-state temperature",
        specificEnthalpy "Prescribed specific enthalpy",
        specificEnthalpySS "Steady-state specific enthalpy",
        Gibbs "Prescribed Gibbs potential",
        GibbsSS "Steady-state Gibbs potential")
      "Methods of initializing a thermodynamic quantity (material or energy)";

  end Enumerations;
  annotation (Documentation(info="
<html>
  <p><b>Licensed by the Georgia Tech Research Corporation under the Modelica License 2</b><br>
Copyright 2007&ndash;2013, <a href=\"http://www.gtrc.gatech.edu/\">Georgia Tech Research Corporation</a>.</p>

<p><i>This Modelica package is <u>free</u> software and the use is completely at <u>your own risk</u>;
it can be redistributed and/or modified under the terms of the Modelica License 2. For license conditions (including the
disclaimer of warranty) see <a href=\"modelica://FCSys.UsersGuide.License\">
FCSys.UsersGuide.License</a> or visit <a href=\"http://www.modelica.org/licenses/ModelicaLicense2\">
http://www.modelica.org/licenses/ModelicaLicense2</a>.</i></p>
</html>"));
end Species;

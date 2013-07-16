within FCSys;
package Species "Dynamic models of chemical species"
  extends Modelica.Icons.Package;
  package 'C+' "C"
    extends Modelica.Icons.Package;
    package Graphite "<html>C<sup>+</sup> graphite</html>"
      extends Modelica.Icons.Package;
      model Calibrated "Correlations with adjustment factors"
        extends SpeciesSolid(
          redeclare replaceable package Data = Characteristics.'C+'.Graphite,
          nu=k_nu*Data.nu(T, v),
          theta=k_theta*Data.theta(T, v));

        // Note:  In Dymola 7.4,
        // "redeclare replaceable package Data = FCSys.Characteristics.C.Graphite"
        // must be used instead of
        // "redeclare replaceable FCSys.Characteristics.C.Graphite Data" in
        // order for this model to pass its check.  This applies to the other
        // species models too.

        parameter Q.NumberAbsolute k_nu(final nominal=1) = 1
          "<html>Adjustment factor for thermal independity (<i>k</i><sub>&nu;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        parameter Q.NumberAbsolute k_theta(final nominal=1) = 1
          "<html>Adjustment factor for thermal resistivity (<i>k</i><sub>&theta;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="'C+'",
          Documentation(info="<html><p>Please see the documentation of the
    <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"),

          Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                  {100,100}}), graphics),
          Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
                  100,100}}), graphics={Text(
                extent={{-150,90},{-118,52}},
                lineColor={0,0,255},
                textString="%t.test")}));

      end Calibrated;

      model Correlated "Correlated properties"
        extends SpeciesSolid(redeclare replaceable package Data =
              Characteristics.'C+'.Graphite);
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="'C+'",
          Documentation(info=
                "<html><p>Please see the documentation of the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"),

          Icon(graphics));

      end Correlated;

      model Fixed "Fixed properties"

        extends SpeciesSolid(
          redeclare replaceable package Data = Characteristics.'C+'.Graphite (
              n_c=0,
              T_lim_c={0,Modelica.Constants.inf},
              b_c=[935*U.J*Data.m/(U.kg*U.K)],
              B_c=[Data.Deltah0_f - (935*U.J*Data.m/U.kg)*298.15, 154.663*U.J/(
                  U.mol*U.K) - (935*U.J*Data.m/(U.kg*U.K))*ln(298.15*U.K)]),
          redeclare final parameter Q.Mobility mu=0,
          redeclare parameter Q.TimeAbsolute nu=Data.nu(),
          redeclare parameter Q.ResistivityThermal theta=U.m*U.K/(11.1*U.W));

        // In Dymola 7.4, the specific heat capacity must be entered explicitly in
        // B_c (i.e., 935*U.J*Data.m/(U.kg*U.K) instead of Data.b_c[1, 1]).

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
    <li>Mobility is zero (by default).</li>
    </ol></p>

   <p>The default isobaric specific heat capacity (<code>b_c = [935*U.J*Data.m/(U.kg*U.K)]</code>)
   and thermal
   resistivity (<code>theta = U.m*U.K/(11.1*U.W)</code>) are for graphite fiber epoxy (25% vol)
   composite (with heat flow parallel to the fibers) at 300&nbsp;K
   [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, p. 909].
   The integration offset for specific entropy is set such that
   the specific entropy is 154.663&nbsp;J/(mol&middot;K) at 25&nbsp;&deg;C and <i>p</i><sup>o</sup> (1&nbsp;atm).
   This is the value from Table B in [<a href=\"modelica://FCSys.UsersGuide.References\">McBride2002</a>].
   Additional thermal data is listed in <a href=\"#Tab1\">Table 1</a>.</p>

  <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
    <caption align=\"top\" id=\"Tab1\">Table 1: Properties of forms of C [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, p. 909].</caption>
    <tr>
      <th rowspan=3 valign=\"middle\"><code>T<br>/U.K</code></th>
      <th rowspan=1 colspan=2 width=1 valign=\"middle\">Diamond (type IIa)</th>
      <th rowspan=1 colspan=1 width=1 valign=\"middle\">Amorphous<br>carbon</th>
      <th rowspan=1 colspan=3 width=1 valign=\"middle\">Graphite (pyrolytic)</th>
      <th rowspan=1 colspan=3 width=1>Graphite fiber epoxy (25% vol)<br>composite</th>
    </tr>
    <tr>
      <th rowspan=2 valign=\"middle\"><code>c_p*U.kg<br>*U.K<br>/(U.J<br>*Data.m)</code></th>
      <th rowspan=2 valign=\"middle\"><code>theta<br>*U.W<br>/(U.m<br>*U.K)</code></th>
      <th rowspan=2 valign=\"middle\"><code>theta<br>*U.W<br>/(U.m<br>*U.K)</code></th>
      <th rowspan=2 valign=\"middle\"><code>c_p*U.kg<br>*U.W<br>/(U.J<br>*Data.m)</code></th>
      <th rowspan=1 colspan=2><code>theta*U.W/(U.m*U.K)</code></th>
      <th rowspan=2 valign=\"middle\"><code>c_p*U.kg<br>*U.K<br>/(U.J<br>*Data.m)</code></th>
      <th rowspan=1 colspan=2><code>theta*U.W/(U.m*U.K)</code></th>
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

  <p>For more information, see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"));

      end Fixed;

    end Graphite;

  end 'C+';

  package 'C19HF37O5S-'
    "<html>C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S<sup>-</sup></html>"
    extends Modelica.Icons.Package;
    package Ionomer
      "<html>C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S ionomer</html>"
      extends Modelica.Icons.Package;
      model Calibrated "Correlations with adjustment factors"
        extends SpeciesSolid(
          redeclare replaceable package Data =
              Characteristics.'C19HF37O5S-'.Ionomer,
          nu=k_nu*Data.nu(T, v),
          theta=k_theta*Data.theta(T, v));

        parameter Q.NumberAbsolute k_nu(final nominal=1) = 1
          "<html>Adjustment factor for thermal independity (<i>k</i><sub>&nu;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        parameter Q.NumberAbsolute k_theta(final nominal=1) = 1
          "<html>Adjustment factor for thermal resistivity (<i>k</i><sub>&theta;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="'C19HF37O5S-'",
          Documentation(info=
                "<html><p>Please see the documentation of the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"),

          Diagram(graphics));

      end Calibrated;

      model Correlated "Correlated properties"
        extends SpeciesSolid(redeclare replaceable package Data =
              Characteristics.'C19HF37O5S-'.Ionomer);
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="'C19HF37O5S-'",
          Documentation(info=
                "<html><p>Please see the documentation of the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"));

      end Correlated;

      model Fixed "Fixed properties"
        extends SpeciesSolid(
          redeclare replaceable package Data =
              Characteristics.'C19HF37O5S-'.Ionomer,
          redeclare parameter Q.TimeAbsolute nu=Data.nu(),
          redeclare parameter Q.ResistivityThermal theta=U.m*U.K/(0.16*U.W));

        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="'C19HF37O5S-'",
          Documentation(info="<html><p>Assumptions:
    <ol>
    <li>The thermal independity and thermal resistivity are fixed (e.g., independent of temperature)</li>
    </ol></p>

    <p>The default thermal resistivity (<code>theta=U.m*U.K/(0.16*U.W)</code>) is of dry
  Nafion 115 [<a href=\"modelica://FCSys.UsersGuide.References\">Kandlikar2009</a>, p. 1277].</p>

<p>For more information, see the
    <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"));

      end Fixed;

    end Ionomer;

  end 'C19HF37O5S-';

  package 'e-' "<html>e<sup>-</sup></html>"
    extends Modelica.Icons.Package;
    package Graphite "<html>e<sup>-</sup> in graphite</html>"
      extends Modelica.Icons.Package;
      model Calibrated "Correlations with adjustment factors"
        extends SpeciesIsochoric(
          redeclare replaceable package Data = Characteristics.'e-'.Graphite,
          final tauprime=0,
          final mu=sigma*Data.v_Tp(T, p),
          nu=k_nu*Data.nu(T, v),
          beta=k_beta*Data.beta(T, v),
          zeta=k_zeta*Data.zeta(T, v),
          theta=Modelica.Constants.inf,
          initMaterial=InitScalar.Density);

        Q.ConductivityElectrical sigma=k_sigma*Data.mu(T, p)/Data.v_Tp(T, p)
          "<html>Electrical conductivity (&sigma;)</html>"
          annotation (Dialog(group="Material properties"));

        parameter Q.NumberAbsolute k_nu(final nominal=1) = 1
          "<html>Adjustment factor for thermal independity (<i>k</i><sub>&nu;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        parameter Q.NumberAbsolute k_eta(final nominal=1) = 1
          "<html>Adjustment factor for material resistivity (<i>k</i><sub>&eta;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        parameter Q.NumberAbsolute k_beta(final nominal=1) = 1
          "<html>Adjustment factor for dynamic compressibility (<i>k</i><sub>&beta;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        parameter Q.NumberAbsolute k_zeta(final nominal=1) = 1
          "<html>Adjustment factor for fluidity (<i>k</i><sub>&zeta;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        parameter Q.NumberAbsolute k_sigma(final nominal=1) = 1
          "<html>Adjustment factor for electrical conductivity (<i>k</i><sub>&sigma;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="'e-'",
          Documentation(info="<html>    <p>Assumptions:<ol>
    <li>The density is equal to that of C<sup>+</sup> as graphite.</li>
          <li>The thermal resistivity is infinite.  All of the thermal conductance is attributed to the substrate
          (e.g., <a href=\"modelica://FCSys.Species.'C+'.Graphite\">C+</a>).<li>
          <li>The phase change interval (&tau;&prime;) is zero.  The rate of phase change is
          governed by other configurations.</li>
          <li>The conductivity is mapped to the mobility of the electrons by assuming that
          the mobility of the substrate (e.g., 
          <a href=\"modelica://FCSys.Species.'C+'.Graphite\">C+</a>) is zero.</li>
    </ol></p>

    <p>For more information, see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"));

      end Calibrated;

      model Correlated "Correlated properties"
        extends SpeciesIsochoric(
          redeclare replaceable package Data = Characteristics.'e-'.Graphite,
          final tauprime=0,
          final mu=sigma*Data.v_Tp(T, p),
          initMaterial=InitScalar.Density);

        Q.ConductivityElectrical sigma=Data.mu(T, p)/Data.v_Tp(T, p)
          "<html>Electrical conductivity (&sigma;)</html>"
          annotation (Dialog(group="Material properties"));

        // **Move assumption to Characteristics.'e-'.Graphite, do the same for other e-.
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="'e-'",
          Documentation(info="<html>    <p>Assumptions:<ol>
    <li>The density is equal to that of C<sup>+</sup> as graphite.</li>
          <li>The thermal resistivity is infinite.  All of the thermal conductance is attributed to the substrate
          (e.g., <a href=\"modelica://FCSys.Species.'C+'.Graphite\">C+</a>).<li>
          <li>The phase change interval (&tau;&prime;) is zero.  The rate of phase change is
          governed by other configurations.</li>
          <li>The conductivity is mapped to the mobility of the electrons by assuming that
          the mobility of the substrate (e.g., 
          <a href=\"modelica://FCSys.Species.'C+'.Graphite\">C+</a>) is zero.</li>
    </ol></p>

    <p>For more information, see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"));

      end Correlated;

      model Fixed "Fixed properties"
        extends Species(
          redeclare replaceable package Data = Characteristics.'e-'.Gas (n_v={0,
                  0}, b_v=FCSys.Characteristics.'C+'.Graphite.b_v),
          final tauprime=0,
          final mu=sigma*v,
          redeclare parameter Q.TimeAbsolute nu=Data.nu(),
          redeclare parameter Q.Mobility eta=Data.eta(),
          redeclare parameter Q.Fluidity beta=Data.beta(),
          redeclare parameter Q.Fluidity zeta=Data.zeta(),
          final theta=Modelica.Constants.inf);
        //invertEOS=false
        //initMaterial=InitScalar.Volume
        //
        //consMaterial=Conservation.IC,
        //            redeclare parameter Q.Fluidity beta=1e-5*Data.beta(),
        //redeclare final parameter Q.Fluidity beta=0,

        parameter Q.ConductivityElectrical sigma=1e-10*Data.mu()/Data.v_Tp()
          "<html>Electrical conductivity (&sigma;)</html>"
          annotation (Dialog(group="Material properties"));
        // **temp factor on sigma
        //    redeclare final parameter Q.Fluidity beta=Data.beta(),
        // **recreate Graphite e- Characteristics model, use it here
        //    redeclare parameter Q.Mobility eta=Data.eta(),
        //    redeclare parameter Q.Fluidity beta=1e-5*Data.beta(),
        // add eta to other e- and H+ models.
        annotation (
          group="Material properties",
          defaultComponentPrefixes="replaceable",
          defaultComponentName="'e-'",
          Documentation(info="<html>

    <p>Assumptions:<ol>
    <li>The density is equal to that of C<sup>+</sup> as graphite.</li>
          <li>The thermal resistivity is infinite.  All of the thermal conductance is attributed to 
          the substrate
          (e.g., <a href=\"modelica://FCSys.Species.'C+'.Graphite\">C+</a>).<li>
          <li>The phase change interval (&tau;&prime;) is zero.  The rate of phase change is
          governed by other configurations.</li>
          <li>The conductivity is mapped to the mobility of the electrons by assuming that
          the mobility of the substrate (e.g., 
          <a href=\"modelica://FCSys.Species.'C+'.Graphite\">C+</a>) is zero.</li>
    </ol></p>

    <p>For more information, see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"),

          Diagram(graphics));

      end Fixed;

    end Graphite;

  end 'e-';

  package 'H+' "<html>H<sup>+</sup></html>"
    extends Modelica.Icons.Package;
    package Ionomer "<html>H<sup>+</sup> in ionomer</html>"
      extends Modelica.Icons.Package;
      model Calibrated "Correlations with adjustment factors"
        extends SpeciesIsochoric(
          redeclare replaceable package Data = Characteristics.'H+'.Ionomer,
          initMaterial=InitScalar.Density,
          final tauprime=0,
          mu=k_mu*Data.mu(T, v),
          nu=k_nu*Data.nu(T, v),
          beta=k_beta*Data.beta(T, v),
          zeta=k_zeta*Data.zeta(T, v),
          theta=k_theta*Data.theta(T, v));

        parameter Q.NumberAbsolute k_mu(final nominal=1) = 1
          "<html>Adjustment factor for mobility (<i>k</i><sub>&mu;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        parameter Q.NumberAbsolute k_nu(final nominal=1) = 1
          "<html>Adjustment factor for thermal independity (<i>k</i><sub>&nu;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        parameter Q.NumberAbsolute k_beta(final nominal=1) = 1
          "<html>Adjustment factor for dynamic compressibility (<i>k</i><sub>&beta;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        parameter Q.NumberAbsolute k_zeta(final nominal=1) = 1
          "<html>Adjustment factor for fluidity (<i>k</i><sub>&zeta;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        parameter Q.NumberAbsolute k_theta(final nominal=1) = 1
          "<html>Adjustment factor for thermal resistivity (<i>k</i><sub>&theta;</i></sub>)</html>"
          annotation (Dialog(group="Material properties"));
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="'H+'",
          Documentation(info="<html><p>Assumptions:<ol>
    <li>The density of H<sup>+</sup> is equal to that of
  C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S<sup>-</sup> or approximately 1.912 M.  Note that
  this is greater than that measured by Spry and Fayer (see
  <a href=\"modelica://FCSys.Characteristics.'H+'.Ionomer\">Characteristics.'H+'.Ionomer</a>), but it
  simplifies the model by requiring only C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S<sup>-</sup>
  (not C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S) for charge neutrality.</li>
          <li>The phase change interval (&tau;&prime;) is zero.  The rate of phase change is
          governed by other configurations.</li>
    </ol></p>

  <p>For more information, see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"));

      end Calibrated;

      model Correlated "Correlated properties"
        extends SpeciesIsochoric(
          redeclare replaceable package Data = Characteristics.'H+'.Ionomer,
          final tauprime=0,
          initMaterial=InitScalar.Density);

        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="'H+'",
          Documentation(info="<html><p>Assumptions:<ol>
    <li>The density of H<sup>+</sup> is equal to that of
  C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S<sup>-</sup> or approximately 1.912 M.  Note that
  this is greater than that measured by Spry and Fayer (see
  <a href=\"modelica://FCSys.Characteristics.'H+'.Ionomer\">Characteristics.'H+'.Ionomer</a>), but it
  simplifies the model by requiring only C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S<sup>-</sup>
  (not C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S) for charge neutrality.</li>
          <li>The phase change interval (&tau;&prime;) is zero.  The rate of phase change is governed by other configurations.</li>
    </ol></p>

    <p>For more information, see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"));

      end Correlated;

      model Fixed "Fixed properties"
        extends Species(
          redeclare replaceable package Data = Characteristics.'H+'.Gas,
          final tauprime=0,
          redeclare parameter Q.Mobility mu=Data.mu(),
          redeclare parameter Q.TimeAbsolute nu=Data.nu(),
          redeclare parameter Q.Mobility eta=Data.eta(),
          redeclare parameter Q.Fluidity beta=Data.beta(),
          redeclare parameter Q.Fluidity zeta=1/(5.3e-6*U.Pa*U.s),
          redeclare parameter Q.ResistivityThermal theta=U.m*U.K/(0.1661*U.W));
        //consMaterial=Conservation.IC,
        //invertEOS=false
        //initMaterial=InitScalar.Volume
        //(n_v={0,0},          b_v=FCSys.Characteristics.'C19HF37O5S-'.Ionomer.b_v),

        //        final theta=Modelica.Constants.inf,

        // **recreate Ionomer H+ Characteristics model, use it here

        // **set rho_IC, do the same for Correlated and Calibrated.

        // See the documentation for a table of values.
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="'H+'",
          Documentation(info="<html><p>The initial density corresponds to the measurement
  by Spry and Fayer (0.95 M) in Nafion<sup>&reg;</sup> at
  &lambda; = 12, where &lambda; is the number of
  H<sub>2</sub>O molecules to SO<sub>3</sub>H
  endgroups.  At &lambda; = 22, the concentration was measured at 0.54 M
  [<a href=\"modelica://FCSys.UsersGuide.References\">Spry2009</a>]. **Implement this, copy to Correlated and Calibrated</p>

<p>Assumptions:<ol>
    <li>The generalized resistivities (&beta;, &zeta;, &theta;) are fixed (e.g., independent of temperature).</li>
    <li>The density of H<sup>+</sup> is equal to that of
  C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S<sup>-</sup> or approximately 1.912 M.  Note that
  this is greater than that measured by Spry and Fayer (see
  <a href=\"modelica://FCSys.Characteristics.'H+'.Ionomer\">Characteristics.'H+'.Ionomer</a>), but it
  simplifies the model by requiring only C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S<sup>-</sup>
  (not C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S) for charge neutrality.</li>
          <li>The phase change interval (&tau;&prime;) is zero.  The rate of phase change is governed by other configurations.</li>
    </ol></p>

  <p>The default resistivities (<code>zeta=1/(5.3e-6*U.Pa*U.s)</code> and <code>theta=U.m*U.K/(0.1661*U.W)</code>) are of H gas
  (rather than H<sup>+</sup>) at 300&nbsp;K from [<a href=\"modelica://FCSys.UsersGuide.References\">Schetz1996</a>, p. 139].
  <a href=\"#Tab1\">Table 1</a> lists the properties at other temperatures.</p>

    <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
  <caption align=\"top\" id=\"Tab1\">Table 1: Properties of H gas (not H<sup>+</sup>) [<a href=\"modelica://FCSys.UsersGuide.References\">Schetz1996</a>, p. 139]</caption>
<tr>
      <th valign=\"middle\"><code>T<br>/U.K</code></th>
      <th width=1><code>zeta<br>*U.Pa*U.s</code></th>
      <th width=1><code>theta*U.W<br>/(U.m*U.K)</code></th>
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

<p>For more information, see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"));

      end Fixed;

    end Ionomer;

  end 'H+';

  package H2 "<html>H<sub>2</sub></html>"
    extends Modelica.Icons.Package;
    package Gas "<html>H<sub>2</sub> gas</html>"
      extends Modelica.Icons.Package;
      model Calibrated "Correlations with adjustment factors"
        extends Species(
          redeclare replaceable package Data = FCSys.Characteristics.H2.Gas (
                b_v=[1],n_v={-1,0}),
          tauprime=k_tauprime*Data.tauprime(T, v),
          mu=k_mu*Data.mu(T, v),
          nu=k_nu*Data.nu(T, v),
          eta=k_eta*Data.eta(T, v),
          beta=k_beta*Data.beta(T, v),
          zeta=k_zeta*Data.zeta(T, v),
          theta=k_theta*Data.theta(T, v));
        parameter Q.NumberAbsolute k_tauprime(final nominal=1) = 1
          "<html>Adjustment factor for the phase change interval (<i>k</i><sub>&tau;&prime;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        parameter Q.NumberAbsolute k_mu(final nominal=1) = 1
          "<html>Adjustment factor for mobility (<i>k</i><sub>&mu;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        parameter Q.NumberAbsolute k_nu(final nominal=1) = 1
          "<html>Adjustment factor for thermal independity (<i>k</i><sub>&nu;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        parameter Q.NumberAbsolute k_eta(final nominal=1) = 1
          "<html>Adjustment factor for material resistivity (<i>k</i><sub>&eta;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        parameter Q.NumberAbsolute k_beta(final nominal=1) = 1
          "<html>Adjustment factor for dynamic compressibility (<i>k</i><sub>&beta;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        parameter Q.NumberAbsolute k_zeta(final nominal=1) = 1
          "<html>Adjustment factor for fluidity (<i>k</i><sub>&zeta;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        parameter Q.NumberAbsolute k_theta(final nominal=1) = 1
          "<html>Adjustment factor for thermal resistivity (<i>k</i><sub>&theta;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="H2",
          Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol></p>

<p>For more information, see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"));

      end Calibrated;

      model Correlated "Correlated properties"
        extends Species(redeclare replaceable package Data =
              Characteristics.H2.Gas (b_v=[1], n_v={-1,0}));
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="H2",
          Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol></p>

<p>For more information, see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"));

      end Correlated;

      model Fixed "Fixed properties"
        extends Species(
          redeclare replaceable package Data = FCSys.Characteristics.H2.Gas (
                b_v=[1],n_v={-1,0}),
          redeclare parameter Q.TimeAbsolute tauprime=Data.tauprime(),
          redeclare parameter Q.Mobility mu=Data.mu(),
          redeclare parameter Q.TimeAbsolute nu=Data.nu(),
          redeclare parameter Q.ResistivityMaterial eta=Data.eta(),
          redeclare parameter Q.Fluidity beta=Data.beta(),
          redeclare parameter Q.Fluidity zeta=1/(8.96e-6*U.Pa*U.s),
          redeclare parameter Q.ResistivityThermal theta=U.m*U.K/(0.183*U.W));
        // See the documentation for a table of values.
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="H2",
          Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
        <li>The generalized resistivities (&beta;, &zeta;, &theta;) are fixed (e.g., independent of temperature).</li>

    </ol></p>

<p>Additional notes:<ul>
<li>
  The specific heat capacity is not fixed because it would
  affect the chemical potential and result in an incorrect cell
  potential.</li></ul></p>

<p>The default resistivities (<code>zeta=1/(89.6e-7*U.Pa*U.s)</code>
and <code>theta=U.m*U.K/(183e-3*U.W)</code>) are based on data of H<sub>2</sub> gas at 1&nbsp;atm and
  300&nbsp;K from Incropera and DeWitt [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 919&ndash;920].
  <a href=\"#Tab1\">Table 1</a> lists the properties at other temperatures.</p>

  <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
    <caption align=\"top\" id=\"Tab1\">Table 1: Properties of H<sub>2</sub> gas at 1&nbsp;atm [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 919&ndash;920].</caption>
  <tr>
      <th valign=\"middle\"><code>T<br>/U.K</code></th>
      <th width=1><code>c_p*U.kg*U.K<br>/(U.J*Data.m)</code></th>
      <th width=1 ><code>zeta<br>*U.Pa*U.s</code></th>
      <th width=1 ><code>theta*U.W<br>/(U.m*U.K)</code></th>
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
<p>For more information, see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"),

          Diagram(graphics));

      end Fixed;

    end Gas;

  end H2;

  package H2O "<html>H<sub>2</sub>O</html>"
    extends Modelica.Icons.Package;
    package Gas "<html>H<sub>2</sub>O gas</html>"
      extends Modelica.Icons.Package;
      model Calibrated "Correlations with adjustment factors"
        extends Species(
          redeclare replaceable package Data = FCSys.Characteristics.H2O.Gas (
                b_v=[1],n_v={-1,0}),
          tauprime=k_tauprime*Data.tauprime(T, v),
          mu=k_mu*Data.mu(T, v),
          nu=k_nu*Data.nu(T, v),
          eta=k_eta*Data.eta(T, v),
          beta=k_beta*Data.beta(T, v),
          zeta=k_zeta*Data.zeta(T, v),
          theta=k_theta*Data.theta(T, v));
        parameter Q.NumberAbsolute k_tauprime(final nominal=1) = 1
          "<html>Adjustment factor for the phase change interval (<i>k</i><sub>&tau;&prime;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        parameter Q.NumberAbsolute k_mu(final nominal=1) = 1
          "<html>Adjustment factor for mobility (<i>k</i><sub>&mu;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        parameter Q.NumberAbsolute k_nu(final nominal=1) = 1
          "<html>Adjustment factor for thermal independity (<i>k</i><sub>&nu;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        parameter Q.NumberAbsolute k_eta(final nominal=1) = 1
          "<html>Adjustment factor for material resistivity (<i>k</i><sub>&eta;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        parameter Q.NumberAbsolute k_beta(final nominal=1) = 1
          "<html>Adjustment factor for dynamic compressibility (<i>k</i><sub>&beta;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        parameter Q.NumberAbsolute k_zeta(final nominal=1) = 1
          "<html>Adjustment factor for fluidity (<i>k</i><sub>&zeta;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        parameter Q.NumberAbsolute k_theta(final nominal=1) = 1
          "<html>Adjustment factor for thermal resistivity (<i>k</i><sub>&theta;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="H2O",
          Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol></p>

<p>For more information, see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"));

      end Calibrated;

      model Correlated "Correlated properties"
        extends Species(redeclare replaceable package Data =
              Characteristics.H2O.Gas (b_v=[1], n_v={-1,0}));
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="H2O",
          Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol></p>

<p>For more information, see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"));

      end Correlated;

      model Fixed "Fixed properties"
        extends Species(
          redeclare replaceable package Data = FCSys.Characteristics.H2O.Gas (
                b_v=[1],n_v={-1,0}),
          redeclare parameter Q.TimeAbsolute tauprime=1e6*Data.tauprime(),
          redeclare parameter Q.Mobility mu=Data.mu(),
          redeclare parameter Q.TimeAbsolute nu=Data.nu(),
          redeclare parameter Q.ResistivityMaterial eta=Data.eta(),
          redeclare parameter Q.Fluidity beta=Data.beta(),
          redeclare parameter Q.Fluidity zeta=1/(9.09e-6*U.Pa*U.s),
          redeclare parameter Q.ResistivityThermal theta=U.m*U.K/(19.6e-3*U.W));
        // **temp factor on tauprime
        // See the documentation for tables of values.
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="H2O",
          Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
        <li>The generalized resistivities (&beta;, &zeta;, &theta;) are fixed (e.g., independent of temperature).</li>
    </ol></p>

  <p>Notes:<ul>
  <li>The specific heat capacity is not fixed because it would
  affect the chemical potential and result in an incorrect cell
  potential.</li></ul></p>

<p>The default resistivities (<code>zeta=1/(9.09e-6*U.Pa*U.s)</code>
and <code>theta=U.m*U.K/(19.6e-3*U.W)</code>) are of H<sub>2</sub>O gas at saturation pressure and
  300&nbsp;K from Incropera and DeWitt [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, p. 921].  <a href=\"#Tab1\">Table 1</a> lists the properties at
  saturation pressure and other temperatures.  <a href=\"#Tab2\">Table 2</a> lists the properties of H<sub>2</sub>O gas at 1&nbsp;atm.
  See also
  <a href=\"http://www.engineeringtoolbox.com/water-dynamic-kinematic-viscosity-d_596.html\">http://www.engineeringtoolbox.com/water-dynamic-kinematic-viscosity-d_596.html</a>.</p>

  <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
    <caption align=\"top\" id=\"Tab1\">Table 1: Properties of H<sub>2</sub>O gas at saturation pressure [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 924&ndash;925].</caption>
  <tr>
      <th valign=\"middle\"><code>T<br>/U.K</code></th>
      <th width=1><code>c_p*U.kg*U.K<br>/(U.J*Data.m)</code></th>
      <th width=1><code>zeta<br>*U.Pa*U.s</code></th>
      <th width=1><code>theta*U.W<br>/(U.m*U.K)</code></th>
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
      <th valign=\"middle\"><code>T<br>/U.K</code></th>
      <th width=1><code>c_p*U.kg*U.K<br>/(U.J*Data.m)</code></th>
      <th width=1 ><code>zeta<br>*U.Pa*U.s</code></th>
      <th width=1 ><code>theta*U.W<br>/(U.m*U.K)</code></th>
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

<p>For more information, see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"));

      end Fixed;

    end Gas;

    package Ionomer "<html>H<sub>2</sub>O in ionomer</html>"
      extends Modelica.Icons.Package;
      model Calibrated "Correlations with adjustment factors"
        extends Species(
          redeclare replaceable package Data = Characteristics.H2O.Ionomer,
          final tauprime=0,
          mu=k_mu*Data.mu(T, v),
          nu=k_nu*Data.nu(T, v),
          eta=k_eta*Data.eta(T, v),
          beta=k_beta*Data.beta(T, v),
          zeta=k_zeta*Data.zeta(T, v),
          theta=k_theta*Data.theta(T, v));

        parameter Q.NumberAbsolute k_mu(final nominal=1) = 1
          "<html>Adjustment factor for mobility (<i>k</i><sub>&mu;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        parameter Q.NumberAbsolute k_nu(final nominal=1) = 1
          "<html>Adjustment factor for thermal independity (<i>k</i><sub>&nu;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        parameter Q.NumberAbsolute k_eta(final nominal=1) = 1
          "<html>Adjustment factor for material resistivity (<i>k</i><sub>&eta;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        parameter Q.NumberAbsolute k_beta(final nominal=1) = 1
          "<html>Adjustment factor for dynamic compressibility (<i>k</i><sub>&beta;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        parameter Q.NumberAbsolute k_zeta(final nominal=1) = 1
          "<html>Adjustment factor for fluidity (<i>k</i><sub>&zeta;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        parameter Q.NumberAbsolute k_theta(final nominal=1) = 1
          "<html>Adjustment factor for thermal resistivity (<i>k</i><sub>&theta;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="H2O",
          Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
        <li>The phase change interval (&tau;&prime;) is zero.  The rate of phase change is governed by the
        other configurations (e.g., gas).</li>
          </ol></p>

<p>For more information, see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"));

      end Calibrated;

      model Correlated "Correlated properties"
        extends Species(redeclare replaceable package Data =
              Characteristics.H2O.Ionomer (b_v=[1], n_v={-1,0}), final tauprime
            =0);
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="H2O",
          Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
    <li>The phase change interval (&tau;&prime;) is zero.  The rate of phase change is governed by the
   other configurations (e.g., gas).</li>
          </ol></p>

<p>For more information, see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"));

      end Correlated;

      model Fixed "Fixed properties"
        extends Species(
          redeclare replaceable package Data = Characteristics.H2O.Ionomer (b_v
                =[1], n_v={-1,0}),
          final tauprime=0,
          redeclare parameter Q.Mobility mu=Data.mu(),
          redeclare parameter Q.TimeAbsolute nu=Data.nu(),
          redeclare parameter Q.ResistivityMaterial eta=Data.eta(),
          redeclare parameter Q.Fluidity beta=Data.beta(),
          redeclare parameter Q.Fluidity zeta=Data.zeta(),
          redeclare parameter Q.ResistivityThermal theta=Data.theta());
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="H2O",
          Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
        <li>The generalized resistivities (&beta;, &zeta;, &theta;) are fixed (e.g., independent of temperature).</li>
        <li>The phase change interval (&tau;&prime;) is zero.  The rate of phase change is governed by the
        other configurations (e.g., gas).</li>
    </ol></p></p>

<p>For more information, see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"));

      end Fixed;

    end Ionomer;

    package Liquid "<html>H<sub>2</sub>O liquid</html>"
      extends Modelica.Icons.Package;
      model Calibrated "Correlations with adjustment factors"
        extends SpeciesIsochoric(
          redeclare replaceable package Data = Characteristics.H2O.Liquid,
          final tauprime=0,
          mu=k_mu*Data.mu(T, v),
          nu=k_nu*Data.nu(T, v),
          beta=k_beta*Data.beta(T, v),
          zeta=k_zeta*Data.zeta(T, v),
          theta=k_theta*Data.theta(T, v));

        parameter Q.NumberAbsolute k_mu(final nominal=1) = 1
          "<html>Adjustment factor for mobility (<i>k</i><sub>&mu;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        parameter Q.NumberAbsolute k_nu(final nominal=1) = 1
          "<html>Adjustment factor for thermal independity (<i>k</i><sub>&nu;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        parameter Q.NumberAbsolute k_beta(final nominal=1) = 1
          "<html>Adjustment factor for dynamic compressibility (<i>k</i><sub>&beta;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        parameter Q.NumberAbsolute k_zeta(final nominal=1) = 1
          "<html>Adjustment factor for fluidity (<i>k</i><sub>&zeta;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        parameter Q.NumberAbsolute k_theta(final nominal=1) = 1
          "<html>Adjustment factor for thermal resistivity (<i>k</i><sub>&theta;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="H2O",
          Documentation(info="<html><p>Assumptions:<ol>
        <li>The phase change interval (&tau;&prime;) is zero.  The rate of phase change is governed by the
        other configurations (e.g., gas).</li>
    </ol></p>

<p>For more information, see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"));

      end Calibrated;

      model Correlated "Correlated properties"
        extends SpeciesIsochoric(redeclare replaceable package Data =
              Characteristics.H2O.Liquid, final tauprime=0);
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="H2O",
          Documentation(info="<html><p>Assumptions:<ol>
        <li>The phase change interval (&tau;&prime;) is zero.  The rate of phase change is governed by the
        other configurations (e.g., gas).</li>
    </ol></p>

<p>For more information, see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"));

      end Correlated;

      model Fixed "Fixed properties"
        extends SpeciesIsochoric(
          redeclare replaceable package Data = Characteristics.H2O.Liquid,
          final tauprime=0,
          redeclare parameter Q.Mobility mu=Data.mu(),
          redeclare parameter Q.TimeAbsolute nu=Data.nu(),
          redeclare parameter Q.Fluidity beta=Data.beta(),
          redeclare parameter Q.Fluidity zeta=1/(855e-6*U.Pa*U.s),
          redeclare parameter Q.ResistivityThermal theta=U.m*U.K/(0.613*U.W));

        // See the documentation for tables of values.
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="H2O",
          Documentation(info="<html><p>Assumptions:<ol>
        <li>The generalized resistivities (&beta;, &zeta;, &theta;) are fixed (e.g., independent of temperature).</li>
        <li>The phase change interval (&tau;&prime;) is zero.  The rate of phase change is governed by the
        other configurations (e.g., gas).</li>
    </ol></p>

          <p>Notes:<ul>
  <li>
  The specific heat capacity is not fixed because it would
  affect the chemical potential and result in an incorrect cell
  potential.</li></ul></p>

<p>The default resistivities (<code>zeta=1/(855e-6*U.Pa*U.s)</code>
and <code>theta=U.m*U.K/(613e-3*U.W)</code>) are of H<sub>2</sub>O liquid at saturation pressure and
  300&nbsp;K from Incropera and DeWitt [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, p. 921].  <a href=\"#Tab1\">Table 1</a> lists the properties at
  saturation pressure and other temperatures.
  See also
  <a href=\"http://www.engineeringtoolbox.com/water-dynamic-kinematic-viscosity-d_596.html\">http://www.engineeringtoolbox.com/water-dynamic-kinematic-viscosity-d_596.html</a>.</p>

  <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
    <caption align=\"top\" id=\"Tab1\">Table 1: Properties of H<sub>2</sub>O liquid at saturation pressure [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 924&ndash;925].</caption>
  <tr>
      <th valign=\"middle\"><code>T<br>/U.K</code></th>
      <th width=1><code>c_p*U.kg*U.K<br>/(U.J*Data.m)</code></th>
      <th width=1><code>zeta<br>*U.Pa*U.s</code></th>
      <th width=1><code>theta*U.W<br>/(U.m*U.K)</code></th>
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

  <p>For more information, see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"));

      end Fixed;

    end Liquid;

  end H2O;

  package N2 "<html>N<sub>2</sub></html>"
    extends Modelica.Icons.Package;
    package Gas "<html>N<sub>2</sub> gas</html>"
      extends Modelica.Icons.Package;
      model Calibrated "Correlations with adjustment factors"
        extends Species(
          redeclare replaceable package Data = FCSys.Characteristics.N2.Gas (
                b_v=[1],n_v={-1,0}),
          tauprime=k_tauprime*Data.tauprime(T, v),
          mu=k_mu*Data.mu(T, v),
          nu=k_nu*Data.nu(T, v),
          eta=k_eta*Data.eta(T, v),
          beta=k_beta*Data.beta(T, v),
          zeta=k_zeta*Data.zeta(T, v),
          theta=k_theta*Data.theta(T, v));

        parameter Q.NumberAbsolute k_tauprime(final nominal=1) = 1
          "<html>Adjustment factor for the phase change interval (<i>k</i><sub>&tau;&prime;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        parameter Q.NumberAbsolute k_mu(final nominal=1) = 1
          "<html>Adjustment factor for mobility (<i>k</i><sub>&mu;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        parameter Q.NumberAbsolute k_nu(final nominal=1) = 1
          "<html>Adjustment factor for thermal independity (<i>k</i><sub>&nu;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        parameter Q.NumberAbsolute k_eta(final nominal=1) = 1
          "<html>Adjustment factor for material resistivity (<i>k</i><sub>&eta;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        parameter Q.NumberAbsolute k_beta(final nominal=1) = 1
          "<html>Adjustment factor for dynamic compressibility (<i>k</i><sub>&beta;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        parameter Q.NumberAbsolute k_zeta(final nominal=1) = 1
          "<html>Adjustment factor for fluidity (<i>k</i><sub>&zeta;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        parameter Q.NumberAbsolute k_theta(final nominal=1) = 1
          "<html>Adjustment factor for thermal resistivity (<i>k</i><sub>&theta;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="N2",
          Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol></p>

<p>For more information, see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"));

      end Calibrated;

      model Correlated "Correlated properties"
        extends Species(redeclare replaceable package Data =
              Characteristics.N2.Gas (b_v=[1], n_v={-1,0}));
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="N2",
          Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol></p>

<p>For more information, see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"));

      end Correlated;

      model Fixed "Fixed properties"
        import FCSys.BaseClasses.Utilities.Polynomial;

        extends Species(
          redeclare replaceable package Data = FCSys.Characteristics.N2.Gas (
              b_v=[1],
              n_v={-1,0},
              n_c=0,
              T_lim_c={0,Modelica.Constants.inf},
              b_c=[1041*U.J*Data.m/(U.kg*U.K)],
              B_c=[Data.Deltah0_f - (1041*U.J*Data.m/U.kg)*298.15, 191.610*U.J/
                  (U.mol*U.K) - (1041*U.J*Data.m/(U.kg*U.K))*ln(298.15*U.K)]),
          redeclare parameter Q.TimeAbsolute tauprime=Data.tauprime(),
          redeclare parameter Q.Mobility mu=Data.mu(),
          redeclare parameter Q.TimeAbsolute nu=Data.nu(),
          redeclare parameter Q.Fluidity beta=Data.beta(),
          redeclare parameter Q.ResistivityMaterial eta=Data.eta(),
          redeclare parameter Q.Fluidity zeta=1/(17.82e-6*U.Pa*U.s),
          redeclare parameter Q.ResistivityThermal theta=U.m*U.K/(25.9e-3*U.W));

        // **    redeclare final parameter Q.Fluidity beta=0,

        // **temp factor beta for PipeFlow

        // See the documentation for a table of values.
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="N2",
          Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
    <li>Fixed specific heat capacity (independent of temperature)</ol>
        <li>The generalized resistivities (&beta;, &zeta;, &theta;) are fixed (e.g., independent of temperature).</li>

    </ol></p>

<p>The default specific heat capacity (<code>b_c=[1.041e3*U.J*Data.m/(U.kg*U.K)]</code>) and resistivities
(<code>zeta=1/(17.82e-6*U.Pa*U.s)</code> and <code>theta=U.m*U.K/(25.9e-3*U.W))</code>) are based on data of gas at 1&nbsp;atm and
  300&nbsp;K from Incropera and DeWitt [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, p. 920].
   The integration offset for specific entropy is set such that
   the specific entropy is 191.610&nbsp;J/(mol&middot;K) at 25&nbsp;&deg;C and <i>p</i><sup>o</sup> (1&nbsp;bar).
   This is the value from Table B in [<a href=\"modelica://FCSys.UsersGuide.References\">McBride2002</a>].
   Additional thermal data is listed in <a href=\"#Tab1\">Table 1</a>.  <a href=\"#Tab2\">Table 2</a> lists
  values of the material resistivity or self diffusion coefficient.</p>

  <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
  <caption align=\"top\" id=\"Tab1\">Table 1: Properties of N<sub>2</sub> gas at 1&nbsp;atm [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, p. 920]</caption>
  <tr>
      <th valign=\"middle\"><code>T<br>/U.K</code></th>
      <th width=1><code>c_p*U.kg*U.K<br>/(U.J*Data.m)</code></th>
      <th width=1><code>zeta<br>*U.Pa*U.s</code></th>
      <th width=1><code>theta*U.W<br>/(U.m*U.K)</code></th>
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

<br>

  <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
  <caption align=\"top\" id=\"Tab2\">Table 2: Material resistivity of N<sub>2</sub> gas at 1&nbsp;atm
  [<a href=\"modelica://FCSys.UsersGuide.References\">Present1958</a>, p. 263]</caption>
  <tr>
      <th valign=\"middle\"><code>T<br>/U.K</code></th>
      <th width=1><code>eta*U.s<br>/U.cm^2</code></th>
    </tr>
<tr><td>77.7</td><td>0.0168</td></tr>
<tr><td>194.7</td><td>0.104</td></tr>
<tr><td>273.2</td><td>0.185</td></tr>
<tr><td>353.2</td><td>0.287</td></tr>
  </table></p>

  <p>The fluidity of air at 15.0&nbsp;&deg;C and 1&nbsp;atm is given by
       <code>zeta=1/(17.8e-6*U.Pa*U.s)</code>
   (<a href=\"http://en.wikipedia.org/wiki/Viscosity\">http://en.wikipedia.org/wiki/Viscosity</a>).</p>

<p>For more information, see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"));

      end Fixed;

    end Gas;

  end N2;

  package O2 "<html>O<sub>2</sub></html>"
    extends Modelica.Icons.Package;
    package Gas "<html>O<sub>2</sub> gas</html>"
      extends Modelica.Icons.Package;
      model Calibrated "Correlations with adjustment factors"
        extends Species(
          redeclare replaceable package Data = FCSys.Characteristics.O2.Gas (
                b_v=[1],n_v={-1,0}),
          tauprime=k_tauprime*Data.tauprime(T, v),
          mu=k_mu*Data.mu(T, v),
          nu=k_nu*Data.nu(T, v),
          eta=k_eta*Data.eta(T, v),
          beta=k_beta*Data.beta(T, v),
          zeta=k_zeta*Data.zeta(T, v),
          theta=k_theta*Data.theta(T, v));

        parameter Q.NumberAbsolute k_tauprime(final nominal=1) = 1
          "<html>Adjustment factor for the phase change interval (<i>k</i><sub>&tau;&prime;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        parameter Q.NumberAbsolute k_mu(final nominal=1) = 1
          "<html>Adjustment factor for mobility (<i>k</i><sub>&mu;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        parameter Q.NumberAbsolute k_nu(final nominal=1) = 1
          "<html>Adjustment factor for thermal independity (<i>k</i><sub>&nu;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        parameter Q.NumberAbsolute k_eta(final nominal=1) = 1
          "<html>Adjustment factor for material resistivity (<i>k</i><sub>&eta;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        parameter Q.NumberAbsolute k_beta(final nominal=1) = 1
          "<html>Adjustment factor for dynamic compressibility (<i>k</i><sub>&beta;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        parameter Q.NumberAbsolute k_zeta(final nominal=1) = 1
          "<html>Adjustment factor for fluidity (<i>k</i><sub>&zeta;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        parameter Q.NumberAbsolute k_theta(final nominal=1) = 1
          "<html>Adjustment factor for thermal resistivity (<i>k</i><sub>&theta;</sub>)</html>"
          annotation (Dialog(group="Material properties"));
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="O2",
          Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol></p>

<p>For more information, see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"));

      end Calibrated;

      model Correlated "Correlated properties"
        extends Species(redeclare replaceable package Data =
              Characteristics.O2.Gas (b_v=[1], n_v={-1,0}));
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="O2",
          Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
          </ol></p>

<p>For more information, see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"));

      end Correlated;

      model Fixed "Fixed properties"
        extends Species(
          redeclare replaceable package Data = FCSys.Characteristics.O2.Gas (
                b_v=[1],n_v={-1,0}),
          redeclare parameter Q.TimeAbsolute tauprime=Data.tauprime(),
          redeclare parameter Q.Mobility mu=Data.mu(),
          redeclare parameter Q.TimeAbsolute nu=Data.nu(),
          redeclare parameter Q.ResistivityMaterial eta=Data.eta(),
          redeclare parameter Q.Fluidity beta=Data.beta(),
          redeclare parameter Q.Fluidity zeta=1/(20.72e-6*U.Pa*U.s),
          redeclare parameter Q.ResistivityThermal theta=U.m*U.K/(26.8e-3*U.W));

        // See the documentation for a table of values.
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="O2",
          Documentation(info="<html><p>Assumptions:<ol>
    <li>Ideal gas</li>
        <li>The generalized resistivities (&beta;, &zeta;, &theta;) are fixed (e.g., independent of temperature).</li>

    </ol></p>

<p>Additional notes:
<ul>
          <li>
  The specific heat capacity is not fixed because it would
  affect the chemical potential and result in an incorrect cell
  potential.</li></ul></p>

  <p>The default resistivities (<code>zeta=1/(207.2e-7*U.Pa*U.s)</code> and <code>theta=U.m*U.K/(26.8e-3*U.W)</code>) are based on data of gas at 1&nbsp;atm and
  300&nbsp;K from Incropera and DeWitt [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 920&ndash;921].
  <a href=\"#Tab1\">Table 1</a> lists the properties at other temperatures. <a href=\"#Tab2\">Table 2</a> lists
  values of the material resistivity or self diffusion coefficient.</p>

  <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
  <caption align=\"top\" id=\"Tab1\">Table 1: Properties of O<sub>2</sub> gas at 1&nbsp;atm
  [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 920&ndash;921]</caption>
  <tr>
      <th valign=\"middle\"><code>T<br>/U.K</code></th>
      <th width=1><code>c_p*U.kg*U.K<br>/(U.J*Data.m)</code></th>
      <th width=1><code>zeta<br>*U.Pa*U.s</code></th>
      <th width=1><code>theta*U.W<br>/(U.m*U.K)</code></th>
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

<br>

  <table border=\"1\" cellspacing=0 cellpadding=2 style=\"border-collapse:collapse;\">
  <caption align=\"top\" id=\"Tab2\">Table 2: Material resistivity of O<sub>2</sub> gas at 1&nbsp;atm
  [<a href=\"modelica://FCSys.UsersGuide.References\">Present1958</a>, p. 263]</caption>
  <tr>
      <th valign=\"middle\"><code>T<br>/U.K</code></th>
      <th width=1><code>eta*U.s<br>/U.cm^2</code></th>
    </tr>
<tr><td>77.7</td><td>0.0153</td></tr>
<tr><td>194.7</td><td>0.104</td></tr>
<tr><td>273.2</td><td>0.187</td></tr>
<tr><td>353.2</td><td>0.301</td></tr>
  </table></p>

<p>For more information, see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"),

          Icon(graphics));

      end Fixed;

    end Gas;

  end O2;

  model SpeciesSolid
    "<html><a href=\"modelica://FCSys.Species.Species\">Species</a> model for a solid (inert and zero velocity)</html>"
    extends SpeciesIsochoric(
      final upstreamX=false,
      final upstreamY=false,
      final upstreamZ=false,
      final phi_IC=zeros(3),
      final I_IC,
      final consMaterial=Conservation.IC,
      final consTransX=Conservation.IC,
      final consTransY=Conservation.IC,
      final consTransZ=Conservation.IC,
      final tauprime=0,
      final beta=1,
      final zeta=1);
    // Note:  beta and zeta don't matter as long as they are nonzero.
    annotation (
      defaultComponentPrefixes="replaceable",
      defaultComponentName="species",
      Documentation(info="<html><p>Assumptions:<ol>
  <li>Zero dynamic compressibility (&rArr; uniform velocity in the axial direction)</li>
  <li>Zero fluidity (&rArr; no shearing)</li></ol></p>

  <p>For more information, see the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"));

  end SpeciesSolid;

  model SpeciesIsochoric
    "<html><a href=\"modelica://FCSys.Species.Species\">Species</a> model for an isochoric material</html>"
    extends Species(
      invertEOS=false,
      initMaterial=InitScalar.Volume,
      final eta=1);

    // Note:  Pressure, which is the default material IC for the base model,
    // can't be used to initialize an incompressible species.
    annotation (
      defaultComponentPrefixes="replaceable",
      defaultComponentName="species",
      Documentation(info="<html>
  <p>Please see the documentation of the <a href=\"modelica://FCSys.Species.Species\">Species</a> model.</p></html>"));

  end SpeciesIsochoric;

  model Species
    "Model to exchange, transport, and store the material, momentum, and energy of one species"
    import FCSys.BaseClasses.Utilities.cartWrap;
    import FCSys.BaseClasses.Utilities.inSign;
    import FCSys.BaseClasses.Utilities.Delta;
    import FCSys.BaseClasses.Utilities.Sigma;
    import assert = FCSys.BaseClasses.Utilities.assertEval;
    extends FCSys.BaseClasses.Icons.Names.Top4;

    // Geometry
    parameter Integer n_faces(
      min=1,
      max=3) = 1
      "<html>Number of pairs of faces (<i>n</i><sub>faces</sub>)</html>"
      annotation (Dialog(group="Geometry"), HideResult=true);
    // Note:  This can't be an outer parameter in Dymola 7.4.

    // Material properties
    replaceable package Data = Characteristics.BaseClasses.Characteristic
      constrainedby Characteristics.BaseClasses.Characteristic
      "Characteristic data" annotation (
      Dialog(group="Material properties"),
      choicesAllMatching=true,
      __Dymola_choicesFromPackage=true);
    Q.TimeAbsolute tauprime(nominal=1e-6*U.s) = Data.tauprime(T, v)
      "<html>Phase change interval (&tau;&prime;)</html>"
      annotation (Dialog(group="Material properties"));
    Q.Mobility mu(nominal=0.1*U.C*U.s/U.kg) = Data.mu(T, v)
      "<html>Mobility (&mu;)</html>"
      annotation (Dialog(group="Material properties"));
    Q.TimeAbsolute nu(nominal=1e-9*U.s) = Data.nu(T, v)
      "<html>Thermal independity (&nu;)</html>"
      annotation (Dialog(group="Material properties"));
    Q.ResistivityMaterial eta(nominal=10e-6*U.s/U.m^2) = Data.eta(T, v)
      "<html>Material resistivity (&eta;)</html>"
      annotation (Dialog(group="Material properties"));
    Q.Fluidity beta(nominal=10*U.cm*U.s/U.g) = Data.beta(T, v)
      "<html>Dynamic compressibility (&beta;)</html>"
      annotation (Dialog(group="Material properties"));
    Q.Fluidity zeta(nominal=10*U.cm*U.s/U.g) = Data.zeta(T, v)
      "<html>Fluidity (&zeta;)</html>"
      annotation (Dialog(group="Material properties"));
    Q.ResistivityThermal theta(nominal=10*U.cm/U.A) = Data.theta(T, v)
      "<html>Thermal resistivity (&theta;)</html>"
      annotation (Dialog(group="Material properties"));

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
    parameter BaseClasses.Conservation consMaterial=Conservation.Dynamic
      "Material" annotation (Dialog(
        Evaluate=true,
        tab="Assumptions",
        group="Formulation of conservation equations"));
    parameter BaseClasses.Conservation consTransX=Conservation.Dynamic
      "X-axis translational momentum" annotation (Evaluate=true, Dialog(
        tab="Assumptions",
        group="Formulation of conservation equations",
        enable=inclTrans[1]));
    parameter BaseClasses.Conservation consTransY=Conservation.Dynamic
      "Y-axis translational momentum" annotation (Evaluate=true, Dialog(
        tab="Assumptions",
        group="Formulation of conservation equations",
        enable=inclTrans[2]));
    parameter BaseClasses.Conservation consTransZ=Conservation.Dynamic
      "Z-axis translational momentum" annotation (Evaluate=true, Dialog(
        tab="Assumptions",
        group="Formulation of conservation equations",
        enable=inclTrans[3]));
    parameter BaseClasses.Conservation consEnergy=Conservation.Dynamic "Energy"
      annotation (Evaluate=true, Dialog(tab="Assumptions", group=
            "Formulation of conservation equations"));
    //
    // Flow conditions
    parameter Q.NumberAbsolute Nu_Phi[Axis]={4,4,4}
      "<html>Translational Nusselt numbers (<b><i>Nu</i><sub>&Phi;</sub></b>)</html>"
      annotation (Dialog(tab="Assumptions", group="Flow conditions"));
    parameter Q.NumberAbsolute Nu_Q=1
      "<html>Thermal Nusselt number (<i>Nu</i><sub><i>Q</i></sub>)</html>"
      annotation (Dialog(tab="Assumptions", group="Flow conditions"));

    // Initialization parameters
    // -------------------------
    // Scalar properties
    parameter BaseClasses.InitScalar initMaterial=InitScalar.Pressure
      "Method of initializing the material state" annotation (Evaluate=true,
        Dialog(tab="Initialization", group="Material and energy"));
    parameter BaseClasses.InitScalar initEnergy=InitScalar.Temperature
      "Method of initializing the thermal state" annotation (Evaluate=true,
        Dialog(tab="Initialization", group="Material and energy"));
    parameter Q.Amount N_IC(start=V_IC*rho_IC)
      "<html>Initial particle number (<i>N</i><sub>IC</sub>)</html>"
      annotation (Dialog(tab="Initialization", group="Material and energy"));
    // Note:  This parameter is left enabled even it isn't used to
    // explicitly initialize any states, since it's used as a guess value.
    // Similar notes apply to some other initial conditions below.
    parameter Q.Density rho_IC(start=1/Data.v_Tp(T_IC, p_IC))
      "<html>Initial density (&rho;<sub>IC</sub>)</html>"
      annotation (Dialog(tab="Initialization", group="Material and energy"));
    parameter Q.Volume V_IC(start=product(L))
      "<html>Initial volume (<i>V</i><sub>IC</sub>)</html>"
      annotation (Dialog(tab="Initialization", group="Material and energy"));
    parameter Q.PressureAbsolute p_IC(start=environment.p)
      "<html>Initial pressure (<i>p</i><sub>IC</sub>)</html>"
      annotation (Dialog(tab="Initialization", group="Material and energy"));
    parameter Q.TemperatureAbsolute T_IC(start=environment.T)
      "<html>Initial temperature (<i>T</i><sub>IC</sub>)</html>"
      annotation (Dialog(tab="Initialization", group="Material and energy"));
    parameter Q.Potential h_IC(start=Data.h(T_IC, p_IC), displayUnit="kJ/mol")
      "<html>Initial specific enthalpy (<i>h</i><sub>IC</sub>)</html>"
      annotation (Dialog(tab="Initialization", group="Material and energy"));
    parameter Q.Potential g_IC(start=Data.g(T_IC, p_IC), displayUnit="kJ/mol")
      "<html>Initial Gibbs potential (<i>g</i><sub>IC</sub>)</html>"
      annotation (Dialog(tab="Initialization", group="Material and energy"));
    //
    // Velocity
    parameter BaseClasses.InitTranslational initTransX=InitTranslational.Velocity
      "Method of initializing the x-axis state" annotation (Evaluate=true,
        Dialog(
        tab="Initialization",
        group="Translational momentum",
        enable=inclTrans[1]));
    parameter BaseClasses.InitTranslational initTransY=InitTranslational.Velocity
      "Method of initializing the y-axis state" annotation (Evaluate=true,
        Dialog(
        tab="Initialization",
        group="Translational momentum",
        enable=inclTrans[2]));
    parameter BaseClasses.InitTranslational initTransZ=InitTranslational.Velocity
      "Method of initializing the z-axis state" annotation (Evaluate=true,
        Dialog(
        tab="Initialization",
        group="Translational momentum",
        enable=inclTrans[3]));
    // Note:  Dymola 7.4 doesn't provide pull-down lists for arrays of
    // enumerations; therefore, a parameter is used for each axis.
    parameter Q.Velocity phi_IC[Axis]={0,0,0}
      "<html>Initial velocity (<b>&phi;</b><sub>IC</sub>)</html>"
      annotation (Dialog(tab="Initialization", group="Translational momentum"));
    parameter Q.Current I_IC[Axis]={0,0,0}
      "<html>Initial current (<i><b>I</b></i><sub>IC</sub>)</html>"
      annotation (Dialog(tab="Initialization", group="Translational momentum"));

    // Advanced parameters
    parameter Boolean invertEOS=true "Invert the equation of state" annotation
      (Dialog(
        Evaluate=true,
        tab="Advanced",
        compact=true), choices(__Dymola_checkBox=true));

    // Preferred states
    // Note:  The start value for this variable (and others below) isn't fixed
    // because the related initial condition is applied in the initial
    // equation section.
    Q.Amount N(
      min=Modelica.Constants.small,
      nominal=4*U.C,
      final start=N_IC,
      final fixed=false,
      stateSelect=StateSelect.prefer) "Particle number";
    Q.Velocity phi[n_trans](
      each nominal=10*U.cm/U.s,
      final start=phi_IC[cartTrans],
      each final fixed=false,
      each stateSelect=StateSelect.prefer) "Velocity";
    Q.TemperatureAbsolute T(
      nominal=300*U.K,
      final start=T_IC,
      final fixed=false,
      stateSelect=StateSelect.prefer) "Temperature";

    // Aliases (for common terms)
    Q.PressureAbsolute p(
      nominal=U.atm,
      final start=p_IC,
      final fixed=false,
      stateSelect=StateSelect.never) "Pressure";
    // StateSelect.never avoids dynamic state selection of this variable and
    // others below in Dymola 7.4.
    Q.Volume V(
      nominal=U.cc,
      final start=V_IC,
      final fixed=false,
      stateSelect=StateSelect.never) "Volume";
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
    Q.Current I[n_trans](
      each nominal=U.A,
      final start=I_IC[cartTrans],
      each final fixed=false,
      stateSelect=StateSelect.never) "Current";
    Q.Current Ndot_faces[n_faces, Side](nominal=U.A, final start=outerProduct(
          I_IC[cartFaces], {1,-1}))
      "Total current into the faces (advection and diffusion)";
    Q.PressureAbsolute p_faces[n_faces, Side](each nominal=U.atm, each start=
          p_IC) "Thermodynamic pressures at the faces";

    // Auxiliary variables (for analysis)
    // ----------------------------------
    // Misc. properties and conditions
    output Q.Density rho(stateSelect=StateSelect.never) = 1/v if environment.analysis
      "Density";
    output Q.MassVolumic mrho(stateSelect=StateSelect.never) = Data.m*rho if
      environment.analysis "Volumic mass";
    output Q.Potential g(stateSelect=StateSelect.never) = chemical.mu if
      environment.analysis "Electrochemical potential";
    output Q.Amount S(stateSelect=StateSelect.never) = N*s if environment.analysis
      "Entropy";
    output Q.PressureAbsolute q[n_trans](each stateSelect=StateSelect.never) =
      Data.m*phi .* I ./ (2*A[cartTrans]) if environment.analysis
      "Dynamic pressure";
    output Q.CapacityThermalSpecific c_p(stateSelect=StateSelect.never) =
      Data.c_p(T, p) if environment.analysis "Isobaric specific heat capacity";
    output Q.CapacityThermalSpecific c_v(stateSelect=StateSelect.never) =
      Data.c_v(T, p) if environment.analysis "Isochoric specific heat capacity";
    output Q.PressureReciprocal kappa(stateSelect=StateSelect.never) =
      Data.kappa(T, p) if environment.analysis "Isothermal compressibility";
    //
    // Time constants (only for the axes with translational momentum included;
    // others are infinite)
    output Q.TimeAbsolute tau_NE(
      stateSelect=StateSelect.never,
      start=U.s) = kappa*tauprime*exp((chemical.mu - g0)/T)*T/v if environment.analysis
      "Time constant for phase change";
    output Q.TimeAbsolute tau_PhiE(
      stateSelect=StateSelect.never,
      start=U.s) = Data.m*mu if environment.analysis
      "Time constant for translational exchange";
    output Q.TimeAbsolute tau_QE(
      stateSelect=StateSelect.never,
      start=U.s) = c_p*nu if environment.analysis
      "Time constant for thermal exchange";
    output Q.TimeAbsolute tau_NT[n_faces](
      each stateSelect=StateSelect.never,
      each start=U.s) = fill(V*eta/2, n_faces) ./ Lprime[cartFaces] if
      environment.analysis "Time constants for material transport";
    output Q.TimeAbsolute tau_PhiT_perp[n_faces](
      each stateSelect=StateSelect.never,
      each start=U.s) = fill(M*beta/2, n_faces) ./ Lprime[cartFaces] if
      environment.analysis "Time constants for normal translational transport";
    output Q.TimeAbsolute tau_PhiT_para[n_faces](
      each stateSelect=StateSelect.never,
      each start=U.s) = fill(M*zeta/2, n_faces) ./ Lprime[cartFaces] if
      environment.analysis
      "Time constants for transverse translational transport";
    output Q.TimeAbsolute tau_QT[n_faces](
      each stateSelect=StateSelect.never,
      each start=U.s) = fill(N*c_v*theta/2, n_faces) ./ Lprime[cartFaces] if
      environment.analysis "Time constants for thermal transport";
    //
    // Peclet numbers (only for the axes with translational momentum included;
    // others are zero)
    output Q.Number Pe_N[n_trans](each stateSelect=StateSelect.never) = eta*v*I
       ./ Lprime[cartTrans] if environment.analysis "Material Peclet numbers";
    output Q.Number Pe_Phi_perp[n_trans](each stateSelect=StateSelect.never) =
      beta*Data.m*I ./ Lprime[cartTrans] if environment.analysis
      "Normal translational Peclet numbers";
    output Q.Number Pe_Phi_para[n_trans](each stateSelect=StateSelect.never) =
      zeta*Data.m*I ./ Lprime[cartTrans] if environment.analysis
      "Transverse translational Peclet numbers";
    output Q.Number Pe_Q[n_trans](each stateSelect=StateSelect.never) = theta*
      Data.c_v(T, p)*I ./ Lprime[cartTrans] if environment.analysis
      "Thermal Peclet numbers";
    //
    // Bulk flow rates
    output Q.Force mphiI[n_trans, n_trans](each stateSelect=StateSelect.never)
       = outerProduct(I, Data.m*phi) if environment.analysis
      "Bulk rate of translational advection (1st index: transport axis, 2nd index: translational component)";
    output Q.Power hI[n_trans](each stateSelect=StateSelect.never) = h*I if
      environment.analysis "Bulk enthalpy flow rate";
    //
    // Translational momentum balance
    output Q.Force Ma[n_trans](each stateSelect=StateSelect.never) = M*(der(phi)
      /U.s + environment.a[cartTrans]) + N*Data.z*environment.E[cartTrans] if
      environment.analysis
      "Acceleration force (including acceleration due to body forces)";
    output Q.Force f_thermo[n_trans](each stateSelect=StateSelect.never) = {(
      if inclFaces[cartTrans[i]] then -Delta(p_faces[facesCart[cartTrans[i]], :])
      *A[cartTrans[i]] else 0) for i in 1:n_trans} "Thermodynamic force";
    output Q.Force f_AE[n_trans](each stateSelect=StateSelect.never) = Data.m*(
      (actualStream(chemical.phi) - phi) .* chemical.Ndot + (actualStream(
      physical.phi) - phi) .* physical.Ndot) if environment.analysis
      "Acceleration force due to advective exchange";
    output Q.Force f_DE[n_trans](each stateSelect=StateSelect.never) = inert.translational.mPhidot
       + dalton.mPhidot if environment.analysis
      "Friction from other configurations (diffusive exchange)";
    output Q.Force f_AT[n_trans](each stateSelect=StateSelect.never) = {sum((
      faces[j, :].phi[cartWrap(cartTrans[i] - cartFaces[j] + 1)] - {phi[i],phi[
      i]})*Ndot_faces[j, :]*Data.m for j in 1:n_faces) for i in 1:n_trans} if
      environment.analysis "Acceleration force due to advective transport";
    output Q.Force f_DT[n_trans](each stateSelect=StateSelect.never) = {sum(
      Sigma(faces[j, :].mPhidot[cartWrap(cartTrans[i] - cartFaces[j] + 1)])
      for j in 1:n_faces) for i in 1:n_trans} if environment.analysis
      "Friction from other subregions (diffusive transport, including bulk viscosity)";
    //
    // Energy balance
    output Q.Power Ndere(stateSelect=StateSelect.never) = (N*T*der(s) + M*phi*
      der(phi))/U.s if environment.analysis
      "Rate of energy storage (internal and kinetic) and boundary work at constant mass";
    // Note that T*der(s) = der(u) + p*der(v).
    output Q.Power Edot_AE(stateSelect=StateSelect.never) = (chemical.mu +
      actualStream(chemical.sT) - h + (actualStream(chemical.phi)*actualStream(
      chemical.phi) - phi*phi)*Data.m/2)*chemical.Ndot + (physical.mu +
      actualStream(physical.sT) - h + (actualStream(physical.phi)*actualStream(
      physical.phi) - phi*phi)*Data.m/2)*physical.Ndot if environment.analysis
      "Relative rate of energy (internal, flow, and kinetic) due to phase change and reaction";
    output Q.Power Edot_DE(stateSelect=StateSelect.never) = inert.translational.phi
      *inert.translational.mPhidot + inert.thermal.Qdot + dalton.phi*dalton.mPhidot
       + dalton.Qdot if environment.analysis
      "Rate of diffusion of energy from other configurations";
    output Q.Power Edot_AT(stateSelect=StateSelect.never) = sum((Data.h(faces[j,
      :].T, p_faces[j, :]) - {h,h})*Ndot_faces[j, :] + sum((faces[j, :].phi[
      cartWrap(cartTrans[i] - cartFaces[j] + 1)] .^ 2 - fill(phi[i]^2, 2))*
      Ndot_faces[j, :]*Data.m/2 for i in 1:n_trans) for j in 1:n_faces) if
      environment.analysis
      "Relative rate of energy (internal, flow, and kinetic) due to advective transport";
    output Q.Power Edot_DT(stateSelect=StateSelect.never) = sum(sum(faces[j, :].phi[
      cartWrap(cartTrans[i] - cartFaces[j] + 1)]*faces[j, :].mPhidot[cartWrap(
      cartTrans[i] - cartFaces[j] + 1)] for i in 1:n_trans) for j in 1:n_faces)
       + sum(faces.Qdot) if environment.analysis
      "Rate of diffusion of energy from other subregions";
    // Note:  The structure of the problem should not change if these
    // auxiliary variables are included (hence StateSelect.never).

    Connectors.Chemical chemical(
      final n_trans=n_trans,
      mu(start=g_IC, final fixed=false),
      phi(start=phi_IC[cartTrans], each final fixed=false),
      sT(start=h_IC - g_IC, final fixed=false)) "Connector for reactions"
      annotation (Placement(transformation(extent={{-10,10},{10,30}}),
          iconTransformation(extent={{-29,60},{-49,80}})));
    Connectors.Physical physical(
      final formula=Data.formula,
      final n_trans=n_trans,
      mu(start=g_IC, fixed=false),
      phi(final start=phi_IC[cartTrans], each final fixed=false),
      sT(final start=h_IC - g_IC, final fixed=false))
      "Connector for phase change" annotation (Placement(transformation(extent=
              {{-30,-10},{-10,10}}),iconTransformation(extent={{-59,28},{-79,48}})));
    Connectors.Inert inert(
      final n_trans=n_trans,
      translational(phi(final start=phi_IC[cartTrans], each final fixed=false)),

      thermal(T(final start=T_IC, final fixed=false)))
      "Connector to directly couple velocity or temperature with other species"
      annotation (Placement(transformation(extent={{10,-10},{30,10}}),
          iconTransformation(extent={{59,-28},{79,-48}})));

    FCSys.Connectors.Dalton dalton(
      final n_trans=n_trans,
      V(
        min=0,
        final start=V_IC,
        final fixed=false),
      p(final start=p_IC, final fixed=false),
      phi(start=phi_IC[cartTrans]),
      T(start=T_IC))
      "Connector for translational and thermal diffusive exchange, with additivity of pressure"
      annotation (Placement(transformation(extent={{-10,-30},{10,-10}}),
          iconTransformation(extent={{28,-59},{48,-79}})));
    Connectors.Face faces[n_faces, Side](
      rho(each start=rho_IC),
      Ndot(start=outerProduct(I_IC[cartFaces], {1,-1})),
      phi(start={fill({phi_IC[cartWrap(cartFaces[i] + orientation - 1)] for
            orientation in Orientation}, 2) for i in 1:n_faces}),
      mPhidot(each start=0),
      T(each start=T_IC),
      Qdot(each start=0))
      "Connectors to transport material, translational momentum, and thermal energy through the boundaries"
      annotation (Placement(transformation(extent={{-10,-10},{10,10}}),
          iconTransformation(extent={{-10,-10},{10,10}})));

    // Geometric parameters

  protected
    outer parameter Q.Length L[Axis] "Lengths" annotation (missingInnerMessage="This model should be used within a subregion model.
");
    outer parameter Q.Area A[Axis] "Cross-sectional areas" annotation (
        missingInnerMessage="This model should be used within a subregion model.
");
    outer parameter Q.Length Lprime[Axis]
      "Effective cross-sectional area per length" annotation (
        missingInnerMessage="This model should be used within a phase model.
      ");
    outer parameter Boolean inclTrans[Axis]
      "true, if each component of translational momentum is included"
      annotation (missingInnerMessage="This model should be used within a subregion model.
");
    outer parameter Boolean inclFaces[Axis]
      "true, if each pair of faces is included" annotation (missingInnerMessage
        ="This model should be used within a subregion model.
");
    outer parameter Boolean inclRot[3]
      "true, if each axis of rotation has all its tangential faces included"
      annotation (missingInnerMessage="This model should be used within a subregion model.
");
    // Note:  The size is also Axis, but it can't be specified here due to
    // an error in Dymola 7.4 (failure in check of Phase models).
    outer parameter Integer n_trans
      "Number of components of translational momentum" annotation (
        missingInnerMessage="This model should be used within a subregion model.
");
    outer parameter Integer cartTrans[:]
      "Cartesian-axis indices of the components of translational momentum"
      annotation (missingInnerMessage="This model should be used within a subregion model.
");
    outer parameter Integer cartFaces[:]
      "Cartesian-axis indices of the pairs of faces" annotation (
        missingInnerMessage="This model should be used within a subregion model.
");
    outer parameter Integer cartRot[:]
      "Cartesian-axis indices of the components of rotational momentum"
      annotation (missingInnerMessage="This model should be used within a subregion model.
");
    // Note:  The size of cartTrans, cartFaces, and cartRot is n_trans,
    // but it can't be specified here due to an error in Dymola 7.4.
    outer parameter Integer transCart[3]
      "Translational-momentum-component indices of the Cartesian axes"
      annotation (missingInnerMessage="This model should be used within a subregion model.
");
    // Note:  The size is also Axis, but it can't be specified here due to
    // an error in Dymola 7.4 (failure in check of Phase models).
    outer parameter Integer facesCart[Axis]
      "Face-pair indices of the Cartesian axes" annotation (missingInnerMessage
        ="This model should be used within a subregion model.
");
    outer parameter Q.NumberAbsolute k_E "Coupling factor for exchange"
      annotation (missingInnerMessage="This model should be used within a phase model.
      ");
    final parameter Boolean upstream[Axis]={upstreamX,upstreamY,upstreamZ}
      "true, if each Cartesian axis uses upstream discretization"
      annotation (HideResult=true);
    final parameter BaseClasses.Conservation consTrans[Axis]={consTransX,
        consTransY,consTransZ}
      "Formulation of the translational conservation equations"
      annotation (HideResult=true);
    final parameter BaseClasses.InitTranslational initTrans[Axis]={initTransX,
        initTransY,initTransZ}
      "Initialization methods for translational momentum"
      annotation (HideResult=true);

    // Additional aliases (for common terms)
    Q.Potential g0(
      nominal=U.V,
      final start=g_IC,
      final fixed=false,
      stateSelect=StateSelect.never) "Gibbs potential at reference pressure";
    Q.Force faces_mPhidot[n_faces, Side, 2] "Directly calculated shear forces";
    Q.Velocity phi_actual_chemical[n_trans] "Velocity of the chemical stream";
    Q.Velocity phi_actual_physical[n_trans] "Velocity of the physical stream";
    // Note:  Dymola 7.4 can't individually index the components of a
    // stream variable (e.g., actualStream(chemical.phi[i])), so these
    // variables are necessary.

    outer Conditions.Environment environment "Environmental conditions";

  initial equation
    // Check the initial conditions.
    assert(initMaterial <> initEnergy or initMaterial == InitScalar.None or
      consMaterial == Conservation.Steady or consEnergy == Conservation.Steady,
      "The initialization methods for material and energy must be different (unless None).");

    // Material
    if consMaterial == Conservation.IC then
      // Ensure that a condition is selected since the state is prescribed.
      assert(initMaterial <> InitScalar.None, "The material state of " + Data.formula
         + " is prescribed, yet its condition is not defined.
Choose any condition besides None.");
    elseif consMaterial == Conservation.Dynamic then
      // Initialize since there's a time-varying state.
      if initMaterial == InitScalar.Amount then
        N = N_IC;
      elseif initMaterial == InitScalar.AmountSS then
        der(N) = 0;
      elseif initMaterial == InitScalar.Density then
        1/v = rho_IC;
        assert(Data.isCompressible or Data.hasThermalExpansion, Data.formula +
          " is isochoric, yet its material initial condition is based on density.");
      elseif initMaterial == InitScalar.DensitySS then
        der(1/v) = 0;
        assert(Data.isCompressible or Data.hasThermalExpansion, Data.formula +
          " is isochoric, yet its material initial condition is based on density.");
      elseif initMaterial == InitScalar.Volume then
        V = V_IC;
      elseif initMaterial == InitScalar.VolumeSS then
        der(V) = 0;
      elseif initMaterial == InitScalar.Pressure then
        assert(Data.isCompressible, Data.formula +
          " is incompressible, yet its material initial condition is based on pressure.");
        p = p_IC;
        assert(Data.isCompressible, Data.formula +
          " is incompressible, yet its material initial condition is based on pressure.");
      elseif initMaterial == InitScalar.PressureSS then
        der(p) = 0;
      elseif initMaterial == InitScalar.Temperature then
        T = T_IC;
      elseif initMaterial == InitScalar.TemperatureSS then
        der(T) = 0;
      elseif initMaterial == InitScalar.SpecificEnthalpy then
        h = h_IC;
      elseif initMaterial == InitScalar.SpecificEnthalpySS then
        der(h) = 0;
      elseif initMaterial == InitScalar.PotentialGibbs then
        chemical.mu = g_IC;
      elseif initMaterial == InitScalar.PotentialGibbsSS then
        der(chemical.mu) = 0;
        // Else there's no initial equation since
        // initMaterial == InitScalar.None or
        // consMaterial == Conservation.Steady.
      end if;
    end if;

    // Translational momentum
    for i in 1:n_trans loop
      if consTrans[cartTrans[i]] == Conservation.IC then
        // Ensure that a condition is selected since the state is
        // prescribed.
        assert(initTrans[cartTrans[i]] <> InitTranslational.None,
          "The state for the " + (if cartTrans[i] == Axis.x then "x" else if
          cartTrans[i] == Axis.y then "y" else "z") +
          "-axis component of translational momentum of " + Data.formula + " is prescribed, yet its condition is not defined.
Choose any condition besides None.");
      elseif consTrans[cartTrans[i]] == Conservation.Dynamic then
        // Initialize since there's a time-varying state.
        if initTrans[cartTrans[i]] == InitTranslational.Velocity then
          phi[i] = phi_IC[cartTrans[i]];
        elseif initTrans[cartTrans[i]] == InitTranslational.VelocitySS then
          der(phi[i]) = 0;
        elseif initTrans[cartTrans[i]] == InitTranslational.Current then
          I[i] = I_IC[cartTrans[i]];
        elseif initTrans[cartTrans[i]] == InitTranslational.CurrentSS then
          der(I[i]) = 0;
          // Else there's no initial equation since
          // initTrans[cartTrans[i]] == InitTranslational.None or
          // consTrans[cartTrans[i]] == Conservation.Steady.
        end if;
      end if;
    end for;

    // Energy
    if consEnergy == Conservation.IC then
      // Ensure that a condition is selected since the state is prescribed.
      assert(initEnergy <> InitScalar.None, "The energy state of " + Data.formula
         + " is prescribed, yet its condition is not defined.
Choose any condition besides None.");
    elseif consEnergy == Conservation.Dynamic then
      // Initialize since there's a time-varying state.
      if initEnergy == InitScalar.Amount then
        N = N_IC;
      elseif initEnergy == InitScalar.AmountSS then
        der(N) = 0;
      elseif initEnergy == InitScalar.Density then
        1/v = rho_IC;
        assert(Data.isCompressible or Data.hasThermalExpansion, Data.formula +
          " is isochoric, yet its thermal initial condition is based on density.");
      elseif initEnergy == InitScalar.DensitySS then
        der(1/v) = 0;
        assert(Data.isCompressible or Data.hasThermalExpansion, Data.formula +
          " is isochoric, yet its thermal initial condition is based on density.");
      elseif initEnergy == InitScalar.Volume then
        V = V_IC;
      elseif initEnergy == InitScalar.VolumeSS then
        der(V) = 0;
      elseif initEnergy == InitScalar.Pressure then
        p = p_IC;
        assert(Data.isCompressible, Data.formula +
          " is incompressible, yet its thermal initial condition is based on pressure.");
      elseif initEnergy == InitScalar.PressureSS then
        der(p) = 0;
        assert(Data.isCompressible, Data.formula +
          " is incompressible, yet its thermal initial condition is based on pressure.");
      elseif initEnergy == InitScalar.Temperature then
        T = T_IC;
      elseif initEnergy == InitScalar.TemperatureSS then
        der(T) = 0;
      elseif initEnergy == InitScalar.SpecificEnthalpy then
        h = h_IC;
      elseif initEnergy == InitScalar.SpecificEnthalpySS then
        der(h) = 0;
      elseif initEnergy == InitScalar.PotentialGibbs then
        chemical.mu = g_IC;
      elseif initEnergy == InitScalar.PotentialGibbsSS then
        der(chemical.mu) = 0;
        // Else there's no initial equation since
        // initEnergy == InitScalar.None or
        // consEnergy == Conservation.Steady.
      end if;
    end if;

  equation
    // Aliases (only to clarify and simplify other equations)
    T = inert.thermal.T;
    p = dalton.p;
    V = dalton.V;
    v*N = V;
    M = Data.m*N;
    phi = inert.translational.phi;
    I .* L[cartTrans] = N*phi;
    p_faces = {{Data.p_Tv(faces[i, side].T, 1/faces[i, side].rho) for side in
      Side} for i in 1:n_faces};
    Ndot_faces = faces.Ndot + faces.rho .* {{faces[i, Side.n].phi[Orientation.normal],
      -faces[i, Side.p].phi[Orientation.normal]}*A[cartFaces[i]] for i in 1:
      n_faces};
    phi_actual_chemical = actualStream(chemical.phi);
    phi_actual_physical = actualStream(physical.phi);

    // Thermodynamic correlations
    if invertEOS then
      p = Data.p_Tv(T, v);
    else
      v = Data.v_Tp(T, p);
    end if;
    h = Data.h(T, p);
    s = Data.s(T, p);
    g0 = Data.g(T, Data.p0);

    // Diffusive exchange
    chemical.mu = h - chemical.sT "Reaction (rate equation elsewhere)";
    if tauprime > Modelica.Constants.small then
      tauprime*physical.Ndot = k_E*N*(exp((physical.mu - g0)/T) - exp((chemical.mu
         - g0)/T)) "Phase change";
    else
      0 = if N > 0 or physical.mu > chemical.mu then physical.mu - chemical.mu
         else physical.Ndot;
      // The first branch avoids nonlinear equations when tauprime=0.
      // Dymola 7.4 can't derive it symbolically from the previous equation.
    end if;
    mu*dalton.mPhidot = k_E*N*(dalton.phi - phi) "Translational momentum";
    nu*dalton.Qdot = k_E*N*(dalton.T - T) "Thermal energy";

    // Properties upon outflow due to reaction and phase change
    chemical.phi = phi;
    physical.phi = phi;
    chemical.sT = s*T;
    physical.sT = chemical.sT;

    // Diffusive transport
    for i in 1:n_faces loop
      for side in Side loop
        // Material (central difference)
        eta*faces[i, side].Ndot = Lprime[cartFaces[i]]*(faces[i, side].rho - 1/
          v)*2;

        // Translational momentum
        beta*faces[i, side].mPhidot[Orientation.normal] = Lprime[cartFaces[i]]*
          (faces[i, side].phi[Orientation.normal] - (if inclTrans[cartFaces[i]]
           then phi[transCart[cartFaces[i]]] else 0))*(if inclTrans[cartFaces[i]]
           and upstream[cartFaces[i]] then 1 + exp(-inSign(side)*I[transCart[
          cartFaces[i]]]*beta*Data.m/(2*Lprime[cartFaces[i]])) else 2) "Normal";
        zeta*faces_mPhidot[i, side, Orientation.following - 1] = Nu_Phi[
          cartFaces[i]]*Lprime[cartFaces[i]]*(faces[i, side].phi[Orientation.following]
           - (if inclTrans[cartWrap(cartFaces[i] + 1)] then phi[transCart[
          cartWrap(cartFaces[i] + 1)]] else 0))*(if inclTrans[cartFaces[i]]
           and upstream[cartFaces[i]] then 1 + exp(-inSign(side)*I[transCart[
          cartFaces[i]]]*zeta*Data.m/(2*Lprime[cartFaces[i]])) else 2)
          "1st transverse";
        zeta*faces_mPhidot[i, side, Orientation.preceding - 1] = Nu_Phi[
          cartFaces[i]]*Lprime[cartFaces[i]]*(faces[i, side].phi[Orientation.preceding]
           - (if inclTrans[cartWrap(cartFaces[i] - 1)] then phi[transCart[
          cartWrap(cartFaces[i] - 1)]] else 0))*(if inclTrans[cartFaces[i]]
           and upstream[cartFaces[i]] then 1 + exp(-inSign(side)*I[transCart[
          cartFaces[i]]]*zeta*Data.m/(2*Lprime[cartFaces[i]])) else 2)
          "2nd transverse";

        // Thermal energy
        theta*faces[i, side].Qdot = Nu_Q*Lprime[cartFaces[i]]*(faces[i, side].T
           - T)*(if inclTrans[cartFaces[i]] and upstream[cartFaces[i]] then 1
           + exp(-inSign(side)*I[transCart[cartFaces[i]]]*theta*Data.c_v(T, p)/
          (2*Lprime[cartFaces[i]])) else 2);
      end for;

      // Direct mapping of transverse forces (calculated above)
      if not inclRot[cartWrap(cartFaces[i] - 1)] then
        faces[i, :].mPhidot[Orientation.following] = faces_mPhidot[i, :,
          Orientation.following - 1];
        // Else the force must be mapped for zero torque (below).
      end if;
      if not inclRot[cartWrap(cartFaces[i] + 1)] then
        faces[i, :].mPhidot[Orientation.preceding] = faces_mPhidot[i, :,
          Orientation.preceding - 1];
        // Else the force must be mapped for zero torque (below).
      end if;
    end for;

    // Zero-torque mapping of transverse forces
    for axis in cartRot loop
      4*cat(
          1,
          faces[facesCart[cartWrap(axis + 1)], :].mPhidot[Orientation.following],
          faces[facesCart[cartWrap(axis - 1)], :].mPhidot[Orientation.preceding])
        = {{3,1,L[cartWrap(axis - 1)]/L[cartWrap(axis + 1)],-L[cartWrap(axis -
        1)]/L[cartWrap(axis + 1)]},{1,3,-L[cartWrap(axis - 1)]/L[cartWrap(axis
         + 1)],L[cartWrap(axis - 1)]/L[cartWrap(axis + 1)]},{L[cartWrap(axis +
        1)]/L[cartWrap(axis - 1)],-L[cartWrap(axis + 1)]/L[cartWrap(axis - 1)],
        3,1},{-L[cartWrap(axis + 1)]/L[cartWrap(axis - 1)],L[cartWrap(axis + 1)]
        /L[cartWrap(axis - 1)],1,3}}*cat(
          1,
          faces_mPhidot[facesCart[cartWrap(axis + 1)], :, Orientation.following
           - 1],
          faces_mPhidot[facesCart[cartWrap(axis - 1)], :, Orientation.preceding
           - 1]);
    end for;

    // Material dynamics
    if consMaterial == Conservation.IC then
      // Apply the IC forever (material not conserved).
      if initMaterial == InitScalar.Amount then
        N = N_IC;
      elseif initMaterial == InitScalar.AmountSS then
        der(N) = 0;
      elseif initMaterial == InitScalar.Density then
        1/v = rho_IC;
      elseif initMaterial == InitScalar.DensitySS then
        der(1/v) = 0;
      elseif initMaterial == InitScalar.Volume then
        V = V_IC;
      elseif initMaterial == InitScalar.VolumeSS then
        der(V) = 0;
      elseif initMaterial == InitScalar.Pressure then
        p = p_IC;
      elseif initMaterial == InitScalar.PressureSS then
        der(p) = 0;
      elseif initMaterial == InitScalar.Temperature then
        T = T_IC;
      elseif initMaterial == InitScalar.TemperatureSS then
        der(T) = 0;
      elseif initMaterial == InitScalar.SpecificEnthalpy then
        h = h_IC;
      elseif initMaterial == InitScalar.SpecificEnthalpySS then
        der(h) = 0;
      elseif initMaterial == InitScalar.PotentialGibbs then
        chemical.mu = g_IC;
      else
        // if initMaterial == InitScalar.PotentialGibbsSS then
        der(chemical.mu) = 0;
        // Note:  initMaterial == InitScalar.None can't occur due to an
        // assertion.
      end if;
    else
      (if consMaterial == Conservation.Dynamic then der(N)/U.s else 0) =
        chemical.Ndot + physical.Ndot + sum(Ndot_faces) "Material conservation";
    end if;

    // Translational dynamics
    for j in 1:n_trans loop
      if consTrans[cartTrans[j]] == Conservation.IC then
        // Apply the IC forever (translational momentum isn't conserved along
        // this axis).
        if initTrans[cartTrans[j]] == InitTranslational.Velocity then
          phi[j] = phi_IC[cartTrans[j]];
        elseif initTrans[cartTrans[j]] == InitTranslational.VelocitySS then
          der(phi[j]) = 0;
        elseif initTransX == InitTranslational.Current then
          I[j] = I_IC[cartTrans[j]];
        else
          // if initTrans[cartTrans[j]] == InitTranslational.CurrentSS then
          der(I[j]) = 0;
          // Note:  initTrans[cartTrans[j]] == InitTranslational.None can't
          // occur due to an assertion.
        end if;
      else
        M*((if consTrans[cartTrans[j]] == Conservation.Dynamic then der(phi[j])
          /U.s else 0) + environment.a[cartTrans[j]]) + N*Data.z*environment.E[
          cartTrans[j]] + (if inclFaces[cartTrans[j]] then Delta(p_faces[
          facesCart[cartTrans[j]], :])*A[cartTrans[j]] else 0) = Data.m*((
          phi_actual_chemical[j] - phi[j])*chemical.Ndot + (phi_actual_physical[
          j] - phi[j])*physical.Ndot) + inert.translational.mPhidot[j] + dalton.mPhidot[
          j] + sum((faces[i, :].phi[cartWrap(cartTrans[j] - cartFaces[i] + 1)]
           - {phi[j],phi[j]})*Ndot_faces[i, :]*Data.m + Sigma(faces[i, :].mPhidot[
          cartWrap(cartTrans[j] - cartFaces[i] + 1)]) for i in 1:n_faces)
          "Conservation of translational momentum";
        // Note:  Dymola 7.4 (Dassl integrator) runs better with this intensive
        // form of the balance (M*der(phi) = ... rather than der(M*phi) = ...).
      end if;
    end for;

    // Thermal dynamics
    if consEnergy == Conservation.IC then
      // Apply the IC forever (energy not conserved).
      if initEnergy == InitScalar.Amount then
        N = N_IC;
      elseif initEnergy == InitScalar.AmountSS then
        der(N) = 0;
      elseif initMaterial == InitScalar.Density then
        1/v = rho_IC;
      elseif initMaterial == InitScalar.DensitySS then
        der(1/v) = 0;
      elseif initEnergy == InitScalar.Volume then
        V = V_IC;
      elseif initEnergy == InitScalar.VolumeSS then
        der(V) = 0;
      elseif initEnergy == InitScalar.Pressure then
        p = p_IC;
      elseif initEnergy == InitScalar.PressureSS then
        der(p) = 0;
      elseif initEnergy == InitScalar.Temperature then
        T = T_IC;
      elseif initEnergy == InitScalar.TemperatureSS then
        der(T) = 0;
      elseif initEnergy == InitScalar.SpecificEnthalpy then
        h = h_IC;
      elseif initEnergy == InitScalar.SpecificEnthalpySS then
        der(h) = 0;
      elseif initEnergy == InitScalar.PotentialGibbs then
        chemical.mu = g_IC;
      else
        // if initEnergy == InitScalar.PotentialGibbsSS then
        der(chemical.mu) = 0;
        // Note:  initEnergy == InitScalar.None can't occur due to an
        // assertion.
      end if;
    else
      (if consEnergy == Conservation.Dynamic then (N*T*der(s) + M*phi*der(phi))
        /U.s else 0) = (chemical.mu + actualStream(chemical.sT) - h + (
        actualStream(chemical.phi)*actualStream(chemical.phi) - phi*phi)*Data.m
        /2)*chemical.Ndot + (physical.mu + actualStream(physical.sT) - h + (
        actualStream(physical.phi)*actualStream(physical.phi) - phi*phi)*Data.m
        /2)*physical.Ndot + inert.translational.phi*inert.translational.mPhidot
         + inert.thermal.Qdot + dalton.phi*dalton.mPhidot + dalton.Qdot + sum((
        Data.h(faces[i, :].T, p_faces[i, :]) - {h,h})*Ndot_faces[i, :] + sum((
        faces[i, :].phi[cartWrap(cartTrans[j] - cartFaces[i] + 1)] .^ 2 - fill(
        phi[j]^2, 2))*Ndot_faces[i, :]*Data.m/2 + faces[i, :].phi[cartWrap(
        cartTrans[j] - cartFaces[i] + 1)]*faces[i, :].mPhidot[cartWrap(
        cartTrans[j] - cartFaces[i] + 1)] for j in 1:n_trans) for i in 1:
        n_faces) + sum(faces.Qdot) "Conservation of energy";
    end if;
    annotation (
      defaultComponentPrefixes="replaceable",
      Documentation(info="<html>
    <p>This model is based on the following fixed assumptions:
    <ol>
       <li>All faces are rectangular.
       <li>The material is orthorhombic.  This implies that a gradient which induces diffusion
       along an axis does not induce diffusion along axes orthogonal to it
       [<a href=\"modelica://FCSys.UsersGuide.References\">Bejan2006</a>,
       pp. 691&ndash;692].</li>
       <li>The coordinate system (x, y, z) is aligned with the principle
       axes of transport.  For example if the species is stratified, then the
       layers must be parallel to one of the planes in the rectilinear
       grid.</li>
       <li>The factors that may cause anisotropic behavior (<b><i>k</i></b>)
          are common to material, translational, and thermal transport.</li>
       <li>There is no radiative heat transfer (or else it must be linearized).</li>
       <li>Rotational momentum is not exchanged, transported, or stored.</li>
       <li>For the purpose of the material, translational momentum, and energy balances, the
       cross sectional areas of the faces are assumed to be the full cross-sectional
       areas of the subregion.  If multiple phases are present, then the areas are
       actually smaller.</li>
    </ol>
    Other assumptions are optional via the parameters.</p>

    <p><a href=\"#Fig1\">Figure 1</a> shows how instances of
    <a href=\"modelica://FCSys.Species\">Species</a> models (derived from this
    model) are
    connected within a <a href=\"modelica://FCSys.Subregions\">Subregion</a>.  A single species in
    a single phase is called a <i>configuration</i>. The
    generalized resistances (<i>R</i>) affect the phase change rate, forces, and heat flow rates
    associated with differences in activity, velocity, and temperature (respectively) between
    each configuration and a common node.  These exchange processes are diffusive.

    <p>In general, the resistances are included within the
    <a href=\"modelica://FCSys.Species\">Species</a> model.  For reactions, however,
    the rate equation is more complex and is included in the
    <a href=\"modelica://FCSys.Subregions.Reaction\">Reaction</a> model.</p>

    <p>Translational momentum and thermal energy are advected as material is exchanged
    due to phase change or reactions.  This occurs at the velocity (&phi;) and specific entropy-temperature
    product (<i>sT</i>) of the reactants (source configurations), where the reactant/product designation
    depends on the current conditions.</p>

    <p>The advective exchange is modeled using <code>stream</code> connectors
    (<a href=\"modelica://FCSys.Connectors.Physical\">Physical</a> and
    <a href=\"modelica://FCSys.Connectors.Chemical\">Chemical</a>).
  The rate of advection of translational momentum is the
  product of the velocity of the source (&phi;) and the mass flow rate
  (<i>m</i><i>N&#775;</i>).  The rate of thermal advection is the
  specific entropy-temperature product of the source (<i>sT</i>) times the rate of
  material exchange
  (<i>N&#775;</i>).  If there are multiple sources, then
  their contributions are additive.  If there are multiple sinks, then
  translational momentum is split on a mass basis and the thermal stream is split
  on a particle-number basis.</p>

    <p align=center id=\"Fig1\"><img src=\"modelica://FCSys/Resources/Documentation/Subregions/Species/Species/Exchange.png\">
<br>Figure 1:  Exchange of a quantity (material, translational momentum, or thermal energy) among configurations
    (A, B, and C) within a subregion.</p>

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
    total extensive volume of the phase.  Within a <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a>,
    the <a href=\"modelica://FCSys.Phases\">Phases</a> are combined by Amagat's law of partial volumes
    (see the <a href=\"modelica://FCSys.Subregions.Volume\">Volume</a> model), as shown
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
    <ul>
    <li>Here (and in the rest of <a href=\"modelica://FCSys\">FCSys</a>), the <i>specific</i>
    adjective means that the following extensive quantity is divided by particle number.
    (<i>Massic</i> indicates a quantity divided by mass.)</li>
    <li>In general, if material resistivity, dynamic compressibility, fluidity, or thermal resistivity is zero, then
    it should be set as <code>final</code> so that index reduction may be performed.
    If two configurations
    are connected through their <code>dalton</code> connectors or faces
    and both have zero generalized resistivities for a
    quantity, then index reduction [<a href=\"modelica://FCSys.UsersGuide.References\">Mattsson1993B</a>] is necessary.</li>
    <li>Even if an initialization parameter is not selected for explicit use,
    it may be used a guess value.</li>
    <li>The area fill factor (<b><i>k</i></b>) can be used to account for the effects of porosity and tortousity
    on the rate of transport.  It may reflect anisotropic properties, since it is a vector with independent components
    for each axis.
    By default, its components are unity.  The area fill factor should be adjusted directly with effective
    area and inversely with effective length.
    It affects all of the diffusive transport rates (material, translational, and
    thermal) by the same factor.</li>
    <li>If <code>Conservation.IC</code> is used for a state (via
    <code>consMaterial</code>, <code>consTransX</code>, <code>consTransY</code>,
    <code>consTransZ</code>, or <code>consEnergy</code>),
    then the associated initial condition (IC) will be applied forever instead of the
    corresponding conservation equation.
    If <code>consMaterial</code>, <code>consTransX</code>, <code>consTransY</code>, or <code>consTransZ</code> is
    <code>Conservation.IC</code>, then there may be a secondary effect on the energy conservation equation
    and thus temperature.
    In that case, it may help to set <code>consEnergy</code> to <code>Conservation.IC</code> so that
    the energy conservation equation is not imposed.</li>
    <li>If <code>consTransX</code>, <code>consTransY</code>, or <code>consTransZ</code> is
    <code>Conservation.Steady</code>, then the derivative of the corresponding component of velocity
    is treated as zero and removed from the translational momentum balance.  If <code>consEnergy</code> is
    <code>Conservation.Steady</code>, then <code>T*der(s) + M*phi*der(phi)</code> is treated as
    zero and removed from the energy balance.</li>
    <li>If a component of velocity is not included (via the outer <code>inclTrans[:]</code> parameter
    which maps to <code>{inclTransX, inclTransY, inclTransZ}</code> in the
    <a href=\"modelica://FCSys.Subregions.BaseClasses.EmptySubregion\">Subregion</a> model), then it
    is taken to be zero in each translational transport equation.  However, the corresponding forces
    in the <code>faces</code> connector array are not included in the momentum or energy balances.
    If it is necessary to set a component of velocity to zero but still include it in the energy balance, then
    set the corresponding component of <code>phi_IC</code> to zero and <code>consTransX</code>,
    <code>consTransY</code>, or <code>consTransZ</code> to <code>Conservation.IC</code>.</li>
    <li>If a subregion does not contain any compressible species, then pressure must be prescribed.
    Set <code>consMaterial</code> to <code>Conservation.IC</code> and <code>initMaterial</code>
    to <code>InitScalar.Pressure</code> for one of the species.</li>
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
    <li>If <code>invertEOS</code> is <code>true</code>, then the equation of state is implemented with pressure
    as a function of temperature and specific volume.  Otherwise, specific volume is a function of temperature
    and pressure.</li>
    <li>The default thermal Nusselt number is one, which represents pure conduction through the gas.  Use 
    3.66 for internal flow where the boundaries are uniform in temperature or 48/11 or approximately 4.36 
    if the heat flux is uniform [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>].</li></p>

    <p>In the <code>faces</code> connector array, the transverse translational flow (<i>m</i>&Phi;dot) is only the
    force due to diffusion.  Translational advection is calculated from the velocity and the material current.
    The thermal flow (<i>Q&#775;</i>) is only the rate of heat transfer due to diffusion.  The advection of
    thermal energy is determined from the thermodynamic state at the boundary and the material current.</p>


    <p>In evaluating the dynamics of a phase, it is typically assumed that all of the species
    exist at the same velocity and temperature.  The translational and thermal time constants
    are usually much shorter than the time span of interest due to the very small coupling
    resistances.  If this is the case, connect the <code>inert</code>
    connectors of the species.  This will reduce the index of the problem.</p>

    <p>For the variables that relate to transport,
    the first index is the axis and the second index is the side.  The sides
    are ordered from negative to positive, according to the
    <a href=\"modelica://FCSys.BaseClasses.Side\">Side</a> enumeration.
    Velocity and force are additionally indexed by
    the orientation of the momentum with respect to the face.
    The orientations are ordered in Cartesian space starting with the normal axis,
    according to the
    <a href=\"modelica://FCSys.BaseClasses.Orientation\">Orientation</a> enumeration.</p>
    </html>"),
      Diagram(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}},
          initialScale=0.1), graphics),
      Icon(graphics={
          Rectangle(
            extent={{-98,80},{98,120}},
            fillPattern=FillPattern.Solid,
            fillColor={255,255,255},
            pattern=LinePattern.None),
          Ellipse(
            extent={{-80,80},{80,-80}},
            lineColor={127,127,127},
            pattern=LinePattern.Dash,
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-98,80},{98,120}},
            textString="%name",
            lineColor={0,0,0})}));
  end Species;

  package BaseClasses "Base classes (generally not for direct use)"

    extends Modelica.Icons.BasesPackage;
    type Conservation = enumeration(
        IC "Initial condition imposed forever (no conservation)",
        Steady "Steady (conservation with steady state)",
        Dynamic "Dynamic (conservation with storage)")
      "Options for a conservation equation";
    type InitScalar = enumeration(
        None "No initialization",
        Amount "Prescribed amount",
        AmountSS "Steady-state amount",
        Density "Prescribed density",
        DensitySS "Steady-state density",
        Volume "Prescribed volume",
        VolumeSS "Steady-state volume",
        Pressure "Prescribed pressure",
        PressureSS "Steady-state pressure",
        Temperature "Prescribed temperature",
        TemperatureSS "Steady-state temperature",
        SpecificEnthalpy "Prescribed specific enthalpy",
        SpecificEnthalpySS "Steady-state specific enthalpy",
        PotentialGibbs "Prescribed Gibbs potential",
        PotentialGibbsSS "Steady-state Gibbs potential")
      "Methods of initializing scalar quantities (material and energy)";

    type InitTranslational = enumeration(
        None "No initialization",
        Velocity "Prescribed velocity",
        VelocitySS "Steady-state velocity",
        Current "Prescribed advective current",
        CurrentSS "Steady-state advective current")
      "Methods of initializing translational momentum";

  end BaseClasses;

end Species;

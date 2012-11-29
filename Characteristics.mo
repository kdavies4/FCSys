within FCSys;
package Characteristics "Data and functions to correlate physical properties"

  extends Modelica.Icons.MaterialPropertiesPackage;

  import Modelica.Media.IdealGases.Common.FluidData;
  import Modelica.Media.IdealGases.Common.SingleGasesData;

  package Examples "Examples and tests"
    extends Modelica.Icons.ExamplesPackage;

    model TestCorrelations "Test the property correlations"
      extends Modelica.Icons.Example;

      import FCSys.Characteristics;

      // Data
      FCSys.Characteristics.C.Graphite DataC;
      FCSys.Characteristics.C19HF37O5S.Solid DataC19HF37O5S;
      FCSys.Characteristics.'e-'.Graphite 'Datae-';
      FCSys.Characteristics.H2.Gas DataH2;
      FCSys.Characteristics.H2.Gas DataH2IG(b_v=[1], specVolPow={-1,0})
        "H2 as ideal gas";
      FCSys.Characteristics.H2O.Gas DataH2O;
      FCSys.Characteristics.H2O.Liquid DataH2OLiquid;
      FCSys.Characteristics.'H+'.Solid 'DataH+';
      FCSys.Characteristics.N2.Gas DataN2;
      FCSys.Characteristics.O2.Gas DataO2;

      parameter Q.PotentialAbsolute epsilon_tol=1e-10*U.V
        "Tolerance in potential";

      // Conditions
      Q.PressureAbsolute p=1*U.atm + time*10*U.bar/10;
      Q.TemperatureAbsolute T=273.15*U.K + 0*200*time*U.K/10;

      // Results
      // Volumic amount
      output Q.VolumeSpecific v_C=DataC.v_pT(p, T);
      output Q.VolumeSpecific v_C19HF37O5S=DataC19HF37O5S.v_pT(p, T);
      output Q.AmountVolumic 'rho_e-'='Datae-'.v_pT(p, T);
      output Q.VolumeSpecific v_H2=DataH2.v_pT(p, T);
      output Q.VolumeSpecific v_IG=DataH2IG.v_pT(p, T)
        "Volumic amount of ideal gas";
      output Q.VolumeSpecific v_H2O=DataH2O.v_pT(p, T);
      output Q.VolumeSpecific v_H2O_liquid=DataH2OLiquid.v_pT(p, T);
      output Q.AmountVolumic 'rho_H+'='DataH+'.v_pT(p, T);
      output Q.VolumeSpecific v_N2=DataN2.v_pT(p, T);
      output Q.VolumeSpecific v_O2=DataO2.v_pT(p, T);
      // Pressure
      output Q.PressureAbsolute p_C=DataC.p_vT(v_C, T) if DataC.isCompressible;
      output Q.PressureAbsolute p_C19HF37O5S=DataC19HF37O5S.p_vT(v_C19HF37O5S,
          T) if DataC19HF37O5S.isCompressible;
      output Q.PressureAbsolute 'p_e-'='Datae-'.p_vT('v_e-', T) if 'Datae-'.isCompressible;
      output Q.PressureAbsolute p_H2=DataH2.p_vT(v_H2, T) if DataH2.isCompressible;
      output Q.PressureAbsolute p_IG=DataH2IG.p_vT(v_IG, T) if DataH2IG.isCompressible
        "Pressure of ideal gas";
      output Q.PressureAbsolute p_H2O=DataH2O.p_vT(v_H2O, T) if DataH2O.isCompressible;
      output Q.PressureAbsolute p_H2O_liquid=DataH2OLiquid.p_vT(v_H2O_liquid, T)
        if DataH2OLiquid.isCompressible;
      output Q.PressureAbsolute 'p_H+'='DataH+'.p_vT('v_H+', T) if 'DataH+'.isCompressible;
      output Q.PressureAbsolute p_N2=DataN2.p_vT(v_N2, T) if DataN2.isCompressible;
      output Q.PressureAbsolute p_O2=DataO2.p_vT(v_O2, T) if DataO2.isCompressible;
      // Specific enthalpy
      output Q.Potential h_C=DataC.h0_T(T, referenceEnthalpy=ReferenceEnthalpy.EnthalpyOfFormationAt25degC);

      output Q.Potential h_C19HF37O5S=DataC19HF37O5S.h0_T(T, referenceEnthalpy=
          ReferenceEnthalpy.EnthalpyOfFormationAt25degC);

      output Q.Potential 'h_e-'='Datae-'.h0_T(T, referenceEnthalpy=
          ReferenceEnthalpy.EnthalpyOfFormationAt25degC);

      output Q.Potential h_H2=DataH2.h0_T(T, referenceEnthalpy=
          ReferenceEnthalpy.EnthalpyOfFormationAt25degC);

      output Q.Potential h_H2O=DataH2O.h0_T(T, referenceEnthalpy=
          ReferenceEnthalpy.EnthalpyOfFormationAt25degC);

      output Q.Potential h_H2O_liquid=DataH2OLiquid.h0_T(T, referenceEnthalpy=
          ReferenceEnthalpy.EnthalpyOfFormationAt25degC);

      output Q.Potential 'h_H+'='DataH+'.h0_T(T, referenceEnthalpy=
          ReferenceEnthalpy.EnthalpyOfFormationAt25degC);

      output Q.Potential h_N2=DataN2.h0_T(T, referenceEnthalpy=
          ReferenceEnthalpy.EnthalpyOfFormationAt25degC);

      output Q.Potential h_O2=DataO2.h0_T(T, referenceEnthalpy=
          ReferenceEnthalpy.EnthalpyOfFormationAt25degC);

      // Specific entropy
      output Q.Number s_C=DataC.s_pT(p, T);
      output Q.Number s_C19HF37O5S=DataC19HF37O5S.s_pT(p, T);
      output Q.Number 's_e-'='Datae-'.s_pT(p, T);
      output Q.Number s_H2=DataH2.s_pT(p, T);
      output Q.Number s_H2O=DataH2O.s_pT(p, T);
      output Q.Number s_H2O_liquid=DataH2OLiquid.s_pT(p, T);
      output Q.Number 's_H+'='DataH+'.s_pT(p, T);
      output Q.Number s_N2=DataN2.s_pT(p, T);
      output Q.Number s_O2=DataO2.s_pT(p, T);
      // Gibbs potential
      output Q.Potential g_C=DataC.g_pT(
              p,
              T,
              referenceEnthalpy=ReferenceEnthalpy.EnthalpyOfFormationAt25degC);
      output Q.Potential g_C19HF37O5S=DataC19HF37O5S.g_pT(
              p,
              T,
              referenceEnthalpy=ReferenceEnthalpy.EnthalpyOfFormationAt25degC);
      output Q.Potential 'g_e-'='Datae-'.g_pT(
              p,
              T,
              referenceEnthalpy=ReferenceEnthalpy.EnthalpyOfFormationAt25degC);
      output Q.Potential g_H2=DataH2.g_pT(
              p,
              T,
              referenceEnthalpy=ReferenceEnthalpy.EnthalpyOfFormationAt25degC);
      output Q.Potential g_H2O=DataH2O.g_pT(
              p,
              T,
              referenceEnthalpy=ReferenceEnthalpy.EnthalpyOfFormationAt25degC);
      output Q.Potential g_H2O_liquid=DataH2OLiquid.g_pT(
              p,
              T,
              referenceEnthalpy=ReferenceEnthalpy.EnthalpyOfFormationAt25degC);
      output Q.Potential 'g_H+'='DataH+'.g_pT(
              p,
              T,
              referenceEnthalpy=ReferenceEnthalpy.EnthalpyOfFormationAt25degC);
      output Q.Potential g_N2=DataN2.g_pT(
              p,
              T,
              referenceEnthalpy=ReferenceEnthalpy.EnthalpyOfFormationAt25degC);
      output Q.Potential g_O2=DataO2.g_pT(
              p,
              T,
              referenceEnthalpy=ReferenceEnthalpy.EnthalpyOfFormationAt25degC);
      // Gibbs potential (indirectly calculated)
      output Q.Potential g_C_indirect=h_C - T*s_C;
      output Q.Potential g_C19HF37O5S_indirect=h_C19HF37O5S - T*s_C19HF37O5S;
      output Q.Potential 'g_e-_indirect'='h_e-' - T*'s_e-';
      output Q.Potential g_H2_indirect=h_H2 - T*s_H2;
      output Q.Potential g_H2O_indirect=h_H2O - T*s_H2O;
      output Q.Potential g_H2O_liquid_indirect=h_H2O_liquid - T*s_H2O_liquid;
      output Q.Potential 'g_H+_indirect'='h_H+' - T*'s_H+';
      output Q.Potential g_N2_indirect=h_N2 - T*s_N2;
      output Q.Potential g_O2_indirect=h_O2 - T*s_O2;
      // Resistivity
      output Q.Resistivity gamma=DataC.gamma(T);
      output Q.Resistivity gamma_SH2_Phi=DataH2.gamma_Phi(T);
      output Q.Resistivity gamma_SH2_U=DataH2.gamma_S(T);
      output Q.Resistivity gamma_SH2O_Phi=DataH2O.gamma_Phi(T);
      output Q.Resistivity gamma_SH2O_U=DataH2O.gamma_S(T);
      output Q.Resistivity gamma_N2_Phi=DataN2.gamma_Phi(T);
      output Q.Resistivity gamma_N2_U=DataN2.gamma_S(T);
      output Q.Resistivity gamma_SO2_Phi=DataO2.gamma_Phi(T);
      output Q.Resistivity gamma_SO2_U=DataO2.gamma_S(T);

    equation
      assert(abs(g_C - g_C_indirect) < epsilon_tol,
        "The Gibbs potential for C is not consistent with its specific enthalpy and specific entropy.");
      assert(abs(g_C19HF37O5S - g_C19HF37O5S_indirect) < epsilon_tol,
        "The Gibbs potential for C19HF37O5S is not consistent with its specific enthalpy and specific entropy.");
      assert(abs('g_e-' - 'g_e-_indirect') < epsilon_tol,
        "The Gibbs potential for 'e-' is not consistent with its specific enthalpy and specific entropy.");
      assert(abs(g_H2 - g_H2_indirect) < epsilon_tol,
        "The Gibbs potential for H2 is not consistent with its specific enthalpy and specific entropy.");
      assert(abs(g_H2O - g_H2O_indirect) < epsilon_tol,
        "The Gibbs potential for H2O is not consistent with its specific enthalpy and specific entropy.");
      assert(abs(g_H2O_liquid - g_H2O_liquid_indirect) < epsilon_tol,
        "The Gibbs potential for H2O liquid is not consistent with its specific enthalpy and specific entropy.");
      assert(abs('g_H+' - 'g_H+_indirect') < epsilon_tol,
        "The Gibbs potential for 'H+' is not consistent with its specific enthalpy and specific entropy.");
      assert(abs(g_N2 - g_N2_indirect) < epsilon_tol,
        "The Gibbs potential for N2 is not consistent with its specific enthalpy and specific entropy.");
      assert(abs(g_O2 - g_O2_indirect) < epsilon_tol,
        "The Gibbs potential for O2 is not consistent with its specific enthalpy and specific entropy.");

      annotation (
        experiment(StopTime=10),
        experimentSetupOutput,
        Commands(file=
              "resources/scripts/Dymola/Characteristics.Examples.TestCorrelations.mos"));
    end TestCorrelations;
  end Examples;

  package C "C"
    extends Modelica.Icons.MaterialPropertiesPackage;
    record Graphite "C graphite"

      extends BaseClasses.Characteristic(
        final formula="C",
        final phase="graphite",
        p0=U.atm,
        specVolPow={0,0},
        b_v=[U.cm^3*m/(2.2210*U.g)],
        p_min=-Modelica.Constants.inf,
        m=12.0107*U.g/U.mol,
        Deltah0_f=0*U.J/U.mol,
        Deltah0=1053.500*U.J/U.mol,
        T_lim_c0={200.000,600.000,2000.000,6000.000}*U.K,
        b_c0=[1.132856760e5, -1.980421677e3, 1.365384188e1, -4.636096440e-2,
            1.021333011e-4, -1.082893179e-7, 4.472258860e-11; 3.356004410e5, -2.596528368e3,
            6.948841910, -3.484836090e-3, 1.844192445e-6, -5.055205960e-10,
            5.750639010e-14; 2.023105106e5, -1.138235908e3, 3.700279500, -1.833807727e-4,
            6.343683250e-8, -7.068589480e-12, 3.335435980e-16] .* fill({U.K^(3
             - i) for i in 1:7}, size(T_lim_c0, 1) - 1),
        B_c0=[8.943859760e3*U.K, -7.295824740e1; 1.398412456e4*U.K, -4.477183040e1;
            5.848134850e3, -2.350925275e1] - b_c0[:, 2:3]*ln(U.K),
        r=170*U.pico*U.m/U.q);

      annotation (Documentation(info="<html>

     <p>Assumptions:
     <ol>
     <li>Constant specific volume (i.e., incompressible and without
          thermal expansion)</li>
     </ol>
     </p>

     <p>Additional notes:
     <ul>
     <li>The radius is from <a href=\"http://en.wikipedia.org/wiki/Carbon\">http://en.wikipedia.org/wiki/Carbon</a>.  See also
   <a href=\"http://en.wikipedia.org/wiki/Van_der_Waals_radius\">http://en.wikipedia.org/wiki/Van_der_Waals_radius</a>.</li>
     <li>The default specific volume (<code>v=U.cm^3*m/(2.210*U.g)</code>) is of pyrolytic graphite
  at 300 K according to Incropera and DeWitt [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, p. 909].  Other forms
  are (also at 300 K and based on the same reference) are:
  <ul>
       <li>Amorphous carbon:  <code>v=U.cm^3*m/(1.950*U.g)</code></li>
       <li>Diamond (type IIa):  <code>v=U.cm^3*m/(3.500*U.g)</code></li>
       </li>
     </ul>
     </p>
<p>For more information, see the
  <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic\">Characteristic</a> record.</p></html>"));
    end Graphite;
  end C;

  package C19HF37O5S "<html>C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S</html>"
    extends Modelica.Icons.MaterialPropertiesPackage;
    record Solid
      "<html>C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S solid</html>"

      extends BaseClasses.Characteristic(
        final formula="C9HF17O5S",
        final phase="solid",
        p0=U.atm,
        m=1044.214*U.g/U.mol,
        specVolPow={0,0},
        b_v=[U.cm^3*m/(2.00*U.g)],
        p_min=-Modelica.Constants.inf,
        Deltah0_f=0,
        Deltah0=0,
        specHeatCapPow=0,
        T_lim_c0={0,Modelica.Constants.inf},
        b_c0=[2890*U.J*m/(U.kg*U.K)],
        B_c0=[-298.15*U.K*b_c0[1, 1] + Deltah0_f, 0],
        r=(147 + 2259.8/2)*U.pico*U.m/U.q);

      annotation (Documentation(info="<html>
       <p>Assumptions:
     <ol>
     <li>Constant specific volume (i.e., incompressible and without
          thermal expansion)</li>
     </ol>
     </p>

     <p>Additional notes:
     <ul>
     <li>A form of C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S is
     C<sub>7</sub>HF<sub>13</sub>O<sub>5</sub>S.(C<sub>2</sub>F<sub>4</sub>)<sub>6</sub>, which is a typical
   configuration of Nafion sulfonate resin after hydrolysis
   [<a href=\"modelica://FCSys.UsersGuide.References\">Mark1999</a>, p. 234].</li>
     <li>Thermodynamic data for this material is not available from McBride et
   al. [<a href=\"modelica://FCSys.UsersGuide.References\">McBride2002</a>].  The default specific heat capacity (<code>b_c0=[2890*U.J*m/(U.kg*U.K)]</code>) is
   based on data of paraffin from Incropera and DeWitt [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, p. 916].</li>
     <li>According to [<a href=\"modelica://FCSys.UsersGuide.References\">Avogadro1.03</a>], the furthest distance
   between two atoms of C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S is 2259.8 pm and is between fluorines.
   The radius of F is 147 pm (<a href=\"http://en.wikipedia.org/wiki/Fluorine\">http://en.wikipedia.org/wiki/Fluorine</a>).</li>
     <li>From [<a href=\"modelica://FCSys.UsersGuide.References\">Avogadro1.03</a>], the molecular weight of C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S
   is 1044.214 g/mol.  However, the \"11\" in Nafion 11x indicates a
   molecular weight of 1100 g/mol.  According to
   <a href=\"http://en.wikipedia.org/wiki/Nafion\">http://en.wikipedia.org/wiki/Nafion</a>,
       \"the molecular weight of Nafion is uncertain due to differences in
        processing and solution morphology.\"</li>
     <li>The specific volume (<code>v = U.cm^3*m/(2.00*U.g)</code>) is based on
   [<a href=\"modelica://FCSys.UsersGuide.References\">Lin2006</a>, p. A1327].
       </li>
     </ul>
     </p>

<p>For more information, see the
  <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic\">Characteristic</a> record.</p></html>"),
          Icon(graphics));
    end Solid;
  end C19HF37O5S;

  package 'e-' "<html>e<sup>-</sup></html>"
    extends Modelica.Icons.MaterialPropertiesPackage;
    record Gas "<html>e<sup>-</sup> gas</html>"
      import Data = Modelica.Media.IdealGases.Common.SingleGasesData.eminus;

      extends BaseClasses.Characteristic(
        final formula="e-",
        phase="gas",
        m=Data.MM*U.kg/U.mol,
        b_v=[1],
        specVolPow={-1,0},
        Deltah0_f=Data.MM*Data.Hf*U.J/U.mol,
        Deltah0=Data.MM*Data.H0*U.J/U.mol,
        T_lim_c0={200.000,20000.000}*U.K,
        b_c0={Data.alow} .* fill({U.K^(3 - i) for i in 1:size(Data.alow, 1)},
            size(T_lim_c0, 1) - 1),
        B_c0={Data.blow} .* fill({U.K,1}, size(T_lim_c0, 1) - 1) - b_c0[:, 2:3]
            *ln(U.K),
        r=U.k_A/m);

      annotation (Documentation(info="<html>
     <p>Notes:
     <ul>
     <li>The equation for the radius is the classical radius of an electron based on
  <a href=\"http://en.wikipedia.org/wiki/Classical_electron_radius\">http://en.wikipedia.org/wiki/Classical_electron_radius</a>.</li>
  <li>McBride and Gordon [<a href=\"modelica://FCSys.UsersGuide.References\">McBride1996</a>] provide correlations for the transport
  properties of e<sup>-</sup> gas.  However, they are not entered here, since they
  contain only one temperature range (2000 to 5000 K) which is beyond the expected operating range of the model.</li>
  <li>
  The thermodynamic data [<a href=\"modelica://FCSys.UsersGuide.References\">McBride2002</a>] splits the correlations into three
  temperature ranges, but the coefficients are the same for each.
   Therefore, the temperature limits are set here such that the entire
   range is handled without switching.  The lower temperature limit
   in the source data is 298.150 K, but here it is expanded down to 200 K.
   The constants are independent of temperature anyway.
  </li>
</ul>
</p>
<p>For more information, see the
  <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic\">Characteristic</a> record.</p></html>"));
    end Gas;

    record Graphite "<html>e<sup>-</sup> in graphite</html>"
      extends Gas(
        final phase="graphite",
        specVolPow=C.Graphite.specVolPow,
        b_v=C.Graphite.b_v,
        p_min=-Modelica.Constants.inf);
      annotation (Documentation(info="<html>
     <p>Assumptions:
     <ol>
     <li>There is one free (conduction) electron per atom of carbon.</li>
</ol>
</p>
<p>For more information, see the
  <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic\">Characteristic</a> record.</p></html>"));
    end Graphite;
  end 'e-';

  package H2 "<html>H<sub>2</sub></html>"
    extends Modelica.Icons.MaterialPropertiesPackage;
    record Gas "<html>H<sub>2</sub> gas</html>"
      import Data = Modelica.Media.IdealGases.Common.SingleGasesData.H2;

      extends BaseClasses.CharacteristicNASA(
        final formula=Data.name,
        final phase="gas",
        final m=Data.MM*U.kg/U.mol,
        specVolPow={-1,-3},
        b_v={{0,0,0,1},{8.0282e6*U.K^3,-2.6988e5*U.K^2,-129.26*U.K,17.472}*U.cm
            ^3/U.mol},
        p_min=1e-323*p0,
        Deltah0_f=Data.MM*Data.Hf*U.J/U.mol,
        Deltah0=Data.MM*Data.H0*U.J/U.mol,
        T_lim_c0={200.000,Data.Tlimit,6000.000,20000.000}*U.K,
        b_c0={Data.alow,Data.ahigh,{4.966884120e8,-3.147547149e5,7.984121880e1,
            -8.414789210e-3,4.753248350e-7,-1.371873492e-11,1.605461756e-16}}
             .* fill({U.K^(3 - i) for i in 1:size(Data.alow, 1)}, size(T_lim_c0,
            1) - 1),
        B_c0={Data.blow,Data.bhigh,{2.488433516e6,-6.695728110e2}} .* fill({U.K,
            1}, size(T_lim_c0, 1) - 1) - b_c0[:, 2:3]*ln(U.K),
        r=(120 + 100.3/2)*U.pico*U.m/U.q,
        T_lim_beta={200.0,1000.0,5000.0,15000.0}*U.K,
        b_eta={{0.74553182,43.555109*U.K,-0.32579340e4*U.K^2,0.13556243},{
            0.96730605,679.31897*U.K,-0.21025179e6*U.K^2,-1.8251697},{1.0126129,
            0.14973739e4*U.K,-0.14428484e7*U.K^2,-2.3254928}},
        b_lambda={{1.0059461,279.51262*U.K,-0.29792018e5*U.K^2,1.1996252},{
            1.0582450,248.75372*U.K,0.11736907e5*U.K^2,0.82758695},{-0.22364420,
            -0.69650442e4*U.K,-0.77771313e5*U.K^2,13.189369}});

      // Note:  In Dymola 7.4, ln(1e-323) returns a valid result, but ln(1e-324)
      // doesn't.

      annotation (Documentation(info="<html>
            <p>Notes:
     <ul>
  <li>According to [<a href=\"modelica://FCSys.UsersGuide.References\">Avogadro1.03</a>], the (center-to-center)
   bond length of H-H is 100.3 pm.  The radius of H is from
   <a href=\"http://en.wikipedia.org/wiki/Hydrogen\">http://en.wikipedia.org/wiki/Hydrogen</a>.  See also
   <a href=\"http://en.wikipedia.org/wiki/Van_der_Waals_radius\">http://en.wikipedia.org/wiki/Van_der_Waals_radius</a>.</li>
  <li>The virial coefficients are from [<a href=\"modelica://FCSys.UsersGuide.References\">Dymond2002</a>, p. 41].  The
  temperature range of the coefficients is [60, 500] K, but this is not enforced in the functions.</li>
     </ul>
     </p>

<p>For more information, see the
  <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic\">Characteristic</a> record.</p></html>"));
    end Gas;
  end H2;

  package H2O "<html>H<sub>2</sub>O</html>"
    extends Modelica.Icons.MaterialPropertiesPackage;
    record Gas "<html>H<sub>2</sub>O gas</html>"
      import Data = Modelica.Media.IdealGases.Common.SingleGasesData.H2O;

      extends BaseClasses.CharacteristicNASA(
        final formula=Data.name,
        final phase="gas",
        final m=Data.MM*U.kg/U.mol,
        specVolPow={-1,-3},
        b_v={{0,0,0,1},{-5.6932e10*U.K^3,1.818e8*U.K^2,-3.0107e5*U.K,158.83}*U.cm
            ^3/U.mol},
        p_min=1e-323*p0,
        Deltah0_f=Data.MM*Data.Hf*U.J/U.mol,
        Deltah0=Data.MM*Data.H0*U.J/U.mol,
        T_lim_c0={200.000,Data.Tlimit,6000.000}*U.K,
        b_c0={Data.alow,Data.ahigh} .* fill({U.K^(3 - i) for i in 1:size(Data.alow,
            1)}, size(T_lim_c0, 1) - 1),
        B_c0={Data.blow,Data.bhigh} .* fill({U.K,1}, size(T_lim_c0, 1) - 1) -
            b_c0[:, 2:3]*ln(U.K),
        r=(282/2)*U.pico*U.m/U.q,
        T_lim_beta={373.2,1073.2,5000.0,15000.0}*U.K,
        b_eta={{0.50019557,-697.12796*U.K,0.88163892e5*U.K^2,3.0836508},{
            0.58988538,-537.69814*U.K,0.54263513e5*U.K^2,2.3386375},{0.64330087,
            -95.668913*U.K,-0.37742283e6*U.K^2,1.8125190}},
        b_lambda={{1.0966389,-555.13429*U.K,0.10623408e6*U.K^2,-0.24664550},{
            0.39367933,-0.22524226e4*U.K,0.61217458e6*U.K^2,5.8011317},{-0.41858737,
            -0.14096649e5*U.K,0.19179190e8*U.K^2,14.345613}});
      // Note:  In Dymola 7.4, ln(1e-323) returns a valid result, but ln(1e-324)
      // doesn't.

      annotation (Documentation(info="<html>
        <p>Notes:
     <ul>
  <li>The radius of H<sub>2</sub>O is 282 pm
   (<a href=\"http://www.lsbu.ac.uk/water/molecule.html\">http://www.lsbu.ac.uk/water/molecule.html</a>).  Using the radius of H
   from <a href=\"http://en.wikipedia.org/wiki/Hydrogen\">http://en.wikipedia.org/wiki/Hydrogen</a> and the center-to-center
   distance of hydrogen atoms in H<sub>2</sub>O from [<a href=\"modelica://FCSys.UsersGuide.References\">Avogadro1.03</a>],
   156.6 pm, the radius of H<sub>2</sub>O would be (120 + 156.6/2) pm = 198.3 pm.</li>
  <li>The virial coefficients are from [<a href=\"modelica://FCSys.UsersGuide.References\">Dymond2002</a>, p. 4]).  The
  temperature range of the coefficients is [350, 770] K, but this is not enforced in the functions.</li>
     </ul>
     </p>

  <p>For more information, see the
  <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic\">Characteristic</a> record.</p></html>"));
    end Gas;

    record Liquid "<html>H<sub>2</sub>O liquid</html>"

      extends BaseClasses.Characteristic(
        final formula="H2O",
        final phase="liquid",
        final m=0.01801528*U.kg/U.mol,
        p0=U.atm,
        specVolPow={0,0},
        b_v=[U.cm^3*m/(0.99656*U.g)],
        Deltah0_f=-285830.000*U.J/U.mol,
        Deltah0=13278.000*U.J/U.mol,
        T_lim_c0={273.150,373.150,600.000}*U.K,
        b_c0=[1.326371304e9, -2.448295388e7, 1.879428776e5, -7.678995050e2,
            1.761556813, -2.151167128e-3, 1.092570813e-6; 1.263631001e9, -1.680380249e7,
            9.278234790e4, -2.722373950e2, 4.479243760e-1, -3.919397430e-4,
            1.425743266e-7] .* fill({U.K^(3 - i) for i in 1:7}, size(T_lim_c0,
            1) - 1),
        B_c0=[1.101760476e8*U.K, -9.779700970e5; 8.113176880e7*U.K, -5.134418080e5]
             - b_c0[:, 2:3]*ln(U.K),
        r=(282/2)*U.pico*U.m/U.q);

      annotation (Documentation(info="<html>     <p>Assumptions:
     <ol>
     <li>Constant specific volume (i.e., incompressible and without
          thermal expansion)</li>
     </ol>
     </p>
            <p>Additional notes:
     <ul>
     <li>See note in <a href=\"modelica://FCSys.Characteristics.H2O.Gas\">Characteristics.H2O.Gas</a> regarding the radius.</li>
     <li>The default specific volume (<code>b_v=[U.cm^3*m/(0.99656*U.g)]</code>) is at 300 K based on [<a href=\"modelica://FCSys.UsersGuide.References\">Takenaka1990</a>].</li>
     </ul>
     </p>
  <p>For more information, see the
  <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic\">Characteristic</a> record.</p></html>"));
    end Liquid;
  end H2O;

  package 'H+' "<html>H<sup>+</sup></html>"
    extends Modelica.Icons.MaterialPropertiesPackage;
    record Gas "<html>H<sup>+</sup> gas</html>"
      import Data = Modelica.Media.IdealGases.Common.SingleGasesData.Hplus;

      extends BaseClasses.Characteristic(
        final formula="H+",
        phase="gas",
        final m=Data.MM*U.kg/U.mol,
        b_v=[1],
        specVolPow={-1,0},
        Deltah0_f=Data.MM*Data.Hf*U.J/U.mol,
        Deltah0=Data.MM*Data.H0*U.J/U.mol,
        specHeatCapPow=-2,
        T_lim_c0={200.000,20000.000}*U.K,
        b_c0={Data.alow} .* fill({U.K^(3 - i) for i in 1:size(Data.alow, 1)},
            size(T_lim_c0, 1) - 1),
        B_c0={Data.blow} .* fill({U.K,1}, size(T_lim_c0, 1) - 1) - b_c0[:, 2:3]
            *ln(U.K),
        r=120*U.pico*U.m/U.q);

      annotation (Documentation(info="<html>
         <p>Assumptions:
     <ol>
  <li>The radius (for the purpose of the rigid-sphere assumption of
  kinetic theory) is that of H
  (<a href=\"http://en.wikipedia.org/wiki/Hydrogen\">http://en.wikipedia.org/wiki/Hydrogen</a>).</li>
     </ol>
     </p>

     <p>Additional notes:
     <ul>
     <li>The source data [<a href=\"modelica://FCSys.UsersGuide.References\">McBride2002</a>] breaks the data into three
   temperature ranges, but the constants are the same for each.
   Therefore, the temperature limits are set here such that the entire
   range is handled without switching.
   The lower temperature limit in the source data
   is 298.150 K, but here it is expanded down to 200 K.  The constants are
   independent of temperature anyway.</li>
   <li>The enthalpy to produce H<sup>+</sup> from H<sub>2</sub>O (H<sub>2</sub>O &#8640; OH<sup>-</sup> + H<sup>+</sup>) is
   -96569.804 J/mol.  The enthalpy to produce H<sup>+</sup> from H<sub>3</sub>O<sup>+</sup>
   (H<sub>3</sub>O<sup>+</sup> &#8640; H<sub>2</sub>O + H<sup>+</sup>) is 839826.0 J/mol.  Based on
   Tissandier et al., the
   enthalpy of formation of aqueous H<sup>+</sup> is 1150.1e3 J/mol [<a href=\"modelica://FCSys.UsersGuide.References\">Tissandier1998</a>].</li>
     </ul>
     </p>
  <p>For more information, see the
  <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic\">Characteristic</a> record.</p></html>"));
    end Gas;

    record Solid "<html>H<sup>+</sup> in solid</html>"
      extends Gas(
        final phase="solid",
        specVolPow={0,0},
        b_v=[1/(0.95*U.M)],
        p_min=-Modelica.Constants.inf,
        Deltah0_f=H2O.Gas.Deltah0_f/4,
        h_offset=-'H+'.Gas.Deltah0_f + H2O.Gas.Deltah0_f/4);

      annotation (Documentation(info="<html><p>The enthalpy of formation (<code>Deltah0_f</code>) and the additional enthalpy
  offset (<code>h_offset</code>) are specified such that the enthalpy in the
  <a href=\"modelica://FCSys.Subregiosn.BaseClasses.PartialSpecies\">PartialSpecies</a> model is referenced to zero at 0 &deg;C.
  Otherwise, there is a very large bias that causes problems with the reaction.  The proper
  values are not known from the NASA CEA thermodynamic data
    [<a href=\"modelica://FCSys.UsersGuide.References\">McBride2002</a>] (it only provides data for H<sup>+</sup> gas).
  Since the same H<sup>+</sup> model is used in the anode and cathode, there is no net effect.</p>

  <p>The specific volume of protons corresponds to the concentration measured
  by Spry and Fayer (0.95 M) in Nafion<sup>&reg;</sup> at
  &lambda; = 12, where &lambda; is the number of
  H<sub>2</sub>O molecules to SO<sub>3</sub>H
  endgroups.  At &lambda; = 22, the concentration was measured at 0.54 M
  [<a href=\"modelica://FCSys.UsersGuide.References\">Spry2009</a>].</p>
</html>"));
    end Solid;
  end 'H+';

  package N2 "<html>N<sub>2</sub></html>"
    extends Modelica.Icons.MaterialPropertiesPackage;
    record Gas "<html>N<sub>2</sub> gas</html>"
      import Data = Modelica.Media.IdealGases.Common.SingleGasesData.N2;

      extends BaseClasses.CharacteristicNASA(
        final formula=Data.name,
        final phase="gas",
        final m=Data.MM*U.kg/U.mol,
        specVolPow={-1,-4},
        b_v={{0,0,0,0,1},{2.7198e9*U.K^4,6.1253e7*U.K^3,-1.4164e6*U.K^2,-9337.8
            *U.K,40.286}*U.cm^3/U.mol},
        p_min=1e-323*p0,
        Deltah0_f=Data.MM*Data.Hf*U.J/U.mol,
        Deltah0=Data.MM*Data.H0*U.J/U.mol,
        T_lim_c0={200.000,Data.Tlimit,6000.000,20000.000}*U.K,
        b_c0={Data.alow,Data.ahigh,{8.310139160e8,-6.420733540e5,2.020264635e2,
            -3.065092046e-2,2.486903333e-6,-9.705954110e-11,1.437538881e-15}}
             .* fill({U.K^(3 - i) for i in 1:size(Data.alow, 1)}, size(T_lim_c0,
            1) - 1),
        B_c0={Data.blow,Data.bhigh,{4.938707040e6,-1.672099740e3}} .* fill({U.K,
            1}, size(T_lim_c0, 1) - 1) - b_c0[:, 2:3]*ln(U.K),
        r=(155 + 145.2/2)*U.pico*U.m/U.q,
        T_lim_beta={200.0,1000.0,5000.0,15000.0}*U.K,
        b_eta={{0.62526577,-31.779652*U.K,-0.16407983e4*U.K^2,1.7454992},{
            0.87395209,561.52222*U.K,-0.17394809e6*U.K^2,-0.39335958},{
            0.88503551,909.02171*U.K,-0.73129061e6*U.K^2,-0.53503838}},
        b_lambda={{0.85439436,105.73224*U.K,-0.12347848e5*U.K^2,0.47793128},{
            0.88407146,133.57293*U.K,-0.11429640e5*U.K^2,0.24417019},{2.4176185,
            0.80477749e4*U.K,0.31055802e7*U.K^2,-14.517761}});
      // Note:  In Dymola 7.4, ln(1e-323) returns a valid result, but ln(1e-324)
      // doesn't.

      annotation (Documentation(info="<html>
                  <p>Notes:
     <ul>
  <li>According to [<a href=\"modelica://FCSys.UsersGuide.References\">Avogadro1.03</a>], the (center-to-center)
   bond length of N-N is 145.2 pm.  The radius of N is from
   <a href=\"http://en.wikipedia.org/wiki/Nitrogen\">http://en.wikipedia.org/wiki/Nitrogen</a>.  See also
   <a href=\"http://en.wikipedia.org/wiki/Van_der_Waals_radius\">http://en.wikipedia.org/wiki/Van_der_Waals_radius</a>.</li>
  <li>The virial coefficients are from Dymond et al. (<a href=\"modelica://FCSys.UsersGuide.References\">2002</a>, p. 69).  The
  temperature range of the coefficients is [75, 745] K, but this is not enforced in the functions.  More precise virial coefficients are available from
  <a href=\"http://www.tpub.com/content/nasa1996/NASA-96-cr4755/NASA-96-cr47550059.htm\">http://www.tpub.com/content/nasa1996/NASA-96-cr4755/NASA-96-cr47550059.htm</a>.</li>
     </ul>
     </p>
  <p>For more information, see the
  <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic\">Characteristic</a> record.</p></html>"));
    end Gas;
  end N2;

  package O2 "<html>O<sub>2</sub></html>"
    extends Modelica.Icons.MaterialPropertiesPackage;
    record Gas "<html>O<sub>2</sub> gas</html>"
      import Data = Modelica.Media.IdealGases.Common.SingleGasesData.O2;

      extends BaseClasses.CharacteristicNASA(
        final formula=Data.name,
        final phase="gas",
        final m=Data.MM*U.kg/U.mol,
        specVolPow={-1,-4},
        b_v={{0,0,0,0,1},{5.0855e9*U.K^4,-1.6393e8*U.K^3,5.2007e5*U.K^2,-1.7696e4
            *U.K,42.859}*U.cm^3/U.mol},
        p_min=1e-323*p0,
        Deltah0_f=Data.MM*Data.Hf*U.J/U.mol,
        Deltah0=Data.MM*Data.H0*U.J/U.mol,
        T_lim_c0={200.000,Data.Tlimit,6000.000,20000.000}*U.K,
        b_c0={Data.alow,Data.ahigh,{4.975294300e8,-2.866106874e5,6.690352250e1,
            -6.169959020e-3,3.016396027e-7,-7.421416600e-12,7.278175770e-17}}
             .* fill({U.K^(3 - i) for i in 1:size(Data.alow, 1)}, size(T_lim_c0,
            1) - 1),
        B_c0={Data.blow,Data.bhigh,{2.293554027e6,-5.530621610e2}} .* fill({U.K,
            1}, size(T_lim_c0, 1) - 1) - b_c0[:, 2:3]*ln(U.K),
        r=(152 + 128.2/2)*U.pico*U.m/U.q,
        T_lim_beta={200.0,1000.0,5000.0,15000.0}*U.K,
        b_eta={{0.60916180,-52.244847*U.K,-599.74009*U.K^2,2.0410801},{
            0.72216486,175.50839*U.K,-0.57974816e5*U.K^2,1.0901044},{0.73981127,
            391.94906*U.K,-0.37833168e6*U.K^2,0.90931780}},
        b_lambda={{0.77229167,6.8463210*U.K,-0.58933377e4*U.K^2,1.2210365},{
            0.90917351,291.24182*U.K,-0.79650171e5*U.K^2,0.064851631},{-1.1218262,
            -0.19286378e5*U.K,0.23295011e8*U.K^2,20.342043}});
      // Note:  In Dymola 7.4, ln(1e-323) returns a valid result, but ln(1e-324)
      // doesn't.

      annotation (Documentation(info="<html><p>Notes:<ul>
  <li>According to [<a href=\"modelica://FCSys.UsersGuide.References\">Avogadro1.03</a>], the (center-to-center)
   bond length of O-O is 128.2 pm.  The radius of O is from
   <a href=\"http://en.wikipedia.org/wiki/Oxygen\">http://en.wikipedia.org/wiki/Oxygen</a>.  See also
   <a href=\"http://en.wikipedia.org/wiki/Van_der_Waals_radius\">http://en.wikipedia.org/wiki/Van_der_Waals_radius</a>.</li>
  <li>The virial coefficients are from Dymond et al. (<a href=\"modelica://FCSys.UsersGuide.References\">2002</a>, p. 69).  The
  temperature range of the coefficients is [70, 495] K, but this is not enforced in the functions.</li>
     </ul>
     </p>
<p>For more information, see the
  <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic\">Characteristic</a> record.</p></html>"));
    end Gas;
  end O2;

  package BaseClasses "Base classes (not for direct use)"
    extends Modelica.Icons.BasesPackage;
    record CharacteristicNASA
      "Thermodynamic record with transport properties based on NASA CEA"

      extends FCSys.Characteristics.BaseClasses.Characteristic;

      constant Q.TemperatureAbsolute T_lim_beta[:]={0,Modelica.Constants.inf}
        "Temperature limits for the rows of b_eta and b_lambda";
      constant Real b_eta[size(T_lim_beta, 1) - 1, 4]
        "Constants in the NASA CEA correlation for viscosity";
      constant Real b_lambda[size(T_lim_beta, 1) - 1, 4]
        "Constants in the NASA CEA correlation for thermal conductivity";

      function gamma_Phi
        "<html>Resistivity to transverse transport of linear momentum as a function of temperature (&gamma;<sub>&Phi;</sub>)</html>"

        extends Modelica.Icons.Function;

        input Q.TemperatureAbsolute T "Temperature";
        output Q.Resistivity gamma_Phi "Transverse resistivity";

      protected
        function b_eta_adj
          "Return unit-adjusted NASA CEA constants for viscosity"
          // See note in b_eta_adj() in gamma_Phi().
          output Real b_eta_adj[size(T_lim_beta, 1) - 1, 4]
            "Unit-adjusted NASA CEA constants for viscosity";
        algorithm
          b_eta_adj := transpose({b_eta[:, 1],b_eta[:, 2],b_eta[:, 3],b_eta[:,
            4] - b_eta[:, 1]*ln(U.K) + fill(ln(1e-6*U.g/(U.cm*U.s*m)), size(
            T_lim_beta, 1) - 1)}) annotation (Inline=true);
        end b_eta_adj;

      algorithm
        /*
  assert(T_lim_beta[1] <= T and T <= T_lim_beta[size(T_lim_beta, 1)], "Temperature "
     + String(T/(U.K)) + " K is out of range for the diffusion properties ([" +
    String(T_lim_beta[1]/U.K) + ", " + String(T_lim_beta[size(T_lim_beta, 1)]
    /U.K) + "] K).");
  */
        // Note:  This is commented out so that the function can be inlined.
        // Note:  In Dymola 7.4, T_lim_beta[end] can't be used instead of
        // T_lim_beta[size(T_lim_beta, 1)] due to:
        //     "Error, not all "end" could be expanded."

        gamma_Phi := smooth(0, exp(-sum(if (T_lim_beta[i] <= T or i == 1) and (
          T <= T_lim_beta[i + 1] or i == size(T_lim_beta, 1) - 1) then (
          b_eta_adj())[i, 1]*ln(T) + ((b_eta_adj())[i, 2] + (b_eta_adj())[i, 3]
          /T)/T + (b_eta_adj())[i, 4] else 0 for i in 1:size(T_lim_beta, 1) - 1)))
          annotation (Inline=true, smoothOrder=2);
        // Note:  The annotation is set assuming that the values of the constants
        // result in a function that is first-order continuous.

        annotation (Documentation(info="<html><p>This function is based on based on NASA CEA
  [<a href=\"modelica://FCSys.UsersGuide.References\">McBride1996</a>, <a href=\"modelica://FCSys.UsersGuide.References\">Svehla1995</a>]</p></html>"));
      end gamma_Phi;

      function gamma_S
        "<html>Thermal transport resistivity as a function of temperature (&gamma;<sub><i>S</i></sub>)</html>"

        extends Modelica.Icons.Function;

        input Q.TemperatureAbsolute T "Temperature";
        output Q.Resistivity gamma_S "Thermal resistivity";

      protected
        function b_lambda_adj
          "Return unit-adjusted NASA CEA constants for thermal conductivity"
          // Note:  If b_lambda were defined as a local constant instead of a
          // function, it would prevent the gamma_Phi function from
          // being inlined.  If it were defined as a global final constant, it
          // would need to be updated manually when b_eta is changed.
          output Real b_lambda_adj[size(T_lim_beta, 1) - 1, 4]
            "Unit-adjusted NASA CEA constants for thermal conductivity";
        algorithm
          b_lambda_adj := transpose({b_lambda[:, 1],b_lambda[:, 2],b_lambda[:,
            3],b_lambda[:, 4] - b_lambda[:, 1]*ln(U.K) + fill(ln(1e-6*U.W/(U.cm
            *U.K)), size(T_lim_beta, 1) - 1)}) annotation (Inline=true);
        end b_lambda_adj;

      algorithm
        /*
  assert(T_lim_beta[1] <= T and T <= T_lim_beta[size(T_lim_beta, 1)], "Temperature "
     + String(T/(U.K)) + " K is out of range for the diffusion properties ([" +
    String(T_lim_beta[1]/U.K) + ", " + String(T_lim_beta[size(T_lim_beta, 1)]
    /U.K) + "] K).");
  */
        // Note:  This is commented out so that the function can be inlined.
        // Note:  In Dymola 7.4, T_lim_beta[end] can't be used instead of
        // T_lim_beta[size(T_lim_beta, 1)] due to:
        //     "Error, not all "end" could be expanded."

        gamma_S := smooth(0, exp(-sum(if (T_lim_beta[i] <= T or i == 1) and (T
           <= T_lim_beta[i + 1] or i == size(T_lim_beta, 1) - 1) then (
          b_lambda_adj())[i, 1]*ln(T) + ((b_lambda_adj())[i, 2] + (b_lambda_adj())
          [i, 3]/T)/T + (b_lambda_adj())[i, 4] else 0 for i in 1:size(
          T_lim_beta, 1) - 1))) annotation (Inline=true, smoothOrder=2);
        // Note:  The annotation is set assuming that the values of the constants
        // result in a function that is first-order continuous.
        annotation (Documentation(info="<html><p>This function is based on based on NASA CEA
  [<a href=\"modelica://FCSys.UsersGuide.References\">McBride1996</a>, <a href=\"modelica://FCSys.UsersGuide.References\">Svehla1995</a>]</p></html>"));
      end gamma_S;

      annotation (defaultComponentPrefixes="replaceable",Documentation(info="<html><p>The correlations for transport properties are available in
  [<a href=\"modelica://FCSys.UsersGuide.References\">McBride1996</a>,
  <a href=\"modelica://FCSys.UsersGuide.References\">McBride2002</a>]. For more information, see the
  <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic\">Characteristic</a> record.</p></html>"));
    end CharacteristicNASA;

    record Characteristic "Record for thermodynamic and resistive properties"
      extends FCSys.BaseClasses.Icons.Record;

      constant String formula "Chemical formula";
      constant String phase "Phase (e.g., \"gas\")";
      constant Q.MassSpecific m "Specific mass";
      constant Q.LengthSpecific r "Specific radius";
      final constant Integer z=Chemistry.charge(formula) "Charge number";
      constant Q.PressureAbsolute p0=U.bar "Reference pressure";
      constant Q.Pressure p_min=-Modelica.Constants.inf
        "Minimum pressure for numerical protection";
      constant Real b_v[:, :]=[1]
        "Coefficients of specific volume as a polynomial in p/T and T"
        annotation (HideResult=false);
      constant Integer specVolPow[2]={-1,0}
        "Powers of p/T and T for 1st row and column of b_v, respectively"
        annotation (HideResult=false);
      final constant Boolean isCompressible=sum(countTrue({b_v[i, j] <> 0 and
          specVolPow[1] + i - 1 <> 0 for i in 1:size1b_v}) for j in 1:size2b_v)
           > 0 "true, if specific volume is a function of pressure";
      final constant Boolean hasThermalExpansion=sum(countTrue({b_v[i, j] <> 0
           and specVolPow[2] + j - specVolPow[1] - i <> 0 for i in 1:size1b_v})
          for j in 1:size2b_v) > 0
        "true, if specific volume is a function of temperature";
      constant Q.PotentialChemical Deltah0_f
        "<html>Enthalpy of formation at 298.15 K, <i>p</i>&deg; (&Delta;<i>h</i>&deg;<sub>f</sub>)</html>";
      constant Q.PotentialChemical Deltah0
        "<html><i>h</i>&deg;(298.15 K) - <i>h</i>&deg;(0 K) (&Delta;<i>h</i>&deg;)</html>";
      constant Q.PotentialChemical h_offset=0 "Additional enthalpy offset";
      constant Integer specHeatCapPow=-2 "Power of T for 1st column of b_c0"
        annotation (HideResult=false);
      constant Q.TemperatureAbsolute T_lim_c0[:]={0,Modelica.Constants.inf}
        "Temperature limits for the rows of b_c and B_c";
      constant Real b_c0[size(T_lim_c0, 1) - 1, :]
        "<html>Coefficients of specific heat capacity at <i>p</i>&deg; as a polynomial in <i>T</i> (<i>b</i><sub><i>c</i>&deg;</sub>)</html>"
        annotation (HideResult=false);
      constant Real B_c0[size(T_lim_c0, 1) - 1, 2]
        "<html>Integration constants for specific enthalpy and entropy (<i>B</i><sub><i>c</i>&deg;</sub>)</html>"
        annotation (HideResult=false);

    protected
      final constant Integer pressPow[2]={specVolPow[1] - size1b_v + 1,
          specVolPow[2] + 1}
        "Powers of v and T for 1st row and column of b_p, respectively";
      final constant Real b_p[size1b_v, size2b_v]={(if specVolPow[1] + i == 0
           or specVolPow[1] + i == 1 then b_v[i, :] else (if specVolPow[1] + i
           == 2 and specVolPow[1] <= 0 then b_v[i, :] + b_v[i - 1, :] .^ 2
           else (if specVolPow[1] + i == 3 and specVolPow[1] <= 0 then b_v[i, :]
           + b_v[i - 2, :] .* (b_v[i - 2, :] .^ 2 + 3*b_v[i - 1, :]) else zeros(
          size2b_v)))) for i in size1b_v:-1:1}
        "Coefficients of p as a polynomial in v and T";
      // Note:  pressPow and b_p are only used in h_pT(); however, in Dymola 7.4
      // they must be defined here (global to h_pT()) so that they are updated
      // properly when the size of b_v is changed from its default.
      final constant Integer size1b_v=size(b_v, 1) "Number of rows in b_v";
      final constant Integer size2b_v=size(b_v, 2) "Number of columns in b_v";

    public
      partial function gamma
        "<html>Ideal transport resistivity as a function of temperature (&gamma;)</html>"
        extends Modelica.Icons.Function;

        input Q.TemperatureAbsolute T=298.15*U.K "Temperature";
        output Q.Resistivity gamma "Resistivity";

      algorithm
        gamma := U.pi*r^2*U.q*sqrt(U.pi*m/T)*6 annotation (
          Inline=true,
          smoothOrder=999,
          Documentation(info="<html>
  <p>This function is based on kinetic theory of gases with the rigid-sphere (\"billiard-ball\")
  assumption [<a href=\"modelica://FCSys.UsersGuide.References\">Present1958</a>].</p></html>"));
      end gamma;

      function c0_T
        "<html>Specific heat capacity at reference pressure as a function of temperature (<i>c</i>&deg;)</html>"
        extends Modelica.Icons.Function;

        input Q.TemperatureAbsolute T=298.15*U.K "Temperature";
        output Q.CapacityThermalSpecific c0
          "Specific heat capacity at reference pressure";

      algorithm
        /*
  assert(T_lim_c0[1] <= T and T <= T_lim_c0[size(T_lim_c0, 1)], "Temperature " +
    String(T/U.K) + " K is out of range for " + name + " ([" + String(T_lim_c0[1]
    /U.K) + ", " + String(T_lim_c0[size(T_lim_c0, 1)]/U.K) + "] K).");
  */
        // Note:  This is commented out so that the function can be inlined.
        // Note:  In Dymola 7.4, T_lim_c0[end] can't be used instead of
        // T_lim_c0[size(T_lim_c0, 1)] due to:
        //    "Error, not all 'end' could be expanded."

        c0 := smooth(0, sum(if (T_lim_c0[i] <= T or i == 1) and (T < T_lim_c0[i
           + 1] or i == size(T_lim_c0, 1) - 1) then poly(
                T,
                b_c0[i, :],
                specHeatCapPow) else 0 for i in 1:size(T_lim_c0, 1) - 1))
          annotation (
          InlineNoEvent=true,
          Inline=true,
          smoothOrder=0);

      end c0_T;

      function dp
        "<html>Derivative of pressure as it is defined by <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.p_vT\">p_vT</a>()</html>"
        extends Modelica.Icons.Function;

        input Q.VolumeSpecific v=U.atm/(298.15*U.K) "Specific volume";
        input Q.TemperatureAbsolute T=298.15*U.K "Temperature";
        input Q.VolumeSpecific dv=0 "Derivative of specific volume";
        input Q.Temperature dT=0 "Derivative of temperature";
        output Q.Pressure dp "Derivative of pressure";

      algorithm
        dp := if isCompressible then poly(
                v,
                {(pressPow[1] + i - 1)*poly(
                  T,
                  b_p[i, :],
                  pressPow[2]) for i in 1:size1b_v},
                pressPow[1] - 1)*dv + poly(
                T,
                {(pressPow[2] + i - 1)*poly(
                  v,
                  b_p[:, i],
                  pressPow[1]) for i in 1:size2b_v},
                pressPow[2] - 1)*dT else 0
          annotation (Inline=true, smoothOrder=999);
      end dp;

      function g_pT "Gibbs potential as a function of pressure and temperature"

        extends Modelica.Icons.Function;

        input Q.PressureAbsolute p=1*U.atm "Pressure";
        input Q.TemperatureAbsolute T=298.15*U.K "Temperature";
        input ReferenceEnthalpy referenceEnthalpy=ReferenceEnthalpy.EnthalpyOfFormationAt25degC
          "Choice of enthalpy reference";
        output Q.Potential g "Gibbs potential";

      protected
        function g0_i
          "Return g0 as a function of T using one of the temperature ranges, with enthalpy of formation at 25 degC"
          input Q.TemperatureAbsolute T "Temperature";
          input Integer i "index of the temperature range";
          output Q.Potential g0_i "g0";

        algorithm
          g0_i := poly(
                    T,
                    {b_c0[i, j]*((if specHeatCapPow + j == 0 then ln(T) else 1/
              (specHeatCapPow + j)) - (if specHeatCapPow + j == 1 then ln(T)
               else 1/(specHeatCapPow + j - 1))) - (if specHeatCapPow + j == 1
               then B_c0[i, 2] else 0) for j in 1:size(b_c0, 2)},
                    specHeatCapPow + 1) + B_c0[i, 1] annotation (Inline=true);
        end g0_i;

      algorithm
        /*
    assert(T_lim_c0[1] <= T and T <= T_lim_c0[size(T_lim_c0, 1)], "Temperature " +
    String(T/U.K) + " K is out of range for " + name + " ([" + String(T_lim_c0[1]
    /U.K) + ", " + String(T_lim_c0[size(T_lim_c0, 1)]/U.K) + "] K).");
    */
        // Note:  This is commented out so that the function can be inlined.
        // Note:  In Dymola 7.4, T_lim_c0[end] can't be used instead of
        // T_lim_c0[size(T_lim_c0, 1)] due to:
        //    "Error, not all 'end' could be expanded."

        g := smooth(1, sum(if (T_lim_c0[i] <= T or i == 1) and (T < T_lim_c0[i
           + 1] or i == size(T_lim_c0, 1) - 1) then g0_i(T, i) else 0 for i in
          1:size(T_lim_c0, 1) - 1) + (if referenceEnthalpy == ReferenceEnthalpy.ZeroAt0K
           then Deltah0 else 0) - (if referenceEnthalpy == ReferenceEnthalpy.ZeroAt25degC
           then Deltah0_f else 0) + h_offset + sum((if specVolPow[1] + i == 0
           then ln((if p_min > 0 then max(p, p_min) else p)/p0) else (p^(
          specVolPow[1] + i) - p0^(specVolPow[1] + i))/(specVolPow[1] + i))*
          poly( T,
                b_v[i, :],
                specVolPow[2] - specVolPow[1] - i + 1) for i in 1:size1b_v))
          annotation (
          InlineNoEvent=true,
          Inline=true,
          smoothOrder=1);
        // The first term is the integral of c0*dT up to T with the reference
        // enthalpy at the lower bound [McBride2002, p. 2] plus T times the
        // integral of (c0/T)*dT up to T with absolute entropy at the lower bound.
        // Both of these integrals are taken at p0.  The second polynomial is the
        // integral of v*dP from p0 to p (at T).
      end g_pT;

      function h0_T
        "<html>Specific enthalpy as a function temperature at reference pressure (<i>h</i>&deg;)</html>"

        extends Modelica.Icons.Function;

        input Q.TemperatureAbsolute T=298.15*U.K "Temperature";
        input ReferenceEnthalpy referenceEnthalpy=ReferenceEnthalpy.EnthalpyOfFormationAt25degC
          "Choice of enthalpy reference";
        output Q.Potential h0 "Specific enthalpy at reference pressure";

      protected
        function h0_i
          "Return h0 as a function of T using one of the temperature ranges, with enthalpy of formation at 25 degC"
          input Q.TemperatureAbsolute T "Temperature";
          input Integer i "index of the temperature range";
          output Q.Potential h0_i "h0";

        algorithm
          h0_i := poly(
                    T,
                    b_c0[i, :] .* {if specHeatCapPow + j == 0 then ln(T) else 1
              /(specHeatCapPow + j) for j in 1:size(b_c0, 2)},
                    specHeatCapPow + 1) + B_c0[i, 1] annotation (Inline=true);
        end h0_i;

      algorithm
        /*
    assert(T_lim_c0[1] <= T and T <= T_lim_c0[size(T_lim_c0, 1)], "Temperature " +
    String(T/U.K) + " K is out of range for " + name + " ([" + String(T_lim_c0[1]
    /U.K) + ", " + String(T_lim_c0[size(T_lim_c0, 1)]/U.K) + "] K).");
    */
        // Note:  This is commented out so that the function can be inlined.
        // Note:  In Dymola 7.4, T_lim_c0[end] can't be used instead of
        // T_lim_c0[size(T_lim_c0, 1)] due to:
        //    "Error, not all 'end' could be expanded."

        h0 := smooth(1, sum(if (T_lim_c0[i] <= T or i == 1) and (T < T_lim_c0[i
           + 1] or i == size(T_lim_c0, 1) - 1) then h0_i(T, i) else 0 for i in
          1:size(T_lim_c0, 1) - 1) + (if referenceEnthalpy == ReferenceEnthalpy.ZeroAt0K
           then Deltah0 else 0) - (if referenceEnthalpy == ReferenceEnthalpy.ZeroAt25degC
           then Deltah0_f else 0) + h_offset)
          annotation (
          InlineNoEvent=true,
          Inline=true,
          smoothOrder=1);
        // The first term is the integral of c0*dT up to T at p0 with the
        // reference enthalpy at the lower bound [McBride2002, p. 2].
      end h0_T;

      function p_vT "Pressure as a function of specific volume and temperature"
        extends Modelica.Icons.Function;

        input Q.VolumeSpecific v=U.atm/(298.15*U.K) "Specific volume";
        input Q.TemperatureAbsolute T=298.15*U.K "Temperature";
        output Q.PressureAbsolute p "Pressure";

      algorithm
        // assert(isCompressible,
        //  "The pressure is undefined since the material is incompressible.",
        //  AssertionLevel.warning);
        // Note:  In Dymola 7.4, the assertion level can't be set, although it has
        // been defined as an argument to assert() since Modelica 3.0.

        p := if isCompressible then poly(
                v,
                {poly(
                  T,
                  b_p[i, :],
                  pressPow[2]) for i in 1:size1b_v},
                pressPow[1]) else 0 annotation (
          Inline=true,
          smoothOrder=999,
          inverse(v=v_pT(p, T)),
          derivative=dp);
        annotation (Documentation(info="<html><p>If the species is incompressible, then <i>p</i>(<i>v</i>, <i>T</i> ) is undefined,
  and the function will return a value of zero.</p>
  <p>The derivative of this function is <a href=\"modelica://FCSys.Characteristics.BaseClasses.Characteristic.dp\">dp</a>().</p></html>"));
      end p_vT;

      function s_pT
        "Specific entropy as a function of pressure and temperature"
        extends Modelica.Icons.Function;

        input Q.PressureAbsolute p=1*U.atm "Pressure";
        input Q.TemperatureAbsolute T=298.15*U.K "Temperature";
        output Q.NumberAbsolute s "Specific entropy";

      protected
        function s0_i
          "Return s0 as a function of T using one of the temperature ranges"
          input Q.TemperatureAbsolute T "Temperature";
          input Integer i "index of the temperature range";
          output Q.NumberAbsolute s0_i "s0";

        algorithm
          s0_i := poly(
                    T,
                    b_c0[i, :] .* {if specHeatCapPow + j == 1 then ln(T) else 1
              /(specHeatCapPow + j - 1) for j in 1:size(b_c0, 2)},
                    specHeatCapPow) + B_c0[i, 2] annotation (Inline=true);
        end s0_i;

      algorithm
        /*
  assert(T_lim_c0[1] <= T and T <= T_lim_c0[size(T_lim_c0, 1)], "Temperature " +
    String(T/U.K) + " K is out of range for " + name + " ([" + String(T_lim_c0[1]
    /U.K) + ", " + String(T_lim_c0[size(T_lim_c0, 1)]/U.K) + "] K).");
  */
        // Note:  This is commented out so that the function can be inlined.
        // Note:  In Dymola 7.4, T_lim_c0[end] can't be used instead of
        // T_lim_c0[size(T_lim_c0, 1)] due to:
        //    "Error, not all 'end' could be expanded."

        s := smooth(1, sum(if (T_lim_c0[i] <= T or i == 1) and (T < T_lim_c0[i
           + 1] or i == size(T_lim_c0, 1) - 1) then s0_i(T, i) else 0 for i in
          1:size(T_lim_c0, 1) - 1) - sum((if specVolPow[1] + i == 0 then ln((
          if p_min > 0 then max(p, p_min) else p)/p0) else (p^(specVolPow[1] +
          i) - p0^(specVolPow[1] + i))/(specVolPow[1] + i))*poly(
                T,
                b_v[i, :],
                specVolPow[2] - specVolPow[1] - i) for i in 1:size1b_v))
          annotation (
          InlineNoEvent=true,
          Inline=true,
          smoothOrder=1);
        // The first term is the integral of c0/T*dT up to T at p0 with the
        // absolute entropy at the lower bound [McBride2002, p. 2].  The second
        // polynomial is the integral of v*dP from p0 to p (at T).
      end s_pT;

      function v_pT "Specific volume as a function of pressure and temperature"
        extends Modelica.Icons.Function;

        input Q.PressureAbsolute p=1*U.atm "Pressure";
        input Q.TemperatureAbsolute T=298.15*U.K "Temperature";
        output Q.VolumeSpecific v "Specific volume";

      algorithm
        v := poly(
                p/T,
                {poly(
                  T,
                  b_v[i, :],
                  specVolPow[2]) for i in 1:size1b_v},
                specVolPow[1]) annotation (
          Inline=true,
          smoothOrder=999,
          inverse(p=p_vT(v, T)));
      end v_pT;

      annotation (defaultComponentPrefixes="replaceable",Documentation(info="<html>
    <p>This package is compatible with NASA CEA thermodynamic data
    [<a href=\"modelica://FCSys.UsersGuide.References\">McBride2002</a>] and the generalized virial equation of state
    [<a href=\"modelica://FCSys.UsersGuide.References\">Dymond2002</a>].  It is compatible with&mdash;but not limited to&mdash;the
    assumptions of ideal gas and constant specific volume.</p>

    <p>Assumptions:
    <ol><li>Specific mass is constant.</li></ol>
    </p>
    <p>The following notes apply to the constants:
    <ul>
    <li>Currently, <code>formula</code> may not contain parentheses or brackets.</li>
    <li><code>r</code> is the Van der Waals radius or the radius for the
    rigid-sphere (\"billiard-ball\") approximation of the kinetic theory of gases.</li>
    <li><code>b_v</code>: The powers of <i>p</i>/<i>T</i> increase by row.  The powers of
    <i>T</i> increase by column.  If <code>specVolPow[1] = -1</code>, then the rows
    of <code>b_v</code> correspond to 1, <i>T</i><i>B</i><sup>*</sup>(<i>T</i> ), <i>T</i><sup> 2</sup><i>C</i><sup>*</sup>(<i>T</i> ), <i>T</i><sup> 3</sup><i>D</i><sup>*</sup>(<i>T</i> ), &hellip;
    in [<a href=\"modelica://FCSys.UsersGuide.References\">Dymond2002</a>].</li>
    <li>The defaults for <code>b_v</code> and <code>specVolPow</code> represent ideal gas.</li>
    <li><code>b_c</code>: The rows give the coefficients for different temperature ranges&mdash;bounded
    by the values in <code>T_lim_c0</code>.
    The powers of <i>T</i> increase
    by column.
    By default,
    the powers of <i>T</i> for the 1st column are each -2, which corresponds to [<a href=\"modelica://FCSys.UsersGuide.References\">McBride2002</a>].
    In that case, the dimensionalities of the coefficients are {M2.L4/(N2.T4), M.L2/(N.T2), 1, &hellip;}
    for each row, where L is length, M is mass, N is number, and T is time. (In <a href=\"modelica://FCSys\">FCSys</a>,
    temperature is a potential with dimension M.L2/(N.T2); see
    the <a href=\"modelica://FCSys.Units\">Units</a> package.)</li>
    <li><code>B_c</code>: As in <code>b_c</code>, the rows correspond to different
    temperature ranges.  The first column is for specific enthalpy and has dimensionality
    M.L2/(N.T2).  The second column is for specific entropy and is dimensionless.
    The integration constants for enthalpy are defined such that the enthalpy at
    25 &deg;C is the specific enthalpy of formation at that temperature [<a href=\"modelica://FCSys.UsersGuide.References\">McBride2002</a>, p. 2].
    The integration constants for specific entropy are defined such that specific entropy is absolute.</li>
    <li><code>T_lim_c0</code>: The first and last entries are the minimum and
    maximum valid temperatures.  The intermediate entries are the thresholds
    between rows of <code>b_c</code> (and <code>B_c</code>).  Therefore, if there are <i>n</i> temperature ranges
    (and rows in <code>b_c</code> and <code>B_c</code>), then <code>T_lim_c0</code> must
    have <i>n</i> + 1 entries.</li>
    <li><code>p0</code> is the reference pressure.   In the
    NASA CEA data [<a href=\"modelica://FCSys.UsersGuide.References\">McBride2002</a>], it is 1 bar for gases or 1 atm for condensed
    species.</li>
    <li>If the <i>p</i>-<i>v</i>-<i>T</i> equation of state includes an ideal gas term, then the
    correlation for specific entropy (<i>s</i>) will involve the natural logarithm of partial pressure.
    The <code>p_min</code> parameter is used to guard against the logarithm of a non-positive pressure.
    To disable this protection (and simplify the translated code), set <code>p_min</code> to zero or a
    negative pressure.</li>
    </ul>
    </p></html>"));
    end Characteristic;

    type ReferenceEnthalpy = enumeration(
        ZeroAt0K "Enthalpy at 0 K is 0 (if no additional offset)",
        ZeroAt25degC "Enthalpy at 25 degC is 0 (if no additional offset)",
        EnthalpyOfFormationAt25degC
          "Enthalpy at 25 degC is the enthalpy of formation at 25 degC (if no additional offset)")
      "Enumeration defining the reference enthalpy of a species";
  end BaseClasses;
  annotation (Documentation(info="<html>
  <p>Additional materials may be included as needed.  The thermodynamic data for
  materials that are condensed at standard conditions is available in
  [<a href=\"modelica://FCSys.UsersGuide.References\">McBride2002</a>].
  The thermodynamic data for materials
  that are gases at standard conditions is available in
  <a href=\"modelica://Modelica.Media.IdealGases.Common.SingleGasesData\">Modelica.Media.IdealGases.Common.SingleGasesData</a>
  (and [<a href=\"modelica://FCSys.UsersGuide.References\">McBride2002</a>]). Virial coefficients are available in
  [<a href=\"modelica://FCSys.UsersGuide.References\">Dymond2002</a>].  Transport characteristics are available in
  [<a href=\"modelica://FCSys.UsersGuide.References\">McBride1996</a>,
  <a href=\"modelica://FCSys.UsersGuide.References\">McBride2002</a>].
  </p>
  <p>
<b>Licensed by the Georgia Tech Research Corporation under the Modelica License 2</b><br>
Copyright 2007&ndash;2012, Georgia Tech Research Corporation.
</p>
<p>
<i>This Modelica package is <u>free</u> software and the use is completely at <u>your own risk</u>;
it can be redistributed and/or modified under the terms of the Modelica License 2. For license conditions (including the
disclaimer of warranty) see <a href=\"modelica://FCSys.UsersGuide.ModelicaLicense2\">
FCSys.UsersGuide.ModelicaLicense2</a> or visit <a href=\"http://www.modelica.org/licenses/ModelicaLicense2\">
http://www.modelica.org/licenses/ModelicaLicense2</a>.</i>
</p></html>"));
end Characteristics;

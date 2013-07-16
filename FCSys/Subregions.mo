within FCSys;
package Subregions
  "Control volumes with multi-species transport, exchange, and storage"
  extends Modelica.Icons.Package;
  package Examples "Examples"
    model AirColumn
      "<html>Test a one-dimensional array of subregions with an initial pressure difference (C<sup>+</sup> and N<sub>2</sub> included by default)</html>"
      import FCSys.BaseClasses.Utilities.round;
      extends Modelica.Icons.Example;

      parameter Integer n_y=3
        "<html>Number of discrete subregions along the y axis, besides the 2 boundary subregions (<i>n<i><sub>y</sub>)</html>";
      parameter Q.Pressure Deltap_IC=0
        "<html>Initial pressure difference (&Delta;<i>p</i><sub>IC</sub>)</html>";

      output Q.Pressure Deltap=subregion2.gas.N2.p - subregion1.gas.N2.p
        "Measured pressure difference";
      output Q.Pressure Deltap_ex=-(n_y + 1)*10*U.m*environment.a[Axis.y]*
          Characteristics.N2.Gas.m*subregions[round(n_y/2)].gas.N2.rho
        "Expected pressure difference";

      inner Conditions.Environment environment(analysis=true, p=U.atm - U.kPa)
        annotation (Placement(transformation(extent={{30,60},{50,80}})));
      FCSys.Subregions.Subregion subregion1(
        L={10,10,10}*U.m,
        inclFacesZ=false,
        inclTransZ=false,
        inclTransX=false,
        inclTransY=true,
        inclFacesX=false,
        inclFacesY=true,
        gas(
          inclH2O=false,
          inclN2=true,
          N2(
            p_IC=environment.p - Deltap_IC/2,
            consTransY=Conservation.IC,
            upstreamY=false),
          k_E=1e-8),
        graphite('inclC+'=true, 'C+'(V_IC=U.mm^3)))
        annotation (Placement(transformation(extent={{-10,-40},{10,-20}})));

      FCSys.Subregions.Subregion subregions[n_y](
        each L={10,10,10}*U.m,
        each inclTransZ=false,
        each inclFacesZ=false,
        gas(
          each inclH2O=false,
          each inclN2=true,
          N2(p_IC={environment.p - Deltap_IC/2 - i*Deltap_IC/(n_y + 1) for i
                 in 1:n_y}, each upstreamY=false),
          each k_E=1e-8),
        graphite(each 'inclC+'=true, each 'C+'(V_IC=U.mm^3)),
        each inclTransX=false,
        each inclTransY=true,
        each inclFacesX=false,
        each inclFacesY=true) if n_y > 0
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

      FCSys.Subregions.Subregion subregion2(
        L={10,10,10}*U.m,
        inclFacesZ=false,
        inclTransZ=false,
        inclTransX=false,
        inclTransY=true,
        inclFacesX=false,
        inclFacesY=true,
        gas(
          inclH2O=false,
          inclN2=true,
          N2(
            p_IC=environment.p + Deltap_IC/2,
            consTransY=Conservation.IC,
            upstreamY=false),
          k_E=1e-8),
        graphite('inclC+'=true, 'C+'(V_IC=U.mm^3)))
        annotation (Placement(transformation(extent={{-10,20},{10,40}})));

    equation
      connect(subregion1.yPositive, subregions[1].yNegative) annotation (Line(
          points={{6.10623e-16,-20},{0,-18},{1.22125e-15,-16},{6.10623e-16,-16},
              {6.10623e-16,-10}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      for i in 1:n_y - 1 loop
        connect(subregions[1:n_y - 1].yPositive, subregions[2:n_y].yNegative)
          "Not shown in the diagram";
      end for;
      if n_y == 0 then
        connect(subregion1.yPositive, subregion2.yNegative)
          "Not shown in the diagram";
      end if;
      connect(subregions[n_y].yPositive, subregion2.yNegative) annotation (Line(
          points={{6.10623e-16,10},{0,14},{1.22125e-15,16},{6.10623e-16,16},{
              6.10623e-16,20}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      annotation (
        Placement(transformation(extent={{70,70},{90,90}})),
        experiment(
          StopTime=1.5,
          Tolerance=1e-06,
          Algorithm="Dassl"),
        Documentation(info="<html><p>This is a model of a vertical column of 10&nbsp;&times;&nbsp;10&nbsp;&times;&nbsp;10&nbsp;m
    regions with N<sub>2</sub> gas.  The upper and lower boundary regions are held with zero velocity.  
    The initial pressure difference is zero, but a pressure difference
    develops due to gravity.  There are some oscillations due to the pressure/translational dynamics.
    After about 1.5&nbsp;s, the pressure difference settles to 
      <i>L</i><sub>y</sub> <i>a</i><sub>y</sub> <i>m</i> &rho;
      as expected.</p>
      
      <p>A temperature gradient is created due to the thermodynamics of the expanding and contracting 
      gases.  It takes much longer (about a day) for the temperatures to equalize since the gas has a
      relatively low thermal conductivity.</p>
      
      <p>Assumptions:
      <ol>
      <li>Graphite is included as a solid, stationary species with small relative volume and very slight 
      friction to dampen the oscillations.</li>
      <li>The central difference scheme is used (no upstream discretization).</ol>
      </p></html>"),
        Commands(file=
              "Resources/Scripts/Dymola/Subregions.Examples.AirColumn.mos"
            "Subregions.Examples.AirColumn.mos"),
        experimentSetupOutput);
    end AirColumn;
    extends Modelica.Icons.ExamplesPackage;
    // TODO: In the documentation of each model, insert the sample plots or link
    // to the sample results of the User's Guide.
    // TODO: In the documentation of each model, add discussion from the dissertation.
    model Echo
      "Two regions of gas with initial pressure difference, no dampening"
      extends Subregions('inclC+'=false);
      annotation (experiment(StopTime=0.0003, Tolerance=1e-06), Commands(file(
              ensureTranslated=true) =
            "Resources/Scripts/Dymola/Subregions.Examples.Echo.mos"
            "Subregions.Examples.Echo.mos"));

    end Echo;

    model EchoCentral
      "Two regions of gas with initial pressure difference (no dampening), with central difference scheme"
      extends Subregions(
        'inclC+'=false,
        subregion1(gas(
            H2(upstreamX=false),
            H2O(upstreamX=false),
            N2(upstreamX=false),
            O2(upstreamX=false))),
        subregions(gas(
            each H2(upstreamX=false),
            each H2O(upstreamX=false),
            each N2(upstreamX=false),
            each O2(upstreamX=false))),
        subregion2(gas(
            H2(upstreamX=false),
            H2O(upstreamX=false),
            N2(upstreamX=false),
            O2(upstreamX=false))));
      annotation (experiment(StopTime=0.0003, Tolerance=1e-06), Commands(file(
              ensureTranslated=true) =
            "Resources/Scripts/Dymola/Subregions.Examples.EchoCentral.mos"
            "Subregions.Examples.EchoCentral.mos"));
    end EchoCentral;

    model ElectricalConduction
      "<html>Test a one-dimensional array of subregions with C<sup>+</sup> and e<sup>-</sup></html>"

      output Q.Potential w=(subregion.graphite.'e-'.faces[1, Side.p].mPhidot[
          Orientation.normal]/subregion.graphite.'e-'.faces[1, Side.p].rho -
          subregion.graphite.'e-'.faces[1, Side.n].mPhidot[Orientation.normal]/
          subregion.graphite.'e-'.faces[1, Side.n].rho)/subregion.A[Axis.x]
        "Potential";
      output Q.Current zI=-subregion.graphite.'e-'.Ndot_faces[1, Side.n]
        "Electrical current";
      output Q.ResistanceElectrical R=w/zI "Measured electrical resistance";
      output Q.ResistanceElectrical R_ex=subregion.graphite.'e-'.v*subregion.L[
          Axis.x]/(subregion.graphite.'e-'.mu*subregion.A[Axis.x])
        "Expected electrical resistance";
      output Q.Power P=subregion.graphite.'e-'.Edot_DT
        "Measured rate of heat generation";
      output Q.Power P_ex=zI^2*R "Expected rate of heat generation";
      output Q.TemperatureAbsolute T=subregion.graphite.'C+'.T
        "Measured temperature";
      output Q.TemperatureAbsolute T_ex=environment.T + subregion.graphite.'C+'.theta
          *U.cm*P/(4*U.mm^2) "Expected temperature";

      extends Examples.Subregion(
        'inclC+'=true,
        'incle-'=true,
        inclH2=false,
        subregion(
          L={U.cm,U.mm,U.mm},
          void=true,
          graphite(
            reduceTemp=true,
            'C+'(initMaterial=InitScalar.Pressure),
            'e-'(
              initMaterial=InitScalar.Volume,
              consMaterial=Conservation.IC,
              initTransX=InitTranslational.None,
              initEnergy=InitScalar.None,
              sigma=1e2*U.S/U.m))));

      FCSys.Conditions.ByConnector.FaceBus.Single.FaceBusFlows BC1(graphite(
          final 'incle-'='incle-',
          'inclC+'=true,
          'e-'(redeclare function normalSpec =
                Conditions.ByConnector.Face.Single.TranslationalNormal.currentDensity,
              normalSource(y=-5*U.A/U.cm^2)),
          'C+'(redeclare function thermalSpec =
                Conditions.ByConnector.Face.Single.Thermal.temperature,
              thermalSource(y=environment.T)))) annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={-24,0})));

      FCSys.Conditions.ByConnector.FaceBus.Single.FaceBusFlows BC2(graphite(
          'inclC+'=true,
          final 'incle-'='incle-',
          'C+'(redeclare function thermalSpec =
                Conditions.ByConnector.Face.Single.Thermal.temperature,
              thermalSource(y=environment.T)))) annotation (Placement(
            transformation(
            extent={{-10,10},{10,-10}},
            rotation=90,
            origin={24,0})));

    equation
      connect(BC1.face, subregion.xNegative) annotation (Line(
          points={{-20,3.65701e-16},{-16,3.65701e-16},{-16,6.10623e-16},{-10,
              6.10623e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      connect(subregion.xPositive, BC2.face) annotation (Line(
          points={{10,6.10623e-16},{16,6.10623e-16},{16,-2.54679e-16},{20,-2.54679e-16}},

          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      annotation (
        experiment(StopTime=100, Tolerance=1e-06),
        Commands(file(ensureTranslated=true) =
            "Resources/Scripts/Dymola/Subregions.Examples.ElectricalConduction.mos"
            "Subregions.Examples.ElectricalConduction.mos"),
        Documentation(info="<html><p>This is an example of Ohm's law.  The 1&nbsp;cm &times; 1&nbsp;mm &times; 1&nbsp;mm subregion contains
    carbon (C<sup>+</sup>) and electrons.  An electrical current of 50&nbsp;mA (5&nbsp;A/cm<sup>2</sup>) is delivered into
    the negative boundary and exits the positive boundary.  Due to the finite mobility of
    the electrons a force is required to support the current; this maps directly to
    electrical potential.  The example shows that the measured resistance is
    <i>R</i>&nbsp;=&nbsp;<i>L</i>/(<i>A</i>&nbsp;<i>&rho;</i>&nbsp;&mu;) as expected, where &rho; is electronic
    density and &mu; is electronic mobility.</p>

    <p>The measured rate of heat generation (<code>subregion.graphite.'e-'.Edot_DT</code>)
    is equal to <i>P</i> = (<i>zI</i>)<sup>2</sup> <i>R</i> as expected, where
    <i>zI</i> = 50&nbsp;mA is the electrical current.  This heat is conducted through the carbon
    to the boundaries, which are held at 25&nbsp;&deg;C.  The measured steady state temperature
    is <i>T</i>&nbsp;=&nbsp;<i>T</i><sub>0</sub>&nbsp;+&nbsp;&theta;&nbsp;<i>L</i>&nbsp;<i>P</i>/(4&nbsp;<i>A</i>) as expected, where
    <i>T</i><sub>0</sub>&nbsp;=&nbsp;25&nbsp;&deg;C is the boundary temperature and
    &theta; is the thermal resistance.  The factor of one fourth
    is due to the boundary conditions; the conduction length is half of the total length
    and the heat is rejected to both sides.  There is no thermal convection or 
    radiation&mdash;only conduction to the sides.</p>
    </html>"));
    end ElectricalConduction;

    model Evaporation "Test a subregion with evaporation and condensation"

      output Q.Pressure p_sat=
          Modelica.Media.Air.MoistAir.saturationPressureLiquid(subregion.gas.H2O.T
          /U.K)*U.Pa "Saturation pressure via Modelica.Media";

      extends Examples.Subregion(
        inclH2O=true,
        inclH2=false,
        subregion(gas(H2O(p_IC=U.kPa, consEnergy=Conservation.IC)), liquid(
              inclH2O=inclH2O, H2O(consEnergy=Conservation.IC))));
      FCSys.Conditions.ByConnector.FaceBus.Single.FaceBusFlows BC1(gas(inclH2O=
              true, H2O(redeclare Modelica.Blocks.Sources.Pulse materialSource(
              amplitude=-U.A,
              width=0.01,
              period=10000,
              startTime=1), redeclare function normalSpec =
                FCSys.Conditions.ByConnector.Face.Single.Translational.velocity)),
          liquid(inclH2O=true, H2O(redeclare function normalSpec =
                FCSys.Conditions.ByConnector.Face.Single.Translational.velocity)))
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={-24,0})));

      FCSys.Conditions.ByConnector.FaceBus.Single.FaceBusFlows BC2(gas(inclH2O=
              true, H2O(redeclare function normalSpec =
                FCSys.Conditions.ByConnector.Face.Single.Translational.velocity)),
          liquid(inclH2O=true, H2O(redeclare function normalSpec =
                FCSys.Conditions.ByConnector.Face.Single.Translational.velocity)))
        annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=90,
            origin={24,0})));

    equation
      connect(BC1.face, subregion.xNegative) annotation (Line(
          points={{-20,3.65701e-16},{-16,-3.36456e-22},{-16,6.10623e-16},{-10,
              6.10623e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      connect(subregion.xPositive, BC2.face) annotation (Line(
          points={{10,6.10623e-16},{16,6.10623e-16},{16,-2.54679e-16},{20,-2.54679e-16}},

          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      annotation (
        Documentation(info="<html><p>Initially, the water vapor is below saturation and liquid water is present.
  Some of the liquid evaporates until saturation is reached.  Then from 1&nbsp;to&nbsp;2&nbsp;s additional water vapor is injected from the negative boundary
  and some of it condenses.  All of this occurs at a fixed temperature (25&nbsp;&deg;C).</p></html>"),

        experiment(StopTime=3, Tolerance=1e-06),
        Commands(file(ensureTranslated=true) =
            "Resources/Scripts/Dymola/Subregions.Examples.Evaporation.mos"
            "Subregions.Examples.Evaporation.mos"),
        Diagram(graphics));
    end Evaporation;

    model HOR "Test the oxygen reduction reaction in one subregion"

      output Q.Potential w=positiveBC.ionomer.'H+'.face.mPhidot[1]/(positiveBC.ionomer.
          'H+'.face.rho*subregion.A[Axis.x]) - negativeBC.graphite.'e-'.face.mPhidot[
          1]/(negativeBC.graphite.'e-'.face.rho*subregion.A[Axis.x])
        "Electrical potential";
      /*
  output Q.Current zI=subregion.graphite.'e-'.chemical.Ndot + subregion.graphite.
      'e-'.faces[1, 2].Ndot "Electrical current due to reaction";
  output Q.Number zJ_Apercm2=zI*U.cm^2/(subregion.A[Axis.x]*U.A) 
    "Electrical current density, in A/cm2";
  output Q.Power Qdot=subregion.graphite.'C+'.Edot_DE "Rate of heat generation";
  output Q.Power P=w*zI "Electrical power";
  output Q.NumberAbsolute eta=P/(P + Qdot) "Efficiency";
  // **Fix Qdot, P, eta
  // **Update plots
  */
      extends Examples.Subregion(
        'inclC+'=true,
        'inclC19HF37O5S-'=true,
        'incle-'=true,
        'inclH+'=true,
        inclH2=true,
        inclH2O=true,
        subregion(L={0.287*U.mm,10*U.cm,10*U.cm}, gas(H2(consTransX=
                  Conservation.IC))));

      Conditions.ByConnector.FaceBus.Single.FaceBusEfforts negativeBC(gas(
          inclH2=true,
          inclH2O=true,
          H2(materialSource(y=(environment.p - environment.p_H2O)/environment.T)),

          H2O(materialSource(y=environment.p_H2O/environment.T))), graphite(
            'incle-'=true, 'e-'(redeclare function normalSpec =
                Conditions.ByConnector.Face.Single.TranslationalNormal.currentDensity,
              redeclare Modelica.Blocks.Sources.Ramp normalSource(
              height=-U.A/U.cm^2,
              duration=100.1,
              startTime=0.1)))) annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=270,
            origin={-24,0})));

      Conditions.ByConnector.FaceBus.Single.FaceBusEfforts positiveBC(ionomer(
            'inclH+'=true, 'H+'(redeclare function normalSpec =
                Conditions.ByConnector.Face.Single.TranslationalNormal.force)))
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={24,0})));

    equation
      connect(negativeBC.face, subregion.xNegative) annotation (Line(
          points={{-20,-1.34539e-15},{-16,-1.34539e-15},{-16,6.10623e-16},{-10,
              6.10623e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      connect(subregion.xPositive, positiveBC.face) annotation (Line(
          points={{10,6.10623e-16},{14,6.10623e-16},{14,1.23436e-15},{20,
              1.23436e-15}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      annotation (
        experiment(StopTime=110, Tolerance=1e-06),
        Commands(file="Resources/Scripts/Dymola/Subregions.Examples.HOR.mos"
            "Subregions.Examples.HOR.mos"),
        experimentSetupOutput);
    end HOR;

    model Hydration
      "<html>Test absorption and desorption of H<sub>2</sub>O between the gas and ionomer</html>"

      extends Examples.Subregion(
        'inclC19HF37O5S-'=true,
        inclH2O=true,
        inclH2=false,
        subregion(gas(H2O(p_IC=1.001*U.atm)), ionomer(inclH2O=inclH2O)));
      // In Dymola 7.4, p_IC=1.1*environment.p has no effect on the
      // initial pressure, but p_IC=1.1*U.atm does.
      extends Modelica.Icons.UnderConstruction;
      // TODO: Update the EOS for H2O in ionomer and recheck the result.
      annotation (experiment(StopTime=0.003, Tolerance=1e-06), Commands(file(
              ensureTranslated=true) =
            "Resources/Scripts/Dymola/Subregions.Examples.Hydration.mos"
            "Subregions.Examples.Hydration.mos"));

    end Hydration;

    model InternalFlow "Internal, laminar flow of liquid water"
      import FCSys.BaseClasses.Utilities.Delta;

      final parameter Q.Area A=subregion.A[Axis.x] "Cross-sectional area";

      // Conditions
      parameter Q.VolumeRate Vdot=U.L/U.s "Prescribed volume flow rate";

      // Measurements
      output Q.Pressure Deltap=Delta(subregion.liquid.H2O.p_faces[1, :]) - sum(
          subregion.liquid.H2O.faces[1, :].mPhidot[Orientation.normal])/A
        "Measured pressure difference (thermodynamic + nonequilibrium)";
      output Q.Length D=2*A/(subregion.L[Axis.y] + subregion.L[Axis.z]);
      output Q.Number Re=subregion.liquid.H2O.phi[Axis.x]*D*subregion.liquid.H2O.zeta
          *subregion.liquid.H2O.Data.m*subregion.liquid.H2O.rho
        "Reynolds number";
      output Q.Pressure Deltap_Poiseuille=-32*subregion.L[Axis.x]*subregion.liquid.H2O.phi[
          Axis.x]/(D^2*subregion.liquid.H2O.zeta)
        "Pressure difference according to Poiseuille's law";

      extends Examples.Subregion(inclH2=false, subregion(
          L={U.m,U.cm,U.cm},
          void=true,
          inclFacesY=true,
          inclFacesZ=true,
          liquid(inclH2O=true,H2O(
              final V_IC=subregion.V,
              final beta=0,
              initTransX=InitTranslational.None))));

      Conditions.ByConnector.FaceBus.Single.FaceBusFlows BC1(liquid(inclH2O=
              true, H2O(redeclare function normalSpec =
                FCSys.Conditions.ByConnector.Face.Single.TranslationalNormal.velocity,
              redeclare Modelica.Blocks.Sources.Sine normalSource(
              amplitude=0.05*Vdot/A,
              offset=Vdot/A,
              freqHz=0.01)))) annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=270,
            origin={-24,0})));

      Conditions.ByConnector.FaceBus.Single.FaceBusFlows BC2(liquid(inclH2O=
              true)) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={24,0})));

      Conditions.ByConnector.FaceBus.Single.FaceBusEfforts BC3(liquid(inclH2O=
              true, H2O(
            redeclare function materialSpec =
                Conditions.ByConnector.Face.Single.Material.current,
            materialSource(y=0),
            redeclare function normalSpec =
                Conditions.ByConnector.Face.Single.TranslationalNormal.force)))
        annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=0,
            origin={0,-24})));

      Conditions.ByConnector.FaceBus.Single.FaceBusEfforts BC4(liquid(inclH2O=
              true, H2O(
            redeclare function materialSpec =
                Conditions.ByConnector.Face.Single.Material.current,
            materialSource(y=0),
            redeclare function normalSpec =
                Conditions.ByConnector.Face.Single.TranslationalNormal.force)))
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={0,24})));

      Conditions.ByConnector.FaceBus.Single.FaceBusEfforts BC5(liquid(inclH2O=
              true, H2O(
            redeclare function materialSpec =
                Conditions.ByConnector.Face.Single.Material.current,
            materialSource(y=0),
            redeclare function normalSpec =
                Conditions.ByConnector.Face.Single.TranslationalNormal.force)))
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=315,
            origin={24,24})));

      Conditions.ByConnector.FaceBus.Single.FaceBusEfforts BC6(liquid(inclH2O=
              true, H2O(
            redeclare function materialSpec =
                Conditions.ByConnector.Face.Single.Material.current,
            materialSource(y=0),
            redeclare function normalSpec =
                Conditions.ByConnector.Face.Single.TranslationalNormal.force)))
        annotation (Documentation(info="<html><p>Note that the temperature increases due to viscous dissipation.  
        However, the temperature rise is limited because the walls are held at constant temperature.</p></html>"),
          Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=315,
            origin={-24,-24})));

    equation
      connect(BC1.face, subregion.xNegative) annotation (Line(
          points={{-20,-1.34539e-15},{-16,-1.34539e-15},{-16,6.10623e-16},{-10,
              6.10623e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      connect(BC2.face, subregion.xPositive) annotation (Line(
          points={{20,1.23436e-15},{16,1.23436e-15},{16,6.10623e-16},{10,
              6.10623e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      connect(subregion.yNegative, BC3.face) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-12.5},{6.10623e-16,-12.5},{
              6.10623e-16,-15},{6.10623e-16,-20},{6.10623e-16,-20}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));
      connect(BC4.face, subregion.yPositive) annotation (Line(
          points={{6.10623e-16,20},{-3.36456e-22,16},{6.10623e-16,16},{
              6.10623e-16,10}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      connect(subregion.zPositive, BC6.face) annotation (Line(
          points={{-5,-5},{-21.1716,-21.1716}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      connect(subregion.zNegative, BC5.face) annotation (Line(
          points={{5,5},{21.1716,21.1716}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));
      annotation (
        experiment(
          StopTime=800,
          Tolerance=1e-08,
          Algorithm="Dassl"),
        Commands(file(ensureTranslated=true) =
            "Resources/Scripts/Dymola/Subregions.Examples.InternalFlow.mos"
            "Subregions.Examples.InternalFlow.mos"),
        experimentSetupOutput);
    end InternalFlow;

    model ORR
      "<html>Test a subregion with the oxygen reduction reaction and essential species (C<sup>+</sup>, C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S<sup>-</sup>, e<sup>-</sup>, H<sup>+</sup>, O<sub>2</sub>, and H<sub>2</sub>O)</html>"

      output Q.Potential w=negativeBC.ionomer.'H+'.face.mPhidot[1]/(negativeBC.ionomer.
          'H+'.face.rho*subregion.A[Axis.x]) - positiveBC.graphite.'e-'.face.mPhidot[
          1]/(positiveBC.graphite.'e-'.face.rho*subregion.A[Axis.x])
        "Electrical potential";

      /* **
  output Q.Current zI=-subregion.graphite.'e-'.chemical.Ndot 
    "Electrical current due to reaction";
  output Q.Number zJ_Apercm2=zI*U.cm^2/(subregion.A[Axis.x]*U.A) 
    "Electrical current density, in A/cm2";
  output Q.Power Qdot=subregion.graphite.'C+'.Edot_DE "Rate of heat generation";
  output Q.Power P=w*zI "Electrical power";
  output Q.NumberAbsolute eta=P/(P + Qdot) "Efficiency";
*/
      // **Fix Qdot, P, eta

      extends Examples.Subregion(
        'inclC+'=true,
        'inclC19HF37O5S-'=true,
        'incle-'=true,
        'inclH+'=true,
        inclH2=false,
        inclH2O=true,
        inclO2=true,
        subregion(L={0.287*U.mm,10*U.cm,10*U.cm}, gas(
            reduceTemp=true,
            H2O(consTransX=Conservation.IC,p_IC=environment.p_H2O),
            O2(
              consTransX=Conservation.IC,
              p_IC=environment.p*environment.n_O2,
              initEnergy=InitScalar.None))));

      Conditions.ByConnector.FaceBus.Single.FaceBusEfforts negativeBC(ionomer(
            'inclH+'=true, 'H+'(redeclare function normalSpec =
                Conditions.ByConnector.Face.Single.TranslationalNormal.force)))
        annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=270,
            origin={-24,0})));

      Conditions.ByConnector.FaceBus.Single.FaceBusEfforts positiveBC(gas(
          inclH2O=true,
          inclO2=true,
          H2O(materialSource(y=environment.p_H2O/environment.T)),
          O2(materialSource(y=environment.p_O2/environment.T))), graphite(
            'incle-'=true, 'e-'(redeclare function normalSpec =
                Conditions.ByConnector.Face.Single.TranslationalNormal.currentDensity,
              redeclare Modelica.Blocks.Sources.Ramp normalSource(
              height=-U.A/U.cm^2,
              duration=100.1,
              startTime=0.1)))) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={24,0})));

    equation
      connect(negativeBC.face, subregion.xNegative) annotation (Line(
          points={{-20,-1.34539e-15},{-16,-1.34539e-15},{-16,6.10623e-16},{-10,
              6.10623e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      connect(subregion.xPositive, positiveBC.face) annotation (Line(
          points={{10,6.10623e-16},{14,6.10623e-16},{14,1.23436e-15},{20,
              1.23436e-15}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      annotation (
        experiment(StopTime=110, Tolerance=1e-06),
        Commands(file(ensureSimulated=true) =
            "Resources/Scripts/Dymola/Subregions.Examples.ORR.mos"
            "Subregions.Examples.ORR.mos"),
        experimentSetupOutput,
        Diagram(graphics));
    end ORR;

    model Reaction
      "Test an electrochemical reaction driven by a single charge carrier (positive or negative)"
      extends Modelica.Icons.Example;
      import FCSys.BaseClasses.Utilities.cartWrap;
      import FCSys.BaseClasses.Utilities.countTrue;
      import FCSys.BaseClasses.Utilities.enumerate;
      import FCSys.BaseClasses.Utilities.index;

      // Geometry
      inner parameter Q.Length L[Axis](each min=Modelica.Constants.small) = {10
        *U.cm,10*U.cm,0.1*U.mm} "Lengths";
      final inner parameter Q.Area A[Axis]={L[cartWrap(axis + 1)]*L[cartWrap(
          axis + 2)] for axis in Axis} "Cross-sectional areas";

      // Assumptions
      // -----------
      // Included components of translational momentum
      parameter Boolean inclTransX=true "X" annotation (choices(
            __Dymola_checkBox=true), Dialog(
          tab="Assumptions",
          group="Axes with translational momentum included",
          compact=true));
      parameter Boolean inclTransY=false "Y" annotation (choices(
            __Dymola_checkBox=true), Dialog(
          tab="Assumptions",
          group="Axes with translational momentum included",
          compact=true));
      parameter Boolean inclTransZ=false "Z" annotation (choices(
            __Dymola_checkBox=true), Dialog(
          tab="Assumptions",
          group="Axes with translational momentum included",
          compact=true));

      // Auxiliary variables (for analysis)
      /* **
  output Q.Potential w(stateSelect=StateSelect.never) = reaction.chemical.mu 
    "Electrochemical overpotential";
  output Q.Current zI(stateSelect=StateSelect.never) = -reaction.chemical.Ndot 
    "Electrical current due to reaction";
  output Q.Number zJ_Apercm2(stateSelect=StateSelect.never) = zI*U.cm^2/(
    reaction.A*U.A) "Electrical current density, in A/cm2";
  output Q.Number Qdot[3](each stateSelect=StateSelect.never) = {potential[1].sT_actual
    *potential[1].chemical.Ndot,potential[2].sT_actual*potential[2].chemical.Ndot,
    current.sT_actual*current.chemical.Ndot} "Rates of thermal advection";
  */

      inner Conditions.Environment environment(analysis=true, T=360*U.K)
        annotation (Placement(transformation(extent={{50,70},{70,90}})));

      FCSys.Subregions.DepletionLayer depletion(
        axis=Axis.x,
        minoritySide=Side.n,
        transSubstrate=false,
        redeclare FCSys.Connectors.Reaction electrical) "Depletion layer"
        annotation (Placement(transformation(extent={{-10,10},{10,30}})));

      FCSys.Conditions.ByConnector.Inert.InertEfforts substrate(
        final inclTransX=inclTransX,
        final inclTransY=inclTransY,
        final inclTransZ=inclTransZ,
        thermalSource(y=environment.T))
        annotation (Placement(transformation(extent={{-10,6},{10,-14}})));

      Conditions.ByConnector.Face.Single.FaceEfforts majorityBC(materialSource(
            y=3.5*U.C/U.cc)) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={-44,20})));
      Conditions.ByConnector.Face.Single.FaceEfforts minorityBC(redeclare
          function normalSpec =
            FCSys.Conditions.ByConnector.Face.Single.TranslationalNormal.force)
        annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=90,
            origin={44,20})));
      replaceable Conditions.ByConnector.Chemical.Potential species1(
        final inclTransX=inclTransX,
        final inclTransY=inclTransY,
        final inclTransZ=inclTransZ,
        sT=U.V)
        annotation (Placement(transformation(extent={{-80,-20},{-60,0}})));
      replaceable Conditions.ByConnector.Chemical.Current species2(
        final inclTransX=inclTransX,
        final inclTransY=inclTransY,
        final inclTransZ=inclTransZ,
        sT=3000*U.K,
        redeclare Modelica.Blocks.Sources.Ramp source(duration=100, height=100*
              U.A)) constrainedby Conditions.ByConnector.Chemical.Current(sT=
            3000*U.K,redeclare Modelica.Blocks.Sources.Ramp source(duration=100,
            height=-100*U.A))
        annotation (Placement(transformation(extent={{-80,-50},{-60,-30}})));
      FCSys.Conditions.Adapters.ChemicalReaction exchange1(
        n=-2,
        m=U.g/U.mol,
        redeclare FCSys.Connectors.Reaction electrochem)
        annotation (Placement(transformation(extent={{-50,-30},{-30,-10}})));
      FCSys.Conditions.Adapters.ChemicalReaction exchange2(
        n=1,
        m=U.g/U.mol,
        redeclare FCSys.Connectors.Reaction electrochem)
        annotation (Placement(transformation(extent={{-50,-60},{-30,-40}})));

    protected
      final inner parameter Boolean inclTrans[Axis]={inclTransX,inclTransY,
          inclTransZ}
        "true, if each component of translational momentum is included";
      final inner parameter Integer n_trans=countTrue(inclTrans)
        "Number of components of translational momentum";
      final inner parameter Integer cartTrans[:]=index(inclTrans)
        "Cartesian-axis indices of the components of translational momentum";

    equation
      connect(depletion.inert, substrate.inert) annotation (Line(
          points={{6.10623e-16,16},{-3.36456e-22,8},{6.10623e-16,8},{
              6.10623e-16,6}},
          color={47,107,251},
          smooth=Smooth.None));

      connect(species2.chemical, exchange2.species) annotation (Line(
          points={{-70,-44},{-70,-50},{-44,-50}},
          color={255,195,38},
          smooth=Smooth.None));
      connect(exchange2.electrochem, depletion.electrical) annotation (Line(
          points={{-36,-50},{-20,-50},{-20,12},{-8,12},{0,20},{6.10623e-16,20}},

          color={255,195,38},
          smooth=Smooth.None));

      connect(species1.chemical, exchange1.species) annotation (Line(
          points={{-70,-14},{-70,-20},{-44,-20}},
          color={255,195,38},
          smooth=Smooth.None));
      connect(exchange2.electrochem, exchange1.electrochem) annotation (Line(
          points={{-36,-50},{-20,-50},{-20,-20},{-36,-20}},
          color={255,195,38},
          smooth=Smooth.None));
      connect(majorityBC.face, depletion.majority) annotation (Line(
          points={{-40,20},{-10,20}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(depletion.minority, minorityBC.face) annotation (Line(
          points={{10,20},{40,20}},
          color={127,127,127},
          smooth=Smooth.None));
      annotation (
        experiment(StopTime=120),
        Commands(file=
              "Resources/Scripts/Dymola/Subregions.Examples.Reaction.mos"
            "Subregions.Examples.Reaction.mos"),
        Diagram(graphics));
    end Reaction;

    model SaturationPressure
      "Evaluate the saturation pressure curve by varying temperature"
      import saturationPressureSI =
        Modelica.Media.Air.MoistAir.saturationPressureLiquid;

      output Q.Pressure p_sat=saturationPressureSI(subregion.gas.H2O.T/U.K)*U.Pa
        "Saturation pressure via Modelica.Media";
      output Q.Number T_degC=U.to_degC(subregion.gas.H2O.T)
        "Temperature in degree Celcius";

      extends Examples.Subregion(
        inclH2O=true,
        inclH2=false,
        subregion(liquid(inclH2O=inclH2O, H2O(consTransX=Conservation.IC)), gas(
              H2O(
              p_IC=saturationPressureSI(environment.T/U.K)*U.Pa,
              redeclare package Data = FCSys.Characteristics.H2O.Gas,
              consTransX=Conservation.IC,
              tauprime=1e-4*U.s))),
        environment(T=274.15*U.K));

      FCSys.Conditions.ByConnector.FaceBus.Single.FaceBusFlows BC1(liquid(
            inclH2O=true, H2O(
            redeclare function normalSpec =
                FCSys.Conditions.ByConnector.Face.Single.TranslationalNormal.velocity,

            redeclare function thermalSpec =
                FCSys.Conditions.ByConnector.Face.Single.Thermal.temperature,
            redeclare Modelica.Blocks.Sources.Ramp thermalSource(
              height=99*U.K,
              duration=3600,
              offset=environment.T,
              y))), gas(inclH2O=true, H2O(
            redeclare function normalSpec =
                FCSys.Conditions.ByConnector.Face.Single.TranslationalNormal.velocity,

            redeclare function thermalSpec =
                FCSys.Conditions.ByConnector.Face.Single.Thermal.temperature,
            redeclare Modelica.Blocks.Sources.Ramp thermalSource(
              height=99*U.K,
              duration=1000,
              offset=environment.T,
              y)))) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={-24,0})));

      FCSys.Conditions.ByConnector.FaceBus.Single.FaceBusFlows BC2(gas(inclH2O=
              true, H2O(
            redeclare function normalSpec =
                FCSys.Conditions.ByConnector.Face.Single.TranslationalNormal.velocity,

            redeclare function thermalSpec =
                FCSys.Conditions.ByConnector.Face.Single.Thermal.temperature,
            redeclare Modelica.Blocks.Sources.Ramp thermalSource(
              height=99*U.K,
              duration=1000,
              offset=environment.T,
              y))), liquid(inclH2O=true, H2O(
            redeclare function normalSpec =
                FCSys.Conditions.ByConnector.Face.Single.TranslationalNormal.velocity,

            redeclare function thermalSpec =
                FCSys.Conditions.ByConnector.Face.Single.Thermal.temperature,
            redeclare Modelica.Blocks.Sources.Ramp thermalSource(
              height=99*U.K,
              duration=1000,
              offset=environment.T,
              y)))) annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=90,
            origin={24,0})));

    equation
      connect(BC1.face, subregion.xNegative) annotation (Line(
          points={{-20,3.65701e-16},{-16,-3.36456e-22},{-16,6.10623e-16},{-10,
              6.10623e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      connect(subregion.xPositive, BC2.face) annotation (Line(
          points={{10,6.10623e-16},{16,6.10623e-16},{16,-2.54679e-16},{20,-2.54679e-16}},

          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      annotation (
        experiment(StopTime=3600, Tolerance=1e-08),
        Commands(file(ensureTranslated=true) =
            "Resources/Scripts/Dymola/Subregions.Examples.SaturationPressure.mos"
            "Subregions.Examples.SaturationPressure.mos"),
        experimentSetupOutput);
    end SaturationPressure;

    model Subregion
      "<html>Evaluate a single subregion, with H<sub>2</sub> by default</html>"

      extends Modelica.Icons.Example;

      parameter Boolean 'inclC+'=false
        "<html>Carbon plus (C<sup>+</sup>)</html>" annotation (choices(
            __Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean 'inclC19HF37O5S-'=false
        "<html>Nafion sulfonate (C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S<sup>-</sup>)</html>"
        annotation (choices(__Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean 'incle-'=false "<html>Electrons (e<sup>-</sup>)</html>"
        annotation (choices(__Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean 'inclH+'=false "<html>Protons (H<sup>+</sup>)</html>"
        annotation (choices(__Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean inclH2=true "<html>Hydrogen (H<sub>2</sub>)</html>"
        annotation (choices(__Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean inclH2O=false
        "<html>Water vapor (H<sub>2</sub>O)</html>" annotation (choices(
            __Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean inclN2=false "<html>Nitrogen (N<sub>2</sub>)</html>"
        annotation (choices(__Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean inclO2=false "<html>Oxygen (O<sub>2</sub>)</html>"
        annotation (choices(__Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));

      inner Conditions.Environment environment(analysis=true)
        annotation (Placement(transformation(extent={{30,32},{50,52}})));
      FCSys.Subregions.Subregion subregion(
        L={1,1,1}*U.cm,
        inclTransY=false,
        inclTransZ=false,
        inclFacesY=false,
        inclFacesZ=false,
        inclFacesX=true,
        graphite(
          final 'inclC+'='inclC+',
          final 'incle-'='incle-',
          'C+'(V_IC=subregion.V/4)),
        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2),
        liquid(H2O(V_IC=subregion.V/4)),
        ionomer(
          final 'inclC19HF37O5S-'='inclC19HF37O5S-',
          final 'inclH+'='inclH+',
          'C19HF37O5S-'(V_IC=subregion.V/4)))
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
      annotation (
        Documentation(info="<html><p>This model is boring.  It just sets up a
  single subregion with H<sub>2</sub> by default.  There are no boundary conditions
  other than those implied by the open connectors (no diffusion current, no forces, 
  no thermal conduction).  Other examples in this package are extended from this one.</p>
  </html>"),
        Placement(transformation(extent={{70,70},{90,90}})),
        experiment(StopTime=10),
        Commands(file(ensureTranslated=true) =
            "Resources/Scripts/Dymola/Subregions.Examples.Subregion.mos"
            "Subregions.Examples.Subregion.mos"));

    end Subregion;

    model Subregions
      "<html>Test a one-dimensional array of subregions with an initial pressure difference (C<sup>+</sup> and H<sub>2</sub> included by default)</html>"
      extends Modelica.Icons.Example;

      parameter Integer n_x=0
        "<html>Number of discrete subregions along the x axis, besides the 2 side subregions (<i>n<i><sub>x</sub>)</html>";
      parameter Q.Pressure Deltap_IC=100*U.Pa
        "<html>Initial pressure difference (&Delta;<i>p</i><sub>IC</sub>)</html>";
      parameter Boolean 'inclC+'=true
        "<html>Carbon plus (C<sup>+</sup>)</html>" annotation (choices(
            __Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean 'inclC19HF37O5S-'=false
        "<html>Nafion sulfonate (C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S)</html>"
        annotation (choices(__Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean 'incle-'=false "<html>Electrons (e<sup>-</sup>)</html>"
        annotation (choices(__Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean 'inclH+'=false "<html>Protons (H<sup>+</sup>)</html>"
        annotation (choices(__Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean inclH2=true "<html>Hydrogen (H<sub>2</sub>)</html>"
        annotation (choices(__Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean inclH2O=false
        "<html>Water vapor (H<sub>2</sub>O)</html>" annotation (choices(
            __Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean inclN2=false "<html>Nitrogen (N<sub>2</sub>)</html>"
        annotation (choices(__Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean inclO2=false "<html>Oxygen (O<sub>2</sub>)</html>"
        annotation (choices(__Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));

      inner Conditions.Environment environment(analysis=true)
        annotation (Placement(transformation(extent={{30,32},{50,52}})));
      FCSys.Subregions.Subregion subregion1(
        L={1,1,1}*U.cm,
        inclFacesY=false,
        inclFacesZ=false,
        inclTransY=false,
        inclTransZ=false,
        gas(
          reduceVel=true,
          reduceTemp=true,
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2,
          H2O(p_IC=environment.p - Deltap_IC/2),
          N2(p_IC=environment.p - Deltap_IC/2),
          O2(p_IC=environment.p - Deltap_IC/2),
          H2(p_IC=environment.p - Deltap_IC/2)),
        graphite(
          reduceTemp=true,
          final 'inclC+'='inclC+',
          final 'incle-'='incle-',
          'C+'(V_IC=subregion1.V/1000)),
        ionomer(
          reduceTemp=true,
          final 'inclC19HF37O5S-'='inclC19HF37O5S-',
          final 'inclH+'='inclH+',
          'C19HF37O5S-'(V_IC=subregion1.V/1000)),
        liquid(H2O(V_IC=subregion1.V/4)))
        annotation (Placement(transformation(extent={{-40,-10},{-20,10}})));

      FCSys.Subregions.Subregion subregions[n_x](
        each L={1,1,1}*U.cm,
        each inclTransY=false,
        each inclTransZ=false,
        each inclFacesY=false,
        each inclFacesZ=false,
        graphite(
          each reduceTemp=true,
          each final 'inclC+'='inclC+',
          each final 'incle-'='incle-',
          'C+'(each V_IC=subregions[1].V/1000)),
        ionomer(
          each reduceTemp=true,
          each final 'inclC19HF37O5S-'='inclC19HF37O5S-',
          each final 'inclH+'='inclH+',
          'C19HF37O5S-'(each V_IC=subregions[1].V/1000)),
        gas(
          each reduceTemp=true,
          each reduceVel=true,
          each final inclH2=inclH2,
          each final inclH2O=inclH2O,
          each final inclN2=inclN2,
          each final inclO2=inclO2,
          H2(p_IC={environment.p - Deltap_IC/2 - i*Deltap_IC/(n_x + 1) for i
                 in 1:n_x}),
          H2O(p_IC={environment.p - Deltap_IC/2 - i*Deltap_IC/(n_x + 1) for i
                 in 1:n_x}),
          N2(p_IC={environment.p - Deltap_IC/2 - i*Deltap_IC/(n_x + 1) for i
                 in 1:n_x}),
          O2(p_IC={environment.p - Deltap_IC/2 - i*Deltap_IC/(n_x + 1) for i
                 in 1:n_x})),
        liquid(H2O(each V_IC=subregions[1].V/4))) if n_x > 0
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

      FCSys.Subregions.Subregion subregion2(
        L={1,1,1}*U.cm,
        inclFacesY=false,
        inclFacesZ=false,
        inclTransY=false,
        inclTransZ=false,
        gas(
          reduceVel=true,
          reduceTemp=true,
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2,
          H2O(p_IC=environment.p + Deltap_IC/2),
          N2(p_IC=environment.p + Deltap_IC/2),
          O2(p_IC=environment.p + Deltap_IC/2),
          H2(p_IC=environment.p + Deltap_IC/2)),
        graphite(
          reduceTemp=true,
          final 'inclC+'='inclC+',
          final 'incle-'='incle-',
          'C+'(V_IC=subregion2.V/1000)),
        ionomer(
          reduceTemp=true,
          final 'inclC19HF37O5S-'='inclC19HF37O5S-',
          final 'inclH+'='inclH+',
          'C19HF37O5S-'(V_IC=subregion2.V/1000)))
        annotation (Placement(transformation(extent={{20,-10},{40,10}})));

      FCSys.Conditions.ByConnector.FaceBus.Single.FaceBusFlows BC1(
        ionomer(
          'inclC19HF37O5S-'=false,
          final 'inclH+'='inclH+',
          'H+'(redeclare function normalSpec =
                FCSys.Conditions.ByConnector.Face.Single.TranslationalNormal.velocity)),

        graphite(
          'inclC+'=false,
          final 'incle-'='incle-',
          'e-'(redeclare function normalSpec =
                FCSys.Conditions.ByConnector.Face.Single.TranslationalNormal.velocity)),

        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2,
          H2(redeclare function normalSpec =
                FCSys.Conditions.ByConnector.Face.Single.TranslationalNormal.velocity),

          H2O(redeclare function normalSpec =
                FCSys.Conditions.ByConnector.Face.Single.TranslationalNormal.velocity),

          N2(redeclare function normalSpec =
                FCSys.Conditions.ByConnector.Face.Single.TranslationalNormal.velocity),

          O2(redeclare function normalSpec =
                FCSys.Conditions.ByConnector.Face.Single.TranslationalNormal.velocity)))
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={-56,0})));

      FCSys.Conditions.ByConnector.FaceBus.Single.FaceBusFlows BC2(
        ionomer(
          'inclC19HF37O5S-'=false,
          final 'inclH+'='inclH+',
          'H+'(redeclare function normalSpec =
                FCSys.Conditions.ByConnector.Face.Single.TranslationalNormal.velocity)),

        graphite(
          'inclC+'=false,
          final 'incle-'='incle-',
          'e-'(redeclare function normalSpec =
                FCSys.Conditions.ByConnector.Face.Single.TranslationalNormal.velocity)),

        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2,
          H2(redeclare function normalSpec =
                FCSys.Conditions.ByConnector.Face.Single.TranslationalNormal.velocity),

          H2O(redeclare function normalSpec =
                FCSys.Conditions.ByConnector.Face.Single.TranslationalNormal.velocity),

          N2(redeclare function normalSpec =
                FCSys.Conditions.ByConnector.Face.Single.TranslationalNormal.velocity),

          O2(redeclare function normalSpec =
                FCSys.Conditions.ByConnector.Face.Single.TranslationalNormal.velocity)))
        annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=90,
            origin={56,0})));

    equation
      connect(subregion1.xPositive, subregions[1].xNegative) annotation (Line(
          points={{-20,6.10623e-16},{-16,-3.36456e-22},{-16,6.10623e-16},{-10,
              6.10623e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      for i in 1:n_x - 1 loop
        connect(subregions[i].xPositive, subregions[i + 1].xNegative)
          "Not shown in the diagram";
      end for;
      if n_x == 0 then
        connect(subregion1.xPositive, subregion2.xNegative)
          "Not shown in the diagram";
      end if;
      connect(subregions[n_x].xPositive, subregion2.xNegative) annotation (Line(
          points={{10,6.10623e-16},{20,-3.36456e-22},{20,6.10623e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      connect(BC1.face, subregion1.xNegative) annotation (Line(
          points={{-52,3.65701e-16},{-46,3.65701e-16},{-46,6.10623e-16},{-40,
              6.10623e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      connect(subregion2.xPositive, BC2.face) annotation (Line(
          points={{40,6.10623e-16},{46,6.10623e-16},{46,-2.54679e-16},{52,-2.54679e-16}},

          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      annotation (
        Placement(transformation(extent={{70,70},{90,90}})),
        experiment(
          StopTime=2,
          Tolerance=1e-06,
          Algorithm="Dassl"),
        Commands(file(ensureTranslated=true) =
            "Resources/Scripts/Dymola/Subregions.Examples.Subregions.mos"
            "Subregions.Examples.Subregions.mos"));
    end Subregions;

    model ThermalConduction "Test thermal conduction (through solid)"
      extends Examples.Subregions(
        n_x=8,
        'inclC+'=true,
        inclH2=false,
        subregion1(graphite('C+'(T_IC=environment.T + 30*U.K, initMaterial=
                  InitScalar.Pressure))),
        subregions(graphite('C+'(each initMaterial=InitScalar.Pressure))),
        subregion2(graphite('C+'(initMaterial=InitScalar.Pressure))));

      annotation (
        Commands(file=
              "Resources/Scripts/Dymola/Subregions.Examples.ThermalConduction.mos"
            "Subregions.Examples.ThermalConduction.mos"),
        experiment(StopTime=500, Algorithm="Dassl"),
        experimentSetupOutput);

    end ThermalConduction;

    model ThermalConductionConvection
      "Test combined thermal conduction and convection"

      extends Examples.Subregions(
        n_x=8,
        'inclC+'=true,
        inclH2=false,
        inclN2=true,
        subregion1(gas(N2(
              T_IC=environment.T + 30*U.K,
              p_IC=environment.p,
              T(displayUnit="degC"),
              phi(displayUnit="mm/s"))), graphite('C+'(T_IC=environment.T + 30*
                  U.K, V_IC=0.5*subregion1.V))),
        subregions(gas(N2(
              each p_IC=environment.p,
              each T(displayUnit="degC"),
              phi(each displayUnit="mm/s"))), graphite('C+'(each V_IC=0.5*
                  subregions[1].V))),
        subregion2(gas(N2(p_IC=environment.p, phi(displayUnit="mm/s"))),
            graphite('C+'(V_IC=0.5*subregion2.V, T(displayUnit="degC")))));

      annotation (
        Commands(file=
              "Resources/Scripts/Dymola/Subregions.Examples.ThermalConductionConvection.mos"
            "Subregions.Examples.ThermalConductionConvection.mos"),
        experiment(
          StopTime=140,
          NumberOfIntervals=5000,
          Tolerance=1e-06,
          Algorithm="Dassl"),
        experimentSetupOutput);

    end ThermalConductionConvection;

    model ChargeLayer
      import FCSys.BaseClasses.Utilities.cartWrap;
      import FCSys.BaseClasses.Utilities.countTrue;
      import FCSys.BaseClasses.Utilities.enumerate;
      import FCSys.BaseClasses.Utilities.index;

      // Geometry
      inner parameter Q.Length L[Axis](each min=Modelica.Constants.small) = {10
        *U.cm,10*U.cm,0.1*U.mm} "Lengths";
      final inner parameter Q.Area A[Axis]={L[cartWrap(axis + 1)]*L[cartWrap(
          axis + 2)] for axis in Axis} "Cross-sectional areas";

      // Assumptions
      // -----------
      // Included components of translational momentum
      parameter Boolean inclTransX=true "X" annotation (choices(
            __Dymola_checkBox=true), Dialog(
          tab="Assumptions",
          group="Axes with translational momentum included",
          compact=true));
      parameter Boolean inclTransY=false "Y" annotation (choices(
            __Dymola_checkBox=true), Dialog(
          tab="Assumptions",
          group="Axes with translational momentum included",
          compact=true));
      parameter Boolean inclTransZ=false "Z" annotation (choices(
            __Dymola_checkBox=true), Dialog(
          tab="Assumptions",
          group="Axes with translational momentum included",
          compact=true));

      FCSys.Conditions.Adapters.ChemicalFace chargeLayer(redeclare
          FCSys.Connectors.Reaction electrical)
        annotation (Placement(transformation(extent={{-6,10},{14,30}})));
      FCSys.Conditions.ByConnector.Reaction.ReactionEfforts electrochem(
        redeclare FCSys.Connectors.Reaction electrochem,
        final inclTransX=inclTransX,
        final inclTransY=inclTransY,
        final inclTransZ=inclTransZ)
        annotation (Placement(transformation(extent={{22,10},{42,30}})));
      Conditions.ByConnector.Face.Single.FaceEfforts face(redeclare function
          normalSpec =
            FCSys.Conditions.ByConnector.Face.Single.TranslationalNormal.currentDensity)
        annotation (Placement(transformation(extent={{-36,14},{-16,34}})));
    protected
      final inner parameter Boolean inclTrans[Axis]={inclTransX,inclTransY,
          inclTransZ}
        "true, if each component of translational momentum is included";
      final inner parameter Integer n_trans=countTrue(inclTrans)
        "Number of components of translational momentum";
    equation
      connect(chargeLayer.electrical, reaction.electrochem) annotation (Line(
          points={{8,20},{18,20},{18,10},{32,10}},
          color={255,195,38},
          smooth=Smooth.None));
      connect(face.face, chargeLayer.face) annotation (Line(
          points={{-26,20},{-12.2,20},{-12.2,20},{6.66134e-16,20}},
          color={127,127,127},
          smooth=Smooth.None));
      annotation (Diagram(graphics));
    end ChargeLayer;
  end Examples;

  model Subregion "Subregion with all phases"
    extends FCSys.Subregions.BaseClasses.EmptySubregion(final hasSpecies=gas.n_spec
           + graphite.n_spec + ionomer.n_spec + liquid.n_spec > 0);

    FCSys.Phases.Gas gas(
      inclH2O=true,
      final n_faces=n_faces,
      final inclHOR=inclHOR,
      final inclORR=inclORR) "Gas" annotation (Dialog(group=
            "Phases (click to edit)"), Placement(transformation(extent={{-10,-10},
              {10,10}})));

    FCSys.Phases.Graphite graphite(
      final n_faces=n_faces,
      final inclHOR=inclHOR,
      final inclORR=inclORR) "Graphite" annotation (Dialog(group=
            "Phases (click to edit)"), Placement(transformation(extent={{-10,-10},
              {10,10}})));

    FCSys.Phases.Ionomer ionomer(
      final n_faces=n_faces,
      final inclHOR=inclHOR,
      final inclORR=inclORR) "Ionomer" annotation (Dialog(group=
            "Phases (click to edit)"), Placement(transformation(extent={{-10,-10},
              {10,10}})));

    FCSys.Phases.Liquid liquid(final n_faces=n_faces) "Liquid" annotation (
        Dialog(group="Phases (click to edit)"), Placement(transformation(extent
            ={{-10,-10},{10,10}})));

  protected
    final parameter Boolean inclHOR=graphite.'incle-' and ionomer.'inclH+' and
        gas.inclH2 "Include the hydrogen oxidation reaction";
    final parameter Boolean inclORR=graphite.'incle-' and ionomer.'inclH+' and
        gas.inclO2 and gas.inclH2O "Include the oxygen reduction reaction";

    FCSys.Conditions.Adapters.AmagatDalton gasDA(final n_trans=n_trans) if gas.n_spec
       > 0 and not void
      "Inerface between Dalton and Amagat representations of gas"
      annotation (Placement(transformation(extent={{22,-10},{2,-30}})));
    FCSys.Conditions.Adapters.AmagatDalton graphiteDA(final n_trans=n_trans)
      if graphite.n_spec > 0 and not void
      "Inerface between Dalton and Amagat representations of graphite"
      annotation (Placement(transformation(extent={{22,-10},{2,-30}})));
    FCSys.Conditions.Adapters.AmagatDalton ionomerDA(final n_trans=n_trans) if
      ionomer.n_spec > 0 and not void
      "Inerface between Dalton and Amagat representations of ionomer"
      annotation (Placement(transformation(extent={{22,-10},{2,-30}})));
    FCSys.Conditions.Adapters.AmagatDalton liquidDA(final n_trans=n_trans) if
      liquid.n_spec > 0 and not void
      "Inerface between Dalton and Amagat representations of liquid"
      annotation (Placement(transformation(extent={{22,-10},{2,-30}})));

    Connectors.PhysicalBusInternal physical
      "Internal connector for phase change" annotation (Placement(
          transformation(extent={{-24,4},{-4,24}}), iconTransformation(extent={
              {-80,28},{-60,48}})));

  equation
    // Phase change
    connect(gas.physical, physical) annotation (Line(
        points={{-5.1,4.9},{-14,14}},
        color={38,196,52},
        smooth=Smooth.None,
        thickness=0.5));
    connect(ionomer.physical, physical) annotation (Line(
        points={{-5.1,4.9},{-14,14}},
        color={38,196,52},
        smooth=Smooth.None,
        thickness=0.5));
    connect(liquid.physical, physical) annotation (Line(
        points={{-5.1,4.9},{-14,14}},
        color={38,196,52},
        smooth=Smooth.None,
        thickness=0.5));

    // Reactions
    connect(gas.chemical, graphite.chemical) annotation (Line(
        points={{6.10623e-16,-4.4},{-3.36456e-22,-4},{6.10623e-16,-4},{
            6.10623e-16,-4.4}},
        color={239,142,1},
        smooth=Smooth.None));

    connect(graphite.chemical, ionomer.chemical) annotation (Line(
        points={{6.10623e-16,-4.4},{6.10623e-16,-4.4},{6.10623e-16,-4.4}},
        color={239,142,1},
        smooth=Smooth.None));

    // Phases
    // ------
    // Gas
    connect(gas.dalton, gasDA.dalton) annotation (Line(
        points={{8,-8},{8,-20},{8,-20},{8,-20}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(gasDA.amagat, volume.amagat) annotation (Line(
        points={{16,-20},{16,-13},{16,-6},{16,-6}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(gas.xNegative, xNegative.gas) annotation (Line(
        points={{-8,6.10623e-16},{-40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.xPositive, xPositive.gas) annotation (Line(
        points={{8,6.10623e-16},{40,6.10623e-16},{40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    connect(gas.yNegative, yNegative.gas) annotation (Line(
        points={{6.10623e-16,-8.4},{-4.87687e-22,-8.4},{-4.87687e-22,-40},{
            5.55112e-16,-40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.yPositive, yPositive.gas) annotation (Line(
        points={{6.10623e-16,10},{6.10623e-16,40},{5.55112e-16,40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.zNegative, zNegative.gas) annotation (Line(
        points={{5,5},{20,20}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.zPositive, zPositive.gas) annotation (Line(
        points={{-8,-8},{-20,-20}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    // Graphite
    connect(graphite.dalton, graphiteDA.dalton) annotation (Line(
        points={{8,-8},{8,-20},{8,-20},{8,-20}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(graphiteDA.amagat, volume.amagat) annotation (Line(
        points={{16,-20},{16,-13},{16,-6},{16,-6}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(graphite.xNegative, xNegative.graphite) annotation (Line(
        points={{-8,6.10623e-16},{-40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.xPositive, xPositive.graphite) annotation (Line(
        points={{8,6.10623e-16},{40,6.10623e-16},{40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    connect(graphite.yNegative, yNegative.graphite) annotation (Line(
        points={{6.10623e-16,-8.4},{-4.87687e-22,-8.4},{-4.87687e-22,-40},{
            5.55112e-16,-40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.yPositive, yPositive.graphite) annotation (Line(
        points={{6.10623e-16,10},{6.10623e-16,40},{5.55112e-16,40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.zNegative, zNegative.graphite) annotation (Line(
        points={{5,5},{20,20}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.zPositive, zPositive.graphite) annotation (Line(
        points={{-8,-8},{-20,-20}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    // Ionomer
    connect(ionomer.dalton, ionomerDA.dalton) annotation (Line(
        points={{8,-8},{8,-20},{8,-20},{8,-20}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(ionomerDA.amagat, volume.amagat) annotation (Line(
        points={{16,-20},{16,-13},{16,-6},{16,-6}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(ionomer.xNegative, xNegative.ionomer) annotation (Line(
        points={{-8,6.10623e-16},{-40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(ionomer.xPositive, xPositive.ionomer) annotation (Line(
        points={{8,6.10623e-16},{40,6.10623e-16},{40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    connect(ionomer.yNegative, yNegative.ionomer) annotation (Line(
        points={{6.10623e-16,-8.4},{-4.87687e-22,-8.4},{-4.87687e-22,-40},{
            5.55112e-16,-40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(ionomer.yPositive, yPositive.ionomer) annotation (Line(
        points={{6.10623e-16,10},{5.55112e-16,40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(ionomer.zNegative, zNegative.ionomer) annotation (Line(
        points={{5,5},{20,20}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(ionomer.zPositive, zPositive.ionomer) annotation (Line(
        points={{-8,-8},{-20,-20}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    // Liquid
    connect(liquid.dalton, liquidDA.dalton) annotation (Line(
        points={{8,-8},{8,-14},{8,-14},{8,-20}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(liquidDA.amagat, volume.amagat) annotation (Line(
        points={{16,-20},{16,-6}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(liquid.xNegative, xNegative.liquid) annotation (Line(
        points={{-8,6.10623e-16},{-40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    connect(liquid.xPositive, xPositive.liquid) annotation (Line(
        points={{8,6.10623e-16},{40,6.10623e-16},{40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    connect(liquid.yNegative, yNegative.liquid) annotation (Line(
        points={{6.10623e-16,-8.4},{6.10623e-16,-40},{5.55112e-16,-40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(liquid.yPositive, yPositive.liquid) annotation (Line(
        points={{6.10623e-16,10},{-4.87687e-22,40},{5.55112e-16,40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    connect(liquid.zNegative, zNegative.liquid) annotation (Line(
        points={{5,5},{20,20}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(liquid.zPositive, zPositive.liquid) annotation (Line(
        points={{-8,-8},{-20,-20}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    annotation (Documentation(info="<html>
  <p>At least one configuration must be included.  H<sub>2</sub>O gas is included by default.</p>
     
   <p>Please see the documentation of the
   <a href=\"modelica://FCSys.Subregions.BaseClasses.EmptySubregion\">EmptySubregion</a> model.</p></html>"),
        Diagram(graphics));
  end Subregion;

  model SubregionIonomerOnly "Subregion with only the ionomer phase"
    extends FCSys.Subregions.BaseClasses.EmptySubregion(final hasSpecies=
          ionomer.n_spec > 0);

    FCSys.Phases.Ionomer ionomer(
      final n_faces=n_faces,
      inclH2O=true,
      final inclHOR=false,
      final inclORR=false) "Ionomer" annotation (Dialog(group=
            "Phases (click to edit)"), Placement(transformation(extent={{-10,-10},
              {10,10}})));

  protected
    FCSys.Conditions.Adapters.AmagatDalton ionomerDA(final n_trans=n_trans) if
      ionomer.n_spec > 0 and not void
      "Inerface between Dalton and Amagat representations of ionomer"
      annotation (Placement(transformation(extent={{22,-10},{2,-30}})));
  equation
    // Phases
    // ------
    // Ionomer
    connect(ionomer.dalton, ionomerDA.dalton) annotation (Line(
        points={{8,-8},{8,-20},{8,-20}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(ionomerDA.amagat, volume.amagat) annotation (Line(
        points={{16,-20},{16,-6},{16,-6}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(ionomer.xNegative, xNegative.ionomer) annotation (Line(
        points={{-8,6.10623e-16},{-40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(ionomer.xPositive, xPositive.ionomer) annotation (Line(
        points={{8,6.10623e-16},{8,-4.87687e-22},{40,-4.87687e-22},{40,
            5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    connect(ionomer.yNegative, yNegative.ionomer) annotation (Line(
        points={{6.10623e-16,-8.4},{6.10623e-16,-40},{5.55112e-16,-40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(ionomer.yPositive, yPositive.ionomer) annotation (Line(
        points={{6.10623e-16,10},{5.55112e-16,40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(ionomer.zNegative, zNegative.ionomer) annotation (Line(
        points={{5,5},{20,20}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(ionomer.zPositive, zPositive.ionomer) annotation (Line(
        points={{-8,-8},{-20,-20}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(ionomer.dalton, ionomerDA.dalton) annotation (Line(
        points={{8,-8},{8,-14},{8,-20},{8,-20}},
        color={47,107,251},
        smooth=Smooth.None));

    connect(ionomerDA.amagat, volume.amagat) annotation (Line(
        points={{16,-20},{16,-13},{16,-6},{16,-6}},
        color={47,107,251},
        smooth=Smooth.None));
    annotation (
      defaultComponentName="subregion",
      Documentation(info="<html>
  <p>At least one configuration must be included.  H<sub>2</sub>O is included by default.</p>

<p>Please see the documentation of the
   <a href=\"modelica://FCSys.Subregions.BaseClasses.EmptySubregion\">EmptySubregion</a> model.</p></html>"),

      Diagram(graphics));
  end SubregionIonomerOnly;

  model SubregionNoIonomer "Subregion with all phases except ionomer"
    extends FCSys.Subregions.BaseClasses.EmptySubregion(final hasSpecies=gas.n_spec
           + graphite.n_spec + liquid.n_spec > 0);

    FCSys.Phases.Gas gas(
      inclH2O=true,
      final n_faces=n_faces,
      final inclHOR=false,
      final inclORR=false) "Gas" annotation (Dialog(group=
            "Phases (click to edit)"), Placement(transformation(extent={{-10,-10},
              {10,10}})));
    FCSys.Phases.Graphite graphite(
      final n_faces=n_faces,
      final inclHOR=false,
      final inclORR=false) "Graphite" annotation (Dialog(group=
            "Phases (click to edit)"), Placement(transformation(extent={{-10,-10},
              {10,10}})));
    FCSys.Phases.Liquid liquid(final n_faces=n_faces) "Liquid" annotation (
        Dialog(group="Phases (click to edit)"), Placement(transformation(extent
            ={{-10,-10},{10,10}})));

  protected
    Connectors.PhysicalBusInternal physical "Connector for phase change"
      annotation (Placement(transformation(extent={{-30,10},{-10,30}}),
          iconTransformation(extent={{-76,16},{-56,36}})));

    FCSys.Conditions.Adapters.AmagatDalton gasDA(final n_trans=n_trans) if gas.n_spec
       > 0 and not void
      "Inerface between Dalton and Amagat representations of gas"
      annotation (Placement(transformation(extent={{22,-10},{2,-30}})));
    FCSys.Conditions.Adapters.AmagatDalton graphiteDA(final n_trans=n_trans)
      if graphite.n_spec > 0 and not void
      "Inerface between Dalton and Amagat representations of graphite"
      annotation (Placement(transformation(extent={{22,-10},{2,-30}})));
    FCSys.Conditions.Adapters.AmagatDalton liquidDA(final n_trans=n_trans) if
      liquid.n_spec > 0 and not void
      "Inerface between Dalton and Amagat representations of liquid"
      annotation (Placement(transformation(extent={{22,-10},{2,-30}})));
  equation
    // Phase change
    connect(physical, gas.physical);
    connect(physical, liquid.physical) annotation (Line(
        points={{-20,20},{-5.1,4.9}},
        color={38,196,52},
        smooth=Smooth.None,
        thickness=0.5));

    // Phases
    // ------
    // Gas
    connect(gas.dalton, gasDA.dalton) annotation (Line(
        points={{8,-8},{8,-20},{8,-20}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(gasDA.amagat, volume.amagat) annotation (Line(
        points={{16,-20},{16,-6},{16,-6}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(gas.xNegative, xNegative.gas) annotation (Line(
        points={{-8,6.10623e-16},{-40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.xPositive, xPositive.gas) annotation (Line(
        points={{8,6.10623e-16},{8,-4.87687e-22},{40,-4.87687e-22},{40,
            5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    connect(gas.yNegative, yNegative.gas) annotation (Line(
        points={{6.10623e-16,-8.4},{6.10623e-16,-40},{5.55112e-16,-40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.yPositive, yPositive.gas) annotation (Line(
        points={{6.10623e-16,10},{0,10},{0,10},{-4.87687e-22,40},{-4.87687e-22,
            40},{5.55112e-16,40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.zNegative, zNegative.gas) annotation (Line(
        points={{5,5},{20,20}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.zPositive, zPositive.gas) annotation (Line(
        points={{-8,-8},{-20,-20}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    // Graphite
    connect(graphite.dalton, graphiteDA.dalton) annotation (Line(
        points={{8,-8},{8,-20},{8,-20}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(graphiteDA.amagat, volume.amagat) annotation (Line(
        points={{16,-20},{16,-6},{16,-6}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(graphite.xNegative, xNegative.graphite) annotation (Line(
        points={{-8,6.10623e-16},{-40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.xPositive, xPositive.graphite) annotation (Line(
        points={{8,6.10623e-16},{8,-4.87687e-22},{40,-4.87687e-22},{40,
            5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    connect(graphite.yNegative, yNegative.graphite) annotation (Line(
        points={{6.10623e-16,-8.4},{6.10623e-16,-40},{5.55112e-16,-40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.yPositive, yPositive.graphite) annotation (Line(
        points={{6.10623e-16,10},{0,10},{0,10},{-4.87687e-22,40},{-4.87687e-22,
            40},{5.55112e-16,40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.zNegative, zNegative.graphite) annotation (Line(
        points={{5,5},{20,20}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.zPositive, zPositive.graphite) annotation (Line(
        points={{-8,-8},{-20,-20}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    // Liquid
    connect(liquid.dalton, liquidDA.dalton) annotation (Line(
        points={{8,-8},{8,-20},{8,-20}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(liquidDA.amagat, volume.amagat) annotation (Line(
        points={{16,-20},{16,-6},{16,-6}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(liquid.xNegative, xNegative.liquid) annotation (Line(
        points={{-8,6.10623e-16},{-24,0},{-40,5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    connect(liquid.xPositive, xPositive.liquid) annotation (Line(
        points={{8,6.10623e-16},{8,-4.87687e-22},{40,-4.87687e-22},{40,
            5.55112e-16}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    connect(liquid.yNegative, yNegative.liquid) annotation (Line(
        points={{6.10623e-16,-8.4},{-4.87687e-22,-8.4},{-4.87687e-22,-40},{
            5.55112e-16,-40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(liquid.yPositive, yPositive.liquid) annotation (Line(
        points={{6.10623e-16,10},{0,10},{-4.87687e-22,40},{5.55112e-16,40}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    connect(liquid.zNegative, zNegative.liquid) annotation (Line(
        points={{5,5},{20,20}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(liquid.zPositive, zPositive.liquid) annotation (Line(
        points={{-8,-8},{-20,-20}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gasDA.amagat, volume.amagat) annotation (Line(
        points={{16,-20},{16,-13},{16,-6},{16,-6}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(graphite.dalton, graphiteDA.dalton) annotation (Line(
        points={{8,-8},{8,-14},{8,-20},{8,-20}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(graphiteDA.amagat, volume.amagat) annotation (Line(
        points={{16,-20},{16,-13},{16,-6},{16,-6}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(liquid.dalton, liquidDA.dalton) annotation (Line(
        points={{8,-8},{8,-20},{8,-20}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(liquidDA.amagat, volume.amagat) annotation (Line(
        points={{16,-20},{16,-13},{16,-6},{16,-6}},
        color={47,107,251},
        smooth=Smooth.None));
    annotation (defaultComponentName="subregion", Documentation(info="<html>
  <p>At least one configuration must be included.  H<sub>2</sub>O gas is included by default.</p>

<p>Please see the documentation of the
   <a href=\"modelica://FCSys.Subregions.BaseClasses.EmptySubregion\">EmptySubregion</a> model.</p></html>"));
  end SubregionNoIonomer;

  package BaseClasses "Base classes (generally not for direct use)"

    extends Modelica.Icons.BasesPackage;

    partial model EmptySubregion
      "Base model for multi-dimensional, multi-species storage, transport, and exchange"
      import FCSys.BaseClasses.Utilities.cartWrap;
      import FCSys.BaseClasses.Utilities.countTrue;
      import FCSys.BaseClasses.Utilities.enumerate;
      import FCSys.BaseClasses.Utilities.index;
      // extends FCSys.BaseClasses.Icons.Names.Top3;

      // Geometric parameters
      inner parameter Q.Length L[Axis](each min=Modelica.Constants.small) = {U.cm,
        U.cm,U.cm} "<html>Lengths (<b>L</b>)</html>"
        annotation (Dialog(group="Geometry"));
      final inner parameter Q.Area A[Axis]={L[cartWrap(axis + 1)]*L[cartWrap(
          axis + 2)] for axis in Axis} "Cross-sectional areas";
      final inner parameter Q.Volume V=product(L) "Volume";

      // Assumptions
      // -----------
      // Included components of translational momentum
      parameter Boolean inclTransX=true "X" annotation (choices(
            __Dymola_checkBox=true), Dialog(
          tab="Assumptions",
          group="Axes with translational momentum included",
          compact=true));
      parameter Boolean inclTransY=true "Y" annotation (choices(
            __Dymola_checkBox=true), Dialog(
          tab="Assumptions",
          group="Axes with translational momentum included",
          compact=true));
      parameter Boolean inclTransZ=true "Z" annotation (choices(
            __Dymola_checkBox=true), Dialog(
          tab="Assumptions",
          group="Axes with translational momentum included",
          compact=true));
      //
      // Included faces
      parameter Boolean inclFacesX=true "X" annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          tab="Assumptions",
          group="Axes with faces included",
          compact=true));
      parameter Boolean inclFacesY=true "Y" annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          tab="Assumptions",
          group="Axes with faces included",
          compact=true));
      parameter Boolean inclFacesZ=true "Z" annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          tab="Assumptions",
          group="Axes with faces included",
          compact=true));
      //
      parameter Boolean void=false "Volume not imposed (zero pressure instead)"
        annotation (Dialog(tab="Assumptions", compact=true), choices(
            __Dymola_checkBox=true));

      Connectors.FaceBus xNegative if inclFacesX
        "Negative face along the x axis" annotation (Placement(transformation(
              extent={{-50,-10},{-30,10}}), iconTransformation(extent={{-110,-10},
                {-90,10}})));
      Connectors.FaceBus xPositive if inclFacesX
        "Positive face along the x axis" annotation (Placement(transformation(
              extent={{30,-10},{50,10}}), iconTransformation(extent={{90,-10},{
                110,10}})));
      Connectors.FaceBus yNegative if inclFacesY
        "Negative face along the y axis" annotation (Placement(transformation(
              extent={{-10,-50},{10,-30}}), iconTransformation(extent={{-10,-110},
                {10,-90}})));
      Connectors.FaceBus yPositive if inclFacesY
        "Positive face along the y axis" annotation (Placement(transformation(
              extent={{-10,30},{10,50}}), iconTransformation(extent={{-10,90},{
                10,110}})));
      Connectors.FaceBus zNegative if inclFacesZ
        "Negative face along the z axis" annotation (Placement(transformation(
              extent={{10,10},{30,30}}), iconTransformation(extent={{40,40},{60,
                60}})));
      Connectors.FaceBus zPositive if inclFacesZ
        "Positive face along the z axis" annotation (Placement(transformation(
              extent={{-30,-30},{-10,-10}}), iconTransformation(extent={{-60,-60},
                {-40,-40}})));

    protected
      parameter Boolean hasSpecies "true, if any species are included";
      final inner parameter Boolean inclTrans[Axis]={inclTransX,inclTransY,
          inclTransZ}
        "true, if each component of translational momentum is included";
      final inner parameter Boolean inclFaces[Axis]={inclFacesX,inclFacesY,
          inclFacesZ} "true, if each pairs of faces is included";
      final inner parameter Boolean inclRot[Axis]={inclFacesY and inclFacesZ,
          inclFacesZ and inclFacesX,inclFacesX and inclFacesY}
        "true, if each axis of rotation has all its tangential faces included";
      final inner parameter Integer n_trans=countTrue(inclTrans)
        "Number of components of translational momentum";
      final inner parameter Integer n_faces=countTrue(inclFaces)
        "Number of pairs of faces";
      final inner parameter Integer cartTrans[n_trans]=index(inclTrans)
        "Cartesian-axis indices of the components of translational momentum";
      final inner parameter Integer cartFaces[n_faces]=index(inclFaces)
        "Cartesian-axis indices of the pairs of faces";
      final inner parameter Integer cartRot[:]=index(inclRot)
        "Cartesian-axis indices of the components of rotational momentum";
      final inner parameter Integer transCart[Axis]=enumerate(inclTrans)
        "Translational-momentum-component indices of the Cartesian axes";
      final inner parameter Integer facesCart[Axis]=enumerate(inclFaces)
        "Face-pair indices of the Cartesian axes";

      FCSys.Conditions.ByConnector.Amagat.Volume volume(V=V, final n_trans=
            n_trans) if hasSpecies and not void
        "Model to establish a fixed total volume"
        annotation (Placement(transformation(extent={{-16,-16},{16,16}})));
      annotation (
        defaultComponentPrefixes="replaceable",
        defaultComponentName="subregion",
        Documentation(info="<html>
  <p>At least one component of translational momentum must be included.
  All of the components are included by default.</p>

    <p>At least one pair of faces must be included.
  All of the faces are included by default.</p>

    <p>If <code>void</code> is <code>true</code>, then the total pressure of the species in each
    phase sums to zero.  The total volume is then typically determined by one of the species in each phase (e.g., the solid).
    may be at the species level.  This is useful to when all the species are incompressible.</p>

  <p>This model should be extended to include the appropriate phases and reactions.</p>
  </html>"),
        Icon(graphics={Line(
                  points={{-100,0},{-40,0}},
                  color={127,127,127},
                  thickness=0.5,
                  visible=inclFacesX,
                  smooth=Smooth.None),Line(
                  points={{0,-40},{0,-100}},
                  color={127,127,127},
                  thickness=0.5,
                  visible=inclFacesY,
                  smooth=Smooth.None),Line(
                  points={{40,40},{50,50}},
                  color={127,127,127},
                  thickness=0.5,
                  visible=inclFacesZ,
                  smooth=Smooth.None),Polygon(
                  points={{-40,16},{-16,40},{40,40},{40,-16},{16,-40},{-40,-40},
                {-40,16}},
                  lineColor={127,127,127},
                  smooth=Smooth.None,
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid),Line(
                  points={{-40,-40},{-16,-16}},
                  color={127,127,127},
                  smooth=Smooth.None,
                  pattern=LinePattern.Dash),Line(
                  points={{-16,40},{-16,-16},{40,-16}},
                  color={127,127,127},
                  smooth=Smooth.None,
                  pattern=LinePattern.Dash),Line(
                  points={{-40,0},{28,0}},
                  color={210,210,210},
                  visible=inclFacesX,
                  smooth=Smooth.None,
                  thickness=0.5),Line(
                  points={{0,28},{0,-40}},
                  color={210,210,210},
                  visible=inclFacesY,
                  smooth=Smooth.None,
                  thickness=0.5),Line(
                  points={{28,0},{100,0}},
                  color={127,127,127},
                  thickness=0.5,
                  visible=inclFacesX,
                  smooth=Smooth.None),Line(
                  points={{0,100},{0,28}},
                  color={127,127,127},
                  thickness=0.5,
                  visible=inclFacesY,
                  smooth=Smooth.None),Line(
                  points={{-12,-12},{40,40}},
                  color={210,210,210},
                  visible=inclFacesZ,
                  smooth=Smooth.None,
                  thickness=0.5),Line(
                  points={{-40,16},{16,16},{16,-40}},
                  color={127,127,127},
                  smooth=Smooth.None),Line(
                  points={{-50,-50},{-12,-12}},
                  color={127,127,127},
                  thickness=0.5,
                  visible=inclFacesZ,
                  smooth=Smooth.None),Polygon(
                  points={{-40,16},{-16,40},{40,40},{40,-16},{16,-40},{-40,-40},
                {-40,16}},
                  lineColor={127,127,127},
                  smooth=Smooth.None),Line(
                  points={{40,40},{16,16}},
                  color={127,127,127},
                  smooth=Smooth.None),Text(
                  extent={{-100,56},{100,96}},
                  textString="%name",
                  lineColor={0,0,0})}));

    end EmptySubregion;

  end BaseClasses;

  annotation (Documentation(info="<html>
<p><b>Licensed by the Georgia Tech Research Corporation under the Modelica License 2</b><br>
Copyright 2007&ndash;2012, Georgia Tech Research Corporation.</p>

<p><i>This Modelica package is <u>free</u> software and the use is completely at <u>your own risk</u>;
it can be redistributed and/or modified under the terms of the Modelica License 2. For license conditions (including the
disclaimer of warranty) see <a href=\"modelica://FCSys.UsersGuide.ModelicaLicense2\">
FCSys.UsersGuide.ModelicaLicense2</a> or visit <a href=\"http://www.modelica.org/licenses/ModelicaLicense2\">
http://www.modelica.org/licenses/ModelicaLicense2</a>.</i></p>
</html>"));
end Subregions;

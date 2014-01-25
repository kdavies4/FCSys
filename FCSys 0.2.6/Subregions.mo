within FCSys;
package Subregions "Control volumes with multi-species transfer and storage"
  package Examples "Examples"
    extends Modelica.Icons.ExamplesPackage;

    package PhaseChange "Examples of phase change"
      extends Modelica.Icons.ExamplesPackage;
      model Evaporation "<html>Evaporation of H<sub>2</sub>O</html>"

        output Q.Pressure p_sat=Characteristics.H2O.p_sat(subregion.gas.H2O.T)
          "Saturation pressure via Modelica.Media";

        extends Examples.Subregion(
          inclH2O=true,
          inclH2=false,
          'inclC+'=true,
          subregion(liquid(inclH2O=inclH2O, H2O(epsilon_IC=0.001)), gas(H2O(
                  p_IC=U.kPa))));

        annotation (
          Documentation(info="<html><p>Initially, the water vapor is below saturation and a small amount of liquid water is present (1/1000 of the total volume).
  Some of the liquid evaporates until saturation is reached. The boundaries are adiabatic; therefore, the temperature of the liquid and the gas

  decreases due to the enthalpy of formation.</p>

  <p>See also <a href=\"modelica://FCSys.Characteristics.Examples.SaturationPressure\">Characteristics.Examples.SaturationPressure</a>.

  </html>"),
          experiment(StopTime=0.002),
          Commands(file=
                "Resources/Scripts/Dymola/Subregions.Examples.PhaseChange.Evaporation.mos"
              "Subregions.Examples.PhaseChange.Evaporation.mos"));

      end Evaporation;

      model Hydration
        "<html>Test absorption of H<sub>2</sub>O vapor into the ionomer</html>"
        extends Examples.Subregion(
          'inclSO3-'=true,
          inclN2=true,
          inclH2=false,
          subregion(
            gas(N2(
                consMaterial=ConsThermo.IC,
                consEnergy=ConsThermo.IC,
                p_IC=environment.p)),
            liquid(inclH2O=true, H2O(
                consMaterial=ConsThermo.IC,
                consEnergy=ConsThermo.IC,
                N(stateSelect=StateSelect.always),
                epsilon_IC=0.3)),
            ionomer(
              inclH2O=true,
              'SO3-'(consEnergy=ConsThermo.IC,epsilon=0.3),
              H2O(lambda_IC=8,initEnergy=Init.none))),
          environment(T=333.15*U.K, RH=1));

        annotation (
          Documentation(info="<html><p>The water vapor is held at saturation pressure at the environmental temperature
  Water is supplied as necessary to maintain this condition.  The ionomer begins with hydration of &lambda; = 8 and
  comes to equilibrium at approximately &lambda; &asymp; 14 in about a half an hour.</p>

  <p>See also <a href=\"modelica://FCSys.Characteristics.Examples.HydrationLevel\">Characteristics.Examples.HydrationLevel</a>.</p>

</html>"),
          experiment(StopTime=120),
          Commands(file(ensureTranslated=true) =
              "Resources/Scripts/Dymola/Subregions.Examples.PhaseChange.Hydration.mos"
              "Subregions.Examples.PhaseChange.Hydration.mos"),
          Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                  {100,100}}), graphics));

      end Hydration;

    end PhaseChange;

    package Reactions "Examples of phase change"
      extends Modelica.Icons.ExamplesPackage;

      model HOR "Test the hydrogen oxidation reaction in one subregion"

        output Q.Potential h0=0.5*Characteristics.H2.Gas.Deltah0_f + 0.25*
            Characteristics.O2.Gas.Deltah0_f - 0.5*Characteristics.H2O.Gas.Deltah0_f;

        output Q.Potential w=subregion.graphite.'e-Transfer'.Deltag
          "Overpotential";
        output Q.Current zI=subregion.graphite.'e-Transfer'.I "Reaction rate";
        output Q.Number J_Apercm2=zI*U.cm^2/(subregion.A[Axis.x]*U.A)
          "Electrical current density of the reaction, in A/cm2";
        output Q.Power Qdot=-subregion.graphite.'e-Transfer'.inert.Qdot
          "Rate of heat generation due to reaction";

        extends Examples.Subregion(
          'inclC+'=true,
          'inclSO3-'=true,
          'incle-'=true,
          'inclH+'=true,
          inclH2=true,
          subregion(L={0.287*U.mm,10*U.cm,10*U.cm}),
          environment(T=333.15*U.K, RH=0.8));

        Conditions.ByConnector.BoundaryBus.Single.Sink anBC(graphite(
            'incle-'=true,
            'inclC+'=true,
            redeclare
              Conditions.ByConnector.ThermalDiffusive.Single.Temperature 'C+'(
                set(y=environment.T)),
            'e-'(redeclare function materialSpec =
                  Conditions.ByConnector.Boundary.Single.Material.current,
                materialSet(y=currentSet.y))), gas(inclH2=true, H2(
              materialSet(y=environment.p_dry),
              redeclare function thermalSpec =
                  Conditions.ByConnector.Boundary.Single.Thermal.temperature,
              thermalSet(y=environment.T)))) annotation (Placement(
              transformation(
              extent={{-10,10},{10,-10}},
              rotation=270,
              origin={-24,0})));

        Conditions.ByConnector.BoundaryBus.Single.Sink caBC(ionomer('inclH+'=
                true, 'H+'(redeclare function materialSpec =
                  Conditions.ByConnector.Boundary.Single.Material.potential (
                    redeclare package Data = FCSys.Characteristics.'H+'.Ionomer),
                materialSet(y=0)))) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={24,0})));

        Modelica.Blocks.Sources.Ramp currentSet(
          height=200*U.A,
          duration=20,
          startTime=5)
          annotation (Placement(transformation(extent={{-10,-60},{10,-40}})));

      equation
        connect(subregion.xPositive, caBC.boundary) annotation (Line(
            points={{10,0},{20,6.66134e-016}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect(anBC.boundary, subregion.xNegative) annotation (Line(
            points={{-20,-8.88178e-016},{-10,-8.88178e-016},{-10,0}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));
        annotation (
          experiment(StopTime=30),
          Commands(file=
                "Resources/Scripts/Dymola/Subregions.Examples.Reactions.HOR.mos"
              "Subregions.Examples.Reactions.HOR.mos"),
          Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                  {100,100}}), graphics));
      end HOR;

      model ORR "Test the oxygen reduction reaction in one subregion"

        output Q.Potential w=-subregion.graphite.'e-Transfer'.Deltag
          "Overpotential";
        output Q.Current zI=-subregion.graphite.'e-Transfer'.I "Reaction rate";
        output Q.Number J_Apercm2=zI*U.cm^2/(subregion.A[Axis.x]*U.A)
          "Electrical current density, in A/cm2";
        output Q.Power Qdot=-subregion.graphite.'e-Transfer'.inert.Qdot
          "Rate of heat generation due to reaction";

        extends Examples.Subregion(
          'inclC+'=true,
          'inclSO3-'=true,
          'incle-'=true,
          'inclH+'=true,
          inclH2=false,
          inclH2O=true,
          inclO2=true,
          subregion(L={0.287*U.mm,10*U.cm,10*U.cm}, gas(O2(initEnergy=Init.none))),

          environment(T=333.15*U.K, RH=0.6));

        Conditions.ByConnector.BoundaryBus.Single.Source anBC(ionomer('inclH+'=
                true, 'H+'(
              thermalSet(y=environment.T),
              materialSet(y=0),
              redeclare function materialSpec =
                  Conditions.ByConnector.Boundary.Single.Material.potential (
                    redeclare package Data = FCSys.Characteristics.'H+'.Ionomer))))
          annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=270,
              origin={-24,0})));

        Conditions.ByConnector.BoundaryBus.Single.Source caBC(gas(
            inclH2O=true,
            inclO2=true,
            H2O(
              redeclare function materialSpec =
                  Conditions.ByConnector.Boundary.Single.Material.pressure,
              materialSet(y=environment.p_H2O),
              thermalSet(y=environment.T)),
            O2(
              redeclare function materialSpec =
                  Conditions.ByConnector.Boundary.Single.Material.pressure,
              materialSet(y=environment.p_O2),
              thermalSet(y=environment.T))), graphite(
            'incle-'=true,
            'e-'(materialSet(y=currentSet.y), thermalSet(y=environment.T)),
            'inclC+'=true,
            'C+'(set(y=environment.T)))) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={24,0})));

        Modelica.Blocks.Sources.Ramp currentSet(
          duration=20,
          height=-200*U.A,
          startTime=5)
          annotation (Placement(transformation(extent={{-10,-60},{10,-40}})));
        // offset=-U.mA,

      equation
        connect(subregion.xPositive, caBC.boundary) annotation (Line(
            points={{10,0},{10,8.88178e-016},{20,8.88178e-016}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect(anBC.boundary, subregion.xNegative) annotation (Line(
            points={{-20,-8.88178e-016},{-10,-8.88178e-016},{-10,0}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));
        annotation (
          experiment(StopTime=30),
          Commands(file=
                "Resources/Scripts/Dymola/Subregions.Examples.Reactions.ORR.mos"
              "Subregions.Examples.Reactions.ORR.mos"),
          Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                  {100,100}}), graphics));
      end ORR;

    end Reactions;

    model AirColumn
      "<html>Gas in a vertical array of subregions, affected by gravity</html>"
      import FCSys.Utilities.round;
      extends Modelica.Icons.Example;

      parameter Integer n_y=3
        "Number of discrete subregions along the y axis, besides the 2 boundary subregions"
        annotation (Dialog(__Dymola_label="<html><i>n<i><sub>y</sub></html>"));
      parameter Q.Pressure Deltap_IC=0 "Initial pressure difference"
        annotation (Dialog(__Dymola_label=
              "<html>&Delta;<i>p</i><sub>IC</sub></html>"));

      output Q.Pressure Deltap=subregion2.gas.N2.p - subregion1.gas.N2.p
        "Measured pressure difference";
      output Q.Pressure Deltap_ex=-(n_y + 1)*10*U.m*environment.a[Axis.y]*
          Characteristics.N2.Gas.m*subregions[round(n_y/2)].gas.N2.rho
        "Expected pressure difference";

      parameter Q.NumberAbsolute k=5 "Damping factor";

      inner Conditions.Environment environment(p=U.bar)
        annotation (Placement(transformation(extent={{20,-40},{40,-20}})));
      FCSys.Subregions.Subregion subregion1(
        L={10,10,10}*U.m,
        inclTransX=false,
        inclTransY=true,
        inclTransZ=false,
        gas(inclN2=true, N2(
            zeta=k*FCSys.Characteristics.H2.Gas.zeta(),
            p_IC=environment.p - Deltap_IC/2,
            upstreamY=false)))
        annotation (Placement(transformation(extent={{-10,-40},{10,-20}})));

      FCSys.Subregions.Subregion subregions[n_y](
        each L={10,10,10}*U.m,
        each inclTransZ=false,
        gas(each inclN2=true, N2(
            each zeta=k*FCSys.Characteristics.H2.Gas.zeta(),
            p_IC={environment.p - Deltap_IC/2 - i*Deltap_IC/(n_y + 1) for i in
                1:n_y},
            each upstreamY=false)),
        each inclTransX=false,
        each inclTransY=true) if n_y > 0
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

      FCSys.Subregions.Subregion subregion2(
        L={10,10,10}*U.m,
        inclTransZ=false,
        inclTransX=false,
        inclTransY=true,
        gas(inclN2=true, N2(
            zeta=k*FCSys.Characteristics.H2.Gas.zeta(),
            p_IC=environment.p + Deltap_IC/2,
            upstreamY=false)))
        annotation (Placement(transformation(extent={{-10,20},{10,40}})));

      Conditions.ByConnector.BoundaryBus.Single.Sink BC(gas(inclN2=true, N2(
              materialSet(y=environment.p))))
        annotation (Placement(transformation(extent={{-10,46},{10,66}})));

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

      connect(BC.boundary, subregion2.yPositive) annotation (Line(
          points={{0,52},{0,40}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));
      annotation (
        experiment(StopTime=3, __Dymola_Algorithm="Dassl"),
        Documentation(info="<html><p>This is a model of a vertical column of 10&nbsp;&times;&nbsp;10&nbsp;&times;&nbsp;10&nbsp;m
    regions with N<sub>2</sub> gas.  The upper boundary is held at 1 bar and the lower boundary has zero velocity.

    The initial pressure difference is zero, but a gas enters the upper boundary and a pressure difference
    develops due to gravity.  There are oscillations due to the inertance and compression of the gas.
    After about 1.5&nbsp;s, the pressure difference settles to

      <i>L</i><sub>y</sub> <i>a</i><sub>y</sub> <i>m</i> &rho;
      as expected.</p>

      <p>A temperature gradient is created due to the thermodynamics of the expanding and contracting

      gases.  It takes much longer (over a year!) for the temperatures to equalize due to the size of the system and the low
      thermal conductivity of the gas.  With a stiff solver, the model should simulate at this time scale as well.</p>

      <p>The damping factor (<i>k</i>) can be used to scale the continuity (&zeta;) of the gas in the regions.

      The oscillations are dampened considerably at <i>k</i> = 100.  However, with high values of the factor, the boundary pressures
      are decoupled from the pressures in the region because the nonequilibrium force is considerable.</p>

      <p>Assumptions:</p>
      <ol>
      <li>The central difference scheme is used (no upstream discretization).</ol></html>"),

        Commands(file(ensureTranslated=true) =
            "Resources/Scripts/Dymola/Subregions.Examples.AirColumn.mos"
            "Subregions.Examples.AirColumn.mos"),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics));
    end AirColumn;

    model BinaryDiffusion
      "<html>Pressure-driven flow of H<sub>2</sub>, dragging H<sub>2</sub>O</html>"
      extends Subregion(
        inclH2O=true,
        'inclC+'=true,
        subregion(common(k_Phi={1e-5,1e-5,1e-5}),gas(
            common(k_Phi={1e4,1e4,1e4}),
            H2(T(stateSelect=StateSelect.always), phi(each stateSelect=
                    StateSelect.always)),
            H2O(initEnergy=Init.none, boundaries(Ndot(each stateSelect=
                      StateSelect.always))))),
        environment(RH=0.6));

      Conditions.ByConnector.BoundaryBus.Single.Source source(gas(
          inclH2=true,
          inclH2O=true,
          H2(
            thermalSet(y=environment.T),
            redeclare function materialSpec =
                Conditions.ByConnector.Boundary.Single.Material.pressure,
            redeclare Modelica.Blocks.Sources.Ramp materialSet(
              height=-U.Pa,
              offset=environment.p_dry,
              duration=10)),
          H2O(
            redeclare function materialSpec =
                Conditions.ByConnector.Boundary.Single.Material.pressure,
            materialSet(y=environment.p_H2O),
            thermalSet(y=environment.T)))) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={-24,0})));

      Conditions.ByConnector.BoundaryBus.Single.Source sink(gas(
          inclH2=true,
          inclH2O=true,
          H2(materialSet(y=-source.gas.H2.boundary.Ndot), thermalSet(y=
                  environment.T)),
          H2O(
            redeclare function materialSpec =
                Conditions.ByConnector.Boundary.Single.Material.pressure,
            materialSet(y=environment.p_H2O),
            thermalSet(y=environment.T)))) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={24,0})));

    equation
      connect(source.boundary, subregion.xNegative) annotation (Line(
          points={{-20,-4.44089e-016},{-16,-4.44089e-016},{-16,0},{-10,0}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));
      connect(subregion.xPositive, sink.boundary) annotation (Line(
          points={{10,0},{20,0}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));
      annotation (
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics),
        experiment(StopTime=20, __Dymola_Algorithm="Dassl"),
        Commands(file=
              "Resources/Scripts/Dymola/Subregions.Examples.BinaryDiffusion.mos"
            "Subregions.Examples.BinaryDiffusion.mos"));
    end BinaryDiffusion;

    model Echo "Two regions of gas with an initial pressure difference"
      parameter Q.NumberAbsolute k=1 "Damping factor";

      extends Subregions(subregion1(gas(H2(zeta=k*Characteristics.H2.Gas.zeta()))),
          subregion2(gas(H2(zeta=k*FCSys.Characteristics.H2.Gas.zeta()))));
      annotation (experiment(StopTime=5e-005), Commands(file(ensureTranslated=
                true) = "Resources/Scripts/Dymola/Subregions.Examples.Echo.mos"
            "Subregions.Examples.Echo.mos"));

    end Echo;

    model EchoCentral
      "Two regions of gas with an initial pressure difference, with central difference scheme"
      extends Echo(
        subregion1(gas(
            H2(upstreamX=false),
            H2O(upstreamX=false),
            N2(upstreamX=false),
            O2(upstreamX=false))),
        subregions(gas(
            H2(each upstreamX=false),
            H2O(each upstreamX=false),
            N2(each upstreamX=false),
            O2(each upstreamX=false))),
        subregion2(gas(
            H2(upstreamX=false),
            H2O(upstreamX=false),
            N2(upstreamX=false),
            O2(upstreamX=false))));
      annotation (experiment(StopTime=5e-005), Commands(file(ensureTranslated=
                true) =
            "Resources/Scripts/Dymola/Subregions.Examples.EchoCentral.mos"
            "Subregions.Examples.EchoCentral.mos"));

    end EchoCentral;

    model ElectricalConduction
      "<html>Migration of e<sup>-</sup> through C<sup>+</sup></html>"

      output Q.Potential w=subregion.graphite.'e-'.Deltag[1]
        "Electrical potential";
      output Q.Current zI=-subregion.graphite.'e-'.I[1] "Electrical current";
      output Q.ResistanceElectrical R=w/zI "Measured electrical resistance";
      output Q.ResistanceElectrical R_ex=subregion.L[Axis.x]/(subregion.graphite.
          'e-'.sigma*subregion.A[Axis.x]) "Expected electrical resistance";
      output Q.Power P=subregion.graphite.'e-'.Edot_AT
        "Measured rate of heat generation";
      output Q.Power P_ex=zI^2*R "Expected rate of heat generation";
      output Q.TemperatureAbsolute T=subregion.graphite.'C+'.T
        "Measured temperature";
      output Q.TemperatureAbsolute T_ex=environment.T + subregion.graphite.'C+'.theta
          *U.cm*P/(4*subregion.A[Axis.x]) "Expected temperature";

      extends Examples.Subregion(
        'inclC+'=true,
        'incle-'=true,
        inclH2=false,
        subregion(L={U.cm,U.mm,U.mm}, graphite('C+'(epsilon=1),'e-'(sigma=1e2*U.S
                  /U.m))));

      Conditions.ByConnector.BoundaryBus.Single.Source BC1(graphite(
          'incle-'='incle-',
          'inclC+'=true,
          'C+'(set(y=environment.T)),
          'e-'(
            materialSet(y=-0.05*U.A),
            redeclare function thermalSpec =
                Conditions.ByConnector.Boundary.Single.Thermal.heatRate,
            thermalSet(y=0)))) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={-24,0})));

      Conditions.ByConnector.BoundaryBus.Single.Sink BC2(graphite(
          'inclC+'=true,
          'incle-'='incle-',
          redeclare Conditions.ByConnector.ThermalDiffusive.Single.Temperature
            'C+'(set(y=environment.T)))) annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=90,
            origin={24,0})));

    equation
      connect(BC1.boundary, subregion.xNegative) annotation (Line(
          points={{-20,3.65701e-16},{-16,3.65701e-16},{-16,6.10623e-16},{-10,
              6.10623e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      connect(subregion.xPositive, BC2.boundary) annotation (Line(
          points={{10,6.10623e-16},{16,6.10623e-16},{16,-2.54679e-16},{20,-2.54679e-16}},

          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      annotation (
        experiment(StopTime=40),
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
    </html>"),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics));
    end ElectricalConduction;

    model InternalFlow "Internal, laminar flow of liquid water"
      import FCSys.Utilities.Delta;
      final parameter Q.Area A=subregion.A[Axis.x] "Cross-sectional area"
        annotation (Dialog(__Dymola_label="<html><i>A</i></html>"));

      // Conditions
      parameter Q.VolumeRate Vdot_large=-1.5*U.cc/U.s
        "Prescribed large signal volumetric flow rate";

      // Measurements
      output Q.Pressure Deltap=Delta(subregion.liquid.H2O.boundaries[1, :].p)
        "Measured pressure difference";
      output Q.Length D=2*A/(subregion.L[Axis.y] + subregion.L[Axis.z]);
      output Q.Number Re=subregion.liquid.H2O.phi[Axis.x]*D*subregion.liquid.H2O.eta
          *subregion.liquid.H2O.Data.m*subregion.liquid.H2O.rho if environment.analysis
        "Reynolds number";
      output Q.Pressure Deltap_Poiseuille=-32*subregion.L[Axis.x]*subregion.liquid.H2O.phi[
          Axis.x]/(D^2*subregion.liquid.H2O.eta)
        "Pressure difference according to Poiseuille's law";
      Real x=Deltap_Poiseuille/Deltap;

      output Q.Power Qdot_gen_Poiseuille=-Deltap_Poiseuille*Vdot
        "Rate of heat generation according to Poiseuille's law";
      output Q.VolumeRate Vdot=BC1.liquid.H2O.materialSet.y
        "Total volumetric flow rate";

      extends Examples.Subregion(inclH2=false, subregion(
          L={U.m,U.mm,U.mm},
          inclTransY=true,
          inclTransZ=true,
          liquid(inclH2O=true, H2O(initMaterial=Init.none))));

      Conditions.ByConnector.BoundaryBus.Single.Source BC1(liquid(inclH2O=true,
            H2O(
            redeclare Modelica.Blocks.Sources.Sine materialSet(
              amplitude=0.2*Vdot_large,
              offset=Vdot_large,
              freqHz=1),
            redeclare function materialSpec =
                Conditions.ByConnector.Boundary.Single.Material.volumeRate (
                  redeclare package Data = FCSys.Characteristics.H2O.Liquid),
            thermalSet(y=environment.T)))) annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=270,
            origin={-24,0})));

      Conditions.ByConnector.BoundaryBus.Single.Sink BC2(liquid(inclH2O=true,
            H2O(materialSet(y=U.atm)))) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={24,0})));

      Conditions.ByConnector.BoundaryBus.Single.Source BC3(liquid(inclH2O=true,
            H2O(redeclare function thermalSpec =
                Conditions.ByConnector.Boundary.Single.Thermal.heatRate,
              thermalSet(y=0)))) annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=0,
            origin={0,-24})));

      Conditions.ByConnector.BoundaryBus.Single.Source BC4(liquid(inclH2O=true,
            H2O(redeclare function thermalSpec =
                Conditions.ByConnector.Boundary.Single.Thermal.heatRate,
              thermalSet(y=0)))) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={0,24})));

      Conditions.ByConnector.BoundaryBus.Single.Source BC5(liquid(inclH2O=true,
            H2O(redeclare function thermalSpec =
                Conditions.ByConnector.Boundary.Single.Thermal.heatRate,
              thermalSet(y=0)))) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=315,
            origin={24,24})));

      Conditions.ByConnector.BoundaryBus.Single.Source BC6(liquid(inclH2O=true,
            H2O(redeclare function thermalSpec =
                Conditions.ByConnector.Boundary.Single.Thermal.heatRate,
              thermalSet(y=0)))) annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=315,
            origin={-24,-24})));

    equation
      connect(BC1.boundary, subregion.xNegative) annotation (Line(
          points={{-20,-1.34539e-15},{-16,-1.34539e-15},{-16,6.10623e-16},{-10,
              6.10623e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      connect(BC2.boundary, subregion.xPositive) annotation (Line(
          points={{20,1.23436e-15},{16,1.23436e-15},{16,6.10623e-16},{10,
              6.10623e-16}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      connect(subregion.yNegative, BC3.boundary) annotation (Line(
          points={{6.10623e-016,-10},{6.10623e-016,-20}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));
      connect(BC4.boundary, subregion.yPositive) annotation (Line(
          points={{6.10623e-16,20},{-3.36456e-22,16},{6.10623e-16,16},{
              6.10623e-16,10}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      connect(subregion.zPositive, BC6.boundary) annotation (Line(
          points={{-5,-5},{-21.1716,-21.1716}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      connect(subregion.zNegative, BC5.boundary) annotation (Line(
          points={{5,5},{21.1716,21.1716}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));
      annotation (
        Documentation(info="<html><p>A small-signal variation is added to the time-average flow rate in order to demonstrate the effects of inertance.</p>

    <p>Note that the temperature increases due to viscous dissipation, but the increase is limited by thermal convection.  The walls are adiabatic.</p></html>"),

        experiment(
          StopTime=5,
          Tolerance=1e-008,
          __Dymola_Algorithm="Dassl"),
        Commands(file=
              "Resources/Scripts/Dymola/Subregions.Examples.InternalFlow.mos"
            "Subregions.Examples.InternalFlow.mos"),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics));
    end InternalFlow;

    model Subregion
      "<html>Single subregion, with H<sub>2</sub> by default</html>"

      extends Modelica.Icons.Example;

      parameter Boolean 'inclC+'=false
        "<html>Carbon plus (C<sup>+</sup>)</html>" annotation (choices(
            __Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean 'inclSO3-'=false
        "<html>Nafion sulfonate (C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S<sup>-</sup>, abbreviated as SO<sub>3</sub><sup>-</sup>)</html>"
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

      inner Conditions.Environment environment(RH=0)
        annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

      FCSys.Subregions.Subregion subregion(
        L={1,1,1}*U.cm,
        inclTransY=false,
        inclTransZ=false,
        graphite(final 'inclC+'='inclC+', final 'incle-'='incle-'),
        ionomer(final 'inclSO3-'='inclSO3-', final 'inclH+'='inclH+'),
        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2))
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

      annotation (
        Documentation(info="<html><p>This model is boring.  It just establishes a
  single subregion with H<sub>2</sub> by default.  There are no boundary conditions
  other than those implied by the open connectors (no diffusion current, no forces,

  no thermal conduction).  Other examples in this package extend from this one.</p>
  </html>"),
        experiment(StopTime=10),
        Commands(file(ensureTranslated=true) =
            "Resources/Scripts/Dymola/Subregions.Examples.Subregion.mos"
            "Subregions.Examples.Subregion.mos"));

    end Subregion;

    model Subregions
      "<html>Horizontal array of subregions with an initial pressure difference (H<sub>2</sub> included by default)</html>"
      extends Modelica.Icons.Example;

      parameter Integer n_x=0
        "Number of discrete subregions along the x axis, besides the 2 side subregions"
        annotation (Dialog(__Dymola_label="<html><i>n<i><sub>x</sub></html>"));
      parameter Q.Pressure Deltap_IC=-100*U.Pa "Initial pressure difference"
        annotation (Dialog(__Dymola_label=
              "<html>&Delta;<i>p</i><sub>IC</sub></html>"));
      parameter Boolean 'inclC+'=false
        "<html>Carbon plus (C<sup>+</sup>)</html>" annotation (choices(
            __Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));
      parameter Boolean 'inclSO3-'=false
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

      inner Conditions.Environment environment
        annotation (Placement(transformation(extent={{-10,-30},{10,-10}})));
      FCSys.Subregions.Subregion subregion1(
        L={1,1,1}*U.cm,
        inclTransY=false,
        inclTransZ=false,
        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2,
          H2(p_IC=environment.p - Deltap_IC/2, phi(each stateSelect=StateSelect.always,
                each fixed=true)),
          H2O(p_IC=environment.p - Deltap_IC/2, phi(each stateSelect=
                  StateSelect.always, each fixed=true)),
          N2(p_IC=environment.p - Deltap_IC/2, phi(each stateSelect=StateSelect.always,
                each fixed=true)),
          O2(p_IC=environment.p - Deltap_IC/2, phi(each stateSelect=StateSelect.always,
                each fixed=true))),
        graphite(final 'inclC+'='inclC+', final 'incle-'='incle-'),
        ionomer(final 'inclSO3-'='inclSO3-', final 'inclH+'='inclH+'),
        liquid(H2O(epsilon_IC=0.25)))
        annotation (Placement(transformation(extent={{-30,-10},{-10,10}})));

      FCSys.Subregions.Subregion subregions[n_x](
        each L={1,1,1}*U.cm,
        each inclTransY=false,
        each inclTransZ=false,
        gas(
          each final inclH2=inclH2,
          each final inclH2O=inclH2O,
          each final inclN2=inclN2,
          each final inclO2=inclO2,
          H2(p_IC={environment.p - Deltap_IC/2 - i*Deltap_IC/(n_x + 1) for i
                 in 1:n_x}, each phi(each stateSelect=StateSelect.always, each
                fixed=true)),
          H2O(p_IC={environment.p - Deltap_IC/2 - i*Deltap_IC/(n_x + 1) for i
                 in 1:n_x}, each phi(each stateSelect=StateSelect.always, each
                fixed=true)),
          N2(p_IC={environment.p - Deltap_IC/2 - i*Deltap_IC/(n_x + 1) for i
                 in 1:n_x}, each phi(each stateSelect=StateSelect.always, each
                fixed=true)),
          O2(p_IC={environment.p - Deltap_IC/2 - i*Deltap_IC/(n_x + 1) for i
                 in 1:n_x}, each phi(each stateSelect=StateSelect.always, each
                fixed=true))),
        graphite(each final 'inclC+'='inclC+', each final 'incle-'='incle-'),
        ionomer(each final 'inclSO3-'='inclSO3-', each final 'inclH+'='inclH+'),

        liquid(H2O(each epsilon_IC=0.25))) if n_x > 0
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

      FCSys.Subregions.Subregion subregion2(
        L={1,1,1}*U.cm,
        inclTransY=false,
        inclTransZ=false,
        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2,
          H2(p_IC=environment.p + Deltap_IC/2),
          H2O(p_IC=environment.p + Deltap_IC/2),
          N2(p_IC=environment.p + Deltap_IC/2),
          O2(p_IC=environment.p + Deltap_IC/2)),
        graphite(final 'inclC+'='inclC+', final 'incle-'='incle-'),
        ionomer(final 'inclSO3-'='inclSO3-', final 'inclH+'='inclH+'),
        liquid(H2O(each epsilon_IC=0.25)))
        annotation (Placement(transformation(extent={{10,-10},{30,10}})));

    equation
      connect(subregion1.xPositive, subregions[1].xNegative) annotation (Line(
          points={{-10,6.10623e-016},{-10,6.10623e-016}},
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
          points={{10,6.10623e-016},{10,6.10623e-016}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));

      annotation (experiment(
          StopTime=2,
          Tolerance=1e-006,
          __Dymola_Algorithm="Dassl"), Commands(file(ensureTranslated=true) =
            "Resources/Scripts/Dymola/Subregions.Examples.Subregions.mos"
            "Subregions.Examples.Subregions.mos"));
    end Subregions;

    model ThermalConduction "Thermal conduction (through solid)"
      extends Examples.Subregions(
        n_x=6,
        'inclC+'=true,
        inclH2=false,
        subregion1(graphite('C+'(T_IC=environment.T + 30*U.K, epsilon=1))),
        subregions(each graphite('C+'(epsilon=1))),
        subregion2(graphite('C+'(epsilon=1))));

      annotation (Commands(file=
              "Resources/Scripts/Dymola/Subregions.Examples.ThermalConduction.mos"
            "Subregions.Examples.ThermalConduction.mos"), experiment(StopTime=
              500, __Dymola_Algorithm="Dassl"));

    end ThermalConduction;

    model ThermalConductionConvection
      "Thermal conduction through a solid alongside conduction and convection through a gas"

      extends ThermalConduction(
        inclN2=true,
        subregion1(gas(N2(
              T_IC=subregion1.graphite.'C+'.T_IC,
              phi(stateSelect=StateSelect.always, displayUnit="mm/s"),
              phi_boundaries(each displayUnit="mm/s"))), graphite('C+'(epsilon=
                  0.5))),
        subregions(
          common(each k_Phi={10,10,10}),
          gas(N2(each phi(each stateSelect=StateSelect.always, displayUnit=
                    "mm/s"), each phi_boundaries(each displayUnit="mm/s"))),
          each graphite('C+'(epsilon=0.5))),
        subregion2(gas(N2(phi(displayUnit="mm/s"), phi_boundaries(each
                  displayUnit="mm/s"))), graphite('C+'(epsilon=0.5))),
        environment(psi_O2_dry=0, RH=0));

      annotation (Commands(file=
              "Resources/Scripts/Dymola/Subregions.Examples.ThermalConductionConvection.mos"
            "Subregions.Examples.ThermalConductionConvection.mos"), experiment(
            StopTime=400, __Dymola_Algorithm="Dassl"));

    end ThermalConductionConvection;

  end Examples;
  extends Modelica.Icons.Package;

  model Subregion "Subregion with all phases"
    import Modelica.Constants.inf;

    extends PartialSubregion(final n_spec=gas.n_spec + graphite.n_spec +
          ionomer.n_spec + liquid.n_spec);

    FCSys.Phases.Gas gas(
      n_inter=2,
      final n_trans=n_trans,
      final k_inter_Phi={common.k_Phi[cartTrans],gasLiq.k_Phi[cartTrans]},
      final k_inter_Q={common.k_Q,gasLiq.k_Q}) "Gas" annotation (Dialog(group=
            "Phases (click to edit)"),Placement(transformation(extent={{-30,-22},
              {-10,-2}})));

    FCSys.Phases.Graphite graphite(
      n_inter=1,
      final n_trans=n_trans,
      'incle-Transfer'=inclHOR or inclORR,
      final k_inter_Phi={common.k_Phi[cartTrans]},
      final k_inter_Q={common.k_Q}) "Graphite" annotation (Dialog(group=
            "Phases (click to edit)"), Placement(transformation(extent={{10,-22},
              {30,-2}})));

    FCSys.Phases.Ionomer ionomer(
      n_inter=1,
      final n_trans=n_trans,
      final k_inter_Phi={common.k_Phi[cartTrans]},
      final k_inter_Q={common.k_Q}) "Ionomer" annotation (Dialog(group=
            "Phases (click to edit)"), Placement(transformation(extent={{50,-22},
              {70,-2}})));

    FCSys.Phases.Liquid liquid(
      n_inter=2,
      final n_trans=n_trans,
      final k_inter_Phi={common.k_Phi[cartTrans],gasLiq.k_Phi[cartTrans]},
      final k_inter_Q={common.k_Q,gasLiq.k_Q}) "Liquid" annotation (Dialog(
          group="Phases (click to edit)"), Placement(transformation(extent={{-70,
              -22},{-50,-2}})));

    Chemistry.HOR HOR(final n_trans=n_trans) if inclHOR
      "Hydrogen oxidation reaction"
      annotation (Placement(transformation(extent={{88,-30},{108,-10}})));
    Chemistry.ORR ORR(final n_trans=n_trans) if inclORR
      "Oxygen reduction reaction"
      annotation (Placement(transformation(extent={{88,-46},{108,-26}})));

    // Independence factors
    Phases.ExchangeParams common(k_Phi={1.7,1.7,1.7}) "Among all phases"
      annotation (Dialog(group="Independence factors"));
    Phases.ExchangeParams gasLiq(k_Phi={inf,inf,inf}) "Between gas and liquid"
      annotation (Dialog(group="Independence factors"));

    Connectors.BoundaryBus xNegative if inclTransX
      "Negative boundary along the x axis" annotation (Placement(transformation(
            extent={{-120,-40},{-100,-20}}), iconTransformation(extent={{-110,-10},
              {-90,10}})));
    Connectors.BoundaryBus yNegative if inclTransY
      "Negative boundary along the y axis" annotation (Placement(transformation(
            extent={{-96,-64},{-76,-44}}), iconTransformation(extent={{-10,-110},
              {10,-90}})));
    Connectors.BoundaryBus zNegative if inclTransZ
      "Negative boundary along the z axis" annotation (Placement(transformation(
            extent={{88,8},{108,28}}), iconTransformation(extent={{40,40},{60,
              60}})));
    Connectors.BoundaryBus xPositive if inclTransX
      "Positive boundary along the x axis" annotation (Placement(transformation(
            extent={{100,-4},{120,16}}), iconTransformation(extent={{90,-10},{
              110,10}})));
    Connectors.BoundaryBus yPositive if inclTransY
      "Positive boundary along the y axis" annotation (Placement(transformation(
            extent={{76,20},{96,40}}), iconTransformation(extent={{-10,90},{10,
              110}})));
    Connectors.BoundaryBus zPositive if inclTransZ
      "Positive boundary along the z axis" annotation (Placement(transformation(
            extent={{-108,-52},{-88,-32}}), iconTransformation(extent={{-60,-60},
              {-40,-40}})));

    Chemistry.CapillaryVolume volume(
      final V=V,
      final inclGas=gas.n_spec > 0,
      final inclLiquid=liquid.n_spec > 0,
      final inclSolid=graphite.n_spec + ionomer.n_spec > 0)
      "Volume with capillary pressure included" annotation (Dialog, Placement(
          transformation(extent={{-24,-80},{-4,-60}})));

  protected
    final parameter Boolean inclHOR=graphite.'incle-' and ionomer.'inclH+' and
        gas.inclH2 "Include the hydrogen oxidation reaction";
    final parameter Boolean inclORR=graphite.'incle-' and ionomer.'inclH+' and
        gas.inclO2 and (gas.inclH2O or liquid.inclH2O)
      "Include the oxygen reduction reaction";

    outer Conditions.Environment environment "Environmental conditions";
    Connectors.InertNode exchCommon "Connector for exchange among all species"
      annotation (HideResult=true, Placement(transformation(extent={{76,32},{96,
              52}}), iconTransformation(extent={{100,18},{120,38}})));
    Connectors.InertNode exchGasLiq
      "Connector for exchange between gas and liquid" annotation (HideResult=
          true, Placement(transformation(extent={{76,44},{96,64}}),
          iconTransformation(extent={{100,18},{120,38}})));

  equation
    // Boundaries
    // ----------
    // Gas
    connect(gas.yNegative, yNegative.gas) annotation (Line(
        points={{-20,-20.4},{-20,-54},{-86,-54}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.zPositive, zPositive.gas) annotation (Line(
        points={{-28,-20},{-28,-42},{-98,-42}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.xNegative, xNegative.gas) annotation (Line(
        points={{-28,-12},{-36,-12},{-36,-30},{-110,-30}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.yPositive, yPositive.gas) annotation (Line(
        points={{-20,-2},{-20,30},{86,30}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.zNegative, zNegative.gas) annotation (Line(
        points={{-15,-7},{-12,-4},{-12,18},{98,18}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.xPositive, xPositive.gas) annotation (Line(
        points={{-12,-12},{-4,-12},{-4,6},{110,6}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    // Graphite
    connect(graphite.yNegative, yNegative.graphite) annotation (Line(
        points={{20,-20.4},{20,-54},{-86,-54}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.yPositive, yPositive.graphite) annotation (Line(
        points={{20,-2},{20,30},{86,30}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.zNegative, zNegative.graphite) annotation (Line(
        points={{25,-7},{28,-4},{28,18},{98,18}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.zPositive, zPositive.graphite) annotation (Line(
        points={{12,-20},{12,-42},{-98,-42}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.xNegative, xNegative.graphite) annotation (Line(
        points={{12,-12},{4,-12},{4,-30},{-110,-30}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.xPositive, xPositive.graphite) annotation (Line(
        points={{28,-12},{36,-12},{36,6},{110,6}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    // Ionomer
    connect(ionomer.yNegative, yNegative.ionomer) annotation (Line(
        points={{60,-20.4},{60,-54},{-86,-54}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(ionomer.yPositive, yPositive.ionomer) annotation (Line(
        points={{60,-2},{60,30},{86,30}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(ionomer.zNegative, zNegative.ionomer) annotation (Line(
        points={{65,-7},{68,-4},{68,18},{98,18}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(ionomer.zPositive, zPositive.ionomer) annotation (Line(
        points={{52,-20},{52,-42},{-98,-42}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(ionomer.xNegative, xNegative.ionomer) annotation (Line(
        points={{52,-12},{44,-12},{44,-30},{-110,-30}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(ionomer.xPositive, xPositive.ionomer) annotation (Line(
        points={{68,-12},{76,-12},{76,6},{110,6}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    // Liquid
    connect(liquid.yNegative, yNegative.liquid) annotation (Line(
        points={{-60,-20.4},{-60,-54},{-86,-54}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(liquid.yPositive, yPositive.liquid) annotation (Line(
        points={{-60,-2},{-60,30},{86,30}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(liquid.zNegative, zNegative.liquid) annotation (Line(
        points={{-55,-7},{-52,-4},{-52,18},{98,18}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(liquid.zPositive, zPositive.liquid) annotation (Line(
        points={{-68,-20},{-68,-42},{-98,-42}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(liquid.xNegative, xNegative.liquid) annotation (Line(
        points={{-68,-12},{-76,-12},{-76,-30},{-110,-30}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(liquid.xPositive, xPositive.liquid) annotation (Line(
        points={{-52,-12},{-44,-12},{-44,6},{110,6}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    // Mixing
    connect(gas.dalton, volume.gas) annotation (Line(
        points={{-12,-20},{-12,-68}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(liquid.amagat, volume.liquid) annotation (Line(
        points={{-52,-20},{-52,-71},{-15,-71}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(volume.solid, graphite.amagat) annotation (Line(
        points={{-8,-72},{28,-72},{28,-20}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(volume.solid, ionomer.amagat) annotation (Line(
        points={{-8,-72},{68,-72},{68,-20}},
        color={47,107,251},
        smooth=Smooth.None));

    // Inert exchange
    // --------------
    // Common
    connect(gas.inter[1], exchCommon.node) annotation (Line(
        points={{-25,-6.5},{-25,42},{86,42}},
        color={221,23,47},
        smooth=Smooth.None));
    connect(graphite.inter[1], exchCommon.node) annotation (Line(
        points={{15,-7},{15,42},{86,42}},
        color={221,23,47},
        smooth=Smooth.None));
    connect(ionomer.inter[1], exchCommon.node) annotation (Line(
        points={{55,-7},{55,42},{86,42}},
        color={221,23,47},
        smooth=Smooth.None));
    connect(liquid.inter[1], exchCommon.node) annotation (Line(
        points={{-65,-6.5},{-65,42},{86,42}},
        color={221,23,47},
        smooth=Smooth.None));
    // Gas-liquid
    connect(gas.inter[2], exchGasLiq.node) annotation (Line(
        points={{-25,-7.5},{-25,54},{86,54}},
        color={221,23,47},
        smooth=Smooth.None));
    connect(liquid.inter[2], exchGasLiq.node) annotation (Line(
        points={{-65,-7.5},{-65,54},{86,54}},
        color={221,23,47},
        smooth=Smooth.None));

    // Reactions and phase change (not shown in diagram)
    // -------------------------------------------------
    connect(gas.chemH2[1], HOR.chemH2);
    connect(gas.chemO2[1], ORR.chemO2);
    connect(liquid.chemH2O[1], ORR.chemH2O);
    if liquid.inclH2O then
      if gas.inclH2O then
        connect(liquid.chemH2O[2], gas.chemH2O[2]);
      end if;
      if ionomer.inclH2O then
        connect(liquid.chemH2O[3], ionomer.chemH2O[1]);
      end if;
    elseif gas.inclH2O then
      connect(gas.chemH2O[1], ORR.chemH2O);
    end if;
    connect(graphite.'cheme-'[1], HOR.'cheme-');
    connect(graphite.'cheme-'[1], ORR.'cheme-');
    connect(ionomer.'chemH+'[1], HOR.'chemH+');
    connect(ionomer.'chemH+'[1], ORR.'chemH+');
    annotation (Documentation(info="<html><p>Assumptions:</p><ol>
<li>The oxygen reduction reaction generates liquid water if it is included; otherwise,
it generates H<sub>2</sub>O vapor.  Since phase change is a dynamic, nonequilibrium
process, there is a difference.</li></ol>

   <p>Please see the documentation of the
   <a href=\"modelica://FCSys.Subregions.PartialSubregion\">PartialSubregion</a> model.</p></html>"),
        Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-120,-80},{
              120,60}}), graphics={Text(
              extent={{78,-44},{118,-50}},
              lineColor={127,127,127},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid,
              textString="(connections not shown
in diagram)")}));
  end Subregion;

  model SubregionIonomer "Subregion with only the ionomer phase"

    extends PartialSubregion(final n_spec=ionomer.n_spec);

    FCSys.Phases.Ionomer ionomer(n_inter=0, final n_trans=n_trans) "Ionomer"
      annotation (Dialog(group="Phases (click to edit)"), Placement(
          transformation(extent={{-20,0},{0,20}})));

    Connectors.BoundaryBus xNegative if inclTransX
      "Negative boundary along the x axis" annotation (Placement(transformation(
            extent={{-60,0},{-40,20}}), iconTransformation(extent={{-110,-10},{
              -90,10}})));
    Connectors.BoundaryBus yNegative if inclTransY
      "Negative boundary along the y axis" annotation (Placement(transformation(
            extent={{-20,-40},{0,-20}}), iconTransformation(extent={{-10,-110},
              {10,-90}})));
    Connectors.BoundaryBus zNegative if inclTransZ
      "Negative boundary along the z axis" annotation (Placement(transformation(
            extent={{0,20},{20,40}}), iconTransformation(extent={{40,40},{60,60}})));
    Connectors.BoundaryBus xPositive if inclTransX
      "Positive boundary along the x axis" annotation (Placement(transformation(
            extent={{20,0},{40,20}}), iconTransformation(extent={{90,-10},{110,
              10}})));
    Connectors.BoundaryBus yPositive if inclTransY
      "Positive boundary along the y axis" annotation (Placement(transformation(
            extent={{-20,40},{0,60}}), iconTransformation(extent={{-10,90},{10,
              110}})));
    Connectors.BoundaryBus zPositive if inclTransZ
      "Positive boundary along the z axis" annotation (Placement(transformation(
            extent={{-40,-20},{-20,0}}), iconTransformation(extent={{-60,-60},{
              -40,-40}})));

    Chemistry.CapillaryVolume volume(
      final V=V,
      final inclSolid=ionomer.n_spec > 0,
      final inclGas=false,
      final inclLiquid=false) "Volume with capillary pressure included"
      annotation (Dialog, Placement(transformation(extent={{0,-20},{20,0}})));

  equation
    // Boundaries
    // ----------
    // Ionomer

    connect(ionomer.yNegative, yNegative.ionomer) annotation (Line(
        points={{-10,1.6},{-10,-30}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    connect(ionomer.yPositive, yPositive.ionomer) annotation (Line(
        points={{-10,20},{-10,50}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(ionomer.zNegative, zNegative.ionomer) annotation (Line(
        points={{-5,15},{10,30}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(ionomer.zPositive, zPositive.ionomer) annotation (Line(
        points={{-18,2},{-30,-10}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    connect(ionomer.xNegative, xNegative.ionomer) annotation (Line(
        points={{-18,10},{-50,10}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(ionomer.xPositive, xPositive.ionomer) annotation (Line(
        points={{-2,10},{30,10}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    connect(volume.solid, ionomer.amagat) annotation (Line(
        points={{16,-12},{20,-12},{20,2},{-2,2}},
        color={47,107,251},
        smooth=Smooth.None));
    annotation (
      defaultComponentName="subregion",
      Documentation(info="<html>
   <p>Please see the documentation of the
   <a href=\"modelica://FCSys.Subregions.PartialSubregion\">PartialSubregion</a> model.</p></html>"),

      Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-60,-40},{40,
              60}}), graphics));
  end SubregionIonomer;

  model SubregionNoIonomer "Subregion with all phases except ionomer"
    import Modelica.Constants.inf;

    extends PartialSubregion(final n_spec=gas.n_spec + graphite.n_spec + liquid.n_spec);

    FCSys.Phases.Gas gas(
      n_inter=2,
      final n_trans=n_trans,
      final k_inter_Phi={common.k_Phi[cartTrans],gasLiq.k_Phi[cartTrans]},
      final k_inter_Q={common.k_Q,gasLiq.k_Q}) "Gas" annotation (Dialog(group=
            "Phases (click to edit)"),Placement(transformation(extent={{-10,-22},
              {10,-2}})));

    FCSys.Phases.Graphite graphite(
      n_inter=1,
      final n_trans=n_trans,
      final k_inter_Phi={common.k_Phi[cartTrans]},
      final k_inter_Q={common.k_Q}) "Graphite" annotation (Dialog(group=
            "Phases (click to edit)"), Placement(transformation(extent={{30,-22},
              {50,-2}})));

    FCSys.Phases.Liquid liquid(
      n_inter=2,
      final n_trans=n_trans,
      final k_inter_Phi={common.k_Phi[cartTrans],gasLiq.k_Phi[cartTrans]},
      final k_inter_Q={common.k_Q,gasLiq.k_Q}) "Liquid" annotation (Dialog(
          group="Phases (click to edit)"), Placement(transformation(extent={{-50,
              -22},{-30,-2}})));

    // Independence factors
    Phases.ExchangeParams common(k_Phi={1.7,1.7,1.7}) "Among all phases"
      annotation (Dialog(group="Independence factors"));
    Phases.ExchangeParams gasLiq(k_Phi={inf,inf,inf}) "Between gas and liquid"
      annotation (Dialog(group="Independence factors"));

    Connectors.BoundaryBus xNegative if inclTransX
      "Negative boundary along the x axis" annotation (Placement(transformation(
            extent={{-100,-40},{-80,-20}}), iconTransformation(extent={{-110,-10},
              {-90,10}})));
    Connectors.BoundaryBus yNegative if inclTransY
      "Negative boundary along the y axis" annotation (Placement(transformation(
            extent={{-76,-64},{-56,-44}}), iconTransformation(extent={{-10,-110},
              {10,-90}})));
    Connectors.BoundaryBus zNegative if inclTransZ
      "Negative boundary along the z axis" annotation (Placement(transformation(
            extent={{68,8},{88,28}}), iconTransformation(extent={{40,40},{60,60}})));
    Connectors.BoundaryBus xPositive if inclTransX
      "Positive boundary along the x axis" annotation (Placement(transformation(
            extent={{80,-4},{100,16}}), iconTransformation(extent={{90,-10},{
              110,10}})));
    Connectors.BoundaryBus yPositive if inclTransY
      "Positive boundary along the y axis" annotation (Placement(transformation(
            extent={{56,20},{76,40}}), iconTransformation(extent={{-10,90},{10,
              110}})));
    Connectors.BoundaryBus zPositive if inclTransZ
      "Positive boundary along the z axis" annotation (Placement(transformation(
            extent={{-88,-52},{-68,-32}}), iconTransformation(extent={{-60,-60},
              {-40,-40}})));

    Chemistry.CapillaryVolume volume(
      final V=V,
      final inclGas=gas.n_spec > 0,
      final inclLiquid=liquid.n_spec > 0,
      final inclSolid=graphite.n_spec > 0)
      "Volume with capillary pressure included"
      annotation (Placement(transformation(extent={{-4,-80},{16,-60}})));

  protected
    outer Conditions.Environment environment "Environmental conditions";

    // Exchange
    Connectors.InertNode exchCommon "Among all phases" annotation (HideResult=
          true,Placement(transformation(extent={{56,32},{76,52}}),
          iconTransformation(extent={{100,18},{120,38}})));
    Connectors.InertNode exchGasLiq "Between gas and liquid" annotation (
        HideResult=true, Placement(transformation(extent={{56,44},{76,64}}),
          iconTransformation(extent={{100,18},{120,38}})));

  equation
    // Boundaries
    // ----------
    // Gas
    connect(gas.yNegative, yNegative.gas) annotation (Line(
        points={{0,-20.4},{0,-54},{-66,-54}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.zPositive, zPositive.gas) annotation (Line(
        points={{-8,-20},{-8,-42},{-78,-42}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.xNegative, xNegative.gas) annotation (Line(
        points={{-8,-12},{-16,-12},{-16,-30},{-90,-30}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.yPositive, yPositive.gas) annotation (Line(
        points={{0,-2},{0,30},{66,30}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.zNegative, zNegative.gas) annotation (Line(
        points={{5,-7},{8,-4},{8,18},{78,18}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.xPositive, xPositive.gas) annotation (Line(
        points={{8,-12},{16,-12},{16,6},{90,6}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    // Graphite
    connect(graphite.yNegative, yNegative.graphite) annotation (Line(
        points={{40,-20.4},{40,-54},{-66,-54}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.yPositive, yPositive.graphite) annotation (Line(
        points={{40,-2},{40,30},{66,30}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.zNegative, zNegative.graphite) annotation (Line(
        points={{45,-7},{48,-4},{48,18},{78,18}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.zPositive, zPositive.graphite) annotation (Line(
        points={{32,-20},{32,-42},{-78,-42}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.xNegative, xNegative.graphite) annotation (Line(
        points={{32,-12},{24,-12},{24,-30},{-90,-30}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.xPositive, xPositive.graphite) annotation (Line(
        points={{48,-12},{56,-12},{56,6},{90,6}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    // Liquid
    connect(liquid.yNegative, yNegative.liquid) annotation (Line(
        points={{-40,-20.4},{-40,-54},{-66,-54}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(liquid.yPositive, yPositive.liquid) annotation (Line(
        points={{-40,-2},{-40,30},{66,30}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(liquid.zNegative, zNegative.liquid) annotation (Line(
        points={{-35,-7},{-32,-4},{-32,18},{78,18}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(liquid.zPositive, zPositive.liquid) annotation (Line(
        points={{-48,-20},{-48,-42},{-78,-42}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(liquid.xNegative, xNegative.liquid) annotation (Line(
        points={{-48,-12},{-56,-12},{-56,-30},{-90,-30}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(liquid.xPositive, xPositive.liquid) annotation (Line(
        points={{-32,-12},{-24,-12},{-24,6},{90,6}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    // Mixing
    connect(liquid.amagat, volume.liquid) annotation (Line(
        points={{-32,-20},{-32,-71},{5,-71}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(gas.dalton, volume.gas) annotation (Line(
        points={{8,-20},{8,-68}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(volume.solid, graphite.amagat) annotation (Line(
        points={{12,-72},{48,-72},{48,-20}},
        color={47,107,251},
        smooth=Smooth.None));

    // Inert exchange
    // --------------
    // Common
    connect(gas.inter[1], exchCommon.node) annotation (Line(
        points={{-5,-6.5},{-5,42},{66,42}},
        color={221,23,47},
        smooth=Smooth.None));
    connect(graphite.inter[1], exchCommon.node) annotation (Line(
        points={{35,-7},{35,42},{66,42}},
        color={221,23,47},
        smooth=Smooth.None));
    connect(liquid.inter[1], exchCommon.node) annotation (Line(
        points={{-45,-6.5},{-45,42},{66,42}},
        color={221,23,47},
        smooth=Smooth.None));
    // Gas-liquid
    connect(gas.inter[2], exchGasLiq.node) annotation (Line(
        points={{-5,-7.5},{-5,54},{66,54}},
        color={221,23,47},
        smooth=Smooth.None));
    connect(liquid.inter[2], exchGasLiq.node) annotation (Line(
        points={{-45,-7.5},{-45,54},{66,54}},
        color={221,23,47},
        smooth=Smooth.None));

    // Phase change (not shown in diagram)
    // -----------------------------------
    if gas.inclH2O and liquid.inclH2O then
      connect(gas.chemH2O[2], liquid.chemH2O[2]);
    end if;

    annotation (
      defaultComponentName="subregion",
      Documentation(info="<html>
   <p>Please see the documentation of the
   <a href=\"modelica://FCSys.Subregions.PartialSubregion\">PartialSubregion</a> model.</p></html>"),

      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-80},{
              100,60}}), graphics));
  end SubregionNoIonomer;

  partial model PartialSubregion
    "Base model for multi-dimensional, multi-species storage, transport, and exchange"
    import Modelica.Math.BooleanVectors.countTrue;
    import Modelica.Math.BooleanVectors.enumerate;
    import Modelica.Math.BooleanVectors.index;
    // extends FCSys.Icons.Names.Top3;

    // Geometric parameters
    inner parameter Q.Length L[Axis](each min=Modelica.Constants.small) = {U.cm,
      U.cm,U.cm} "Lengths" annotation (Dialog(group="Geometry", __Dymola_label=
            "<html><b><i>L</i></b></html>"));

    // Assumptions
    // -----------
    // Included boundaries
    parameter Boolean inclTransX=true "X" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        tab="Assumptions",
        group="Included transport axes",
        compact=true));
    parameter Boolean inclTransY=true "Y" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        tab="Assumptions",
        group="Included transport axes",
        compact=true));
    parameter Boolean inclTransZ=true "Z" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        tab="Assumptions",
        group="Included transport axes",
        compact=true));

    // Auxiliary variables (for analysis)
    final inner parameter Q.Volume V=product(L) "Volume";
    final parameter Q.Area A[Axis]=fill(V, 3) ./ L "Cross-sectional areas";

  protected
    parameter Integer n_spec(start=0) "Number of species"
      annotation (HideResult=true);
    final inner parameter Boolean inclTrans[Axis]={inclTransX,inclTransY,
        inclTransZ} "true, if each pair of boundaries is included";
    final inner parameter Boolean inclRot[Axis]={inclTransY and inclTransZ,
        inclTransZ and inclTransX,inclTransX and inclTransY}
      "true, if each axis of rotation has all its tangential boundaries included";
    final parameter Integer n_trans=countTrue(inclTrans)
      "Number of transport axes";
    final inner parameter Integer cartRot[:]=index(inclRot)
      "Cartesian-axis indices of the components of rotational momentum";
    final inner parameter Integer cartTrans[n_trans]=index(inclTrans)
      "Cartesian-axis indices of the transport axes";
    final inner parameter Integer transCart[Axis]=enumerate(inclTrans)
      "Transport-axis indices of the Cartesian axes";

    annotation (
      defaultComponentPrefixes="replaceable",
      defaultComponentName="subregion",
      Documentation(info="<html>
  <p>At least one component of translational momentum must be included.
  All of the components are included by default.</p>

  <p>This model should be extended to include the appropriate phases, reactions, etc.</p>
  </html>"),
      Icon(graphics={Line(
              points={{-100,0},{-40,0}},
              color={127,127,127},
              thickness=0.5,
              visible=inclTransX,
              smooth=Smooth.None),Line(
              points={{0,-40},{0,-100}},
              color={127,127,127},
              thickness=0.5,
              visible=inclTransY,
              smooth=Smooth.None),Line(
              points={{40,40},{50,50}},
              color={127,127,127},
              thickness=0.5,
              visible=inclTransZ,
              smooth=Smooth.None),Polygon(
              points={{-40,16},{-16,40},{40,40},{40,-16},{16,-40},{-40,-40},{-40,
              16}},
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
              visible=inclTransX,
              smooth=Smooth.None,
              thickness=0.5),Line(
              points={{0,28},{0,-40}},
              color={210,210,210},
              visible=inclTransY,
              smooth=Smooth.None,
              thickness=0.5),Line(
              points={{28,0},{100,0}},
              color={127,127,127},
              thickness=0.5,
              visible=inclTransX,
              smooth=Smooth.None),Line(
              points={{0,100},{0,28}},
              color={127,127,127},
              thickness=0.5,
              visible=inclTransY,
              smooth=Smooth.None),Line(
              points={{-12,-12},{40,40}},
              color={210,210,210},
              visible=inclTransZ,
              smooth=Smooth.None,
              thickness=0.5),Line(
              points={{-40,16},{16,16},{16,-40}},
              color={127,127,127},
              smooth=Smooth.None),Line(
              points={{-50,-50},{-12,-12}},
              color={127,127,127},
              thickness=0.5,
              visible=inclTransZ,
              smooth=Smooth.None),Polygon(
              points={{-40,16},{-16,40},{40,40},{40,-16},{16,-40},{-40,-40},{-40,
              16}},
              lineColor={127,127,127},
              smooth=Smooth.None),Line(
              points={{40,40},{16,16}},
              color={127,127,127},
              smooth=Smooth.None),Text(
              extent={{-100,56},{100,96}},
              textString="%name",
              lineColor={0,0,0})}),
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
              100,100}}), graphics));

  end PartialSubregion;
  annotation (Documentation(info="
<html>
  <p><b>Licensed by the Hawaii Natural Energy Institute under the Modelica License 2</b><br>
Copyright &copy; 2007&ndash;2014, <a href=\"http://www.hnei.hawaii.edu/\">Hawaii Natural Energy Institute</a> and <a href=\"http://www.gtrc.gatech.edu/\">Georgia Tech Research Corporation</a>.</p>

<p><i>This Modelica package is <u>free</u> software and the use is completely at <u>your own risk</u>;
it can be redistributed and/or modified under the terms of the Modelica License 2. For license conditions (including the
disclaimer of warranty) see <a href=\"modelica://FCSys.UsersGuide.License\">
FCSys.UsersGuide.License</a> or visit <a href=\"http://www.modelica.org/licenses/ModelicaLicense2\">
http://www.modelica.org/licenses/ModelicaLicense2</a>.</i></p>
</html>"));

end Subregions;

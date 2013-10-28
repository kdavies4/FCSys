within FCSys;
package Subregions
  "Control volumes with multi-species transport, exchange, and storage"
  extends Modelica.Icons.Package;

  package Examples "Examples"
    extends Modelica.Icons.ExamplesPackage;
    // TODO: In the documentation of each model:
    //   1. Insert the sample plots or link to the sample results of the User's Guide.
    //   2. Add discussion from the dissertation.
    package PhaseChange "Examples of phase change"
      extends Modelica.Icons.ExamplesPackage;
      model Evaporation
        "<html>Evaporation of H<sub>2</sub>O liquid into sub-saturated vapor</html>"

        output Q.Pressure p_sat=
            Modelica.Media.Air.MoistAir.saturationPressureLiquid(subregion.gas.H2O.T
            /U.K)*U.Pa "Saturation pressure via Modelica.Media";

        extends Examples.Subregion(
          inclH2O=true,
          inclH2=false,
          subregion(gas(H2O(
                p_IC=U.kPa,
                consEnergy=ConsThermo.dynamic,
                N(stateSelect=StateSelect.always))), liquid(inclH2O=inclH2O,
                H2O(
                consEnergy=ConsThermo.dynamic,
                N(stateSelect=StateSelect.always),
                epsilon_IC=0.001))));

        annotation (
          Documentation(info="<html><p>Initially, the water vapor is below saturation and a small amount of liquid water is present (1/1000 of the total volume).
  Some of the liquid evaporates until saturation is reached. The boundaries are adiabatic; therefore, the temperature of the liquid and the gas 
  decreases due to the enthalpy of formation.</p></html>"),
          experiment(StopTime=120),
          Commands(file=
                "Resources/Scripts/Dymola/Subregions.Examples.PhaseChange.Evaporation.mos"
              "Subregions.Examples.PhaseChange.Evaporation.mos"),
          Diagram(graphics),
          __Dymola_experimentSetupOutput);

      end Evaporation;

      model Hydration
        "<html>Test absorption and desorption of H<sub>2</sub>O between the gas and ionomer</html>"

        extends Examples.Subregion(
          'inclSO3-'=true,
          inclH2O=true,
          inclH2=false,
          subregion(gas(H2O(
                initEnergy=Init.temperature,
                p_IC=1.001*U.atm,
                N(stateSelect=StateSelect.always))), ionomer(
              reduceThermal=true,
              inclH2O=inclH2O,
              H2O(N(stateSelect=StateSelect.always)))));
        // In Dymola 7.4, p_IC=1.1*environment.p has no effect on the
        // initial pressure, but p_IC=1.1*U.atm does.
        annotation (experiment(StopTime=0.003, Tolerance=1e-06), Commands(file(
                ensureTranslated=true) =
              "Resources/Scripts/Dymola/Subregions.Examples.PhaseChange.Hydration.mos"
              "Subregions.Examples.PhaseChange.Hydration.mos"));

      end Hydration;

      model SaturationPressure
        "<html>Evaluate the saturation pressure curve of H<sub>2</sub>O by varying temperature</html>"
        extends SaturationPressureIdeal(subregion(gas(H2O(redeclare package
                  Data = FCSys.Characteristics.H2O.Gas))));

        annotation (
          Commands(file=
                "Resources/Scripts/Dymola/Subregions.Examples.PhaseChange.SaturationPressure.mos"
              "Subregions.Examples.PhaseChange.SaturationPressure.mos"),
          experiment(StopTime=3600, Tolerance=1e-006),
          __Dymola_experimentSetupOutput);

      end SaturationPressure;

      model SaturationPressureIdeal
        "<html>Evaluate the saturation pressure curve of H<sub>2</sub>O by varying temperature, assuming ideal gas</html>"
        import saturationPressureSI =
          Modelica.Media.Air.MoistAir.saturationPressureLiquid;

        output Q.Pressure p_sat=saturationPressureSI(subregion.gas.H2O.T/U.K)*U.Pa
          "Saturation pressure via Modelica.Media";
        output Q.Number T_degC=U.to_degC(subregion.gas.H2O.T)
          "Temperature in degree Celsius";

        extends Examples.Subregion(
          inclH2O=true,
          inclH2=false,
          subregion(gas(H2O(N(stateSelect=StateSelect.always),p_IC=
                    saturationPressureSI(environment.T/U.K)*U.Pa)), liquid(H2O(
                  N(stateSelect=StateSelect.always), epsilon_IC=0.25), inclH2O=
                  true)),
          environment(T=274.15*U.K));

        FCSys.Conditions.ByConnector.FaceBus.Single.Flows BC1(liquid(inclH2O=
                true, H2O(redeclare function thermalSpec =
                  FCSys.Conditions.ByConnector.Face.Single.Thermal.temperature,
                redeclare Modelica.Blocks.Sources.Ramp thermalSet(
                height=99*U.K,
                duration=3600,
                offset=environment.T,
                y))), gas(inclH2O=inclH2O, H2O(redeclare function thermalSpec
                = FCSys.Conditions.ByConnector.Face.Single.Thermal.temperature,
                redeclare Modelica.Blocks.Sources.Ramp thermalSet(
                height=99*U.K,
                duration=3600,
                offset=environment.T,
                y)))) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=90,
              origin={-24,0})));

        FCSys.Conditions.ByConnector.FaceBus.Single.Flows BC2(gas(inclH2O=
                inclH2O, H2O(redeclare function thermalSpec =
                  FCSys.Conditions.ByConnector.Face.Single.Thermal.temperature,
                redeclare Modelica.Blocks.Sources.Ramp thermalSet(
                height=99*U.K,
                duration=3600,
                offset=environment.T,
                y))), liquid(inclH2O=true, H2O(redeclare function thermalSpec
                = FCSys.Conditions.ByConnector.Face.Single.Thermal.temperature,
                redeclare Modelica.Blocks.Sources.Ramp thermalSet(
                height=99*U.K,
                duration=3600,
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
          experiment(StopTime=3600, Tolerance=1e-006),
          Commands(file(ensureTranslated=true) =
              "Resources/Scripts/Dymola/Subregions.Examples.PhaseChange.SaturationPressureIdeal.mos"
              "Subregions.Examples.PhaseChange.SaturationPressureIdeal.mos"),
          __Dymola_experimentSetupOutput);
      end SaturationPressureIdeal;
    end PhaseChange;

    package Reactions "Examples of phase change"
      extends Modelica.Icons.ExamplesPackage;

      model HOR "Test the hydrogen oxidation reaction in one subregion"

        output Q.Potential wprime=subregion.graphite.'e-'.wprime
          "Overpotential";
        output Q.Current zI=-subregion.graphite.'e-'.faces[1, 1].Ndot
          "Electrical current";
        output Q.Number zJ_Apercm2=zI*U.cm^2/(subregion.A[Axis.x]*U.A)
          "Electrical current density, in A/cm2";
        output Q.Power Qdot=-subregion.graphite.'e-'.Edot_DE
          "Rate of heat generation";
        output Q.Power P=wprime*zI "Electrical power";

        extends Examples.Subregion(
          'inclC+'=true,
          'inclSO3-'=true,
          'incle-'=true,
          'inclH+'=true,
          inclH2=true,
          subregion(
            L={0.287*U.mm,10*U.cm,10*U.cm},
            gas(
              reduceThermal=true,
              T(stateSelect=StateSelect.always),
              H2(phi(each fixed=true, each stateSelect=StateSelect.always))),
            graphite(T(stateSelect=StateSelect.always)),
            dielectric(final A=0),
            ionomer(
              reduceThermal=true,
              T(stateSelect=StateSelect.always),
              'H+'(consTransX=ConsMom.steady))));

        // Note:  The double layer capacitance introduces nonlinear equations in Dymola 2014.
        // A = 0 eliminates it.

        Conditions.ByConnector.FaceBus.Single.Efforts negativeBC(graphite(
              'incle-'=true, 'e-'(redeclare function materialSpec =
                  FCSys.Conditions.ByConnector.Face.Single.Material.current,
                redeclare Modelica.Blocks.Sources.Ramp materialSet(height=200*U.A,
                  duration=20))), gas(inclH2=true)) annotation (Placement(
              transformation(
              extent={{-10,10},{10,-10}},
              rotation=270,
              origin={-24,0})));

        Conditions.ByConnector.FaceBus.Single.Efforts positiveBC(ionomer(
              'inclH+'=true)) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={24,0})));

      equation
        connect(subregion.xPositive, positiveBC.face) annotation (Line(
            points={{10,0},{20,6.66134e-016}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect(negativeBC.face, subregion.xNegative) annotation (Line(
            points={{-20,-8.88178e-016},{-10,-8.88178e-016},{-10,0}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));
        annotation (
          experiment(StopTime=30, Tolerance=1e-006),
          Commands(file=
                "Resources/Scripts/Dymola/Subregions.Examples.Reactions.HOR.mos"
              "Subregions.Examples.Reactions.HOR.mos"),
          __Dymola_experimentSetupOutput,
          Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                  {100,100}}), graphics));
      end HOR;

      model ORR "Test the oxygen reduction reaction in one subregion"

        output Q.Potential wprime=-subregion.graphite.'e-'.wprime
          "Overpotential";
        output Q.Current zI=subregion.graphite.'e-'.faces[1, 2].Ndot
          "Electrical current";
        output Q.Number zJ_Apercm2=zI*U.cm^2/(subregion.A[Axis.x]*U.A)
          "Electrical current density, in A/cm2";
        output Q.Power Qdot=-subregion.graphite.'e-'.Edot_DE
          "Rate of heat generation";
        output Q.Power P=wprime*zI "Electrical power";

        extends Examples.Subregion(
          'inclC+'=true,
          'inclSO3-'=true,
          'incle-'=true,
          'inclH+'=true,
          inclH2=false,
          inclH2O=true,
          inclO2=true,
          subregion(
            L={0.287*U.mm,10*U.cm,10*U.cm},
            gas(
              reduceThermal=true,
              T(stateSelect=StateSelect.always),
              H2O(phi(each fixed=true, each stateSelect=StateSelect.always)),
              O2(phi(each fixed=true, each stateSelect=StateSelect.always),
                  initEnergy=Init.none)),
            graphite(T(stateSelect=StateSelect.always)),
            dielectric(final A=0),
            ionomer(
              reduceThermal=true,
              T(stateSelect=StateSelect.always),
              'H+'(consTransX=ConsMom.steady))));

        // Note:  The double layer capacitance introduces nonlinear equations in Dymola 2014.
        // A = 0 eliminates it.

        Conditions.ByConnector.FaceBus.Single.Efforts negativeBC(ionomer(
              'inclH+'=true)) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={24,0})));
        Conditions.ByConnector.FaceBus.Single.Efforts positiveBC(graphite(
              'incle-'=true, 'e-'(redeclare function materialSpec =
                  FCSys.Conditions.ByConnector.Face.Single.Material.current,
                redeclare Modelica.Blocks.Sources.Ramp materialSet(height=-200*
                    U.A, duration=20))), gas(inclH2O=true, inclO2=true))
          annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=270,
              origin={-24,0})));

      equation
        connect(subregion.xPositive, positiveBC.face) annotation (Line(
            points={{10,0},{-20,6.66134e-016}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));
        connect(negativeBC.face, subregion.xNegative) annotation (Line(
            points={{20,-8.88178e-016},{-10,-8.88178e-016},{-10,0}},
            color={127,127,127},
            thickness=0.5,
            smooth=Smooth.None));
        annotation (
          experiment(StopTime=30, Tolerance=1e-006),
          Commands(file=
                "Resources/Scripts/Dymola/Subregions.Examples.Reactions.ORR.mos"
              "Subregions.Examples.Reactions.ORR.mos"),
          __Dymola_experimentSetupOutput,
          Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                  {100,100}}), graphics));
      end ORR;
    end Reactions;

    model AirColumn
      "<html>Vertical array of subregions with an initial pressure difference (C<sup>+</sup> and N<sub>2</sub> included by default)</html>"
      import FCSys.Utilities.round;
      extends Modelica.Icons.Example;

      parameter Integer n_y=3
        "Number of discrete subregions along the y axis, besides the 2 boundary subregions"
        annotation (Dialog(__Dymola_label="<html><i>n<i><sub>y</sub></html>"));
      parameter Q.Pressure Deltap_IC=0 "Initial pressure difference"
        annotation (Dialog(__Dymola_label=
              "<html>&Delta;<i>p</i><sub>IC</sub></html>"));
      parameter Boolean 'inclC+'=true
        "<html>Carbon plus (C<sup>+</sup>)</html>" annotation (choices(
            __Dymola_checkBox=true), Dialog(group="Species",
            __Dymola_descriptionLabel=true));

      output Q.Pressure Deltap=subregion2.gas.N2.p - subregion1.gas.N2.p
        "Measured pressure difference";
      output Q.Pressure Deltap_ex=-(n_y + 1)*10*U.m*environment.a[Axis.y]*
          Characteristics.N2.Gas.m*subregions[round(n_y/2)].gas.N2.rho
        "Expected pressure difference";

      inner Conditions.Environment environment(p=U.bar)
        annotation (Placement(transformation(extent={{20,-40},{40,-20}})));
      FCSys.Subregions.Subregion subregion1(
        L={10,10,10}*U.m,
        inclTransX=false,
        inclTransY=true,
        inclTransZ=false,
        k_common=1e-7,
        gas(inclN2=true, N2(
            p_IC=environment.p - Deltap_IC/2,
            N(stateSelect=StateSelect.always),
            upstreamY=false,
            phi(each stateSelect=StateSelect.always, each fixed=true))),
        graphite('inclC+'='inclC+', 'C+'(epsilon=1e-9, T(stateSelect=
                  StateSelect.always))))
        annotation (Placement(transformation(extent={{-10,-40},{10,-20}})));

      FCSys.Subregions.Subregion subregions[n_y](
        each L={10,10,10}*U.m,
        each inclTransZ=false,
        each k_common=1e-7,
        gas(each inclN2=true, N2(
            p_IC={environment.p - Deltap_IC/2 - i*Deltap_IC/(n_y + 1) for i in
                1:n_y},
            each N(stateSelect=StateSelect.always),
            each upstreamY=false,
            each phi(each stateSelect=StateSelect.always,each fixed=true))),
        graphite(each 'inclC+'='inclC+', each 'C+'(epsilon=1e-9,T(stateSelect=
                  StateSelect.always))),
        each inclTransX=false,
        each inclTransY=true) if n_y > 0
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

      FCSys.Subregions.Subregion subregion2(
        L={10,10,10}*U.m,
        inclTransZ=false,
        inclTransX=false,
        inclTransY=true,
        each k_common=1e-7,
        graphite('inclC+'='inclC+', 'C+'(epsilon=1e-9,T(stateSelect=StateSelect.always))),

        gas(inclN2=true, N2(
            p_IC=environment.p + Deltap_IC/2,
            N(stateSelect=StateSelect.always),
            upstreamY=false,
            phi(each stateSelect=StateSelect.always,each fixed=true))))
        annotation (Placement(transformation(extent={{-10,20},{10,40}})));

      Conditions.ByConnector.FaceBus.Single.Efforts BC1(gas(inclN2=true, N2(
            redeclare function thermalSpec =
                FCSys.Conditions.ByConnector.Face.Single.Thermal.heatRate,
            materialSet(y=environment.p),
            thermalSet(y=0))))
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

      connect(BC1.face, subregion2.yPositive) annotation (Line(
          points={{0,52},{0,40}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));
      annotation (
        experiment(
          StopTime=8,
          Tolerance=1e-005,
          __Dymola_Algorithm="Dassl"),
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
        Commands(file(ensureTranslated=true) =
            "Resources/Scripts/Dymola/Subregions.Examples.AirColumn.mos"
            "Subregions.Examples.AirColumn.mos"),
        __Dymola_experimentSetupOutput,
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics));
    end AirColumn;

    model Echo
      "Two regions of gas with initial pressure difference, no dampening"
      extends Subregions(inclH2=true, 'inclC+'=false);
      annotation (
        experiment(StopTime=0.0001),
        Commands(file(ensureTranslated=true) =
            "Resources/Scripts/Dymola/Subregions.Examples.Echo.mos"
            "Subregions.Examples.Echo.mos"),
        __Dymola_experimentSetupOutput);

    end Echo;

    model EchoCentral
      "Two regions of gas with initial pressure difference, no dampening, with central difference scheme"
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
      annotation (
        experiment(StopTime=0.0001),
        Commands(file(ensureTranslated=true) =
            "Resources/Scripts/Dymola/Subregions.Examples.EchoCentral.mos"
            "Subregions.Examples.EchoCentral.mos"),
        __Dymola_experimentSetupOutput);
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
          *U.cm*P/(4*U.mm^2) "Expected temperature";

      extends Examples.Subregion(
        'inclC+'=true,
        'incle-'=true,
        inclH2=false,
        subregion(L={U.cm,U.mm,U.mm}, graphite('C+'(epsilon=1),'e-'(sigma=1e2*U.S
                  /U.m))));

      FCSys.Conditions.ByConnector.FaceBus.Single.Flows BC1(graphite(
          final 'incle-'='incle-',
          'inclC+'=true,
          'e-'(materialSet(y=-0.05*U.A)),
          redeclare FCSys.Conditions.ByConnector.ThermalDiffusion.Temperature
            'C+'(source(y=environment.T)))) annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={-24,0})));

      FCSys.Conditions.ByConnector.FaceBus.Single.Flows BC2(graphite(
          'inclC+'=true,
          final 'incle-'='incle-',
          'e-'(redeclare function materialSpec =
                FCSys.Conditions.ByConnector.Face.Single.Material.pressure),
          redeclare FCSys.Conditions.ByConnector.ThermalDiffusion.Temperature
            'C+'(source(y=environment.T)))) annotation (Placement(
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
    </html>"),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics));
    end ElectricalConduction;

    model InternalFlow "Internal, laminar flow of liquid water"
      import FCSys.Utilities.Delta;

      final parameter Q.Area A=subregion.A[Axis.x] "Cross-sectional area"
        annotation (Dialog(__Dymola_label="<html><i>A</i></html>"));

      // Conditions
      parameter Q.VolumeRate Vdot=-0.1*U.L/U.s
        "Prescribed large signal volumetric flow rate"
        annotation (Dialog(__Dymola_label="<html><i>V&#775;</i></html>"));

      // Measurements
      output Q.Pressure Deltap=Delta(subregion.liquid.H2O.faces[1, :].p)
        "Measured pressure difference";
      output Q.Length D=2*A/(subregion.L[Axis.y] + subregion.L[Axis.z]);
      output Q.Number Re=subregion.liquid.H2O.phi[Axis.x]*D*subregion.liquid.H2O.eta
          *subregion.liquid.H2O.Data.m*subregion.liquid.H2O.rho if environment.analysis
        "Reynolds number";
      output Q.Pressure Deltap_Poiseuille=-32*subregion.L[Axis.x]*subregion.liquid.H2O.phi[
          Axis.x]/(D^2*subregion.liquid.H2O.eta)
        "Pressure difference according to Poiseuille's law";
      output Q.Power Qdot_gen=subregion.liquid.H2O.Edot_AT
        "Rate of heat generation";
      output Q.Power Qdot_gen_Poiseuille=-Deltap_Poiseuille*subregion.liquid.H2O.Vdot[
          1] "Rate of heat generation according to Poiseuille's law";

      extends Examples.Subregion(inclH2=false, subregion(
          L={U.m,U.cm,U.cm},
          inclTransY=true,
          inclTransZ=true,
          liquid(inclH2O=true, H2O(final epsilon_IC=1))));

      Conditions.ByConnector.FaceBus.Single.Flows BC1(liquid(inclH2O=true, H2O(
              redeclare Modelica.Blocks.Sources.Sine materialSet(
              amplitude=0.2*Vdot,
              offset=Vdot,
              freqHz=0.01), redeclare function materialSpec =
                FCSys.Conditions.ByConnector.Face.Single.Material.volumeRate (
                  redeclare package Data = FCSys.Characteristics.H2O.Liquid))))
        annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=270,
            origin={-24,0})));

      Conditions.ByConnector.FaceBus.Single.Flows BC2(liquid(inclH2O=true, H2O(
              redeclare function materialSpec =
                FCSys.Conditions.ByConnector.Face.Single.Material.pressure,
              materialSet(y=U.atm)))) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={24,0})));

      Conditions.ByConnector.FaceBus.Single.Flows BC3(liquid(inclH2O=true, H2O(
            redeclare function afterSpec =
                FCSys.Conditions.ByConnector.Face.Single.Translational.velocity,

            redeclare function beforeSpec =
                FCSys.Conditions.ByConnector.Face.Single.Translational.velocity,

            redeclare function thermalSpec =
                FCSys.Conditions.ByConnector.Face.Single.Thermal.temperature,
            thermalSet(y=environment.T)))) annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=0,
            origin={0,-24})));

      Conditions.ByConnector.FaceBus.Single.Flows BC4(liquid(inclH2O=true, H2O(
            redeclare function afterSpec =
                FCSys.Conditions.ByConnector.Face.Single.Translational.velocity,

            redeclare function beforeSpec =
                FCSys.Conditions.ByConnector.Face.Single.Translational.velocity,

            redeclare function thermalSpec =
                FCSys.Conditions.ByConnector.Face.Single.Thermal.temperature,
            thermalSet(y=environment.T)))) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={0,24})));

      Conditions.ByConnector.FaceBus.Single.Flows BC5(liquid(inclH2O=true, H2O(
            redeclare function afterSpec =
                FCSys.Conditions.ByConnector.Face.Single.Translational.velocity,

            redeclare function beforeSpec =
                FCSys.Conditions.ByConnector.Face.Single.Translational.velocity,

            redeclare function thermalSpec =
                FCSys.Conditions.ByConnector.Face.Single.Thermal.temperature,
            thermalSet(y=environment.T)))) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=315,
            origin={24,24})));

      Conditions.ByConnector.FaceBus.Single.Flows BC6(liquid(inclH2O=true, H2O(
            redeclare function afterSpec =
                FCSys.Conditions.ByConnector.Face.Single.Translational.velocity,

            redeclare function beforeSpec =
                FCSys.Conditions.ByConnector.Face.Single.Translational.velocity,

            redeclare function thermalSpec =
                FCSys.Conditions.ByConnector.Face.Single.Thermal.temperature,
            thermalSet(y=environment.T)))) annotation (Placement(transformation(
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
          points={{6.10623e-016,-10},{6.10623e-016,-20}},
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
        Documentation(info="<html><p>Note that the temperature increases due to viscous dissipation.  
        However, the temperature rise is limited because the walls are held at constant temperature.</p></html>"),

        experiment(
          StopTime=800,
          Tolerance=1e-008,
          __Dymola_Algorithm="Dassl"),
        Commands(file=
              "Resources/Scripts/Dymola/Subregions.Examples.InternalFlow.mos"
            "Subregions.Examples.InternalFlow.mos"),
        __Dymola_experimentSetupOutput,
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

      inner Conditions.Environment environment
        annotation (Placement(transformation(extent={{20,20},{40,40}})));

      FCSys.Subregions.Subregion subregion(
        L={1,1,1}*U.cm,
        inclTransY=false,
        inclTransZ=false,
        inclTransX=true,
        graphite(
          final 'inclC+'='inclC+',
          final 'incle-'='incle-',
          'C+'(epsilon=0.25)),
        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2,
          H2(N(stateSelect=StateSelect.always))),
        liquid(H2O(epsilon_IC=0.25)),
        ionomer(
          final 'inclSO3-'='inclSO3-',
          final 'inclH+'='inclH+',
          'SO3-'(epsilon=0.25)))
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

      annotation (
        Documentation(info="<html><p>This model is boring.  It just sets up a
  single subregion with H<sub>2</sub> by default.  There are no boundary conditions
  other than those implied by the open connectors (no diffusion current, no forces, 
  no thermal conduction).  Other examples in this package are extended from this one.</p>
  </html>"),
        experiment(StopTime=10),
        Commands(file(ensureTranslated=true) =
            "Resources/Scripts/Dymola/Subregions.Examples.Subregion.mos"
            "Subregions.Examples.Subregion.mos"));

    end Subregion;

    model Subregions
      "<html>Horizontal array of subregions with an initial pressure difference (C<sup>+</sup> and H<sub>2</sub> included by default)</html>"
      extends Modelica.Icons.Example;

      parameter Integer n_x=0
        "Number of discrete subregions along the x axis, besides the 2 side subregions"
        annotation (Dialog(__Dymola_label="<html><i>n<i><sub>x</sub></html>"));
      parameter Q.Pressure Deltap_IC=-100*U.Pa "Initial pressure difference"
        annotation (Dialog(__Dymola_label=
              "<html>&Delta;<i>p</i><sub>IC</sub></html>"));
      parameter Boolean 'inclC+'=true
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
        annotation (Placement(transformation(extent={{20,20},{40,40}})));
      FCSys.Subregions.Subregion subregion1(
        L={1,1,1}*U.cm,
        inclTransY=false,
        inclTransZ=false,
        gas(
          reduceTrans=true,
          reduceThermal=true,
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2,
          H2O(
            p_IC=environment.p - Deltap_IC/2,
            N(stateSelect=StateSelect.always),
            T(stateSelect=StateSelect.always),
            phi(each stateSelect=StateSelect.always,each fixed=true)),
          N2(
            p_IC=environment.p - Deltap_IC/2,
            N(stateSelect=StateSelect.always),
            T(stateSelect=StateSelect.always),
            phi(each stateSelect=StateSelect.always,each fixed=true)),
          O2(
            p_IC=environment.p - Deltap_IC/2,
            N(stateSelect=StateSelect.always),
            T(stateSelect=StateSelect.always),
            phi(each stateSelect=StateSelect.always,each fixed=true)),
          H2(
            p_IC=environment.p - Deltap_IC/2,
            N(stateSelect=StateSelect.always),
            T(stateSelect=StateSelect.always),
            phi(each stateSelect=StateSelect.always,each fixed=true))),
        graphite(
          final 'inclC+'='inclC+',
          final 'incle-'='incle-',
          'C+'(epsilon=0.001,T(stateSelect=StateSelect.always))),
        ionomer(
          final 'inclSO3-'='inclSO3-',
          final 'inclH+'='inclH+',
          'SO3-'(epsilon=0.001)),
        liquid(H2O(epsilon_IC=0.25)))
        annotation (Placement(transformation(extent={{-40,-10},{-20,10}})));

      FCSys.Subregions.Subregion subregions[n_x](
        each L={1,1,1}*U.cm,
        each inclTransY=false,
        each inclTransZ=false,
        graphite(
          each final 'inclC+'='inclC+',
          each final 'incle-'='incle-',
          'C+'(each epsilon=0.001)),
        ionomer(
          each final 'inclSO3-'='inclSO3-',
          each final 'inclH+'='inclH+',
          'SO3-'(each epsilon=0.001)),
        gas(
          each reduceThermal=true,
          each reduceTrans=true,
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
        liquid(H2O(each epsilon_IC=0.25))) if n_x > 0
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

      FCSys.Subregions.Subregion subregion2(
        L={1,1,1}*U.cm,
        inclTransY=false,
        inclTransZ=false,
        graphite(
          final 'inclC+'='inclC+',
          final 'incle-'='incle-',
          'C+'(epsilon=0.001,T(stateSelect=StateSelect.always))),
        ionomer(
          final 'inclSO3-'=false,
          final 'inclH+'='inclH+',
          'SO3-'(epsilon=0.001)),
        gas(
          reduceTrans=true,
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2,
          H2O(p_IC=environment.p + Deltap_IC/2),
          N2(p_IC=environment.p + Deltap_IC/2),
          O2(p_IC=environment.p + Deltap_IC/2),
          H2(
            p_IC=environment.p + Deltap_IC/2,
            N(stateSelect=StateSelect.always),
            T(stateSelect=StateSelect.always))))
        annotation (Placement(transformation(extent={{20,-10},{40,10}})));

      FCSys.Conditions.ByConnector.FaceBus.Single.Flows BC1(
        ionomer(final 'inclSO3-'='inclSO3-', final 'inclH+'='inclH+'),
        graphite(final 'inclC+'='inclC+', final 'incle-'='incle-'),
        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2)) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={-56,0})));

      FCSys.Conditions.ByConnector.FaceBus.Single.Flows BC2(
        ionomer(final 'inclSO3-'='inclSO3-', final 'inclH+'='inclH+'),
        gas(
          final inclH2=inclH2,
          final inclH2O=inclH2O,
          final inclN2=inclN2,
          final inclO2=inclO2),
        graphite(final 'inclC+'='inclC+', final 'incle-'='incle-')) annotation
        (Placement(transformation(
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
        experiment(
          StopTime=4,
          Tolerance=1e-006,
          __Dymola_Algorithm="Dassl"),
        Commands(file(ensureTranslated=true) =
            "Resources/Scripts/Dymola/Subregions.Examples.Subregions.mos"
            "Subregions.Examples.Subregions.mos"),
        __Dymola_experimentSetupOutput);
    end Subregions;

    model ThermalConduction "Thermal conduction (through solid)"
      extends Examples.Subregions(
        n_x=6,
        'inclC+'=true,
        inclH2=false,
        subregion1(graphite('C+'(T_IC=environment.T + 30*U.K))),
        subregions(graphite('C+'(each T(stateSelect=StateSelect.always)))),
        subregion2(graphite('C+'(s(stateSelect=StateSelect.never)))));

      annotation (
        Commands(file=
              "Resources/Scripts/Dymola/Subregions.Examples.ThermalConduction.mos"
            "Subregions.Examples.ThermalConduction.mos"),
        experiment(StopTime=500, Algorithm="Dassl"),
        __Dymola_experimentSetupOutput);

    end ThermalConduction;

    model ThermalConductionConvection
      "Thermal conduction through a solid in parallel with conduction and convection through a gas"

      extends Examples.Subregions(
        n_x=6,
        'inclC+'=true,
        inclH2=false,
        inclN2=true,
        subregion1(gas(N2(
              T_IC=environment.T + 30*U.K,
              p_IC=environment.p,
              N(stateSelect=StateSelect.always),
              T(stateSelect=StateSelect.always),
              phi(stateSelect=StateSelect.always,displayUnit="mm/s"))),
            graphite('C+'(
              T_IC=environment.T + 30*U.K,
              epsilon=0.5,
              T(stateSelect=StateSelect.always)))),
        subregions(each gas(N2(
              p_IC=environment.p,
              N(stateSelect=StateSelect.always),
              T(stateSelect=StateSelect.always),
              phi(each stateSelect=StateSelect.always,displayUnit="mm/s"))),
            each graphite('C+'(epsilon=0.5, T(stateSelect=StateSelect.always)))),

        subregion2(gas(N2(
              p_IC=environment.p,
              phi(displayUnit="mm/s"),
              N(stateSelect=StateSelect.always),
              T(stateSelect=StateSelect.always))), graphite('C+'(epsilon=0.5, T(
                  stateSelect=StateSelect.always)))));

      annotation (
        Commands(file=
              "Resources/Scripts/Dymola/Subregions.Examples.ThermalConductionConvection.mos"
            "Subregions.Examples.ThermalConductionConvection.mos"),
        experiment(
          StopTime=500,
          Tolerance=1e-005,
          __Dymola_Algorithm="Dassl"),
        __Dymola_experimentSetupOutput);

    end ThermalConductionConvection;

  end Examples;

  model Subregion "Subregion with all phases"
    import Modelica.Constants.eps;
    extends Partial(final n_spec=gas.n_spec + graphite.n_spec + ionomer.n_spec
           + liquid.n_spec);

    parameter Q.NumberAbsolute k_gas_liq=eps
      "Additional coupling factor between gas and liquid" annotation (Dialog(
          group="Geometry",__Dymola_label=
            "<html><i>k</i><sub>gas liq</sub></html>"));
    parameter Q.NumberAbsolute k_graphite_liq=eps
      "Additional coupling factor between graphite and liquid" annotation (
        Dialog(group="Geometry", __Dymola_label=
            "<html><i>k</i><sub>graphite liq</sub></html>"));

    FCSys.Phases.Gas gas(
      n_inter=2,
      final n_trans=n_trans,
      k_inter={k_common,k_gas_liq}) "Gas" annotation (Dialog(group=
            "Phases (click to edit)"), Placement(transformation(extent={{-90,-22},
              {-70,-2}})));

    FCSys.Phases.Graphite graphite(
      n_inter=2,
      final n_trans=n_trans,
      k_inter={k_common,k_graphite_liq}) "Graphite" annotation (Dialog(group=
            "Phases (click to edit)"), Placement(transformation(extent={{-50,-22},
              {-30,-2}})));

    FCSys.Phases.Ionomer ionomer(
      n_inter=1,
      final n_trans=n_trans,
      k_inter={k_common}) "Ionomer" annotation (Dialog(group=
            "Phases (click to edit)"), Placement(transformation(extent={{10,-22},
              {30,-2}})));

    FCSys.Phases.Liquid liquid(
      n_inter=3,
      final n_trans=n_trans,
      k_inter={k_common,k_gas_liq,k_graphite_liq}) "Liquid" annotation (Dialog(
          group="Phases (click to edit)"), Placement(transformation(extent={{50,
              -22},{70,-2}})));

    inner Reactions.HOR HOR(n_trans=n_trans) if graphite.'incle-' and ionomer.
      'inclH+' and gas.inclH2 "Hydrogen oxidation reaction"
      annotation (Placement(transformation(extent={{80,-30},{100,-10}})));
    inner Reactions.ORR ORR(final n_trans=n_trans) if graphite.'incle-' and
      ionomer.'inclH+' and gas.inclO2 and gas.inclH2O
      "Oxygen reduction reaction"
      annotation (Placement(transformation(extent={{80,-46},{100,-26}})));

    Connectors.FaceBus xNegative if inclTransX "Negative face along the x axis"
      annotation (Placement(transformation(extent={{-140,-40},{-120,-20}}),
          iconTransformation(extent={{-110,-10},{-90,10}})));
    Connectors.FaceBus yNegative if inclTransY "Negative face along the y axis"
      annotation (Placement(transformation(extent={{-116,-64},{-96,-44}}),
          iconTransformation(extent={{-10,-110},{10,-90}})));
    Connectors.FaceBus zNegative if inclTransZ "Negative face along the z axis"
      annotation (Placement(transformation(extent={{88,8},{108,28}}),
          iconTransformation(extent={{40,40},{60,60}})));
    Connectors.FaceBus xPositive if inclTransX "Positive face along the x axis"
      annotation (Placement(transformation(extent={{100,-4},{120,16}}),
          iconTransformation(extent={{90,-10},{110,10}})));
    Connectors.FaceBus yPositive if inclTransY "Positive face along the y axis"
      annotation (Placement(transformation(extent={{76,20},{96,40}}),
          iconTransformation(extent={{-10,90},{10,110}})));
    Connectors.FaceBus zPositive if inclTransZ "Positive face along the z axis"
      annotation (Placement(transformation(extent={{-128,-52},{-108,-32}}),
          iconTransformation(extent={{-60,-60},{-40,-40}})));

    Phases.Dielectric dielectric if graphite.'incle-' and ionomer.'inclH+'
      "Gap of the double layer"
      annotation (Placement(transformation(extent={{-20,-26},{0,-6}})));

  protected
    Conditions.ByConnector.Amagat.Volume2 volume(final V=V, final setVolume=gas.n_spec
           > 0 or liquid.n_spec > 0) if n_spec > 0
      "Model to establish a fixed total volume"
      annotation (Placement(transformation(extent={{76,-76},{96,-56}})));
    Connectors.InertNode common
      "Connector for translational and thermal exchange among all species"
      annotation (Placement(transformation(extent={{76,32},{96,52}}),
          iconTransformation(extent={{100,18},{120,38}})));
    Connectors.InertNode gasLiq
      "Connector for translational and thermal exchange between gas and liquid"
      annotation (Placement(transformation(extent={{76,56},{96,76}}),
          iconTransformation(extent={{100,18},{120,38}})));
    Connectors.InertNode graphiteLiq
      "Connector for translational and thermal exchange between graphite and liquid"
      annotation (Placement(transformation(extent={{76,44},{96,64}}),
          iconTransformation(extent={{100,18},{120,38}})));

  equation
    // Boundaries and mixing
    // ---------------------
    // Gas
    connect(gas.amagat, volume.amagat) annotation (Line(
        points={{-72,-20},{-72,-66},{86,-66}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(gas.yNegative, yNegative.gas) annotation (Line(
        points={{-80,-20.4},{-80,-54},{-106,-54}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.zPositive, zPositive.gas) annotation (Line(
        points={{-88,-20},{-88,-42},{-118,-42}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.xNegative, xNegative.gas) annotation (Line(
        points={{-88,-12},{-96,-12},{-96,-30},{-130,-30}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.yPositive, yPositive.gas) annotation (Line(
        points={{-80,-2},{-80,30},{86,30}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.zNegative, zNegative.gas) annotation (Line(
        points={{-75,-7},{-74,-7},{-72,-4},{-72,18},{98,18}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.xPositive, xPositive.gas) annotation (Line(
        points={{-72,-12},{-64,-12},{-64,6},{110,6}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    // Graphite
    connect(graphite.amagat, volume.amagat) annotation (Line(
        points={{-32,-20},{-32,-66},{86,-66}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(graphite.yNegative, yNegative.graphite) annotation (Line(
        points={{-40,-20.4},{-40,-54},{-106,-54}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.yPositive, yPositive.graphite) annotation (Line(
        points={{-40,-2},{-40,30},{86,30}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.zNegative, zNegative.graphite) annotation (Line(
        points={{-35,-7},{-34,-7},{-32,-4},{-32,18},{98,18}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.zPositive, zPositive.graphite) annotation (Line(
        points={{-48,-20},{-48,-42},{-118,-42}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.xNegative, xNegative.graphite) annotation (Line(
        points={{-48,-12},{-56,-12},{-56,-30},{-130,-30}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.xPositive, xPositive.graphite) annotation (Line(
        points={{-32,-12},{-24,-12},{-24,6},{110,6}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    // Ionomer
    connect(ionomer.amagat, volume.amagat) annotation (Line(
        points={{28,-20},{28,-66},{86,-66}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(ionomer.yNegative, yNegative.ionomer) annotation (Line(
        points={{20,-20.4},{20,-54},{-106,-54}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(ionomer.yPositive, yPositive.ionomer) annotation (Line(
        points={{20,-2},{20,30},{86,30}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(ionomer.zNegative, zNegative.ionomer) annotation (Line(
        points={{25,-7},{28,-4},{28,18},{98,18}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(ionomer.zPositive, zPositive.ionomer) annotation (Line(
        points={{12,-20},{12,-42},{-118,-42}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(ionomer.xNegative, xNegative.ionomer) annotation (Line(
        points={{12,-12},{4,-12},{4,-30},{-130,-30}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(ionomer.xPositive, xPositive.ionomer) annotation (Line(
        points={{28,-12},{36,-12},{36,6},{110,6}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    // Liquid
    connect(liquid.amagat, volume.amagat) annotation (Line(
        points={{68,-20},{68,-66},{86,-66}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(liquid.yNegative, yNegative.liquid) annotation (Line(
        points={{60,-20.4},{60,-54},{-106,-54}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(liquid.yPositive, yPositive.liquid) annotation (Line(
        points={{60,-2},{60,30},{86,30}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(liquid.zNegative, zNegative.liquid) annotation (Line(
        points={{65,-7},{68,-4},{68,18},{98,18}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(liquid.zPositive, zPositive.liquid) annotation (Line(
        points={{52,-20},{52,-42},{-118,-42}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(liquid.xNegative, xNegative.liquid) annotation (Line(
        points={{52,-12},{44,-12},{44,-30},{-130,-30}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(liquid.xPositive, xPositive.liquid) annotation (Line(
        points={{68,-12},{76,-12},{76,6},{110,6}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    // Dielectric
    connect(dielectric.amagat, volume.amagat) annotation (Line(
        points={{-10,-16},{-10,-66},{86,-66}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(graphite.electrostatic, dielectric.negative) annotation (Line(
        points={{-30,-16},{-16,-16}},
        color={255,195,38},
        smooth=Smooth.None));
    connect(dielectric.positive, ionomer.electrostatic) annotation (Line(
        points={{-4,-16},{10,-16}},
        color={255,195,38},
        smooth=Smooth.None));

    // Inert exchange
    // --------------
    // Common
    connect(gas.inter[1], common.exchange) annotation (Line(
        points={{-85,-6.5},{-85,42},{86,42}},
        color={38,196,52},
        smooth=Smooth.None));
    connect(graphite.inter[1], common.exchange) annotation (Line(
        points={{-45,-6.5},{-45,42},{86,42}},
        color={38,196,52},
        smooth=Smooth.None));
    connect(ionomer.inter[1], common.exchange) annotation (Line(
        points={{15,-7},{15,42},{86,42}},
        color={38,196,52},
        smooth=Smooth.None));
    connect(liquid.inter[1], common.exchange) annotation (Line(
        points={{55,-6.33333},{55,42},{86,42}},
        color={38,196,52},
        smooth=Smooth.None));
    // Gas-liquid
    connect(gas.inter[2], gasLiq.exchange) annotation (Line(
        points={{-85,-7.5},{-85,66},{86,66}},
        color={38,196,52},
        smooth=Smooth.None));
    connect(liquid.inter[2], gasLiq.exchange) annotation (Line(
        points={{55,-7},{55,66},{86,66}},
        color={38,196,52},
        smooth=Smooth.None));
    // Graphite-liquid
    connect(graphite.inter[2], graphiteLiq.exchange) annotation (Line(
        points={{-45,-7.5},{-45,54},{86,54}},
        color={38,196,52},
        smooth=Smooth.None));
    connect(liquid.inter[3], graphiteLiq.exchange) annotation (Line(
        points={{55,-7.66667},{55,54},{86,54}},
        color={38,196,52},
        smooth=Smooth.None));

    // Reactions and phase change (not shown in diagram)
    // -------------------------------------------------
    connect(gas.connH2, HOR.connH2);
    connect(gas.connH2O, ORR.connH2O);
    connect(gas.connO2, ORR.connO2);
    connect(gas.connH2O, ionomer.connH2O);
    connect(gas.connH2O, liquid.connH2O);
    connect(graphite.'conne-', HOR.'conne-');
    connect(graphite.'conne-', ORR.'conne-');
    connect(ionomer.'connH+', HOR.'connH+');
    connect(ionomer.'connH+', ORR.'connH+');
    annotation (Documentation(info="<html>
   <p>Please see the documentation of the
   <a href=\"modelica://FCSys.Subregions.BaseClasses.EmptySubregion\">EmptySubregion</a> model.</p></html>"),
        Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-140,-80},{
              120,80}}), graphics={Text(
            extent={{70,-44},{110,-50}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            textString="(connections not shown 
on diagram)")}));
  end Subregion;

  model SubregionIonomerOnly "Subregion with only the ionomer phase"
    import Modelica.Constants.eps;
    extends Partial(final n_spec=ionomer.n_spec);

    FCSys.Phases.Ionomer ionomer(n_inter=0, final n_trans=n_trans) "Ionomer"
      annotation (Dialog(group="Phases (click to edit)"), Placement(
          transformation(extent={{-20,0},{0,20}})));

    Connectors.FaceBus xNegative if inclTransX "Negative face along the x axis"
      annotation (Placement(transformation(extent={{-60,0},{-40,20}}),
          iconTransformation(extent={{-110,-10},{-90,10}})));
    Connectors.FaceBus yNegative if inclTransY "Negative face along the y axis"
      annotation (Placement(transformation(extent={{-20,-40},{0,-20}}),
          iconTransformation(extent={{-10,-110},{10,-90}})));
    Connectors.FaceBus zNegative if inclTransZ "Negative face along the z axis"
      annotation (Placement(transformation(extent={{0,20},{20,40}}),
          iconTransformation(extent={{40,40},{60,60}})));
    Connectors.FaceBus xPositive if inclTransX "Positive face along the x axis"
      annotation (Placement(transformation(extent={{20,0},{40,20}}),
          iconTransformation(extent={{90,-10},{110,10}})));
    Connectors.FaceBus yPositive if inclTransY "Positive face along the y axis"
      annotation (Placement(transformation(extent={{-20,40},{0,60}}),
          iconTransformation(extent={{-10,90},{10,110}})));
    Connectors.FaceBus zPositive if inclTransZ "Positive face along the z axis"
      annotation (Placement(transformation(extent={{-40,-20},{-20,0}}),
          iconTransformation(extent={{-60,-60},{-40,-40}})));

  protected
    Conditions.ByConnector.Amagat.Volume2 volume(V=V,final setVolume=false) if
      n_spec > 0 "Model to establish a fixed total volume"
      annotation (Placement(transformation(extent={{10,-30},{30,-10}})));
  equation
    // Boundaries and mixing
    // ---------------------
    // Ionomer
    connect(ionomer.amagat, volume.amagat) annotation (Line(
        points={{-2,2},{20,-20}},
        color={47,107,251},
        smooth=Smooth.None));

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

    annotation (Documentation(info="<html>
   <p>Please see the documentation of the
   <a href=\"modelica://FCSys.Subregions.BaseClasses.EmptySubregion\">EmptySubregion</a> model.</p></html>"),
        Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-60,-40},{
              40,60}}), graphics));
  end SubregionIonomerOnly;

  model SubregionNoIonomer "Subregion with all phases except ionomer"
    import Modelica.Constants.eps;
    extends Partial(final n_spec=gas.n_spec + graphite.n_spec + liquid.n_spec);

    parameter Q.NumberAbsolute k_gas_liq=eps
      "Additional coupling factor between gas and liquid" annotation (Dialog(
          group="Geometry",__Dymola_label=
            "<html><i>k</i><sub>gas liq</sub></html>"));
    parameter Q.NumberAbsolute k_graphite_liq=eps
      "Additional coupling factor between graphite and liquid" annotation (
        Dialog(group="Geometry", __Dymola_label=
            "<html><i>k</i><sub>graphite liq</sub></html>"));

    FCSys.Phases.Gas gas(
      n_inter=2,
      final n_trans=n_trans,
      k_inter={k_common,k_gas_liq}) "Gas" annotation (Dialog(group=
            "Phases (click to edit)"), Placement(transformation(extent={{-50,-18},
              {-30,2}})));

    FCSys.Phases.Graphite graphite(
      n_inter=2,
      final n_trans=n_trans,
      k_inter={k_common,k_graphite_liq}) "Graphite" annotation (Dialog(group=
            "Phases (click to edit)"), Placement(transformation(extent={{-10,-18},
              {10,2}})));

    FCSys.Phases.Liquid liquid(
      n_inter=3,
      final n_trans=n_trans,
      k_inter={k_common,k_gas_liq,k_graphite_liq}) "Liquid" annotation (Dialog(
          group="Phases (click to edit)"), Placement(transformation(extent={{30,
              -18},{50,2}})));

    Connectors.FaceBus xNegative if inclTransX "Negative face along the x axis"
      annotation (Placement(transformation(extent={{-100,-36},{-80,-16}}),
          iconTransformation(extent={{-110,-10},{-90,10}})));
    Connectors.FaceBus yNegative if inclTransY "Negative face along the y axis"
      annotation (Placement(transformation(extent={{-76,-60},{-56,-40}}),
          iconTransformation(extent={{-10,-110},{10,-90}})));
    Connectors.FaceBus zNegative if inclTransZ "Negative face along the z axis"
      annotation (Placement(transformation(extent={{68,12},{88,32}}),
          iconTransformation(extent={{40,40},{60,60}})));
    Connectors.FaceBus xPositive if inclTransX "Positive face along the x axis"
      annotation (Placement(transformation(extent={{80,0},{100,20}}),
          iconTransformation(extent={{90,-10},{110,10}})));
    Connectors.FaceBus yPositive if inclTransY "Positive face along the y axis"
      annotation (Placement(transformation(extent={{56,24},{76,44}}),
          iconTransformation(extent={{-10,90},{10,110}})));
    Connectors.FaceBus zPositive if inclTransZ "Positive face along the z axis"
      annotation (Placement(transformation(extent={{-88,-48},{-68,-28}}),
          iconTransformation(extent={{-60,-60},{-40,-40}})));

  protected
    Conditions.ByConnector.Amagat.Volume2 volume(V=V, final setVolume=gas.n_spec
           > 0 or liquid.n_spec > 0) if n_spec > 0
      "Model to establish a fixed total volume"
      annotation (Placement(transformation(extent={{56,-72},{76,-52}})));
    Connectors.InertNode common
      "Connector for translational and thermal exchange among all species"
      annotation (Placement(transformation(extent={{56,36},{76,56}}),
          iconTransformation(extent={{100,18},{120,38}})));
    Connectors.InertNode gasLiq
      "Connector for translational and thermal exchange between gas and liquid"
      annotation (Placement(transformation(extent={{56,60},{76,80}}),
          iconTransformation(extent={{100,18},{120,38}})));
    Connectors.InertNode graphiteLiq
      "Connector for translational and thermal exchange between graphite and liquid"
      annotation (Placement(transformation(extent={{56,48},{76,68}}),
          iconTransformation(extent={{100,18},{120,38}})));

  equation
    // Boundaries and mixing
    // ---------------------
    // Gas
    connect(gas.amagat, volume.amagat) annotation (Line(
        points={{-32,-16},{-32,-62},{66,-62}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(gas.yNegative, yNegative.gas) annotation (Line(
        points={{-40,-16.4},{-40,-50},{-66,-50}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.zPositive, zPositive.gas) annotation (Line(
        points={{-48,-16},{-48,-38},{-78,-38}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.xNegative, xNegative.gas) annotation (Line(
        points={{-48,-8},{-56,-8},{-56,-26},{-90,-26}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.yPositive, yPositive.gas) annotation (Line(
        points={{-40,2},{-40,34},{66,34}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(gas.zNegative, zNegative.gas) annotation (Line(
        points={{-35,-3},{-32,0},{-32,22},{78,22}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    connect(gas.xPositive, xPositive.gas) annotation (Line(
        points={{-32,-8},{-24,-8},{-24,10},{90,10}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    // Graphite
    connect(graphite.amagat, volume.amagat) annotation (Line(
        points={{8,-16},{8,-62},{66,-62}},
        color={47,107,251},
        smooth=Smooth.None));

    connect(graphite.yNegative, yNegative.graphite) annotation (Line(
        points={{0,-16.4},{0,-50},{-66,-50}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.yPositive, yPositive.graphite) annotation (Line(
        points={{0,2},{0,34},{66,34}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.zNegative, zNegative.graphite) annotation (Line(
        points={{5,-3},{8,0},{8,22},{78,22}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.zPositive, zPositive.graphite) annotation (Line(
        points={{-8,-16},{-8,-38},{-78,-38}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.xNegative, xNegative.graphite) annotation (Line(
        points={{-8,-8},{-16,-8},{-16,-26},{-90,-26}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(graphite.xPositive, xPositive.graphite) annotation (Line(
        points={{8,-8},{16,-8},{16,10},{90,10}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    // Liquid
    connect(liquid.amagat, volume.amagat) annotation (Line(
        points={{48,-16},{48,-62},{66,-62}},
        color={47,107,251},
        smooth=Smooth.None));

    connect(liquid.yNegative, yNegative.liquid) annotation (Line(
        points={{40,-16.4},{40,-50},{-66,-50}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(liquid.yPositive, yPositive.liquid) annotation (Line(
        points={{40,2},{40,34},{66,34}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(liquid.zNegative, zNegative.liquid) annotation (Line(
        points={{45,-3},{48,0},{48,22},{78,22}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));
    connect(liquid.zPositive, zPositive.liquid) annotation (Line(
        points={{32,-16},{32,-38},{-78,-38}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    connect(liquid.xNegative, xNegative.liquid) annotation (Line(
        points={{32,-8},{24,-8},{24,-26},{-90,-26}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    connect(liquid.xPositive, xPositive.liquid) annotation (Line(
        points={{48,-8},{56,-8},{56,10},{90,10}},
        color={127,127,127},
        pattern=LinePattern.None,
        thickness=0.5,
        smooth=Smooth.None));

    // Inert exchange
    // --------------
    // Common
    connect(gas.inter[1], common.exchange) annotation (Line(
        points={{-45,-2.5},{-45,46},{66,46}},
        color={38,196,52},
        smooth=Smooth.None));
    connect(graphite.inter[1], common.exchange) annotation (Line(
        points={{-5,-2.5},{-5,46},{66,46}},
        color={38,196,52},
        smooth=Smooth.None));
    connect(liquid.inter[1], common.exchange) annotation (Line(
        points={{35,-2.33333},{35,46},{66,46}},
        color={38,196,52},
        smooth=Smooth.None));
    // Gas-liquid
    connect(gas.inter[2], gasLiq.exchange) annotation (Line(
        points={{-45,-3.5},{-45,70},{66,70}},
        color={38,196,52},
        smooth=Smooth.None));
    connect(liquid.inter[2], gasLiq.exchange) annotation (Line(
        points={{35,-3},{35,70},{66,70}},
        color={38,196,52},
        smooth=Smooth.None));
    // Graphite-liquid
    connect(graphite.inter[2], graphiteLiq.exchange) annotation (Line(
        points={{-5,-3.5},{-5,58},{66,58}},
        color={38,196,52},
        smooth=Smooth.None));
    connect(liquid.inter[3], graphiteLiq.exchange) annotation (Line(
        points={{35,-3.66667},{35,58},{66,58}},
        color={38,196,52},
        smooth=Smooth.None));

    // Phase change (not shown in diagram)
    // -----------------------------------
    connect(gas.connH2O, liquid.connH2O);

    annotation (Documentation(info="<html>
   <p>Please see the documentation of the
   <a href=\"modelica://FCSys.Subregions.BaseClasses.EmptySubregion\">EmptySubregion</a> model.</p></html>"),
        Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-100,-80},{
              100,80}}), graphics));
  end SubregionNoIonomer;

  partial model Partial
    "Base model for multi-dimensional, multi-species storage, transport, and exchange"
    import FCSys.Utilities.Coordinates.cartWrap;
    import Modelica.Math.BooleanVectors.countTrue;
    import Modelica.Math.BooleanVectors.enumerate;
    import Modelica.Math.BooleanVectors.index;
    // extends FCSys.Icons.Names.Top3;

    // Geometric parameters
    inner parameter Q.Length L[Axis](each min=Modelica.Constants.small) = {U.cm,
      U.cm,U.cm} "Lengths" annotation (Dialog(group="Geometry", __Dymola_label=
            "<html><b><i>L</i></b></html>"));
    final inner parameter Q.Area A[Axis]={L[cartWrap(axis + 1)]*L[cartWrap(axis
         + 2)] for axis in Axis} "Cross-sectional areas";
    final inner parameter Q.Volume V=product(L) "Volume";
    parameter Q.NumberAbsolute k_common=1
      "Coupling factor for exchange among all phases" annotation (Dialog(group=
            "Geometry", __Dymola_label="<html><i>k</i><sub>common</sub></html>"));

    // Assumptions
    // -----------
    // Included faces
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

  protected
    parameter Integer n_spec(start=0) "Number of species"
      annotation (HideResult=true);
    final inner parameter Boolean inclTrans[Axis]={inclTransX,inclTransY,
        inclTransZ} "true, if each pairs of faces is included";
    final inner parameter Boolean inclRot[Axis]={inclTransY and inclTransZ,
        inclTransZ and inclTransX,inclTransX and inclTransY}
      "true, if each axis of rotation has all its tangential faces included";
    final inner parameter Integer n_trans=countTrue(inclTrans)
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

    <p>At least one transport axis must be included.
  All of the faces are included by default.</p>

  <p>This model should be extended to include the appropriate phases and reactions.</p>
  </html>"),
      Icon(graphics={
          Line(
            points={{-100,0},{-40,0}},
            color={127,127,127},
            thickness=0.5,
            visible=inclTransX,
            smooth=Smooth.None),
          Line(
            points={{0,-40},{0,-100}},
            color={127,127,127},
            thickness=0.5,
            visible=inclTransY,
            smooth=Smooth.None),
          Line(
            points={{40,40},{50,50}},
            color={127,127,127},
            thickness=0.5,
            visible=inclTransZ,
            smooth=Smooth.None),
          Polygon(
            points={{-40,16},{-16,40},{40,40},{40,-16},{16,-40},{-40,-40},{-40,
                16}},
            lineColor={127,127,127},
            smooth=Smooth.None,
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Line(
            points={{-40,-40},{-16,-16}},
            color={127,127,127},
            smooth=Smooth.None,
            pattern=LinePattern.Dash),
          Line(
            points={{-16,40},{-16,-16},{40,-16}},
            color={127,127,127},
            smooth=Smooth.None,
            pattern=LinePattern.Dash),
          Line(
            points={{-40,0},{28,0}},
            color={210,210,210},
            visible=inclTransX,
            smooth=Smooth.None,
            thickness=0.5),
          Line(
            points={{0,28},{0,-40}},
            color={210,210,210},
            visible=inclTransY,
            smooth=Smooth.None,
            thickness=0.5),
          Line(
            points={{28,0},{100,0}},
            color={127,127,127},
            thickness=0.5,
            visible=inclTransX,
            smooth=Smooth.None),
          Line(
            points={{0,100},{0,28}},
            color={127,127,127},
            thickness=0.5,
            visible=inclTransY,
            smooth=Smooth.None),
          Line(
            points={{-12,-12},{40,40}},
            color={210,210,210},
            visible=inclTransZ,
            smooth=Smooth.None,
            thickness=0.5),
          Line(
            points={{-40,16},{16,16},{16,-40}},
            color={127,127,127},
            smooth=Smooth.None),
          Line(
            points={{-50,-50},{-12,-12}},
            color={127,127,127},
            thickness=0.5,
            visible=inclTransZ,
            smooth=Smooth.None),
          Polygon(
            points={{-40,16},{-16,40},{40,40},{40,-16},{16,-40},{-40,-40},{-40,
                16}},
            lineColor={127,127,127},
            smooth=Smooth.None),
          Line(
            points={{40,40},{16,16}},
            color={127,127,127},
            smooth=Smooth.None),
          Text(
            extent={{-100,56},{100,96}},
            textString="%name",
            lineColor={0,0,0})}),
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
              100,100}}), graphics));

  end Partial;
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
end Subregions;

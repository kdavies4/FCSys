within FCSys;
package Regions "3D arrays of discrete, interconnected subregions"
  extends Modelica.Icons.Package;

  package Examples "Examples"
    model AnFP "Test the anode flow plate"

      extends Modelica.Icons.Example;

      parameter Q.NumberAbsolute psi_H2O=environment.psi_H2O
        "Mole fraction of H2O at the inlet";
      parameter Q.NumberAbsolute psi_H2=environment.psi_dry
        "Mole fraction of H2 at the inlet";

      output Q.Potential w=anFP.subregions[1, 1, 1].graphite.'e-'.Deltag[1] if
        environment.analysis "Electrical potential";
      output Q.ResistanceElectrical R=w/zI if environment.analysis
        "Measured electrical resistance";
      output Q.ResistanceElectrical R_ex=anFP.L[Axis.x]/(anFP.subregions[1, 1,
          1].graphite.'e-'.sigma*anFP.A[Axis.x]*anFP.subregions[1, 1, 1].graphite.epsilon
          ^1.5) if environment.analysis "Expected electrical resistance";

      AnFPs.AnFP anFP
        annotation (Placement(transformation(extent={{-70,-10},{-50,10}})));

      // Conditions
      Conditions.ByConnector.BoundaryBus.Single.Sink anBC[anFP.n_y, anFP.n_z](
          each graphite(
          'inclC+'=true,
          'incle-'=true,
          redeclare Conditions.ByConnector.ThermalDiffusive.Single.Temperature
            'C+'(set(y=environment.T)),
          'e-'(materialSet(y=0)))) annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=90,
            origin={-84,0})));

      Conditions.ByConnector.BoundaryBus.Single.Source caBC[anFP.n_y, anFP.n_z]
        (
        each gas(
          inclH2=true,
          inclH2O=true,
          H2(
            materialSet(y=zI/2),
            redeclare function thermalSpec =
                Conditions.ByConnector.Boundary.Single.Thermal.heatRate,
            thermalSet(y=0)),
          H2O(
            redeclare function materialSpec =
                Conditions.ByConnector.Boundary.Single.Material.pressure,
            materialSet(y=environment.p_H2O),
            redeclare function thermalSpec =
                Conditions.ByConnector.Boundary.Single.Thermal.heatRate,
            thermalSet(y=0))),
        each graphite('incle-'=true, 'e-'(materialSet(y=-zI), thermalSet(y=
                  environment.T))),
        each liquid(inclH2O=true, H2O(redeclare function thermalSpec =
                Conditions.ByConnector.Boundary.Single.Thermal.heatRate,
              thermalSet(y=0)))) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-36,0})));

      Conditions.ByConnector.BoundaryBus.Single.Source anSource[anFP.n_x, anFP.n_z]
        (each gas(
          inclH2=true,
          inclH2O=true,
          H2(materialSet(y=-testConditions.Ndot_H2), thermalSet(y=environment.T)),

          H2O(materialSet(y=-testConditions.Ndot_H2O_an), thermalSet(y=
                  environment.T)))) annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=180,
            origin={-60,-24})));

      Conditions.ByConnector.BoundaryBus.Single.Sink anSink[anFP.n_x, anFP.n_z]
        (gas(
          each inclH2=true,
          each inclH2O=true,
          H2O(materialSet(y=fill(
                      environment.p,
                      anFP.n_x,
                      anFP.n_z) - anSink.gas.H2.p)),
          H2(materialSet(y=anSink.gas.H2O.boundary.Ndot .* anFP.subregions[:,
                  anFP.n_y, :].gas.H2O.v ./ anFP.subregions[:, anFP.n_y, :].gas.H2.v),
              redeclare each function materialSpec =
                Conditions.ByConnector.Boundary.Single.Material.current)),
          each liquid(inclH2O=true, H2O(materialSet(y=environment.p))))
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={-60,24})));

      inner Conditions.Environment environment(
        analysis=true,
        T=333.15*U.K,
        p=U.from_kPag(48.3),
        RH=0.8) "Environmental conditions"
        annotation (Placement(transformation(extent={{-30,30},{-10,50}})));

      Modelica.Blocks.Sources.Ramp currentSet(
        height=100*U.A,
        duration=300,
        offset=U.mA)
        annotation (Placement(transformation(extent={{-20,-60},{0,-40}})));

    protected
      Connectors.RealOutputInternal zI(unit="N/T") "Electrical current"
        annotation (Placement(transformation(extent={{0,-60},{20,-40}})));
    public
      Assemblies.Cells.Examples.TestConditions testConditions(I_ca=2*zI, I_an=
            1.5*zI) "Test conditions" annotation (Dialog,Placement(
            transformation(extent={{10,30},{30,50}})));
    equation
      connect(anBC.boundary, anFP.xNegative) annotation (Line(
          points={{-80,0},{-70,0}},
          color={127,127,127},
          pattern=LinePattern.None,
          thickness=0.5,
          smooth=Smooth.None));

      connect(caBC.boundary, anFP.xPositive) annotation (Line(
          points={{-40,0},{-50,0}},
          color={240,0,0},
          thickness=0.5,
          smooth=Smooth.None));

      connect(anSource.boundary, anFP.yNegative) annotation (Line(
          points={{-60,-20},{-60,-10}},
          color={240,0,0},
          thickness=0.5,
          smooth=Smooth.None));
      connect(anSink.boundary, anFP.yPositive) annotation (Line(
          points={{-60,20},{-60,10}},
          color={240,0,0},
          thickness=0.5,
          smooth=Smooth.None));
      connect(currentSet.y, zI) annotation (Line(
          points={{1,-50},{10,-50}},
          color={0,0,127},
          smooth=Smooth.None));
      annotation (
        experiment(StopTime=350),
        Commands(file="Resources/Scripts/Dymola/Regions.Examples.AnFP.mos"
            "Regions.Examples.AnFP.mos", file=
              "Resources/Scripts/Dymola/Regions.Examples.anFP-states.mos"
            "Regions.Examples.anFP-states.mos"),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics));
    end AnFP;
    extends Modelica.Icons.ExamplesPackage;

    model AnGDL "Test the anode gas diffusion layer"

      extends Modelica.Icons.Example;

      output Q.Potential w=anGDL.subregions[1, 1, 1].graphite.'e-'.Deltag[1]
        if environment.analysis "Electrical potential";
      output Q.Current zI=-sum(anGDL.subregions[1, :, :].graphite.'e-'.I[1])
        if environment.analysis "Electrical current";
      output Q.ResistanceElectrical R=w/zI if environment.analysis
        "Measured electrical resistance";
      output Q.ResistanceElectrical R_ex=anGDL.L[Axis.x]/(anGDL.subregions[1, 1,
          1].graphite.'e-'.sigma*anGDL.A[Axis.x]*anGDL.subregions[1, 1, 1].graphite.epsilon
          ^1.5) if environment.analysis "Expected electrical resistance";

      AnGDLs.AnGDL anGDL
        annotation (Placement(transformation(extent={{-50,-10},{-30,10}})));

      // Conditions
      Conditions.ByConnector.BoundaryBus.Single.Sink anBC[anGDL.n_y, anGDL.n_z]
        (
        each gas(
          inclH2=true,
          inclH2O=true,
          H2(materialSet(y=environment.p_dry)),
          H2O(materialSet(y=environment.p_H2O))),
        each graphite(
          'inclC+'=true,
          'incle-'=true,
          'e-'(materialSet(y=0))),
        each liquid(inclH2O=true, H2O(materialSet(y=environment.p))))
        annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=90,
            origin={-64,0})));

      Conditions.ByConnector.BoundaryBus.Single.Source caBC[anGDL.n_y, anGDL.n_z]
        (
        each graphite(
          'inclC+'=true,
          'incle-'=true,
          'C+'(set(y=environment.T)),
          'e-'(materialSet(y=-currentSet.y), thermalSet(y=environment.T))),
        each liquid(inclH2O=true, H2O(
            thermalSet(y=environment.T),
            redeclare function afterSpec =
                FCSys.Conditions.ByConnector.Boundary.Single.Translational.force,

            redeclare function beforeSpec =
                FCSys.Conditions.ByConnector.Boundary.Single.Translational.force)),

        each gas(
          inclH2=true,
          inclH2O=true,
          H2O(
            thermalSet(y=environment.T),
            redeclare function afterSpec =
                FCSys.Conditions.ByConnector.Boundary.Single.Translational.force,

            redeclare function beforeSpec =
                FCSys.Conditions.ByConnector.Boundary.Single.Translational.force),

          H2(
            thermalSet(y=environment.T),
            redeclare function afterSpec =
                FCSys.Conditions.ByConnector.Boundary.Single.Translational.force,

            redeclare function beforeSpec =
                FCSys.Conditions.ByConnector.Boundary.Single.Translational.force)))
        annotation (Placement(transformation(
            extent={{10,10},{-10,-10}},
            rotation=90,
            origin={-16,0})));

      inner Conditions.Environment environment(
        analysis=true,
        T=333.15*U.K,
        p=U.from_kPag(48.3),
        RH=0.8) "Environmental conditions"
        annotation (Placement(transformation(extent={{-10,30},{10,50}})));

      Modelica.Blocks.Sources.Ramp currentSet(
        duration=20,
        height=100*U.A,
        offset=U.mA)
        annotation (Placement(transformation(extent={{-10,-60},{10,-40}})));

    equation
      connect(anGDL.xPositive, caBC.boundary) annotation (Line(
          points={{-30,0},{-20,0}},
          color={240,0,0},
          thickness=0.5,
          smooth=Smooth.None));
      connect(anBC.boundary, anGDL.xNegative) annotation (Line(
          points={{-60,0},{-50,0}},
          color={240,0,0},
          thickness=0.5,
          smooth=Smooth.None));
      annotation (
        experiment(StopTime=30),
        Commands(file="Resources/Scripts/Dymola/Regions.Examples.AnGDL.mos"
            "Regions.Examples.AnGDL.mos"),
        __Dymola_experimentSetupOutput);
    end AnGDL;

    model AnCL "Test the anode catalyst layer"

      extends Modelica.Icons.Example;

      output Q.Current zI=sum(anCL.subregions[1, :, :].graphite.'e-Transfer'.I)
        if environment.analysis "Reaction rate";
      output Q.Potential w=anCL.subregions[1, 1, 1].graphite.'e-Transfer'.Deltag
        "Overpotential";
      output Q.Number J_Apercm2=zI*U.cm^2/(anCL.A[Axis.x]*U.A)
        "Electrical current density, in A/cm2";

      AnCLs.AnCL anCL
        annotation (Placement(transformation(extent={{-30,-10},{-10,10}})));

      // Conditions
      Conditions.ByConnector.BoundaryBus.Single.Sink anBC[anCL.n_y, anCL.n_z](
        each gas(
          inclH2=true,
          inclH2O=true,
          H2(materialSet(y=environment.p_dry)),
          H2O(materialSet(y=environment.p_H2O))),
        each graphite(
          'inclC+'=true,
          'incle-'=true,
          'e-'(
            redeclare function materialSpec =
                Conditions.ByConnector.Boundary.Single.Material.current,
            materialSet(y=currentSet.y),
            redeclare function thermalSpec =
                Conditions.ByConnector.Boundary.Single.Thermal.temperature,
            thermalSet(y=environment.T))),
        each liquid(inclH2O=true, H2O(materialSet(y=environment.p - U.kPa))))
        annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=270,
            origin={-44,0})));

      Conditions.ByConnector.BoundaryBus.Single.Sink caBC[anCL.n_y, anCL.n_z](
          each ionomer(
          'inclH+'=true,
          inclH2O=false,
          'inclSO3-'=true,
          'H+'(materialSet(y=0)))) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={4,0})));

      inner Conditions.Environment environment(
        analysis=true,
        T=333.15*U.K,
        p=U.from_kPag(48.3),
        RH=0.8) "Environmental conditions"
        annotation (Placement(transformation(extent={{-10,30},{10,50}})));
      Modelica.Blocks.Sources.Ramp currentSet(
        duration=20,
        height=100*U.A,
        offset=U.mA)
        annotation (Placement(transformation(extent={{-10,-60},{10,-40}})));
    equation
      connect(anBC.boundary, anCL.xNegative) annotation (Line(
          points={{-40,0},{-30,0}},
          color={240,0,0},
          thickness=0.5,
          smooth=Smooth.None));

      connect(anCL.xPositive, caBC.boundary) annotation (Line(
          points={{-10,0},{0,0}},
          color={240,0,0},
          thickness=0.5,
          smooth=Smooth.None));

      annotation (
        Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                {100,100}}), graphics),
        experiment(StopTime=30),
        Commands(file="Resources/Scripts/Dymola/Regions.Examples.AnCL.mos"
            "Regions.Examples.AnCL.mos"),
        __Dymola_experimentSetupOutput);
    end AnCL;

    model PEM "Test the proton exchange membrane"

      extends Modelica.Icons.Example;

      output Q.Potential w=environment.T*log(PEM.subregions[1, 1, 1].ionomer.
          'H+'.boundaries[1, Side.p].p/PEM.subregions[1, 1, 1].ionomer.'H+'.boundaries[
          1, Side.n].p) if environment.analysis
        "Isothermal electrical potential";
      output Q.Current zI=-sum(PEM.subregions[1, :, :].ionomer.'H+'.I[1]) if
        environment.analysis "Electrical current";
      output Q.ConductanceElectrical G=zI/w if environment.analysis
        "Measured electrical conductance";
      output Q.ConductanceElectrical G_ex=PEM.subregions[1, 1, 1].ionomer.'H+'.sigma
          *PEM.A[Axis.x]/PEM.L[Axis.x] if environment.analysis
        "Expected electrical conductance";

      PEMs.PEM PEM
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

      // Conditions
      Conditions.ByConnector.BoundaryBus.Single.Source anBC[PEM.n_y, PEM.n_z](
          each ionomer(
          'inclSO3-'=true,
          'inclH+'=true,
          inclH2O=true,
          'H+'(materialSet(y=-currentSet.y), thermalSet(y=environment.T)),
          H2O(thermalSet(y=environment.T)),
          'SO3-'(set(y=environment.T)))) annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=90,
            origin={-24,0})));

      Conditions.ByConnector.BoundaryBus.Single.Sink caBC[PEM.n_y, PEM.n_z](
          each ionomer(
          'inclSO3-'=true,
          'inclH+'=true,
          inclH2O=false,
          'H+'(materialSet(y=U.atm)))) annotation (Placement(transformation(
            extent={{10,10},{-10,-10}},
            rotation=90,
            origin={24,0})));

      inner Conditions.Environment environment(
        analysis=true,
        T=333.15*U.K,
        p=U.from_kPag(48.3),
        RH=0.7) "Environmental conditions"
        annotation (Placement(transformation(extent={{-10,30},{10,50}})));

      Modelica.Blocks.Sources.Ramp currentSet(
        duration=20,
        height=100*U.A,
        offset=U.mA)
        annotation (Placement(transformation(extent={{-10,-60},{10,-40}})));
    equation
      connect(anBC.boundary, PEM.xNegative) annotation (Line(
          points={{-20,0},{-10,0}},
          color={240,0,0},
          thickness=0.5,
          smooth=Smooth.None));

      connect(PEM.xPositive, caBC.boundary) annotation (Line(
          points={{10,0},{20,0}},
          color={0,0,240},
          thickness=0.5,
          smooth=Smooth.None));

      annotation (
        Documentation(info="<html><p>The protonic resistance is slightly higher than expected
  due solely to electrical resistance.  Electro-osmotic drag applies additional resistance to the migration of protons.</p>
  </html>"),
        experiment(StopTime=30),
        Commands(file="Resources/Scripts/Dymola/Regions.Examples.PEM.mos"
            "Regions.Examples.PEM.mos"),
        __Dymola_experimentSetupOutput);
    end PEM;

    model CaCL "Test the cathode catalyst layer"

      extends Modelica.Icons.Example;

      output Q.Current zI=-sum(caCL.subregions[1, :, :].graphite.'e-Transfer'.I)
        if environment.analysis "Reaction rate";
      output Q.Potential w=-caCL.subregions[1, 1, 1].graphite.'e-Transfer'.Deltag
        "Overpotential";
      output Q.Number J_Apercm2=zI*U.cm^2/(caCL.A[Axis.x]*U.A)
        "Electrical current density, in A/cm2";

      CaCLs.CaCL caCL
        annotation (Placement(transformation(extent={{10,-10},{30,10}})));

      // Conditions
      Conditions.ByConnector.BoundaryBus.Single.Source anBC[caCL.n_y, caCL.n_z]
        (each ionomer(
          'inclH+'=true,
          inclH2O=false,
          'inclSO3-'=true,
          'H+'(
            redeclare function materialSpec =
                Conditions.ByConnector.Boundary.Single.Material.pressure,
            materialSet(y=0),
            thermalSet(y=environment.T)),
          'SO3-'(set(y=environment.T)))) annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=270,
            origin={-4,0})));
      Conditions.ByConnector.BoundaryBus.Single.Source caBC[caCL.n_y, caCL.n_z]
        (
        each gas(
          inclH2O=true,
          inclN2=true,
          inclO2=true,
          H2O(
            redeclare function materialSpec =
                Conditions.ByConnector.Boundary.Single.Material.pressure,
            materialSet(y=environment.p_H2O),
            thermalSet(y=environment.T)),
          N2(
            redeclare function materialSpec =
                Conditions.ByConnector.Boundary.Single.Material.pressure,
            materialSet(y=environment.p_dry - environment.p_O2),
            thermalSet(y=environment.T)),
          O2(
            redeclare function materialSpec =
                Conditions.ByConnector.Boundary.Single.Material.pressure,
            materialSet(y=environment.p_O2),
            thermalSet(y=environment.T))),
        each graphite(
          'inclC+'=false,
          'incle-'=true,
          'C+'(set(y=environment.T)),
          'e-'(
            redeclare function materialSpec =
                Conditions.ByConnector.Boundary.Single.Material.current,
            materialSet(y=-currentSet.y),
            thermalSet(y=environment.T))),
        each liquid(inclH2O=true, H2O(
            redeclare function materialSpec =
                Conditions.ByConnector.Boundary.Single.Material.pressure,
            materialSet(y=environment.p),
            thermalSet(y=environment.T)))) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={44,0})));

      inner Conditions.Environment environment(
        analysis=true,
        T=333.15*U.K,
        p=U.from_kPag(48.3),
        RH=0.6) "Environmental conditions"
        annotation (Placement(transformation(extent={{-10,30},{10,50}})));
      Modelica.Blocks.Sources.Ramp currentSet(
        duration=20,
        height=100*U.A,
        offset=U.mA)
        annotation (Placement(transformation(extent={{-10,-60},{10,-40}})));

    equation
      connect(anBC.boundary, caCL.xNegative) annotation (Line(
          points={{0,0},{10,0}},
          color={0,0,240},
          thickness=0.5,
          smooth=Smooth.None));

      connect(caCL.xPositive, caBC.boundary) annotation (Line(
          points={{30,0},{40,0}},
          color={0,0,240},
          thickness=0.5,
          smooth=Smooth.None));

      annotation (
        Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},
                {100,100}}), graphics),
        experiment(StopTime=50),
        Commands(file="Resources/Scripts/Dymola/Regions.Examples.CaCL.mos"
            "Regions.Examples.CaCL.mos"),
        __Dymola_experimentSetupOutput);
    end CaCL;

    model CaGDL "Test the cathode gas diffusion layer"

      extends Modelica.Icons.Example;

      output Q.Potential w=caGDL.subregions[1, 1, 1].graphite.'e-'.Deltag[1]
        if environment.analysis "Electrical potential";
      output Q.Current zI=-sum(caGDL.subregions[1, :, :].graphite.'e-'.I[1])
        if environment.analysis "Electrical current";
      output Q.ResistanceElectrical R=w/zI if environment.analysis
        "Measured electrical resistance";
      output Q.ResistanceElectrical R_ex=caGDL.L[Axis.x]/(caGDL.subregions[1, 1,
          1].graphite.'e-'.sigma*caGDL.A[Axis.x]*caGDL.subregions[1, 1, 1].graphite.epsilon
          ^1.5) if environment.analysis "Expected electrical resistance";

      CaGDLs.CaGDL caGDL
        annotation (Placement(transformation(extent={{30,-10},{50,10}})));

      // Conditions
      Conditions.ByConnector.BoundaryBus.Single.Sink anBC[caGDL.n_y, caGDL.n_z]
        (
        each gas(
          inclH2O=true,
          inclN2=true,
          inclO2=true,
          H2O(materialSet(y=environment.p_H2O)),
          N2(materialSet(y=environment.p_dry - environment.p_O2)),
          O2(materialSet(y=environment.p_O2))),
        each graphite(
          'inclC+'=true,
          'incle-'=true,
          redeclare Conditions.ByConnector.ThermalDiffusive.Single.Temperature
            'C+'(set(y=environment.T)),
          'e-'(redeclare function thermalSpec =
                Conditions.ByConnector.Boundary.Single.Thermal.temperature,
              thermalSet(y=environment.T))),
        each liquid(inclH2O=true, H2O(materialSet(y=environment.p))))
        annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=90,
            origin={16,0})));

      Conditions.ByConnector.BoundaryBus.Single.Source caBC[caGDL.n_y, caGDL.n_z]
        (
        each gas(
          inclH2O=true,
          inclN2=true,
          inclO2=true,
          H2O(
            thermalSet(y=environment.T),
            redeclare function afterSpec =
                FCSys.Conditions.ByConnector.Boundary.Single.Translational.force,

            redeclare function beforeSpec =
                FCSys.Conditions.ByConnector.Boundary.Single.Translational.force),

          N2(
            thermalSet(y=environment.T),
            redeclare function afterSpec =
                FCSys.Conditions.ByConnector.Boundary.Single.Translational.force,

            redeclare function beforeSpec =
                FCSys.Conditions.ByConnector.Boundary.Single.Translational.force),

          O2(
            thermalSet(y=environment.T),
            redeclare function afterSpec =
                FCSys.Conditions.ByConnector.Boundary.Single.Translational.force,

            redeclare function beforeSpec =
                FCSys.Conditions.ByConnector.Boundary.Single.Translational.force)),

        each graphite(
          'inclC+'=true,
          'incle-'=true,
          'C+'(set(y=environment.T)),
          'e-'(materialSet(y=-currentSet.y), thermalSet(y=environment.T))),
        each liquid(inclH2O=true, H2O(
            thermalSet(y=environment.T),
            redeclare function afterSpec =
                FCSys.Conditions.ByConnector.Boundary.Single.Translational.force,

            redeclare function beforeSpec =
                FCSys.Conditions.ByConnector.Boundary.Single.Translational.force)))
        annotation (Placement(transformation(
            extent={{10,10},{-10,-10}},
            rotation=90,
            origin={64,0})));

      inner Conditions.Environment environment(
        analysis=true,
        T=333.15*U.K,
        p=U.from_kPag(48.3),
        RH=0.6) "Environmental conditions"
        annotation (Placement(transformation(extent={{-10,30},{10,50}})));

      Modelica.Blocks.Sources.Ramp currentSet(
        duration=20,
        height=100*U.A,
        offset=U.mA)
        annotation (Placement(transformation(extent={{-10,-60},{10,-40}})));

    equation
      connect(anBC.boundary, caGDL.xNegative) annotation (Line(
          points={{20,0},{30,0}},
          color={0,0,240},
          thickness=0.5,
          smooth=Smooth.None));
      connect(caGDL.xPositive, caBC.boundary) annotation (Line(
          points={{50,0},{60,0}},
          color={0,0,240},
          thickness=0.5,
          smooth=Smooth.None));
      annotation (
        experiment(StopTime=30),
        Commands(file="Resources/Scripts/Dymola/Regions.Examples.CaGDL.mos"
            "Regions.Examples.CaGDL.mos"),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics));
    end CaGDL;

    model CLtoCL "Test one catalyst layer to the other"

      extends Modelica.Icons.Example;

      output Q.Potential w=anCL.subregions[1, 1, 1].graphite.'e-'.g_boundaries[
          1, Side.n] - caCL.subregions[1, 1, 1].graphite.'e-'.g_boundaries[1,
          Side.p] if environment.analysis "Electrical potential";
      output Q.Current zI=-sum(anCL.subregions[1, :, :].graphite.'e-'.boundaries[
          1, Side.n].Ndot) if environment.analysis "Electrical current";
      output Q.Number J_Apercm2=zI*U.cm^2/(anCL.A[Axis.x]*U.A) if environment.analysis
        "Electrical current density, in A/cm2";

      parameter Q.Length L_y[:]={8}*U.cm "Lengths in the y direction";
      parameter Q.Length L_z[:]={6.25}*U.cm "Lengths in the z direction";

      AnCLs.AnCL anCL(final L_y=L_y, final L_z=L_z)
        annotation (Placement(transformation(extent={{-30,-10},{-10,10}})));

      PEMs.PEM PEM(final L_y=L_y, final L_z=L_z)
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
      CaCLs.CaCL caCL(final L_y=L_y, final L_z=L_z)
        annotation (Placement(transformation(extent={{10,-10},{30,10}})));

      // Conditions
      Conditions.ByConnector.BoundaryBus.Single.Sink anBC[anCL.n_y, anCL.n_z](
        each gas(
          inclH2=true,
          inclH2O=true,
          H2(
            materialSet(y=environment.p_dry),
            redeclare function thermalSpec =
                Conditions.ByConnector.Boundary.Single.Thermal.temperature,
            thermalSet(y=environment.T)),
          H2O(
            redeclare function materialSpec =
                Conditions.ByConnector.Boundary.Single.Material.pressure,
            materialSet(y=environment.p_H2O),
            redeclare function thermalSpec =
                Conditions.ByConnector.Boundary.Single.Thermal.temperature,
            thermalSet(y=environment.T))),
        each graphite(
          'inclC+'=true,
          'incle-'=true,
          redeclare Conditions.ByConnector.ThermalDiffusive.Single.Temperature
            'C+'(set(y=environment.T)),
          'e-'(redeclare function thermalSpec =
                Conditions.ByConnector.Boundary.Single.Thermal.temperature,
              thermalSet(y=environment.T))),
        each liquid(inclH2O=true, H2O(materialSet(y=environment.p - U.kPa))))
        annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=270,
            origin={-44,0})));

      Conditions.ByConnector.BoundaryBus.Single.Source caBC[caCL.n_y, caCL.n_z]
        (
        each gas(
          inclH2O=true,
          inclO2=true,
          H2O(
            redeclare function materialSpec =
                Conditions.ByConnector.Boundary.Single.Material.pressure,
            materialSet(y=environment.p_H2O),
            redeclare function thermalSpec =
                Conditions.ByConnector.Boundary.Single.Thermal.heatRate,
            thermalSet(y=0)),
          O2(
            redeclare function materialSpec =
                Conditions.ByConnector.Boundary.Single.Material.pressure,
            materialSet(y=environment.p_O2),
            thermalSet(y=environment.T))),
        each graphite(
          'inclC+'=true,
          'incle-'=true,
          'C+'(set(y=environment.T)),
          'e-'(
            redeclare function materialSpec =
                Conditions.ByConnector.Boundary.Single.Material.current,
            materialSet(y=-currentSet.y),
            thermalSet(y=environment.T))),
        each liquid(inclH2O=true, H2O(
            redeclare function materialSpec =
                Conditions.ByConnector.Boundary.Single.Material.pressure,
            materialSet(y=environment.p),
            thermalSet(y=environment.T)))) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={44,0})));

      inner Conditions.Environment environment(
        T=333.15*U.K,
        p=U.from_kPag(48.3),
        analysis=true,
        RH=0.7) "Environmental conditions"
        annotation (Placement(transformation(extent={{-10,30},{10,50}})));

      Modelica.Blocks.Sources.Ramp currentSet(duration=300, height=100*U.A)
        annotation (Placement(transformation(extent={{-10,-60},{10,-40}})));

    equation
      connect(anCL.xPositive, PEM.xNegative) annotation (Line(
          points={{-10,0},{-10,0}},
          color={240,0,0},
          thickness=0.5,
          smooth=Smooth.None));
      connect(PEM.xPositive, caCL.xNegative) annotation (Line(
          points={{10,0},{10,0}},
          color={0,0,240},
          thickness=0.5,
          smooth=Smooth.None));
      connect(anBC.boundary, anCL.xNegative) annotation (Line(
          points={{-40,0},{-30,0}},
          color={240,0,0},
          thickness=0.5,
          smooth=Smooth.None));
      connect(caCL.xPositive, caBC.boundary) annotation (Line(
          points={{30,0},{40,0}},
          color={0,0,240},
          thickness=0.5,
          smooth=Smooth.None));
      annotation (
        Commands(file="Resources/Scripts/Dymola/Regions.Examples.CLtoCL.mos"
            "Regions.Examples.CLtoCL.mos"),
        experiment(StopTime=350, __Dymola_Algorithm="Dassl"),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics));
    end CLtoCL;

    model GDLtoGDL "Test one gas diffusion layer to the other"

      extends Modelica.Icons.Example;

      output Q.Potential w=anCL.subregions[1, 1, 1].graphite.'e-'.g_boundaries[
          1, Side.n] - caCL.subregions[1, 1, 1].graphite.'e-'.g_boundaries[1,
          Side.p] if environment.analysis "Electrical potential";
      output Q.Current zI=-sum(anCL.subregions[1, :, :].graphite.'e-'.boundaries[
          1, Side.n].Ndot) if environment.analysis "Electrical current";
      output Q.Number J_Apercm2=zI*U.cm^2/(anCL.A[Axis.x]*U.A)
        "Electrical current density, in A/cm2";

      parameter Q.Length L_y[:]={8}*U.cm "Lengths in the y direction";
      parameter Q.Length L_z[:]={6.25}*U.cm "Lengths in the z direction";

      AnGDLs.AnGDL anGDL
        annotation (Placement(transformation(extent={{-50,-10},{-30,10}})));
      AnCLs.AnCL anCL(final L_y=L_y, final L_z=L_z)
        annotation (Placement(transformation(extent={{-30,-10},{-10,10}})));
      PEMs.PEM PEM(final L_y=L_y, final L_z=L_z)
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

      CaCLs.CaCL caCL(final L_y=L_y, final L_z=L_z)
        annotation (Placement(transformation(extent={{10,-10},{30,10}})));
      CaGDLs.CaGDL caGDL(
        final L_y=L_y,
        final L_z=L_z,
        subregions(liquid(H2O(each phi(each stateSelect=StateSelect.default,
                  each fixed=false)))))
        annotation (Placement(transformation(extent={{30,-10},{50,10}})));

      // Conditions
      Conditions.ByConnector.BoundaryBus.Single.Sink anBC[anCL.n_y, anCL.n_z](
        each gas(
          inclH2=true,
          inclH2O=true,
          H2(
            materialSet(y=environment.p_dry),
            redeclare function thermalSpec =
                Conditions.ByConnector.Boundary.Single.Thermal.temperature,
            thermalSet(y=environment.T)),
          H2O(
            materialSet(y=environment.p_H2O),
            redeclare function thermalSpec =
                Conditions.ByConnector.Boundary.Single.Thermal.temperature,
            thermalSet(y=environment.T))),
        each graphite(
          'inclC+'=true,
          'incle-'=true,
          redeclare Conditions.ByConnector.ThermalDiffusive.Single.Temperature
            'C+'(set(y=environment.T)),
          'e-'(
            materialSet(y=0),
            redeclare function thermalSpec =
                Conditions.ByConnector.Boundary.Single.Thermal.temperature,
            thermalSet(y=environment.T))),
        each liquid(inclH2O=true, H2O(
            materialSet(y=environment.p),
            redeclare function thermalSpec =
                Conditions.ByConnector.Boundary.Single.Thermal.temperature,
            thermalSet(y=environment.T)))) annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=270,
            origin={-64,0})));

      Conditions.ByConnector.BoundaryBus.Single.Source caBC[caGDL.n_y, caGDL.n_z]
        (
        each gas(
          inclH2O=true,
          inclN2=true,
          inclO2=true,
          H2O(
            redeclare function materialSpec =
                Conditions.ByConnector.Boundary.Single.Material.pressure,
            materialSet(y=environment.p_H2O),
            thermalSet(y=environment.T),
            redeclare function afterSpec =
                FCSys.Conditions.ByConnector.Boundary.Single.Translational.force,

            redeclare function beforeSpec =
                FCSys.Conditions.ByConnector.Boundary.Single.Translational.force),

          N2(
            redeclare function materialSpec =
                Conditions.ByConnector.Boundary.Single.Material.pressure,
            materialSet(y=environment.p*(1 - environment.psi_O2_dry)),
            thermalSet(y=environment.T),
            redeclare function afterSpec =
                FCSys.Conditions.ByConnector.Boundary.Single.Translational.force,

            redeclare function beforeSpec =
                FCSys.Conditions.ByConnector.Boundary.Single.Translational.force),

          O2(
            redeclare function materialSpec =
                Conditions.ByConnector.Boundary.Single.Material.pressure,
            materialSet(y=environment.p_O2),
            thermalSet(y=environment.T),
            redeclare function afterSpec =
                FCSys.Conditions.ByConnector.Boundary.Single.Translational.force,

            redeclare function beforeSpec =
                FCSys.Conditions.ByConnector.Boundary.Single.Translational.force)),

        each graphite(
          'inclC+'=true,
          'incle-'=true,
          'C+'(set(y=environment.T)),
          'e-'(
            redeclare function materialSpec =
                Conditions.ByConnector.Boundary.Single.Material.current,
            materialSet(y=-currentSet.y),
            thermalSet(y=environment.T))),
        each liquid(inclH2O=true,H2O(
            materialSet(y=0),
            thermalSet(y=environment.T),
            redeclare function afterSpec =
                FCSys.Conditions.ByConnector.Boundary.Single.Translational.force,

            redeclare function beforeSpec =
                FCSys.Conditions.ByConnector.Boundary.Single.Translational.force)))
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={64,0})));

      inner Conditions.Environment environment(
        analysis=true,
        T=333.15*U.K,
        p=U.from_kPag(48.3),
        RH=0.7) "Environmental conditions"
        annotation (Placement(transformation(extent={{-10,30},{10,50}})));

      Modelica.Blocks.Sources.Ramp currentSet(duration=300, height=100*U.A)
        annotation (Placement(transformation(extent={{-10,-60},{10,-40}})));

    equation
      connect(anCL.xPositive, PEM.xNegative) annotation (Line(
          points={{-10,0},{-10,0}},
          color={240,0,0},
          thickness=0.5,
          smooth=Smooth.None));
      connect(PEM.xPositive, caCL.xNegative) annotation (Line(
          points={{10,0},{10,0}},
          color={0,0,240},
          thickness=0.5,
          smooth=Smooth.None));
      connect(caGDL.xPositive, caBC.boundary) annotation (Line(
          points={{50,0},{60,0}},
          color={0,0,240},
          thickness=0.5,
          smooth=Smooth.None));
      connect(caCL.xPositive, caGDL.xNegative) annotation (Line(
          points={{30,0},{30,0}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));
      connect(anBC.boundary, anGDL.xNegative) annotation (Line(
          points={{-60,0},{-50,0}},
          color={240,0,0},
          thickness=0.5,
          smooth=Smooth.None));
      connect(anGDL.xPositive, anCL.xNegative) annotation (Line(
          points={{-30,0},{-30,0}},
          color={240,0,0},
          thickness=0.5,
          smooth=Smooth.None));
      annotation (
        Commands(file="Resources/Scripts/Dymola/Regions.Examples.GDLtoGDL.mos"
            "Regions.Examples.GDLtoGDL.mos", file=
              "Resources/Scripts/Dymola/Regions.Examples.GDLtoGDL-states.mos"
            "Regions.Examples.GDLtoGDL-states.mos"),
        experiment(StopTime=350, __Dymola_Algorithm="Dassl"),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics));
    end GDLtoGDL;

    model CaFP "Test the cathode flow plate"

      extends Modelica.Icons.Example;

      parameter Q.NumberAbsolute psi_H2O=environment.psi_H2O
        "Mole fraction of H2O at the inlet";
      parameter Q.NumberAbsolute psi_N2=(1 - environment.psi_O2_dry)*
          environment.psi_dry "Mole fraction of N2 at the inlet";
      parameter Q.NumberAbsolute psi_O2=environment.psi_O2_dry*environment.psi_dry
        "Mole fraction of O2 at the inlet";

      output Q.Potential w=caFP.subregions[1, 1, 1].graphite.'e-'.Deltag[1] if
        environment.analysis "Electrical potential";
      output Q.ResistanceElectrical R=w/zI if environment.analysis
        "Measured electrical resistance";
      output Q.ResistanceElectrical R_ex=caFP.L[Axis.x]/(caFP.subregions[1, 1,
          1].graphite.'e-'.sigma*caFP.A[Axis.x]*caFP.subregions[1, 1, 1].graphite.epsilon
          ^1.5) if environment.analysis "Expected electrical resistance";

      CaFPs.CaFP caFP
        annotation (Placement(transformation(extent={{50,-10},{70,10}})));

      // Conditions
      Conditions.ByConnector.BoundaryBus.Single.Sink caBC[caFP.n_y, caFP.n_z](
          each graphite(
          'inclC+'=true,
          'incle-'=true,
          redeclare Conditions.ByConnector.ThermalDiffusive.Single.Temperature
            'C+'(set(y=environment.T)),
          'e-'(materialSet(y=0)))) annotation (Placement(transformation(
            extent={{10,10},{-10,-10}},
            rotation=90,
            origin={84,0})));

      Conditions.ByConnector.BoundaryBus.Single.Source anBC[caFP.n_y, caFP.n_z]
        (
        each gas(
          inclH2O=true,
          inclN2=false,
          inclO2=true,
          H2O(
            redeclare function materialSpec =
                Conditions.ByConnector.Boundary.Single.Material.pressure,
            materialSet(y=environment.p_H2O),
            redeclare function thermalSpec =
                Conditions.ByConnector.Boundary.Single.Thermal.heatRate,
            thermalSet(y=0)),
          N2(
            redeclare function materialSpec =
                Conditions.ByConnector.Boundary.Single.Material.pressure,
            materialSet(y=environment.p*psi_N2),
            redeclare function thermalSpec =
                Conditions.ByConnector.Boundary.Single.Thermal.heatRate,
            thermalSet(y=0)),
          O2(
            materialSet(y=zI/2),
            redeclare function thermalSpec =
                Conditions.ByConnector.Boundary.Single.Thermal.heatRate,
            thermalSet(y=0))),
        each graphite('incle-'=true, 'e-'(materialSet(y=-zI), thermalSet(y=
                  environment.T))),
        each liquid(inclH2O=true, H2O(redeclare function thermalSpec =
                Conditions.ByConnector.Boundary.Single.Thermal.heatRate,
              thermalSet(y=0)))) annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=270,
            origin={36,0})));

      Conditions.ByConnector.BoundaryBus.Single.Source caSource[caFP.n_x, caFP.n_z]
        (each gas(
          inclH2O=true,
          inclN2=false,
          inclO2=true,
          H2O(materialSet(y=-testConditions.Ndot_H2O_ca), thermalSet(y=
                  environment.T)),
          N2(materialSet(y=-testConditions.Ndot_N2), thermalSet(y=environment.T)),

          O2(materialSet(y=-testConditions.Ndot_O2), thermalSet(y=environment.T))))
        annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=180,
            origin={60,-24})));

      Conditions.ByConnector.BoundaryBus.Single.Sink caSink[caFP.n_x, caFP.n_z]
        (gas(
          each inclH2O=true,
          each inclN2=true,
          each inclO2=true,
          H2O(materialSet(y=fill(
                      environment.p,
                      caFP.n_x,
                      caFP.n_z) - caSink.gas.N2.p - caSink.gas.O2.p)),
          N2(materialSet(y=caSink.gas.H2O.boundary.Ndot .* caFP.subregions[:,
                  caFP.n_y, :].gas.H2O.v ./ caFP.subregions[:, caFP.n_y, :].gas.N2.v),
              redeclare each function materialSpec =
                Conditions.ByConnector.Boundary.Single.Material.current),
          O2(materialSet(y=caSink.gas.H2O.boundary.Ndot .* caFP.subregions[:,
                  caFP.n_y, :].gas.H2O.v ./ caFP.subregions[:, caFP.n_y, :].gas.O2.v),
              redeclare each function materialSpec =
                Conditions.ByConnector.Boundary.Single.Material.current)),
          each liquid(inclH2O=true, H2O(materialSet(y=environment.p))))
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={60,24})));

      inner Conditions.Environment environment(
        analysis=true,
        T=333.15*U.K,
        p=U.from_kPag(48.3),
        RH=0.8) "Environmental conditions"
        annotation (Placement(transformation(extent={{-30,30},{-10,50}})));

      Modelica.Blocks.Sources.Ramp currentSet(
        height=100*U.A,
        duration=300,
        offset=U.mA)
        annotation (Placement(transformation(extent={{-20,-60},{0,-40}})));

    protected
      Connectors.RealOutputInternal zI(unit="N/T") "Electrical current"
        annotation (Placement(transformation(extent={{0,-60},{20,-40}})));
    public
      Assemblies.Cells.Examples.TestConditions testConditions(I_ca=2*zI, I_an=
            1.5*zI) "Test conditions" annotation (Dialog,Placement(
            transformation(extent={{10,30},{30,50}})));
    equation

      connect(caSource.boundary, caFP.yNegative) annotation (Line(
          points={{60,-20},{60,-10}},
          color={0,0,240},
          thickness=0.5,
          smooth=Smooth.None));
      connect(caSink.boundary, caFP.yPositive) annotation (Line(
          points={{60,20},{60,10}},
          color={0,0,240},
          thickness=0.5,
          smooth=Smooth.None));
      connect(currentSet.y, zI) annotation (Line(
          points={{1,-50},{10,-50}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(anBC.boundary, caFP.xNegative) annotation (Line(
          points={{40,0},{50,0}},
          color={0,0,240},
          thickness=0.5,
          smooth=Smooth.None));
      connect(caFP.xPositive, caBC.boundary) annotation (Line(
          points={{70,0},{80,0}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));
      annotation (
        experiment(StopTime=350),
        Commands(file="Resources/Scripts/Dymola/Regions.Examples.CaFP.mos"
            "Regions.Examples.CaFP.mos", file=
              "Resources/Scripts/Dymola/Regions.Examples.CaFP-states.mos"
            "Regions.Examples.CaFP-states.mos"),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics));
    end CaFP;

    model FPtoFP "Test one flow plate to the other"
      extends Modelica.Icons.Example;

      parameter Boolean inclLiq=true "Include liquid H2O";
      parameter Q.NumberAbsolute psi_H2O=environment.psi_H2O
        "Mole fraction of H2O at the inlet";
      parameter Q.NumberAbsolute psi_H2=environment.psi_dry
        "Mole fraction of H2 at the inlet";
      parameter Q.NumberAbsolute psi_O2=environment.psi_O2_dry*environment.psi_dry
        "Mole fraction of O2 at the inlet";
      parameter Q.NumberAbsolute psi_N2=(1 - environment.psi_O2_dry)*
          environment.psi_dry "Mole fraction of N2 at the inlet";
      output Q.Number J_Apercm2=zI*U.cm^2/(caFP.A[Axis.x]*U.A)
        "Electrical current density, in A/cm2";
      output Q.Potential w=anFP.subregions[1, 1, 1].graphite.'e-'.g_boundaries[
          1, Side.n] - caFP.subregions[end, 1, 1].graphite.'e-'.g_boundaries[1,
          Side.p] if environment.analysis "Electrical potential";

      parameter Q.Length L_y[:]={8}*U.cm "Lengths in the y direction";
      parameter Q.Length L_z[:]={6.25}*U.cm "Lengths in the z direction";

      // Layers
      AnFPs.AnFP anFP(
        final L_y=L_y,
        final L_z=L_z,
        subregions(each liquid(inclH2O=inclLiq)))
        annotation (Placement(transformation(extent={{-70,-10},{-50,10}})));
      AnGDLs.AnGDL anGDL(
        final L_y=L_y,
        final L_z=L_z,
        subregions(each liquid(inclH2O=inclLiq)))
        annotation (Placement(transformation(extent={{-50,-10},{-30,10}})));
      AnCLs.AnCL anCL(
        final L_y=L_y,
        final L_z=L_z,
        subregions(each liquid(inclH2O=inclLiq)))
        annotation (Placement(transformation(extent={{-30,-10},{-10,10}})));
      PEMs.PEM PEM(final L_y=L_y, final L_z=L_z)
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
      CaCLs.CaCL caCL(
        final L_y=L_y,
        final L_z=L_z,
        subregions(each liquid(inclH2O=inclLiq)))
        annotation (Placement(transformation(extent={{10,-10},{30,10}})));
      CaGDLs.CaGDL caGDL(
        final L_y=L_y,
        final L_z=L_z,
        subregions(each liquid(inclH2O=inclLiq)))
        annotation (Placement(transformation(extent={{30,-10},{50,10}})));
      CaFPs.CaFP caFP(
        final L_y=L_y,
        final L_z=L_z,
        subregions(each liquid(inclH2O=inclLiq)))
        annotation (Placement(transformation(extent={{50,-10},{70,10}})));

      // Conditions
      Conditions.ByConnector.BoundaryBus.Single.Sink anBC[anFP.n_y, anFP.n_z](
          each graphite(
          'inclC+'=true,
          'incle-'=true,
          redeclare Conditions.ByConnector.ThermalDiffusive.Single.Temperature
            'C+'(set(y=environment.T)),
          'e-'(materialSet(y=0)))) annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=90,
            origin={-84,0})));

      Conditions.ByConnector.BoundaryBus.Single.Source anSource[anFP.n_x, anFP.n_z]
        (each gas(
          inclH2=true,
          inclH2O=true,
          H2(materialSet(y=-testConditions.Ndot_H2), thermalSet(y=environment.T)),

          H2O(materialSet(y=-testConditions.Ndot_H2O_an), thermalSet(y=
                  environment.T)))) annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=180,
            origin={-60,-24})));

      Conditions.ByConnector.BoundaryBus.Single.Sink anSink[anFP.n_x, anFP.n_z]
        (gas(
          each inclH2=true,
          each inclH2O=true,
          H2(redeclare each function materialSpec =
                Conditions.ByConnector.Boundary.Single.Material.current,
              materialSet(y=anSink.gas.H2O.boundary.Ndot .* anFP.subregions[:,
                  anFP.n_y, :].gas.H2O.v ./ anFP.subregions[:, anFP.n_y, :].gas.H2.v)),

          H2O(materialSet(y=fill(
                      environment.p,
                      anFP.n_x,
                      anFP.n_z) - anSink.gas.H2.p))), each liquid(inclH2O=
              inclLiq, H2O(materialSet(y=environment.p)))) annotation (
          Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={-60,24})));

      Conditions.ByConnector.BoundaryBus.Single.Source caBC[caFP.n_y, caFP.n_z]
        (each graphite(
          'inclC+'=true,
          'incle-'=true,
          'C+'(set(y=environment.T)),
          'e-'(
            redeclare function materialSpec =
                Conditions.ByConnector.Boundary.Single.Material.current,
            materialSet(y=-zI),
            thermalSet(y=environment.T)))) annotation (Placement(transformation(
            extent={{10,10},{-10,-10}},
            rotation=90,
            origin={84,0})));

      Conditions.ByConnector.BoundaryBus.Single.Source caSource[caFP.n_x, caFP.n_z]
        (each gas(
          inclH2O=true,
          inclN2=true,
          inclO2=true,
          H2O(materialSet(y=-testConditions.Ndot_H2O_ca), thermalSet(y=
                  environment.T)),
          N2(materialSet(y=-testConditions.Ndot_N2), thermalSet(y=environment.T)),

          O2(materialSet(y=-testConditions.Ndot_O2), thermalSet(y=environment.T))))
        annotation (Placement(transformation(
            extent={{10,-10},{-10,10}},
            rotation=180,
            origin={60,-24})));

      Conditions.ByConnector.BoundaryBus.Single.Sink caSink[caFP.n_x, caFP.n_z]
        (gas(
          each inclH2O=true,
          each inclN2=true,
          each inclO2=true,
          H2O(materialSet(y=fill(
                      environment.p,
                      caFP.n_x,
                      caFP.n_z) - caSink.gas.N2.p - caSink.gas.O2.p)),
          N2(redeclare function materialSpec =
                Conditions.ByConnector.Boundary.Single.Material.current,
              materialSet(y=caSink.gas.H2O.boundary.Ndot .* caFP.subregions[:,
                  caFP.n_y, :].gas.H2O.v ./ caFP.subregions[:, caFP.n_y, :].gas.N2.v)),

          O2(redeclare function materialSpec =
                Conditions.ByConnector.Boundary.Single.Material.current,
              materialSet(y=caSink.gas.H2O.boundary.Ndot .* caFP.subregions[:,
                  caFP.n_y, :].gas.H2O.v ./ caFP.subregions[:, caFP.n_y, :].gas.O2.v))),
          each liquid(inclH2O=inclLiq,H2O(materialSet(y=environment.p))))
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=0,
            origin={60,24})));

      Modelica.Blocks.Sources.Ramp currentSet(
        offset=U.mA,
        startTime=50,
        duration=600,
        height=100*U.A)
        annotation (Placement(transformation(extent={{-20,-60},{0,-40}})));

    protected
      Connectors.RealOutputInternal zI(unit="N/T") "Electrical current"
        annotation (Placement(transformation(extent={{0,-60},{20,-40}})));

      inner Conditions.Environment environment(
        a={0,0,0},
        analysis=true,
        T=333.15*U.K,
        p=U.from_kPag(48.3),
        RH=0.7) "Environmental conditions"
        annotation (Placement(transformation(extent={{-30,30},{-10,50}})));

    public
      Assemblies.Cells.Examples.TestConditions testConditions(I_ca=2*zI, I_an=
            1.5*zI) "Test conditions" annotation (Dialog,Placement(
            transformation(extent={{10,30},{30,50}})));
    equation
      connect(anCL.xPositive, PEM.xNegative) annotation (Line(
          points={{-10,0},{-10,0}},
          color={240,0,0},
          thickness=0.5,
          smooth=Smooth.None));
      connect(PEM.xPositive, caCL.xNegative) annotation (Line(
          points={{10,0},{10,0}},
          color={0,0,240},
          thickness=0.5,
          smooth=Smooth.None));
      connect(caCL.xPositive, caGDL.xNegative) annotation (Line(
          points={{30,0},{30,0}},
          color={0,0,240},
          thickness=0.5,
          smooth=Smooth.None));
      connect(anCL.xNegative, anGDL.xPositive) annotation (Line(
          points={{-30,0},{-30,0}},
          color={240,0,0},
          thickness=0.5,
          smooth=Smooth.None));
      connect(anBC.boundary, anFP.xNegative) annotation (Line(
          points={{-80,0},{-70,0}},
          color={127,127,127},
          pattern=LinePattern.None,
          thickness=0.5,
          smooth=Smooth.None));
      connect(anSource.boundary, anFP.yNegative) annotation (Line(
          points={{-60,-20},{-60,-10}},
          color={240,0,0},
          thickness=0.5,
          smooth=Smooth.None));
      connect(anSink.boundary, anFP.yPositive) annotation (Line(
          points={{-60,20},{-60,10}},
          color={240,0,0},
          thickness=0.5,
          smooth=Smooth.None));
      connect(currentSet.y, zI) annotation (Line(
          points={{1,-50},{10,-50}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(caBC.boundary, caFP.xPositive) annotation (Line(
          points={{80,0},{70,0}},
          color={127,127,127},
          thickness=0.5,
          smooth=Smooth.None));
      connect(anGDL.xNegative, anFP.xPositive) annotation (Line(
          points={{-50,0},{-50,0}},
          color={240,0,0},
          thickness=0.5,
          smooth=Smooth.None));
      connect(caGDL.xPositive, caFP.xNegative) annotation (Line(
          points={{50,0},{50,0}},
          color={0,0,240},
          thickness=0.5,
          smooth=Smooth.None));
      connect(caFP.yNegative, caSource.boundary) annotation (Line(
          points={{60,-10},{60,-20}},
          color={0,0,240},
          thickness=0.5,
          smooth=Smooth.None));
      connect(caSink.boundary, caFP.yPositive) annotation (Line(
          points={{60,20},{60,10}},
          color={0,0,240},
          thickness=0.5,
          smooth=Smooth.None));
      annotation (
        Commands(file="Resources/Scripts/Dymola/Regions.Examples.FPtoFP.mos"
            "Regions.Examples.FPtoFP.mos", file=
              "Resources/Scripts/Dymola/Regions.Examples.FPtoFP-states.mos"
            "Regions.Examples.FPtoFP-states.mos"),
        experiment(StopTime=650, __Dymola_Algorithm="Dassl"),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics),
        __Dymola_experimentSetupOutput);
    end FPtoFP;

  end Examples;

  package AnFPs "Anode flow plates"
    extends Modelica.Icons.Package;

    model AnFP "Anode flow plate"
      import Modelica.Constants.inf;
      // extends FCSys.Icons.Names.Top4;

      extends Region(
        L_x={8}*U.mm,
        L_y={8}*U.cm,
        L_z={6.25}*U.cm,
        final inclTransX=true,
        final inclTransY=true,
        inclTransZ=false,
        redeclare replaceable model Subregion =
            FCSys.Subregions.SubregionNoIonomer (
            common(k_Phi={1e7,inf,1e7},k_Q=1e5),
            gasLiq(k_Phi={inf,1e6,inf},k_Q=inf),
            gas(
              common(k_Phi={inf,inf,inf}),
              k={epsilon/2,11,1/11},
              inclH2=true,
              inclH2O=true,
              H2(
                upstreamX=false,
                Nu_Phi={4,16*A[Axis.z]*epsilon/D^2,4},
                zeta=100*Characteristics.H2.Gas.zeta(),
                T(stateSelect=StateSelect.always)),
              H2O(
                upstreamX=false,
                Nu_Phi={4,16*A[Axis.z]*epsilon/D^2,4},
                zeta=100*Characteristics.H2O.Gas.zeta(),
                I(each stateSelect=StateSelect.always, each fixed=true),
                initEnergy=Init.none)),
            graphite(
              'inclC+'=true,
              'incle-'=true,
              'C+'(theta=U.m*U.K/(95*U.W),epsilon=1 - epsilon),
              'e-'(sigma=U.S/(1.470e-3*U.cm))),
            liquid(
              k={epsilon/2,11,1/11},
              inclH2O=true,
              H2O(
                upstreamX=false,
                Nu_Phi={4,16*A[Axis.z]*epsilon/D^2,4},
                epsilon_IC=1e-5,
                N0=0.1*U.C)),
            volume(inclCapillary=false))) annotation (IconMap(primitivesVisible
            =false));

      parameter Q.NumberAbsolute epsilon(nominal=1) = 0.0625
        "Fraction of volume for the fluid" annotation (Dialog(group="Geometry",
            __Dymola_label="<html>&epsilon;</html>"));

      parameter Q.Length D=1.5*U.mm "Hydraulic diameter of the channel"
        annotation (Dialog(__Dymola_label="<html><i>D</i></html>"));

    protected
      Q.Velocity phi_states_H2[:, :, :](
        each stateSelect=StateSelect.always,
        each start=0,
        each fixed=true) = subregions[:, 2:n_y, :].gas.H2.phi[2] if n_y > 1
        "Forced states for H2";
      // Note:  This avoids dynamic state selection in Dymola 2014.

      outer Conditions.Environment environment "Environmental conditions";

      // Thermal resistivity of some other flow plate materials [Incropera2002, pp. 905 & 907]:
      //                                Stainless steel
      //                                ---------------------------------------------------------------------------------
      //             Aluminium (pure)   AISI 302             AISI 304              AISI 316              AISI 347
      //             ----------------   ------------------   ------------------    -----------------    -----------------
      //                        theta              theta                theta                 theta                theta
      //             c_p*U.kg   *U.W    c_p*U.kg   *U.W      c_p*U.kg   *U.W       c_p*U.kg   *U.W      c_p*U.kg   *U.W
      //             *U.K       /(U.m   *U.K       /(U.m     *U.K       /(U.m      *U.K       /(U.m     *U.K       /(U.m
      //     T/K     /(U.J*m)   *U.K)   /(U.J*m)   *U.K)     /(U.J*m)   *U.K)      /(U.J*m)   *U.K)     /(U.J*m)   *U.K)
      //     ----    --------   -----   --------   ------    --------   -------    --------   ------    --------   ------
      //     100     482        1/302                        272        1/9.2
      //     200     798        1/237                        402        1/12.6
      //     300     903        1/237   480        1/15.1    477        1/14.9     468        1/13.4    480        1/14.2
      //     400     949        1/240   512        1/17.3    515        1/16.6     504        1/15.2    513        1/15.8
      //     600     1033       1/231   559        1/20.0    557        1/19.8     550        1/18.3    559        1/18.9
      //     800     1146       1/218   585        1/22.8    582        1/22.6     576        1/21.3    585        1/21.9
      //     1000                       606        1/25.4    611        1/25.4     602        1/24.2    606        1/24.7

      // Electrical resistivities:
      //     Aluminium
      //     (http://en.wikipedia.org/wiki/Electrical_resistivity, 2008):
      //         2.82e-8 ohm.m
      //     Graphite (http://hypertextbook.com/facts/2004/AfricaBelgrave.shtml):
      //         7.837e-6 to 41e-6 ohm.m
      //     Copper (http://en.wikipedia.org/wiki/Electrical_resistivity, 2008):
      //         1.72e-8 ohm.m
      //     Stainless steel AISI 304
      //     (http://hypertextbook.com/facts/2006/UmranUgur.shtml, 2008):
      //         6.897e-7 ohm.m

      annotation (
        Documentation(info="<html>
<p>This model represents the anode flow plate of a PEMFC.
The x axis extends from the anode to the cathode.  Fluid is considered to travel

in the y direction, with the associated length factor (<i>k</i><sub>y</sub>) greater than one (by default)
to represent a serpentine channel.
The model is
bidirectional, meaning that either <code>yNegative</code> or <code>yPositive</code> can be
used as the inlet.  By default, the cross-sectional area in the yz plane is 50 cm<sup>2</sup>.</p>

<p>The solid and the fluid phases are assumed to exist in the same subregions, 
even though a 
typical flow plate is impermeable
to the fluid (except for the channel).  
In theory, it is possible
to discretize the flow plate into smaller subregions for the bulk solid, lands, and valleys.  However,
this would significantly increase the mathematical size of the model.  Currently, that level of detail is best left to
computational fluid dynamics.</p>
 
<p>The x axis-component of the transport factor (<i>k</i><sub>x</sub>) for the gas and the liquid should
generally be less than one because the transport distance into/out of the GDL is less that half the 
thickness of the flow plate. It is equal to the product of two ratios:<ol>
<li>the depth of the channels to the thickness of the flow plate</li>
<li>the product of the total area of the flow plate in the yz plane (land + valleys) and the fraction of
the total volume available for the fluid (&epsilon;) to the area of the valleys in the yz plane</li>
</ol> The default is &epsilon;/2.</p>

<p>See <a href=\"modelica://FCSys.Species.'C+'.Graphite.Fixed\">Species.'C+'.Graphite.Fixed</a>
regarding the default specific heat capacity.  The default thermal resistivity
of the carbon (&theta; = <code>U.m*U.K/(95*U.W)</code>) and the
electrical conductivity (&sigma; = <code>U.S/(1.470e-3*U.cm)</code>)
are that of Entegris/Poco Graphite AXF-5Q
[<a href=\"modelica://FCSys.UsersGuide.References.Entegris2012\">Entegris2012</a>].
There is additional data in the
text layer of this model.</p>

<p>For more information, please see the
    <a href=\"modelica://FCSys.Regions.Region\">Region</a> model.</p>
</html>"),
        Icon(coordinateSystem(
            preserveAspectRatio=true,
            extent={{-100,-100},{100,100}},
            initialScale=0.1), graphics={Rectangle(
                  extent={{-98,62},{98,98}},
                  fillColor={255,255,255},
                  visible=not inclTransY,
                  fillPattern=FillPattern.Solid,
                  pattern=LinePattern.None),Rectangle(
                  extent={{-76.648,66.211},{-119.073,52.0689}},
                  fillPattern=FillPattern.HorizontalCylinder,
                  rotation=45,
                  fillColor={135,135,135},
                  origin={111.017,77.3801},
                  pattern=LinePattern.None,
                  lineColor={95,95,95}),Rectangle(
                  extent={{-20,40},{0,-60}},
                  lineColor={95,95,95},
                  fillPattern=FillPattern.VerticalCylinder,
                  fillColor={135,135,135}),Polygon(
                  points={{20,0},{42,0},{42,80},{-42,80},{-42,0},{-20,0},{-20,
                40},{0,60},{20,60},{20,0}},
                  smooth=Smooth.None,
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid,
                  pattern=LinePattern.None),Polygon(
                  points={{20,0},{42,0},{42,-80},{-42,-80},{-42,0},{-20,0},{-20,
                -60},{0,-60},{20,-40},{20,0}},
                  smooth=Smooth.None,
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid,
                  pattern=LinePattern.None),Rectangle(extent={{-20,40},{0,-60}},
              lineColor={0,0,0}),Polygon(
                  points={{-20,40},{0,60},{20,60},{0,40},{-20,40}},
                  lineColor={0,0,0},
                  smooth=Smooth.None),Polygon(
                  points={{20,60},{0,40},{0,-60},{20,-40},{20,60}},
                  lineColor={0,0,0},
                  fillColor={95,95,95},
                  fillPattern=FillPattern.Solid),Polygon(
                  points={{16,48},{4,36},{4,32},{14,42},{14,36},{4,26},{4,12},{
                16,24},{16,28},{6,18},{6,24},{16,34},{16,48}},
                  smooth=Smooth.None,
                  fillColor={0,0,0},
                  fillPattern=FillPattern.Solid,
                  pattern=LinePattern.None,
                  lineColor={0,0,0}),Polygon(
                  points={{16,28},{4,16},{4,12},{14,22},{14,16},{4,6},{4,-8},{
                16,4},{16,8},{6,-2},{6,4},{16,14},{16,28}},
                  smooth=Smooth.None,
                  fillColor={0,0,0},
                  fillPattern=FillPattern.Solid,
                  pattern=LinePattern.None,
                  lineColor={0,0,0}),Polygon(
                  points={{16,8},{4,-4},{4,-8},{14,2},{14,-4},{4,-14},{4,-28},{
                16,-16},{16,-12},{6,-22},{6,-16},{16,-6},{16,8}},
                  smooth=Smooth.None,
                  fillColor={0,0,0},
                  fillPattern=FillPattern.Solid,
                  pattern=LinePattern.None,
                  lineColor={0,0,0}),Polygon(
                  points={{16,-12},{4,-24},{4,-28},{14,-18},{14,-24},{4,-34},{4,
                -48},{16,-36},{16,-32},{6,-42},{6,-36},{16,-26},{16,-12}},
                  smooth=Smooth.None,
                  fillColor={0,0,0},
                  fillPattern=FillPattern.Solid,
                  pattern=LinePattern.None,
                  lineColor={0,0,0}),Line(
                  points={{10,0},{100,0}},
                  color={240,0,0},
                  visible=inclTransX,
                  thickness=0.5),Line(
                  points={{-20,0},{-100,0}},
                  color={127,127,127},
                  visible=inclTransX,
                  thickness=0.5),Line(
                  points={{0,-60},{0,-100}},
                  color={240,0,0},
                  visible=inclTransY,
                  smooth=Smooth.None,
                  thickness=0.5),Line(
                  points={{20,20},{50,50}},
                  color={253,52,56},
                  visible=inclTransZ,
                  smooth=Smooth.None,
                  thickness=0.5),Line(
                  points={{-50,-50},{-10,-10}},
                  color={253,52,56},
                  visible=inclTransZ,
                  smooth=Smooth.None,
                  thickness=0.5),Ellipse(
                  extent={{-4,52},{4,48}},
                  lineColor={135,135,135},
                  fillColor={240,0,0},
                  visible=inclTransY,
                  fillPattern=FillPattern.Sphere),Line(
                  points={{0,100},{0,50}},
                  color={240,0,0},
                  visible=inclTransY,
                  smooth=Smooth.None,
                  thickness=0.5),Text(
                  extent={{-100,60},{100,100}},
                  textString="%name",
                  visible=not inclTransY,
                  lineColor={0,0,0})}),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-60,-60},{
                40,40}}), graphics));
    end AnFP;

  end AnFPs;

  package AnGDLs "Anode gas diffusion layers"
    extends Modelica.Icons.Package;

    model AnGDL "Anode gas diffusion layer"
      import Modelica.Constants.inf;
      // extends FCSys.Icons.Names.Top4;

      // Note:  Extensions of AnGDL should be placed directly in the AnGDLs
      // package rather than subpackages (e.g., by manufacturer) so that
      // __Dymola_choicesFromPackage can be used.  Dymola 7.4 launches the
      // parameter dialogs too slowly when __Dymola_choicesAllMatching is
      // used.

      extends Region(
        L_x={0.3}*U.mm,
        L_y={8}*U.cm,
        L_z={6.25}*U.cm,
        final inclTransX=true,
        inclTransY=false,
        inclTransZ=false,
        redeclare replaceable model Subregion =
            FCSys.Subregions.SubregionNoIonomer (
            common(k_Q=0),
            gasLiq(k_Phi={2,2,2}, k_Q=inf),
            gas(
              common(k_Q=inf),
              k=fill(epsilon^(-0.5), 3),
              inclH2=true,
              inclH2O=true,
              H2(
                initEnergy=Init.none,
                upstreamX=false,
                phi(each stateSelect=StateSelect.always, each fixed=true),
                zeta=100*Characteristics.H2.Gas.zeta(),
                final eta=0),
              H2O(
                initEnergy=Init.none,
                upstreamX=false,
                phi(each stateSelect=StateSelect.always, each fixed=true),
                zeta=100*Characteristics.H2O.Gas.zeta(),
                final eta=0)),
            graphite(
              k=fill((1 - epsilon)^(-0.5), 3),
              'inclC+'=true,
              'incle-'=true,
              'C+'(theta=U.m*U.K/(1.18*U.W)*(0.5*(1 - epsilon))^1.5,final
                  epsilon=1 - epsilon),
              'e-'(sigma=40*U.S/(12*U.cm)/(1 - epsilon)^1.5)),
            liquid(
              k=fill(epsilon^(-0.5), 3),
              inclH2O=true,
              H2O(
                upstreamX=false,
                epsilon_IC=1e-5,
                final eta=0,
                mu=100*Characteristics.H2O.Liquid.mu(),
                phi(each stateSelect=StateSelect.always, each fixed=true),
                each initEnergy=Init.none,
                T(each stateSelect=StateSelect.default))))) annotation (IconMap(
            primitivesVisible=false));

      // Note:  The fluid species have zero fluidity (eta=0) so that the transverse
      // velocity is zero at the interface with the flow plate.  That condition
      // is necessary to produce the appropriate pressure loss down the channel.

      // See the documentation layer of Phases.PartialPhase regarding the
      // settings of k for each phase.

      parameter Q.NumberAbsolute epsilon(nominal=1) = 0.88 "Porosity"
        annotation (Dialog(group="Geometry", __Dymola_label=
              "<html>&epsilon;</html>"));

    protected
      outer Conditions.Environment environment "Environmental conditions";
      annotation (Documentation(info="<html>
<p>This model represents the anode gas diffusion layer of a PEMFC.
The x axis extends from the anode to the cathode.

By default, the cross-sectional area in the yz plane is 50 cm<sup>2</sup>.</p>

<p>The default porosity (&epsilon; = 0.88) is that of SGL Carbon Group Sigracet&reg; 10 BA and 25 BA GDLs.

The porosity of a GDL may be lower than specified due to compression (e.g., 0.4 according to

[<a href=\"modelica://FCSys.UsersGuide.References.Bernardi1992\">Bernardi1992</a>, p. 2483], although
that reference may be outdated).
  The default thermal conductivity of the carbon (&theta; = <code>U.m*U.K/(1.18*U.W)</code>)
   represents a compressed Sigracet&reg; 10 BA gas diffusion layer
  [<a href=\"modelica://FCSys.UsersGuide.References.Nitta2008\">Nitta2008</a>].  The default
  electrical conductivity
  is also for Sigracet&reg; 10 BA [<a href=\"modelica://FCSys.UsersGuide.References.SGL2007\">SGL2007</a>].</p>

<p>Default assumptions (may be adjusted):
<ol>
<li>All of the species have the same temperature, even in different phases.</li> 
</ol></p>

<p>For more information, please see the
    <a href=\"modelica://FCSys.Regions.Region\">Region</a> model.</p></html>"),
          Icon(coordinateSystem(
            preserveAspectRatio=true,
            extent={{-100,-100},{100,100}},
            initialScale=0.1), graphics={Rectangle(
                  extent={{-98,62},{98,98}},
                  fillColor={255,255,255},
                  visible=not inclTransY,
                  fillPattern=FillPattern.Solid,
                  pattern=LinePattern.None),Rectangle(
                  extent={{-78.7855,18.6813},{-50.5004,-23.7455}},
                  lineColor={64,64,64},
                  fillColor={127,127,127},
                  rotation=-45,
                  fillPattern=FillPattern.VerticalCylinder,
                  origin={42.5001,11.0805}),Rectangle(
                  extent={{-40,40},{0,-60}},
                  lineColor={64,64,64},
                  fillColor={127,127,127},
                  fillPattern=FillPattern.VerticalCylinder),Polygon(
                  points={{20,0},{42,0},{42,80},{-42,80},{-42,0},{-20,0},{-20,
                40},{0,60},{20,60},{20,0}},
                  smooth=Smooth.None,
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid,
                  pattern=LinePattern.None),Polygon(
                  points={{20,0},{42,0},{42,-80},{-42,-80},{-42,0},{-20,0},{-20,
                -60},{0,-60},{20,-40},{20,0}},
                  smooth=Smooth.None,
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid,
                  pattern=LinePattern.None),Polygon(
                  points={{0,40},{20,60},{20,-40},{0,-60},{0,40}},
                  lineColor={0,0,0},
                  smooth=Smooth.None,
                  fillPattern=FillPattern.Solid,
                  fillColor={64,64,64}),Rectangle(extent={{-20,40},{0,-60}},
              lineColor={0,0,0}),Polygon(
                  points={{0,60},{20,60},{0,40},{-20,40},{0,60}},
                  lineColor={0,0,0},
                  smooth=Smooth.None),Line(
                  points={{-20,0},{-100,0}},
                  color={240,0,0},
                  visible=inclTransX,
                  thickness=0.5),Line(
                  points={{10,0},{100,0}},
                  color={240,0,0},
                  visible=inclTransX,
                  thickness=0.5),Line(
                  points={{0,-60},{0,-100}},
                  color={253,52,56},
                  visible=inclTransY,
                  smooth=Smooth.None,
                  thickness=0.5),Line(
                  points={{0,100},{0,50}},
                  color={253,52,56},
                  visible=inclTransY,
                  smooth=Smooth.None,
                  thickness=0.5),Line(
                  points={{-50,-50},{-10,-10}},
                  color={253,52,56},
                  visible=inclTransZ,
                  smooth=Smooth.None,
                  thickness=0.5),Line(
                  points={{20,20},{50,50}},
                  color={253,52,56},
                  visible=inclTransZ,
                  smooth=Smooth.None,
                  thickness=0.5),Text(
                  extent={{-100,60},{100,100}},
                  textString="%name",
                  visible=not inclTransY,
                  lineColor={0,0,0})}));

    end AnGDL;

    model Sigracet10BA "<html>SGL Carbon Group Sigracet&reg; 10 BA</html>"
      extends AnGDL(
        L_x={0.400}*U.mm,
        epsilon=0.88,
        Subregion(graphite('e-'(sigma=40*U.S/(12*U.cm)))));
      // Additional properties not incorporated [SGL2007]:
      //     Diffusivity:  L = 0.400 mm, p = 0.85 m/s (for air) => D = P*L = 340 mm2/s
      //     Density:  (85 g/m2)/(0.400 mm)/0.88 = 212.5 kg/m3

      annotation (defaultComponentName="anGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References.SGL2007\">SGL2007</a>].</p>

  <p>For more information, please see the
    <a href=\"modelica://FCSys.Regions.AnGDLs.AnGDL\">AnGDL</a> model.</p></html>"));

    end Sigracet10BA;

    model Sigracet10BB "<html>SGL Carbon Group Sigracet&reg; 10 BB</html>"
      extends AnGDL(
        L_x={0.420}*U.mm,
        epsilon=0.84,
        Subregion(graphite('e-'(sigma=42*U.S/(15*U.cm)))));
      // Additional properties not incorporated [SGL2007]:
      //     Diffusivity:  L = 0.420 mm, p = 0.03 m/s (for air) => D = P*L = 12.6 mm2/s
      //     Density:  (125 g/m2)/(0.420 mm) = 297.62 kg/m3
      annotation (defaultComponentName="anGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References.SGL2007\">SGL2007</a>].</p>

  <p>For more information, please see the
    <a href=\"modelica://FCSys.Regions.AnGDLs.AnGDL\">AnGDL</a> model.</p></html>"));

    end Sigracet10BB;

    model Sigracet10BC "<html>SGL Carbon Group Sigracet&reg; 10 BC</html>"
      extends AnGDL(
        L_x={0.420}*U.mm,
        epsilon=0.82,
        Subregion(graphite('e-'(sigma=42*U.S/(16*U.cm)))));
      // Additional properties not incorporated [SGL2007]:
      //     Diffusivity:  L = 0.420 mm, p = 0.0145 m/s (for air) => D = P*L = 6.09 mm2/s
      //     Density:  (135 g/m2)/(0.420 mm) = 321.43 kg/m3
      annotation (defaultComponentName="anGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References.SGL2007\">SGL2007</a>].</p>

  <p>For more information, please see the
    <a href=\"modelica://FCSys.Regions.AnGDLs.AnGDL\">AnGDL</a> model.</p></html>"));

    end Sigracet10BC;

    model Sigracet24BA "<html>SGL Carbon Group Sigracet&reg; 24 BA</html>"
      extends AnGDL(
        L_x={0.190}*U.mm,
        epsilon=0.84,
        Subregion(graphite('e-'(sigma=19*U.S/(10*U.cm)))));
      // Additional properties not incorporated [SGL2004]:
      //     Diffusivity:  L = 0.190 mm, p = 0.30 m/s (for air) => D = P*L = 57 mm2/s
      //     Density:  (54 g/m2)/(0.190 mm) = 284.21 kg/m3
      annotation (defaultComponentName="anGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References.SGL2004\">SGL2004</a>].</p>

  <p>For more information, please see the
    <a href=\"modelica://FCSys.Regions.AnGDLs.AnGDL\">AnGDL</a> model.</p></html>"));

    end Sigracet24BA;

    model Sigracet24BC "<html>SGL Carbon Group Sigracet&reg; 24 BC</html>"
      extends AnGDL(
        L_x={0.235}*U.mm,
        epsilon=0.76,
        Subregion(graphite('e-'(sigma=23.5*U.S/(11*U.cm)))));
      // Additional properties not incorporated [SGL2004]:
      //     Diffusivity:  L = 0.235 mm, p = 0.0045 m/s (for air) => D = P*L = 1.0575 mm2/s
      //     Density:  (100 g/m2)/(0.235 mm) = 425.53 kg/m3
      annotation (defaultComponentName="anGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References.SGL2004\">SGL2004</a>].</p>

  <p>For more information, please see the
    <a href=\"modelica://FCSys.Regions.AnGDLs.AnGDL\">AnGDL</a> model.</p></html>"));

    end Sigracet24BC;

    model Sigracet25BA "<html>SGL Carbon Group Sigracet&reg; 25 BA</html>"
      extends AnGDL(
        L_x={0.190}*U.mm,
        epsilon=0.88,
        Subregion(graphite('e-'(sigma=19*U.S/(10*U.cm)))));
      // Additional properties not incorporated [SGL2004]:
      //     Diffusivity:  L = 0.190 mm, p = 0.90 m/s (for air) => D = P*L = 171 mm2/s
      //     Density:  (40 g/m2)/(0.190 mm) = 210.53 kg/m3
      annotation (defaultComponentName="anGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References.SGL2004\">SGL2004</a>].</p>

  <p>For more information, please see the
    <a href=\"modelica://FCSys.Regions.AnGDLs.AnGDL\">AnGDL</a> model.</p></html>"));

    end Sigracet25BA;

    model Sigracet25BC "<html>SGL Carbon Group Sigracet&reg; 25 BC</html>"
      extends AnGDL(
        L_x={0.235}*U.mm,
        epsilon=0.80,
        Subregion(graphite('e-'(sigma=23.5*U.S/(12*U.cm)))));
      // Additional properties not incorporated [SGL2004]:
      //     Diffusivity:  L = 0.235 mm, p = 0.008 m/s (for air) => D = P*L = 1.88 mm2/s
      //     Density:  (86 g/m2)/(0.235 mm) = 365.96 kg/m3
      annotation (defaultComponentName="anGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References.SGL2004\">SGL2004</a>].</p>

  <p>For more information, please see the
    <a href=\"modelica://FCSys.Regions.AnGDLs.AnGDL\">AnGDL</a> model.</p></html>"));

    end Sigracet25BC;

    model TorayTGPH030 "Toray Industries TGP-H-030"
      extends AnGDL(
        L_x={0.11}*U.mm,
        epsilon=0.80,
        Subregion(graphite('C+'(theta=U.m*U.K/(1.7*U.W)),'e-'(sigma=U.S/(0.80*U.mm)))));
      // Additional properties not incorporated [Toray2010]:
      //     Diffusivity:  L = 0.110 mm, P/p = 2500 ml.mm/(cm2.hr.mmAq) = 0.70814e-3 m/(s.kPa)
      //         => D = P*L = 7.89e-6 m2/s, assuming p = 101.325 kPa
      //     Bulk density:  0.44 g/cm3
      annotation (defaultComponentName="anGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References.Toray2010\">Toray2010</a>].  The electrical conductivity (&sigma;) is based on the
  through-plane value of resistivity.  The thermal conductivity is not listed but is
  taken to be the same as for TGP-H-060, TGP-H-090, and TGP-H-120.</p>

  <p>For more information, please see the
    <a href=\"modelica://FCSys.Regions.AnGDLs.AnGDL\">AnGDL</a> model.</p></html>"));

    end TorayTGPH030;

    model TorayTGPH060 "Toray Industries TGP-H-060"
      extends AnGDL(
        L_x={0.19}*U.mm,
        epsilon=0.78,
        Subregion(graphite('C+'(theta=U.m*U.K/(1.7*U.W)),'e-'(sigma=U.S/(0.80*U.mm)))));

      // Additional properties not incorporated [Toray2010]:
      //     Diffusivity: L = 0.190 mm, P/p = 1900 ml.mm/(cm2.hr.mmAq) = 0.53818e-3 m/(s.kPa)
      //         => D = P*L = 10.36e-6 m2/s, assuming p = 101.325 kPa
      //     Bulk density:  0.44 g/cm3
      annotation (defaultComponentName="anGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References.Toray2010\">Toray2010</a>].  The electrical conductivity (&sigma;)  is based on the
  through-plane value of resistivity.  The thermal conductivity is through the plane at
  room temperature.</p>

  <p>For more information, please see the
    <a href=\"modelica://FCSys.Regions.AnGDLs.AnGDL\">AnGDL</a> model.</p></html>"));

    end TorayTGPH060;

    model TorayTGPH090 "Toray Industries TGP-H-090"
      extends AnGDL(
        L_x={0.28}*U.mm,
        epsilon=0.78,
        Subregion(graphite('C+'(theta=U.m*U.K/(1.7*U.W)),'e-'(sigma=U.S/(0.80*U.mm)))));

      // Additional properties not incorporated [Toray2010]:
      //     Diffusivity: L = 0.280 mm, P/p = 1700 ml.mm/(cm2.hr.mmAq) = 0.48153e-3 m/(s.kPa)
      //         => D = P*L = 13.66e-6 m2/s, assuming p = 101.325 kPa
      //     Bulk density:  0.44 g/cm3
      annotation (defaultComponentName="anGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References.Toray2010\">Toray2010</a>].  The electrical conductivity  (&sigma;) is based on the
  through-plane value of resistivity.  The thermal conductivity is through the plane at
  room temperature.</p>

  <p>For more information, please see the
    <a href=\"modelica://FCSys.Regions.AnGDLs.AnGDL\">AnGDL</a> model.</p></html>"));

    end TorayTGPH090;

    model TorayTGPH120 "Toray Industries TGP-H-030"
      extends AnGDL(
        L_x={0.37}*U.mm,
        epsilon=0.78,
        Subregion(graphite('C+'(theta=U.m*U.K/(1.7*U.W)),'e-'(sigma=U.S/(0.80*U.mm)))));

      // Additional properties not incorporated [Toray2010]:
      //     Diffusivity: L = 0.370 mm, P/p = 1500 ml.mm/(cm2.hr.mmAq) = 0.42488e-3 m/(s.kPa)
      //         => D = P*L = 15.93e-6 m2/s, assuming p = 101.325 kPa
      //     Bulk density:  0.44 g/cm3
      annotation (defaultComponentName="anGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References.Toray2010\">Toray2010</a>].  The electrical conductivity  (&sigma;) is based on the
  through-plane value of resistivity.  The thermal conductivity is through the plane at
  room temperature.</p>

  <p>For more information, please see the
    <a href=\"modelica://FCSys.Regions.AnGDLs.AnGDL\">AnGDL</a> model.</p></html>"));

    end TorayTGPH120;

  end AnGDLs;

  package AnCLs "Anode catalyst layers"
    extends Modelica.Icons.Package;

    model AnCL "Anode catalyst layer"
      import Modelica.Constants.inf;
      // extends FCSys.Icons.Names.Top4;

      extends Region(
        L_x={28.7}*U.um,
        L_y={8}*U.cm,
        L_z={6.25}*U.cm,
        final inclTransX=true,
        inclTransY=false,
        inclTransZ=false,
        redeclare replaceable model Subregion = Subregions.Subregion (
            common(k_Q=0),
            gasLiq(k_Phi={2,2,2}, k_Q=inf),
            gas(
              common(k_Q=inf),
              k=fill(epsilon^(-0.5), 3),
              inclH2=true,
              inclH2O=true,
              H2(
                initEnergy=Init.none,
                upstreamX=false,
                zeta=100*Characteristics.H2.Gas.zeta(),
                phi(each stateSelect=StateSelect.always, each fixed=true)),
              H2O(
                initEnergy=Init.none,
                upstreamX=false,
                zeta=100*Characteristics.H2O.Gas.zeta())),
            graphite(
              k=fill((0.5*(1 - epsilon))^(-0.5), 3),
              'inclC+'=true,
              'incle-'=true,
              'incle-Transfer'=true,
              'C+'(theta=U.m*U.K/(1.18*U.W)*(0.5*(1 - epsilon))^1.5, epsilon=
                    0.5*(1 - epsilon)),
              'e-'(sigma=40*U.S/(12*U.cm)/(0.5*(1 - epsilon))^1.5)),
            ionomer(
              k=fill((0.5*(1 - epsilon))^(-0.5), 3),
              'inclSO3-'=true,
              'inclH+'=true,
              inclH2O=true,
              'H+'(initEnergy=Init.none, sigma=40*U.S/U.m),
              'SO3-'(epsilon=0.5*(1 - epsilon),T(each fixed=false, each
                    stateSelect=StateSelect.default)),
              H2O(initEnergy=Init.none, phi(each stateSelect=StateSelect.always,
                    each fixed=true))),
            liquid(
              k=fill(epsilon^(-0.5), 3),
              inclH2O=true,
              H2O(
                upstreamX=false,
                epsilon_IC=1e-5,
                phi(each stateSelect=StateSelect.always, each fixed=true),
                each initEnergy=Init.none,
                mu=100*Characteristics.H2O.Liquid.mu(),
                N0=0.1*U.C,
                T(each stateSelect=StateSelect.default))),
            volume(kappa=6.46e-6*U.mm^2)),
        subregions(graphite('e-Transfer'(final I0=J0*subregions.A[Axis.x]))))
        annotation (IconMap(primitivesVisible=false));

      // See the documentation layer of Phases.PartialPhase regarding the
      // settings of k for each phase.

      parameter Q.NumberAbsolute epsilon(nominal=1) = 0.4 "Porosity"
        annotation (Dialog(group="Geometry", __Dymola_label=
              "<html>&epsilon;</html>"));

      parameter Q.CurrentAreic J0(min=0) = U.A/U.cm^2
        "Exchange current density @ 300 K"
        annotation (Dialog(__Dymola_label="<html><i>J</i><sup>o</sup></html>"));

    protected
      outer Conditions.Environment environment "Environmental conditions";

      annotation (
        defaultComponentPrefixes="replaceable",
        Documentation(info="<html>
<p>This model represents the anode catalyst layer of a PEMFC.
The x axis extends from the anode to the cathode.

By default, the cross-sectional area in the yz plane is 50 cm<sup>2</sup>.</p>

<p>The default thickness (<i>L</i><sub>x</sub> = <code>{28.7*U.um}</code>) is from [<a href=\"modelica://FCSys.UsersGuide.References.Gurau1998\">Gurau1998</a>].
The default thermal conductivity of the carbon (&theta; = <code>U.m*U.K/(1.18*U.W)</code>)
   represents a compressed SGL Sigracet 10 BA gas diffusion layer
  [<a href=\"modelica://FCSys.UsersGuide.References.Nitta2008\">Nitta2008</a>].  The default
  electronic conductivity (&sigma; = <code>40*U.S/(12*U.cm)</code>)
  is for SGL Carbon Group Sigracet&reg; 10 BA (see <a href=\"modelica://FCSys.Regions.AnGDLs.Sigracet10BA\">AnGDLs.Sigracet10BA</a>).</p>

<p>Default assumptions (may be adjusted):
<ol>
<li>All of the species have the same temperature, even in different phases.</li> 
<li>Half of the solid is graphite and half is ionomer (by volume).</li>
</ol></p>

<p>For more information, please see the
    <a href=\"modelica://FCSys.Regions.Region\">Region</a> model.</p>
</html>"),
        Icon(coordinateSystem(
            preserveAspectRatio=true,
            extent={{-100,-100},{100,100}},
            initialScale=0.1), graphics={Rectangle(
                  extent={{-98,62},{98,98}},
                  fillColor={255,255,255},
                  visible=not inclTransY,
                  fillPattern=FillPattern.Solid,
                  pattern=LinePattern.None),Rectangle(
                  extent={{-22.6085,-62.7355},{-56.551,-79.7156}},
                  lineColor={64,64,64},
                  rotation=45,
                  fillColor={127,127,127},
                  fillPattern=FillPattern.HorizontalCylinder,
                  origin={-34.3741,132.347}),Rectangle(
                  extent={{-105.39,79.1846},{-139.33,70.6991}},
                  lineColor={0,0,0},
                  fillColor={200,200,200},
                  rotation=45,
                  fillPattern=FillPattern.HorizontalCylinder,
                  origin={148.513,82.5291}),Polygon(
                  points={{-14,40},{6,60},{14,60},{-6,40},{-14,40}},
                  fillPattern=FillPattern.HorizontalCylinder,
                  smooth=Smooth.None,
                  fillColor={0,0,0},
                  pattern=LinePattern.None),Rectangle(
                  extent={{-6,40},{6,-60}},
                  lineColor={0,0,0},
                  fillColor={200,200,200},
                  fillPattern=FillPattern.VerticalCylinder),Rectangle(
                  extent={{-38,40},{-14,-60}},
                  lineColor={64,64,64},
                  fillColor={127,127,127},
                  fillPattern=FillPattern.VerticalCylinder),Rectangle(
                  extent={{-14,40},{-6,-60}},
                  fillPattern=FillPattern.Solid,
                  fillColor={0,0,0},
                  pattern=LinePattern.None),Polygon(
                  points={{-20,0},{-20,40},{0,60},{20,60},{20,0},{42,0},{42,80},
                {-42,80},{-42,0},{-20,0}},
                  fillPattern=FillPattern.HorizontalCylinder,
                  smooth=Smooth.None,
                  fillColor={255,255,255},
                  pattern=LinePattern.None),Polygon(
                  points={{-20,0},{-20,-60},{0,-60},{20,-40},{20,0},{42,0},{42,
                -80},{-42,-80},{-42,0},{-20,0}},
                  fillPattern=FillPattern.HorizontalCylinder,
                  smooth=Smooth.None,
                  fillColor={255,255,255},
                  pattern=LinePattern.None),Polygon(points={{0,60},{20,60},{0,
              40},{-20,40},{0,60}}, lineColor={0,0,0}),Rectangle(
                  extent={{-20,40},{0,-60}},
                  pattern=LinePattern.None,
                  lineColor={0,0,0}),Polygon(
                  points={{20,60},{0,40},{0,-60},{20,-40},{20,60}},
                  lineColor={0,0,0},
                  fillColor={200,200,200},
                  fillPattern=FillPattern.Solid),Line(
                  points={{-20,0},{-100,0}},
                  color={240,0,0},
                  visible=inclTransX,
                  thickness=0.5),Line(
                  points={{10,0},{100,0}},
                  color={240,0,0},
                  visible=inclTransX,
                  thickness=0.5),Line(
                  points={{0,-60},{0,-100}},
                  color={253,52,56},
                  visible=inclTransY,
                  smooth=Smooth.None,
                  thickness=0.5),Line(
                  points={{0,100},{0,50}},
                  color={253,52,56},
                  visible=inclTransY,
                  smooth=Smooth.None,
                  thickness=0.5),Line(
                  points={{-50,-50},{-10,-10}},
                  color={253,52,56},
                  visible=inclTransZ,
                  smooth=Smooth.None,
                  thickness=0.5),Line(
                  points={{20,20},{50,50}},
                  color={253,52,56},
                  visible=inclTransZ,
                  smooth=Smooth.None,
                  thickness=0.5),Text(
                  extent={{-100,60},{100,100}},
                  textString="%name",
                  visible=not inclTransY,
                  lineColor={0,0,0})}));

    end AnCL;

    model AnCGDL "Integrated anode catalyst/gas diffusion layer"
      extends AnCLs.AnCL(
        L_x={(28.7*U.um + 0.3*U.mm)},
        epsilon=0.5,
        subregions(
          common(each k_Phi={0.3,0.3,0.3}),
          gas(H2(each final eta=0),H2O(each final eta=0)),
          ionomer('H+'(each sigma=650*U.S/U.m)),
          liquid(H2O(each final eta=0))));

      annotation (Documentation(info="<html><p>The default thickness is the total thickness of
  <a href=\"modelica://FCSys.Regions.AnCLs.AnCL\">AnCL</a> and
  <a href=\"modelica://FCSys.Regions.AnGDLs.AnGDL\">AnGDL</a>.</p>

  <p>For more information, please see the
    <a href=\"modelica://FCSys.Regions.AnCLs.AnCL\">AnCL</a> model.</p></html>"));

    end AnCGDL;

  end AnCLs;

  package PEMs "Proton exchange membranes"
    extends Modelica.Icons.Package;
    model PEM "Proton exchange membrane"
      // extends FCSys.Icons.Names.Top4;

      // Note:  Extensions of PEM should be placed directly in the PEMs
      // package rather than subpackages (e.g., by manufacturer) so that
      // __Dymola_choicesFromPackage can be used.  Dymola 7.4 launches the
      // parameter dialogs too slowly when __Dymola_choicesAllMatching is
      // used.

      extends Region(
        L_x={100}*U.um,
        L_y={8}*U.cm,
        L_z={6.25}*U.cm,
        final inclTransX=true,
        inclTransY=false,
        inclTransZ=false,
        redeclare replaceable model Subregion =
            FCSys.Subregions.SubregionIonomer (ionomer(
              'inclSO3-'=true,
              'inclH+'=true,
              inclH2O=true,
              'SO3-'(final mu=0,final epsilon=1),
              'H+'(initEnergy=Init.none, sigma=40*U.S/U.m),
              H2O(initEnergy=Init.none,upstreamX=false))),
        subregions(ionomer('H+'(consTransX={{{if x > 1 or (y == 1 and z == 1)
                   then ConsTrans.steady else ConsTrans.dynamic for z in 1:n_z}
                  for y in 1:n_y} for x in 1:n_x})))) annotation (IconMap(
            primitivesVisible=false));

    protected
      Q.Current I_states_1[:](
        each stateSelect=StateSelect.always,
        each start=0,
        each fixed=true) = subregions[1, 2:n_y, 1].ionomer.'H+'.I[1] if n_y > 1
        "Forced states, set #1";
      Q.Current I_states_2[:, :](
        each stateSelect=StateSelect.always,
        each start=0,
        each fixed=true) = subregions[1, :, 2:n_z].ionomer.'H+'.I[1] if n_z > 1
        "Forced states, set #2";
      // Note:  These variables avoid dynamic state selection in Dymola 2014.

      outer Conditions.Environment environment "Environmental conditions";

      annotation (
        defaultComponentName="PEM",
        Documentation(info="<html>
<p>This model represents the proton exchange membrane of a PEMFC.
The x axis extends from the anode to the cathode.
By default, the cross-sectional area in the yz plane is 50 cm<sup>2</sup>.</p>

<p>The transport of protons includes inertance or inductance (i.e., translational momentum is stored) in the

x direction in the subregions with
index (<i>x</i> = 1, <i>y</i> > 1, <i>z</i> > 1).  This is for numerical reasons,
although in reality there is inductance.</p>

<p>Assumptions:<ol>
<li>There are no pores in the PEM.  All H<sub>2</sub>O is absorbed into the ionomer itself.</li>
<li>There is no cross-over of the reactant gases.</li>
</ul></p>

<p>For more information, please see the
    <a href=\"modelica://FCSys.Regions.Region\">Region</a> model.</p>
</html>"),
        Icon(coordinateSystem(
            preserveAspectRatio=true,
            extent={{-100,-100},{100,100}},
            initialScale=0.1), graphics={Rectangle(
                  extent={{-98,62},{98,98}},
                  fillColor={255,255,255},
                  visible=not inclTransY,
                  fillPattern=FillPattern.Solid,
                  pattern=LinePattern.None),Rectangle(
                  extent={{-99.092,-21.1179},{-84.9489,-63.5448}},
                  lineColor={200,200,200},
                  fillColor={255,255,255},
                  rotation=-45,
                  fillPattern=FillPattern.VerticalCylinder,
                  origin={95.001,14.864}),Rectangle(
                  extent={{-20,40},{0,-60}},
                  lineColor={200,200,200},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.VerticalCylinder),Polygon(
                  points={{20,0},{42,0},{42,80},{-42,80},{-42,0},{-20,0},{-20,
                40},{0,60},{20,60},{20,0}},
                  smooth=Smooth.None,
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid,
                  pattern=LinePattern.None),Polygon(
                  points={{20,0},{42,0},{42,-80},{-42,-80},{-42,0},{-20,0},{-20,
                -60},{0,-60},{20,-40},{20,0}},
                  smooth=Smooth.None,
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid,
                  pattern=LinePattern.None),Polygon(
                  points={{0,40},{20,60},{20,-40},{0,-60},{0,40}},
                  lineColor={0,0,0},
                  smooth=Smooth.None,
                  fillPattern=FillPattern.Solid,
                  fillColor={200,200,200}),Rectangle(extent={{-20,40},{0,-60}},
              lineColor={0,0,0}),Polygon(
                  points={{0,60},{20,60},{0,40},{-20,40},{0,60}},
                  lineColor={0,0,0},
                  smooth=Smooth.None),Line(
                  points={{-20,0},{-100,0}},
                  color={240,0,0},
                  visible=inclTransX,
                  thickness=0.5),Line(
                  points={{10,0},{100,0}},
                  color={0,0,240},
                  visible=inclTransX,
                  thickness=0.5),Line(
                  points={{0,-60},{0,-100}},
                  color={127,127,127},
                  visible=inclTransY,
                  smooth=Smooth.None,
                  thickness=0.5),Line(
                  points={{0,100},{0,50}},
                  color={127,127,127},
                  visible=inclTransY,
                  smooth=Smooth.None,
                  thickness=0.5),Line(
                  points={{-50,-50},{-10,-10}},
                  color={127,127,127},
                  visible=inclTransZ,
                  smooth=Smooth.None,
                  thickness=0.5),Line(
                  points={{20,20},{50,50}},
                  color={127,127,127},
                  visible=inclTransZ,
                  smooth=Smooth.None,
                  thickness=0.5),Text(
                  extent={{-100,60},{100,100}},
                  textString="%name",
                  visible=not inclTransY,
                  lineColor={0,0,0})}));

    end PEM;

    model DuPontN112 "<html>DuPont<sup>TM</sup> Nafion&reg; N-112</html>"
      extends PEM(L_x={51*U.um},Subregion(ionomer('H+'(sigma=0.083*U.S/U.cm))));
      // Additional properties not yet incorporated [DuPont2004N]:
      //     Density:  (100 g/m2)/(51 um) = 1.9608 g/cm3
      annotation (defaultComponentName="PEM", Documentation(info="<html>

    <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References.DuPont2004N\">DuPont2004N</a>].</p>

  <p>For more information, please see the
          <a href=\"modelica://FCSys.Regions.PEMs.PEM\">PEM</a> model.</p></html>"));

    end DuPontN112;

    model DuPontN115 "<html>DuPont<sup>TM</sup> Nafion&reg; N-115</html>"
      extends PEM(L_x={127*U.um},Subregion(ionomer('H+'(sigma=0.083*U.S/U.cm))));
      // Additional properties not yet incorporated [DuPont2004N] and [DuPont2005]:
      //     Density:  (250 g/m2)/(127 um) = 1.9685 g/cm3
      annotation (defaultComponentName="PEM", Documentation(info="<html>
    <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References\">DuPont2004N</a>].</p>

  <p>For more information, please see the
          <a href=\"modelica://FCSys.Regions.PEMs.PEM\">PEM</a> model.</p></html>"));

    end DuPontN115;

    model DuPontN117 "<html>DuPont<sup>TM</sup> Nafion&reg; N-117</html>"
      extends PEM(L_x={183*U.um},Subregion(ionomer('H+'(sigma=0.083*U.S/U.cm))));
      // Additional properties not yet incorporated [DuPont2004N] and [DuPont2005]:
      //     Density:  (360 g/m2)/(183 um) = 1.9672 g/cm3
      annotation (defaultComponentName="PEM", Documentation(info="<html>
    <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References.DuPont2004N\">DuPont2004N</a>].</p>

  <p>For more information, please see the
          <a href=\"modelica://FCSys.Regions.PEMs.PEM\">PEM</a> model.</p></html>"));

    end DuPontN117;

    model DuPontNE1110 "<html>DuPont<sup>TM</sup> Nafion&reg; NE-1110</html>"
      extends PEM(L_x={254*U.um},Subregion(ionomer('H+'(sigma=0.083*U.S/U.cm))));
      // Additional properties not yet incorporated [DuPont2004N] and [DuPont2005]:
      //     Density:  (500 g/m2)/(254 um) = 1.9685 g/cm3
      annotation (defaultComponentName="PEM", Documentation(info="<html>
    <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References.DuPont2004N\">DuPont2004N</a>].</p>

  <p>For more information, please see the
          <a href=\"modelica://FCSys.Regions.PEMs.PEM\">PEM</a> model.</p></html>"));

    end DuPontNE1110;

    model DuPontNE1135 "<html>DuPont<sup>TM</sup> Nafion&reg; NE-1135</html>"
      extends PEM(L_x={89*U.um},Subregion(ionomer('H+'(sigma=0.083*U.S/U.cm))));
      // Additional properties not yet incorporated [DuPont2004N] and [DuPont2005]:
      //     Density:  (190 g/m2)/(89 um) = 2.1348 g/cm3
      annotation (defaultComponentName="PEM", Documentation(info="<html>
    <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References.DuPont2004N\">DuPont2004N</a>].</p>

  <p>For more information, please see the
          <a href=\"modelica://FCSys.Regions.PEMs.PEM\">PEM</a> model.</p></html>"));

    end DuPontNE1135;

    model DuPontNRE211 "<html>DuPont<sup>TM</sup> Nafion&reg; NRE-1110</html>"
      extends PEM(L_x={25.4*U.um},Subregion(ionomer('H+'(sigma=0.083*U.S/U.cm))));
      // Additional properties not yet incorporated [DuPont2004NRE]:
      //     Density:  (50 g/m2)/(25.4 um) = 1.9685 g/cm3
      annotation (defaultComponentName="PEM", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References.DuPont2004N\">DuPont2004N</a>], except that

  the default value of protonic conductivity (&sigma; = <code>0.083*U.S/U.cm</code>)
  is for the

  DuPont<sup>TM</sup> Nafion&reg; N and NE series
  [<a href=\"modelica://FCSys.UsersGuide.References.DuPont2004N\">DuPont2004N</a>].

  It is not listed for DuPont<sup>TM</sup> Nafion&reg; NRE-211
  in [<a href=\"modelica://FCSys.UsersGuide.References.DuPont2004NRE\">DuPont2004NRE</a>].</p>

  <p>For more information, please see the
          <a href=\"modelica://FCSys.Regions.PEMs.PEM\">PEM</a> model.</p></html>"));

    end DuPontNRE211;

    model DuPontNRE212 "<html>DuPont<sup>TM</sup> Nafion&reg; NRE-1110</html>"
      extends PEM(L_x={50.8*U.um},Subregion(ionomer('H+'(sigma=0.083*U.S/U.cm))));
      // Additional properties not yet incorporated [DuPont2004NRE]:
      //     Density:  (100 g/m2)/(50.8 um) = 1.9685 g/cm3
      annotation (defaultComponentName="PEM", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References.DuPont2004N\">DuPont2004N</a>], except that

  the default value of protonic conductivity (&sigma; = <code>0.083*U.S/U.cm</code>)
  is for the
  DuPont<sup>TM</sup> Nafion&reg; N and NE series
  [<a href=\"modelica://FCSys.UsersGuide.References.DuPont2004N\">DuPont2004N</a>].
    It is not listed for DuPont<sup>TM</sup> Nafion&reg; NRE-212
  in [<a href=\"modelica://FCSys.UsersGuide.References.DuPont2004NRE\">DuPont2004NRE</a>].</p></p>

  <p>For more information, please see the
          <a href=\"modelica://FCSys.Regions.PEMs.PEM\">PEM</a> model.</p></html>"));

    end DuPontNRE212;

  end PEMs;

  package CaCLs "Cathode catalyst layers"
    extends Modelica.Icons.Package;

    model CaCL "Cathode catalyst layer"
      import Modelica.Constants.inf;
      // extends FCSys.Icons.Names.Top4;

      extends Region(
        L_x={28.7}*U.um,
        L_y={8}*U.cm,
        L_z={6.25}*U.cm,
        final inclTransX=true,
        inclTransY=false,
        inclTransZ=false,
        redeclare replaceable model Subregion = Subregions.Subregion (
            common(k_Q=0),
            gasLiq(k_Phi={2,2,2}, k_Q=inf),
            gas(
              common(k_Q=inf),
              k=fill(epsilon^(-0.5), 3),
              inclH2O=true,
              inclN2=true,
              inclO2=true,
              H2O(
                initEnergy=Init.none,
                upstreamX=false,
                phi(each stateSelect=StateSelect.always, each fixed=true),
                zeta=100*Characteristics.H2O.Gas.zeta()),
              N2(
                initEnergy=Init.none,
                upstreamX=false,
                zeta=100*Characteristics.N2.Gas.zeta()),
              O2(
                initEnergy=Init.none,
                upstreamX=false,
                I(each stateSelect=StateSelect.always, each fixed=true),
                zeta=100*Characteristics.O2.Gas.zeta(),
                p_stop=12*U.Pa)),
            graphite(
              k=fill((0.5*(1 - epsilon))^(-0.5), 3),
              'inclC+'=true,
              'incle-'=true,
              'incle-Transfer'=true,
              'C+'(theta=U.m*U.K/(1.18*U.W)*(0.5*(1 - epsilon))^1.5, epsilon=
                    0.5*(1 - epsilon)),
              'e-'(sigma=40*U.S/(12*U.cm)/(0.5*(1 - epsilon))^1.5),
              'e-Transfer'(E_A=0.46*U.V)),
            ionomer(
              k=fill((0.5*(1 - epsilon))^(-0.5), 3),
              'inclSO3-'=true,
              'inclH+'=true,
              inclH2O=true,
              'SO3-'(epsilon=0.5*(1 - epsilon), T(each fixed=false, each
                    stateSelect=StateSelect.default)),
              'H+'(initEnergy=Init.none, sigma=40*U.S/U.m),
              H2O(initEnergy=Init.none, phi(each stateSelect=StateSelect.always,
                    each fixed=true))),
            liquid(
              k=fill(epsilon^(-0.5), 3),
              inclH2O=true,
              H2O(
                upstreamX=false,
                epsilon_IC=1e-5,
                phi(each stateSelect=StateSelect.always, each fixed=true),
                each initEnergy=Init.none,
                mu=100*Characteristics.H2O.Liquid.mu(),
                N0=0.1*U.C,
                T(each stateSelect=StateSelect.default))),
            volume(kappa=6.46e-6*U.mm^2)),
        subregions(graphite('e-Transfer'(final I0=J0*subregions.A[Axis.x]))))
        annotation (IconMap(primitivesVisible=false));

      // See the documentation layer of Phases.PartialPhase regarding the
      // settings of k for each phase.

      parameter Q.NumberAbsolute epsilon(nominal=1) = 0.4 "Porosity"
        annotation (Dialog(group="Geometry", __Dymola_label=
              "<html>&epsilon;</html>"));

      parameter Q.CurrentAreic J0(min=0) = 0.11*U.mA/U.cm^2
        "Exchange current density @ 300 K"
        annotation (Dialog(__Dymola_label="<html><i>J</i><sup>o</sup></html>"));

      // Auxiliary variables (for analysis)
      output Q.Potential Deltaw[n_y, n_z]=-subregions[n_x, :, :].graphite.'e-'.g_boundaries[
          1, Side.p] - subregions[1, :, :].ionomer.'H+'.g_boundaries[1, Side.n]
        if environment.analysis
        "Electrical potential differences over the yz plane";

    protected
      outer Conditions.Environment environment "Environmental conditions";

    initial equation

      annotation (Documentation(info="<html>
<p>This model represents the cathode catalyst layer of a PEMFC.
The x axis extends from the anode to the cathode.

By default, the cross-sectional area in the yz plane is 50 cm<sup>2</sup>.</p>

<p>The default thickness (<i>L</i><sub>x</sub> = <code>{28.7*U.um}</code>) is from [<a href=\"modelica://FCSys.UsersGuide.References.Gurau1998\">Gurau1998</a>].
The default thermal conductivity of the carbon (&theta; = <code>U.m*U.K/(1.18*U.W)</code>)
   represents a compressed SGL Sigracet 10 BA gas diffusion layer
  [<a href=\"modelica://FCSys.UsersGuide.References.Nitta2008\">Nitta2008</a>].  The default
  electronic mobility (&sigma; = <code>40*U.S/(12*U.cm)</code>)
  is for SGL Carbon Group Sigracet&reg; 10 BA (see <a href=\"modelica://FCSys.Regions.CaGDLs.Sigracet10BA\">CAGDLs.Sigracet10BA</a>).  The

  default activation energy for the oxygen reduction reaction (<code>E_A = 73.2e3*U.J/U.mol</code>) is from

  [<a href=\"modelica://FCSys.UsersGuide.References.Sivertsen2005\">Sivertsen2005</a>].</p>

<p>Default assumptions (may be adjusted):
<ol>
<li>All of the species have the same temperature, even in different phases.</li> 
<li>Half of the solid is graphite and half is ionomer (by volume).</li>
</ol></p>

<p>For more information, please see the
    <a href=\"modelica://FCSys.Regions.Region\">Region</a> model.</p>
</html>"), Icon(coordinateSystem(
            preserveAspectRatio=true,
            extent={{-100,-100},{100,100}},
            initialScale=0.1), graphics={Rectangle(
                  extent={{-98,62},{98,98}},
                  fillColor={255,255,255},
                  visible=not inclTransY,
                  fillPattern=FillPattern.Solid,
                  pattern=LinePattern.None),Rectangle(
                  extent={{-21.6329,-68.4511},{-58.4038,-85.4311}},
                  lineColor={64,64,64},
                  rotation=45,
                  fillColor={127,127,127},
                  fillPattern=FillPattern.HorizontalCylinder,
                  origin={-15.1055,127.699}),Rectangle(
                  extent={{-105.385,79.1805},{-139.323,70.6948}},
                  lineColor={0,0,0},
                  fillColor={200,200,200},
                  rotation=45,
                  fillPattern=FillPattern.HorizontalCylinder,
                  origin={130.507,84.5292}),Polygon(
                  points={{-14,40},{6,60},{14,60},{-6,40},{-14,40}},
                  fillPattern=FillPattern.HorizontalCylinder,
                  smooth=Smooth.None,
                  fillColor={0,0,0},
                  pattern=LinePattern.None),Rectangle(
                  extent={{-26,40},{-14,-60}},
                  lineColor={0,0,0},
                  fillColor={200,200,200},
                  fillPattern=FillPattern.VerticalCylinder),Rectangle(
                  extent={{-6,40},{18,-60}},
                  lineColor={64,64,64},
                  fillColor={127,127,127},
                  fillPattern=FillPattern.VerticalCylinder),Rectangle(
                  extent={{-14,40},{-6,-60}},
                  fillPattern=FillPattern.Solid,
                  fillColor={0,0,0},
                  pattern=LinePattern.None),Polygon(
                  points={{-20,0},{-20,40},{0,60},{20,60},{20,0},{42,0},{42,80},
                {-42,80},{-42,0},{-20,0}},
                  fillPattern=FillPattern.HorizontalCylinder,
                  smooth=Smooth.None,
                  fillColor={255,255,255},
                  pattern=LinePattern.None),Polygon(
                  points={{-20,0},{-20,-60},{0,-60},{20,-40},{20,0},{42,0},{42,
                -80},{-42,-80},{-42,0},{-20,0}},
                  fillPattern=FillPattern.HorizontalCylinder,
                  smooth=Smooth.None,
                  fillColor={255,255,255},
                  pattern=LinePattern.None),Polygon(points={{0,60},{20,60},{0,
              40},{-20,40},{0,60}}, lineColor={0,0,0}),Rectangle(
                  extent={{-20,40},{0,-60}},
                  pattern=LinePattern.None,
                  lineColor={0,0,0}),Polygon(
                  points={{20,60},{0,40},{0,-60},{20,-40},{20,60}},
                  lineColor={0,0,0},
                  fillColor={120,120,120},
                  fillPattern=FillPattern.Solid),Line(
                  points={{-20,0},{-100,0}},
                  color={0,0,240},
                  visible=inclTransX,
                  thickness=0.5),Line(
                  points={{10,0},{100,0}},
                  color={0,0,240},
                  visible=inclTransX,
                  thickness=0.5),Line(
                  points={{0,-60},{0,-98}},
                  color={0,0,240},
                  visible=inclTransY,
                  smooth=Smooth.None,
                  thickness=0.5),Line(
                  points={{0,100},{0,50}},
                  color={0,0,240},
                  visible=inclTransY,
                  smooth=Smooth.None,
                  thickness=0.5),Line(
                  points={{-50,-50},{-10,-10}},
                  color={0,0,240},
                  visible=inclTransZ,
                  smooth=Smooth.None,
                  thickness=0.5),Line(
                  points={{20,20},{50,50}},
                  color={0,0,240},
                  visible=inclTransZ,
                  smooth=Smooth.None,
                  thickness=0.5),Text(
                  extent={{-100,60},{100,100}},
                  textString="%name",
                  visible=not inclTransY,
                  lineColor={0,0,0})}));
    end CaCL;

    model CaCGDL "Integrated cathode catalyst/gas diffusion layer"
      extends CaCLs.CaCL(
        L_x={(28.7*U.um + 0.3*U.mm)},
        epsilon=0.5,
        subregions(
          common(each k_Phi={0.3,0.3,0.3}),
          gas(H2(each final eta=0),H2O(each final eta=0)),
          ionomer('H+'(each sigma=650*U.S/U.m)),
          liquid(H2O(each final eta=0))));
      annotation (Documentation(info="<html><p>The default thickness is the total thickness of
  <a href=\"modelica://FCSys.Regions.CaCLs.CaCL\">CaCL</a> and
  <a href=\"modelica://FCSys.Regions.CaGDLs.CaGDL\">CaGDL</a>.</p>

  <p>For more information, please see the
    <a href=\"modelica://FCSys.Regions.CaCLs.CaCL\">CaCL</a> model.</p></html>"));

    end CaCGDL;

  end CaCLs;

  package CaGDLs "Cathode gas diffusion layers"
    model CaGDL "Cathode gas diffusion layer"
      import Modelica.Constants.inf;
      // extends FCSys.Icons.Names.Top4;

      // Note:  Extensions of CaGDL should be placed directly in the CaGDLs
      // package rather than subpackages (e.g., by manufacturer) so that
      // __Dymola_choicesFromPackage can be used.  Dymola 7.4 launches the
      // parameter dialogs too slowly when __Dymola_choicesAllMatching is
      // used.

      extends Region(
        L_x={0.3}*U.mm,
        L_y={8}*U.cm,
        L_z={6.25}*U.cm,
        final inclTransX=true,
        inclTransY=false,
        inclTransZ=false,
        redeclare replaceable model Subregion =
            FCSys.Subregions.SubregionNoIonomer (
            common(k_Phi={200,200,200},k_Q=0),
            gasLiq(k_Phi={2,2,2}, k_Q=inf),
            gas(
              common(k_Phi={0.6,0.6,0.6},k_Q=inf),
              k=fill(epsilon^(-0.5), 3),
              inclH2O=true,
              inclN2=true,
              inclO2=true,
              H2O(
                initEnergy=Init.none,
                upstreamX=false,
                phi(each stateSelect=StateSelect.always, each fixed=true),
                zeta=100*Characteristics.H2O.Gas.zeta(),
                final eta=0),
              N2(
                initEnergy=Init.none,
                upstreamX=false,
                phi(each stateSelect=StateSelect.always, each fixed=true),
                zeta=100*Characteristics.N2.Gas.zeta(),
                final eta=0),
              O2(
                initEnergy=Init.none,
                upstreamX=false,
                I(each stateSelect=StateSelect.always, each fixed=true),
                zeta=100*Characteristics.O2.Gas.zeta(),
                final eta=0)),
            graphite(
              k=fill((1 - epsilon)^(-0.5), 3),
              'inclC+'=true,
              'incle-'=true,
              'C+'(theta=U.m*U.K/(1.18*U.W)*(1 - epsilon)^1.5,final epsilon=1
                     - epsilon),
              'e-'(sigma=40*U.S/(12*U.cm)/(1 - epsilon)^1.5)),
            liquid(
              k=fill(epsilon^(-0.5), 3),
              inclH2O=true,
              H2O(
                upstreamX=false,
                epsilon_IC=1e-5,
                final eta=0,
                phi(each stateSelect=StateSelect.always, each fixed=true),
                each initEnergy=Init.none,
                mu=100*Characteristics.H2O.Liquid.mu(),
                T(each stateSelect=StateSelect.default))))) annotation (IconMap(
            primitivesVisible=false));

      // Note:  The fluid species have zero fluidity (eta=0) so that the transverse
      // velocity is zero at the interface with the flow plate.  That condition
      // is necessary to produce the appropriate pressure loss down the channel.

      // See the documentation layer of Phases.PartialPhase regarding the
      // settings of k for each phase.

      parameter Q.NumberAbsolute epsilon(nominal=1) = 0.88 "Porosity"
        annotation (Dialog(group="Geometry", __Dymola_label=
              "<html>&epsilon;</html>"));

    protected
      outer Conditions.Environment environment "Environmental conditions";
      annotation (Documentation(info="<html>
<p>This model represents the cathode gas diffusion layer of a PEMFC.
The x axis extends from the anode to the cathode.

By default, the cross-sectional area in the yz plane is 50 cm<sup>2</sup>.</p>

<p>The default porosity (&epsilon; = 0.88) is that of SGL Carbon Group Sigracet&reg; 10 BA and 25 BA GDLs.

The porosity of a GDL may be lower than specified due to compression (e.g., 0.4 according to

[<a href=\"modelica://FCSys.UsersGuide.References.Bernardi1992\">Bernardi1992</a>, p. 2483], although
that reference may be outdated).
  The default thermal conductivity of the carbon (&theta; = <code>U.m*U.K/(1.18*U.W)</code>)
   represents a compressed Sigracet&reg; 10 BA gas diffusion layer
  [<a href=\"modelica://FCSys.UsersGuide.References.Nitta2008\">Nitta2008</a>].  The default
  electrical conductivity
  is also for Sigracet&reg; 10 BA [<a href=\"modelica://FCSys.UsersGuide.References.SGL2007\">SGL2007</a>].</p>

<p>Default assumptions (may be adjusted):
<ol>
<li>All of the species have the same temperature, even in different phases.</li> 
</ol></p>

<p>For more information, please see the
    <a href=\"modelica://FCSys.Regions.Region\">Region</a> model.</p>
</html>"), Icon(coordinateSystem(
            preserveAspectRatio=true,
            extent={{-100,-100},{100,100}},
            initialScale=0.1), graphics={Rectangle(
                  extent={{-98,62},{98,98}},
                  fillColor={255,255,255},
                  visible=not inclTransY,
                  fillPattern=FillPattern.Solid,
                  pattern=LinePattern.None),Rectangle(
                  extent={{-78.7855,18.6813},{-50.5004,-23.7455}},
                  lineColor={64,64,64},
                  fillColor={127,127,127},
                  rotation=-45,
                  fillPattern=FillPattern.VerticalCylinder,
                  origin={52.5001,1.0805}),Rectangle(
                  extent={{-20,40},{20,-60}},
                  lineColor={64,64,64},
                  fillColor={127,127,127},
                  fillPattern=FillPattern.VerticalCylinder),Polygon(
                  points={{20,0},{42,0},{42,80},{-42,80},{-42,0},{-20,0},{-20,
                40},{0,60},{20,60},{20,0}},
                  smooth=Smooth.None,
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid,
                  pattern=LinePattern.None),Polygon(
                  points={{20,0},{42,0},{42,-80},{-42,-80},{-42,0},{-20,0},{-20,
                -60},{0,-60},{20,-40},{20,0}},
                  smooth=Smooth.None,
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid,
                  pattern=LinePattern.None),Polygon(
                  points={{0,40},{20,60},{20,-40},{0,-60},{0,40}},
                  lineColor={0,0,0},
                  smooth=Smooth.None,
                  fillPattern=FillPattern.Solid,
                  fillColor={127,127,127}),Rectangle(extent={{-20,40},{0,-60}},
              lineColor={0,0,0}),Polygon(
                  points={{0,60},{20,60},{0,40},{-20,40},{0,60}},
                  lineColor={0,0,0},
                  smooth=Smooth.None),Line(
                  points={{-20,0},{-100,0}},
                  color={0,0,240},
                  visible=inclTransX,
                  thickness=0.5),Line(
                  points={{10,0},{100,0}},
                  color={0,0,240},
                  thickness=0.5),Line(
                  points={{0,-60},{0,-100}},
                  color={0,0,240},
                  visible=inclTransY,
                  smooth=Smooth.None,
                  thickness=0.5),Line(
                  points={{0,100},{0,50}},
                  color={0,0,240},
                  visible=inclTransY,
                  smooth=Smooth.None,
                  thickness=0.5),Line(
                  points={{-50,-50},{-10,-10}},
                  color={0,0,240},
                  visible=inclTransZ,
                  smooth=Smooth.None,
                  thickness=0.5),Line(
                  points={{20,20},{50,50}},
                  color={0,0,240},
                  visible=inclTransZ,
                  smooth=Smooth.None,
                  thickness=0.5),Text(
                  extent={{-100,60},{100,100}},
                  textString="%name",
                  visible=not inclTransY,
                  lineColor={0,0,0})}));

    end CaGDL;
    extends Modelica.Icons.Package;

    // Note:  Extensions of CaGDL should be placed directly in the CaGDLs
    // package rather than nested packages (e.g., by manufacturer) so that
    // __Dymola_choicesFromPackage can be used.  In Dymola 7.4 the
    // parameter dialogs launch too slowly when __Dymola_choicesAllMatching
    // is used.

    model Sigracet10BA "<html>SGL Carbon Group Sigracet&reg; 10 BA</html>"
      extends CaGDL(
        L_x={0.400}*U.mm,
        epsilon=0.88,
        Subregion(graphite('e-'(sigma=U.S/(0.012*U.cm)))));
      // Additional properties not incorporated [SGL2007]:
      //     Diffusivity:  L = 0.400 mm, p = 0.85 m/s (for air) => D = P*L = 340 mm2/s
      //     Density:  (85 g/m2)/(0.400 mm)/0.88 = 212.5 kg/m3
      annotation (defaultComponentName="caGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References.SGL2007\">SGL2007</a>].</p>

  <p>For more information, please see the
    <a href=\"modelica://FCSys.Regions.CaGDLs.CaGDL\">CaGDL</a> model.</p></html>"));

    end Sigracet10BA;

    model Sigracet10BB "<html>SGL Carbon Group Sigracet&reg; 10 BB</html>"
      extends CaGDL(
        L_x={0.420}*U.mm,
        epsilon=0.84,
        Subregion(graphite('e-'(sigma=42*U.S/(15*U.cm)))));
      // Additional properties not incorporated [SGL2007]:
      //     Diffusivity:  L = 0.420 mm, p = 0.03 m/s (for air) => D = P*L = 12.6 mm2/s
      //     Density:  (125 g/m2)/(0.420 mm) = 297.62 kg/m3
      annotation (defaultComponentName="caGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References.SGL2007\">SGL2007</a>].</p>

  <p>For more information, please see the
    <a href=\"modelica://FCSys.Regions.CaGDLs.CaGDL\">CaGDL</a> model.</p></html>"));

    end Sigracet10BB;

    model Sigracet10BC "<html>SGL Carbon Group Sigracet&reg; 10 BC</html>"
      extends CaGDL(
        L_x={0.420}*U.mm,
        epsilon=0.82,
        Subregion(graphite('e-'(sigma=42*U.S/(16*U.cm)))));
      // Additional properties not incorporated [SGL2007]:
      //     Diffusivity:  L = 0.420 mm, p = 0.0145 m/s (for air) => D = P*L = 6.09 mm2/s
      //     Density:  (135 g/m2)/(0.420 mm) = 321.43 kg/m3
      annotation (defaultComponentName="caGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References.SGL2007\">SGL2007</a>].</p>

  <p>For more information, please see the
    <a href=\"modelica://FCSys.Regions.CaGDLs.CaGDL\">CaGDL</a> model.</p></html>"));

    end Sigracet10BC;

    model Sigracet24BA "<html>SGL Carbon Group Sigracet&reg; 24 BA</html>"
      extends CaGDL(
        L_x={0.190}*U.mm,
        epsilon=0.84,
        Subregion(graphite('e-'(sigma=19*U.S/(10*U.cm)))));
      // Additional properties not incorporated [SGL2004]:
      //     Diffusivity:  L = 0.190 mm, p = 0.30 m/s (for air) => D = P*L = 57 mm2/s
      //     Density:  (54 g/m2)/(0.190 mm) = 284.21 kg/m3
      annotation (defaultComponentName="caGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References.SGL2004\">SGL2004</a>].</p>

  <p>For more information, please see the
    <a href=\"modelica://FCSys.Regions.CaGDLs.CaGDL\">CaGDL</a> model.</p></html>"));

    end Sigracet24BA;

    model Sigracet24BC "<html>SGL Carbon Group Sigracet&reg; 24 BC</html>"
      extends CaGDL(
        L_x={0.235}*U.mm,
        epsilon=0.76,
        Subregion(graphite('e-'(sigma=23.5*U.S/(11*U.cm)))));
      // Additional properties not incorporated [SGL2004]:
      //     Diffusivity:  L = 0.235 mm, p = 0.0045 m/s (for air) => D = P*L = 1.0575 mm2/s
      //     Density:  (100 g/m2)/(0.235 mm) = 425.53 kg/m3
      annotation (defaultComponentName="caGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References.SGL2004\">SGL2004</a>].</p>

  <p>For more information, please see the
    <a href=\"modelica://FCSys.Regions.CaGDLs.CaGDL\">CaGDL</a> model.</p></html>"));

    end Sigracet24BC;

    model Sigracet25BA "<html>SGL Carbon Group Sigracet&reg; 25 BA</html>"
      extends CaGDL(
        L_x={0.190}*U.mm,
        epsilon=0.88,
        Subregion(graphite('e-'(sigma=19*U.S/(10*U.cm)))));
      // Additional properties not incorporated [SGL2004]:
      //     Diffusivity:  L = 0.190 mm, p = 0.90 m/s (for air) => D = P*L = 171 mm2/s
      //     Density:  (40 g/m2)/(0.190 mm) = 210.53 kg/m3
      annotation (defaultComponentName="caGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References.SGL2004\">SGL2004</a>].</p>

  <p>For more information, please see the
    <a href=\"modelica://FCSys.Regions.CaGDLs.CaGDL\">CaGDL</a> model.</p></html>"));

    end Sigracet25BA;

    model Sigracet25BC "<html>SGL Carbon Group Sigracet&reg; 25 BC</html>"
      extends CaGDL(
        L_x={0.235}*U.mm,
        epsilon=0.80,
        Subregion(graphite('e-'(sigma=23.5*U.S/(12*U.cm)))));
      // Additional properties not incorporated [SGL2004]:
      //     Diffusivity:  L = 0.235 mm, p = 0.008 m/s (for air) => D = P*L = 1.88 mm2/s
      //     Density:  (86 g/m2)/(0.235 mm) = 365.96 kg/m3
      annotation (defaultComponentName="caGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References.SGL2004\">SGL2004</a>].</p>

  <p>For more information, please see the
    <a href=\"modelica://FCSys.Regions.CaGDLs.CaGDL\">CaGDL</a> model.</p></html>"));

    end Sigracet25BC;

    model TorayTGPH030 "Toray Industries TGP-H-030"
      extends CaGDL(
        L_x={0.11}*U.mm,
        epsilon=0.80,
        Subregion(graphite('C+'(theta=U.m*U.K/(1.7*U.W)),'e-'(sigma=U.S/(0.80*U.mm)))));
      // Additional properties not incorporated [Toray2010]:
      //     Diffusivity:  L = 0.110 mm, P/p = 2500 ml.mm/(cm2.hr.mmAq) = 0.70814e-3 m/(s.kPa)
      //         => D = P*L = 7.89e-6 m2/s, assuming p = 101.325 kPa
      //     Bulk density:  0.44 g/cm3
      annotation (defaultComponentName="caGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References.Toray2010\">Toray2010</a>].  The electronic mobility (&mu;) is based on the
  through-plane value of resistivity.  The thermal conductivity is not listed but is
  taken to be the same as for TGP-H-060, TGP-H-090, and TGP-H-120.</p>

  <p>For more information, please see the
    <a href=\"modelica://FCSys.Regions.CaGDLs.CaGDL\">CaGDL</a> model.</p></html>"));

    end TorayTGPH030;

    model TorayTGPH060 "Toray Industries TGP-H-060"
      extends CaGDL(
        L_x={0.19}*U.mm,
        epsilon=0.78,
        Subregion(graphite('C+'(theta=U.m*U.K/(1.7*U.W)),'e-'(sigma=U.S/(0.80*U.mm)))));

      // Additional properties not incorporated [Toray2010]:
      //     Diffusivity: L = 0.190 mm, P/p = 1900 ml.mm/(cm2.hr.mmAq) = 0.53818e-3 m/(s.kPa)
      //         => D = P*L = 10.36e-6 m2/s, assuming p = 101.325 kPa
      //     Bulk density:  0.44 g/cm3
      annotation (defaultComponentName="caGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References.Toray2010\">Toray2010</a>].  The electronic mobility (&mu;) is based on the
  through-plane value of resistivity.  The thermal conductivity is through the plane at
  room temperature.</p>

  <p>For more information, please see the
    <a href=\"modelica://FCSys.Regions.CaGDLs.CaGDL\">CaGDL</a> model.</p></html>"));

    end TorayTGPH060;

    model TorayTGPH090 "Toray Industries TGP-H-090"
      extends CaGDL(
        L_x={0.28}*U.mm,
        epsilon=0.78,
        Subregion(graphite('C+'(theta=U.m*U.K/(1.7*U.W)),'e-'(sigma=U.S/(0.80*U.mm)))));

      // Additional properties not incorporated [Toray2010]:
      //     Diffusivity: L = 0.280 mm, P/p = 1700 ml.mm/(cm2.hr.mmAq) = 0.48153e-3 m/(s.kPa)
      //         => D = P*L = 13.66e-6 m2/s, assuming p = 101.325 kPa
      //     Bulk density:  0.44 g/cm3
      annotation (defaultComponentName="caGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References.Toray2010\">Toray2010</a>].  The electronic mobility (&mu;) is based on the
  through-plane value of resistivity.  The thermal conductivity is through the plane at
  room temperature.</p>

  <p>For more information, please see the
    <a href=\"modelica://FCSys.Regions.CaGDLs.CaGDL\">CaGDL</a> model.</p></html>"));

    end TorayTGPH090;

    model TorayTGPH120 "Toray Industries TGP-H-030"
      extends CaGDL(
        L_x={0.37}*U.mm,
        epsilon=0.78,
        Subregion(graphite('C+'(theta=U.m*U.K/(1.7*U.W)),'e-'(sigma=U.S/(0.80*U.mm)))));

      // Additional properties not incorporated [Toray2010]:
      //     Diffusivity: L = 0.370 mm, P/p = 1500 ml.mm/(cm2.hr.mmAq) = 0.42488e-3 m/(s.kPa)
      //         => D = P*L = 15.93e-6 m2/s, assuming p = 101.325 kPa
      //     Bulk density:  0.44 g/cm3
      annotation (defaultComponentName="caGDL", Documentation(info="<html>
  <p>The default properties are based on [<a href=\"modelica://FCSys.UsersGuide.References.Toray2010\">Toray2010</a>].  The electronic mobility (&mu;) is based on the
  through-plane value of resistivity.  The thermal conductivity is through the plane at
  room temperature.</p>

  <p>For more information, please see the
    <a href=\"modelica://FCSys.Regions.CaGDLs.CaGDL\">CaGDL</a> model.</p></html>"));

    end TorayTGPH120;

  end CaGDLs;

  package CaFPs "Cathode flow plates"
    extends Modelica.Icons.Package;

    model CaFP "Cathode flow plate"
      import Modelica.Constants.inf;
      // extends FCSys.Icons.Names.Top4;

      extends Region(
        L_x={8}*U.mm,
        L_y={8}*U.cm,
        L_z={6.25}*U.cm,
        final inclTransX=true,
        final inclTransY=true,
        inclTransZ=false,
        redeclare replaceable model Subregion =
            FCSys.Subregions.SubregionNoIonomer (
            common(k_Phi={1e7,inf,1e7},k_Q=1e5),
            gasLiq(k_Phi={inf,1e6,inf},k_Q=inf),
            gas(
              common(k_Phi={inf,inf,inf}),
              k={epsilon/2,11,1/11},
              inclH2O=true,
              inclN2=true,
              inclO2=true,
              H2O(
                upstreamX=false,
                Nu_Phi={4,16*A[Axis.z]*epsilon/D^2,4},
                zeta=100*Characteristics.H2O.Gas.zeta(),
                T(stateSelect=StateSelect.always)),
              N2(
                upstreamX=false,
                Nu_Phi={4,16*A[Axis.z]*epsilon/D^2,4},
                zeta=100*Characteristics.N2.Gas.zeta(),
                initEnergy=Init.none,
                I(each stateSelect=StateSelect.always, each fixed=true)),
              O2(
                upstreamX=false,
                Nu_Phi={4,16*A[Axis.z]*epsilon/D^2,4},
                zeta=100*Characteristics.O2.Gas.zeta(),
                initEnergy=Init.none)),
            graphite(
              'inclC+'=true,
              'incle-'=true,
              'C+'(theta=U.m*U.K/(95*U.W),epsilon=1 - epsilon),
              'e-'(sigma=U.S/(1.470e-3*U.cm))),
            liquid(
              k={epsilon/2,11,1/11},
              inclH2O=true,
              H2O(
                upstreamX=false,
                Nu_Phi={4,16*A[Axis.z]*epsilon/D^2,4},
                epsilon_IC=1e-5,
                N0=0.1*U.C)),
            volume(inclCapillary=false))) annotation (IconMap(primitivesVisible
            =false));

      parameter Q.NumberAbsolute epsilon(nominal=1) = 0.0625
        "Fraction of volume for the fluid" annotation (Dialog(group="Geometry",
            __Dymola_label="<html>&epsilon;</html>"));

      parameter Q.Length D=1.5*U.mm "Hydraulic diameter of the channel"
        annotation (Dialog(__Dymola_label="<html><i>D</i></html>"));

    protected
      Q.Velocity phi_states_H2O[:, :, :](
        each stateSelect=StateSelect.always,
        each start=0,
        each fixed=true) = subregions[:, 2:n_y, :].gas.H2O.phi[2] if n_y > 1
        "Forced states for H2O";
      Q.Velocity phi_states_O2[:, :, :](
        each stateSelect=StateSelect.always,
        each start=0,
        each fixed=true) = subregions[:, 2:n_y, :].gas.O2.phi[2] if n_y > 1
        "Forced states for O2";
      // Note:  These variables avoid dynamic state selection in Dymola 2014.

      outer Conditions.Environment environment "Environmental conditions";

      // See AnFPs.AnFP for data on additional materials.
      annotation (Documentation(info="<html>
<p>This model represents the cathode flow plate of a PEMFC.
The x axis extends from the anode to the cathode.
Fluid is considered to travel
in the y direction, with the associated length factor (<i>k</i><sub>y</sub>) greater than one (by default)
to represent a serpentine channel.
The model is
bidirectional, meaning that either <code>yNegative</code> or <code>yPositive</code> can be
used as the inlet.  By default, the cross-sectional area in the yz plane is 50 cm<sup>2</sup>.</p>

<p>The solid and the fluid phases are assumed to exist in the same subregions, 
even though a 
typical flow plate is impermeable
to the fluid (except for the channel).  
In theory, it is possible
to discretize the flow plate into smaller subregions for the bulk solid, lands, and valleys.  However,
this would significantly increase the mathematical size of the model.  Currently, that level of detail is best left to
computational fluid dynamics.</p>
 
<p>The x axis-component of the transport factor (<i>k</i><sub>x</sub>) for the gas and the liquid should
generally be less than one because the transport distance into/out of the GDL is less that half the 
thickness of the flow plate. It is equal to the product of two ratios:<ol>
<li>the depth of the channels to the thickness of the flow plate</li>
<li>the product of the total area of the flow plate in the yz plane (land + valleys) and the fraction of
the total volume available for the fluid (&epsilon;) to the area of the valleys in the yz plane</li>
</ol> The default is &epsilon;/2.</p>

<p>See <a href=\"modelica://FCSys.Species.'C+'.Graphite.Fixed\">Species.'C+'.Graphite.Fixed</a>
regarding the default specific heat capacity.  The default thermal resistivity
of the carbon (&theta; = <code>U.m*U.K/(95*U.W)</code>) and the
electrical conductivity (&sigma; = <code>U.S/(1.470e-3*U.cm)</code>)
are that of Entegris/Poco Graphite AXF-5Q
[<a href=\"modelica://FCSys.UsersGuide.References.Entegris2012\">Entegris2012</a>].
  There is additional data in the
text layer of the <a href=\"modelica://FCSys.Regions.AnFPs.AnFP\">AnFP</a> model.</p>

<p>For more information, please see the
    <a href=\"modelica://FCSys.Regions.Region\">Region</a> model.</p>
</html>"), Icon(coordinateSystem(
            preserveAspectRatio=true,
            extent={{-100,-100},{100,100}},
            initialScale=0.1), graphics={Rectangle(
                  extent={{-98,62},{98,98}},
                  fillColor={255,255,255},
                  visible=not inclTransY,
                  fillPattern=FillPattern.Solid,
                  pattern=LinePattern.None),Rectangle(
                  extent={{-76.648,66.211},{-119.073,52.0689}},
                  fillPattern=FillPattern.HorizontalCylinder,
                  rotation=45,
                  fillColor={135,135,135},
                  origin={111.017,77.3801},
                  pattern=LinePattern.None,
                  lineColor={95,95,95}),Rectangle(
                  extent={{-20,40},{0,-60}},
                  lineColor={95,95,95},
                  fillPattern=FillPattern.VerticalCylinder,
                  fillColor={135,135,135}),Polygon(
                  points={{20,0},{42,0},{42,80},{-42,80},{-42,0},{-20,0},{-20,
                40},{0,60},{20,60},{20,0}},
                  smooth=Smooth.None,
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid,
                  pattern=LinePattern.None),Polygon(
                  points={{20,0},{42,0},{42,-80},{-42,-80},{-42,0},{-20,0},{-20,
                -60},{0,-60},{20,-40},{20,0}},
                  smooth=Smooth.None,
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid,
                  pattern=LinePattern.None),Rectangle(extent={{-20,40},{0,-60}},
              lineColor={0,0,0}),Polygon(
                  points={{-20,40},{0,60},{20,60},{0,40},{-20,40}},
                  lineColor={0,0,0},
                  smooth=Smooth.None),Polygon(
                  points={{20,60},{0,40},{0,-60},{20,-40},{20,60}},
                  lineColor={0,0,0},
                  fillColor={95,95,95},
                  fillPattern=FillPattern.Solid),Line(
                  points={{-20,0},{-100,0}},
                  color={0,0,240},
                  visible=inclTransX,
                  thickness=0.5),Line(
                  points={{10,0},{100,0}},
                  color={127,127,127},
                  visible=inclTransX,
                  thickness=0.5),Ellipse(
                  extent={{-4,52},{4,48}},
                  lineColor={135,135,135},
                  fillColor={0,0,240},
                  visible=inclTransY,
                  fillPattern=FillPattern.Sphere),Line(
                  points={{0,-60},{0,-100}},
                  color={0,0,240},
                  visible=inclTransY,
                  smooth=Smooth.None,
                  thickness=0.5),Line(
                  points={{0,100},{0,50}},
                  color={0,0,240},
                  visible=inclTransY,
                  smooth=Smooth.None,
                  thickness=0.5),Line(
                  points={{-50,-50},{-10,-10}},
                  color={0,0,240},
                  visible=inclTransZ,
                  smooth=Smooth.None,
                  thickness=0.5),Line(
                  points={{20,20},{50,50}},
                  color={0,0,240},
                  visible=inclTransZ,
                  smooth=Smooth.None,
                  thickness=0.5),Text(
                  extent={{-100,60},{100,100}},
                  textString="%name",
                  visible=not inclTransY,
                  lineColor={0,0,0})}));

    end CaFP;

  end CaFPs;

  model Region "Base model for a 3D array of subregions"
    import FCSys.Utilities.Coordinates.cartWrap;
    // extends FCSys.Icons.Names.Top3;
    // extends FCSys.Icons.Names.Top6;

    // Geometric parameters
    parameter Q.Length L_x[:]={U.cm} "Lengths along the x axis" annotation (
        Dialog(group="Geometry",__Dymola_label=
            "<html><i>L</i><sub>x</sub></html>"));
    parameter Q.Length L_y[:]={U.cm} "Lengths along the y axis" annotation (
        Dialog(group="Geometry",__Dymola_label=
            "<html><i>L</i><sub>y</sub></html>"));
    parameter Q.Length L_z[:]={U.cm} "Lengths across the z axis" annotation (
        Dialog(group="Geometry",__Dymola_label=
            "<html><i>L</i><sub>z</sub></html>"));

    final parameter Integer n_x=size(L_x, 1)
      "Number of sets of subregions along the x axis"
      annotation (HideResult=true);
    final parameter Integer n_y=size(L_y, 1)
      "Number of sets of subregions along the y axis"
      annotation (HideResult=true);
    final parameter Integer n_z=size(L_z, 1)
      "Number of sets of subregions along the z axis"
      annotation (HideResult=true);

    // Assumptions
    // -----------
    // Included transport axes
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

    // Auxiliary parameters (for analysis only)
    final parameter Q.Length L[Axis]={sum(L_x),sum(L_y),sum(L_z)} if
      hasSubregions "Length";
    final parameter Q.Area A[Axis]={L[cartWrap(axis + 1)]*L[cartWrap(axis + 2)]
        for axis in Axis} if hasSubregions "Cross-sectional areas";
    final parameter Q.Volume V=product(L) if hasSubregions "Volume";

    replaceable model Subregion = Subregions.Subregion constrainedby
      FCSys.Subregions.PartialSubregion(
      final inclTransX=inclTransX,
      final inclTransY=inclTransY,
      final inclTransZ=inclTransZ) "Base subregion model" annotation (
        __Dymola_choicesAllMatching=true);
    // Note:  In Dymola 2014, the inclTransX, inclTransY, and inclTransZ parameters
    // still appear in the parameter dialog even though they are final.

    Subregion subregions[n_x, n_y, n_z](final L={{L_x[i_x],L_y[i_y],L_z[i_z]}
          for i_z in 1:n_z, i_y in 1:n_y, i_x in 1:n_x}) if hasSubregions
      "Instances of the subregion model"
      annotation (Placement(transformation(extent={{-20,-20},{0,0}})));
    Connectors.BoundaryBus xNegative[n_y, n_z] if inclTransX
      "Negative boundary along the x axis" annotation (Placement(transformation(
            extent={{-60,-20},{-40,0}}), iconTransformation(extent={{-110,-10},
              {-90,10}})));
    Connectors.BoundaryBus xPositive[n_y, n_z] if inclTransX
      "Positive boundary along the x axis" annotation (Placement(transformation(
            extent={{20,-20},{40,0}}), iconTransformation(extent={{90,-10},{110,
              10}})));
    Connectors.BoundaryBus yNegative[n_x, n_z] if inclTransY
      "Negative boundary along the y axis" annotation (Placement(transformation(
            extent={{-20,-60},{0,-40}}), iconTransformation(extent={{-10,-110},
              {10,-90}})));
    Connectors.BoundaryBus yPositive[n_x, n_z] if inclTransY
      "Positive boundary along the y axis" annotation (Placement(transformation(
            extent={{-20,20},{0,40}}), iconTransformation(extent={{-10,90},{10,
              110}})));
    Connectors.BoundaryBus zNegative[n_x, n_y] if inclTransZ
      "Negative boundary along the z axis" annotation (Placement(transformation(
            extent={{0,0},{20,20}}), iconTransformation(extent={{40,40},{60,60}})));
    Connectors.BoundaryBus zPositive[n_x, n_y] if inclTransZ
      "Positive boundary along the z axis" annotation (Placement(transformation(
            extent={{-40,-40},{-20,-20}}), iconTransformation(extent={{-60,-60},
              {-40,-40}})));

  protected
    final parameter Boolean hasSubregions=n_x > 0 and n_y > 0 and n_z > 0
      "true, if there are any subregions";

  equation
    // X axis
    connect(xNegative, subregions[1, :, :].xNegative) annotation (Line(
        points={{-50,-10},{-20,-10}},
        color={127,127,127},
        smooth=Smooth.None,
        thickness=0.5));
    if inclTransX then
      for i in 1:n_x - 1 loop
        connect(subregions[i, :, :].xPositive, subregions[i + 1, :, :].xNegative)
          "Connection b/w neighboring subregions (not shown in the diagram)";
      end for;
    end if;
    connect(subregions[n_x, :, :].xPositive, xPositive) annotation (Line(
        points={{0,-10},{30,-10}},
        color={127,127,127},
        smooth=Smooth.None,
        thickness=0.5));
    if n_x == 0 then
      connect(xNegative, xPositive)
        "Direct pass-through (not shown in the diagram)";
    end if;

    // Y axis
    connect(yNegative, subregions[:, 1, :].yNegative) annotation (Line(
        points={{-10,-50},{-10,-20}},
        color={127,127,127},
        smooth=Smooth.None,
        thickness=0.5));
    if inclTransY then
      for i in 1:n_y - 1 loop
        connect(subregions[:, i, :].yPositive, subregions[:, i + 1, :].yNegative)
          "Connection b/w neighboring subregions (not shown in the diagram)";
      end for;
    end if;
    connect(subregions[:, n_y, :].yPositive, yPositive) annotation (Line(
        points={{-10,0},{-10,30}},
        color={127,127,127},
        smooth=Smooth.None,
        thickness=0.5));

    if n_y == 0 then
      connect(yNegative, yPositive)
        "Direct pass-through (not shown in the diagram)";
    end if;

    // Z axis
    connect(zNegative, subregions[:, :, 1].zNegative) annotation (Line(
        points={{10,10},{-5,-5}},
        color={127,127,127},
        smooth=Smooth.None,
        thickness=0.5));
    if inclTransZ then
      for i in 1:n_z - 1 loop
        connect(subregions[:, :, i].zPositive, subregions[:, :, i + 1].zNegative)
          "Connection b/w neighboring subregions (not shown in the diagram)";
      end for;
    end if;
    connect(zPositive, subregions[:, :, n_z].zPositive) annotation (Line(
        points={{-30,-30},{-15,-15}},
        color={127,127,127},
        smooth=Smooth.None,
        thickness=0.5));
    if n_z == 0 then
      connect(zNegative, zPositive)
        "Direct pass-through (not shown in the diagram)";
    end if;
    // TODO:  Once primitivesVisible is supported by Modelica tools (not
    // supported as of Dymola 2014), complete the icon of this model.
    // Until then, the icon should be blank so that the layer models (AnFP, AnGDL, etc.)
    // are not affected.
    annotation (
      Documentation(info="<html>
  <p>If <i>L</i><sub>x</sub> is an empty vector (e.g., <code>zeros(0)</code>,
  <code>ones(0)</code>, or <code>fill(1, 0)</code>), then there are no
  subregions along the x axis and the boundaries along the x axis are
  directly connected.  The same applies to the other axes.</li></p>
</html>"),
      Diagram(coordinateSystem(preserveAspectRatio=true,extent={{-60,-60},{40,
              40}}), graphics),
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
              100}}), graphics={Text(
            extent={{-100,120},{100,160}},
            textString="%name",
            visible=inclTransY,
            lineColor={0,0,0}), Text(
            extent={{-100,56},{100,96}},
            textString="%name",
            visible=not inclTransY,
            lineColor={0,0,0})}));
  end Region;
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

end Regions;

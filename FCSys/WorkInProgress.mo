within FCSys;
package WorkInProgress "Incomplete classes under development"
  extends Modelica.Icons.Package;


  model CellMSL
    "<html>Single-cell PEMFC with interfaces from the <a href=\"modelica://Modelica\">Modelica Standard Library</a></html>"
    extends FCSys.Icons.Cell;

    extends Modelica.Icons.UnderConstruction;
    Assemblies.Cells.Examples.Cell cell(anFP(redeclare
          FCSys.Subregions.Subregion subregions(
          each final inclX=true,
          each inclY=true,
          each graphite('incle-'=true, 'e-'(perfectMaterialDiff={{{{true,false}}}})),

          each gas(inclH2=true, inclH2O=true))))
      annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

    inner FCSys.Conditions.Environment environment(analysis=false)
      annotation (Placement(transformation(extent={{40,60},{60,80}})));
    Conditions.Adapters.Phases.Graphite caModelicaAdapt(A=cell.L_y[1]*cell.L_z[
          1]) annotation (Placement(transformation(extent={{20,-10},{40,10}})));
    Conditions.Adapters.Phases.Graphite anModelicaAdapt(A=cell.L_y[1]*cell.L_z[
          1])
      annotation (Placement(transformation(extent={{-20,-10},{-40,10}})));
    FCSys.WorkInProgress.TanConduct tanConduct
      annotation (Placement(transformation(extent={{10,40},{-10,60}})));
    Modelica.Blocks.Sources.Ramp loadSweep(duration=1000)
      "This is the arctangent of conductance."
      annotation (Placement(transformation(extent={{-40,60},{-20,80}})));

  equation
    connect(anModelicaAdapt.normal, cell.anFPX[1, 1]) annotation (Line(
        points={{-20,6.10623e-16},{-16,6.10623e-16},{-16,5.55112e-16},{-10,
            5.55112e-16}},
        color={0,0,0},
        thickness=0.5,
        smooth=Smooth.None));

    connect(caModelicaAdapt.normal, cell.caFPX[1, 1]) annotation (Line(
        points={{20,6.10623e-16},{16,6.10623e-16},{16,5.55112e-16},{10,
            5.55112e-16}},
        color={0,0,0},
        thickness=0.5,
        smooth=Smooth.None));

    connect(loadSweep.y, tanConduct.atanGstar) annotation (Line(
        points={{-19,70},{-6.66134e-16,70},{-6.66134e-16,61}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(tanConduct.p, caModelicaAdapt.pin) annotation (Line(
        points={{10,50},{60,50},{60,4},{40,4}},
        color={0,0,255},
        smooth=Smooth.None));
    connect(tanConduct.n, anModelicaAdapt.pin) annotation (Line(
        points={{-10,50},{-60,50},{-60,4},{-40,4}},
        color={0,0,255},
        smooth=Smooth.None));
    annotation (defaultComponentName="cell",experiment(StopTime=1000));
  end CellMSL;

  function plot "Create plots using FCRes"
    extends Modelica.Icons.Function;
    extends Modelica.Icons.UnderConstruction;

  algorithm
    Modelica.Utilities.System.command("loadres");

  end plot;

  model ElectroOsmoticDrag
    "<html>Example to calibrate the coupling between H<sup>+</sup> and H<sub>2</sub>O in the PEM</html>"

    extends Regions.Examples.CLtoCL(anCL(redeclare model Subregion =
            Subregions.Subregion (ionomer(redeclare FCSys.Species.H2O.Gas.Fixed
                H2O(
                redeclare package Data = FCSys.Characteristics.'H+'.Ionomer,
                p_IC=65536*U.kPa,
                consMaterial=ConsThermo.IC)))), caCL(redeclare model Subregion
          = Subregions.Subregion (ionomer(redeclare FCSys.Species.H2O.Gas.Fixed
                H2O(
                redeclare package Data = FCSys.Characteristics.'H+'.Ionomer,
                p_IC=65536*U.kPa,
                consMaterial=ConsThermo.IC)))));

    extends Modelica.Icons.UnderConstruction;

  end ElectroOsmoticDrag;





  model RegionsExamplesCLtoCLVoltage
    "Test one catalyst layer to the other, with prescribed voltage"

    extends Modelica.Icons.Example;
    extends Modelica.Icons.UnderConstruction;
    output Q.Potential w=anCL.subregions[1, 1, 1].graphite.'e-'.g_boundaries[1,
        Side.n] - caCL.subregions[1, 1, 1].graphite.'e-'.g_boundaries[1, Side.p]
      if environment.analysis "Electrical potential";
    output Q.Current zI=-sum(anCL.subregions[1, :, :].graphite.'e-'.boundaries[
        1, Side.n].Ndot) if environment.analysis "Electrical current";
    output Q.Number J_Apercm2=zI*U.cm^2/(anCL.A[Axis.x]*U.A)
      "Electrical current density, in A/cm2";

    parameter Q.Length L_y[:]={8}*U.cm "**Lengths in the y direction";
    parameter Q.Length L_z[:]={6.25}*U.cm "**Lengths in the z direction";
    Regions.AnCLs.AnCL anCL(
      final L_y=L_y,
      final L_z=L_z,
      subregions(graphite(each inclDL=true, 'e-Transfer'(each fromI=false))))
      annotation (Placement(transformation(extent={{-30,30},{-10,50}})));

    Regions.PEMs.PEM PEM(
      final L_y=L_y,
      final L_z=L_z,
      subregions(ionomer('H+'(each consTransX=ConsTrans.dynamic))))
      annotation (Placement(transformation(extent={{-10,30},{10,50}})));
    Regions.CaCLs.CaCL caCL(
      final L_y=L_y,
      final L_z=L_z,
      subregions(graphite(each inclDL=true, 'e-Transfer'(each fromI=false)),
          each ORR('e-'(reaction(Ndot(stateSelect=StateSelect.always))))))
      annotation (Placement(transformation(extent={{10,30},{30,50}})));

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
              Conditions.ByConnector.Boundary.Single.Thermal.heatRate,
          thermalSet(y=0))), each graphite(
        'inclC+'=true,
        'incle-'=true,
        redeclare Conditions.ByConnector.ThermalDiffusive.Single.Temperature
          'C+'(set(y=environment.T)),
        'e-'(
          redeclare function thermalSpec =
              Conditions.ByConnector.Boundary.Single.Thermal.temperature,
          thermalSet(y=environment.T),
          redeclare function materialMeas =
              Conditions.ByConnector.Boundary.Single.Material.potential (
                redeclare package Data = FCSys.Characteristics.'e-'.Graphite))))
      annotation (Placement(transformation(
          extent={{-10,10},{10,-10}},
          rotation=270,
          origin={-44,40})));

    Conditions.ByConnector.BoundaryBus.Single.Source caBC[caCL.n_y, caCL.n_z](
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
          thermalSet(y=environment.T))), graphite(
        each 'inclC+'=true,
        each 'incle-'=true,
        each 'C+'(set(y=environment.T)),
        'e-'(
          redeclare each function materialSpec =
              Conditions.ByConnector.Boundary.Single.Material.potential (
                redeclare package Data = FCSys.Characteristics.'e-'.Graphite),
          materialSet(y=anBC.graphite.'e-'.materialOut.y + fill(
                  -voltageSet.y,
                  caCL.n_y,
                  caCL.n_z)),
          each thermalSet(y=environment.T)))) annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={44,40})));

    Modelica.Blocks.Sources.Ramp voltageSet(
      duration=300,
      offset=1.19997*U.V,
      height=-0.8*U.V)
      annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
    inner Conditions.Environment environment(
      analysis=true,
      T=333.15*U.K,
      p=U.from_kPag(48.3),
      RH=0.7) "Environmental conditions"
      annotation (Placement(transformation(extent={{-10,70},{10,90}})));
  equation
    connect(anCL.xPositive, PEM.xNegative) annotation (Line(
        points={{-10,40},{-10,40}},
        color={240,0,0},
        thickness=0.5,
        smooth=Smooth.None));
    connect(PEM.xPositive, caCL.xNegative) annotation (Line(
        points={{10,40},{10,40}},
        color={0,0,240},
        thickness=0.5,
        smooth=Smooth.None));
    connect(anBC.boundary, anCL.xNegative) annotation (Line(
        points={{-40,40},{-30,40}},
        color={240,0,0},
        thickness=0.5,
        smooth=Smooth.None));
    connect(caCL.xPositive, caBC.boundary) annotation (Line(
        points={{30,40},{40,40}},
        color={0,0,240},
        thickness=0.5,
        smooth=Smooth.None));
    annotation (
      Commands(file=
            "Resources/Scripts/Dymola/Regions.Examples.CLtoCLVoltage.mos"
          "Regions.Examples.CLtoCLVoltage.mos"),
      experiment(
        StopTime=600,
        Tolerance=1e-007,
        __Dymola_Algorithm="Dassl"),
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
              100,100}}), graphics));
  end RegionsExamplesCLtoCLVoltage;



  annotation (Commands(
      file="../../units.mos"
        "Establish the constants and units in the workspace (first translate a model besides Units.Evaluate).",

      file="test/check.mos" "Check all of FCSys using Dymola's check function.",

      file="../../../LaTeX/Dissertation/Results/Cell/Simulation/sim.mos"), Icon(
        graphics={Polygon(
          points={{-80,-72},{0,68},{80,-72},{-80,-72}},
          lineColor={255,0,0},
          lineThickness=0.5)}));

end WorkInProgress;

within FCSys;
package Sensors "Models to measure conditions"
  extends Modelica.Icons.SensorsPackage;

  package ChemicalBus
    "<html>Sensors for the <a href=\"modelica://FCSys.Connectors.ChemicalBus\">ChemicalBus</a> connector</html>"
    extends Modelica.Icons.Package;

    model Gas "Sensor for gas"

      extends BaseClasses.NullPhase;

      // Conditionally include species.
      parameter Boolean inclH2=false "<html>Hydrogen (H<sub>2</sub>)</html>"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      Chemical.Species H2(final n_lin=n_lin) if inclH2 "Model"
        annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

      parameter Boolean inclH2O=false "<html>Water (H<sub>2</sub>O)</html>"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      Chemical.Species H2O(final n_lin=n_lin) if inclH2O "Model"
        annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

      parameter Boolean inclN2=false "<html>Nitrogen (N<sub>2</sub>)</html>"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      Chemical.Species N2(final n_lin=n_lin) if inclN2 "Model" annotation (
          Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          enable=inclN2), Placement(transformation(extent={{-10,-50},{10,-30}})));

      parameter Boolean inclO2=false "<html>Oxygen (O<sub>2</sub>)</html>"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      Chemical.Species O2(final n_lin=n_lin) if inclO2 "Model"
        annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

    equation
      // H2
      connect(H2.chemical, chemical.H2) annotation (Line(
          points={{6.10623e-16,-40},{5.55112e-16,-40},{5.55112e-16,5.55112e-16}},

          color={208,104,0},
          smooth=Smooth.None));

      connect(H2.muPerT, y.H2.muPerT) annotation (Line(
          points={{-4,-50},{-4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(H2.phi, y.H2.phi) annotation (Line(
          points={{6.10623e-16,-50},{5.55112e-16,-50},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(H2.hbar, y.H2.hbar) annotation (Line(
          points={{4,-50},{4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      // H2O
      connect(H2O.chemical, chemical.H2O) annotation (Line(
          points={{6.10623e-16,-40},{5.55112e-16,-40},{5.55112e-16,5.55112e-16}},

          color={208,104,0},
          smooth=Smooth.None));

      connect(H2O.muPerT, y.H2O.muPerT) annotation (Line(
          points={{-4,-50},{-4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(H2O.phi, y.H2O.phi) annotation (Line(
          points={{6.10623e-16,-50},{5.55112e-16,-50},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(H2O.hbar, y.H2O.hbar) annotation (Line(
          points={{4,-50},{4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      // N2
      connect(N2.chemical, chemical.N2) annotation (Line(
          points={{6.10623e-16,-40},{5.55112e-16,-40},{5.55112e-16,5.55112e-16}},

          color={208,104,0},
          smooth=Smooth.None));

      connect(N2.muPerT, y.N2.muPerT) annotation (Line(
          points={{-4,-50},{-4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(N2.phi, y.N2.phi) annotation (Line(
          points={{6.10623e-16,-50},{-4.87687e-22,-50},{-4.87687e-22,-100},{
              5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(N2.hbar, y.N2.hbar) annotation (Line(
          points={{4,-50},{4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      // O2
      connect(O2.chemical, chemical.O2) annotation (Line(
          points={{6.10623e-16,-40},{5.55112e-16,-40},{5.55112e-16,5.55112e-16}},

          color={208,104,0},
          smooth=Smooth.None));

      connect(O2.muPerT, y.O2.muPerT) annotation (Line(
          points={{-4,-50},{-4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(O2.phi, y.O2.phi) annotation (Line(
          points={{6.10623e-16,-50},{-4.87687e-22,-50},{-4.87687e-22,-100},{
              5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(O2.hbar, y.O2.hbar) annotation (Line(
          points={{4,-50},{4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      annotation (defaultComponentPrefixes="replaceable", defaultComponentName=
            "phaseChemicalSensor");
    end Gas;

    model Graphite "Sensor for graphite"

      extends BaseClasses.NullPhase;

      // Conditionally include species.
      parameter Boolean inclC=false "Carbon (C)" annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      replaceable Chemical.Species C(final n_lin=n_lin) if inclC "Model"
        annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

      parameter Boolean 'incle-'=false "<html>Electrons (e<sup>-</sup>)</html>"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      Chemical.Species 'e-'(final n_lin=n_lin) if 'incle-' "Model"
        annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

    equation
      // C
      connect(C.chemical, chemical.C) annotation (Line(
          points={{6.10623e-16,-40},{5.55112e-16,5.55112e-16}},
          color={208,104,0},
          smooth=Smooth.None));
      connect(C.muPerT, y.C.muPerT) annotation (Line(
          points={{-4,-50},{-4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(C.phi, y.C.phi) annotation (Line(
          points={{6.10623e-16,-50},{5.55112e-16,-50},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(C.hbar, y.C.hbar) annotation (Line(
          points={{4,-50},{4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      // e-
      connect('e-'.chemical, chemical.'e-') annotation (Line(
          points={{6.10623e-16,-40},{5.55112e-16,-40},{5.55112e-16,5.55112e-16}},

          color={208,104,0},
          smooth=Smooth.None));

      connect('e-'.muPerT, y.'e-'.muPerT) annotation (Line(
          points={{-4,-50},{-4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));
      connect('e-'.phi, y.'e-'.phi) annotation (Line(
          points={{6.10623e-16,-50},{5.55112e-16,-50},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      connect('e-'.hbar, y.'e-'.hbar) annotation (Line(
          points={{4,-50},{4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      annotation (defaultComponentPrefixes="replaceable", defaultComponentName=
            "phaseChemicalSensor");
    end Graphite;

    model Ionomer "Sensor for ionomer"

      extends BaseClasses.NullPhase;

      // Conditionally include species.
      parameter Boolean inclC19HF37O5S=false
        "<html>Nafion sulfonate (C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S)</html>"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      Chemical.Species C19HF37O5S(final n_lin=n_lin) if inclC19HF37O5S "Model"
        annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

      parameter Boolean inclH2O=false "<html>Water (H<sub>2</sub>O)</html>"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      Chemical.Species H2O(final n_lin=n_lin) if inclH2O "Model"
        annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

      parameter Boolean 'inclH+'=false "<html>Protons (H<sup>+</sup>)</html>"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      Chemical.Species 'H+'(final n_lin=n_lin) if 'inclH+' "Model"
        annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

    equation
      // C19HF37O5S
      connect(C19HF37O5S.chemical, chemical.C19HF37O5S) annotation (Line(
          points={{6.10623e-16,-40},{5.55112e-16,-40},{5.55112e-16,5.55112e-16}},

          color={208,104,0},
          smooth=Smooth.None));

      connect(C19HF37O5S.muPerT, y.C19HF37O5S.muPerT) annotation (Line(
          points={{-4,-50},{-4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(C19HF37O5S.phi, y.C19HF37O5S.phi) annotation (Line(
          points={{6.10623e-16,-50},{5.55112e-16,-50},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(C19HF37O5S.hbar, y.C19HF37O5S.hbar) annotation (Line(
          points={{4,-50},{4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      // H+
      connect('H+'.chemical, chemical.'H+') annotation (Line(
          points={{6.10623e-16,-40},{5.55112e-16,-40},{5.55112e-16,5.55112e-16}},

          color={208,104,0},
          smooth=Smooth.None));

      connect('H+'.muPerT, y.'H+'.muPerT) annotation (Line(
          points={{-4,-50},{-4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));
      connect('H+'.phi, y.'H+'.phi) annotation (Line(
          points={{6.10623e-16,-50},{5.55112e-16,-50},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      connect('H+'.hbar, y.'H+'.hbar) annotation (Line(
          points={{4,-50},{4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      // H2O
      connect(H2O.chemical, chemical.H2O) annotation (Line(
          points={{6.10623e-16,-40},{5.55112e-16,-40},{5.55112e-16,5.55112e-16}},

          color={208,104,0},
          smooth=Smooth.None));

      connect(H2O.muPerT, y.H2O.muPerT) annotation (Line(
          points={{-4,-50},{-4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(H2O.phi, y.H2O.phi) annotation (Line(
          points={{6.10623e-16,-50},{5.55112e-16,-50},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(H2O.hbar, y.H2O.hbar) annotation (Line(
          points={{4,-50},{4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      annotation (
        defaultComponentPrefixes="replaceable",
        defaultComponentName="phaseChemicalSensor",
        Diagram(graphics));
    end Ionomer;

    model Liquid "Sensor for liquid"

      extends BaseClasses.NullPhase;

      // Conditionally include species.
      parameter Boolean inclH2O=false "<html>Water (H<sub>2</sub>O)</html>"
        annotation (
        Evaluate=true,
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          group="Species",
          __Dymola_descriptionLabel=true,
          __Dymola_joinNext=true));
      Chemical.Species H2O(final n_lin=n_lin) if inclH2O "Model"
        annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

    equation
      // H2O
      connect(H2O.chemical, chemical.H2O) annotation (Line(
          points={{6.10623e-16,-40},{5.55112e-16,-40},{5.55112e-16,5.55112e-16}},

          color={208,104,0},
          smooth=Smooth.None));

      connect(H2O.muPerT, y.H2O.muPerT) annotation (Line(
          points={{-4,-50},{-4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(H2O.phi, y.H2O.phi) annotation (Line(
          points={{6.10623e-16,-50},{5.55112e-16,-50},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(H2O.hbar, y.H2O.hbar) annotation (Line(
          points={{4,-50},{4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      annotation (defaultComponentPrefixes="replaceable", defaultComponentName=
            "phaseChemicalSensor");
    end Liquid;

    package BaseClasses "Base classes (not for direct use)"
      extends Modelica.Icons.BasesPackage;
      model NullPhase "Empty sensor for a phase (no species)"
        extends FCSys.BaseClasses.Icons.Sensor;
        parameter Integer n_lin(
          final min=1,
          final max=3) = 1
          "<html>Number of components of linear momentum (<i>n</i><sub>lin</sub>)</html>"
          annotation (HideResult=true);

        FCSys.Connectors.ChemicalBus chemical
          "Multi-species connector for material" annotation (Placement(
              transformation(extent={{-10,-10},{10,10}}), iconTransformation(
                extent={{-10,-10},{10,10}})));

        FCSys.Connectors.RealOutputBus y "Bus of measurements" annotation (
            Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,-100}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,-100})));
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="phaseChemicalSensor",
          Icon(graphics));
      end NullPhase;
    end BaseClasses;
    annotation (Documentation(info="<html><p>Since the connectors in
<a href=\"modelica://FCSys\">FCSys</a> are hierarchical
(see the <a href=\"modelica://FCSys.Connectors\">Connectors</a> package),
the models for the sensors must be as well.  A
<a href=\"modelica://FCSys.Connectors.ChemicalOutput\">ChemicalOutput</a>
connector
is used in <a href=\"modelica://FCSys.Subregions.Species\">Species</a> models,
and there is a corresponding <a href=\"modelica://FCSys.Sensors.Chemical.Species\">Species
sensor</a> model in this package. The
<a href=\"modelica://FCSys.Connectors.ChemicalBus\">ChemicalBus</a>
connector is used in <a href=\"modelica://FCSys.Subregions.Phases\">Phase</a> models,
and there are corresponding <a href=\"modelica://FCSys.Sensors.ChemicalBus\">ChemicalBus
sensor</a> models.
</p></html>"));
  end ChemicalBus;

  package Chemical
    "<html>Sensors for the <a href=\"modelica://FCSys.Connectors.ChemicalInput\">ChemicalOutput</a> connector</html>"

    extends Modelica.Icons.Package;
    model Species
      "<html>Sensors for the <a href=\"modelica://FCSys.Connectors.ChemicalInput\">ChemicalInput</a> connector</html>"

      extends FCSys.BaseClasses.Icons.Sensor;

      parameter Integer n_lin(
        final min=1,
        final max=3) = 1
        "<html>Number of components of linear momentum (<i>n</i><sub>lin</sub>)</html>"
        annotation (HideResult=true);

      FCSys.Connectors.ChemicalInput chemical(final n_lin=n_lin)
        "Single-species connector for material" annotation (Placement(
            transformation(extent={{-10,-10},{10,10}}), iconTransformation(
              extent={{-10,-10},{10,10}})));

      FCSys.Connectors.RealOutput muPerT(final unit="1")
        "Internal signal for quotient of electrochemical potential and temperature"
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-30,-100}),iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-40,-100})));
      FCSys.Connectors.RealOutput phi[n_lin](each final unit="l/T")
        "Internal signal for velocity of the source" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,-100}),iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,-100})));
      FCSys.Connectors.RealOutput hbar(final unit="l2/T2")
        "Internal signal for massic enthalpy of the source" annotation (
          Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={30,-100}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={40,-100})));

    equation
      // Measurements
      muPerT = chemical.muPerT;
      phi = chemical.phi;
      hbar = chemical.hbar;

      // Conservation (no storage)
      0 = chemical.Ndot "Material";
      zeros(n_lin) = chemical.mPhidot "Linear momentum";
      0 = chemical.Hdot "Energy";
      annotation (
        defaultComponentName="speciesChemicalSensor",
        Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                {100,100}}), graphics),
        Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
                100,100}}), graphics={Line(points={{-40,-100},{-40,-58}}, color
              ={0,0,127}),Line(points={{40,-100},{40,-58}}, color={0,0,127}),
              Text(
                  extent={{-100,-20},{100,-50}},
                  lineColor={127,127,127},
                  textString="mu/T   phi   T     ")}));
    end Species;
  end Chemical;

  package InertAmagat
    "<html>Sensors for the <a href=\"modelica://FCSys.Connectors.InertAmagat\">InertAmagat</a> connector</html>"
    extends Modelica.Icons.Package;

    model Phase
      "<html>Sensor for the <a href=\"modelica://FCSys.Connectors.InertAmagat\">InertAmagat</a> connector, e.g., of a <a href=\"modelica://FCSys.Subregions.Phases\">Phase</a> model</html>"
      extends FCSys.BaseClasses.Icons.Sensor;

      parameter Integer n_lin(
        final min=0,
        final max=3) = 1
        "<html>Number of components of linear momentum (<i>n</i><sub>lin</sub>)</html>";

      Pressure pressure(final n_lin=n_lin) "Sensor for pressure"
        annotation (Placement(transformation(extent={{-50,-60},{-30,-40}})));

      Velocity velocity(final n_lin=n_lin) if n_lin > 1 "Sensor for velocity"
        annotation (Placement(transformation(extent={{-10,-60},{10,-40}})));

      Temperature temperature(final n_lin=n_lin) "Sensor for temperature"
        annotation (Placement(transformation(extent={{30,-60},{50,-40}})));
      FCSys.Connectors.InertAmagat inert(final n_lin=n_lin)
        "Single-species connector for linear momentum and heat, with additivity of volume"
        annotation (Placement(transformation(extent={{-10,-10},{10,10}}),
            iconTransformation(extent={{-10,-10},{10,10}})));
      FCSys.Connectors.RealOutputBus y "Input bus for external signal sources"
        annotation (HideResult=not (internalPress or internalLin1 or
            internalLin2 or internalLin3 or internalEntropy), Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,-100}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,-100})));

    equation
      connect(pressure.inert, inert) annotation (Line(
          points={{-40,-50},{-40,-30},{5.55112e-16,-30},{5.55112e-16,
              5.55112e-16}},
          color={72,90,180},
          smooth=Smooth.None));

      connect(pressure.y, y.p) annotation (Line(
          points={{-40,-60},{-40,-70},{5.55112e-16,-70},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%second",
          index=1,
          extent={{6,3},{6,3}}));
      connect(velocity.inert, inert) annotation (Line(
          points={{6.10623e-16,-50},{6.10623e-16,-30},{0,-30},{0,5.55112e-16},{
              5.55112e-16,5.55112e-16}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(velocity.y, y.phi) annotation (Line(
          points={{6.10623e-16,-60},{6.10623e-16,-70},{5.55112e-16,-70},{
              5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%second",
          index=1,
          extent={{6,3},{6,3}}));
      connect(temperature.inert, inert) annotation (Line(
          points={{40,-50},{40,-30},{5.55112e-16,-30},{5.55112e-16,5.55112e-16}},

          color={72,90,180},
          smooth=Smooth.None));

      connect(temperature.y, y.T) annotation (Line(
          points={{40,-60},{40,-70},{5.55112e-16,-70},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%second",
          index=1,
          extent={{6,3},{6,3}}));
      annotation (defaultComponentName="phaseInertSensor");
    end Phase;

    model Pressure "Prescribed pressure"
      extends BaseClasses.PartialSensor(redeclare FCSys.Connectors.RealOutput y(
            final unit="m/(l.T2)"));
    equation
      y = inert.p;
      annotation (
        defaultComponentPrefixes="replaceable",
        defaultComponentName="volSensor",
        Icon(graphics={Text(
                  extent={{-100,-20},{100,-50}},
                  lineColor={127,127,127},
                  textString="p")}));
    end Pressure;

    model Velocity "Measured velocity"
      extends FCSys.BaseClasses.Icons.Sensor;

      parameter Integer n_lin(
        final min=0,
        final max=3) = 1
        "<html>Number of components of linear momentum (<i>n</i><sub>lin</sub>)</html>"
        annotation (HideResult=true);

      FCSys.Connectors.InertDalton inert(final n_lin=n_lin)
        "Connector for linear momentum and heat, with additivity of pressure"
        annotation (HideResult=true, Placement(transformation(extent={{-10,-10},
                {10,10}}), iconTransformation(extent={{-10,-10},{10,10}})));
      FCSys.Connectors.RealOutput y[n_lin](final unit="l/T") "Measurement"
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,-100})));

    equation
      // Measurement
      y = inert.phi;

      // Conservation (no storage)
      0 = inert.V "No (additional) volume";
      zeros(n_lin) = inert.mPhidot "No force";
      0 = inert.Qdot "Adiabatic";

      annotation (
        defaultComponentPrefixes="replaceable",
        defaultComponentName="momSensor",
        Icon(graphics={Text(
                  extent={{-100,-20},{100,-50}},
                  lineColor={127,127,127},
                  textString="phi")}));
    end Velocity;

    model Temperature "Measured temperature"
      extends BaseClasses.PartialSensor(redeclare FCSys.Connectors.RealOutput y(
            final unit="l2.m/(N.T2)", displayUnit="K"));
    equation
      y = inert.T;
      annotation (
        defaultComponentPrefixes="replaceable",
        defaultComponentName="entropySensor",
        Icon(graphics={Text(
                  extent={{-100,-20},{100,-50}},
                  lineColor={127,127,127},
                  textString="T")}));
    end Temperature;

    package BaseClasses "Base classes (not for direct use)"
      extends Modelica.Icons.BasesPackage;
      partial model PartialSensor "Partial model for a sensor"
        extends FCSys.BaseClasses.Icons.Sensor;

        parameter Integer n_lin(
          final min=1,
          final max=3) = 1
          "<html>Number of components of linear momentum (<i>n</i><sub>lin</sub>)</html>"
          annotation (HideResult=true);

        FCSys.Connectors.InertAmagat inert(final n_lin=n_lin)
          "Connector for linear momentum and heat, with additivity of volume"
          annotation (HideResult=true, Placement(transformation(extent={{-10,-10},
                  {10,10}}), iconTransformation(extent={{-10,-10},{10,10}})));
        FCSys.Connectors.RealOutput y "Measurement" annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,-100}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,-100})));
      equation
        // Conservation (no storage)
        0 = inert.V "No (additional) volume";
        zeros(n_lin) = inert.mPhidot "No force";
        0 = inert.Qdot "Adiabatic";
        annotation (defaultComponentName="phaseInertSensor");
      end PartialSensor;
    end BaseClasses;
  end InertAmagat;

  package InertDalton
    "<html>Sensors for the <a href=\"modelica://FCSys.Connectors.InertDalton\">InertDalton</a> connector</html>"
    extends Modelica.Icons.Package;

    model Species
      "<html>Sensor for the <a href=\"modelica://FCSys.Connectors.InertDalton\">InertDalton</a> connector, e.g., of a <a href=\"modelica://FCSys.Subregions.Species\">Species</a> model</html>"

      extends FCSys.BaseClasses.Icons.Sensor;

      parameter Integer n_lin(
        final min=0,
        final max=3) = 1
        "<html>Number of components of linear momentum (<i>n</i><sub>lin</sub>)</html>";

      Volume volume(final n_lin=n_lin) "Sensor for volume"
        annotation (Placement(transformation(extent={{-50,-60},{-30,-40}})));

      Velocity velocity(final n_lin=n_lin) if n_lin > 1 "Sensor for velocity"
        annotation (Placement(transformation(extent={{-10,-60},{10,-40}})));

      Temperature temperature(final n_lin=n_lin) "Sensor for temperature"
        annotation (Placement(transformation(extent={{30,-60},{50,-40}})));

      FCSys.Connectors.InertDalton inert(final n_lin=n_lin)
        "Single-species connector for linear momentum and heat, with additivity of pressure"
        annotation (Placement(transformation(extent={{-10,-10},{10,10}}),
            iconTransformation(extent={{-10,-10},{10,10}})));
      FCSys.Connectors.RealOutputBus y "Input bus for external signal sources"
        annotation (HideResult=not (internalPress or internalLin1 or
            internalLin2 or internalLin3 or internalEntropy), Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,-100}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,-100})));

    equation
      connect(volume.inert, inert) annotation (Line(
          points={{-40,-50},{-40,-30},{5.55112e-16,-30},{5.55112e-16,
              5.55112e-16}},
          color={72,90,180},
          smooth=Smooth.None));

      connect(volume.y, y.V) annotation (Line(
          points={{-40,-60},{-40,-70},{5.55112e-16,-70},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%second",
          index=1,
          extent={{6,3},{6,3}}));
      connect(velocity.inert, inert) annotation (Line(
          points={{6.10623e-16,-50},{6.10623e-16,-30},{0,-30},{0,5.55112e-16},{
              5.55112e-16,5.55112e-16}},
          color={72,90,180},
          smooth=Smooth.None));
      connect(velocity.y, y.phi) annotation (Line(
          points={{6.10623e-16,-60},{6.10623e-16,-70},{5.55112e-16,-70},{
              5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%second",
          index=1,
          extent={{6,3},{6,3}}));
      connect(temperature.inert, inert) annotation (Line(
          points={{40,-50},{40,-30},{5.55112e-16,-30},{5.55112e-16,5.55112e-16}},

          color={72,90,180},
          smooth=Smooth.None));

      connect(temperature.y, y.T) annotation (Line(
          points={{40,-60},{40,-70},{5.55112e-16,-70},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%second",
          index=1,
          extent={{6,3},{6,3}}));
      annotation (defaultComponentName="speciesInertSensor", Diagram(graphics));
    end Species;

    model Volume "Measured volume"
      extends BaseClasses.PartialSensor(redeclare FCSys.Connectors.RealOutput y(
            final unit="l3"));

    equation
      y = inert.V;

      annotation (Icon(graphics={Text(
                  extent={{-100,-20},{100,-50}},
                  lineColor={127,127,127},
                  textString="V")}), Diagram(graphics));
    end Volume;

    model Velocity "Measured velocity"
      extends FCSys.BaseClasses.Icons.Sensor;

      parameter Integer n_lin(
        final min=0,
        final max=3) = 1
        "<html>Number of components of linear momentum (<i>n</i><sub>lin</sub>)</html>"
        annotation (HideResult=true);

      FCSys.Connectors.InertDalton inert(final n_lin=n_lin)
        "Connector for linear momentum and heat, with additivity of pressure"
        annotation (HideResult=true, Placement(transformation(extent={{-10,-10},
                {10,10}}), iconTransformation(extent={{-10,-10},{10,10}})));
      FCSys.Connectors.RealOutput y[n_lin](final unit="l/T") "Measurement"
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,-100})));

    equation
      // Measurement
      y = inert.phi;

      // Conservation (no storage)
      0 = inert.p "No (additional) pressure";
      zeros(n_lin) = inert.mPhidot "No force";
      0 = inert.Qdot "Adiabatic";
      annotation (Icon(graphics={Text(
                  extent={{-100,-20},{100,-50}},
                  lineColor={127,127,127},
                  textString="phi")}), Diagram(graphics));
    end Velocity;

    model Temperature "Measured temperature"

      extends BaseClasses.PartialSensor(redeclare FCSys.Connectors.RealOutput y(
            final unit="l2.m/(N.T2)", displayUnit="K"));

    equation
      y = inert.T;

      annotation (Icon(graphics={Text(
                  extent={{-100,-20},{100,-50}},
                  lineColor={127,127,127},
                  textString="T")}), Diagram(graphics));
    end Temperature;

    package BaseClasses "Base classes (not for direct use)"
      extends Modelica.Icons.BasesPackage;
      partial model PartialSensor "Partial model for a sensor"
        extends FCSys.Sensors.BaseClasses.PartialSensor;

        parameter Integer n_lin(
          final min=0,
          final max=3) = 1
          "<html>Number of components of linear momentum (<i>n</i><sub>lin</sub>)</html>"
          annotation (HideResult=true);

        FCSys.Connectors.InertDalton inert(final n_lin=n_lin)
          "Connector for linear momentum and heat, with additivity of pressure"
          annotation (HideResult=true, Placement(transformation(extent={{-10,-10},
                  {10,10}}), iconTransformation(extent={{-10,-10},{10,10}})));

      equation
        // Conservation (no storage)
        0 = inert.p "No (additional) pressure";
        zeros(n_lin) = inert.mPhidot "No force";
        0 = inert.Qdot "Adiabatic";

        annotation (defaultComponentName="speciesInertSensor");
      end PartialSensor;
    end BaseClasses;
  end InertDalton;

  package FaceBus
    "<html>Sensors for a single <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connector</html>"

    extends Modelica.Icons.Package;
    model Subregion
      "<html>Sensor for a face of a <a href=\"modelica://FCSys.Subregions.Subregion\">Region</a> or <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a> model</html>"

      extends FCSys.BaseClasses.Icons.Sensor;

      parameter Axis axis=Axis.x "Axis normal to the face";

      Phases.Gas gas(final axis=axis) "Gas" annotation (Dialog(group="Phases",
            __Dymola_descriptionLabel=true), Placement(transformation(extent={{
                -10,-50},{10,-30}})));

      Phases.Graphite graphite(final axis=axis) "Graphite" annotation (Dialog(
            group="Phases", __Dymola_descriptionLabel=true), Placement(
            transformation(extent={{-10,-50},{10,-30}})));

      Phases.Ionomer ionomer(final axis=axis) "Ionomer" annotation (Dialog(
            group="Phases", __Dymola_descriptionLabel=true), Placement(
            transformation(extent={{-10,-50},{10,-30}})));

      Phases.Liquid liquid(final axis=axis) "Liquid" annotation (Dialog(group=
              "Phases", __Dymola_descriptionLabel=true), Placement(
            transformation(extent={{-10,-50},{10,-30}})));

      FCSys.Connectors.FaceBus face
        "Connector for linear momentum and heat of multiple species"
        annotation (Placement(transformation(extent={{-10,-10},{10,10}}),
            iconTransformation(extent={{-10,-10},{10,10}})));
      FCSys.Connectors.RealOutputBus y "Bus of measurements" annotation (
          Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,-100})));
    equation
      connect(gas.face, face.gas) annotation (Line(
          points={{6.10623e-16,-40},{6.10623e-16,5.55112e-16},{5.55112e-16,
              5.55112e-16}},
          color={127,127,127},
          pattern=LinePattern.None,
          thickness=0.5,
          smooth=Smooth.None));

      connect(graphite.face, face.graphite) annotation (Line(
          points={{6.10623e-16,-40},{6.10623e-16,5.55112e-16},{5.55112e-16,
              5.55112e-16}},
          color={127,127,127},
          pattern=LinePattern.None,
          thickness=0.5,
          smooth=Smooth.None));

      connect(ionomer.face, face.ionomer) annotation (Line(
          points={{6.10623e-16,-40},{6.10623e-16,5.55112e-16},{5.55112e-16,
              5.55112e-16}},
          color={127,127,127},
          pattern=LinePattern.None,
          thickness=0.5,
          smooth=Smooth.None));

      connect(y.gas, gas.y) annotation (Line(
          points={{5.55112e-16,-100},{6.10623e-16,-50}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None));

      connect(y.graphite, graphite.y) annotation (Line(
          points={{5.55112e-16,-100},{6.10623e-16,-50}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None));

      connect(y.ionomer, ionomer.y) annotation (Line(
          points={{5.55112e-16,-100},{6.10623e-16,-50}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None));

      annotation (
        defaultComponentPrefixes="replaceable",
        defaultComponentName="subregionFaceSensor",
        Diagram(graphics));
    end Subregion;

    package Phases
      "<html>Sensors for the <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connector, e.g., of a <a href=\"modelica://FCSys.Subregions.Phase\">Phase</a> model (multi-species)</html>"
      extends Modelica.Icons.Package;

      model Gas "Sensor for gas"

        extends BaseClasses.NullPhase;

        // Conditionally include species.
        parameter Boolean inclH2=false "<html>Hydrogen (H<sub>2</sub>)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        Face.Species H2(final axis=axis) if inclH2 "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclH2), Placement(transformation(extent={{-10,-40},{10,-20}})));
        parameter Boolean inclH2O=false "<html>Water (H<sub>2</sub>O)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        Face.Species H2O(final axis=axis) if inclH2O "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclH2O), Placement(transformation(extent={{-10,-40},{10,-20}})));
        parameter Boolean inclN2=false "<html>Nitrogen (N<sub>2</sub>)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        Face.Species N2(final axis=axis) if inclN2 "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclN2), Placement(transformation(extent={{-10,-40},{10,-20}})));
        parameter Boolean inclO2=false "<html>Oxygen (O<sub>2</sub>)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        Face.Species O2(final axis=axis) if inclO2 "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclO2), Placement(transformation(extent={{-10,-40},{10,-20}})));
      equation
        // H2
        connect(H2.face.normal, face.H2.normal) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.face.thermal, face.H2.thermal) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.face.transverseX, face.H2.transverseX) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,-30},{5.55112e-16,
                5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(H2.face.transverseY, face.H2.transverseY) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,-30},{5.55112e-16,
                5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(H2.face.transverseZ, face.H2.transverseZ) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,-30},{5.55112e-16,
                5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(y.H2, H2.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-40}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // H2O
        connect(H2O.face.normal, face.H2O.normal) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.face.thermal, face.H2O.thermal) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.face.transverseX, face.H2O.transverseX) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,-30},{5.55112e-16,
                5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(H2O.face.transverseY, face.H2O.transverseY) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,-30},{5.55112e-16,
                5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(H2O.face.transverseZ, face.H2O.transverseZ) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,-30},{5.55112e-16,
                5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(y.H2O, H2O.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-40}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // N2
        connect(N2.face.normal, face.N2.normal) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.face.thermal, face.N2.thermal) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.face.transverseX, face.N2.transverseX) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,-30},{5.55112e-16,
                5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(N2.face.transverseY, face.N2.transverseY) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,-30},{5.55112e-16,
                5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(N2.face.transverseZ, face.N2.transverseZ) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,-30},{5.55112e-16,
                5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(y.N2, N2.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-40}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // O2
        connect(O2.face.normal, face.O2.normal) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.face.thermal, face.O2.thermal) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.face.transverseX, face.O2.transverseX) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,-30},{5.55112e-16,
                5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(O2.face.transverseY, face.O2.transverseY) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,-30},{5.55112e-16,
                5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(O2.face.transverseZ, face.O2.transverseZ) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,-30},{5.55112e-16,
                5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(y.O2, O2.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-40}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));
      end Gas;

      model Graphite "Sensor for graphite"

        extends BaseClasses.NullPhase;

        // Conditionally include species.
        parameter Boolean inclC=false "Carbon (C)" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        Face.Species C(final axis=axis) if inclC "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclC), Placement(transformation(extent={{-10,-40},{10,-20}})));
        parameter Boolean 'incle-'=false
          "<html>Electrons (e<sup>-</sup>)</html>" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        Face.Species 'e-'(final axis=axis) if 'incle-' "Model" annotation (
            Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable='incle-'), Placement(transformation(extent={{-10,-40},{10,-20}})));

      equation
        // C
        connect(C.face.normal, face.C.normal) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.face.thermal, face.C.thermal) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.face.transverseX, face.C.transverseX) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,-30},{5.55112e-16,
                5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(C.face.transverseY, face.C.transverseY) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,-30},{5.55112e-16,
                5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(C.face.transverseZ, face.C.transverseZ) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,-30},{5.55112e-16,
                5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(y.C, C.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-40}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // e-
        connect('e-'.face.normal, face.'e-'.normal) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.face.thermal, face.'e-'.thermal) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.face.transverseX, face.'e-'.transverseX) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,-30},{5.55112e-16,
                5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect('e-'.face.transverseY, face.'e-'.transverseY) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,-30},{5.55112e-16,
                5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect('e-'.face.transverseZ, face.'e-'.transverseZ) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,-30},{5.55112e-16,
                5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(y.'e-', 'e-'.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-40}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));
      end Graphite;

      model Ionomer "Sensor for ionomer"

        extends BaseClasses.NullPhase;

        // Conditionally include species.
        parameter Boolean inclC19HF37O5S=false
          "<html>Nafion sulfonate (C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        Face.Species C19HF37O5S(final axis=axis) if inclC19HF37O5S "Model"
          annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclC19HF37O5S), Placement(transformation(extent={{-10,-40},
                  {10,-20}})));
        parameter Boolean 'inclH+'=false "<html>Protons (H<sup>+</sup>)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        Face.Species 'H+'(final axis=axis) if 'inclH+' "Model" annotation (
            Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable='inclH+'), Placement(transformation(extent={{-10,-40},{10,-20}})));
        parameter Boolean inclH2O=false "<html>Water (H<sub>2</sub>O)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        Face.Species H2O(final axis=axis) if inclH2O "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclH2O), Placement(transformation(extent={{-10,-40},{10,-20}})));

      equation
        // C19HF37O5S
        connect(C19HF37O5S.face.normal, face.C19HF37O5S.normal) annotation (
            Line(
            points={{6.10623e-16,-30},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.face.thermal, face.C19HF37O5S.thermal) annotation (
            Line(
            points={{6.10623e-16,-30},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.face.transverseX, face.C19HF37O5S.transverseX)
          annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,-30},{5.55112e-16,
                5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(C19HF37O5S.face.transverseY, face.C19HF37O5S.transverseY)
          annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,-30},{5.55112e-16,
                5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(C19HF37O5S.face.transverseZ, face.C19HF37O5S.transverseZ)
          annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,-30},{5.55112e-16,
                5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(y.C19HF37O5S, C19HF37O5S.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-40}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // H+
        connect('H+'.face.normal, face.'H+'.normal) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.face.thermal, face.'H+'.thermal) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.face.transverseX, face.'H+'.transverseX) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,-30},{5.55112e-16,
                5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect('H+'.face.transverseY, face.'H+'.transverseY) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,-30},{5.55112e-16,
                5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect('H+'.face.transverseZ, face.'H+'.transverseZ) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,-30},{5.55112e-16,
                5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(y.'H+', 'H+'.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-40}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // H2O
        connect(H2O.face.normal, face.H2O.normal) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.face.thermal, face.H2O.thermal) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.face.transverseX, face.H2O.transverseX) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,-30},{5.55112e-16,
                5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(H2O.face.transverseY, face.H2O.transverseY) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,-30},{5.55112e-16,
                5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(H2O.face.transverseZ, face.H2O.transverseZ) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,-30},{5.55112e-16,
                5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(y.H2O, H2O.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-40}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));
        annotation (Diagram(graphics));
      end Ionomer;

      model Liquid "Sensor for liquid"

        extends BaseClasses.NullPhase;

        // Conditionally include species.
        parameter Boolean inclH2O=false "<html>Water (H<sub>2</sub>O)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        Face.Species H2O(final axis=axis) if inclH2O "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclH2O), Placement(transformation(extent={{-10,-40},{10,-20}})));

      equation
        // H2O
        connect(H2O.face.normal, face.H2O.normal) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.face.thermal, face.H2O.thermal) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.face.transverseX, face.H2O.transverseX) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,-30},{5.55112e-16,
                5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(H2O.face.transverseY, face.H2O.transverseY) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,-30},{5.55112e-16,
                5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(H2O.face.transverseZ, face.H2O.transverseZ) annotation (Line(
            points={{6.10623e-16,-30},{5.55112e-16,-30},{5.55112e-16,
                5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(y.H2O, H2O.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-40}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));
      end Liquid;

      package BaseClasses "Base classes (not for direct use)"
        extends Modelica.Icons.BasesPackage;
        model NullPhase "Empty sensor for a phase (no species)"
          extends FCSys.BaseClasses.Icons.Sensor;

          parameter Axis axis=Axis.x "Axis normal to the face";

          FCSys.Connectors.FaceBus face
            "Multi-species connector for linear momentum and heat" annotation (
              Placement(transformation(extent={{-10,-10},{10,10}}),
                iconTransformation(extent={{-10,-10},{10,10}})));
          FCSys.Connectors.RealOutputBus y "Bus of measurements" annotation (
              Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={0,-100}), iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={0,-100})));

          annotation (defaultComponentName="phaseFaceSensor");
        end NullPhase;
      end BaseClasses;
    end Phases;

    annotation (Documentation(info="<html><p>Since the connectors in
<a href=\"modelica://FCSys\">FCSys</a> are hierarchical
(see the <a href=\"modelica://FCSys.Connectors\">Connectors</a> package),
the models for the sensors must be as well.  A
<a href=\"modelica://FCSys.Connectors.Face\">Face</a>,
<a href=\"modelica://FCSys.Connectors.FaceX\">FaceX</a>,
<a href=\"modelica://FCSys.Connectors.FaceY\">FaceY</a>, or
<a href=\"modelica://FCSys.Connectors.FaceZ\">FaceZ</a> connector
is used in <a href=\"modelica://FCSys.Subregions.Species\">Species</a> models,
and there is a corresponding <a href=\"modelica://FCSys.Sensors.Face.Species\">Species
sensor</a> model. The
<a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a>
connector is used in <a href=\"modelica://FCSys.Subregions.Phases\">Phase</a> models,
and there are corresponding <a href=\"modelica://FCSys.Sensors.Face.Phases\">Phase
sensor</a> models.  The
<a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connector is nested once more
in models such as the <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a>,
and there is a corresponding <a href=\"modelica://FCSys.Sensors.Face.Subregion\">Subregion
sensor</a> model.
</p></html>"));
  end FaceBus;

  package Face
    "<html>Sensors for a single <a href=\"modelica://FCSys.Connectors.Face\">Face</a> connector, e.g., of a <a href=\"modelica://FCSys.Subregions.Species\">Species</a> model (single-species)</html>"
    extends Modelica.Icons.Package;

    model Species
      "<html>Sensor for a face of a <a href=\"modelica://FCSys.Subregions.Species\">Species</a> model (single-species)</html>"

      extends BaseClasses.PartialSpecies;

      parameter Axis axis "Axis normal to the face";

      // X-axis linear momentum
      parameter Boolean inviscidX=true "Inviscid" annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          compact=true,
          group="X-axis linear momentum",
          enable=axis <> 1,
          __Dymola_descriptionLabel=true));
      // Note:  Dymola 7.4 doesn't recognize enumerations in the dialog enable
      // option, e.g.,
      //     enable=axis <> Axis.x.
      // Therefore, the values of the enumerations are specified numerically.
      Velocity velocityX if axis <> Axis.x and not inviscidX "Type of sensor"
        annotation (Placement(transformation(extent={{-40,-10},{-20,10}})));

      // Y-axis linear momentum
      parameter Boolean inviscidY=true "Inviscid" annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          compact=true,
          group="Y-axis linear momentum",
          enable=axis <> 2,
          __Dymola_descriptionLabel=true));
      Velocity velocityY if axis <> Axis.y and not inviscidY "Type of sensor"
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

      // Z-axis linear momentum
      parameter Boolean inviscidZ=true "Inviscid" annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          compact=true,
          group="Z-axis linear momentum",
          enable=axis <> 3,
          __Dymola_descriptionLabel=true));
      Velocity velocityZ if axis <> Axis.z and not inviscidZ "Type of sensor"
        annotation (Placement(transformation(extent={{20,-10},{40,10}})));

      FCSys.Connectors.Face face(
        final isobaric=isobaric,
        final inviscidX=inviscidX,
        final inviscidY=inviscidY,
        final inviscidZ=inviscidZ)
        "Single-species connector for linear momentum and heat" annotation (
          Placement(transformation(extent={{-10,10},{10,30}}),
            iconTransformation(extent={{-10,-10},{10,10}})));

    equation
      // Density
      connect(density.normal, face.normal) annotation (Line(
          points={{-60,6.10623e-16},{-60,20},{5.55112e-16,20}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));

      // X-axis velocity
      connect(velocityX.mechanical, face.transverseX) annotation (Line(
          points={{-30,6.10623e-16},{-30,20},{5.55112e-16,20}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(velocityX.y, y.phi_x) annotation (Line(
          points={{-30,-10},{-30,-40},{0,-40},{0,-100},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%second",
          index=1,
          extent={{-8,-3},{-8,-3}}));

      // Y-axis velocity
      connect(velocityY.mechanical, face.transverseY) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{0,0},{0,20},{5.55112e-16,20}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(velocityY.y, y.phi_y) annotation (Line(
          points={{6.10623e-16,-10},{5.55112e-16,-10},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%second",
          index=1,
          extent={{-8,-3},{-8,-3}}));

      // Z-axis velocity
      connect(velocityZ.mechanical, face.transverseZ) annotation (Line(
          points={{30,6.10623e-16},{30,20},{5.55112e-16,20}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(velocityZ.y, y.phi_z) annotation (Line(
          points={{30,-10},{30,-40},{5.55112e-16,-40},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%second",
          index=1,
          extent={{-8,-3},{-8,-3}}));

      // Temperature
      connect(temperature.thermal, face.thermal) annotation (Line(
          points={{60,6.10623e-16},{60,20},{5.55112e-16,20}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));

      annotation (defaultComponentName="speciesFaceSensor", Diagram(graphics));
    end Species;

    model SpeciesX
      "<html>Sensor for an x-axis face of a <a href=\"modelica://FCSys.Subregions.Species\">Species</a> model (single-species)</html>"

      extends BaseClasses.PartialSpecies;

      // Y-axis linear momentum
      parameter Boolean inviscidY=true "Inviscid" annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          compact=true,
          group="Y-axis linear momentum",
          __Dymola_descriptionLabel=true));
      Velocity velocityY if not inviscidY "Type of sensor"
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

      // Z-axis linear momentum
      parameter Boolean inviscidZ=true "Inviscid" annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          compact=true,
          group="Z-axis linear momentum",
          __Dymola_descriptionLabel=true));
      Velocity velocityZ if not inviscidZ "Type of sensor"
        annotation (Placement(transformation(extent={{20,-10},{40,10}})));

      FCSys.Connectors.FaceX face(
        final isobaric=isobaric,
        final inviscidY=inviscidY,
        final inviscidZ=inviscidZ)
        "Single-species connector for linear momentum and heat" annotation (
          Placement(transformation(extent={{-10,10},{10,30}}),
            iconTransformation(extent={{-10,-10},{10,10}})));

    equation
      // Density
      connect(density.normal, face.normal) annotation (Line(
          points={{-60,6.10623e-16},{-60,20},{5.55112e-16,20}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));

      // Y-axis velocity
      connect(velocityY.mechanical, face.transverseY) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{0,0},{0,20},{5.55112e-16,20}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(velocityY.y, y.phi_y) annotation (Line(
          points={{6.10623e-16,-10},{5.55112e-16,-10},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%second",
          index=1,
          extent={{-8,-3},{-8,-3}}));

      // Z-axis velocity
      connect(velocityZ.mechanical, face.transverseZ) annotation (Line(
          points={{30,6.10623e-16},{30,20},{5.55112e-16,20}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(velocityZ.y, y.phi_z) annotation (Line(
          points={{30,-10},{30,-40},{5.55112e-16,-40},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%second",
          index=1,
          extent={{-8,-3},{-8,-3}}));

      // Temperature
      connect(temperature.thermal, face.thermal) annotation (Line(
          points={{60,6.10623e-16},{60,20},{5.55112e-16,20}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));

      annotation (defaultComponentName="speciesFaceSensor", Diagram(graphics));
    end SpeciesX;

    model SpeciesY
      "<html>Sensor for a y-axis face of a <a href=\"modelica://FCSys.Subregions.Species\">Species</a> model (single-species)</html>"

      extends BaseClasses.PartialSpecies;

      // Z-axis linear momentum
      parameter Boolean inviscidZ=true "Inviscid" annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          compact=true,
          group="Z-axis linear momentum",
          __Dymola_descriptionLabel=true));
      Velocity velocityZ if not inviscidZ "Type of sensor"
        annotation (Placement(transformation(extent={{20,-10},{40,10}})));

      // X-axis linear momentum
      parameter Boolean inviscidX=true "Inviscid" annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          compact=true,
          group="X-axis linear momentum",
          __Dymola_descriptionLabel=true));
      Velocity velocityX if not inviscidX "Type of sensor"
        annotation (Placement(transformation(extent={{-40,-10},{-20,10}})));

      FCSys.Connectors.FaceY face(
        final isobaric=isobaric,
        final inviscidX=inviscidX,
        final inviscidZ=inviscidZ)
        "Single-species connector for linear momentum and heat" annotation (
          Placement(transformation(extent={{-10,10},{10,30}}),
            iconTransformation(extent={{-10,-10},{10,10}})));

    equation
      // Density
      connect(density.normal, face.normal) annotation (Line(
          points={{-60,6.10623e-16},{-60,20},{5.55112e-16,20}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));

      // Z-axis velocity
      connect(velocityZ.mechanical, face.transverseZ) annotation (Line(
          points={{30,6.10623e-16},{30,20},{5.55112e-16,20}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(velocityZ.y, y.phi_z) annotation (Line(
          points={{30,-10},{30,-40},{5.55112e-16,-40},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%second",
          index=1,
          extent={{-8,-3},{-8,-3}}));

      // X-axis velocity
      connect(velocityX.mechanical, face.transverseX) annotation (Line(
          points={{-30,6.10623e-16},{-30,20},{5.55112e-16,20}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(velocityX.y, y.phi_x) annotation (Line(
          points={{-30,-10},{-30,-40},{0,-40},{0,-100},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%second",
          index=1,
          extent={{-8,-3},{-8,-3}}));

      // Temperature
      connect(temperature.thermal, face.thermal) annotation (Line(
          points={{60,6.10623e-16},{60,20},{5.55112e-16,20}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));

      annotation (defaultComponentName="speciesFaceSensor", Diagram(graphics));
    end SpeciesY;

    model SpeciesZ
      "<html>Sensor for a z-axis face of a <a href=\"modelica://FCSys.Subregions.Species\">Species</a> model (single-species)</html>"

      extends BaseClasses.PartialSpecies;

      // X-axis linear momentum
      parameter Boolean inviscidX=true "Inviscid" annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          compact=true,
          group="X-axis linear momentum",
          __Dymola_descriptionLabel=true));
      Velocity velocityX if not inviscidX "Type of sensor"
        annotation (Placement(transformation(extent={{-40,-10},{-20,10}})));

      // Y-axis linear momentum
      parameter Boolean inviscidY=true "Inviscid" annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          compact=true,
          group="Y-axis linear momentum",
          __Dymola_descriptionLabel=true));
      Velocity velocityY if not inviscidY "Type of sensor"
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

      FCSys.Connectors.FaceZ face(
        final isobaric=isobaric,
        final inviscidX=inviscidX,
        final inviscidY=inviscidY)
        "Single-species connector for linear momentum and heat" annotation (
          Placement(transformation(extent={{-10,10},{10,30}}),
            iconTransformation(extent={{-10,-10},{10,10}})));

    equation
      // Density
      connect(density.normal, face.normal) annotation (Line(
          points={{-60,6.10623e-16},{-60,20},{5.55112e-16,20}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));

      // X-axis velocity
      connect(velocityX.mechanical, face.transverseX) annotation (Line(
          points={{-30,6.10623e-16},{-30,20},{5.55112e-16,20}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(velocityX.y, y.phi_x) annotation (Line(
          points={{-30,-10},{-30,-40},{0,-40},{0,-100},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%second",
          index=1,
          extent={{-8,-3},{-8,-3}}));

      // Y-axis velocity
      connect(velocityY.mechanical, face.transverseY) annotation (Line(
          points={{6.10623e-16,6.10623e-16},{0,0},{0,20},{5.55112e-16,20}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(velocityY.y, y.phi_y) annotation (Line(
          points={{6.10623e-16,-10},{5.55112e-16,-10},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None), Text(
          string="%second",
          index=1,
          extent={{-8,-3},{-8,-3}}));

      // Temperature
      connect(temperature.thermal, face.thermal) annotation (Line(
          points={{60,6.10623e-16},{60,20},{5.55112e-16,20}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));

      annotation (defaultComponentName="speciesFaceSensor", Diagram(graphics));
    end SpeciesZ;

    model Density "Sensor for density"

      extends FCSys.Sensors.BaseClasses.PartialSensor(redeclare
          FCSys.Connectors.RealOutput y(final unit="N/l3"));

      FCSys.Connectors.Normal material
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
    equation
      y = material.rho;
      0 = material.Ndot "Conservation of material (no storage)";
      annotation (Icon(graphics={Text(
                  extent={{-100,-20},{100,-50}},
                  lineColor={127,127,127},
                  textString="rho")}));
    end Density;

    model Velocity "Sensor for velocity"
      extends FCSys.Sensors.BaseClasses.PartialSensor(redeclare
          FCSys.Connectors.RealOutput y(final unit="l/T"));
      FCSys.Connectors.Transverse mechanical "Mechanical subconnector"
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

    equation
      y = mechanical.phi "Measurement";
      0 = mechanical.mPhidot "Conservation of linear momentum (no storage)";
      annotation (Icon(graphics={Text(
                  extent={{-100,-20},{100,-50}},
                  lineColor={127,127,127},
                  textString="phi")}));
    end Velocity;

    model Temperature "Sensor for temperature"
      extends FCSys.Sensors.BaseClasses.PartialSensor(redeclare
          FCSys.Connectors.RealOutput y(final unit="l2.m/(N.T2)", displayUnit=
              "K"));
      FCSys.Connectors.Thermal thermal "Thermal subconnector"
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
    equation
      y = thermal.T "Measurement";
      0 = thermal.Qdot "Conservation of energy (no storage)";
      annotation (Icon(graphics={Text(
                  extent={{-100,-20},{100,-50}},
                  lineColor={127,127,127},
                  textString="T")}));
    end Temperature;

    package BaseClasses "Base classes (not for direct use)"
      extends Modelica.Icons.BasesPackage;

      model PartialSpecies
        "<html>Partial sensor for a face of a <a href=\"modelica://FCSys.Subregions.Species\">Species</a> model (single-species)</html>"

        extends FCSys.BaseClasses.Icons.Sensor;

        parameter ThermoOpt isobaric=true
          "Options for material and thermal subconnectors"
          annotation (Dialog(compact=true));

        // Material
        Density density if not isobaric "Type of sensor"
          annotation (Placement(transformation(extent={{-70,-10},{-50,10}})));

        // Heat
        Temperature temperature if thermoOpt <> ThermoOpt.ClosedAdiabatic
          "Type of sensor"
          annotation (Placement(transformation(extent={{50,-10},{70,10}})));

        FCSys.Connectors.RealOutputBus y "Output bus for measurements"
          annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,-100}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,-100})));

      equation
        // Density
        connect(density.y, y.rho) annotation (Line(
            points={{-60,-10},{-60,-40},{0,-40},{0,-100},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%second",
            index=1,
            extent={{-8,-3},{-8,-3}}));

        // Temperature
        connect(temperature.y, y.T) annotation (Line(
            points={{60,-10},{60,-40},{0,-40},{0,-100},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%second",
            index=1,
            extent={{-8,-3},{-8,-3}}));

        annotation (defaultComponentName="speciesFaceSensor", Diagram(graphics));
      end PartialSpecies;
    end BaseClasses;
  end Face;

  package FaceBusDifferential
    "<html>Sensors for a pair of <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connectors</html>"

    extends Modelica.Icons.Package;

    model Subregion
      "<html>Sensor for a face of a <a href=\"modelica://FCSys.Subregions.Subregion\">Region</a> or <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a> model</html>"
      extends FCSys.BaseClasses.Icons.Sensor;

      parameter Axis axis=Axis.x "Axis normal to the face";

      Phases.Gas gas(final axis=axis) "Gas" annotation (Dialog(group="Phases",
            __Dymola_descriptionLabel=true), Placement(transformation(extent={{
                -10,-10},{10,10}})));

      Phases.Graphite graphite(final axis=axis) "Graphite" annotation (Dialog(
            group="Phases", __Dymola_descriptionLabel=true), Placement(
            transformation(extent={{-10,-10},{10,10}})));

      Phases.Ionomer ionomer(final axis=axis) "Ionomer" annotation (Dialog(
            group="Phases", __Dymola_descriptionLabel=true), Placement(
            transformation(extent={{-10,-10},{10,10}})));

      Phases.Liquid liquid(final axis=axis) "Gas" annotation (Dialog(group=
              "Phases", __Dymola_descriptionLabel=true), Placement(
            transformation(extent={{-10,-10},{10,10}})));

      FCSys.Connectors.RealOutputBus y "Bus of measurements" annotation (
          Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,-100})));
      FCSys.Connectors.FaceBus negative
        "Negative-side connector for linear momentum and heat" annotation (
          Placement(transformation(extent={{-110,-10},{-90,10}}),
            iconTransformation(extent={{-110,-10},{-90,10}})));
      FCSys.Connectors.FaceBus positive
        "Positive-side connector for linear momentum and heat" annotation (
          Placement(transformation(extent={{90,-10},{110,10}}),
            iconTransformation(extent={{90,-10},{110,10}})));
    equation
      connect(gas.negative, negative.gas) annotation (Line(
          points={{-10,6.10623e-16},{-10,5.55112e-16},{-100,5.55112e-16}},
          color={127,127,127},
          pattern=LinePattern.None,
          thickness=0.5,
          smooth=Smooth.None));

      connect(gas.positive, positive.gas) annotation (Line(
          points={{10,6.10623e-16},{10,5.55112e-16},{100,5.55112e-16}},
          color={127,127,127},
          pattern=LinePattern.None,
          thickness=0.5,
          smooth=Smooth.None));

      connect(graphite.negative, negative.graphite) annotation (Line(
          points={{-10,6.10623e-16},{-10,5.55112e-16},{-100,5.55112e-16}},
          color={127,127,127},
          pattern=LinePattern.None,
          thickness=0.5,
          smooth=Smooth.None));

      connect(graphite.positive, positive.graphite) annotation (Line(
          points={{10,6.10623e-16},{10,5.55112e-16},{100,5.55112e-16}},
          color={127,127,127},
          pattern=LinePattern.None,
          thickness=0.5,
          smooth=Smooth.None));

      connect(ionomer.negative, negative.ionomer) annotation (Line(
          points={{-10,6.10623e-16},{-10,5.55112e-16},{-100,5.55112e-16}},
          color={127,127,127},
          pattern=LinePattern.None,
          thickness=0.5,
          smooth=Smooth.None));

      connect(ionomer.positive, positive.ionomer) annotation (Line(
          points={{10,6.10623e-16},{10,5.55112e-16},{100,5.55112e-16}},
          color={127,127,127},
          pattern=LinePattern.None,
          thickness=0.5,
          smooth=Smooth.None));

      connect(y.gas, gas.y) annotation (Line(
          points={{5.55112e-16,-100},{6.10623e-16,-10}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None));

      connect(y.graphite, graphite.y) annotation (Line(
          points={{5.55112e-16,-100},{6.10623e-16,-10}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None));

      connect(y.ionomer, ionomer.y) annotation (Line(
          points={{5.55112e-16,-100},{6.10623e-16,-10}},
          color={0,0,127},
          thickness=0.5,
          smooth=Smooth.None));

      annotation (
        defaultComponentPrefixes="replaceable",
        defaultComponentName="subregionFaceSensor",
        Icon(graphics={Line(
                  points={{-70,0},{-100,0}},
                  color={127,127,127},
                  smooth=Smooth.None),Line(
                  points={{100,0},{70,0}},
                  color={127,127,127},
                  smooth=Smooth.None)}),
        Diagram(graphics));
    end Subregion;

    package Phases
      "<html>Sensors for the <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connector, e.g., of a <a href=\"modelica://FCSys.Subregions.Phase\">Phase</a> model (multi-species)</html>"
      extends Modelica.Icons.Package;

      model Gas "Sensor for gas"

        extends BaseClasses.NullPhase;

        // Conditionally include species.
        parameter Boolean inclH2=false "<html>Hydrogen (H<sub>2</sub>)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        Face.Species H2(final axis=axis) if inclH2 "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclH2), Placement(transformation(extent={{-10,-10},{10,10}})));
        parameter Boolean inclH2O=false "<html>Water (H<sub>2</sub>O)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        Face.Species H2O(final axis=axis) if inclH2O "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclH2O), Placement(transformation(extent={{-10,-10},{10,10}})));
        parameter Boolean inclN2=false "<html>Nitrogen (N<sub>2</sub>)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        Face.Species N2(final axis=axis) if inclN2 "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclN2), Placement(transformation(extent={{-10,-10},{10,10}})));
        parameter Boolean inclO2=false "<html>Oxygen (O<sub>2</sub>)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        Face.Species O2(final axis=axis) if inclO2 "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclO2), Placement(transformation(extent={{-10,-10},{10,10}})));

      equation
        // H2
        connect(H2.negative, negative.H2) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.positive, positive.H2) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(y.H2, H2.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-10}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // H2O
        connect(H2O.negative, negative.H2O) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive, positive.H2O) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(y.H2O, H2O.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-10}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // N2
        connect(N2.negative, negative.N2) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.positive, positive.N2) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(y.N2, N2.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-10}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // O2
        connect(O2.negative, negative.O2) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.positive, positive.O2) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(y.O2, O2.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-10}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));
        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="phaseFaceSensor");
      end Gas;

      model Graphite "Sensor for graphite"

        extends BaseClasses.NullPhase;

        // Conditionally include species.
        parameter Boolean inclC=false "Carbon (C)" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        Face.Species C(final axis=axis) if inclC "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclC), Placement(transformation(extent={{-10,-10},{10,10}})));
        parameter Boolean 'incle-'=false
          "<html>Electrons (e<sup>-</sup>)</html>" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        Face.Species 'e-'(final axis=axis) if 'incle-' "Model" annotation (
            Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable='incle-'), Placement(transformation(extent={{-10,-10},{10,10}})));

      equation
        // C
        connect(C.negative, negative.C) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.positive, positive.C) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(y.C, C.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-10}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // e-
        connect('e-'.negative, negative.'e-') annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.positive, positive.'e-') annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(y.'e-', 'e-'.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-10}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="phaseFaceSensor");
      end Graphite;

      model Ionomer "Sensor for ionomer"

        extends BaseClasses.NullPhase;

        // Conditionally include species.
        parameter Boolean inclC19HF37O5S=false
          "<html>Nafion sulfonate (C<sub>19</sub>HF<sub>37</sub>O<sub>5</sub>S)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FaceDifferential.Species C19HF37O5S(final axis=axis) if inclC19HF37O5S
          "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclC19HF37O5S), Placement(transformation(extent={{-10,-10},
                  {10,10}})));
        parameter Boolean 'inclH+'=false "<html>Protons (H<sup>+</sup>)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FaceDifferential.Species 'H+'(final axis=axis) if 'inclH+' "Model"
          annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable='inclH+'), Placement(transformation(extent={{-10,-10},{10,10}})));
        parameter Boolean inclH2O=false "<html>Water (H<sub>2</sub>O)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FaceDifferential.Species H2O(final axis=axis) if inclH2O "Model"
          annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclH2O), Placement(transformation(extent={{-10,-10},{10,10}})));

      equation
        connect(C19HF37O5S.negative, negative.C19HF37O5S) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.positive, positive.C19HF37O5S) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(y.C19HF37O5S, C19HF37O5S.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-9.8}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // H+
        connect('H+'.negative, negative.'H+') annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.positive, positive.'H+') annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(y.'H+', 'H+'.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-9.8}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // H2O
        connect(H2O.negative, negative.H2O) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive, positive.H2O) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(y.H2O, H2O.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-9.8}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="phaseFaceSensor",
          Diagram(graphics));
      end Ionomer;

      model Liquid "Sensor for liquid"

        extends FaceBusDifferential.Phases.BaseClasses.NullPhase;

        // Conditionally include species.
        parameter Boolean inclH2O=false "<html>Water (H<sub>2</sub>O)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FaceDifferential.Species H2O(final axis=axis) if inclH2O "Model"
          annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclH2O), Placement(transformation(extent={{-10,-10},{10,10}})));

      equation
        // H2O
        connect(H2O.negative, negative.H2O) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive, positive.H2O) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(y.H2O, H2O.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-9.8}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="phaseFaceSensor",
          Diagram(graphics));
      end Liquid;

      package BaseClasses "Base classes (not for direct use)"
        extends Modelica.Icons.BasesPackage;

        model NullPhase "Empty sensor for a phase (no species)"
          extends FCSys.BaseClasses.Icons.Sensor;

          parameter Axis axis=Axis.x "Axis normal to the face";

          FCSys.Connectors.FaceBus negative
            "Negative-side connector for linear momentum and heat" annotation (
              Placement(transformation(extent={{-110,-10},{-90,10}}),
                iconTransformation(extent={{-110,-10},{-90,10}})));
          FCSys.Connectors.FaceBus positive
            "Positive-side connector for linear momentum and heat" annotation (
              Placement(transformation(extent={{90,-10},{110,10}}),
                iconTransformation(extent={{90,-10},{110,10}})));
          FCSys.Connectors.RealOutputBus y "Bus of measurements" annotation (
              Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={0,-100}), iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={0,-100})));

          annotation (defaultComponentName="phaseFaceSensor", Icon(graphics={
                  Line(   points={{-70,0},{-100,0}},
                          color={127,127,127},
                          smooth=Smooth.None),Line(
                          points={{100,0},{70,0}},
                          color={127,127,127},
                          smooth=Smooth.None)}));
        end NullPhase;
      end BaseClasses;
    end Phases;

    annotation (Documentation(info="<html><p>Since the connectors in
<a href=\"modelica://FCSys\">FCSys</a> are hierarchical
(see the <a href=\"modelica://FCSys.Connectors\">Connectors</a> package),
the models for the sensors must be as well.   A
<a href=\"modelica://FCSys.Connectors.Face\">Face</a>,
<a href=\"modelica://FCSys.Connectors.FaceX\">FaceX</a>,
<a href=\"modelica://FCSys.Connectors.FaceY\">FaceY</a>, or
<a href=\"modelica://FCSys.Connectors.FaceZ\">FaceZ</a>
connector
is used in <a href=\"modelica://FCSys.Subregions.Species\">Species</a> models,
and there is a corresponding <a href=\"modelica://FCSys.Sensors.FaceDifferential.Species\">Species
sensor</a> model. The
<a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a>
connector is used in <a href=\"modelica://FCSys.Subregions.Phases\">Phase</a> models,
and there are corresponding <a href=\"modelica://FCSys.Sensors.Face.Phases\">Phase
sensor</a> models.  The
<a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connector is nested once more
in models such as the <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a>,
and there is a corresponding <a href=\"modelica://FCSys.Sensors.Face.Subregion\">Subregion
sensor</a> model.
</p></html>"));
  end FaceBusDifferential;

  package FaceDifferential
    "<html>Sensors for a pair of <a href=\"modelica://FCSys.Connectors.Face\">Face</a> connectors, e.g., of a <a href=\"modelica://FCSys.Subregions.Species\">Species</a> model (single-species)</html>"
    extends Modelica.Icons.Package;

    model Species
      "<html>Sensor for faces of a <a href=\"modelica://FCSys.Subregions.Species\">Species</a> model (single-species)</html>"

      extends BaseClasses.PartialSpecies;

      parameter Axis axis=Axis.x "Axis normal to the face";

      // X-axis mechanical
      parameter Boolean inviscidX=true "Inviscid" annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          compact=true,
          group="X-axis linear momentum",
          enable=axis <> 1,
          __Dymola_descriptionLabel=true));
      // Note:  Dymola 7.4 doesn't recognize enumerations in the dialog enable
      // option, e.g.,
      //     enable=axis <> Axis.x.
      // Therefore, the values of the enumerations are specified numerically.
      replaceable Mechanical.Velocity transverseX if axis <> Axis.x and not
        inviscidX "Condition" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="X-axis linear momentum",
          enable=axis <> 1 and not inviscidX,
          __Dymola_descriptionLabel=true),
        Placement(transformation(extent={{-40,-10},{-20,10}})));

      // Y-axis mechanical
      parameter Boolean inviscidY=true "Inviscid" annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          compact=true,
          group="Y-axis linear momentum",
          enable=axis <> 2,
          __Dymola_descriptionLabel=true));
      replaceable Mechanical.Velocity transverseY if axis <> Axis.y and not
        inviscidY "Condition" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Y-axis linear momentum",
          enable=axis <> 2 and not inviscidY,
          __Dymola_descriptionLabel=true),
        Placement(transformation(extent={{-10,-10},{10,10}})));

      // Z-axis mechanical
      parameter Boolean inviscidZ=true "Inviscid" annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          compact=true,
          group="Z-axis linear momentum",
          enable=axis <> 3,
          __Dymola_descriptionLabel=true));
      replaceable Mechanical.Velocity transverseZ if axis <> Axis.z and not
        inviscidZ "Condition" annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Z-axis linear momentum",
          enable=axis <> 3 and not inviscidZ,
          __Dymola_descriptionLabel=true),
        Placement(transformation(extent={{20,-10},{40,10}})));

      FCSys.Connectors.Face negative(
        final isobaric=isobaric,
        final inviscidX=inviscidX,
        final inviscidY=inviscidY,
        final inviscidZ=inviscidZ)
        "Negative-side connector for linear momentum and heat" annotation (
          Placement(transformation(extent={{-110,-10},{-90,10}}),
            iconTransformation(extent={{-110,-10},{-90,10}})));

      FCSys.Connectors.Face positive(
        final isobaric=isobaric,
        final inviscidX=inviscidX,
        final inviscidY=inviscidY,
        final inviscidZ=inviscidZ)
        "Positive-side connector for linear momentum and heat" annotation (
          Placement(transformation(extent={{90,-2},{110,18}}),
            iconTransformation(extent={{90,-10},{110,10}})));

    equation
      // Material
      connect(material.negative, negative.normal) annotation (Line(
          points={{-70,6.10623e-16},{-70,5.55112e-16},{-100,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(material.positive, positive.normal) annotation (Line(
          points={{-50,6.10623e-16},{90,6.10623e-16},{90,8},{100,8}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));

      // X-axis mechanical
      connect(transverseX.negative, negative.transverseX) annotation (Line(
          points={{-40,6.10623e-16},{-90,6.10623e-16},{-90,5.55112e-16},{-100,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(transverseX.positive, positive.transverseX) annotation (Line(
          points={{-20,6.10623e-16},{90,6.10623e-16},{90,8},{100,8}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(transverseX.y, y.transverseX) annotation (Line(
          points={{-30,-10},{-30,-20},{5.55112e-16,-20},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      // Y-axis mechanical
      connect(transverseY.negative, negative.transverseY) annotation (Line(
          points={{-10,6.10623e-16},{-90,6.10623e-16},{-90,5.55112e-16},{-100,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(transverseY.positive, positive.transverseY) annotation (Line(
          points={{10,6.10623e-16},{90,6.10623e-16},{90,8},{100,8}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(transverseY.y, y.transverseY) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-20},{5.55112e-16,-20},{
              5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      // Z-axis mechanical
      connect(transverseZ.negative, negative.transverseZ) annotation (Line(
          points={{20,6.10623e-16},{-90,6.10623e-16},{-90,5.55112e-16},{-100,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(transverseZ.positive, positive.transverseZ) annotation (Line(
          points={{40,6.10623e-16},{90,6.10623e-16},{90,8},{100,8}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(transverseZ.y, y.transverseZ) annotation (Line(
          points={{30,-10},{30,-20},{5.55112e-16,-20},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      // Thermal
      connect(thermal.negative, negative.thermal) annotation (Line(
          points={{50,6.10623e-16},{-90,6.10623e-16},{-90,5.55112e-16},{-100,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(thermal.positive, positive.thermal) annotation (Line(
          points={{70,6.10623e-16},{90,6.10623e-16},{90,8},{100,8}},
          color={127,127,127},
          smooth=Smooth.None));

      annotation (defaultComponentName="speciesFaceSensor", Icon(graphics={Line(
                  points={{-70,0},{-100,0}},
                  color={0,0,0},
                  smooth=Smooth.None),Line(
                  points={{100,0},{70,0}},
                  color={0,0,0},
                  smooth=Smooth.None)}));
    end Species;

    model SpeciesX
      "<html>Sensor for x-axis faces of a <a href=\"modelica://FCSys.Subregions.Species\">Species</a> model (single-species)</html>"

      extends BaseClasses.PartialSpecies;

      // Y-axis mechanical
      parameter Boolean inviscidY=true "Inviscid" annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          compact=true,
          group="Y-axis linear momentum",
          __Dymola_descriptionLabel=true));
      replaceable Mechanical.Velocity transverseY if not inviscidY "Condition"
        annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Y-axis linear momentum",
          enable=not inviscidY,
          __Dymola_descriptionLabel=true),
        Placement(transformation(extent={{-10,-10},{10,10}})));

      // Z-axis mechanical
      parameter Boolean inviscidZ=true "Inviscid" annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          compact=true,
          group="Z-axis linear momentum",
          __Dymola_descriptionLabel=true));
      replaceable Mechanical.Velocity transverseZ if not inviscidZ "Condition"
        annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Z-axis linear momentum",
          enable=not inviscidZ,
          __Dymola_descriptionLabel=true),
        Placement(transformation(extent={{20,-10},{40,10}})));

      FCSys.Connectors.FaceX negative(
        final isobaric=isobaric,
        final inviscidY=inviscidY,
        final inviscidZ=inviscidZ)
        "Negative-side connector for linear momentum and heat" annotation (
          Placement(transformation(extent={{-110,-10},{-90,10}}),
            iconTransformation(extent={{-110,-10},{-90,10}})));

      FCSys.Connectors.FaceX positive(
        final isobaric=isobaric,
        final inviscidY=inviscidY,
        final inviscidZ=inviscidZ)
        "Positive-side connector for linear momentum and heat" annotation (
          Placement(transformation(extent={{90,-2},{110,18}}),
            iconTransformation(extent={{90,-10},{110,10}})));

    equation
      // Material
      connect(material.negative, negative.normal) annotation (Line(
          points={{-70,6.10623e-16},{-70,5.55112e-16},{-100,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(material.positive, positive.normal) annotation (Line(
          points={{-50,6.10623e-16},{90,6.10623e-16},{90,8},{100,8}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));

      // Y-axis mechanical
      connect(transverseY.negative, negative.transverseY) annotation (Line(
          points={{-10,6.10623e-16},{-90,6.10623e-16},{-90,5.55112e-16},{-100,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(transverseY.positive, positive.transverseY) annotation (Line(
          points={{10,6.10623e-16},{90,6.10623e-16},{90,8},{100,8}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(transverseY.y, y.transverseY) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-20},{5.55112e-16,-20},{
              5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      // Z-axis mechanical
      connect(transverseZ.negative, negative.transverseZ) annotation (Line(
          points={{20,6.10623e-16},{-90,6.10623e-16},{-90,5.55112e-16},{-100,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(transverseZ.positive, positive.transverseZ) annotation (Line(
          points={{40,6.10623e-16},{90,6.10623e-16},{90,8},{100,8}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(transverseZ.y, y.transverseZ) annotation (Line(
          points={{30,-10},{30,-20},{5.55112e-16,-20},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      // Thermal
      connect(thermal.negative, negative.thermal) annotation (Line(
          points={{50,6.10623e-16},{-90,6.10623e-16},{-90,5.55112e-16},{-100,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(thermal.positive, positive.thermal) annotation (Line(
          points={{70,6.10623e-16},{90,6.10623e-16},{90,8},{100,8}},
          color={127,127,127},
          smooth=Smooth.None));

      annotation (defaultComponentName="speciesFaceSensor", Icon(graphics={Line(
                  points={{-70,0},{-100,0}},
                  color={0,0,0},
                  smooth=Smooth.None),Line(
                  points={{100,0},{70,0}},
                  color={0,0,0},
                  smooth=Smooth.None)}));
    end SpeciesX;

    model SpeciesY
      "<html>Sensor for y-axis faces of a <a href=\"modelica://FCSys.Subregions.Species\">Species</a> model (single-species)</html>"

      extends BaseClasses.PartialSpecies;

      // Z-axis mechanical
      parameter Boolean inviscidZ=true "Inviscid" annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          compact=true,
          group="Z-axis linear momentum",
          __Dymola_descriptionLabel=true));
      replaceable Mechanical.Velocity transverseZ if not inviscidZ "Condition"
        annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Z-axis linear momentum",
          enable=not inviscidZ,
          __Dymola_descriptionLabel=true),
        Placement(transformation(extent={{20,-10},{40,10}})));

      // X-axis mechanical
      parameter Boolean inviscidX=true "Inviscid" annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          compact=true,
          group="X-axis linear momentum",
          __Dymola_descriptionLabel=true));
      replaceable Mechanical.Velocity transverseX if not inviscidX "Condition"
        annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="X-axis linear momentum",
          enable=not inviscidX,
          __Dymola_descriptionLabel=true),
        Placement(transformation(extent={{-40,-10},{-20,10}})));

      FCSys.Connectors.FaceY negative(
        final isobaric=isobaric,
        final inviscidX=inviscidX,
        final inviscidZ=inviscidZ)
        "Negative-side connector for linear momentum and heat" annotation (
          Placement(transformation(extent={{-110,-10},{-90,10}}),
            iconTransformation(extent={{-110,-10},{-90,10}})));

      FCSys.Connectors.FaceY positive(
        final isobaric=isobaric,
        final inviscidX=inviscidX,
        final inviscidZ=inviscidZ)
        "Positive-side connector for linear momentum and heat" annotation (
          Placement(transformation(extent={{90,-2},{110,18}}),
            iconTransformation(extent={{90,-10},{110,10}})));

    equation
      // Material
      connect(material.negative, negative.normal) annotation (Line(
          points={{-70,6.10623e-16},{-70,5.55112e-16},{-100,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(material.positive, positive.normal) annotation (Line(
          points={{-50,6.10623e-16},{90,6.10623e-16},{90,8},{100,8}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));

      // Z-axis mechanical
      connect(transverseZ.negative, negative.transverseZ) annotation (Line(
          points={{20,6.10623e-16},{-90,6.10623e-16},{-90,5.55112e-16},{-100,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(transverseZ.positive, positive.transverseZ) annotation (Line(
          points={{40,6.10623e-16},{90,6.10623e-16},{90,8},{100,8}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(transverseZ.y, y.transverseZ) annotation (Line(
          points={{30,-10},{30,-20},{5.55112e-16,-20},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      // X-axis mechanical
      connect(transverseX.negative, negative.transverseX) annotation (Line(
          points={{-40,6.10623e-16},{-90,6.10623e-16},{-90,5.55112e-16},{-100,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(transverseX.positive, positive.transverseX) annotation (Line(
          points={{-20,6.10623e-16},{90,6.10623e-16},{90,8},{100,8}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(transverseX.y, y.transverseX) annotation (Line(
          points={{-30,-10},{-30,-20},{5.55112e-16,-20},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      // Thermal
      connect(thermal.negative, negative.thermal) annotation (Line(
          points={{50,6.10623e-16},{-90,6.10623e-16},{-90,5.55112e-16},{-100,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(thermal.positive, positive.thermal) annotation (Line(
          points={{70,6.10623e-16},{90,6.10623e-16},{90,8},{100,8}},
          color={127,127,127},
          smooth=Smooth.None));

      annotation (defaultComponentName="speciesFaceSensor", Icon(graphics={Line(
                  points={{-70,0},{-100,0}},
                  color={0,0,0},
                  smooth=Smooth.None),Line(
                  points={{100,0},{70,0}},
                  color={0,0,0},
                  smooth=Smooth.None)}));
    end SpeciesY;

    model SpeciesZ
      "<html>Sensor for z-axis faces of a <a href=\"modelica://FCSys.Subregions.Species\">Species</a> model (single-species)</html>"

      extends BaseClasses.PartialSpecies;

      // X-axis mechanical
      parameter Boolean inviscidX=true "Inviscid" annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          compact=true,
          group="X-axis linear momentum",
          __Dymola_descriptionLabel=true));
      replaceable Mechanical.Velocity transverseX if not inviscidX "Condition"
        annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="X-axis linear momentum",
          enable=not inviscidX,
          __Dymola_descriptionLabel=true),
        Placement(transformation(extent={{-40,-10},{-20,10}})));

      // Y-axis mechanical
      parameter Boolean inviscidY=true "Inviscid" annotation (
        HideResult=true,
        choices(__Dymola_checkBox=true),
        Dialog(
          compact=true,
          group="Y-axis linear momentum",
          __Dymola_descriptionLabel=true));
      replaceable Mechanical.Velocity transverseY if not inviscidY "Condition"
        annotation (
        __Dymola_choicesFromPackage=true,
        Dialog(
          group="Y-axis linear momentum",
          enable=not inviscidY,
          __Dymola_descriptionLabel=true),
        Placement(transformation(extent={{-10,-10},{10,10}})));

      FCSys.Connectors.FaceZ negative(
        final isobaric=isobaric,
        final inviscidX=inviscidX,
        final inviscidY=inviscidY)
        "Negative-side connector for linear momentum and heat" annotation (
          Placement(transformation(extent={{-110,-10},{-90,10}}),
            iconTransformation(extent={{-110,-10},{-90,10}})));

      FCSys.Connectors.FaceZ positive(
        final isobaric=isobaric,
        final inviscidX=inviscidX,
        final inviscidY=inviscidY)
        "Positive-side connector for linear momentum and heat" annotation (
          Placement(transformation(extent={{90,-2},{110,18}}),
            iconTransformation(extent={{90,-10},{110,10}})));

    equation
      // Material
      connect(material.negative, negative.normal) annotation (Line(
          points={{-70,6.10623e-16},{-70,5.55112e-16},{-100,5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));
      connect(material.positive, positive.normal) annotation (Line(
          points={{-50,6.10623e-16},{90,6.10623e-16},{90,8},{100,8}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));

      // X-axis mechanical
      connect(transverseX.negative, negative.transverseX) annotation (Line(
          points={{-40,6.10623e-16},{-90,6.10623e-16},{-90,5.55112e-16},{-100,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(transverseX.positive, positive.transverseX) annotation (Line(
          points={{-20,6.10623e-16},{90,6.10623e-16},{90,8},{100,8}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(transverseX.y, y.transverseX) annotation (Line(
          points={{-30,-10},{-30,-20},{5.55112e-16,-20},{5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      // Y-axis mechanical
      connect(transverseY.negative, negative.transverseY) annotation (Line(
          points={{-10,6.10623e-16},{-90,6.10623e-16},{-90,5.55112e-16},{-100,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(transverseY.positive, positive.transverseY) annotation (Line(
          points={{10,6.10623e-16},{90,6.10623e-16},{90,8},{100,8}},
          color={127,127,127},
          pattern=LinePattern.None,
          smooth=Smooth.None));
      connect(transverseY.y, y.transverseY) annotation (Line(
          points={{6.10623e-16,-10},{6.10623e-16,-20},{5.55112e-16,-20},{
              5.55112e-16,-100}},
          color={0,0,127},
          smooth=Smooth.None));

      // Thermal
      connect(thermal.negative, negative.thermal) annotation (Line(
          points={{50,6.10623e-16},{-90,6.10623e-16},{-90,5.55112e-16},{-100,
              5.55112e-16}},
          color={127,127,127},
          smooth=Smooth.None));

      connect(thermal.positive, positive.thermal) annotation (Line(
          points={{70,6.10623e-16},{90,6.10623e-16},{90,8},{100,8}},
          color={127,127,127},
          smooth=Smooth.None));

      annotation (defaultComponentName="speciesFaceSensor", Icon(graphics={Line(
                  points={{-70,0},{-100,0}},
                  color={0,0,0},
                  smooth=Smooth.None),Line(
                  points={{100,0},{70,0}},
                  color={0,0,0},
                  smooth=Smooth.None)}));
    end SpeciesZ;

    package Material "Relative sensors for normal linear momentum"
      extends Modelica.Icons.Package;
      model Density "Sensor for density difference (closed condition)"

        extends BaseClasses.PartialSensor(final sensorType=BaseClasses.SensorType.Density,
            redeclare FCSys.Connectors.RealOutput y(final unit="N/l3"));

      equation
        y = negative.rho - positive.rho "Measurement";
        0 = negative.Ndot "Condition of no current";
        // Note:  In conjunction with the material conservation equation, this
        // means that there's no current into either face.
        annotation (Icon(graphics={Text(
                      extent={{-100,-20},{100,-50}},
                      lineColor={127,127,127},
                      textString="  rho"),Polygon(
                      points={{-26,-26},{-18,-44},{-34,-44},{-26,-26}},
                      lineColor={127,127,127},
                      smooth=Smooth.None,
                      lineThickness=0.5)}));
      end Density;

      model Current "Sensor for current (isochoric condition)"
        extends BaseClasses.PartialSensor(final sensorType=BaseClasses.SensorType.Current,
            redeclare FCSys.Connectors.RealOutput y(final unit="N/T"));

      equation
        y = negative.Ndot "Measurement";
        negative.rho = positive.rho
          "Condition of equal electrochemical potentials";

        annotation (Icon(graphics={Text(
                      extent={{-100,-20},{100,-50}},
                      lineColor={127,127,127},
                      textString="I")}));
      end Current;

      package BaseClasses "Base classes (not for direct use)"
        extends Modelica.Icons.BasesPackage;

        partial model PartialSensor
          "<html>Partial sensor for a pair of <a href=\"modelica://FCSys.Connectors.Normal\">Normal</a> connectors</html>"

          extends FCSys.Sensors.BaseClasses.PartialSensor;

          constant SensorType sensorType "Type of sensor";
          // Note:  This is included so that the type of BC is recorded with the
          // results.

          FCSys.Connectors.Normal negative
            "Material connector for the negative face" annotation (Placement(
                transformation(extent={{-110,-10},{-90,10}}),
                iconTransformation(extent={{-110,-10},{-90,10}})));
          FCSys.Connectors.Normal positive
            "Material connector for the positive face" annotation (Placement(
                transformation(extent={{90,-10},{110,10}}), iconTransformation(
                  extent={{90,-10},{110,10}})));

        equation
          0 = negative.Ndot + positive.Ndot
            "Material rate balance (no storage)";

          annotation (Icon(graphics={Line(
                          points={{-70,0},{-100,0}},
                          color={0,0,0},
                          smooth=Smooth.None),Line(
                          points={{100,0},{70,0}},
                          color={0,0,0},
                          smooth=Smooth.None)}));
        end PartialSensor;

        type SensorType = enumeration(
            Density "Density difference",
            Current "Current") "Types of sensors";
      end BaseClasses;
    end Material;

    package Mechanical "Relative sensors for transverse linear momentum"
      extends Modelica.Icons.Package;
      model Velocity "Sensor for relative velocity (inviscid condition)"

        extends BaseClasses.PartialSensor(final sensorType=BaseClasses.SensorType.Velocity,
            redeclare FCSys.Connectors.RealOutput y(final unit="l/T"));

      equation
        y = negative.phi - positive.phi "Measurement";
        0 = negative.mPhidot "Condition of no shear force";
        // Note:  In conjunction with the conservation of linear momentum, this
        // means that there's no force on either face.
        annotation (
          Documentation(info="<html><p>If <code>ax=1</code>, then the difference in normal velocity
  is measured.  The relative first
  and second transverse components of linear momentum may be measured by
  <code>ax=2</code> and <code>ax=3</code>, but only if <code>inclTrans</code>
   is <code>true</code>.
  </p></html>"),
          Icon(graphics={Text(
                      extent={{-100,-20},{100,-50}},
                      lineColor={127,127,127},
                      textString="  phi"),Polygon(
                      points={{-20,-26},{-12,-44},{-28,-44},{-20,-26}},
                      lineColor={127,127,127},
                      smooth=Smooth.None,
                      lineThickness=0.5)}),
          Diagram(graphics));
      end Velocity;

      model Force "Sensor for shear force (condition of uniform velocity)"
        extends BaseClasses.PartialSensor(final sensorType=BaseClasses.SensorType.Force,
            redeclare FCSys.Connectors.RealOutput y(final unit="l.m/T2"));

      equation
        y = negative.mPhidot "Measurement";
        negative.phi = positive.phi "Condition of uniform velocities";
        annotation (Documentation(info="<html><p>If <code>ax=1</code>, then normal force
  is measured.  The first
  and second transverse components of force may be measured by
  <code>ax=2</code> and <code>ax=3</code>, but only if <code>inclTrans</code>
   is <code>true</code>.
  </p></html>"), Icon(graphics={Text(
                      extent={{-100,-20},{100,-50}},
                      lineColor={127,127,127},
                      textString="mPhi"),Text(
                      extent={{-100,8},{100,-22}},
                      lineColor={127,127,127},
                      textString=".")}));
      end Force;

      package BaseClasses "Base classes (not for direct use)"
        extends Modelica.Icons.BasesPackage;

        partial model PartialSensor
          "<html>Partial sensor for a pair of <a href=\"modelica://FCSys.Connectors.Transverse\">Transverse</a> connectors</html>"

          extends FCSys.Sensors.BaseClasses.PartialSensor;

          constant SensorType sensorType "Type of sensor";
          // Note:  This is included so that the type of BC is recorded with the
          // results.

          FCSys.Connectors.Transverse negative
            "Mechanical connector for the negative face" annotation (Placement(
                transformation(extent={{-110,-10},{-90,10}}),
                iconTransformation(extent={{-110,-10},{-90,10}})));
          FCSys.Connectors.Transverse positive
            "Mechanical connector for the positive face" annotation (Placement(
                transformation(extent={{90,-10},{110,10}}), iconTransformation(
                  extent={{90,-10},{110,10}})));
        equation
          0 = negative.mPhidot + positive.mPhidot
            "Conservation of linear momentum (no storage)";

          annotation (Icon(graphics={Line(
                          points={{-70,0},{-100,0}},
                          color={0,0,0},
                          smooth=Smooth.None),Line(
                          points={{100,0},{70,0}},
                          color={0,0,0},
                          smooth=Smooth.None)}));
        end PartialSensor;

        type SensorType = enumeration(
            Velocity "Shear velocity",
            Force "Shear force") "Types of sensors";
      end BaseClasses;
    end Mechanical;

    package Thermal "Relative thermal sensors"
      extends Modelica.Icons.Package;
      model Temperature
        "Sensor for temperature difference (adiabatic condition)"

        extends BaseClasses.PartialSensor(final sensorType=BaseClasses.SensorType.Temperature,
            redeclare FCSys.Connectors.RealOutput y(final unit="l2.m/(N.T2)",
              displayUnit="K"));

      equation
        y = negative.T - positive.T "Measurement";
        0 = negative.Qdot "Adiabatic condition";
        // Note:  In conjunction with the energy rate balance, this means
        // that both faces are adiabatic.
        annotation (Icon(graphics={Text(
                      extent={{-100,-20},{100,-50}},
                      lineColor={127,127,127},
                      textString="  T"),Polygon(
                      points={{-10,-26},{-2,-44},{-18,-44},{-10,-26}},
                      lineColor={127,127,127},
                      smooth=Smooth.None,
                      lineThickness=0.5)}));
      end Temperature;

      model HeatRate "Sensor for heat flow rate (isothermal condition)"
        extends BaseClasses.PartialSensor(final sensorType=BaseClasses.SensorType.HeatRate,
            redeclare FCSys.Connectors.RealOutput y(final unit="l2.m/T3"));

      equation
        y = negative.Qdot "Measurement";
        negative.T = positive.T "Isothermal condition";
        annotation (Icon(graphics={Text(
                      extent={{-100,-20},{100,-50}},
                      lineColor={127,127,127},
                      textString="Q"),Text(
                      extent={{-100,8},{100,-22}},
                      lineColor={127,127,127},
                      textString=".")}));
      end HeatRate;

      package BaseClasses "Base classes (not for direct use)"
        extends Modelica.Icons.BasesPackage;

        partial model PartialSensor
          "<html>Partial sensor for a pair of <a href=\"modelica://FCSys.Connectors.Thermal\">Thermal</a> connectors</html>"

          extends FCSys.Sensors.BaseClasses.PartialSensor;

          constant SensorType sensorType "Type of sensor";
          // Note:  This is included so that the type of BC is recorded with the
          // results.

          FCSys.Connectors.Thermal negative
            "Thermal connector for the negative face" annotation (Placement(
                transformation(extent={{-110,-10},{-90,10}}),
                iconTransformation(extent={{-110,-10},{-90,10}})));
          FCSys.Connectors.Thermal positive
            "Thermal connector for the positive face" annotation (Placement(
                transformation(extent={{90,-10},{110,10}}), iconTransformation(
                  extent={{90,-10},{110,10}})));
        equation
          0 = negative.Qdot + positive.Qdot
            "Conservation of energy (no storage)";

          annotation (Icon(graphics={Line(
                          points={{-70,0},{-100,0}},
                          color={0,0,0},
                          smooth=Smooth.None),Line(
                          points={{100,0},{70,0}},
                          color={0,0,0},
                          smooth=Smooth.None)}));
        end PartialSensor;

        type SensorType = enumeration(
            Temperature "Temperature difference",
            HeatRate "Heat flow rate") "Types of sensors";
      end BaseClasses;
    end Thermal;

    package BaseClasses "Base classes (not for direct use)"
      extends Modelica.Icons.BasesPackage;

      partial model PartialSpecies
        "<html>Partial sensor for faces of a <a href=\"modelica://FCSys.Subregions.Species\">Species</a> model (single-species)</html>"

        extends FCSys.BaseClasses.Icons.Sensor;

        parameter ThermoOpt isobaric=true
          "Options for material and thermal subconnectors"
          annotation (Dialog(compact=true));

        // Material
        replaceable Material.Density material if not isobaric
          "Type of sensor"
          annotation (Placement(transformation(extent={{-70,-10},{-50,10}})));

        // Thermal
        replaceable Thermal.Temperature thermal if thermoOpt <> ThermoOpt.ClosedAdiabatic
          "Type of condition"
          annotation (Placement(transformation(extent={{50,-10},{70,10}})));

        FCSys.Connectors.RealOutputBus y "Output bus for measurements"
          annotation (HideResult=not (internalMaterial or internalLin1 or
              internalLin1 or internalLin2 or internalEntropy), Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,-100}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,-98})));

      equation
        // Material
        connect(material.y, y.normal) annotation (Line(
            points={{-60,-10},{-60,-70},{5.55112e-16,-70},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));

        // Thermal
        connect(thermal.y, y.thermal) annotation (Line(
            points={{60,-10},{60,-70},{5.55112e-16,-70},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));

        annotation (defaultComponentName="speciesFaceSensor", Icon(graphics={
                Line( points={{-70,0},{-100,0}},
                      color={0,0,0},
                      smooth=Smooth.None),Line(
                      points={{100,0},{70,0}},
                      color={0,0,0},
                      smooth=Smooth.None)}));
      end PartialSpecies;
    end BaseClasses;
  end FaceDifferential;

  package BaseClasses "Base classes (not for direct use)"
    extends Modelica.Icons.BasesPackage;

    partial model PartialSensorNonideal
      "Partial sensor with dynamics, saturation, and noise"

      extends PartialSensor;

      replaceable PartialSensor sensor constrainedby PartialSensor "Sensor"
        annotation (choicesAllMatching=true, Placement(transformation(extent={{
                -10,44},{10,64}})));

      parameter Real y_min=-y_max
        "<html>Lower limit of measurement (<i>y</i><sub>min</sub>)</html>";
      parameter Real y_max=Modelica.Constants.inf
        "<html>Upper limit of measurement (<i>y</i><sub>max</sub>)</html>";
      parameter Real k_noise=0
        "<html>Amplitude of noise signal (<i>k</i><sub>noise</sub>)</html>";

      replaceable Modelica.Blocks.Continuous.FirstOrder dynamics(each initType=
            Modelica.Blocks.Types.Init.SteadyState, each T=Modelica.Constants.eps)
        constrainedby Modelica.Blocks.Interfaces.SISO "Dynamics" annotation (
          choicesAllMatching=true, Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,16})));
      Modelica.Blocks.Nonlinear.Limiter saturation(uMax=y_max, uMin=y_min)
        "Saturation" annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,-40})));
      FCSys.Blocks.Math.AddSkipInclIncl add "Addition" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,-10})));
      replaceable FCSys.Blocks.Continuous.Sources.RandomNormal noise "Noise"
        annotation (Placement(transformation(extent={{-40,-20},{-20,0}})));

    equation
      connect(dynamics.y, add.u_1[1]) annotation (Line(
          points={{-1.40998e-15,5},{2.16493e-15,-1}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(add.y[1], saturation.u) annotation (Line(
          points={{-1.05471e-15,-19},{-1.05471e-15,-27.5},{2.87043e-15,-27.5},{
              2.87043e-15,-28}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(noise.y, add.u_2) annotation (Line(
          points={{-19,-10},{-14,-10},{-14,-10},{-9,-10}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(dynamics.u, sensor.y) annotation (Line(
          points={{2.87043e-15,28},{2.87043e-15,36},{6.10623e-16,36},{
              6.10623e-16,44}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(y, saturation.y) annotation (Line(
          points={{5.55112e-16,-100},{5.55112e-16,-75},{-1.40998e-15,-75},{-1.40998e-15,
              -51}},
          color={0,0,127},
          smooth=Smooth.None));
      annotation (Icon(coordinateSystem(
            preserveAspectRatio=true,
            extent={{-100,-100},{100,100}},
            grid={2,2}), graphics), Diagram(coordinateSystem(
            preserveAspectRatio=true,
            extent={{-100,-100},{100,100}},
            grid={2,2}), graphics));
    end PartialSensorNonideal;

    partial model PartialSensor "Partial model for a sensor"

      extends FCSys.BaseClasses.Icons.Sensor;
      FCSys.Connectors.RealOutput y "Measurement" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,-100})));
    end PartialSensor;
  end BaseClasses;
  annotation (Documentation(info="<html>
  <p>
<b>Licensed by the Georgia Tech Research Corporation under the Modelica License 2</b><br>
Copyright 2007&ndash;2013, Georgia Tech Research Corporation.
</p>
<p>
<i>This Modelica package is <u>free</u> software and the use is completely at <u>your own risk</u>;
it can be redistributed and/or modified under the terms of the Modelica License 2. For license conditions (including the
disclaimer of warranty) see <a href=\"modelica://FCSys.UsersGuide.ModelicaLicense2\">
FCSys.UsersGuide.ModelicaLicense2</a> or visit <a href=\"http://www.modelica.org/licenses/ModelicaLicense2\">
http://www.modelica.org/licenses/ModelicaLicense2</a>.</i>
</p></html>"));
end Sensors;

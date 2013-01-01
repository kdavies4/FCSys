within FCSys;
package Sensors "Models to measure conditions"

  extends Modelica.Icons.SensorsPackage;

  // TODO:  Recheck this package, fix errors and warnings.

  package Chemical
    "<html>Sensors for the <a href=\"modelica://FCSys.Connectors.BaseClasses.PartialChemical\">Chemical</a> and <a href=\"modelica://FCSys.Connectors.ChemicalBus\">ChemicalBus</a> connectors</html>"
    extends Modelica.Icons.Package;

    package Phases
      "<html>Sensors for the <a href=\"modelica://FCSys.Connectors.ChemicalBus\">ChemicalBus</a> connector, e.g., of a <a href=\"modelica://FCSys.Subregions.Phase\">Phase</a> model (multi-species)</html>"
      extends Modelica.Icons.Package;
      model Phase "Sensor for a phase with all species conditionally included"

        extends FCSys.Sensors.Chemical.Phases.NullPhase;

        // Conditionally include species.
        parameter Boolean inclC=false "Carbon (C)" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        replaceable FCSys.Sensors.Chemical.Species C(final n_vel=n_vel) if
          inclC "Model"
          annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

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
        FCSys.Sensors.Chemical.Species C19HF37O5S(final n_vel=n_vel) if
          inclC19HF37O5S "Model"
          annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

        parameter Boolean 'incle-'=false
          "<html>Electrons (e<sup>-</sup>)</html>" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FCSys.Sensors.Chemical.Species 'e-'(final n_vel=n_vel) if 'incle-'
          "Model"
          annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

        parameter Boolean inclH2=false "<html>Hydrogen (H<sub>2</sub>)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FCSys.Sensors.Chemical.Species H2(final n_vel=n_vel) if inclH2 "Model"
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
        FCSys.Sensors.Chemical.Species H2O(final n_vel=n_vel) if inclH2O
          "Model"
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
        FCSys.Sensors.Chemical.Species 'H+'(final n_vel=n_vel) if 'inclH+'
          "Model"
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
        FCSys.Sensors.Chemical.Species N2(final n_vel=n_vel) if inclN2 "Model"
          annotation (Dialog(
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
        FCSys.Sensors.Chemical.Species O2(final n_vel=n_vel) if inclO2 "Model"
          annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

      equation
        // C
        connect(C.chemical, chemical.C) annotation (Line(
            points={{6.10623e-16,-40},{5.55112e-16,5.55112e-16}},
            color={208,104,0},
            smooth=Smooth.None));
        connect(C.mu, y.C.mu) annotation (Line(
            points={{-4,-50},{-4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(C.phi, y.C.phi) annotation (Line(
            points={{6.10623e-16,-50},{5.55112e-16,-50},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(C.s, y.C.s) annotation (Line(
            points={{4,-50},{4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));

        // C19HF37O5S
        connect(C19HF37O5S.chemical, chemical.C19HF37O5S) annotation (Line(
            points={{6.10623e-16,-40},{5.55112e-16,-40},{5.55112e-16,
                5.55112e-16}},
            color={208,104,0},
            smooth=Smooth.None));

        connect(C19HF37O5S.mu, y.C19HF37O5S.mu) annotation (Line(
            points={{-4,-50},{-4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(C19HF37O5S.phi, y.C19HF37O5S.phi) annotation (Line(
            points={{6.10623e-16,-50},{5.55112e-16,-50},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(C19HF37O5S.s, y.C19HF37O5S.s) annotation (Line(
            points={{4,-50},{4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));

        // e-
        connect('e-'.chemical, chemical.'e-') annotation (Line(
            points={{6.10623e-16,-40},{5.55112e-16,-40},{5.55112e-16,
                5.55112e-16}},
            color={208,104,0},
            smooth=Smooth.None));

        connect('e-'.mu, y.'e-'.mu) annotation (Line(
            points={{-4,-50},{-4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));
        connect('e-'.phi, y.'e-'.phi) annotation (Line(
            points={{6.10623e-16,-50},{5.55112e-16,-50},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));

        connect('e-'.s, y.'e-'.s) annotation (Line(
            points={{4,-50},{4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));

        // H2
        connect(H2.chemical, chemical.H2) annotation (Line(
            points={{6.10623e-16,-40},{5.55112e-16,-40},{5.55112e-16,
                5.55112e-16}},
            color={208,104,0},
            smooth=Smooth.None));

        connect(H2.mu, y.H2.mu) annotation (Line(
            points={{-4,-50},{-4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(H2.phi, y.H2.phi) annotation (Line(
            points={{6.10623e-16,-50},{5.55112e-16,-50},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(H2.s, y.H2.s) annotation (Line(
            points={{4,-50},{4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));

        // H2O
        connect(H2O.chemical, chemical.H2O) annotation (Line(
            points={{6.10623e-16,-40},{5.55112e-16,-40},{5.55112e-16,
                5.55112e-16}},
            color={208,104,0},
            smooth=Smooth.None));

        connect(H2O.mu, y.H2O.mu) annotation (Line(
            points={{-4,-50},{-4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(H2O.phi, y.H2O.phi) annotation (Line(
            points={{6.10623e-16,-50},{5.55112e-16,-50},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(H2O.s, y.H2O.s) annotation (Line(
            points={{4,-50},{4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));

        // H+
        connect('H+'.chemical, chemical.'H+') annotation (Line(
            points={{6.10623e-16,-40},{5.55112e-16,-40},{5.55112e-16,
                5.55112e-16}},
            color={208,104,0},
            smooth=Smooth.None));

        connect('H+'.mu, y.'H+'.mu) annotation (Line(
            points={{-4,-50},{-4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));
        connect('H+'.phi, y.'H+'.phi) annotation (Line(
            points={{6.10623e-16,-50},{5.55112e-16,-50},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));

        connect('H+'.s, y.'H+'.s) annotation (Line(
            points={{4,-50},{4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));

        // N2
        connect(N2.chemical, chemical.N2) annotation (Line(
            points={{6.10623e-16,-40},{5.55112e-16,-40},{5.55112e-16,
                5.55112e-16}},
            color={208,104,0},
            smooth=Smooth.None));

        connect(N2.mu, y.N2.mu) annotation (Line(
            points={{-4,-50},{-4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(N2.phi, y.N2.phi) annotation (Line(
            points={{6.10623e-16,-50},{-4.87687e-22,-50},{-4.87687e-22,-100},{
                5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(N2.s, y.N2.s) annotation (Line(
            points={{4,-50},{4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));

        // O2
        connect(O2.chemical, chemical.O2) annotation (Line(
            points={{6.10623e-16,-40},{5.55112e-16,-40},{5.55112e-16,
                5.55112e-16}},
            color={208,104,0},
            smooth=Smooth.None));

        connect(O2.mu, y.O2.mu) annotation (Line(
            points={{-4,-50},{-4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(O2.phi, y.O2.phi) annotation (Line(
            points={{6.10623e-16,-50},{-4.87687e-22,-50},{-4.87687e-22,-100},{
                5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(O2.s, y.O2.s) annotation (Line(
            points={{4,-50},{4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));

        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="phaseChemicalSensor",
          Diagram(graphics),
          Icon(graphics));
      end Phase;

      model Gas "Sensor for gas"

        extends FCSys.Sensors.Chemical.Phases.NullPhase;

        // Conditionally include species.
        parameter Boolean inclC=false "Carbon (C)" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        replaceable FCSys.Sensors.Chemical.Species C(final n_vel=n_vel) if
          inclC "Model"
          annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

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
        FCSys.Sensors.Chemical.Species C19HF37O5S(final n_vel=n_vel) if
          inclC19HF37O5S "Model"
          annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

        parameter Boolean 'incle-'=false
          "<html>Electrons (e<sup>-</sup>)</html>" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FCSys.Sensors.Chemical.Species 'e-'(final n_vel=n_vel) if 'incle-'
          "Model"
          annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

        parameter Boolean inclH2=false "<html>Hydrogen (H<sub>2</sub>)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FCSys.Sensors.Chemical.Species H2(final n_vel=n_vel) if inclH2 "Model"
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
        FCSys.Sensors.Chemical.Species H2O(final n_vel=n_vel) if inclH2O
          "Model"
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
        FCSys.Sensors.Chemical.Species 'H+'(final n_vel=n_vel) if 'inclH+'
          "Model"
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
        FCSys.Sensors.Chemical.Species N2(final n_vel=n_vel) if inclN2 "Model"
          annotation (Dialog(
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
        FCSys.Sensors.Chemical.Species O2(final n_vel=n_vel) if inclO2 "Model"
          annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

      equation
        // H2
        connect(H2.chemical, chemical.H2) annotation (Line(
            points={{6.10623e-16,-40},{5.55112e-16,-40},{5.55112e-16,
                5.55112e-16}},
            color={208,104,0},
            smooth=Smooth.None));

        connect(H2.mu, y.H2.mu) annotation (Line(
            points={{-4,-50},{-4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(H2.phi, y.H2.phi) annotation (Line(
            points={{6.10623e-16,-50},{5.55112e-16,-50},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(H2.s, y.H2.s) annotation (Line(
            points={{4,-50},{4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));

        // H2O
        connect(H2O.chemical, chemical.H2O) annotation (Line(
            points={{6.10623e-16,-40},{5.55112e-16,-40},{5.55112e-16,
                5.55112e-16}},
            color={208,104,0},
            smooth=Smooth.None));

        connect(H2O.mu, y.H2O.mu) annotation (Line(
            points={{-4,-50},{-4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(H2O.phi, y.H2O.phi) annotation (Line(
            points={{6.10623e-16,-50},{5.55112e-16,-50},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(H2O.s, y.H2O.s) annotation (Line(
            points={{4,-50},{4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));

        // N2
        connect(N2.chemical, chemical.N2) annotation (Line(
            points={{6.10623e-16,-40},{5.55112e-16,-40},{5.55112e-16,
                5.55112e-16}},
            color={208,104,0},
            smooth=Smooth.None));

        connect(N2.mu, y.N2.mu) annotation (Line(
            points={{-4,-50},{-4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(N2.phi, y.N2.phi) annotation (Line(
            points={{6.10623e-16,-50},{-4.87687e-22,-50},{-4.87687e-22,-100},{
                5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(N2.s, y.N2.s) annotation (Line(
            points={{4,-50},{4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));

        // O2
        connect(O2.chemical, chemical.O2) annotation (Line(
            points={{6.10623e-16,-40},{5.55112e-16,-40},{5.55112e-16,
                5.55112e-16}},
            color={208,104,0},
            smooth=Smooth.None));

        connect(O2.mu, y.O2.mu) annotation (Line(
            points={{-4,-50},{-4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(O2.phi, y.O2.phi) annotation (Line(
            points={{6.10623e-16,-50},{-4.87687e-22,-50},{-4.87687e-22,-100},{
                5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(O2.s, y.O2.s) annotation (Line(
            points={{4,-50},{4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));

        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="phaseChemicalSensor");
      end Gas;

      model Graphite "Sensor for graphite"

        extends FCSys.Sensors.Chemical.Phases.NullPhase;

        // Conditionally include species.
        parameter Boolean inclC=false "Carbon (C)" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        replaceable FCSys.Sensors.Chemical.Species C(final n_vel=n_vel) if
          inclC "Model"
          annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

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
        FCSys.Sensors.Chemical.Species C19HF37O5S(final n_vel=n_vel) if
          inclC19HF37O5S "Model"
          annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

        parameter Boolean 'incle-'=false
          "<html>Electrons (e<sup>-</sup>)</html>" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FCSys.Sensors.Chemical.Species 'e-'(final n_vel=n_vel) if 'incle-'
          "Model"
          annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

        parameter Boolean inclH2=false "<html>Hydrogen (H<sub>2</sub>)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FCSys.Sensors.Chemical.Species H2(final n_vel=n_vel) if inclH2 "Model"
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
        FCSys.Sensors.Chemical.Species H2O(final n_vel=n_vel) if inclH2O
          "Model"
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
        FCSys.Sensors.Chemical.Species 'H+'(final n_vel=n_vel) if 'inclH+'
          "Model"
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
        FCSys.Sensors.Chemical.Species N2(final n_vel=n_vel) if inclN2 "Model"
          annotation (Dialog(
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
        FCSys.Sensors.Chemical.Species O2(final n_vel=n_vel) if inclO2 "Model"
          annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

      equation
        // C
        connect(C.chemical, chemical.C) annotation (Line(
            points={{6.10623e-16,-40},{5.55112e-16,5.55112e-16}},
            color={208,104,0},
            smooth=Smooth.None));
        connect(C.mu, y.C.mu) annotation (Line(
            points={{-4,-50},{-4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(C.phi, y.C.phi) annotation (Line(
            points={{6.10623e-16,-50},{5.55112e-16,-50},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(C.s, y.C.s) annotation (Line(
            points={{4,-50},{4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));

        // e-
        connect('e-'.chemical, chemical.'e-') annotation (Line(
            points={{6.10623e-16,-40},{5.55112e-16,-40},{5.55112e-16,
                5.55112e-16}},
            color={208,104,0},
            smooth=Smooth.None));

        connect('e-'.mu, y.'e-'.mu) annotation (Line(
            points={{-4,-50},{-4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));
        connect('e-'.phi, y.'e-'.phi) annotation (Line(
            points={{6.10623e-16,-50},{5.55112e-16,-50},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));

        connect('e-'.s, y.'e-'.s) annotation (Line(
            points={{4,-50},{4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));

        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="phaseChemicalSensor");
      end Graphite;

      model Ionomer "Sensor for ionomer"
        extends FCSys.Sensors.Chemical.Phases.NullPhase;

        // Conditionally include species.
        parameter Boolean inclC=false "Carbon (C)" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        replaceable FCSys.Sensors.Chemical.Species C(final n_vel=n_vel) if
          inclC "Model"
          annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

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
        FCSys.Sensors.Chemical.Species C19HF37O5S(final n_vel=n_vel) if
          inclC19HF37O5S "Model"
          annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

        parameter Boolean 'incle-'=false
          "<html>Electrons (e<sup>-</sup>)</html>" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FCSys.Sensors.Chemical.Species 'e-'(final n_vel=n_vel) if 'incle-'
          "Model"
          annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

        parameter Boolean inclH2=false "<html>Hydrogen (H<sub>2</sub>)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FCSys.Sensors.Chemical.Species H2(final n_vel=n_vel) if inclH2 "Model"
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
        FCSys.Sensors.Chemical.Species H2O(final n_vel=n_vel) if inclH2O
          "Model"
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
        FCSys.Sensors.Chemical.Species 'H+'(final n_vel=n_vel) if 'inclH+'
          "Model"
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
        FCSys.Sensors.Chemical.Species N2(final n_vel=n_vel) if inclN2 "Model"
          annotation (Dialog(
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
        FCSys.Sensors.Chemical.Species O2(final n_vel=n_vel) if inclO2 "Model"
          annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));

      equation
        // C19HF37O5S
        connect(C19HF37O5S.chemical, chemical.C19HF37O5S) annotation (Line(
            points={{6.10623e-16,-40},{5.55112e-16,-40},{5.55112e-16,
                5.55112e-16}},
            color={208,104,0},
            smooth=Smooth.None));

        connect(C19HF37O5S.mu, y.C19HF37O5S.mu) annotation (Line(
            points={{-4,-50},{-4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(C19HF37O5S.phi, y.C19HF37O5S.phi) annotation (Line(
            points={{6.10623e-16,-50},{5.55112e-16,-50},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(C19HF37O5S.s, y.C19HF37O5S.s) annotation (Line(
            points={{4,-50},{4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));

        // H+
        connect('H+'.chemical, chemical.'H+') annotation (Line(
            points={{6.10623e-16,-40},{5.55112e-16,-40},{5.55112e-16,
                5.55112e-16}},
            color={208,104,0},
            smooth=Smooth.None));

        connect('H+'.mu, y.'H+'.mu) annotation (Line(
            points={{-4,-50},{-4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));
        connect('H+'.phi, y.'H+'.phi) annotation (Line(
            points={{6.10623e-16,-50},{5.55112e-16,-50},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));

        connect('H+'.s, y.'H+'.s) annotation (Line(
            points={{4,-50},{4,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));

        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="phaseChemicalSensor",
          Diagram(graphics));
      end Ionomer;

      model NullPhase "Empty sensor for a phase (no species)"
        extends FCSys.BaseClasses.Icons.Sensor;
        parameter Integer n_vel(
          final min=1,
          final max=3) = 1
          "<html>Number of components of velocity (<i>n</i><sub>vel</sub>)</html>"
          annotation (HideResult=true);

        FCSys.Connectors.ChemicalBus chemical
          "Multi-species connector for material" annotation (Placement(
              transformation(extent={{-10,-10},{10,10}}), iconTransformation(
                extent={{-10,-10},{10,10}})));

        Connectors.RealOutputBus y "Bus of measurements" annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,-100}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,-100})));
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="phaseChemicalSensor",
          Diagram(graphics),
          Icon(graphics));
      end NullPhase;
    end Phases;

    model Species
      "<html>Sensors for the <a href=\"modelica://FCSys.Connectors.BaseClasses.PartialChemical\">Chemical</a> connector, e.g., of a <a href=\"modelica://FCSys.Subregions.Species\">Species</a> model (single-species)</html>"

      extends FCSys.BaseClasses.Icons.Sensor;

      parameter Integer n_vel(
        final min=1,
        final max=3) = 1
        "<html>Number of components of velocity (<i>n</i><sub>vel</sub>)</html>"
        annotation (HideResult=true);

      Connectors.ChemicalInput chemical(final n_vel=n_vel)
        "Single-species connector for material" annotation (Placement(
            transformation(extent={{-10,-10},{10,10}}), iconTransformation(
              extent={{-10,-10},{10,10}})));

      Connectors.RealOutput mu(final unit="m.l2/(N.T2)")
        "Internal signal for electrochemical potential" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-30,-100}),iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-40,-100})));
      Connectors.RealOutput mphi[n_vel](each final unit="l.m/(N.T)")
        "Internal signal for specific mass times velocity" annotation (
          Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,-100}),iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,-100})));
      Connectors.RealOutput sT(final unit="l2.m/(N.T)")
        "Internal signal for specific entropy times temperature" annotation (
          Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={30,-100}), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={40,-100})));

    equation
      // Measurements
      mu = chemical.mu;
      mphi = chemical.mphi;
      sT = chemical.sT;

      // Conditions
      0 = chemical.Ndot "No current";
      zeros(n_vel) = chemical.mphi "Zero outflow specific mass times velocity";
      0 = chemical.sT "Zero outflow specific entropy times temperature";
      annotation (
        defaultComponentName="speciesChemicalSensor",
        Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                {100,100}}),graphics),
        Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
                100,100}}), graphics={Line(points={{-40,-100},{-40,-58}}, color
              ={0,0,127}),Line(points={{40,-100},{40,-58}}, color={0,0,127}),
              Text(
                  extent={{-100,-20},{100,-50}},
                  lineColor={127,127,127},
                  textString="mu mphi sT ")}));
    end Species;
    annotation (Documentation(info="<html><p>Since the connectors in
<a href=\"modelica://FCSys\">FCSys</a> are hierarchical
(see the <a href=\"modelica://FCSys.Connectors\">Connectors</a> package),
the models for the sensors must be as well.  A
<a href=\"modelica://FCSys.Connectors.BaseClasses.PartialChemical\">Chemical</a>
connector
is used in <a href=\"modelica://FCSys.Subregions.Species\">Species</a> models,
and there is a corresponding <a href=\"modelica://FCSys.Sensors.Chemical.Species.Species\">Species
sensor</a> model in this package. The
<a href=\"modelica://FCSys.Connectors.ChemicalBus\">ChemicalBus</a>
connector is used in <a href=\"modelica://FCSys.Subregions.Phases\">Phase</a> models,
and there are corresponding <a href=\"modelica://FCSys.Sensors.Chemical.Phases\">Phase
sensor</a> models.
</p></html>"));
  end Chemical;

  package InertAmagat
    "<html>Sensors for the <a href=\"modelica://FCSys.Connectors.InertAmagat\">InertAmagat</a> connector</html>"
    extends Modelica.Icons.Package;

    model Phase
      "<html>Sensor for the <a href=\"modelica://FCSys.Connectors.InertAmagat\">InertAmagat</a> connector, e.g., of a <a href=\"modelica://FCSys.Subregions.Phases\">Phase</a> model</html>"
      extends FCSys.BaseClasses.Icons.Sensor;

      parameter Integer n_vel(
        final min=0,
        final max=3) = 1
        "<html>Number of components of velocity (<i>n</i><sub>vel</sub>)</html>";

      Pressure pressure(final n_vel=n_vel) "Sensor for pressure"
        annotation (Placement(transformation(extent={{-50,-60},{-30,-40}})));

      Velocity velocity(final n_vel=n_vel) if n_vel > 1 "Sensor for velocity"
        annotation (Placement(transformation(extent={{-10,-60},{10,-40}})));

      Temperature temperature(final n_vel=n_vel) "Sensor for temperature"
        annotation (Placement(transformation(extent={{30,-60},{50,-40}})));
      Connectors.InertAmagat inert(final n_vel=n_vel)
        "Single-species connector for linear momentum and entropy, with additivity of volume"
        annotation (Placement(transformation(extent={{-10,-10},{10,10}}),
            iconTransformation(extent={{-10,-10},{10,10}})));
      Connectors.RealOutputBus y "Input bus for external signal sources"
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
      annotation (
        defaultComponentName="phaseInertSensor",
        Diagram(graphics),
        Icon(graphics));
    end Phase;

    model Pressure "Prescribed pressure"
      extends BaseClasses.PartialSensor(redeclare Connectors.RealOutput y(
            final unit="m/(l.T2)"));
    equation
      y = inert.p "Measurement";
      annotation (
        defaultComponentPrefixes="replaceable",
        defaultComponentName="volSensor",
        Diagram(graphics),
        Icon(graphics={Text(
                  extent={{-100,-20},{100,-50}},
                  lineColor={127,127,127},
                  textString="P")}));
    end Pressure;

    model Velocity "Measured velocity"
      extends FCSys.BaseClasses.Icons.Sensor;

      parameter Integer n_vel(
        final min=0,
        final max=3) = 1
        "<html>Number of components of velocity (<i>n</i><sub>vel</sub>)</html>"
        annotation (HideResult=true);

      Connectors.InertDalton inert(final n_vel=n_vel)
        "Connector for linear momentum and entropy, with additivity of pressure"
        annotation (HideResult=true, Placement(transformation(extent={{-10,-10},
                {10,10}}), iconTransformation(extent={{-10,-10},{10,10}})));
      FCSys.Connectors.RealOutput y[n_vel](final unit="l/T") "Measurement"
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,-100})));

    equation
      // Measurement
      y = inert.phi;

      // Conditions
      0 = inert.V "No (additional) volume";
      zeros(n_vel) = inert.mPhidot "No force";
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
      extends BaseClasses.PartialSensor(redeclare Connectors.RealOutput y(
            final unit="l2.m/(N.T2)", displayUnit="K"));
    equation
      y = inert.T "Measurement";
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

        parameter Integer n_vel(
          final min=1,
          final max=3) = 1
          "<html>Number of components of velocity (<i>n</i><sub>vel</sub>)</html>"
          annotation (HideResult=true);

        FCSys.Connectors.InertAmagat inert(final n_vel=n_vel)
          "Connector for linear momentum and entropy, with additivity of volume"
          annotation (HideResult=true, Placement(transformation(extent={{-10,-10},
                  {10,10}}), iconTransformation(extent={{-10,-10},{10,10}})));
        Connectors.RealOutput y "Measurement" annotation (Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,-100}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,-100})));
      equation
        // Conditions
        0 = inert.V "No (additional) volume";
        zeros(n_vel) = inert.mPhidot "No force";
        0 = inert.Qdot "Adiabatic";
        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="phaseInertSensor",
          Diagram(graphics),
          Icon(graphics));
      end PartialSensor;
    end BaseClasses;
  end InertAmagat;

  package InertDalton
    "<html>Sensors for the <a href=\"modelica://FCSys.Connectors.InertDalton\">InertDalton</a> connector</html>"
    extends Modelica.Icons.Package;

    model Species
      "<html>Sensor for the <a href=\"modelica://FCSys.Connectors.InertDalton\">InertDalton</a> connector, e.g., of a <a href=\"modelica://FCSys.Subregions.Species\">Species</a> model</html>"

      extends FCSys.BaseClasses.Icons.Sensor;

      parameter Integer n_vel(
        final min=0,
        final max=3) = 1
        "<html>Number of components of velocity (<i>n</i><sub>vel</sub>)</html>";

      Volume volume(final n_vel=n_vel) "Sensor for volume"
        annotation (Placement(transformation(extent={{-50,-60},{-30,-40}})));

      Velocity velocity(final n_vel=n_vel) if n_vel > 1 "Sensor for velocity"
        annotation (Placement(transformation(extent={{-10,-60},{10,-40}})));

      Temperature temperature(final n_vel=n_vel) "Sensor for temperature"
        annotation (Placement(transformation(extent={{30,-60},{50,-40}})));

      Connectors.InertDalton inert(final n_vel=n_vel)
        "Single-species connector for linear momentum and entropy, with additivity of pressure"
        annotation (Placement(transformation(extent={{-10,-10},{10,10}}),
            iconTransformation(extent={{-10,-10},{10,10}})));
      Connectors.RealOutputBus y "Input bus for external signal sources"
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
      annotation (
        defaultComponentName="speciesInertSensor",
        Diagram(graphics),
        Icon(graphics));
    end Species;

    model Volume "Measured volume"
      extends BaseClasses.PartialSensor(redeclare Connectors.RealOutput y(
            final unit="l3"));

    equation
      y = inert.V;

      annotation (Icon(graphics={Text(
                  extent={{-100,-20},{100,-50}},
                  lineColor={127,127,127},
                  textString="V")}));
    end Volume;

    model Velocity "Measured velocity"
      extends FCSys.BaseClasses.Icons.Sensor;

      parameter Integer n_vel(
        final min=0,
        final max=3) = 1
        "<html>Number of components of velocity (<i>n</i><sub>vel</sub>)</html>"
        annotation (HideResult=true);

      Connectors.InertDalton inert(final n_vel=n_vel)
        "Connector for linear momentum and entropy, with additivity of pressure"
        annotation (HideResult=true, Placement(transformation(extent={{-10,-10},
                {10,10}}), iconTransformation(extent={{-10,-10},{10,10}})));
      FCSys.Connectors.RealOutput y[n_vel](final unit="l/T") "Measurement"
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,-100})));

    equation
      // Measurement
      y = inert.phi;

      // Conditions
      inert.p = 0 "No (additional) pressure";
      inert.mPhidot = zeros(n_vel) "No force";
      0 = inert.Qdot "Adiabatic";

      annotation (Icon(graphics={Text(
                  extent={{-100,-20},{100,-50}},
                  lineColor={127,127,127},
                  textString="phi")}));
    end Velocity;

    model Temperature "Measured temperature"

      extends FCSys.Sensors.InertDalton.BaseClasses.PartialSensor(redeclare
          Connectors.RealOutput y(final unit="l2.m/(N.T2)", displayUnit="K"));

    equation
      y = inert.T "Measurement";

      annotation (Icon(graphics={Text(
                  extent={{-100,-20},{100,-50}},
                  lineColor={127,127,127},
                  textString="T")}));
    end Temperature;

    package BaseClasses "Base classes (not for direct use)"
      extends Modelica.Icons.BasesPackage;
      partial model PartialSensor "Partial model for a sensor"
        extends FCSys.Sensors.BaseClasses.PartialSensor;

        parameter Integer n_vel(
          final min=0,
          final max=3) = 1
          "<html>Number of components of velocity (<i>n</i><sub>vel</sub>)</html>"
          annotation (HideResult=true);

        Connectors.InertDalton inert(final n_vel=n_vel)
          "Connector for linear momentum and entropy, with additivity of pressure"
          annotation (HideResult=true, Placement(transformation(extent={{-10,-10},
                  {10,10}}), iconTransformation(extent={{-10,-10},{10,10}})));

      equation
        // Conditions
        0 = inert.p "No (additional) pressure";
        zeros(n_vel) = inert.mPhidot "No force";
        0 = inert.Qdot "Adiabatic";

        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="speciesInertSensor",
          Diagram(graphics),
          Icon(graphics));
      end PartialSensor;
    end BaseClasses;
  end InertDalton;

  package Face
    "<html>Sensors for <a href=\"modelica://FCSys.Connectors.BaseClasses.PartialFace\">Face</a> and <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connectors</html>"

    extends Modelica.Icons.Package;

    model Subregion
      "<html>Sensor for a face of a <a href=\"modelica://FCSys.Subregions.Subregion\">Region</a> or <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a> model</html>"
      import FCSys.BaseClasses.Axis;
      extends FCSys.BaseClasses.Icons.Sensor;

      parameter FCSys.BaseClasses.Axis axis=Axis.x "Axis normal to the face";

      replaceable FCSys.Sensors.Face.Phases.Gas gas(final axis=axis) "Gas"
        annotation (Dialog(group="Phases", __Dymola_descriptionLabel=true),
          Placement(transformation(extent={{-10,-50},{10,-30}})));

      replaceable FCSys.Sensors.Face.Phases.Graphite graphite(final axis=axis)
        "Graphite" annotation (Dialog(group="Phases", __Dymola_descriptionLabel
            =true), Placement(transformation(extent={{-10,-50},{10,-30}})));

      replaceable FCSys.Sensors.Face.Phases.Ionomer ionomer(final axis=axis)
        "Ionomer" annotation (Dialog(group="Phases", __Dymola_descriptionLabel=
              true), Placement(transformation(extent={{-10,-50},{10,-30}})));

      FCSys.Connectors.FaceBus face
        "Connector for material, linear momentum, and entropy of multiple species"
        annotation (Placement(transformation(extent={{-10,-10},{10,10}}),
            iconTransformation(extent={{-10,-10},{10,10}})));
      Connectors.RealOutputBus y "Bus of measurements" annotation (Placement(
            transformation(
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
        Icon(graphics={Text(
                  extent={{-100,-20},{100,-60}},
                  lineColor={127,127,127},
                  fillColor={0,0,127},
                  fillPattern=FillPattern.Solid,
                  textString="mu phi T")}),
        Diagram(graphics));
    end Subregion;

    package Phases
      "<html>Sensors for the <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connector, e.g., of a <a href=\"modelica://FCSys.Subregions.Phase\">Phase</a> model (multi-species)</html>"
      extends Modelica.Icons.Package;
      model Phase "Sensor for a phase with all species conditionally included"

        extends FCSys.BaseClasses.Icons.Sensor;

        extends FCSys.Sensors.Face.Phases.BaseClasses.NullPhase;

        // Conditionally include species.
        parameter Boolean inclC=false "Carbon (C)" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        replaceable Species.Species C if inclC "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclC), Placement(transformation(extent={{-10,-10},{10,10}})));
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
        replaceable Species.Species C19HF37O5S(final axis=axis) if
          inclC19HF37O5S "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclC19HF37O5S), Placement(transformation(extent={{-10,-10},
                  {10,10}})));
        parameter Boolean 'incle-'=false
          "<html>Electrons (e<sup>-</sup>)</html>" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        replaceable Species.Species 'e-'(final axis=axis) if 'incle-' "Model"
          annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable='incle-'), Placement(transformation(extent={{-10,-10},{10,10}})));
        parameter Boolean inclH2=false "<html>Hydrogen (H<sub>2</sub>)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        replaceable Species.Species H2(final axis=axis) if inclH2 "Model"
          annotation (Dialog(
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
        replaceable Species.Species H2O(final axis=axis) if inclH2O "Model"
          annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclH2O), Placement(transformation(extent={{-10,-10},{10,10}})));
        parameter Boolean 'inclH+'=false "<html>Protons (H<sup>+</sup>)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        replaceable Species.Species 'H+'(final axis=axis) if 'inclH+' "Model"
          annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable='inclH+'), Placement(transformation(extent={{-10,-10},{10,10}})));
        parameter Boolean inclN2=false "<html>Nitrogen (N<sub>2</sub>)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        replaceable Species.Species N2(final axis=axis) if inclN2 "Model"
          annotation (Dialog(
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
        replaceable Species.Species O2(final axis=axis) if inclO2 "Model"
          annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclO2), Placement(transformation(extent={{-10,-10},{10,10}})));

      equation
        // Note:  It would be helpful if Modelica allowed elements of expandable
        // connectors to be named by the contents of a string variable and the
        // name of an instance of a model was accessible through a string (like
        // %name is expanded to be the name of the instance of the model).  Then,
        // the connection equations that follow could be generic.

        // C
        connect(C.face.material, face.C.material) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.face.entropy, face.C.entropy) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.face.momentum1, xInt.C.momentumY) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{-20,6.10623e-16},{-20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.face.momentum1, yInt.C.momentumZ) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,6.10623e-16},{
                5.55112e-16,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(C.face.momentum1, zInt.C.momentumX) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{20,6.10623e-16},{20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.face.momentum2, xInt.C.momentumZ) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{-20,6.10623e-16},{-20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.face.momentum2, yInt.C.momentumX) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,6.10623e-16},{
                5.55112e-16,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(C.face.momentum2, zInt.C.momentumY) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{20,6.10623e-16},{20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(y.C, C.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-9.8}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // C19HF37O5S
        connect(C19HF37O5S.face.material, face.C19HF37O5S.material) annotation
          (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.face.entropy, face.C19HF37O5S.entropy) annotation (
            Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.face.momentum1, xInt.C19HF37O5S.momentumY)
          annotation (Line(
            points={{6.10623e-16,6.10623e-16},{-20,6.10623e-16},{-20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.face.momentum1, yInt.C19HF37O5S.momentumZ)
          annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,6.10623e-16},{
                5.55112e-16,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(C19HF37O5S.face.momentum1, zInt.C19HF37O5S.momentumX)
          annotation (Line(
            points={{6.10623e-16,6.10623e-16},{20,6.10623e-16},{20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.face.momentum2, xInt.C19HF37O5S.momentumZ)
          annotation (Line(
            points={{6.10623e-16,6.10623e-16},{-20,6.10623e-16},{-20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.face.momentum2, yInt.C19HF37O5S.momentumX)
          annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,6.10623e-16},{
                5.55112e-16,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(C19HF37O5S.face.momentum2, zInt.C19HF37O5S.momentumY)
          annotation (Line(
            points={{6.10623e-16,6.10623e-16},{20,6.10623e-16},{20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(y.C19HF37O5S, C19HF37O5S.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-9.8}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // e-
        connect('e-'.face.material, face.'e-'.material) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.face.entropy, face.'e-'.entropy) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.face.momentum1, xInt.'e-'.momentumY) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{-20,6.10623e-16},{-20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.face.momentum1, yInt.'e-'.momentumZ) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,6.10623e-16},{
                5.55112e-16,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect('e-'.face.momentum1, zInt.'e-'.momentumX) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{20,6.10623e-16},{20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.face.momentum2, xInt.'e-'.momentumZ) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{-20,6.10623e-16},{-20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.face.momentum2, yInt.'e-'.momentumX) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,6.10623e-16},{
                5.55112e-16,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect('e-'.face.momentum2, zInt.'e-'.momentumY) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{20,6.10623e-16},{20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(y.'e-', 'e-'.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-9.8}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // H2
        connect(H2.face.material, face.H2.material) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.face.entropy, face.H2.entropy) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.face.momentum1, xInt.H2.momentumY) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{-20,6.10623e-16},{-20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.face.momentum1, yInt.H2.momentumZ) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,6.10623e-16},{
                5.55112e-16,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(H2.face.momentum1, zInt.H2.momentumX) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{20,6.10623e-16},{20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.face.momentum2, xInt.H2.momentumZ) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{-20,6.10623e-16},{-20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.face.momentum2, yInt.H2.momentumX) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,6.10623e-16},{
                5.55112e-16,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(H2.face.momentum2, zInt.H2.momentumY) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{20,6.10623e-16},{20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(y.H2, H2.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-9.8}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // H2O
        connect(H2O.face.material, face.H2O.material) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.face.entropy, face.H2O.entropy) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.face.momentum1, xInt.H2O.momentumY) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{-20,6.10623e-16},{-20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.face.momentum1, yInt.H2O.momentumZ) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,6.10623e-16},{
                5.55112e-16,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(H2O.face.momentum1, zInt.H2O.momentumX) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{20,6.10623e-16},{20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.face.momentum2, xInt.H2O.momentumZ) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{-20,6.10623e-16},{-20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.face.momentum2, yInt.H2O.momentumX) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,6.10623e-16},{
                5.55112e-16,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(H2O.face.momentum2, zInt.H2O.momentumY) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{20,6.10623e-16},{20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(y.H2O, H2O.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-9.8}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // H+
        connect('H+'.face.material, face.'H+'.material) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.face.entropy, face.'H+'.entropy) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.face.momentum1, xInt.'H+'.momentumY) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{-20,6.10623e-16},{-20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.face.momentum1, yInt.'H+'.momentumZ) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,6.10623e-16},{
                5.55112e-16,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect('H+'.face.momentum1, zInt.'H+'.momentumX) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{20,6.10623e-16},{20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.face.momentum2, xInt.'H+'.momentumZ) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{-20,6.10623e-16},{-20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.face.momentum2, yInt.'H+'.momentumX) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,6.10623e-16},{
                5.55112e-16,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect('H+'.face.momentum2, zInt.'H+'.momentumY) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{20,6.10623e-16},{20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(y.'H+', 'H+'.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-9.8}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // N2
        connect(N2.face.material, face.N2.material) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.face.entropy, face.N2.entropy) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.face.momentum1, xInt.N2.momentumY) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{-20,6.10623e-16},{-20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.face.momentum1, yInt.N2.momentumZ) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,6.10623e-16},{
                5.55112e-16,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(N2.face.momentum1, zInt.N2.momentumX) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{20,6.10623e-16},{20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.face.momentum2, xInt.N2.momentumZ) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{-20,6.10623e-16},{-20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.face.momentum2, yInt.N2.momentumX) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,6.10623e-16},{
                5.55112e-16,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(N2.face.momentum2, zInt.N2.momentumY) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{20,6.10623e-16},{20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(y.N2, N2.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-9.8}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // O2
        connect(O2.face.material, face.O2.material) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.face.entropy, face.O2.entropy) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.face.momentum1, xInt.O2.momentumY) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{-20,6.10623e-16},{-20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.face.momentum1, yInt.O2.momentumZ) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,6.10623e-16},{
                5.55112e-16,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(O2.face.momentum1, zInt.O2.momentumX) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{20,6.10623e-16},{20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.face.momentum2, xInt.O2.momentumZ) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{-20,6.10623e-16},{-20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.face.momentum2, yInt.O2.momentumX) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,6.10623e-16},{
                5.55112e-16,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(O2.face.momentum2, zInt.O2.momentumY) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{20,6.10623e-16},{20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(y.O2, O2.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-9.8}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="phaseFaceSensor",
          Diagram(graphics));
      end Phase;

      model Gas "Sensor for gas"

        extends FCSys.Sensors.Face.Phases.BaseClasses.NullPhase;

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
        Species.Species H2(final axis=axis) if inclH2 "Model" annotation (
            Dialog(
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
        Species.Species H2O(final axis=axis) if inclH2O "Model" annotation (
            Dialog(
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
        Species.Species N2(final axis=axis) if inclN2 "Model" annotation (
            Dialog(
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
        Species.Species O2(final axis=axis) if inclO2 "Model" annotation (
            Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclO2), Placement(transformation(extent={{-10,-10},{10,10}})));
      equation
        // H2
        connect(H2.face.material, face.H2.material) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.face.entropy, face.H2.entropy) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.face.momentum1, xInt.H2.momentumY) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{-20,6.10623e-16},{-20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.face.momentum1, yInt.H2.momentumZ) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,6.10623e-16},{
                5.55112e-16,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(H2.face.momentum1, zInt.H2.momentumX) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{20,6.10623e-16},{20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.face.momentum2, xInt.H2.momentumZ) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{-20,6.10623e-16},{-20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.face.momentum2, yInt.H2.momentumX) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,6.10623e-16},{
                5.55112e-16,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(H2.face.momentum2, zInt.H2.momentumY) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{20,6.10623e-16},{20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(y.H2, H2.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-9.8}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // H2O
        connect(H2O.face.material, face.H2O.material) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.face.entropy, face.H2O.entropy) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.face.momentum1, xInt.H2O.momentumY) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{-20,6.10623e-16},{-20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.face.momentum1, yInt.H2O.momentumZ) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,6.10623e-16},{
                5.55112e-16,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(H2O.face.momentum1, zInt.H2O.momentumX) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{20,6.10623e-16},{20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.face.momentum2, xInt.H2O.momentumZ) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{-20,6.10623e-16},{-20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.face.momentum2, yInt.H2O.momentumX) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,6.10623e-16},{
                5.55112e-16,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(H2O.face.momentum2, zInt.H2O.momentumY) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{20,6.10623e-16},{20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(y.H2O, H2O.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-9.8}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // N2
        connect(N2.face.material, face.N2.material) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.face.entropy, face.N2.entropy) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.face.momentum1, xInt.N2.momentumY) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{-20,6.10623e-16},{-20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.face.momentum1, yInt.N2.momentumZ) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,6.10623e-16},{
                5.55112e-16,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(N2.face.momentum1, zInt.N2.momentumX) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{20,6.10623e-16},{20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.face.momentum2, xInt.N2.momentumZ) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{-20,6.10623e-16},{-20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.face.momentum2, yInt.N2.momentumX) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,6.10623e-16},{
                5.55112e-16,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(N2.face.momentum2, zInt.N2.momentumY) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{20,6.10623e-16},{20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(y.N2, N2.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-9.8}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // O2
        connect(O2.face.material, face.O2.material) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.face.entropy, face.O2.entropy) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.face.momentum1, xInt.O2.momentumY) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{-20,6.10623e-16},{-20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.face.momentum1, yInt.O2.momentumZ) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,6.10623e-16},{
                5.55112e-16,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(O2.face.momentum1, zInt.O2.momentumX) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{20,6.10623e-16},{20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.face.momentum2, xInt.O2.momentumZ) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{-20,6.10623e-16},{-20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.face.momentum2, yInt.O2.momentumX) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,6.10623e-16},{
                5.55112e-16,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(O2.face.momentum2, zInt.O2.momentumY) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{20,6.10623e-16},{20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(y.O2, O2.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-9.8}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="phaseFaceSensor");
      end Gas;

      model Graphite "Sensor for graphite"

        extends FCSys.Sensors.Face.Phases.BaseClasses.NullPhase;

        // Conditionally include species.
        parameter Boolean inclC=false "Carbon (C)" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        Species.Species C(final axis=axis) if inclC "Model" annotation (Dialog(
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
        Species.Species 'e-'(final axis=axis) if 'incle-' "Model" annotation (
            Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable='incle-'), Placement(transformation(extent={{-10,-10},{10,10}})));

      equation
        // C
        connect(C.face.material, face.C.material) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.face.entropy, face.C.entropy) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.face.momentum1, xInt.C.momentumY) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{-20,6.10623e-16},{-20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.face.momentum1, yInt.C.momentumZ) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,6.10623e-16},{
                5.55112e-16,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(C.face.momentum1, zInt.C.momentumX) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{20,6.10623e-16},{20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.face.momentum2, xInt.C.momentumZ) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{-20,6.10623e-16},{-20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.face.momentum2, yInt.C.momentumX) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,6.10623e-16},{
                5.55112e-16,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(C.face.momentum2, zInt.C.momentumY) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{20,6.10623e-16},{20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(y.C, C.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-9.8}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // e-
        connect('e-'.face.material, face.'e-'.material) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.face.entropy, face.'e-'.entropy) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.face.momentum1, xInt.'e-'.momentumY) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{-20,6.10623e-16},{-20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.face.momentum1, yInt.'e-'.momentumZ) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,6.10623e-16},{
                5.55112e-16,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect('e-'.face.momentum1, zInt.'e-'.momentumX) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{20,6.10623e-16},{20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.face.momentum2, xInt.'e-'.momentumZ) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{-20,6.10623e-16},{-20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.face.momentum2, yInt.'e-'.momentumX) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,6.10623e-16},{
                5.55112e-16,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect('e-'.face.momentum2, zInt.'e-'.momentumY) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{20,6.10623e-16},{20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(y.'e-', 'e-'.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-9.8}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="phaseFaceSensor");
      end Graphite;

      model Ionomer "Sensor for ionomer"

        extends FCSys.Sensors.Face.Phases.BaseClasses.NullPhase;

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
        Species.Species C19HF37O5S(final axis=axis) if inclC19HF37O5S "Model"
          annotation (Dialog(
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
        Species.Species 'H+'(final axis=axis) if 'inclH+' "Model" annotation (
            Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable='inclH+'), Placement(transformation(extent={{-10,-10},{10,10}})));

      equation
        // C19HF37O5S
        connect(C19HF37O5S.face.material, face.C19HF37O5S.material) annotation
          (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.face.entropy, face.C19HF37O5S.entropy) annotation (
            Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.face.momentum1, xInt.C19HF37O5S.momentumY)
          annotation (Line(
            points={{6.10623e-16,6.10623e-16},{-20,6.10623e-16},{-20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.face.momentum1, yInt.C19HF37O5S.momentumZ)
          annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,6.10623e-16},{
                5.55112e-16,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(C19HF37O5S.face.momentum1, zInt.C19HF37O5S.momentumX)
          annotation (Line(
            points={{6.10623e-16,6.10623e-16},{20,6.10623e-16},{20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.face.momentum2, xInt.C19HF37O5S.momentumZ)
          annotation (Line(
            points={{6.10623e-16,6.10623e-16},{-20,6.10623e-16},{-20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.face.momentum2, yInt.C19HF37O5S.momentumX)
          annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,6.10623e-16},{
                5.55112e-16,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(C19HF37O5S.face.momentum2, zInt.C19HF37O5S.momentumY)
          annotation (Line(
            points={{6.10623e-16,6.10623e-16},{20,6.10623e-16},{20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(y.C19HF37O5S, C19HF37O5S.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-9.8}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // H+
        connect('H+'.face.material, face.'H+'.material) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.face.entropy, face.'H+'.entropy) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.face.momentum1, xInt.'H+'.momentumY) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{-20,6.10623e-16},{-20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.face.momentum1, yInt.'H+'.momentumZ) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,6.10623e-16},{
                5.55112e-16,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect('H+'.face.momentum1, zInt.'H+'.momentumX) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{20,6.10623e-16},{20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.face.momentum2, xInt.'H+'.momentumZ) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{-20,6.10623e-16},{-20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.face.momentum2, yInt.'H+'.momentumX) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{5.55112e-16,6.10623e-16},{
                5.55112e-16,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect('H+'.face.momentum2, zInt.'H+'.momentumY) annotation (Line(
            points={{6.10623e-16,6.10623e-16},{20,6.10623e-16},{20,-30}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(y.'H+', 'H+'.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-9.8}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="phaseFaceSensor",
          Diagram(graphics));
      end Ionomer;

      package BaseClasses "Base classes (not for direct use)"
        extends Modelica.Icons.BasesPackage;
        model NullPhase "Empty sensor for a phase (no species)"
          import FCSys.BaseClasses.Axis;
          extends FCSys.BaseClasses.Icons.Sensor;

          parameter FCSys.BaseClasses.Axis axis=Axis.x
            "Axis normal to the face";

          FCSys.Connectors.FaceBus face
            "Multi-species connector for material, linear momentum, and entropy"
            annotation (Placement(transformation(extent={{-10,-10},{10,10}}),
                iconTransformation(extent={{-10,-10},{10,10}})));
          Connectors.RealOutputBus y "Bus of measurements" annotation (
              Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={0,-100}), iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={0,-100})));

        protected
          Connectors.FaceBusInternal xInt if axis == Axis.x
            "Internal connector enabled if x axis" annotation (Placement(
                transformation(extent={{-30,-40},{-10,-20}})));
          Connectors.FaceBusInternal yInt if axis == Axis.y
            "Internal connector enabled if y axis"
            annotation (Placement(transformation(extent={{-10,-40},{10,-20}})));
          Connectors.FaceBusInternal zInt if axis == Axis.z
            "Internal connector enabled if z axis"
            annotation (Placement(transformation(extent={{10,-40},{30,-20}})));
        equation

          connect(xInt, face) annotation (Line(
              points={{-20,-30},{-20,0},{0,0},{0,5.55112e-16},{5.55112e-16,
                  5.55112e-16}},
              color={127,127,127},
              thickness=0.5,
              smooth=Smooth.None));

          connect(yInt, face) annotation (Line(
              points={{5.55112e-16,-30},{5.55112e-16,-22.5},{5.55112e-16,-22.5},
                  {5.55112e-16,-15},{5.55112e-16,5.55112e-16},{5.55112e-16,
                  5.55112e-16}},
              color={127,127,127},
              thickness=0.5,
              smooth=Smooth.None));

          connect(zInt, face) annotation (Line(
              points={{20,-30},{20,0},{0,0},{0,5.55112e-16},{5.55112e-16,
                  5.55112e-16}},
              color={127,127,127},
              thickness=0.5,
              smooth=Smooth.None));

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="phaseFaceSensor",
            Diagram(graphics),
            Icon(graphics));
        end NullPhase;
      end BaseClasses;
    end Phases;

    package Species
      "<html>Sensors for a single <a href=\"modelica://FCSys.Connectors.BaseClasses.BaseClasses.PartialFace\">Face</a> connector, e.g., of a <a href=\"modelica://FCSys.Subregions.Species\">Species</a> model (single-species)</html>"
      extends Modelica.Icons.Package;

      model Species
        "<html>Sensor for a face of a <a href=\"modelica://FCSys.Subregions.Species\">Species</a> model (single-species)</html>"
        import FCSys.BaseClasses.Axis;

        extends FCSys.BaseClasses.Icons.Sensor;

        parameter FCSys.BaseClasses.Axis axis=Axis.x "Axis normal to the face";

        parameter ThermoOpt thermoOpt=ThermoOpt.ClosedAdiabatic
          "Options for material and entropy subconnectors"
          annotation (Dialog(compact=true));

        // Material
        final parameter Boolean open=thermoOpt == ThermoOpt.OpenDiabatic "Open"
          annotation (choices(__Dymola_checkBox=true), Dialog(
            group="Material",
            compact=true,
            __Dymola_label="Included",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FCSys.Sensors.Face.Species.PotentialElectrochemical echemPot if open
          "Type of condition"
          annotation (Placement(transformation(extent={{-70,-6},{-50,14}})));

        // 1st transverse linear momentum
        parameter Boolean slip1=false
          "<html>Viscous (1<sup>st</sup> transverse momentum included)</html>"
          annotation (choices(__Dymola_checkBox=true), Dialog(compact=true));
        FCSys.Sensors.Face.Species.Velocity velocity1 if slip1
          "Type of condition"
          annotation (Placement(transformation(extent={{-30,-6},{-10,14}})));

        // 2nd transverse linear momentum
        parameter Boolean slip2=false
          "<html>Viscous (2<sup>nd</sup> transverse momentum included)</html>"
          annotation (choices(__Dymola_checkBox=true), Dialog(compact=true));
        FCSys.Sensors.Face.Species.Velocity velocity2 if slip2
          "Type of condition"
          annotation (Placement(transformation(extent={{10,-6},{30,14}})));

        // Entropy
        final parameter Boolean diabatic=thermoOpt == ThermoOpt.ClosedDiabatic
             or thermoOpt == ThermoOpt.OpenDiabatic
          "Diabatic (entropy included)"
          annotation (choices(__Dymola_checkBox=true), Dialog(compact=true));
        FCSys.Sensors.Face.Species.Temperature temperature if diabatic
          "Type of condition"
          annotation (Placement(transformation(extent={{50,-6},{70,14}})));

        FCSys.Connectors.FaceGeneric face(
          final thermoOpt=thermoOpt,
          final slip1=slip1,
          final slip2=slip2)
          "Single-species connector for material, linear momentum, and entropy"
          annotation (Placement(transformation(extent={{-10,10},{10,30}}),
              iconTransformation(extent={{-10,-10},{10,10}})));
        Connectors.RealOutputBus y "Output bus for measurements" annotation (
            HideResult=not (internalMaterial or internalLin1 or internalLin1
               or internalLin2 or internalEntropy), Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,-100}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,-98})));

      protected
        Connectors.RealOutputBusInternal xInt if axis == Axis.x
          "Internal bus for measurements if face is normal to x axis"
          annotation (HideResult=not (internalMaterial or internalLin1 or
              internalLin1 or internalLin2 or internalEntropy), Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-10,-20}),iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,-98})));
        Connectors.RealOutputBusInternal yInt if axis == Axis.y
          "Internal bus for measurements if face is normal to y axis"
          annotation (HideResult=not (internalMaterial or internalLin1 or
              internalLin1 or internalLin2 or internalEntropy), Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,-40}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,-98})));
        Connectors.RealOutputBusInternal zInt if axis == Axis.z
          "Internal bus for measurements if face is normal to z axis"
          annotation (HideResult=not (internalMaterial or internalLin1 or
              internalLin1 or internalLin2 or internalEntropy), Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={10,-60}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,-98})));

      equation
        // Electrochemical potential
        connect(echemPot.material, face.material) annotation (Line(
            points={{-60,4},{-48,4},{-48,20},{5.55112e-16,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(y.mu, echemPot.y) annotation (Line(
            points={{5.55112e-16,-100},{0,-100},{0,-80},{-60,-80},{-60,-6}},
            color={0,0,127},
            smooth=Smooth.None));

        // 1st transverse velocity
        connect(velocity1.momentum, face.momentum1) annotation (Line(
            points={{-20,4},{-8,4},{-8,20},{5.55112e-16,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(velocity1.y, xInt.phi_y) annotation (Line(
            points={{-20,-6},{-20,-10},{-10,-10},{-10,-20}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(velocity1.y, yInt.phi_z) annotation (Line(
            points={{-20,-6},{-20,-30},{0,-30},{0,-40},{5.55112e-16,-40}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(velocity1.y, zInt.phi_x) annotation (Line(
            points={{-20,-6},{-20,-50},{10,-50},{10,-60}},
            color={0,0,127},
            smooth=Smooth.None));

        // 2nd transverse velocity
        connect(velocity2.momentum, face.momentum2) annotation (Line(
            points={{20,4},{8,4},{8,20},{5.55112e-16,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(velocity2.y, xInt.phi_z) annotation (Line(
            points={{20,-6},{20,-10},{-10,-10},{-10,-20}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(velocity2.y, yInt.phi_x) annotation (Line(
            points={{20,-6},{20,-30},{0,-30},{0,-40},{5.55112e-16,-40}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(velocity2.y, zInt.phi_y) annotation (Line(
            points={{20,-6},{20,-50},{10,-50},{10,-60}},
            color={0,0,127},
            smooth=Smooth.None));

        // Temperature
        connect(temperature.entropy, face.entropy) annotation (Line(
            points={{60,4},{48,4},{48,20},{5.55112e-16,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(y.T, temperature.y) annotation (Line(
            points={{5.55112e-16,-100},{5.55112e-16,-80},{60,-80},{60,-6}},
            color={0,0,127},
            smooth=Smooth.None), Text(
            string="%first",
            index=-1,
            extent={{-6,3},{-6,3}}));

        // Conditional axes
        connect(xInt, y) annotation (Line(
            points={{-10,-20},{-10,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(yInt, y) annotation (Line(
            points={{5.55112e-16,-40},{-4.87687e-22,-64},{5.55112e-16,-64},{
                5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(zInt, y) annotation (Line(
            points={{10,-60},{10,-70},{0,-70},{0,-100},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));

        annotation (
          defaultComponentName="speciesFaceSensor",
          Diagram(graphics),
          Icon(graphics));
      end Species;

      model PotentialElectrochemical "Sensor for electrochemical potential"

        extends FCSys.Sensors.BaseClasses.PartialSensor(redeclare
            Connectors.RealOutput y(final unit="l2.m/(N.T2)"));

        FCSys.Connectors.MaterialTransport material
          annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
      equation
        y = material.mu;
        0 = material.Ndot "Material rate balance (no storage)";
        annotation (Icon(graphics={Text(
                      extent={{-100,-20},{100,-50}},
                      lineColor={127,127,127},
                      textString="mu")}));
      end PotentialElectrochemical;

      model Velocity "Sensor for velocity"
        extends FCSys.Sensors.BaseClasses.PartialSensor(redeclare
            Connectors.RealOutput y(final unit="l/T"));
        FCSys.Connectors.MechanicalTransport momentum
          annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

      equation
        y = momentum.phi "Measurement";
        0 = momentum.mPhidot "Linear momentum rate balance (no storage)";
        annotation (Icon(graphics={Text(
                      extent={{-100,-20},{100,-50}},
                      lineColor={127,127,127},
                      textString="phi")}));
      end Velocity;

      model Temperature "Sensor for temperature"
        extends FCSys.Sensors.BaseClasses.PartialSensor(redeclare
            Connectors.RealOutput y(final unit="l2.m/(N.T2)", displayUnit="K"));
        FCSys.Connectors.Heat entropy "Entropy connector for the face"
          annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
      equation
        y = entropy.T "Measurement";
        0 = entropy.Qdot "Entropy rate balance (no storage)";
        annotation (Diagram(graphics), Icon(graphics={Text(
                      extent={{-100,-20},{100,-50}},
                      lineColor={127,127,127},
                      textString="T")}));
      end Temperature;
    end Species;
    annotation (Documentation(info="<html><p>Since the connectors in
<a href=\"modelica://FCSys\">FCSys</a> are hierarchical
(see the <a href=\"modelica://FCSys.Connectors\">Connectors</a> package),
the models for the sensors must be as well.  A
<a href=\"modelica://FCSys.Connectors.BaseClasses.PartialFace\">Face</a>
connector (<a href=\"modelica://FCSys.Connectors.FaceX\">FaceX</a>,
<a href=\"modelica://FCSys.Connectors.FaceY\">FaceY</a>,
<a href=\"modelica://FCSys.Connectors.FaceZ\">FaceZ</a>,
or <a href=\"modelica://FCSys.Connectors.FaceGeneric\">FaceGeneric</a>)
is used in <a href=\"modelica://FCSys.Subregions.Species\">Species</a> models,
and there is a corresponding <a href=\"modelica://FCSys.Sensors.Face.Species.Species\">Species
sensor</a> model in this package. The
<a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a>
connector is used in <a href=\"modelica://FCSys.Subregions.Phases\">Phase</a> models,
and there are corresponding <a href=\"modelica://FCSys.Sensors.Face.Phases\">Phase
sensor</a> models.  The
<a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connector is nested once more
in models such as the <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a>,
and there is a corresponding <a href=\"modelica://FCSys.Sensors.Face.Subregion\">Subregion
sensor</a> model.
</p></html>"));
  end Face;

  package FaceDifferential
    "<html>Sensors for pairs of <a href=\"modelica://FCSys.Connectors.BaseClasses.PartialFace\">Face</a> or <a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connectors</html>"

    extends Modelica.Icons.Package;

    model Subregion
      "<html>Sensor for a face of a <a href=\"modelica://FCSys.Subregions.Subregion\">Region</a> or <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a> model</html>"
      import FCSys.BaseClasses.Axis;
      extends FCSys.BaseClasses.Icons.Sensor;

      parameter Boolean effort=true "true, if effort sensor (otherwise, flow)"
        annotation (choices(__Dymola_checkBox=true), Dialog(compact=true));

      parameter FCSys.BaseClasses.Axis axis=Axis.x "Axis normal to the face";

      replaceable Phases.Gas gas(final effort=effort,final axis=axis) "Gas"
        annotation (Dialog(group="Phases", __Dymola_descriptionLabel=true),
          Placement(transformation(extent={{-10,-10},{10,10}})));

      replaceable Phases.Graphite graphite(final effort=effort,final axis=axis)
        "Graphite" annotation (Dialog(group="Phases", __Dymola_descriptionLabel
            =true), Placement(transformation(extent={{-10,-10},{10,10}})));

      replaceable Phases.Ionomer ionomer(final effort=effort,final axis=axis)
        "Ionomer" annotation (Dialog(group="Phases", __Dymola_descriptionLabel=
              true), Placement(transformation(extent={{-10,-10},{10,10}})));

      Connectors.RealOutputBus y "Bus of measurements" annotation (Placement(
            transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,-100})));
      Connectors.FaceBus negative
        "Negative-side connector for material, linear momentum, and entropy"
        annotation (Placement(transformation(extent={{-110,-10},{-90,10}}),
            iconTransformation(extent={{-110,-10},{-90,10}})));
      Connectors.FaceBus positive
        "Positive-side connector for material, linear momentum, and entropy"
        annotation (Placement(transformation(extent={{90,-10},{110,10}}),
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

      model Phase "Sensor for a phase with all species conditionally included"

        extends FCSys.BaseClasses.Icons.Sensor;

        extends FCSys.Sensors.FaceDifferential.Phases.BaseClasses.NullPhase;

        // Conditionally include species.
        parameter Boolean inclC=false "Carbon (C)" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FCSys.Sensors.FaceDifferential.Species.Species C(final effort=effort,
            final axis=axis) if inclC "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclC), Placement(transformation(extent={{-10,-10},{10,10}})));
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
        FCSys.Sensors.FaceDifferential.Species.Species C19HF37O5S(final effort=
              effort, final axis=axis) if inclC19HF37O5S "Model" annotation (
            Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclC19HF37O5S), Placement(transformation(extent={{-10,-10},
                  {10,10}})));
        parameter Boolean 'incle-'=false
          "<html>Electrons (e<sup>-</sup>)</html>" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FCSys.Sensors.FaceDifferential.Species.Species 'e-'(final effort=effort,
            final axis=axis) if 'incle-' "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable='incle-'), Placement(transformation(extent={{-10,-10},{10,10}})));
        parameter Boolean inclH2=false "<html>Hydrogen (H<sub>2</sub>)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FCSys.Sensors.FaceDifferential.Species.Species H2(final effort=effort,
            final axis=axis) if inclH2 "Model" annotation (Dialog(
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
        FCSys.Sensors.FaceDifferential.Species.Species H2O(final effort=effort,
            final axis=axis) if inclH2O "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclH2O), Placement(transformation(extent={{-10,-10},{10,10}})));
        parameter Boolean 'inclH+'=false "<html>Protons (H<sup>+</sup>)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FCSys.Sensors.FaceDifferential.Species.Species 'H+'(final effort=effort,
            final axis=axis) if 'inclH+' "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable='inclH+'), Placement(transformation(extent={{-10,-10},{10,10}})));
        parameter Boolean inclN2=false "<html>Nitrogen (N<sub>2</sub>)</html>"
          annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FCSys.Sensors.FaceDifferential.Species.Species N2(final effort=effort,
            final axis=axis) if inclN2 "Model" annotation (Dialog(
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
        FCSys.Sensors.FaceDifferential.Species.Species O2(final effort=effort,
            final axis=axis) if inclO2 "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclO2), Placement(transformation(extent={{-10,-10},{10,10}})));

      equation
        // Note:  It would be helpful if Modelica allowed elements of expandable
        // connectors to be named by the contents of a string variable and the
        // name of an instance of a model was accessible through a string (like
        // %name is expanded to be the name of the instance of the model).  Then,
        // the connection equations that follow could be generic.

        // C
        connect(C.negative.material, negative.C.material) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.positive.material, positive.C.material) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.negative.momentum1, xNegative.C.momentumY) annotation (Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.positive.momentum1, xPositive.C.momentumY) annotation (Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.negative.momentum1, yNegative.C.momentumZ) annotation (Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.positive.momentum1, yPositive.C.momentumZ) annotation (Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.negative.momentum1, zNegative.C.momentumX) annotation (Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.positive.momentum1, zPositive.C.momentumX) annotation (Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.negative.momentum2, xNegative.C.momentumZ) annotation (Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.positive.momentum2, xPositive.C.momentumZ) annotation (Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.negative.momentum2, yNegative.C.momentumX) annotation (Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.positive.momentum2, yPositive.C.momentumX) annotation (Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.negative.momentum2, zNegative.C.momentumY) annotation (Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.positive.momentum2, zPositive.C.momentumY) annotation (Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.negative.entropy, negative.C.entropy) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.positive.entropy, positive.C.entropy) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(y.C, C.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-9.8}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // C19HF37O5S
        connect(C19HF37O5S.negative.material, negative.C19HF37O5S.material)
          annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.positive.material, positive.C19HF37O5S.material)
          annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.negative.momentum1, xNegative.C19HF37O5S.momentumY)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.positive.momentum1, xPositive.C19HF37O5S.momentumY)
          annotation (Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.negative.momentum1, yNegative.C19HF37O5S.momentumZ)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.positive.momentum1, yPositive.C19HF37O5S.momentumZ)
          annotation (Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.negative.momentum1, zNegative.C19HF37O5S.momentumX)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.positive.momentum1, zPositive.C19HF37O5S.momentumX)
          annotation (Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.negative.momentum2, xNegative.C19HF37O5S.momentumZ)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.positive.momentum2, xPositive.C19HF37O5S.momentumZ)
          annotation (Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.negative.momentum2, yNegative.C19HF37O5S.momentumX)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.positive.momentum2, yPositive.C19HF37O5S.momentumX)
          annotation (Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.negative.momentum2, zNegative.C19HF37O5S.momentumY)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.positive.momentum2, zPositive.C19HF37O5S.momentumY)
          annotation (Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.negative.entropy, negative.C19HF37O5S.entropy)
          annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.positive.entropy, positive.C19HF37O5S.entropy)
          annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(y.C19HF37O5S, C19HF37O5S.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-9.8}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // e-
        connect('e-'.negative.material, negative.'e-'.material) annotation (
            Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.positive.material, positive.'e-'.material) annotation (
            Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.negative.momentum1, xNegative.'e-'.momentumY) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.positive.momentum1, xPositive.'e-'.momentumY) annotation (
            Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.negative.momentum1, yNegative.'e-'.momentumZ) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.positive.momentum1, yPositive.'e-'.momentumZ) annotation (
            Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.negative.momentum1, zNegative.'e-'.momentumX) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.positive.momentum1, zPositive.'e-'.momentumX) annotation (
            Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.negative.momentum2, xNegative.'e-'.momentumZ) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.positive.momentum2, xPositive.'e-'.momentumZ) annotation (
            Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.negative.momentum2, yNegative.'e-'.momentumX) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.positive.momentum2, yPositive.'e-'.momentumX) annotation (
            Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.negative.momentum2, zNegative.'e-'.momentumY) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.positive.momentum2, zPositive.'e-'.momentumY) annotation (
            Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.negative.entropy, negative.'e-'.entropy) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.positive.entropy, positive.'e-'.entropy) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(y.'e-', 'e-'.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-9.8}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // H2
        connect(H2.negative.material, negative.H2.material) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.positive.material, positive.H2.material) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.negative.momentum1, xNegative.H2.momentumY) annotation (Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.positive.momentum1, xPositive.H2.momentumY) annotation (Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.negative.momentum1, yNegative.H2.momentumZ) annotation (Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.positive.momentum1, yPositive.H2.momentumZ) annotation (Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.negative.momentum1, zNegative.H2.momentumX) annotation (Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.positive.momentum1, zPositive.H2.momentumX) annotation (Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.negative.momentum2, xNegative.H2.momentumZ) annotation (Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.positive.momentum2, xPositive.H2.momentumZ) annotation (Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.negative.momentum2, yNegative.H2.momentumX) annotation (Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.positive.momentum2, yPositive.H2.momentumX) annotation (Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.negative.momentum2, zNegative.H2.momentumY) annotation (Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.positive.momentum2, zPositive.H2.momentumY) annotation (Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.negative.entropy, negative.H2.entropy) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.positive.entropy, positive.H2.entropy) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(y.H2, H2.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-9.8}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // H2O
        connect(H2O.negative.material, negative.H2O.material) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive.material, positive.H2O.material) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.negative.momentum1, xNegative.H2O.momentumY) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive.momentum1, xPositive.H2O.momentumY) annotation (
            Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.negative.momentum1, yNegative.H2O.momentumZ) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive.momentum1, yPositive.H2O.momentumZ) annotation (
            Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.negative.momentum1, zNegative.H2O.momentumX) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive.momentum1, zPositive.H2O.momentumX) annotation (
            Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.negative.momentum2, xNegative.H2O.momentumZ) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive.momentum2, xPositive.H2O.momentumZ) annotation (
            Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.negative.momentum2, yNegative.H2O.momentumX) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive.momentum2, yPositive.H2O.momentumX) annotation (
            Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.negative.momentum2, zNegative.H2O.momentumY) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive.momentum2, zPositive.H2O.momentumY) annotation (
            Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.negative.entropy, negative.H2O.entropy) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive.entropy, positive.H2O.entropy) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(y.H2O, H2O.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-9.8}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // H+
        connect('H+'.negative.material, negative.'H+'.material) annotation (
            Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.positive.material, positive.'H+'.material) annotation (
            Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.negative.momentum1, xNegative.'H+'.momentumY) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.positive.momentum1, xPositive.'H+'.momentumY) annotation (
            Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.negative.momentum1, yNegative.'H+'.momentumZ) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.positive.momentum1, yPositive.'H+'.momentumZ) annotation (
            Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.negative.momentum1, zNegative.'H+'.momentumX) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.positive.momentum1, zPositive.'H+'.momentumX) annotation (
            Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.negative.momentum2, xNegative.'H+'.momentumZ) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.positive.momentum2, xPositive.'H+'.momentumZ) annotation (
            Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.negative.momentum2, yNegative.'H+'.momentumX) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.positive.momentum2, yPositive.'H+'.momentumX) annotation (
            Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.negative.momentum2, zNegative.'H+'.momentumY) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.positive.momentum2, zPositive.'H+'.momentumY) annotation (
            Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.negative.entropy, negative.'H+'.entropy) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.positive.entropy, positive.'H+'.entropy) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(y.'H+', 'H+'.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-9.8}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // N2
        connect(N2.negative.material, negative.N2.material) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.positive.material, positive.N2.material) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.negative.momentum1, xNegative.N2.momentumY) annotation (Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.positive.momentum1, xPositive.N2.momentumY) annotation (Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.negative.momentum1, yNegative.N2.momentumZ) annotation (Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.positive.momentum1, yPositive.N2.momentumZ) annotation (Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.negative.momentum1, zNegative.N2.momentumX) annotation (Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.positive.momentum1, zPositive.N2.momentumX) annotation (Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.negative.momentum2, xNegative.N2.momentumZ) annotation (Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.positive.momentum2, xPositive.N2.momentumZ) annotation (Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.negative.momentum2, yNegative.N2.momentumX) annotation (Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.positive.momentum2, yPositive.N2.momentumX) annotation (Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.negative.momentum2, zNegative.N2.momentumY) annotation (Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.positive.momentum2, zPositive.N2.momentumY) annotation (Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.negative.entropy, negative.N2.entropy) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.positive.entropy, positive.N2.entropy) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(y.N2, N2.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-9.8}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // O2
        connect(O2.negative.material, negative.O2.material) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.positive.material, positive.O2.material) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.negative.momentum1, xNegative.O2.momentumY) annotation (Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.positive.momentum1, xPositive.O2.momentumY) annotation (Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.negative.momentum1, yNegative.O2.momentumZ) annotation (Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.positive.momentum1, yPositive.O2.momentumZ) annotation (Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.negative.momentum1, zNegative.O2.momentumX) annotation (Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.positive.momentum1, zPositive.O2.momentumX) annotation (Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.negative.momentum2, xNegative.O2.momentumZ) annotation (Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.positive.momentum2, xPositive.O2.momentumZ) annotation (Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.negative.momentum2, yNegative.O2.momentumX) annotation (Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.positive.momentum2, yPositive.O2.momentumX) annotation (Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.negative.momentum2, zNegative.O2.momentumY) annotation (Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.positive.momentum2, zPositive.O2.momentumY) annotation (Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.negative.entropy, negative.O2.entropy) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.positive.entropy, positive.O2.entropy) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(y.O2, O2.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-9.8}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="phaseFaceSensor",
          Diagram(graphics));
      end Phase;

      model Gas "Sensor for gas"

        extends FCSys.Sensors.FaceDifferential.Phases.BaseClasses.NullPhase;

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
        FCSys.Sensors.FaceDifferential.Species.Species H2(final effort=effort,
            final axis=axis) if inclH2 "Model" annotation (Dialog(
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
        FCSys.Sensors.FaceDifferential.Species.Species H2O(final effort=effort,
            final axis=axis) if inclH2O "Model" annotation (Dialog(
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
        FCSys.Sensors.FaceDifferential.Species.Species N2(final effort=effort,
            final axis=axis) if inclN2 "Model" annotation (Dialog(
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
        FCSys.Sensors.FaceDifferential.Species.Species O2(final effort=effort,
            final axis=axis) if inclO2 "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable=inclO2), Placement(transformation(extent={{-10,-10},{10,10}})));

      equation
        // H2
        connect(H2.negative.material, negative.H2.material) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.positive.material, positive.H2.material) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.negative.momentum1, xNegative.H2.momentumY) annotation (Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.positive.momentum1, xPositive.H2.momentumY) annotation (Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.negative.momentum1, yNegative.H2.momentumZ) annotation (Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.positive.momentum1, yPositive.H2.momentumZ) annotation (Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.negative.momentum1, zNegative.H2.momentumX) annotation (Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.positive.momentum1, zPositive.H2.momentumX) annotation (Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.negative.momentum2, xNegative.H2.momentumZ) annotation (Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.positive.momentum2, xPositive.H2.momentumZ) annotation (Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.negative.momentum2, yNegative.H2.momentumX) annotation (Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.positive.momentum2, yPositive.H2.momentumX) annotation (Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.negative.momentum2, zNegative.H2.momentumY) annotation (Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.positive.momentum2, zPositive.H2.momentumY) annotation (Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.negative.entropy, negative.H2.entropy) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2.positive.entropy, positive.H2.entropy) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(y.H2, H2.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-9.8}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // H2O
        connect(H2O.negative.material, negative.H2O.material) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive.material, positive.H2O.material) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.negative.momentum1, xNegative.H2O.momentumY) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive.momentum1, xPositive.H2O.momentumY) annotation (
            Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.negative.momentum1, yNegative.H2O.momentumZ) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive.momentum1, yPositive.H2O.momentumZ) annotation (
            Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.negative.momentum1, zNegative.H2O.momentumX) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive.momentum1, zPositive.H2O.momentumX) annotation (
            Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.negative.momentum2, xNegative.H2O.momentumZ) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive.momentum2, xPositive.H2O.momentumZ) annotation (
            Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.negative.momentum2, yNegative.H2O.momentumX) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive.momentum2, yPositive.H2O.momentumX) annotation (
            Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.negative.momentum2, zNegative.H2O.momentumY) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive.momentum2, zPositive.H2O.momentumY) annotation (
            Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.negative.entropy, negative.H2O.entropy) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(H2O.positive.entropy, positive.H2O.entropy) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(y.H2O, H2O.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-9.8}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // N2
        connect(N2.negative.material, negative.N2.material) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.positive.material, positive.N2.material) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.negative.momentum1, xNegative.N2.momentumY) annotation (Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.positive.momentum1, xPositive.N2.momentumY) annotation (Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.negative.momentum1, yNegative.N2.momentumZ) annotation (Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.positive.momentum1, yPositive.N2.momentumZ) annotation (Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.negative.momentum1, zNegative.N2.momentumX) annotation (Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.positive.momentum1, zPositive.N2.momentumX) annotation (Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.negative.momentum2, xNegative.N2.momentumZ) annotation (Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.positive.momentum2, xPositive.N2.momentumZ) annotation (Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.negative.momentum2, yNegative.N2.momentumX) annotation (Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.positive.momentum2, yPositive.N2.momentumX) annotation (Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.negative.momentum2, zNegative.N2.momentumY) annotation (Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.positive.momentum2, zPositive.N2.momentumY) annotation (Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.negative.entropy, negative.N2.entropy) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(N2.positive.entropy, positive.N2.entropy) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(y.N2, N2.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-9.8}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // O2
        connect(O2.negative.material, negative.O2.material) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.positive.material, positive.O2.material) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.negative.momentum1, xNegative.O2.momentumY) annotation (Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.positive.momentum1, xPositive.O2.momentumY) annotation (Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.negative.momentum1, yNegative.O2.momentumZ) annotation (Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.positive.momentum1, yPositive.O2.momentumZ) annotation (Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.negative.momentum1, zNegative.O2.momentumX) annotation (Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.positive.momentum1, zPositive.O2.momentumX) annotation (Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.negative.momentum2, xNegative.O2.momentumZ) annotation (Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.positive.momentum2, xPositive.O2.momentumZ) annotation (Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.negative.momentum2, yNegative.O2.momentumX) annotation (Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.positive.momentum2, yPositive.O2.momentumX) annotation (Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.negative.momentum2, zNegative.O2.momentumY) annotation (Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.positive.momentum2, zPositive.O2.momentumY) annotation (Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.negative.entropy, negative.O2.entropy) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(O2.positive.entropy, positive.O2.entropy) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(y.O2, O2.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-9.8}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));
        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="phaseFaceSensor");
      end Gas;

      model Graphite "Sensor for graphite"

        extends FCSys.Sensors.FaceDifferential.Phases.BaseClasses.NullPhase;

        // Conditionally include species.
        parameter Boolean inclC=false "Carbon (C)" annotation (
          Evaluate=true,
          HideResult=true,
          choices(__Dymola_checkBox=true),
          Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        FCSys.Sensors.FaceDifferential.Species.Species C(final effort=effort,
            final axis=axis) if inclC "Model" annotation (Dialog(
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
        FCSys.Sensors.FaceDifferential.Species.Species 'e-'(final effort=effort,
            final axis=axis) if 'incle-' "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable='incle-'), Placement(transformation(extent={{-10,-10},{10,10}})));

      equation
        // C
        connect(C.negative.material, negative.C.material) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.positive.material, positive.C.material) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.negative.momentum1, xNegative.C.momentumY) annotation (Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.positive.momentum1, xPositive.C.momentumY) annotation (Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.negative.momentum1, yNegative.C.momentumZ) annotation (Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.positive.momentum1, yPositive.C.momentumZ) annotation (Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.negative.momentum1, zNegative.C.momentumX) annotation (Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.positive.momentum1, zPositive.C.momentumX) annotation (Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.negative.momentum2, xNegative.C.momentumZ) annotation (Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.positive.momentum2, xPositive.C.momentumZ) annotation (Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.negative.momentum2, yNegative.C.momentumX) annotation (Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.positive.momentum2, yPositive.C.momentumX) annotation (Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.negative.momentum2, zNegative.C.momentumY) annotation (Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.positive.momentum2, zPositive.C.momentumY) annotation (Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.negative.entropy, negative.C.entropy) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C.positive.entropy, positive.C.entropy) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(y.C, C.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-9.8}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        // e-
        connect('e-'.negative.material, negative.'e-'.material) annotation (
            Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.positive.material, positive.'e-'.material) annotation (
            Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.negative.momentum1, xNegative.'e-'.momentumY) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.positive.momentum1, xPositive.'e-'.momentumY) annotation (
            Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.negative.momentum1, yNegative.'e-'.momentumZ) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.positive.momentum1, yPositive.'e-'.momentumZ) annotation (
            Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.negative.momentum1, zNegative.'e-'.momentumX) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.positive.momentum1, zPositive.'e-'.momentumX) annotation (
            Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.negative.momentum2, xNegative.'e-'.momentumZ) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.positive.momentum2, xPositive.'e-'.momentumZ) annotation (
            Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.negative.momentum2, yNegative.'e-'.momentumX) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.positive.momentum2, yPositive.'e-'.momentumX) annotation (
            Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.negative.momentum2, zNegative.'e-'.momentumY) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.positive.momentum2, zPositive.'e-'.momentumY) annotation (
            Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.negative.entropy, negative.'e-'.entropy) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('e-'.positive.entropy, positive.'e-'.entropy) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(y.'e-', 'e-'.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-9.8}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        annotation (defaultComponentPrefixes="replaceable",
            defaultComponentName="phaseFaceSensor");
      end Graphite;

      model Ionomer "Sensor for ionomer"

        extends FCSys.Sensors.FaceDifferential.Phases.BaseClasses.NullPhase;

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
        FCSys.Sensors.FaceDifferential.Species.Species C19HF37O5S(final effort=
              effort, final axis=axis) if inclC19HF37O5S "Model" annotation (
            Dialog(
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
        FCSys.Sensors.FaceDifferential.Species.Species 'H+'(final effort=effort,
            final axis=axis) if 'inclH+' "Model" annotation (Dialog(
            group="Species",
            __Dymola_descriptionLabel=true,
            enable='inclH+'), Placement(transformation(extent={{-10,-10},{10,10}})));

      equation
        // C19HF37O5S
        connect(C19HF37O5S.negative.material, negative.C19HF37O5S.material)
          annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.positive.material, positive.C19HF37O5S.material)
          annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.negative.momentum1, xNegative.C19HF37O5S.momentumY)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.positive.momentum1, xPositive.C19HF37O5S.momentumY)
          annotation (Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.negative.momentum1, yNegative.C19HF37O5S.momentumZ)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.positive.momentum1, yPositive.C19HF37O5S.momentumZ)
          annotation (Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.negative.momentum1, zNegative.C19HF37O5S.momentumX)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.positive.momentum1, zPositive.C19HF37O5S.momentumX)
          annotation (Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.negative.momentum2, xNegative.C19HF37O5S.momentumZ)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.positive.momentum2, xPositive.C19HF37O5S.momentumZ)
          annotation (Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.negative.momentum2, yNegative.C19HF37O5S.momentumX)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.positive.momentum2, yPositive.C19HF37O5S.momentumX)
          annotation (Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.negative.momentum2, zNegative.C19HF37O5S.momentumY)
          annotation (Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.positive.momentum2, zPositive.C19HF37O5S.momentumY)
          annotation (Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.negative.entropy, negative.C19HF37O5S.entropy)
          annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(C19HF37O5S.positive.entropy, positive.C19HF37O5S.entropy)
          annotation (Line(
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
        connect('H+'.negative.material, negative.'H+'.material) annotation (
            Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.positive.material, positive.'H+'.material) annotation (
            Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.negative.momentum1, xNegative.'H+'.momentumY) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.positive.momentum1, xPositive.'H+'.momentumY) annotation (
            Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.negative.momentum1, yNegative.'H+'.momentumZ) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.positive.momentum1, yPositive.'H+'.momentumZ) annotation (
            Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.negative.momentum1, zNegative.'H+'.momentumX) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.positive.momentum1, zPositive.'H+'.momentumX) annotation (
            Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.negative.momentum2, xNegative.'H+'.momentumZ) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.positive.momentum2, xPositive.'H+'.momentumZ) annotation (
            Line(
            points={{10,6.10623e-16},{60,20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.negative.momentum2, yNegative.'H+'.momentumX) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.positive.momentum2, yPositive.'H+'.momentumX) annotation (
            Line(
            points={{10,6.10623e-16},{60,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.negative.momentum2, zNegative.'H+'.momentumY) annotation (
            Line(
            points={{-10,6.10623e-16},{-60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.positive.momentum2, zPositive.'H+'.momentumY) annotation (
            Line(
            points={{10,6.10623e-16},{60,-20}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.negative.entropy, negative.'H+'.entropy) annotation (Line(
            points={{-10,6.10623e-16},{-100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect('H+'.positive.entropy, positive.'H+'.entropy) annotation (Line(
            points={{10,6.10623e-16},{100,5.55112e-16}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(y.'H+', 'H+'.y) annotation (Line(
            points={{5.55112e-16,-100},{6.10623e-16,-9.8}},
            color={0,0,127},
            thickness=0.5,
            smooth=Smooth.None));

        annotation (
          defaultComponentPrefixes="replaceable",
          defaultComponentName="phaseFaceSensor",
          Diagram(graphics));
      end Ionomer;

      package BaseClasses "Base classes (not for direct use)"
        extends Modelica.Icons.BasesPackage;

        model NullPhase "Empty sensor for a phase (no species)"
          import FCSys.BaseClasses.Axis;
          extends FCSys.BaseClasses.Icons.Sensor;

          parameter Boolean effort=true
            "true, if effort sensor (otherwise, flow)"
            annotation (choices(__Dymola_checkBox=true), Dialog(compact=true));

          parameter FCSys.BaseClasses.Axis axis=Axis.x
            "Axis normal to the face";

          FCSys.Connectors.FaceBus negative
            "Negative-side connector for material, linear momentum, and entropy"
            annotation (Placement(transformation(extent={{-110,-10},{-90,10}}),
                iconTransformation(extent={{-110,-10},{-90,10}})));
          FCSys.Connectors.FaceBus positive
            "Positive-side connector for material, linear momentum, and entropy"
            annotation (Placement(transformation(extent={{90,-10},{110,10}}),
                iconTransformation(extent={{90,-10},{110,10}})));
          Connectors.RealOutputBus y "Bus of measurements" annotation (
              Placement(transformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={0,-100}), iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={0,-100})));

        protected
          Connectors.FaceBusInternal xNegative if axis == Axis.x
            "Internal negative-side connector enabled if x axis" annotation (
              Placement(transformation(extent={{-70,10},{-50,30}}),
                iconTransformation(extent={{-30,-40},{-10,-20}})));
          Connectors.FaceBusInternal yNegative if axis == Axis.y
            "Internal negative-side connector enabled if y axis" annotation (
              Placement(transformation(extent={{-70,-10},{-50,10}}),
                iconTransformation(extent={{-10,-40},{10,-20}})));
          Connectors.FaceBusInternal zNegative if axis == Axis.z
            "Internal negative-side connector enabled if z axis" annotation (
              Placement(transformation(extent={{-70,-30},{-50,-10}}),
                iconTransformation(extent={{10,-40},{30,-20}})));
          Connectors.FaceBusInternal xPositive if axis == Axis.x
            "Internal positive-side connector enabled if x axis" annotation (
              Placement(transformation(extent={{50,10},{70,30}}),
                iconTransformation(extent={{-30,-40},{-10,-20}})));
          Connectors.FaceBusInternal yPositive if axis == Axis.y
            "Internal positive-side connector enabled if y axis" annotation (
              Placement(transformation(extent={{50,-10},{70,10}}),
                iconTransformation(extent={{-10,-40},{10,-20}})));
          Connectors.FaceBusInternal zPositive if axis == Axis.z
            "Internal positive-side connector enabled if z axis" annotation (
              Placement(transformation(extent={{50,-30},{70,-10}}),
                iconTransformation(extent={{10,-40},{30,-20}})));
        equation

          connect(xNegative, negative) annotation (Line(
              points={{-60,20},{-100,5.55112e-16}},
              color={127,127,127},
              thickness=0.5,
              smooth=Smooth.None));

          connect(yNegative, negative) annotation (Line(
              points={{-60,5.55112e-16},{-70,5.55112e-16},{-70,5.55112e-16},{-80,
                  5.55112e-16},{-80,5.55112e-16},{-100,5.55112e-16}},
              color={127,127,127},
              thickness=0.5,
              smooth=Smooth.None));

          connect(zNegative, negative) annotation (Line(
              points={{-60,-20},{-100,5.55112e-16}},
              color={127,127,127},
              thickness=0.5,
              smooth=Smooth.None));

          connect(xPositive, positive) annotation (Line(
              points={{60,20},{100,5.55112e-16}},
              color={127,127,127},
              thickness=0.5,
              smooth=Smooth.None));

          connect(yPositive, positive) annotation (Line(
              points={{60,5.55112e-16},{70,5.55112e-16},{70,5.55112e-16},{80,
                  5.55112e-16},{80,5.55112e-16},{100,5.55112e-16}},
              color={127,127,127},
              thickness=0.5,
              smooth=Smooth.None));

          connect(zPositive, positive) annotation (Line(
              points={{60,-20},{100,5.55112e-16}},
              color={127,127,127},
              thickness=0.5,
              smooth=Smooth.None));

          annotation (
            defaultComponentPrefixes="replaceable",
            defaultComponentName="phaseFaceSensor",
            Diagram(graphics),
            Icon(graphics={Line(
                          points={{-70,0},{-100,0}},
                          color={127,127,127},
                          smooth=Smooth.None),Line(
                          points={{100,0},{70,0}},
                          color={127,127,127},
                          smooth=Smooth.None)}));
        end NullPhase;
      end BaseClasses;
    end Phases;

    package Species
      "<html>Sensors for a single <a href=\"modelica://FCSys.Connectors.BaseClasses.BaseClasses.PartialFace\">Face</a> connector, e.g., of a <a href=\"modelica://FCSys.Subregions.Species\">Species</a> model (single-species)</html>"
      extends Modelica.Icons.Package;

      model Species
        "<html>Sensor for faces of a <a href=\"modelica://FCSys.Subregions.Species\">Species</a> model (single-species)</html>"
        import FCSys.BaseClasses.Axis;

        extends FCSys.BaseClasses.Icons.Sensor;

        parameter Boolean effort=true
          "true, if effort sensor (otherwise, flow)"
          annotation (choices(__Dymola_checkBox=true), Dialog(compact=true));

        parameter FCSys.BaseClasses.Axis axis=Axis.x "Axis normal to the face";

        parameter ThermoOpt thermoOpt=ThermoOpt.ClosedAdiabatic
          "Options for material and entropy subconnectors"
          annotation (Dialog(compact=true));

        // Material
        final parameter Boolean open=thermoOpt == ThermoOpt.OpenDiabatic "Open"
          annotation (choices(__Dymola_checkBox=true), Dialog(
            group="Material",
            compact=true,
            __Dymola_label="Included",
            __Dymola_descriptionLabel=true,
            __Dymola_joinNext=true));
        Material.PotentialElectrochemical echemPot if open and effort
          annotation (Placement(transformation(extent={{-90,64},{-70,84}})));
        Material.Current current if open and not effort
          annotation (Placement(transformation(extent={{-70,48},{-50,68}})));

        // 1st transverse linear momentum
        parameter Boolean slip1=false
          "<html>Viscous (1<sup>st</sup> transverse momentum included)</html>"
          annotation (choices(__Dymola_checkBox=true), Dialog(compact=true));
        Momentum.Velocity velocity1 if slip1 and effort
          annotation (Placement(transformation(extent={{-50,34},{-30,54}})));
        Momentum.Force force1 if slip1 and not effort
          annotation (Placement(transformation(extent={{-30,18},{-10,38}})));

        // 2nd transverse linear momentum
        parameter Boolean slip2=false
          "<html>Viscous (2<sup>nd</sup> transverse momentum included)</html>"
          annotation (choices(__Dymola_checkBox=true), Dialog(compact=true));
        Momentum.Velocity velocity2 if slip2 and effort
          annotation (Placement(transformation(extent={{-10,4},{10,24}})));
        Momentum.Force force2 if slip2 and not effort
          annotation (Placement(transformation(extent={{10,-10},{30,10}})));

        // Entropy
        final parameter Boolean diabatic=thermoOpt == ThermoOpt.ClosedDiabatic
             or thermoOpt == ThermoOpt.OpenDiabatic
          "Diabatic (entropy included)"
          annotation (choices(__Dymola_checkBox=true), Dialog(compact=true));
        Entropy.Temperature temperature if diabatic and effort
          annotation (Placement(transformation(extent={{30,-24},{50,-4}})));
        Entropy.HeatRate heatRate if diabatic and not effort
          annotation (Placement(transformation(extent={{70,-52},{90,-32}})));
        Entropy.EntropyRate entropyRate if diabatic and not effort
          annotation (Placement(transformation(extent={{50,-38},{70,-18}})));

        FCSys.Connectors.FaceGeneric negative(
          final thermoOpt=thermoOpt,
          final slip1=slip1,
          final slip2=slip2)
          "Negative-side connector for material, linear momentum, and entropy"
          annotation (Placement(transformation(extent={{-110,-10},{-90,10}}),
              iconTransformation(extent={{-110,-10},{-90,10}})));
        Connectors.RealOutputBus y "Output bus for measurements" annotation (
            HideResult=not (internalMaterial or internalLin1 or internalLin1
               or internalLin2 or internalEntropy), Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,-100}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,-98})));

        FCSys.Connectors.FaceGeneric positive(
          final thermoOpt=thermoOpt,
          final slip1=slip1,
          final slip2=slip2)
          "Positive-side connector for material, linear momentum, and entropy"
          annotation (Placement(transformation(extent={{90,-2},{110,18}}),
              iconTransformation(extent={{90,-10},{110,10}})));
      protected
        Connectors.RealOutputBusInternal xInt if axis == Axis.x
          "Internal bus for measurements if face is normal to x axis"
          annotation (HideResult=not (internalMaterial or internalLin1 or
              internalLin1 or internalLin2 or internalEntropy), Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-34,-60}),iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,-98})));
        Connectors.RealOutputBusInternal yInt if axis == Axis.y
          "Internal bus for measurements if face is normal to y axis"
          annotation (HideResult=not (internalMaterial or internalLin1 or
              internalLin1 or internalLin2 or internalEntropy), Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-20,-60}),iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,-98})));
        Connectors.RealOutputBusInternal zInt if axis == Axis.z
          "Internal bus for measurements if face is normal to z axis"
          annotation (HideResult=not (internalMaterial or internalLin1 or
              internalLin1 or internalLin2 or internalEntropy), Placement(
              transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-8,-60}), iconTransformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,-98})));
      equation
        // Electrochemical potential
        connect(echemPot.negative, negative.material) annotation (Line(
            points={{-90,74},{-90,5.55112e-16},{-100,5.55112e-16}},
            color={127,127,127},
            smooth=Smooth.None));
        connect(echemPot.positive, positive.material) annotation (Line(
            points={{-70,74},{90,74},{90,8},{100,8}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(echemPot.y, y.Deltamu) annotation (Line(
            points={{-80,64},{-80,-70},{5.55112e-16,-70},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(current.negative, negative.material) annotation (Line(
            points={{-70,58},{-90,58},{-90,5.55112e-16},{-100,5.55112e-16}},
            color={127,127,127},
            smooth=Smooth.None));
        connect(current.positive, positive.material) annotation (Line(
            points={{-50,58},{90,58},{90,8},{100,8}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(current.y, y.Ndot) annotation (Line(
            points={{-60,48},{-60,-70},{5.55112e-16,-70},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));

        // 1st transverse velocity
        connect(velocity1.negative, negative.momentum1) annotation (Line(
            points={{-50,44},{-90,44},{-90,5.55112e-16},{-100,5.55112e-16}},
            color={127,127,127},
            smooth=Smooth.None));
        connect(velocity1.positive, positive.momentum1) annotation (Line(
            points={{-30,44},{90,44},{90,8},{100,8}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(velocity1.y, xInt.Deltaphi_y) annotation (Line(
            points={{-40,34},{-40,-20},{-34,-20},{-34,-60}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(velocity1.y, yInt.Deltaphi_z) annotation (Line(
            points={{-40,34},{-40,-34},{-20,-34},{-20,-60}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(velocity1.y, zInt.Deltaphi_x) annotation (Line(
            points={{-40,34},{-40,-48},{-8,-48},{-8,-60}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(force1.negative, negative.momentum1) annotation (Line(
            points={{-30,28},{-90,28},{-90,5.55112e-16},{-100,5.55112e-16}},
            color={127,127,127},
            smooth=Smooth.None));
        connect(force1.positive, positive.momentum1) annotation (Line(
            points={{-10,28},{90,28},{90,8},{100,8}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(force1.y, xInt.mPhidot_y) annotation (Line(
            points={{-20,18},{-20,-20},{-34,-20},{-34,-60}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(force1.y, yInt.mPhidot_z) annotation (Line(
            points={{-20,18},{-20,-60}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(force1.y, zInt.mPhidot_x) annotation (Line(
            points={{-20,18},{-20,-48},{-8,-48},{-8,-60}},
            color={0,0,127},
            smooth=Smooth.None));

        // 2nd transverse velocity
        connect(velocity2.negative, negative.momentum2) annotation (Line(
            points={{-10,14},{-90,14},{-90,5.55112e-16},{-100,5.55112e-16}},
            color={127,127,127},
            smooth=Smooth.None));
        connect(velocity2.positive, positive.momentum2) annotation (Line(
            points={{10,14},{90,14},{90,8},{100,8}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(velocity2.y, xInt.Deltaphi_z) annotation (Line(
            points={{6.10623e-16,4},{0,4},{0,-20},{-34,-20},{-34,-60}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(velocity2.y, yInt.Deltaphi_x) annotation (Line(
            points={{6.10623e-16,4},{0,4},{0,-34},{-20,-34},{-20,-60}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(velocity2.y, zInt.Deltaphi_y) annotation (Line(
            points={{6.10623e-16,4},{0,4},{0,-48},{-8,-48},{-8,-60}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(force2.negative, negative.momentum2) annotation (Line(
            points={{10,6.10623e-16},{-90,6.10623e-16},{-90,5.55112e-16},{-100,
                5.55112e-16}},
            color={127,127,127},
            smooth=Smooth.None));

        connect(force2.positive, positive.momentum2) annotation (Line(
            points={{30,6.10623e-16},{48,6.10623e-16},{48,0},{90,0},{90,8},{100,
                8}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));

        connect(force2.y, xInt.mPhidot_z) annotation (Line(
            points={{20,-10},{20,-20},{-34,-20},{-34,-60}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(force2.y, yInt.mPhidot_x) annotation (Line(
            points={{20,-10},{20,-34},{-20,-34},{-20,-60}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(force2.y, zInt.mPhidot_y) annotation (Line(
            points={{20,-10},{20,-48},{-8,-48},{-8,-60}},
            color={0,0,127},
            smooth=Smooth.None));

        // Conditional axes
        connect(xInt, y) annotation (Line(
            points={{-34,-60},{-34,-70},{5.55112e-16,-70},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(yInt, y) annotation (Line(
            points={{-20,-60},{-20,-70},{5.55112e-16,-70},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(zInt, y) annotation (Line(
            points={{-8,-60},{-8,-70},{5.55112e-16,-70},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));

        // Temperature
        connect(temperature.negative, negative.entropy) annotation (Line(
            points={{30,-14},{-90,-14},{-90,5.55112e-16},{-100,5.55112e-16}},
            color={127,127,127},
            smooth=Smooth.None));

        connect(temperature.positive, positive.entropy) annotation (Line(
            points={{50,-14},{90,-14},{90,8},{100,8}},
            color={127,127,127},
            smooth=Smooth.None));
        connect(temperature.y, y.DeltaT) annotation (Line(
            points={{40,-24},{40,-70},{5.55112e-16,-70},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));

        connect(entropyRate.negative, negative.entropy) annotation (Line(
            points={{50,-28},{-90,-28},{-90,5.55112e-16},{-100,5.55112e-16}},
            color={127,127,127},
            smooth=Smooth.None));
        connect(entropyRate.positive, positive.entropy) annotation (Line(
            points={{70,-28},{90,-28},{90,8},{100,8}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(entropyRate.y, y.Qdot_avg) annotation (Line(
            points={{60,-38},{60,-70},{5.55112e-16,-70},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(heatRate.negative, negative.entropy) annotation (Line(
            points={{70,-42},{-90,-42},{-90,5.55112e-16},{-100,5.55112e-16}},
            color={127,127,127},
            smooth=Smooth.None));
        connect(heatRate.positive, positive.entropy) annotation (Line(
            points={{90,-42},{90,8},{100,8}},
            color={127,127,127},
            pattern=LinePattern.None,
            smooth=Smooth.None));
        connect(heatRate.y, y.Qdot) annotation (Line(
            points={{80,-52},{80,-70},{5.55112e-16,-70},{5.55112e-16,-100}},
            color={0,0,127},
            smooth=Smooth.None));

        annotation (
          defaultComponentName="speciesFaceSensor",
          Diagram(graphics),
          Icon(graphics={Line(
                      points={{-70,0},{-100,0}},
                      color={0,0,0},
                      smooth=Smooth.None),Line(
                      points={{100,0},{70,0}},
                      color={0,0,0},
                      smooth=Smooth.None)}));
      end Species;

      package Material "Relative sensors for material"
        extends Modelica.Icons.Package;
        model PotentialElectrochemical
          "Sensor for difference in electrochemical potential"

          extends BaseClasses.PartialSensor(redeclare Connectors.RealOutput y(
                final unit="l2.m/(N.T2)"));

        equation
          y = negative.mu - positive.mu "Measurement";
          0 = negative.Ndot "Condition of no current";
          // Note:  In conjunction with the momentum rate balance, this means
          // that there's no current into either face.
          annotation (Icon(graphics={Text(
                          extent={{-100,-20},{100,-50}},
                          lineColor={127,127,127},
                          textString="  mu"),Polygon(
                          points={{-26,-26},{-18,-44},{-34,-44},{-26,-26}},
                          lineColor={127,127,127},
                          smooth=Smooth.None,
                          lineThickness=0.5)}), Diagram(graphics));
        end PotentialElectrochemical;

        model Current "Sensor for current"
          extends BaseClasses.PartialSensor(redeclare Connectors.RealOutput y(
                final unit="N/T"));

        equation
          y = negative.Ndot "Measurement";
          negative.mu = positive.mu
            "Condition of equal electrochemical potentials";

          annotation (Icon(graphics={Text(
                          extent={{-100,-20},{100,-50}},
                          lineColor={127,127,127},
                          textString="N"),Text(
                          extent={{-100,8},{100,-22}},
                          lineColor={127,127,127},
                          textString=".")}));
        end Current;

        package BaseClasses "Base classes (not for direct use)"
          extends Modelica.Icons.BasesPackage;

          partial model PartialSensor
            "<html>Partial sensor for a pair of <a href=\"modelica://FCSys.Connectors.Material\">Material</a> connectors</html>"

            extends FCSys.Sensors.BaseClasses.PartialSensor;

            FCSys.Connectors.MaterialTransport negative
              "Material connector for the negative face" annotation (Placement(
                  transformation(extent={{-110,-10},{-90,10}}),
                  iconTransformation(extent={{-110,-10},{-90,10}})));
            FCSys.Connectors.MaterialTransport positive
              "Material connector for the positive face" annotation (Placement(
                  transformation(extent={{90,-10},{110,10}}),
                  iconTransformation(extent={{90,-10},{110,10}})));

          equation
            0 = negative.Ndot + positive.Ndot
              "Material rate balance (no storage)";

            annotation (Icon(graphics={Line(
                              points={{-70,0},{-100,0}},
                              color={0,0,0},
                              smooth=Smooth.None),Line(
                              points={{100,0},{70,0}},
                              color={0,0,0},
                              smooth=Smooth.None)}), Diagram(graphics));
          end PartialSensor;
        end BaseClasses;
      end Material;

      package Momentum "Relative sensors for linear momentum"
        extends Modelica.Icons.Package;
        model Velocity "Sensor for relative velocity"

          extends BaseClasses.PartialSensor(redeclare Connectors.RealOutput y(
                final unit="l/T"));

        equation
          y = negative.phi - positive.phi "Measurement";
          0 = negative.mPhidot "Condition of no compressive force";
          // Note:  In conjunction with the momentum rate balance, this means
          // that there's no force on either face.
          annotation (
            Documentation(info="<html><p>If <code>ax=1</code>, then the difference in normal velocity
  is measured.  The relative first
  and second transverse components of velocity may be measured by
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

        model Force "Sensor for compressive force"
          extends BaseClasses.PartialSensor(redeclare Connectors.RealOutput y(
                final unit="l2.m/T2"));

        equation
          y = negative.mPhidot "Measurement";
          negative.phi = positive.phi "Condition of equal velocities";
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
            "<html>Partial sensor for a pair of <a href=\"modelica://FCSys.Connectors.MomentumLineic\">MomentumLineic</a> connectors</html>"

            extends FCSys.Sensors.BaseClasses.PartialSensor;

            FCSys.Connectors.MechanicalTransport negative
              "Linear momentum connector for the negative face" annotation (
                Placement(transformation(extent={{-110,-10},{-90,10}}),
                  iconTransformation(extent={{-110,-10},{-90,10}})));
            FCSys.Connectors.MechanicalTransport positive
              "Linear momentum connector for the positive face" annotation (
                Placement(transformation(extent={{90,-10},{110,10}}),
                  iconTransformation(extent={{90,-10},{110,10}})));
          equation
            0 = negative.mPhidot + positive.mPhidot
              "Linear momentum rate balance (no storage)";

            annotation (Icon(graphics={Line(
                              points={{-70,0},{-100,0}},
                              color={0,0,0},
                              smooth=Smooth.None),Line(
                              points={{100,0},{70,0}},
                              color={0,0,0},
                              smooth=Smooth.None)}), Diagram(graphics));
          end PartialSensor;
        end BaseClasses;
      end Momentum;

      package Entropy "Relative sensors for entropy"
        extends Modelica.Icons.Package;
        model Temperature "Sensor for temperature difference"

          extends BaseClasses.PartialSensor(redeclare Connectors.RealOutput y(
                final unit="l2.m/(N.T2)", displayUnit="K"));

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
                          lineThickness=0.5)}), Diagram(graphics));
        end Temperature;

        model HeatRate "Sensor for heat flow rate"
          extends BaseClasses.PartialSensor(redeclare Connectors.RealOutput y(
                final unit="l2.m/T3"));

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

        model EntropyRate "Sensor for entropy flow rate"
          extends BaseClasses.PartialSensor(redeclare Connectors.RealOutput y(
                final unit="N/T"));

        equation
          y = (negative.Qdot - positive.Qdot)/2 "Measurement";
          negative.T = positive.T "Isothermal condition";
          annotation (Icon(graphics={Text(
                          extent={{-100,-20},{100,-50}},
                          lineColor={127,127,127},
                          textString="S/2"),Text(
                          extent={{-100,8},{100,-22}},
                          lineColor={127,127,127},
                          textString=".  "),Polygon(
                          points={{-32,-26},{-24,-44},{-40,-44},{-32,-26}},
                          lineColor={127,127,127},
                          smooth=Smooth.None,
                          lineThickness=0.5)}));
        end EntropyRate;

        package BaseClasses "Base classes (not for direct use)"
          extends Modelica.Icons.BasesPackage;

          partial model PartialSensor
            "<html>Partial sensor for a pair of <a href=\"modelica://FCSys.Connectors.Entropy\">Entropy</a> connectors</html>"

            extends FCSys.Sensors.BaseClasses.PartialSensor;

            FCSys.Connectors.Heat negative(T(min=Modelica.Constants.small))
              "Entropy connector for the negative face" annotation (Placement(
                  transformation(extent={{-110,-10},{-90,10}}),
                  iconTransformation(extent={{-110,-10},{-90,10}})));
            FCSys.Connectors.Heat positive(T(min=Modelica.Constants.small))
              "Entropy connector for the positive face" annotation (Placement(
                  transformation(extent={{90,-10},{110,10}}),
                  iconTransformation(extent={{90,-10},{110,10}})));
            // By increasing the minimum from zero to small (1e-60), symbolic
            // simplification becomes possible.
          equation
            0 = negative.Qdot + positive.Qdot
              "Energy rate balance (no storage)";

            annotation (Icon(graphics={Line(
                              points={{-70,0},{-100,0}},
                              color={0,0,0},
                              smooth=Smooth.None),Line(
                              points={{100,0},{70,0}},
                              color={0,0,0},
                              smooth=Smooth.None)}), Diagram(graphics));
          end PartialSensor;
        end BaseClasses;
      end Entropy;
    end Species;

    annotation (Documentation(info="<html><p>Since the connectors in
<a href=\"modelica://FCSys\">FCSys</a> are hierarchical
(see the <a href=\"modelica://FCSys.Connectors\">Connectors</a> package),
the models for the sensors must be as well.  A
<a href=\"modelica://FCSys.Connectors.BaseClasses.PartialFace\">Face</a>
connector (<a href=\"modelica://FCSys.Connectors.FaceX\">FaceX</a>,
<a href=\"modelica://FCSys.Connectors.FaceY\">FaceY</a>,
<a href=\"modelica://FCSys.Connectors.FaceZ\">FaceZ</a>,
or <a href=\"modelica://FCSys.Connectors.FaceGeneric\">FaceGeneric</a>)
is used in <a href=\"modelica://FCSys.Subregions.Species\">Species</a> models,
and there is a corresponding <a href=\"modelica://FCSys.Sensors.Face.Species.Species\">Species
sensor</a> model in this package. The
<a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a>
connector is used in <a href=\"modelica://FCSys.Subregions.Phases\">Phase</a> models,
and there are corresponding <a href=\"modelica://FCSys.Sensors.Face.Phases\">Phase
sensor</a> models.  The
<a href=\"modelica://FCSys.Connectors.FaceBus\">FaceBus</a> connector is nested once more
in models such as the <a href=\"modelica://FCSys.Subregions.Subregion\">Subregion</a>,
and there is a corresponding <a href=\"modelica://FCSys.Sensors.Face.Subregion\">Subregion
sensor</a> model.
</p></html>"));
  end FaceDifferential;

  package BaseClasses "Base classes (not for direct use)"
    extends Modelica.Icons.BasesPackage;

    partial model PartialSensorNonideal
      "Partial sensor with dynamics, saturation, and noise"

      extends PartialSensor;

      replaceable FCSys.Sensors.BaseClasses.PartialSensor sensor constrainedby
        FCSys.Sensors.BaseClasses.PartialSensor "Sensor" annotation (
          choicesAllMatching=true, Placement(transformation(extent={{-10,44},{
                10,64}})));

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
Copyright 2007&ndash;2012, Georgia Tech Research Corporation.
</p>
<p>
<i>This Modelica package is <u>free</u> software and the use is completely at <u>your own risk</u>;
it can be redistributed and/or modified under the terms of the Modelica License 2. For license conditions (including the
disclaimer of warranty) see <a href=\"modelica://FCSys.UsersGuide.ModelicaLicense2\">
FCSys.UsersGuide.ModelicaLicense2</a> or visit <a href=\"http://www.modelica.org/licenses/ModelicaLicense2\">
http://www.modelica.org/licenses/ModelicaLicense2</a>.</i>
</p></html>"));
end Sensors;

within FCSys;
package Phases "Mixtures of species"
  extends Modelica.Icons.Package;
  // **Use __Dymola_label for incl and species tags
  model Gas "Gas phase"
    import Modelica.Math.BooleanVectors.countTrue;
    extends Partial(final n_spec=countTrue({inclH2,inclH2O,inclN2,inclO2}));

    // Conditionally include species.
    parameter Boolean inclH2=false "Include H2" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>Hydrogen (H<sub>2</sub>)</html>",
        __Dymola_joinNext=true));

    replaceable FCSys.Species.H2.Gas.Fixed H2(final n_faces) if inclH2
      constrainedby FCSys.Species.Compressible(
      n_faces=n_faces,
      phi(each stateSelect=if reduceTrans then StateSelect.default else
            StateSelect.prefer),
      T(stateSelect=if reduceThermal then StateSelect.default else StateSelect.prefer))
      "H2 model" annotation (
      __Dymola_choicesFromPackage=true,
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>H<sub>2</sub> model</html>",
        enable=inclH2),
      Placement(transformation(extent={{-80,-10},{-60,10}})));

    parameter Boolean inclH2O=false "Include H2O" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>Water (H<sub>2</sub>O)</html>",
        __Dymola_joinNext=true));

    replaceable FCSys.Species.H2O.Gas.Fixed H2O(final n_faces) if inclH2O
      constrainedby FCSys.Species.Compressible(
      n_faces=n_faces,
      phi(each stateSelect=if reduceTrans then StateSelect.default else
            StateSelect.prefer),
      T(stateSelect=if reduceThermal then StateSelect.default else StateSelect.prefer))
      "H2O model" annotation (
      __Dymola_choicesFromPackage=true,
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>H<sub>2</sub>O model</html>",
        enable=inclH2O),
      Placement(transformation(extent={{-40,-10},{-20,10}})));

    parameter Boolean inclN2=false "Include N2" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>Nitrogen (N<sub>2</sub>)</html>",
        __Dymola_joinNext=true));

    replaceable FCSys.Species.N2.Gas.Fixed N2(final n_faces) if inclN2
      constrainedby FCSys.Species.Compressible(
      n_faces=n_faces,
      phi(each stateSelect=if reduceTrans then StateSelect.default else
            StateSelect.prefer),
      T(stateSelect=if reduceThermal then StateSelect.default else StateSelect.prefer))
      "N2 model" annotation (
      __Dymola_choicesFromPackage=true,
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>N<sub>2</sub> model</html>",
        enable=inclN2),
      Placement(transformation(extent={{0,-10},{20,10}})));

    parameter Boolean inclO2=false "Include O2" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>Oxygen (O<sub>2</sub>)</html>",
        __Dymola_joinNext=true));

    replaceable FCSys.Species.O2.Gas.Fixed O2(final n_faces) if inclO2
      constrainedby FCSys.Species.Compressible(
      n_faces=n_faces,
      phi(each stateSelect=if reduceTrans then StateSelect.default else
            StateSelect.prefer),
      T(stateSelect=if reduceThermal then StateSelect.default else StateSelect.prefer))
      "O2 model" annotation (
      __Dymola_choicesFromPackage=true,
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>O<sub>2</sub> model</html>",
        enable=inclO2),
      Placement(transformation(extent={{40,-10},{60,10}})));

    parameter Boolean inclHOR=false "Hydrogen oxidation" annotation (
      HideResult=true,
      Dialog(group="Included reactions", compact=true),
      choices(__Dymola_checkBox=true));
    parameter Boolean inclORR=false "Oxygen reduction" annotation (
      HideResult=true,
      Dialog(group="Included reactions", compact=true),
      choices(__Dymola_checkBox=true));
    // These can't be outer parameters in Dymola 7.4.

    // Auxiliary variables (for analysis)
    output Q.PressureAbsolute p(stateSelect=StateSelect.never) = amagat.p if
      n_spec > 0 "Total thermodynamic pressure";

    Connectors.FaceBus xNegative if inclFaces[Axis.x]
      "Negative face along the x axis" annotation (Placement(transformation(
            extent={{-130,-10},{-110,10}}), iconTransformation(extent={{-90,-10},
              {-70,10}})));
    Connectors.FaceBus yNegative if inclFaces[Axis.y]
      "Negative face along the y axis" annotation (Placement(transformation(
            extent={{-106,-34},{-86,-14}}), iconTransformation(extent={{-10,-94},
              {10,-74}})));
    Connectors.FaceBus zNegative if inclFaces[Axis.z]
      "Negative face along the z axis" annotation (Placement(transformation(
            extent={{78,2},{98,22}}), iconTransformation(extent={{40,40},{60,60}})));
    Connectors.FaceBus xPositive if inclFaces[Axis.x]
      "Positive face along the x axis" annotation (Placement(transformation(
            extent={{90,-10},{110,10}}),iconTransformation(extent={{70,-10},{90,
              10}})));
    Connectors.FaceBus yPositive if inclFaces[Axis.y]
      "Positive face along the y axis" annotation (Placement(transformation(
            extent={{66,14},{86,34}}), iconTransformation(extent={{-10,90},{10,
              110}})));
    Connectors.FaceBus zPositive if inclFaces[Axis.z]
      "Positive face along the z axis" annotation (Placement(transformation(
            extent={{-118,-22},{-98,-2}}), iconTransformation(extent={{-90,-90},
              {-70,-70}})));
    Connectors.Amagat amagat if n_spec > 0 "Connector for additivity of volume"
      annotation (Placement(transformation(extent={{90,-64},{110,-44}}),
          iconTransformation(extent={{-60,40},{-40,60}})));
    Connectors.Inter inter[n_inter](each final n_trans=n_trans) if n_spec > 0
      "Connector to exchange momentum and energy with other phases" annotation
      (Placement(transformation(extent={{90,-40},{110,-20}}),
          iconTransformation(extent={{80,-60},{100,-80}})));
    Connectors.Stoichiometric chemical(final n_trans=n_trans) if inclHOR or
      inclORR "Connector for an electrochemical reaction" annotation (Placement(
          transformation(extent={{-70,50},{-50,70}}), iconTransformation(extent
            ={{50,-94},{70,-74}})));
    Connectors.PhysicalBus physical if inclH2O annotation (Placement(
          transformation(extent={{-49,50},{-29,70}}), iconTransformation(extent
            ={{20,-94},{40,-74}})));

    // Aliases
    /*
    Q.Velocity phi[n_trans](each stateSelect=StateSelect.prefer) = direct.translational.phi
      if n_spec > 1 and reduceTrans "Velocity";
    Q.Temperature T(each stateSelect=StateSelect.prefer) = direct.thermal.T if 
      n_spec > 1 and reduceThermal "Temperature";
    // These make the selected states more readable.
  */

  protected
    Conditions.Adapters.AmagatDalton DA if n_spec > 0
      "Adapter between additivity of pressure and additivity of volume"
      annotation (Placement(transformation(extent={{90,-64},{70,-44}})));
    Conditions.Adapters.ChemicalReaction HOR(
      final n_trans=n_trans,
      n=-1,
      final m=H2.Data.m) if inclHOR
      "Adapter for the contribution of H2 to the hydrogen oxidation reaction"
      annotation (Placement(transformation(extent={{-81,30},{-61,50}})));
    Conditions.Adapters.ChemicalReaction ORR[2](
      n={-1,2},
      final m={O2.Data.m,H2O.Data.m},
      each final n_trans=n_trans) if inclORR
      "Adapter for the contribution of O2 and H2O to the oxygen reduction reaction"
      annotation (Placement(transformation(extent={{-39,30},{-59,50}})));
    Connectors.DirectNode direct(
      final n_trans=n_trans,
      final inclTrans=reduceTrans,
      final inclThermal=reduceThermal) if n_spec > 0 and (reduceTrans or
      reduceThermal)
      "Connector to directly couple velocities and temperatures within the phase"
      annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={
              -4,-42}), iconTransformation(extent={{-10,-10},{10,10}}, origin={
              26,-52})));

  equation
    // Phase change
    connect(H2O.physical, physical.H2O) annotation (Line(
        points={{-38.6,5},{-39,5},{-39,60}},
        color={38,196,52},
        smooth=Smooth.None));

    // Reactions
    // ---------
    // HOR
    connect(H2.chemical, HOR.chemical) annotation (Line(
        points={{-75,8.8},{-75,40}},
        color={255,195,38},
        smooth=Smooth.None));
    connect(HOR.reaction, chemical) annotation (Line(
        points={{-67,40},{-60,40},{-60,60}},
        color={255,195,38},
        smooth=Smooth.None));

    // ORR
    connect(O2.chemical, ORR[1].chemical) annotation (Line(
        points={{45,8.8},{45,40},{-45,40}},
        color={255,195,38},
        smooth=Smooth.None));
    connect(H2O.chemical, ORR[2].chemical) annotation (Line(
        points={{-35,8.8},{-35,40},{-45,40}},
        color={255,195,38},
        smooth=Smooth.None));
    connect(ORR[1].reaction, chemical) annotation (Line(
        points={{-53,40},{-60,40},{-60,60}},
        color={255,195,38},
        smooth=Smooth.None));
    connect(ORR[2].reaction, chemical) annotation (Line(
        points={{-53,40},{-60,40},{-60,60}},
        color={255,195,38},
        smooth=Smooth.None));

    // Inert exchange
    connect(H2.inter, inter) annotation (Line(
        points={{-60.6,-3},{-58,-6},{-58,-30},{100,-30}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(H2O.inter, inter) annotation (Line(
        points={{-20.6,-3},{-18,-6},{-18,-30},{100,-30}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(O2.inter, inter) annotation (Line(
        points={{59.4,-3},{62,-6},{62,-30},{100,-30}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(N2.inter, inter) annotation (Line(
        points={{19.4,-3},{22,-6},{22,-30},{100,-30}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(H2.direct.translational, direct.translational) annotation (Line(
        points={{-64,-8},{-64,-42},{-4,-42}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(H2.direct.thermal, direct.thermal) annotation (Line(
        points={{-64,-8},{-64,-42},{-4,-42}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(H2O.direct.translational, direct.translational) annotation (Line(
        points={{-24,-8},{-24,-42},{-4,-42}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(H2O.direct.thermal, direct.thermal) annotation (Line(
        points={{-24,-8},{-24,-42},{-4,-42}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(N2.direct.translational, direct.translational) annotation (Line(
        points={{16,-8},{16,-42},{-4,-42}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(N2.direct.thermal, direct.thermal) annotation (Line(
        points={{16,-8},{16,-42},{-4,-42}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(O2.direct.translational, direct.translational) annotation (Line(
        points={{56,-8},{56,-42},{-4,-42}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(O2.direct.thermal, direct.thermal) annotation (Line(
        points={{56,-8},{56,-42},{-4,-42}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(H2.dalton, DA.dalton) annotation (Line(
        points={{-67,-9.4},{-67,-54},{76,-54}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(H2O.dalton, DA.dalton) annotation (Line(
        points={{-27,-9.4},{-27,-54},{76,-54}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(N2.dalton, DA.dalton) annotation (Line(
        points={{13,-9.4},{13,-54},{76,-54}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(O2.dalton, DA.dalton) annotation (Line(
        points={{53,-9.4},{53,-54},{76,-54}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(DA.amagat, amagat) annotation (Line(
        points={{84,-54},{100,-54}},
        color={47,107,251},
        smooth=Smooth.None));

    // Transport
    // ---------
    // H2
    connect(H2.faces[facesCart[Axis.x], Side.n], xNegative.H2) annotation (Line(
        points={{-70,6.10623e-16},{-120,5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2.faces[facesCart[Axis.x], Side.p], xPositive.H2) annotation (Line(
        points={{-70,6.10623e-16},{100,5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2.faces[facesCart[Axis.y], Side.n], yNegative.H2) annotation (Line(
        points={{-70,6.10623e-16},{-70,-24},{-96,-24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2.faces[facesCart[Axis.y], Side.p], yPositive.H2) annotation (Line(
        points={{-70,6.10623e-16},{-70,24},{76,24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2.faces[facesCart[Axis.z], Side.n], zNegative.H2) annotation (Line(
        points={{-70,6.10623e-16},{-58,12},{88,12}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(H2.faces[facesCart[Axis.z], Side.p], zPositive.H2) annotation (Line(
        points={{-70,6.10623e-16},{-82,-12},{-108,-12}},
        color={127,127,127},
        smooth=Smooth.None));
    // H2O
    connect(H2O.faces[facesCart[Axis.x], Side.n], xNegative.H2O) annotation (
        Line(
        points={{-30,6.10623e-16},{-120,5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.faces[facesCart[Axis.x], Side.p], xPositive.H2O) annotation (
        Line(
        points={{-30,6.10623e-16},{100,5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.faces[facesCart[Axis.y], Side.n], yNegative.H2O) annotation (
        Line(
        points={{-30,6.10623e-16},{-30,-24},{-96,-24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.faces[facesCart[Axis.y], Side.p], yPositive.H2O) annotation (
        Line(
        points={{-30,6.10623e-16},{-30,24},{76,24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.faces[facesCart[Axis.z], Side.n], zNegative.H2O) annotation (
        Line(
        points={{-30,6.10623e-16},{-18,12},{88,12}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(H2O.faces[facesCart[Axis.z], Side.p], zPositive.H2O) annotation (
        Line(
        points={{-30,6.10623e-16},{-42,-12},{-108,-12}},
        color={127,127,127},
        smooth=Smooth.None));
    // N2
    connect(N2.faces[facesCart[Axis.x], Side.n], xNegative.N2) annotation (Line(
        points={{10,6.10623e-16},{-120,5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(N2.faces[facesCart[Axis.x], Side.p], xPositive.N2) annotation (Line(
        points={{10,6.10623e-16},{100,5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(N2.faces[facesCart[Axis.y], Side.n], yNegative.N2) annotation (Line(
        points={{10,6.10623e-16},{10,-24},{-96,-24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(N2.faces[facesCart[Axis.y], Side.p], yPositive.N2) annotation (Line(
        points={{10,6.10623e-16},{10,24},{76,24}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(N2.faces[facesCart[Axis.z], Side.n], zNegative.N2) annotation (Line(
        points={{10,6.10623e-16},{10,0},{22,12},{88,12}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(N2.faces[facesCart[Axis.z], Side.p], zPositive.N2) annotation (Line(
        points={{10,6.10623e-16},{-2,-12},{-108,-12}},
        color={127,127,127},
        smooth=Smooth.None));

    // O2
    connect(O2.faces[facesCart[Axis.x], Side.n], xNegative.O2) annotation (Line(
        points={{50,6.10623e-16},{-120,5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(O2.faces[facesCart[Axis.x], Side.p], xPositive.O2) annotation (Line(
        points={{50,6.10623e-16},{100,5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(O2.faces[facesCart[Axis.y], Side.n], yNegative.O2) annotation (Line(
        points={{50,6.10623e-16},{50,-24},{-96,-24}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(O2.faces[facesCart[Axis.y], Side.p], yPositive.O2) annotation (Line(
        points={{50,6.10623e-16},{50,24},{76,24}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(O2.faces[facesCart[Axis.z], Side.n], zNegative.O2) annotation (Line(
        points={{50,6.10623e-16},{62,12},{88,12}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(O2.faces[facesCart[Axis.z], Side.p], zPositive.O2) annotation (Line(
        points={{50,6.10623e-16},{38,-12},{-108,-12}},
        color={127,127,127},
        smooth=Smooth.None));
    annotation (
      Documentation(info="<html>
<p>Please see the documentation of the <a href=\"modelica://FCSys.Phases.BaseClasses.EmptyPhase\">EmptyPhase</a> model.</p></html>"),

      Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-120,-100},{
              100,100}}), graphics),
      Icon(coordinateSystem(preserveAspectRatio=true, extent={{-120,-100},{100,
              100}}), graphics));
  end Gas;

  model Graphite "Graphite phase"
    import assert = FCSys.Utilities.assertEval;
    import Modelica.Math.BooleanVectors.countTrue;
    extends Partial(final n_spec=countTrue({'inclC+','incle-'}));

    // Conditionally include species.
    parameter Boolean 'inclC+'=false "Include C+" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>Carbon plus (C<sup>+</sup>)</html>",
        __Dymola_joinNext=true));

    replaceable FCSys.Species.'C+'.Graphite.Fixed 'C+'(final n_faces) if
      'inclC+' constrainedby FCSys.Species.Isochoric(
      n_faces=n_faces,
      phi(each stateSelect=if reduceTrans then StateSelect.default else
            StateSelect.prefer),
      T(stateSelect=if reduceThermal then StateSelect.default else StateSelect.prefer))
      "C+ model" annotation (
      __Dymola_choicesFromPackage=true,
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>C<sup>+</sup> model</html>",
        enable='inclC+'),
      Placement(transformation(extent={{-30,-10},{-10,10}})));

    parameter Boolean 'incle-'=false "Include e-" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>Electrons (e<sup>-</sup>)</html>",
        __Dymola_joinNext=true));

    replaceable FCSys.Species.'e-'.Graphite.Fixed 'e-'(final n_faces) if
      'incle-' and not (inclHOR or inclORR) constrainedby
      FCSys.Species.Compressible(
      n_faces=n_faces,
      phi(each stateSelect=if reduceTrans then StateSelect.default else
            StateSelect.prefer),
      T(stateSelect=if reduceThermal then StateSelect.default else StateSelect.prefer))
      "e- model" annotation (
      __Dymola_choicesFromPackage=true,
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>'e-' model</html>",
        enable='incle-'),
      Placement(transformation(extent={{10,-10},{30,10}})));

    //initMaterial=if inclVoid or 'inclC+' then else ,

    // Reaction parameters
    // -------------------
    parameter Boolean inclHOR=false "Include hydrogen oxidation" annotation (
      HideResult=true,
      Dialog(group="Reactions", compact=true),
      choices(__Dymola_checkBox=true));
    parameter Boolean inclORR=false "Include oxygen reduction" annotation (
      HideResult=true,
      Dialog(group="Reactions", compact=true),
      choices(__Dymola_checkBox=true));
    parameter Q.CurrentAreicAbsolute J0=1e-7*U.A/U.cm^2 if inclHOR or inclORR
      "Exchange current density" annotation (Dialog(
        group="Reactions",
        enable=inclHOR or inclORR,
        __Dymola_label="<html><i>J</i><sup>o</sup></html>"));
    parameter Q.CurrentAreic J_irr=0 if inclHOR or inclORR
      "Irreversible current of side-reactions" annotation (Dialog(
        group="Reactions",
        enable=inclHOR or inclORR,
        __Dymola_label="<html><i>J</i><sub>irr</sub></html>"));
    // These can't be outer parameters in Dymola 7.4.

    Connectors.FaceBus xNegative if inclFaces[Axis.x]
      "Negative face along the x axis" annotation (Placement(transformation(
            extent={{-80,-10},{-60,10}}), iconTransformation(extent={{-90,-10},
              {-70,10}})));
    Connectors.FaceBus yNegative if inclFaces[Axis.y]
      "Negative face along the y axis" annotation (Placement(transformation(
            extent={{-56,-34},{-36,-14}}),iconTransformation(extent={{-10,-94},
              {10,-74}})));
    Connectors.FaceBus zNegative if inclFaces[Axis.z]
      "Negative face along the z axis" annotation (Placement(transformation(
            extent={{48,2},{68,22}}), iconTransformation(extent={{40,40},{60,60}})));
    Connectors.FaceBus xPositive if inclFaces[Axis.x]
      "Positive face along the x axis" annotation (Placement(transformation(
            extent={{60,-10},{80,10}}), iconTransformation(extent={{70,-10},{90,
              10}})));
    Connectors.FaceBus yPositive if inclFaces[Axis.y]
      "Positive face along the y axis" annotation (Placement(transformation(
            extent={{36,14},{56,34}}), iconTransformation(extent={{-10,90},{10,
              110}})));
    Connectors.FaceBus zPositive if inclFaces[Axis.z]
      "Positive face along the z axis" annotation (Placement(transformation(
            extent={{-68,-22},{-48,-2}}), iconTransformation(extent={{-90,-90},
              {-70,-70}})));
    Connectors.Inter inter[n_inter](each final n_trans=n_trans) if n_spec > 0
      "Connector to exchange momentum and energy with other phases" annotation
      (Placement(transformation(extent={{60,-40},{80,-20}}), iconTransformation(
            extent={{80,-60},{100,-80}})));

    Connectors.Amagat amagat if n_spec > 0 "Connector for additivity of volume"
      annotation (Placement(transformation(extent={{60,-64},{80,-44}}),
          iconTransformation(extent={{-60,40},{-40,60}})));
    Connectors.Stoichiometric chemical(final n_trans=n_trans) if inclHOR or
      inclORR "Connector for an electrochemical reaction" annotation (Placement(
          transformation(extent={{-10,70},{10,90}}), iconTransformation(extent=
              {{50,-94},{70,-74}})));

    Conditions.Adapters.FaceReaction HOR(
      final axis=Axis.x,
      final side=Side.n,
      final A=A[Axis.x],
      final cartTrans=cartTrans,
      redeclare package Data = FCSys.Characteristics.'e-'.Gas,
      n=2) if inclHOR
      "Adapter for the contribution of e- to the hydrogen oxidation reaction"
      annotation (Placement(transformation(extent={{-70,50},{-50,70}})));

    Conditions.Adapters.FaceReaction ORR(
      final axis=Axis.x,
      final side=Side.p,
      final A=A[Axis.x],
      cartTrans=cartTrans,
      redeclare package Data = FCSys.Characteristics.'e-'.Gas,
      n=-4) if inclORR
      "Adapter for the contribution of e- to the oxygen reduction reaction"
      annotation (Placement(transformation(extent={{70,50},{50,70}})));

    // Aliases
    /*
    Q.Velocity phi[n_trans](each stateSelect=StateSelect.prefer) = direct.translational.phi
      if n_spec > 1 and reduceTrans "Velocity";
    Q.Temperature T(each stateSelect=StateSelect.prefer) = direct.thermal.T if 
      n_spec > 1 and reduceThermal "Temperature";
    // These make the selected states more readable.
*/

  protected
    FCSys.Species.Reaction reaction(
      final A=A[Axis.x],
      transSubstrate=true,
      thermalSubstrate=true,
      fromCurrent=true,
      J0=J0,
      J_irr=J_irr) if inclHOR or inclORR "Model to determine the reaction rate"
      annotation (Placement(transformation(extent={{-24,50},{-4,70}})));

    Connectors.DirectNode direct(
      final n_trans=n_trans,
      final inclTrans=reduceTrans,
      final inclThermal=reduceThermal) if n_spec > 0 and (reduceTrans or
      reduceThermal)
      "Connector to directly couple velocities and temperatures within the phase"
      annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={
              6,-42}), iconTransformation(extent={{-10,-10},{10,10}}, origin={
              26,-52})));

  protected
    Conditions.Adapters.AmagatDalton DA if n_spec > 0
      "Adapter between additivity of pressure and additivity of volume"
      annotation (Placement(transformation(extent={{62,-64},{42,-44}})));
  equation
    // Reactions
    // ---------
    connect(reaction.inert, 'C+'.direct) annotation (Line(
        points={{-14,56},{-14,-8}},
        color={47,107,251},
        smooth=Smooth.None));
    // HOR
    connect(HOR.face, xNegative.'e-') annotation (Line(
        points={{-64,60},{-70,60},{-70,5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None));
    // ORR
    connect(ORR.face, xPositive.'e-') annotation (Line(
        points={{64,60},{70,60},{70,5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None));

    // Inert exchange
    connect('C+'.inter, inter) annotation (Line(
        points={{-10.6,-3},{-8,-6},{-8,-30},{70,-30}},
        color={47,107,251},
        smooth=Smooth.None));
    connect('e-'.inter, inter) annotation (Line(
        points={{29.4,-3},{32,-6},{32,-30},{70,-30}},
        color={47,107,251},
        smooth=Smooth.None));
    connect('C+'.direct.translational, direct.translational) annotation (Line(
        points={{-14,-8},{-14,-42},{6,-42}},
        color={47,107,251},
        smooth=Smooth.None));
    connect('C+'.direct.thermal, direct.thermal) annotation (Line(
        points={{-14,-8},{-14,-42},{6,-42}},
        color={47,107,251},
        smooth=Smooth.None));
    connect('e-'.direct.translational, direct.translational) annotation (Line(
        points={{26,-8},{26,-42},{6,-42}},
        color={47,107,251},
        smooth=Smooth.None));
    connect('e-'.direct.thermal, direct.thermal) annotation (Line(
        points={{26,-8},{26,-42},{6,-42}},
        color={47,107,251},
        smooth=Smooth.None));

    // Transport
    // ---------
    // C+
    connect('C+'.faces[facesCart[Axis.x], Side.n], xNegative.'C+') annotation (
        Line(
        points={{-20,6.10623e-16},{-20,5.55112e-16},{-70,5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('C+'.faces[facesCart[Axis.x], Side.p], xPositive.'C+') annotation (
        Line(
        points={{-20,6.10623e-16},{20,0},{20,5.55112e-16},{70,5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('C+'.faces[facesCart[Axis.y], Side.n], yNegative.'C+') annotation (
        Line(
        points={{-20,6.10623e-16},{-20,-24},{-46,-24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('C+'.faces[facesCart[Axis.y], Side.p], yPositive.'C+') annotation (
        Line(
        points={{-20,6.10623e-16},{-20,24},{46,24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('C+'.faces[facesCart[Axis.z], Side.n], zNegative.'C+') annotation (
        Line(
        points={{-20,6.10623e-16},{-8,12},{58,12}},
        color={127,127,127},
        smooth=Smooth.None));
    connect('C+'.faces[facesCart[Axis.z], Side.p], zPositive.'C+') annotation (
        Line(
        points={{-20,6.10623e-16},{-32,-12},{-58,-12}},
        color={127,127,127},
        smooth=Smooth.None));

    // e-
    connect('e-'.faces[facesCart[Axis.x], Side.n], xNegative.'e-') annotation (
        Line(
        points={{20,6.10623e-16},{-20,0},{-20,5.55112e-16},{-70,5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('e-'.faces[facesCart[Axis.x], Side.p], xPositive.'e-') annotation (
        Line(
        points={{20,6.10623e-16},{20,5.55112e-16},{70,5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('e-'.faces[facesCart[Axis.y], Side.n], yNegative.'e-') annotation (
        Line(
        points={{20,6.10623e-16},{20,-24},{-46,-24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('e-'.faces[facesCart[Axis.y], Side.p], yPositive.'e-') annotation (
        Line(
        points={{20,6.10623e-16},{20,24},{46,24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('e-'.faces[facesCart[Axis.z], Side.n], zNegative.'e-') annotation (
        Line(
        points={{20,6.10623e-16},{32,12},{58,12}},
        color={127,127,127},
        smooth=Smooth.None));
    connect('e-'.faces[facesCart[Axis.z], Side.p], zPositive.'e-') annotation (
        Line(
        points={{20,6.10623e-16},{8,-12},{-58,-12}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(HOR.reaction, reaction.reaction) annotation (Line(
        points={{-56,60},{-14,60}},
        color={255,195,38},
        smooth=Smooth.None));

    connect(ORR.reaction, reaction.reaction) annotation (Line(
        points={{56,60},{-14,60}},
        color={255,195,38},
        smooth=Smooth.None));
    connect(chemical, reaction.reaction) annotation (Line(
        points={{5.55112e-16,80},{5.55112e-16,60},{-14,60}},
        color={255,195,38},
        smooth=Smooth.None));
    connect(DA.amagat, amagat) annotation (Line(
        points={{56,-54},{70,-54}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(DA.dalton, 'C+'.dalton) annotation (Line(
        points={{48,-54},{-17,-54},{-17,-9.4}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(DA.dalton, 'e-'.dalton) annotation (Line(
        points={{48,-54},{23,-54},{23,-9.4}},
        color={47,107,251},
        smooth=Smooth.None));
    annotation (
      Documentation(info="<html>
    <p>See <a href=\"modelica://FCSys.Species.'e-'.Graphite.Fixed\">Species.'e-'.Graphite.Fixed</a> for assumptions.
    For more information, see the
 <a href=\"modelica://FCSys.Phases.BaseClasses.EmptyPhase\">EmptyPhase</a> model.</p></html>"),

      Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
              100,100}}), graphics),
      Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
              100}}), graphics));
  end Graphite;

  model Ionomer "Ionomer phase"
    import Modelica.Math.BooleanVectors.countTrue;
    extends Partial(final n_spec=countTrue({'inclSO3-','inclH+',inclH2O}));

    parameter Q.NumberAbsolute k_EOD=1 if 'inclH+' or inclH2O
      "Coupling factor between H+ and H2O" annotation (Dialog(group="Geometry",
          __Dymola_label="<html><i>k</i><sub>EOD</sub></html>"));
    parameter Q.NumberAbsolute k_PEMH2O=18 if 'inclSO3-' or inclH2O
      "Coupling factor between SO3- and H2O" annotation (Dialog(group=
            "Geometry", __Dymola_label=
            "<html><i>k</i><sub>PEM H2O</sub></html>"));

    // Conditionally include species.
    parameter Boolean 'inclSO3-'=false
      "Include C19HF37O5S- (abbreviated as SO3-)" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label=
            "<html>Nafion sulfonate minus (abbreviated as SO<sub>3</sub><sup>-</sup>)</html>",

        __Dymola_joinNext=true));

    replaceable FCSys.Species.'SO3-'.Ionomer.Fixed 'SO3-'(final n_faces) if
      'inclSO3-' constrainedby FCSys.Species.Isochoric(
      n_faces=n_faces,
      V=V,
      n_intra=1,
      k_intra={k_PEMH2O},
      phi(each stateSelect=if reduceTrans then StateSelect.default else
            StateSelect.prefer),
      T(stateSelect=if reduceThermal then StateSelect.default else StateSelect.prefer))
      "SO3- model" annotation (
      __Dymola_choicesFromPackage=true,
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>SO<sub>3</sub><sup>-</sup> model</html>",
        enable='inclSO3-'),
      Placement(transformation(extent={{30,-10},{50,10}})));

    parameter Boolean 'inclH+'=false "Include H+" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>Protons (H<sup>+</sup>)</html>",
        __Dymola_joinNext=true));

    replaceable FCSys.Species.'H+'.Ionomer.Fixed 'H+'(final n_faces) if
      'inclH+' and not (inclHOR or inclORR) constrainedby
      FCSys.Species.Compressible(
      n_faces=n_faces,
      V=V,
      n_intra=1,
      k_intra={k_EOD},
      phi(each stateSelect=if reduceTrans then StateSelect.default else
            StateSelect.prefer),
      T(stateSelect=if reduceThermal then StateSelect.default else StateSelect.prefer))
      "H+ model" annotation (
      __Dymola_choicesFromPackage=true,
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>H<sup>+</sup> model</html>",
        enable='inclH+'),
      Placement(transformation(extent={{-50,-10},{-30,10}})));

    parameter Boolean inclH2O=false "Include H2O" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>Water (H<sub>2</sub>O)</html>",
        __Dymola_joinNext=true));

    replaceable FCSys.Species.H2O.Ionomer.Fixed H2O(final n_faces) if inclH2O
      constrainedby FCSys.Species.Compressible(
      n_faces=n_faces,
      V=V,
      n_intra=2,
      k_intra={k_EOD,k_PEMH2O},
      phi(each stateSelect=if reduceTrans then StateSelect.default else
            StateSelect.prefer),
      T(stateSelect=if reduceThermal then StateSelect.default else StateSelect.prefer))
      "H2O model" annotation (
      __Dymola_choicesFromPackage=true,
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>H<sub>2</sub>O model</html>",
        enable=inclH2O),
      Placement(transformation(extent={{-10,-10},{10,10}})));

    parameter Boolean inclHOR=false "Hydrogen oxidation" annotation (
      HideResult=true,
      Dialog(group="Included reactions", compact=true),
      choices(__Dymola_checkBox=true));
    parameter Boolean inclORR=false "Oxygen reduction" annotation (
      HideResult=true,
      Dialog(group="Included reactions", compact=true),
      choices(__Dymola_checkBox=true));
    // These can't be outer parameters in Dymola 7.4.

    Connectors.FaceBus xNegative if inclFaces[Axis.x]
      "Negative face along the x axis" annotation (Placement(transformation(
            extent={{-100,-10},{-80,10}}),iconTransformation(extent={{-90,-10},
              {-70,10}})));
    Connectors.FaceBus yNegative if inclFaces[Axis.y]
      "Negative face along the y axis" annotation (Placement(transformation(
            extent={{-76,-34},{-56,-14}}),iconTransformation(extent={{-10,-94},
              {10,-74}})));
    Connectors.FaceBus zNegative if inclFaces[Axis.z]
      "Negative face along the z axis" annotation (Placement(transformation(
            extent={{68,2},{88,22}}), iconTransformation(extent={{40,40},{60,60}})));
    Connectors.FaceBus xPositive if inclFaces[Axis.x]
      "Positive face along the x axis" annotation (Placement(transformation(
            extent={{80,-10},{100,10}}),iconTransformation(extent={{70,-10},{90,
              10}})));
    Connectors.FaceBus yPositive if inclFaces[Axis.y]
      "Positive face along the y axis" annotation (Placement(transformation(
            extent={{56,14},{76,34}}), iconTransformation(extent={{-10,90},{10,
              110}})));
    Connectors.FaceBus zPositive if inclFaces[Axis.z]
      "Positive face along the z axis" annotation (Placement(transformation(
            extent={{-88,-22},{-68,-2}}), iconTransformation(extent={{-90,-90},
              {-70,-70}})));
    Connectors.Inter inter[n_inter](each final n_trans=n_trans) if n_spec > 0
      "Connector to exchange momentum and energy with other phases" annotation
      (Placement(transformation(extent={{80,-40},{100,-20}}),
          iconTransformation(extent={{80,-60},{100,-80}})));
    Connectors.Amagat amagat if n_spec > 0 "Connector for additivity of volume"
      annotation (Placement(transformation(extent={{80,-76},{100,-56}}),
          iconTransformation(extent={{-60,40},{-40,60}})));

    Connectors.Stoichiometric chemical(final n_trans=n_trans) if inclHOR or
      inclORR "Connector for an electrochemical reaction" annotation (Placement(
          transformation(extent={{10,70},{30,90}}), iconTransformation(extent={
              {50,-94},{70,-74}})));
    Connectors.PhysicalBus physical if inclH2O annotation (Placement(
          transformation(extent={{-19,70},{1,90}}), iconTransformation(extent={
              {20,-94},{40,-74}})));
    Conditions.Adapters.FaceReaction HOR(
      final axis=Axis.x,
      final side=Side.p,
      final A=A[Axis.x],
      redeclare package Data = FCSys.Characteristics.'H+'.Gas,
      cartTrans=cartTrans,
      n=2) if inclHOR
      "Charge layer or adapter for the contribution of H+ to the hydrogen oxidation reaction"
      annotation (Placement(transformation(extent={{90,50},{70,70}})));

    Conditions.Adapters.FaceReaction ORR(
      final axis=Axis.x,
      final side=Side.n,
      final A=A[Axis.x],
      redeclare package Data = FCSys.Characteristics.'H+'.Gas,
      cartTrans=cartTrans,
      n=-4) if inclORR
      "Charge layer or adapter for the contribution of H+ to the oxygen reduction reaction"
      annotation (Placement(transformation(extent={{-90,50},{-70,70}})));

    // Aliases
    /*
    Q.Velocity phi[n_trans](each stateSelect=StateSelect.prefer) = direct.translational.phi
      if n_spec > 1 and reduceTrans "Velocity";
    Q.Temperature T(each stateSelect=StateSelect.prefer) = direct.thermal.T if 
      n_spec > 1 and reduceThermal "Temperature";
    // These make the selected states more readable.
*/

  protected
    Conditions.Adapters.AmagatDalton DA if n_spec > 0
      "Adapter between additivity of pressure and additivity of volume"
      annotation (Placement(transformation(extent={{80,-76},{60,-56}})));
    Connectors.DirectNode direct(
      final n_trans=n_trans,
      final inclTrans=reduceTrans,
      final inclThermal=reduceThermal) if n_spec > 0 and (reduceTrans or
      reduceThermal)
      "Connector to directly couple velocities and temperatures within the phase"
      annotation (Placement(transformation(extent={{-10,-10},{10,10}}, origin={
              26,-54}), iconTransformation(extent={{-10,-10},{10,10}}, origin={
              26,-52})));

  equation
    // Reactions
    // ---------
    // HOR
    connect(HOR.face, xPositive.'H+') annotation (Line(
        points={{84,60},{90,60},{90,5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None));
    // ORR
    connect(ORR.face, xNegative.'H+') annotation (Line(
        points={{-84,60},{-90,60},{-90,5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None));

    // Inert exchange
    // --------------
    connect('SO3-'.inter, inter) annotation (Line(
        points={{49.4,-3},{52,-6},{52,-30},{90,-30}},
        color={47,107,251},
        smooth=Smooth.None));
    connect('H+'.inter, inter) annotation (Line(
        points={{-30.6,-3},{-28,-6},{-28,-30},{90,-30}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(H2O.inter, inter) annotation (Line(
        points={{9.4,-3},{12,-6},{12,-30},{90,-30}},
        color={47,107,251},
        smooth=Smooth.None));
    connect('H+'.intra[1], H2O.intra[1]) "Electro-osmotic drag (EOD)"
      annotation (Line(
        points={{-32,-5.9},{-31,-7},{-31,-42},{9,-42},{9,-7},{8,-6},{8,-5.9}},
        color={47,107,251},
        smooth=Smooth.None));

    connect('SO3-'.intra[1], H2O.intra[2]) "SO3- and H2O" annotation (Line(
        points={{48,-5.9},{49,-7},{49,-42},{9,-42},{9,-7},{8,-6},{8,-5.9}},
        color={47,107,251},
        smooth=Smooth.None));

    connect('SO3-'.direct.translational, direct.translational) annotation (Line(
        points={{46,-8},{46,-54},{26,-54}},
        color={47,107,251},
        smooth=Smooth.None));
    connect('SO3-'.direct.thermal, direct.thermal) annotation (Line(
        points={{46,-8},{46,-54},{26,-54}},
        color={47,107,251},
        smooth=Smooth.None));
    connect('H+'.direct.translational, direct.translational) annotation (Line(
        points={{-34,-8},{-34,-54},{26,-54}},
        color={47,107,251},
        smooth=Smooth.None));
    connect('H+'.direct.thermal, direct.thermal) annotation (Line(
        points={{-34,-8},{-34,-54},{26,-54}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(H2O.direct.translational, direct.translational) annotation (Line(
        points={{6,-8},{6,-54},{26,-54}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(H2O.direct.thermal, direct.thermal) annotation (Line(
        points={{6,-8},{6,-54},{26,-54}},
        color={47,107,251},
        smooth=Smooth.None));
    connect('SO3-'.dalton, DA.dalton) annotation (Line(
        points={{43,-9.4},{43,-66},{66,-66}},
        color={47,107,251},
        smooth=Smooth.None));

    connect('H+'.dalton, DA.dalton) annotation (Line(
        points={{-37,-9.4},{-37,-66},{66,-66}},
        color={47,107,251},
        smooth=Smooth.None));

    connect(H2O.dalton, DA.dalton) annotation (Line(
        points={{3,-9.4},{3,-66},{66,-66}},
        color={47,107,251},
        smooth=Smooth.None));
    connect(DA.amagat, amagat) annotation (Line(
        points={{74,-66},{90,-66}},
        color={47,107,251},
        smooth=Smooth.None));

    // Transport
    // ---------
    // SO3-
    connect('SO3-'.faces[facesCart[Axis.x], Side.n], xNegative.'SO3-')
      annotation (Line(
        points={{40,6.10623e-16},{-90,5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('SO3-'.faces[facesCart[Axis.x], Side.p], xPositive.'SO3-')
      annotation (Line(
        points={{40,6.10623e-16},{30,5.55112e-16},{90,5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('SO3-'.faces[facesCart[Axis.y], Side.n], yNegative.'SO3-')
      annotation (Line(
        points={{40,6.10623e-16},{40,-24},{-66,-24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('SO3-'.faces[facesCart[Axis.y], Side.p], yPositive.'SO3-')
      annotation (Line(
        points={{40,6.10623e-16},{40,24},{66,24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('SO3-'.faces[facesCart[Axis.z], Side.n], zNegative.'SO3-')
      annotation (Line(
        points={{40,6.10623e-16},{52,12},{78,12}},
        color={127,127,127},
        smooth=Smooth.None));
    connect('SO3-'.faces[facesCart[Axis.z], Side.p], zPositive.'SO3-')
      annotation (Line(
        points={{40,6.10623e-16},{28,-12},{-78,-12}},
        color={127,127,127},
        smooth=Smooth.None));

    // 'H+'
    connect('H+'.faces[facesCart[Axis.x], Side.n], xNegative.'H+') annotation (
        Line(
        points={{-40,6.10623e-16},{-64,0},{-90,0},{-90,5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('H+'.faces[facesCart[Axis.x], Side.p], xPositive.'H+') annotation (
        Line(
        points={{-40,6.10623e-16},{-20,5.55112e-16},{90,5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('H+'.faces[facesCart[Axis.y], Side.n], yNegative.'H+') annotation (
        Line(
        points={{-40,6.10623e-16},{-40,-24},{-66,-24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('H+'.faces[facesCart[Axis.y], Side.p], yPositive.'H+') annotation (
        Line(
        points={{-40,6.10623e-16},{-40,24},{66,24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('H+'.faces[facesCart[Axis.z], Side.n], zNegative.'H+') annotation (
        Line(
        points={{-40,6.10623e-16},{-28,12},{78,12}},
        color={127,127,127},
        smooth=Smooth.None));
    connect('H+'.faces[facesCart[Axis.z], Side.p], zPositive.'H+') annotation (
        Line(
        points={{-40,6.10623e-16},{-52,-12},{-78,-12}},
        color={127,127,127},
        smooth=Smooth.None));

    // H2O
    connect(H2O.faces[facesCart[Axis.x], Side.n], xNegative.H2O) annotation (
        Line(
        points={{6.10623e-16,6.10623e-16},{-28,6.10623e-16},{-28,5.55112e-16},{
            -90,5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.faces[facesCart[Axis.x], Side.p], xPositive.H2O) annotation (
        Line(
        points={{6.10623e-16,6.10623e-16},{90,5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.faces[facesCart[Axis.y], Side.n], yNegative.H2O) annotation (
        Line(
        points={{6.10623e-16,6.10623e-16},{0,-24},{-66,-24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.faces[facesCart[Axis.y], Side.p], yPositive.H2O) annotation (
        Line(
        points={{6.10623e-16,6.10623e-16},{6.10623e-16,24},{66,24}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.faces[facesCart[Axis.z], Side.n], zNegative.H2O) annotation (
        Line(
        points={{6.10623e-16,6.10623e-16},{12,12},{78,12}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(H2O.faces[facesCart[Axis.z], Side.p], zPositive.H2O) annotation (
        Line(
        points={{6.10623e-16,6.10623e-16},{-12,-12},{-78,-12}},
        color={127,127,127},
        smooth=Smooth.None));

    // Phase change
    connect(H2O.physical, physical.H2O) annotation (Line(
        points={{-8.6,5},{-8.6,80.5},{-9,80.5},{-9,80}},
        color={38,196,52},
        smooth=Smooth.None));
    // Reactions
    connect(ORR.reaction, chemical) annotation (Line(
        points={{-76,60},{20,60},{20,80}},
        color={255,195,38},
        smooth=Smooth.None));

    connect(HOR.reaction, chemical) annotation (Line(
        points={{76,60},{20,60},{20,80}},
        color={255,195,38},
        smooth=Smooth.None));
    annotation (
      Documentation(info="<html><p>Assumptions:<ol>
    <li>The water in the ionomer does not participate in the reaction (only the water vapor does).</li>
    </ol</p>
    
    <p>See <a href=\"modelica://FCSys.Species.'H+'.Ionomer.Fixed\">Species.'H+'.Ionomer.Fixed</a> 
    for additional assumptions.
    For more information, see the
 <a href=\"modelica://FCSys.Phases.BaseClasses.EmptyPhase\">EmptyPhase</a> model.</p></html>"),

      Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{
              100,100}}), graphics),
      Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
              100}}), graphics));
  end Ionomer;

  model Liquid "Liquid phase"
    extends Partial(final n_spec=if inclH2O then 1 else 0);

    // Conditionally include species.
    parameter Boolean inclH2O=false "Include H2O" annotation (
      HideResult=true,
      choices(__Dymola_checkBox=true),
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>Water (H<sub>2</sub>O)</html>",
        __Dymola_joinNext=true));

    replaceable FCSys.Species.H2O.Liquid.Fixed H2O(final n_faces) if inclH2O
      constrainedby FCSys.Species.Isochoric(
      n_faces=n_faces,
      phi(each stateSelect=if reduceTrans then StateSelect.default else
            StateSelect.prefer),
      T(stateSelect=if reduceThermal then StateSelect.default else StateSelect.prefer))
      "H2O model" annotation (
      __Dymola_choicesFromPackage=true,
      Dialog(
        group="Species",
        __Dymola_descriptionLabel=true,
        __Dymola_label="<html>H<sub>2</sub>O model</html>",
        enable=inclH2O),
      Placement(transformation(extent={{-10,-10},{10,10}})));

    Connectors.FaceBus xNegative if inclFaces[Axis.x]
      "Negative face along the x axis" annotation (Placement(transformation(
            extent={{-50,-10},{-30,10}}), iconTransformation(extent={{-90,-10},
              {-70,10}})));
    Connectors.FaceBus yNegative if inclFaces[Axis.y]
      "Negative face along the y axis" annotation (Placement(transformation(
            extent={{-10,-50},{10,-30}}), iconTransformation(extent={{-10,-94},
              {10,-74}})));
    Connectors.FaceBus zNegative if inclFaces[Axis.z]
      "Negative face along the z axis" annotation (Placement(transformation(
            extent={{10,10},{30,30}}), iconTransformation(extent={{40,40},{60,
              60}})));
    Connectors.FaceBus xPositive if inclFaces[Axis.x]
      "Positive face along the x axis" annotation (Placement(transformation(
            extent={{30,-10},{50,10}}), iconTransformation(extent={{70,-10},{90,
              10}})));
    Connectors.FaceBus yPositive if inclFaces[Axis.y]
      "Positive face along the y axis" annotation (Placement(transformation(
            extent={{-10,30},{10,50}}), iconTransformation(extent={{-10,90},{10,
              110}})));
    Connectors.FaceBus zPositive if inclFaces[Axis.z]
      "Positive face along the z axis" annotation (Placement(transformation(
            extent={{-30,-30},{-10,-10}}), iconTransformation(extent={{-90,-90},
              {-70,-70}})));
    Connectors.Inter inter[n_inter](each final n_trans=n_trans) if n_spec > 0
      "Connector to exchange momentum and energy with other phases" annotation
      (Placement(transformation(extent={{30,-24},{50,-4}}), iconTransformation(
            extent={{80,-60},{100,-80}})));
    Connectors.Amagat amagat if n_spec > 0 "Connector for additivity of volume"
      annotation (Placement(transformation(extent={{30,-38},{50,-18}}),
          iconTransformation(extent={{-60,40},{-40,60}})));
    Connectors.PhysicalBus physical if inclH2O annotation (Placement(
          transformation(extent={{-19,50},{1,70}}), iconTransformation(extent={
              {20,-94},{40,-74}})));

  equation
    // Inert exchange
    connect(H2O.inter, inter) annotation (Line(
        points={{9.4,-3},{12,-6},{12,-14},{40,-14}},
        color={47,107,251},
        smooth=Smooth.None));

    // Phase change
    connect(H2O.physical, physical.H2O) annotation (Line(
        points={{-8.6,5},{-8.6,60.5},{-9,60.5},{-9,60}},
        color={38,196,52},
        smooth=Smooth.None));

    // H2O
    // ---
    // Transport
    connect(H2O.faces[Axis.x, Side.n], xNegative.H2O) annotation (Line(
        points={{6.10623e-16,6.10623e-16},{6.10623e-16,5.55112e-16},{-40,
            5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.faces[Axis.x, Side.p], xPositive.H2O) annotation (Line(
        points={{6.10623e-16,6.10623e-16},{40,5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.faces[facesCart[Axis.y], Side.n], yNegative.H2O) annotation (
        Line(
        points={{6.10623e-16,6.10623e-16},{5.55112e-16,-40}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.faces[facesCart[Axis.y], Side.p], yPositive.H2O) annotation (
        Line(
        points={{6.10623e-16,6.10623e-16},{-4.87687e-22,20},{5.55112e-16,40}},
        color={127,127,127},
        smooth=Smooth.None));

    connect(H2O.faces[facesCart[Axis.z], Side.n], zNegative.H2O) annotation (
        Line(
        points={{6.10623e-16,6.10623e-16},{20,20}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(H2O.faces[facesCart[Axis.z], Side.p], zPositive.H2O) annotation (
        Line(
        points={{6.10623e-16,6.10623e-16},{-20,-20}},
        color={127,127,127},
        smooth=Smooth.None));
    connect(amagat, H2O.amagat) annotation (Line(
        points={{40,-28},{3,-28},{3,-9.4}},
        color={47,107,251},
        smooth=Smooth.None));
    annotation (
      Documentation(info="<html>
    <p>Please see the documentation of the
 <a href=\"modelica://FCSys.Phases.BaseClasses.EmptyPhase\">EmptyPhase</a> model.</p></html>"),

      Icon(graphics),
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
              100,100}}), graphics));
  end Liquid;


  annotation (Documentation(info="
<html><p>The graphite, ionomer, and
liquid phases can only be used with a compressible phase (gas).</p>

  <p><b>Licensed by the Georgia Tech Research Corporation under the Modelica License 2</b><br>
Copyright 2007&ndash;2013, <a href=\"http://www.gtrc.gatech.edu/\">Georgia Tech Research Corporation</a>.</p>

<p><i>This Modelica package is <u>free</u> software and the use is completely at <u>your own risk</u>;
it can be redistributed and/or modified under the terms of the Modelica License 2. For license conditions (including the
disclaimer of warranty) see <a href=\"modelica://FCSys.UsersGuide.License\">
FCSys.UsersGuide.License</a> or visit <a href=\"http://www.modelica.org/licenses/ModelicaLicense2\">
http://www.modelica.org/licenses/ModelicaLicense2</a>.</i></p>
</html>"));
protected
  partial model Partial "Base model for a phase"
    import Modelica.Math.BooleanVectors.index;
    // extends FCSys.Icons.Names.Middle;

    parameter Integer n_spec(start=0) "Number of species"
      annotation (HideResult=true);
    inner parameter Integer n_inter=0
      "Number of independent couplings with other phases"
      annotation (Dialog(connectorSizing=true),HideResult=n_spec == 0);

    // Geometry
    parameter Integer n_faces(min=1, max=3) "Number of pairs of faces"
      annotation (Dialog(group="Geometry",__Dymola_label=
            "<html><i>n</i><sub>faces</sub></html>"), HideResult=true);
    // This can't be an outer parameter in Dymola 7.4.
    parameter Q.NumberAbsolute k_DT[Axis](
      each min=Modelica.Constants.small,
      each final nominal=1) = {1,1,1} if n_spec > 0
      "Scaling factor for diffusive transport" annotation (Dialog(group=
            "Geometry", __Dymola_label=
            "<html><b><i>k</i><sub>DT</sub></b></html>"));
    inner parameter Q.NumberAbsolute k_inter[n_inter]=ones(n_inter) if n_spec
       > 0 "Coupling factor for exchange with other phases" annotation (Dialog(
          group="Geometry", __Dymola_label=
            "<html><i>k</i><sub>inter</sub></html>"));

    // Assumptions
    parameter Boolean reduceTrans=false "Same velocity for all species"
      annotation (Dialog(tab="Assumptions", enable=n_spec > 1), choices(
          __Dymola_checkBox=true));
    //if n_spec > 0
    parameter Boolean reduceThermal=false "Same temperature for all species"
      annotation (Dialog(tab="Assumptions", enable=n_spec > 1), choices(
          __Dymola_checkBox=true));
    //if n_spec > 0

    outer parameter Q.Length L[Axis] if n_spec > 0 "Length" annotation (
        HideResult=true,missingInnerMessage="This model should be used within a subregion model.
");
    outer parameter Q.Area A[Axis] if n_spec > 0 "Cross-sectional area"
      annotation (HideResult=true,missingInnerMessage="This model should be used within a subregion model.
");
    // Note:  These must be public in Dymola 7.4, so HideResult is used.

  protected
    final inner parameter Q.Length Lprime[Axis]=k_DT ./ L if n_spec > 0
      "**dimension, **Effective cross-sectional area per length **rename";
    outer parameter Integer n_trans
      "Number of components of translational momentum" annotation (
        missingInnerMessage="This model should be used within a subregion model.
");
    outer parameter Integer cartTrans[:]
      "Cartesian-axis indices of the components of translational momentum"
      annotation (missingInnerMessage="This model should be used within a subregion model.
");
    outer parameter Integer facesCart[:]
      "Face-pair indices of the Cartesian axes" annotation (missingInnerMessage
        ="This model should be used within a subregion model.
");
    outer parameter Boolean inclTrans[Axis]
      "true, if each component of translational momentum is included"
      annotation (missingInnerMessage="This model should be used within a subregion model.
");
    outer parameter Boolean inclFaces[Axis]
      "true, if each pairs of faces is included" annotation (
        missingInnerMessage="This model should be used within a subregion model.
");

    outer Conditions.Environment environment "Environmental conditions";
    // This component is conditional to prevent a mathematical singularity
    // when two or more empty phases (without any species included) are
    // connected.

    annotation (
      defaultComponentPrefixes="replaceable",
      defaultComponentName="phase",
      Documentation(info="<html><p>The area fill factor (<b><i>k</i></b>) is a vector which inversely scales all
    the transport coefficients (&beta;, &zeta;, &eta;, and &theta;) of all of the species
    within the phase.  It can be used to introduce minor head loss or the effects of
    porosity or tortousity.  These effects may be anisotropic.</p>

    <p>Porosity is often quoted in material data sheets (e.g.,
    [<a href=\"modelica://FCSys.UsersGuide.References\">Toray2010</a>]) as volumetric porosity.  Using the
    Bruggeman correction factor [<a href=\"modelica://FCSys.UsersGuide.References\">Weber2004</a>, p. 4696],
    the area fill factor for the solid should be set to (1 - &epsilon;)<sup>3/2</sup>
    along each axis, where &epsilon; is the volumetric porosity (or volumetric fill factor
    of the gas).</p>


</html>"),
      Icon(graphics={
          Ellipse(
            extent={{-40,100},{40,20}},
            lineColor={127,127,127},
            startAngle=30,
            endAngle=149,
            pattern=LinePattern.Dash,
            fillPattern=FillPattern.Solid,
            fillColor={225,225,225}),
          Ellipse(
            extent={{20,-4},{100,-84}},
            lineColor={127,127,127},
            startAngle=270,
            endAngle=390,
            pattern=LinePattern.Dash,
            fillPattern=FillPattern.Solid,
            fillColor={225,225,225}),
          Ellipse(
            extent={{-100,-4},{-20,-84}},
            lineColor={127,127,127},
            startAngle=149,
            endAngle=270,
            pattern=LinePattern.Dash,
            fillPattern=FillPattern.Solid,
            fillColor={225,225,225}),
          Polygon(
            points={{60,-84},{-60,-84},{-94.5,-24},{-34.5,80},{34.5,80},{94.5,-24},
                {60,-84}},
            pattern=LinePattern.None,
            fillPattern=FillPattern.Sphere,
            smooth=Smooth.None,
            fillColor={225,225,225},
            lineColor={0,0,0}),
          Line(
            points={{-60,-84.1},{60,-84.1}},
            color={127,127,127},
            pattern=LinePattern.Dash,
            smooth=Smooth.None),
          Line(
            points={{34.5,80},{94.5,-24}},
            color={127,127,127},
            pattern=LinePattern.Dash,
            smooth=Smooth.None),
          Line(
            points={{-34.5,80},{-94.5,-24}},
            color={127,127,127},
            pattern=LinePattern.Dash,
            smooth=Smooth.None),
          Text(
            extent={{-100,-20},{100,20}},
            textString="%name",
            lineColor={0,0,0})}),
      Diagram(graphics));
  end Partial;
end Phases;

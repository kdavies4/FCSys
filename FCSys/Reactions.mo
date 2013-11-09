within FCSys;
package Reactions "Chemical reactions"
  extends Modelica.Icons.Package;

  package Examples "Examples"
    extends Modelica.Icons.ExamplesPackage;

    model Overpotential "**"
      extends Modelica.Icons.Example;
      Electrochemistry.ElectronTransfer rate(n_trans=1, I0=1e-5*U.A)
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
      Conditions.ByConnector.Chemical.Potential potential(
        inclTransY=false,
        inclTransZ=false,
        chemical(redeclare constant Integer n_trans=1)) annotation (Placement(
            transformation(
            extent={{-10,10},{10,-10}},
            rotation=90,
            origin={36,0})));
      Conditions.ByConnector.Chemical.Current current(
        inclTransY=false,
        inclTransZ=false,
        chemical(redeclare constant Integer n_trans=1),
        redeclare Modelica.Blocks.Sources.Ramp source(height=U.A, duration=1))
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={-36,0})));
      Electrochemistry.DoubleLayer capacitor(n_trans=1, useTrans=false)
        annotation (Placement(transformation(extent={{-10,10},{10,30}})));
      Conditions.ByConnector.Direct.Efforts direct(inclTransY=false, inclTransZ
          =false)
        annotation (Placement(transformation(extent={{-10,-34},{10,-54}})));
    equation
      connect(current.chemical, rate.negative) annotation (Line(
          points={{-32,0},{-6,0}},
          color={255,195,38},
          smooth=Smooth.None));
      connect(potential.chemical, rate.positive) annotation (Line(
          points={{32,0},{6,0}},
          color={255,195,38},
          smooth=Smooth.None));
      connect(capacitor.negative, rate.negative) annotation (Line(
          points={{-6,20},{-20,20},{-20,0},{-6,0}},
          color={255,195,38},
          smooth=Smooth.None));
      connect(capacitor.positive, rate.positive) annotation (Line(
          points={{6,20},{20,20},{20,0},{6,0}},
          color={255,195,38},
          smooth=Smooth.None));
      connect(direct.direct, capacitor.direct) annotation (Line(
          points={{0,-34},{0,20}},
          color={221,23,47},
          smooth=Smooth.None));
      connect(rate.direct, direct.direct) annotation (Line(
          points={{0,-4},{0,-34}},
          color={221,23,47},
          smooth=Smooth.None));
      annotation (
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics),
        experiment(StopTime=2),
        __Dymola_experimentSetupOutput);
    end Overpotential;

    model Stoichiometry
      "<html>Test the stoichiometry of the <a href=\"modelica://FCSys.Reactions.HOR\">HOR</a></html>"
      extends Modelica.Icons.Example;

      HOR hOR(n_trans=3)
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
      Conditions.ByConnector.Chemical.Potential H2(sT=1000*U.K, source(y=U.A))
        annotation (Placement(transformation(extent={{-40,-10},{-20,-30}})));
      Conditions.ByConnector.Chemical.Current 'e-'(sT=2000*U.K, redeclare
          Modelica.Blocks.Sources.Ramp source(
          height=100*U.A,
          duration=100,
          startTime=10,
          offset=U.mA))
        annotation (Placement(transformation(extent={{-10,-10},{10,-30}})));
      Conditions.ByConnector.Chemical.Potential 'H+'(sT=3000*U.K)
        annotation (Placement(transformation(extent={{20,-10},{40,-30}})));

      inner Conditions.Environment environment
        annotation (Placement(transformation(extent={{20,20},{40,40}})));
    equation
      connect(hOR.'connH+', 'H+'.chemical) annotation (Line(
          points={{4,0},{4,-8},{30,-8},{30,-16}},
          color={221,23,47},
          smooth=Smooth.None));
      connect(H2.chemical, hOR.connH2) annotation (Line(
          points={{-30,-16},{-30,-8},{-4,-8},{-4,0}},
          color={221,23,47},
          smooth=Smooth.None));
      connect('e-'.chemical, hOR.'conne-') annotation (Line(
          points={{0,-16},{0,0}},
          color={221,23,47},
          smooth=Smooth.None));
      annotation (
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics),
        Commands(file=
              "Resources/Scripts/Dymola/Reactions.Examples.Stoichiometry.mos"
            "Reactions.Examples.Stoichiometry.mos"),
        experiment(StopTime=200, Tolerance=1e-006),
        __Dymola_experimentSetupOutput);
    end Stoichiometry;

  end Examples;

  model HOR "Hydrogen oxidation reaction"
    extends FCSys.Icons.Names.Top2;

    constant Integer n_trans(min=0, max=3)
      "Number of components of translational momentum" annotation (Dialog(
          __Dymola_label="<html><i>n</i><sub>trans</sub></html>"));
    // Note:  This must be a constant rather than a parameter due to errors
    // in Dymola 2014.
    Conditions.Adapters.ChemicalReaction 'e-'(
      final n_trans=n_trans,
      m=Characteristics.'e-'.Gas.m,
      n=-2) annotation (Placement(transformation(
          extent={{-10,10},{10,-10}},
          rotation=180,
          origin={10,10})));
    Conditions.Adapters.ChemicalReaction 'H+'(
      final n_trans=n_trans,
      m=Characteristics.'H+'.Gas.m,
      n=-2) annotation (Placement(transformation(extent={{20,0},{0,-20}})));
    Conditions.Adapters.ChemicalReaction H2(
      final n_trans=n_trans,
      m=Characteristics.H2.Gas.m,
      n=1) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={-10,0})));

    Connectors.Chemical 'conne-'(redeclare final constant Integer n_trans=
          n_trans) "Connector for e-" annotation (Placement(transformation(
            extent={{20,0},{40,20}}), iconTransformation(extent={{-10,-10},{10,
              10}})));
    Connectors.Chemical 'connH+'(redeclare final constant Integer n_trans=
          n_trans) "Connector for H+" annotation (Placement(transformation(
            extent={{20,-20},{40,0}}), iconTransformation(extent={{30,-10},{50,
              10}})));
    Connectors.Chemical connH2(redeclare final constant Integer n_trans=n_trans)
      "Connector for H2" annotation (Placement(transformation(extent={{-40,-10},
              {-20,10}}), iconTransformation(extent={{-50,-10},{-30,10}})));
    // Note:  These redeclarations are necessary due to errors in Dymola 2014.

  equation
    connect(connH2, H2.chemical) annotation (Line(
        points={{-30,0},{-14,0}},
        color={255,195,38},
        smooth=Smooth.None));
    connect('e-'.chemical, 'conne-') annotation (Line(
        points={{14,10},{30,10}},
        color={255,195,38},
        smooth=Smooth.None));
    connect('H+'.chemical, 'connH+') annotation (Line(
        points={{14,-10},{30,-10}},
        color={255,195,38},
        smooth=Smooth.None));
    connect(H2.reaction, 'e-'.reaction) annotation (Line(
        points={{-6,0},{0,0},{0,10},{6,10}},
        color={255,195,38},
        smooth=Smooth.None));
    connect('H+'.reaction, H2.reaction) annotation (Line(
        points={{6,-10},{0,-10},{0,0},{-6,0}},
        color={255,195,38},
        smooth=Smooth.None));
    annotation (
      defaultComponentName="HOR",
      Documentation(info=
            "<html><table border=1><tr><td><font size=100 color=\"gray\">H<sub>2</sub> &#8652; 2e<sup>-</sup> + 2H<sup>+</sup></td></tr></table></div></html>"),

      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
              100}}), graphics={Rectangle(
              extent={{-100,40},{100,-50}},
              pattern=LinePattern.Dash,
              lineColor={127,127,127},
              radius=15,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),Bitmap(extent={{-100,-20},{100,-40}},
            fileName=
            "modelica://FCSys/Resources/Documentation/Reactions/HOR.png")}),
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-40,-20},{40,
              20}}), graphics));
  end HOR;

  model ORR "Oxygen reduction reaction"
    extends FCSys.Icons.Names.Top2;

    constant Integer n_trans(min=0, max=3)
      "Number of components of translational momentum" annotation (Dialog(
          __Dymola_label="<html><i>n</i><sub>trans</sub></html>"));
    // Note:  This must be a constant rather than a parameter due to errors
    // in Dymola 2014.
    Conditions.Adapters.ChemicalReaction 'e-'(
      final n_trans=n_trans,
      m=Characteristics.'e-'.Gas.m,
      n=4) annotation (Placement(transformation(
          extent={{10,10},{-10,-10}},
          rotation=180,
          origin={-10,10})));
    Conditions.Adapters.ChemicalReaction 'H+'(
      final n_trans=n_trans,
      m=Characteristics.'H+'.Gas.m,
      n=4) annotation (Placement(transformation(extent={{-20,-20},{0,0}})));
    Conditions.Adapters.ChemicalReaction O2(
      final n_trans=n_trans,
      m=Characteristics.O2.Gas.m,
      n=1) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={-10,-30})));
    Conditions.Adapters.ChemicalReaction H2O(
      final n_trans=n_trans,
      m=Characteristics.H2O.Gas.m,
      n=-2) annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=0,
          origin={10,-10})));

    Connectors.Chemical 'conne-'(redeclare final constant Integer n_trans=
          n_trans) "Connector for e-" annotation (Placement(transformation(
            extent={{-40,0},{-20,20}}), iconTransformation(extent={{-70,-10},{-50,
              10}})));
    Connectors.Chemical 'connH+'(redeclare final constant Integer n_trans=
          n_trans) "Connector for H+" annotation (Placement(transformation(
            extent={{-40,-20},{-20,0}}), iconTransformation(extent={{-30,-10},{
              -10,10}})));
    Connectors.Chemical connO2(redeclare final constant Integer n_trans=n_trans)
      "Connector for O2" annotation (Placement(transformation(extent={{-40,-40},
              {-20,-20}}), iconTransformation(extent={{10,-10},{30,10}})));
    Connectors.Chemical connH2O(redeclare final constant Integer n_trans=
          n_trans) "Connector for H2O" annotation (Placement(transformation(
            extent={{20,-20},{40,0}}), iconTransformation(extent={{50,-10},{70,
              10}})));
    // Note:  These redeclarations are necessary due to errors in Dymola 2014.

  equation
    connect('H+'.chemical, 'connH+') annotation (Line(
        points={{-14,-10},{-30,-10}},
        color={255,195,38},
        smooth=Smooth.None));
    connect('e-'.chemical, 'conne-') annotation (Line(
        points={{-14,10},{-30,10}},
        color={255,195,38},
        smooth=Smooth.None));
    connect(O2.chemical, connO2) annotation (Line(
        points={{-14,-30},{-30,-30}},
        color={255,195,38},
        smooth=Smooth.None));
    connect(H2O.chemical, connH2O) annotation (Line(
        points={{14,-10},{30,-10}},
        color={255,195,38},
        smooth=Smooth.None));
    connect('e-'.reaction, H2O.reaction) annotation (Line(
        points={{-6,10},{0,10},{0,-10},{6,-10}},
        color={255,195,38},
        smooth=Smooth.None));
    connect('H+'.reaction, H2O.reaction) annotation (Line(
        points={{-6,-10},{6,-10}},
        color={255,195,38},
        smooth=Smooth.None));
    connect(O2.reaction, H2O.reaction) annotation (Line(
        points={{-6,-30},{0,-30},{0,-10},{6,-10}},
        color={255,195,38},
        smooth=Smooth.None));
    annotation (
      defaultComponentName="ORR",
      Documentation(info=
            "<html><table border=1><tr><td><font size=100 color=\"gray\">4e<sup>-</sup> + 4H<sup>+</sup> + O<sub>2</sub> &#8652; 2H<sub>2</sub>O</td></tr></table></div></html>"),

      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
              100}}), graphics={Rectangle(
              extent={{-100,40},{100,-50}},
              pattern=LinePattern.Dash,
              lineColor={127,127,127},
              radius=15,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),Bitmap(extent={{-100,-20},{100,-40}},
            fileName=
            "modelica://FCSys/Resources/Documentation/Reactions/ORR.png")}),
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-40,-40},{40,
              20}}), graphics));
  end ORR;

  package Electrochemistry "Models related to electrochemical reactions"
    extends Modelica.Icons.Package;

    model DoubleLayer "Electrolytic double layer"
      extends FCSys.Icons.Names.Top2;

      parameter Integer n_trans(min=1,max=3)
        "Number of components of translational momentum" annotation (Evaluate=
            true, Dialog(__Dymola_label="<html><i>n</i><sub>trans</sub></html>"));
      parameter Q.Capacitance C=1e11*U.m*U.epsilon_0 "Capacitance";
      replaceable package Data = Characteristics.'e-'.Graphite constrainedby
        Characteristics.BaseClasses.Characteristic "Material properties"
        annotation (
        Evaluate=true,
        choicesAllMatching=true,
        __Dymola_choicesFromPackage=true);

      parameter Boolean useTrans=true
        "<html>Material exits at the velocity of the <code>direct</code> connector</html>"
        annotation (
        Evaluate=true,
        Dialog(tab="Assumptions", compact=true),
        choices(__Dymola_checkBox=true));

      Q.Potential w(stateSelect=StateSelect.always, start=0)
        "Electrical potential";
      Q.Current I "Material current";

      output Q.Amount Z(stateSelect=StateSelect.never) = C*w if environment.analysis
        "Amount of charge shifted in the positive direction";

      Connectors.Chemical negative(final n_trans=n_trans)
        "Material connector on the 1st side" annotation (Placement(
            transformation(extent={{-70,-10},{-50,10}}), iconTransformation(
              extent={{-70,-10},{-50,10}})));
      Connectors.Chemical positive(final n_trans=n_trans)
        "Material connector on the 2nd side" annotation (Placement(
            transformation(extent={{50,-10},{70,10}}), iconTransformation(
              extent={{50,-10},{70,10}})));
      Connectors.Direct direct(final n_trans=n_trans)
        "Interface with the substrate" annotation (Placement(transformation(
              extent={{-10,-10},{10,10}}), iconTransformation(extent={{-10,-10},
                {10,10}})));

      outer Conditions.Environment environment "Environmental conditions";

    equation
      // Streams
      if useTrans then
        negative.phi = direct.trans.phi;
        positive.phi = direct.trans.phi;
      else
        negative.phi = inStream(positive.phi);
        positive.phi = inStream(negative.phi);
      end if;
      negative.sT = inStream(positive.sT);
      positive.sT = inStream(negative.sT);

      // Aliases
      Data.z*w = positive.g - negative.g;
      I = positive.Ndot;

      // Conservation
      0 = negative.Ndot + positive.Ndot "Material (no storage)";
      zeros(n_trans) = Data.m*(actualStream(negative.phi) - actualStream(
        positive.phi))*I + direct.trans.mPhidot
        "Translational momentum (no storage)";
      der(C*w)/U.s = Data.z*I
        "Electrical energy (reversible; simplified using material balance and divided by potential)";
      0 = direct.therm.Qdot + (actualStream(negative.phi)*actualStream(negative.phi)
         - actualStream(positive.phi)*actualStream(positive.phi))*I*Data.m/2 +
        direct.trans.phi*direct.trans.mPhidot
        "Mechanical and thermal energy (no storage)";

      annotation (
        Documentation(info="<html><p>If <code>useTrans</code> is <code>true</code>, then the material exits with the
  velocity of the <code>direct</code> connector.  Typically, that connector should be connected to the stationary solid,
  in which case heat will be generated if material arrives with a nonzero velocity.  That heat is rejected to the <code>direct</code> connector.
  </p></html>"),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics),
        Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
                100,100}}), graphics={Line(
                  points={{-20,30},{-20,-30}},
                  color={255,195,38},
                  smooth=Smooth.None),Line(
                  points={{20,30},{20,-30}},
                  color={255,195,38},
                  smooth=Smooth.None),Line(
                  points={{-50,0},{-20,0}},
                  color={255,195,38},
                  smooth=Smooth.None),Line(
                  points={{20,0},{50,0}},
                  color={255,195,38},
                  smooth=Smooth.None)}));
    end DoubleLayer;

    model ElectronTransfer "Electron transfer"
      import Modelica.Math.asinh;
      extends FCSys.Icons.Names.Top2;

      parameter Integer n_trans(min=1,max=3)
        "Number of components of translational momentum" annotation (Evaluate=
            true, Dialog(__Dymola_label="<html><i>n</i><sub>trans</sub></html>"));

      parameter Integer z=-1 "Charge number";
      parameter Q.Potential E_A=0 "Activation energy" annotation (Dialog(group=
              "Chemical parameters", __Dymola_label=
              "<html><i>E</i><sub>A</sub></html>"));
      parameter Q.NumberAbsolute alpha(max=1) = 0.5
        "Charge transfer coefficient" annotation (Dialog(group=
              "Chemical parameters", __Dymola_label="<html>&alpha;</html>"));
      parameter Q.Temperature T=300*U.K;
      parameter Q.Current I0=U.A
        "Exchange current (including the activation reference)";

      Connectors.Chemical negative(final n_trans=n_trans)
        "Material connector on the 1st side" annotation (Placement(
            transformation(extent={{-70,-10},{-50,10}}), iconTransformation(
              extent={{-70,-10},{-50,10}})));
      Connectors.Chemical positive(final n_trans=n_trans)
        "Material connector on the 2nd side" annotation (Placement(
            transformation(extent={{50,-10},{70,10}}), iconTransformation(
              extent={{50,-10},{70,10}})));

      // Aliases
      Q.Current I "Reaction rate";
      Q.Potential Deltag "Potential difference";

      Connectors.Direct direct(final n_trans=n_trans) annotation (Placement(
            transformation(extent={{-10,-50},{10,-30}}), iconTransformation(
              extent={{-10,-50},{10,-30}})));

    equation
      // Aliases (only to clarify and simplify other equations)
      I = positive.Ndot;
      Deltag = positive.g - negative.g;

      // Streams
      negative.phi = inStream(positive.phi);
      positive.phi = inStream(negative.phi);
      negative.sT = inStream(positive.sT);
      positive.sT = inStream(negative.sT);

      // Reaction rate
      I*exp(E_A/T) = I0*(exp((1 - alpha)*Deltag/T) - exp(-alpha*Deltag/T));

      // Conservation (without storage)
      0 = negative.Ndot + positive.Ndot "Material";
      zeros(n_trans) = direct.trans.mPhidot "Translational momentum";
      0 = Deltag*I + direct.therm.Qdot "Energy";
      // Note:  Energy and momentum cancel among the stream terms.

      annotation (
        defaultComponentName="'e-Transfer'",
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics),
        Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
                100,100}}), graphics={Line(
                  points={{0,-20},{0,-50}},
                  color={221,23,47},
                  smooth=Smooth.None),Line(
                  points={{-50,0},{50,0}},
                  color={255,195,38},
                  smooth=Smooth.None),Rectangle(
                  extent={{-30,20},{32,-20}},
                  lineColor={255,195,38},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid),Line(
                  points={{-20,4},{20,4},{8,12}},
                  color={255,195,38},
                  smooth=Smooth.None),Line(
                  points={{-20,-5},{20,-5},{8,3}},
                  color={255,195,38},
                  smooth=Smooth.None,
                  origin={0,-11},
                  rotation=180)}));
    end ElectronTransfer;
  end Electrochemistry;
end Reactions;

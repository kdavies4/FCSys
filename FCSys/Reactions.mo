within FCSys;
package Reactions "Electrochemical reactions"
  extends Modelica.Icons.Package;

  model HOR "Hydrogen oxidation reaction"
    extends FCSys.Icons.Names.Top2;

    constant Integer n_trans(min=0, max=3)
      "Number of components of translational momentum" annotation (Dialog(
          __Dymola_label="<html><i>n</i><sub>trans</sub></html>"));
    // Note:  This must be a constant rather than a parameter due to errors
    // in Dymola 2014.
    Conditions.Adapters.ElectrochemicalReaction 'e-'(
      final n_trans=n_trans,
      m=Characteristics.'e-'.Gas.m,
      n=-2) annotation (Placement(transformation(
          extent={{-10,10},{10,-10}},
          rotation=180,
          origin={10,10})));
    Conditions.Adapters.ElectrochemicalReaction 'H+'(
      final n_trans=n_trans,
      m=Characteristics.'H+'.Gas.m,
      n=-2) annotation (Placement(transformation(extent={{20,0},{0,-20}})));
    Conditions.Adapters.ElectrochemicalReaction H2(
      final n_trans=n_trans,
      m=Characteristics.H2.Gas.m,
      n=1) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={-10,0})));

    Connectors.Electrochemical 'conne-'(redeclare final constant Integer
        n_trans=n_trans) "Connector for e-" annotation (Placement(
          transformation(extent={{20,0},{40,20}}), iconTransformation(extent={{
              -10,-10},{10,10}})));
    Connectors.Electrochemical 'connH+'(redeclare final constant Integer
        n_trans=n_trans) "Connector for H+" annotation (Placement(
          transformation(extent={{20,-20},{40,0}}), iconTransformation(extent={
              {30,-10},{50,10}})));
    Connectors.Electrochemical connH2(redeclare final constant Integer n_trans=
          n_trans) "Connector for H2" annotation (Placement(transformation(
            extent={{-40,-10},{-20,10}}), iconTransformation(extent={{-50,-10},
              {-30,10}})));
    // Note:  These redeclarations are necessary due to errors in Dymola 2014.

  equation
    connect(connH2, H2.electrochemical) annotation (Line(
        points={{-30,0},{-14,0}},
        color={221,23,47},
        smooth=Smooth.None));
    connect('e-'.electrochemical, 'conne-') annotation (Line(
        points={{14,10},{30,10}},
        color={221,23,47},
        smooth=Smooth.None));
    connect('H+'.electrochemical, 'connH+') annotation (Line(
        points={{14,-10},{30,-10}},
        color={221,23,47},
        smooth=Smooth.None));
    connect(H2.reaction, 'e-'.reaction) annotation (Line(
        points={{-6,0},{0,0},{0,10},{6,10}},
        color={221,23,47},
        smooth=Smooth.None));
    connect('H+'.reaction, H2.reaction) annotation (Line(
        points={{6,-10},{0,-10},{0,0},{-6,0}},
        color={221,23,47},
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
    Conditions.Adapters.ElectrochemicalReaction 'e-'(
      final n_trans=n_trans,
      m=Characteristics.'e-'.Gas.m,
      n=4) annotation (Placement(transformation(
          extent={{10,10},{-10,-10}},
          rotation=180,
          origin={-10,10})));
    Conditions.Adapters.ElectrochemicalReaction 'H+'(
      final n_trans=n_trans,
      m=Characteristics.'H+'.Gas.m,
      n=4) annotation (Placement(transformation(extent={{-20,-20},{0,0}})));
    Conditions.Adapters.ElectrochemicalReaction O2(
      final n_trans=n_trans,
      m=Characteristics.O2.Gas.m,
      n=1) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={-10,-30})));
    Conditions.Adapters.ElectrochemicalReaction H2O(
      final n_trans=n_trans,
      m=Characteristics.H2O.Gas.m,
      n=-2) annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=0,
          origin={10,-10})));

    Connectors.Electrochemical 'conne-'(redeclare final constant Integer
        n_trans=n_trans) "Connector for e-" annotation (Placement(
          transformation(extent={{-40,0},{-20,20}}), iconTransformation(extent=
              {{-70,-10},{-50,10}})));
    Connectors.Electrochemical 'connH+'(redeclare final constant Integer
        n_trans=n_trans) "Connector for H+" annotation (Placement(
          transformation(extent={{-40,-20},{-20,0}}), iconTransformation(extent
            ={{-30,-10},{-10,10}})));
    Connectors.Electrochemical connO2(redeclare final constant Integer n_trans=
          n_trans) "Connector for O2" annotation (Placement(transformation(
            extent={{-40,-40},{-20,-20}}), iconTransformation(extent={{10,-10},
              {30,10}})));
    Connectors.Electrochemical connH2O(redeclare final constant Integer n_trans
        =n_trans) "Connector for H2O" annotation (Placement(transformation(
            extent={{20,-20},{40,0}}), iconTransformation(extent={{50,-10},{70,
              10}})));
    // Note:  These redeclarations are necessary due to errors in Dymola 2014.

  equation
    connect('H+'.electrochemical, 'connH+') annotation (Line(
        points={{-14,-10},{-30,-10}},
        color={221,23,47},
        smooth=Smooth.None));
    connect('e-'.electrochemical, 'conne-') annotation (Line(
        points={{-14,10},{-30,10}},
        color={221,23,47},
        smooth=Smooth.None));
    connect(O2.electrochemical, connO2) annotation (Line(
        points={{-14,-30},{-30,-30}},
        color={221,23,47},
        smooth=Smooth.None));
    connect(H2O.electrochemical, connH2O) annotation (Line(
        points={{14,-10},{30,-10}},
        color={221,23,47},
        smooth=Smooth.None));
    connect('e-'.reaction, H2O.reaction) annotation (Line(
        points={{-6,10},{0,10},{0,-10},{6,-10}},
        color={221,23,47},
        smooth=Smooth.None));
    connect('H+'.reaction, H2O.reaction) annotation (Line(
        points={{-6,-10},{6,-10}},
        color={221,23,47},
        smooth=Smooth.None));
    connect(O2.reaction, H2O.reaction) annotation (Line(
        points={{-6,-30},{0,-30},{0,-10},{6,-10}},
        color={221,23,47},
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

  package Examples "Examples"
    extends Modelica.Icons.ExamplesPackage;

    model Stoichiometry
      "<html>Test the stoichiometry of the <a href=\"modelica://FCSys.Reactions.HOR\">HOR</a></html>"
      extends Modelica.Icons.Example;

      HOR hOR(n_trans=3)
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
      Conditions.ByConnector.Electrochemical.Potential H2(sT=1000*U.K, source(y
            =U.A))
        annotation (Placement(transformation(extent={{-40,-10},{-20,-30}})));
      Conditions.ByConnector.Electrochemical.Current 'e-'(sT=2000*U.K,
          redeclare Modelica.Blocks.Sources.Ramp source(
          height=100*U.A,
          duration=100,
          startTime=10,
          offset=U.mA))
        annotation (Placement(transformation(extent={{-10,-10},{10,-30}})));
      Conditions.ByConnector.Electrochemical.Potential 'H+'(sT=3000*U.K)
        annotation (Placement(transformation(extent={{20,-10},{40,-30}})));

      inner Conditions.Environment environment
        annotation (Placement(transformation(extent={{20,20},{40,40}})));
    equation
      connect(hOR.'connH+', 'H+'.electrochemical) annotation (Line(
          points={{4,0},{4,-8},{30,-8},{30,-16}},
          color={221,23,47},
          smooth=Smooth.None));
      connect(H2.electrochemical, hOR.connH2) annotation (Line(
          points={{-30,-16},{-30,-8},{-4,-8},{-4,0}},
          color={221,23,47},
          smooth=Smooth.None));
      connect('e-'.electrochemical, hOR.'conne-') annotation (Line(
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
end Reactions;

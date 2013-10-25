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
        color={221,23,47},
        smooth=Smooth.None));
    connect('e-'.chemical, 'conne-') annotation (Line(
        points={{14,10},{30,10}},
        color={221,23,47},
        smooth=Smooth.None));
    connect('H+'.chemical, 'connH+') annotation (Line(
        points={{14,-10},{30,-10}},
        color={221,23,47},
        smooth=Smooth.None));
    connect(H2.reaction, 'e-'.reaction) annotation (Line(
        points={{-6,0},{0,0},{0,10},{6,10}},
        color={221,23,47},
        smooth=Smooth.None));
    connect(H2.reaction, 'H+'.reaction) annotation (Line(
        points={{-6,0},{0,0},{0,-10},{6,-10}},
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
            fillPattern=FillPattern.Solid), Bitmap(extent={{-100,-20},{100,-40}},
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
        color={221,23,47},
        smooth=Smooth.None));
    connect('e-'.chemical, 'conne-') annotation (Line(
        points={{-14,10},{-30,10}},
        color={221,23,47},
        smooth=Smooth.None));
    connect(O2.chemical, connO2) annotation (Line(
        points={{-14,-30},{-30,-30}},
        color={221,23,47},
        smooth=Smooth.None));
    connect(H2O.chemical, connH2O) annotation (Line(
        points={{14,-10},{30,-10}},
        color={221,23,47},
        smooth=Smooth.None));
    connect('e-'.reaction, 'H+'.reaction) annotation (Line(
        points={{-6,10},{0,10},{0,-10},{-6,-10}},
        color={221,23,47},
        smooth=Smooth.None));
    connect(H2O.reaction, O2.reaction) annotation (Line(
        points={{6,-10},{0,-10},{0,-30},{-6,-30}},
        color={221,23,47},
        smooth=Smooth.None));
    connect(H2O.reaction, 'e-'.reaction) annotation (Line(
        points={{6,-10},{0,-10},{0,10},{-6,10}},
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
            fillPattern=FillPattern.Solid), Bitmap(extent={{-100,-20},{100,-40}},
              fileName=
                "modelica://FCSys/Resources/Documentation/Reactions/ORR.png")}),

      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-40,-40},{40,
              20}}), graphics));

  end ORR;

end Reactions;

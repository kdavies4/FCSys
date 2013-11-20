within FCSys;
package Chemistry "Models associated with chemical reactions"
  extends Icons.ChemistryPackage;

  package Examples "Examples"
    extends Modelica.Icons.ExamplesPackage;

    model Overpotential "Demonstrate the Butler-Volmer overpotential"
      extends Modelica.Icons.Example;
      extends Modelica.Icons.UnderConstruction;

      output Q.Potential w=-'e-Transfer'.Deltag "Overpotential";
      output Q.Current I_A=-'e-Transfer'.I/U.A if environment.analysis
        "Reaction current in amperes";

      Electrochemistry.ElectronTransfer 'e-Transfer'(
        redeclare constant Integer n_trans=1,
        fromI=false,
        I0=U.mA)
        annotation (Placement(transformation(extent={{-10,10},{10,-10}})));
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
        redeclare Modelica.Blocks.Sources.Sine set(
          freqHz=1,
          amplitude=100*U.A,
          phase=1.5707963267949)) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={-36,0})));
      Electrochemistry.DoubleLayer doubleLayer(
        setVelocity=false,
        inclVolume=false,
        redeclare constant Integer n_trans=1,
        w(fixed=true))
        annotation (Placement(transformation(extent={{-10,20},{10,40}})));
      Conditions.ByConnector.Inter.Efforts substrate(inclTransY=false,
          inclTransZ=false)
        annotation (Placement(transformation(extent={{50,10},{70,-10}})));

      inner Conditions.Environment environment
        annotation (Placement(transformation(extent={{40,40},{60,60}})));
    equation
      connect(current.chemical, 'e-Transfer'.negative) annotation (Line(
          points={{-32,0},{-6,0}},
          color={255,195,38},
          smooth=Smooth.None));
      connect(potential.chemical, 'e-Transfer'.positive) annotation (Line(
          points={{32,0},{6,0}},
          color={255,195,38},
          smooth=Smooth.None));
      connect(doubleLayer.negative, 'e-Transfer'.negative) annotation (Line(
          points={{-6,30},{-20,30},{-20,0},{-6,0}},
          color={255,195,38},
          smooth=Smooth.None));
      connect(doubleLayer.positive, 'e-Transfer'.positive) annotation (Line(
          points={{6,30},{20,30},{20,0},{6,0}},
          color={255,195,38},
          smooth=Smooth.None));
      connect(substrate.inter, doubleLayer.intra) annotation (Line(
          points={{60,10},{60,20},{0,20},{0,26}},
          color={221,23,47},
          smooth=Smooth.None));
      connect(substrate.inter, 'e-Transfer'.intra) annotation (Line(
          points={{60,10},{60,20},{0,20},{0,4}},
          color={221,23,47},
          smooth=Smooth.None));
      annotation (
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics),
        experiment(StopTime=2, __Dymola_NumberOfIntervals=5000),
        __Dymola_experimentSetupOutput,
        Commands(file=
              "Resources/Scripts/Dymola/Chemistry.Examples.Overpotential.mos"
            "Chemistry.Examples.Overpotential.mos"));
    end Overpotential;

    model Stoichiometry
      "<html>Test the stoichiometry of the <a href=\"modelica://FCSys.Chemistry.HOR\">HOR</a></html>"
      extends Modelica.Icons.Example;

      HOR hOR(n_trans=3)
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
      Conditions.ByConnector.Chemical.Potential H2(sT=1000*U.K, set(y=U.A))
        annotation (Placement(transformation(extent={{-40,-10},{-20,-30}})));
      Conditions.ByConnector.Chemical.Current 'e-'(sT=2000*U.K, redeclare
          Modelica.Blocks.Sources.Ramp set(
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
          color={255,195,38},
          smooth=Smooth.None));
      connect(H2.chemical, hOR.connH2) annotation (Line(
          points={{-30,-16},{-30,-8},{-4,-8},{-4,0}},
          color={255,195,38},
          smooth=Smooth.None));
      connect('e-'.chemical, hOR.'conne-') annotation (Line(
          points={{0,-16},{0,0}},
          color={255,195,38},
          smooth=Smooth.None));
      annotation (
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics),
        Commands(file=
              "Resources/Scripts/Dymola/Chemistry.Examples.Stoichiometry.mos"
            "Reactions.Examples.Stoichiometry.mos"),
        experiment(StopTime=200),
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
      n=4,
      reaction(Ndot(stateSelect=StateSelect.prefer))) annotation (Placement(
          transformation(
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
            fillPattern=FillPattern.Solid), Bitmap(extent={{-100,-20},{100,-40}},
              fileName=
                "modelica://FCSys/Resources/Documentation/Reactions/ORR.png")}),

      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-40,-40},{40,
              20}}), graphics));
  end ORR;

public
  model SurfaceTension "Surface tension"

    extends FCSys.Icons.Names.Top2;

    // Geometric parameters
    constant Integer n_trans(min=0, max=3)
      "Number of components of translational momentum" annotation (Dialog(group
          ="Geometry", __Dymola_label="<html><i>n</i><sub>trans</sub></html>"));
    // Note:  This must be a constant rather than a parameter due to errors
    // in Dymola 2014.
    Q.LengthReciprocal overR=1/U.mm "Reciprocal of characteristic radius"
      annotation (Dialog(group="Geometry", __Dymola_label=
            "<html>1/<i>R</i></html>"));

    // Material properties
    parameter Q.VolumeSpecific v_w=Characteristics.H2O.Liquid.v_Tp()
      "Specific volume of the wetting phase" annotation (Dialog(group=
            "Material properties", __Dymola_label=
            "<html><i>v</i><sub>w</sub></html>"));
    parameter Q.SurfaceTension gamma=0.0663*U.N/U.m "Surface tension"
      annotation (Dialog(group="Material properties", __Dymola_label=
            "<html>&gamma;</html>"));

    // Aliases
    Q.Pressure Deltap "Capillary pressure";

    // Auxiliary variables (for analysis only)
    output Q.Potential Deltag=wetting.g - nonwetting.g if environment.analysis
      "Potential of wetting phase w.r.t. nonwetting phase";

    Connectors.Chemical wetting(final n_trans=n_trans)
      "Chemical interface to the wetting phase" annotation (Placement(
          transformation(extent={{-10,-10},{10,10}}), iconTransformation(extent
            ={{-10,-10},{10,10}})));
    Connectors.Chemical nonwetting(final n_trans=n_trans)
      "Chemical interface to the nonwetting phase" annotation (Placement(
          transformation(extent={{50,-10},{70,10}}), iconTransformation(extent=
              {{50,-10},{70,10}})));

  protected
    outer Conditions.Environment environment "Environmental conditions";

  equation
    // Aliases
    Deltap = 2*gamma*overR "Young-Laplace equation";

    // Streams
    wetting.phi = inStream(nonwetting.phi[:]);
    nonwetting.phi = inStream(wetting.phi[:]);
    wetting.sT = inStream(nonwetting.sT);
    nonwetting.sT = inStream(wetting.sT);

    // Potentials
    nonwetting.g = wetting.g + v_w*Deltap "Kelvin equation (see assumption #1)";

    // Conservation (without storage)
    0 = wetting.Ndot + nonwetting.Ndot "Material";

    annotation (
      Documentation(info="<html>
    <p>The characteristic radius (<i>R</i>) is the harmonic mean of the (2) principle radii of the liquid volume.</p>
    
    <p>The default surface tension (&gamma; = 0.0663 N/m) is for saturated water at 60 &deg;C, interpolated from
    [<a href=\"modelica://FCSys.UsersGuide.References\">Incropera2002</a>, pp. 924].  Note that the surface tension in
    [<a href=\"modelica://FCSys.UsersGuide.References\">Wang2001</a>] is incorrect (likely unit conversion error).</p>  
    
    <p>Assumptions:<ol>
    <li>The fluid is incompressible and isothermal across the surface layer.<li></p>
    
    </html>"),
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
              100}}), graphics={Ellipse(
            extent={{30,30},{-30,-30}},
            lineColor={47,107,251},
            fillColor={255,255,255},
            fillPattern=FillPattern.Sphere), Ellipse(extent={{30,30},{-30,-30}},
              lineColor={0,0,0})}),
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
              100,100}}), graphics));
  end SurfaceTension;

  model SurfaceTensionContact "Surface tension, with contact angle"
    extends SurfaceTension(final overR=cos(theta)/(2*sqrt(kappa)));
    // ** *J(V_w/V)

    // Geometry
    parameter Q.Volume V=U.cc "Total volume" annotation (Dialog(group=
            "Geometry", __Dymola_label="<html><i>V</i></html>"));
    Q.Volume V_w=U.cc "Volume of the wetting phase" annotation (Dialog(group=
            "Geometry", __Dymola_label="<html><i>V</i><sub>w</sub></html>"));

    // Material properties
    parameter Q.Angle theta=140*U.degree "Contact angle" annotation (Dialog(
          group="Material properties", __Dymola_label="<html>&theta;</html>"));
    parameter Q.Area kappa=6.46e-5*U.mm^2 "Permeability" annotation (Dialog(
          group="Material properties", __Dymola_label="<html>&kappa;</html>"));

    replaceable function J = FCSys.Characteristics.H2O.J "Leverett J function"
      annotation (choicesAllMatching=true, Dialog(group="Material properties"));

    annotation (
      defaultComponentName="surfaceTension",
      Documentation(info="<html>
    
    **Note modified Leverett correlation
    
    <p>The default permeability (&kappa; = 6.46&times;10<sup>-5</sup> mm<sup>2</sup>) is based on 
    the air permeability of SGL Carbon Group Sigracet&reg; 10 BA 
    [<a href=\"modelica://FCSys.UsersGuide.References\">SGL2007</a>].  
    Wang et al. use &kappa; = 10<sup>-5</sup> mm<sup>2</sup>
    [<a href=\"modelica://FCSys.UsersGuide.References\">Wang2001</a>].</p>
  
    <p>The default contact angle (&theta; = 140&deg;) is typical of the GDL measurments listed at 
    <a href=\"http://www.chem.mtu.edu/cnlm/research/Movement_of_Water-in_Fuel_Cell_Electrodes.htm\">http://www.chem.mtu.edu/cnlm/research/Movement_of_Water-in_Fuel_Cell_Electrodes.htm</a> 
    (accessed Nov. 22, 2103).</p></html>"),
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
              100}}), graphics={
          Rectangle(
            extent={{40,-40},{-40,-20}},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            pattern=LinePattern.None),
          Line(points={{-40,-20},{40,-20}}, color={0,0,0}),
          Line(points={{0,-40},{20,-20}}, color={0,0,0}),
          Line(points={{-20,-40},{0,-20}}, color={0,0,0}),
          Line(points={{20,-40},{40,-20}}, color={0,0,0}),
          Line(points={{-40,-40},{-20,-20}}, color={0,0,0})}));
  end SurfaceTensionContact;

  package Electrochemistry "Models associated with electrochemical reactions"
    extends Modelica.Icons.Package;

    model DoubleLayer "Electrolytic double layer"
      extends FCSys.Icons.Names.Top2;

      parameter Integer n_trans(min=1,max=3)
        "Number of components of translational momentum" annotation (Evaluate=
            true,Dialog(__Dymola_label="<html><i>n</i><sub>trans</sub></html>"));
      parameter Q.Area A=10*U.m^2 "Surface area"
        annotation (Dialog(__Dymola_label="<html><i>A</i></html>"));
      parameter Q.Length L=1e-10*U.m "Length of the gap"
        annotation (Dialog(__Dymola_label="<html><i>L</i></html>"));
      parameter Q.Permittivity epsilon=U.epsilon_0
        "Permittivity of the dielectric"
        annotation (Dialog(__Dymola_label="<html>&epsilon;</html>"));

      final parameter Q.Capacitance C=epsilon*A/L "Capacitance";
      replaceable package Data = Characteristics.'e-'.Graphite constrainedby
        Characteristics.BaseClasses.Characteristic "Material properties"
        annotation (
        Evaluate=true,
        choicesAllMatching=true,
        __Dymola_choicesFromPackage=true);

      parameter Boolean setVelocity=true
        "<html>Material exits at the velocity of the <code>direct</code> connector</html>"
        annotation (
        Evaluate=true,
        Dialog(tab="Assumptions", compact=true),
        choices(__Dymola_checkBox=true));
      parameter Boolean inclVolume=true
        "<html>Include the <code>amagat</code> connector to occupy volume</html>"
        annotation (
        Evaluate=true,
        Dialog(tab="Assumptions", compact=true),
        choices(__Dymola_checkBox=true));

      // Aliases
      Q.Potential w(stateSelect=StateSelect.always, start=0)
        "Electrical potential";
      Q.Current I "Material current";

      output Q.Amount Z(stateSelect=StateSelect.never) = C*w if environment.analysis
        "Amount of charge shifted in the positive direction";

      Connectors.Chemical negative(final n_trans=n_trans)
        "Chemical connector on the 1st side" annotation (Placement(
            transformation(extent={{-70,-10},{-50,10}}), iconTransformation(
              extent={{-70,-10},{-50,10}})));
      Connectors.Chemical positive(final n_trans=n_trans)
        "Chemical connector on the 2nd side" annotation (Placement(
            transformation(extent={{50,-10},{70,10}}), iconTransformation(
              extent={{50,-10},{70,10}})));
      Connectors.Intra intra(final n_trans=n_trans)
        "Translational and thermal interface with the substrate" annotation (
          Placement(transformation(extent={{-10,-50},{10,-30}}),
            iconTransformation(extent={{-10,-50},{10,-30}})));

      Connectors.Amagat amagat(final V=-A*L) if inclVolume
        "Connector for additivity of volume"
        annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

    protected
      outer Conditions.Environment environment "Environmental conditions";

    equation
      // Aliases
      Data.z*w = positive.g - negative.g;
      I = positive.Ndot;

      // Streams
      if setVelocity then
        negative.phi = intra.phi;
        positive.phi = intra.phi;
      else
        negative.phi = inStream(positive.phi);
        positive.phi = inStream(negative.phi);
      end if;
      negative.sT = inStream(positive.sT);
      positive.sT = inStream(negative.sT);

      // Conservation
      0 = negative.Ndot + positive.Ndot "Material (no storage)";
      zeros(n_trans) = Data.m*(actualStream(negative.phi) - actualStream(
        positive.phi))*I + intra.mPhidot "Translational momentum (no storage)";
      der(C*w)/U.s = Data.z*I
        "Electrical energy (reversible; simplified using material conservation and divided by potential)";
      0 = intra.Qdot + (actualStream(negative.phi)*actualStream(negative.phi)
         - actualStream(positive.phi)*actualStream(positive.phi))*I*Data.m/2 +
        intra.phi*intra.mPhidot "Mechanical and thermal energy (no storage)";

      annotation (
        Documentation(info="<html><p>The capacitance (<i>C</i>) is calculated from the surface area (<i>A</i>), length of the gap (<i>L</i>), and the permittivity (&epsilon;) assuming that the 
  charges are uniformly distributed over (infinite) parallel planes.</p>
  
  <p>If <code>setVelocity</code> is <code>true</code>, then the material exits with the
  velocity of the <code>direct</code> connector.  Typically, that connector should be connected to the stationary solid,
  in which case heat will be generated if material arrives with a nonzero velocity.  That heat is rejected to the same connector.
  </p>
  
  </html>"),
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
            true,Dialog(__Dymola_label="<html><i>n</i><sub>trans</sub></html>"));

      parameter Integer z=-1 "Charge number";
      parameter Q.Potential E_A=0 "Activation energy" annotation (Dialog(group=
              "Chemical parameters", __Dymola_label=
              "<html><i>E</i><sub>A</sub></html>"));
      parameter Q.NumberAbsolute alpha(max=1) = 0.5
        "Charge transfer coefficient" annotation (Dialog(group=
              "Chemical parameters", __Dymola_label="<html>&alpha;</html>"));
      parameter Q.Current I0=U.A "Exchange current @ 300 K"
        annotation (Dialog(__Dymola_label="<html><i>I</i><sup>o</sup></html>"));
      parameter Boolean fromI=true
        "<html>Invert the Butler-Volmer equation, if &alpha;=&frac12;</html>"
        annotation (Dialog(tab="Advanced", compact=true), choices(
            __Dymola_checkBox=true));

      Connectors.Chemical negative(final n_trans=n_trans)
        "Chemical connector on the 1st side" annotation (Placement(
            transformation(extent={{-70,-10},{-50,10}}), iconTransformation(
              extent={{-70,-10},{-50,10}})));
      Connectors.Chemical positive(final n_trans=n_trans)
        "Chemical connector on the 2nd side" annotation (Placement(
            transformation(extent={{50,-10},{70,10}}), iconTransformation(
              extent={{50,-10},{70,10}})));

      // Aliases
      Q.TemperatureAbsolute T(start=300*U.K) "Reaction rate";
      Q.Current I(start=0) "Reaction rate";
      Q.Potential Deltag(start=0) "Potential difference";

      Connectors.Intra intra(final n_trans=n_trans)
        "Translational and thermal interface with the substrate" annotation (
          Placement(transformation(extent={{-10,-50},{10,-30}}),
            iconTransformation(extent={{-10,-50},{10,-30}})));

    equation
      // Aliases
      I = positive.Ndot;
      Deltag = positive.g - negative.g;
      T = intra.T;

      // Streams
      negative.phi = inStream(positive.phi);
      positive.phi = inStream(negative.phi);
      negative.sT = inStream(positive.sT);
      positive.sT = inStream(negative.sT);

      // Reaction rate
      if abs(alpha - 0.5) < Modelica.Constants.eps and fromI then
        Deltag = 2*T*asinh(0.5*exp(E_A*(1/T - 1/(300*U.K)))*I/I0);
      else
        I*exp(E_A*(1/T - 1/(300*U.K))) = I0*(exp((1 - alpha)*Deltag/T) - exp(-
          alpha*Deltag/T));
      end if;

      // Conservation (without storage)
      0 = negative.Ndot + positive.Ndot "Material";
      zeros(n_trans) = intra.mPhidot "Translational momentum";
      0 = Deltag*I + intra.Qdot "Energy";
      // Note:  Energy and momentum cancel among the stream terms.

      annotation (
        defaultComponentName="'e-Transfer'",
        Documentation(info="<html><p><code>fromI</code> may help to eliminate nonlinear systems of equations if the 
    <a href=\"modelica://FCSys.Chemistry.Electrochemistry.DoubleLayer\">double layer capacitance</a> is not included.</p></html>"),

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
  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
            {100,100}}), graphics));
end Chemistry;

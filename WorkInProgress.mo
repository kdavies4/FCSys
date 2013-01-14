within FCSys;
package WorkInProgress "Incomplete classes under development"
  extends Modelica.Icons.Package;

  model BCsExamplesClosedVolume
    "<html>Test the <code>ClosedVolume</code> model</html>"

    extends Modelica.Icons.Example;
    extends Modelica.Icons.UnderConstruction;

    FCSys.WorkInProgress.ClosedVolume volume(
      use_portsData=false,
      redeclare Modelica.Media.IdealGases.SingleGases.H2 Medium,
      V=1e-6)
      annotation (Placement(transformation(extent={{-40,-10},{-20,10}})));
    inner Modelica.Fluid.System system
      annotation (Placement(transformation(extent={{40,70},{60,90}})));
    FCSys.Subregions.Subregion subregion(
      L={1,1,1}*U.cm,
      inclReact=false,
      inclYFaces=false,
      inclZFaces=false,
      gas(
        inclH2=true,
        inclH2O=false,
        H2(
          xNegative(isobaric=false),
          xPositive(final isobaric=true),
          setVelX=true)))
      annotation (Placement(transformation(extent={{20,-10},{40,10}})));

    inner BCs.Environment environment(analysis=true, T=293.15*U.K)
      annotation (Placement(transformation(extent={{70,70},{90,90}})));

  protected
    FCSys.Connectors.FaceBusInternal xNegative "Positive face along the x axis"
      annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
  equation

    connect(xNegative.H2.normal, volume.face.normal) annotation (Line(
        points={{5.55112e-16,5.55112e-16},{-6,5.55112e-16},{-6,0},{-10,0},{-10,
            6.10623e-16},{-20,6.10623e-16}},
        color={127,127,127},
        smooth=Smooth.None,
        thickness=0.5));
    connect(xNegative.H2.thermal, volume.face.thermal) annotation (Line(
        points={{5.55112e-16,5.55112e-16},{-6,5.55112e-16},{-6,0},{-10,0},{-10,
            6.10623e-16},{-20,6.10623e-16}},
        color={127,127,127},
        smooth=Smooth.None,
        thickness=0.5));
    connect(xNegative, subregion.xNegative) annotation (Line(
        points={{5.55112e-16,5.55112e-16},{6,5.55112e-16},{6,0},{10,0},{10,
            6.10623e-16},{20,6.10623e-16}},
        color={127,127,127},
        thickness=0.5,
        smooth=Smooth.None));
    annotation (
      experiment(StopTime=10),
      experimentSetupOutput,
      Commands(file="resources/scripts/Dymola/BCs.Examples.FluidAdapt.mos"));
  end BCsExamplesClosedVolume;

  partial model BCsBaseClassesPartialLumpedVessel
    "Lumped volume with a vector of fluid ports and replaceable heat transfer model"
    // Copied and modified from Modelica.Fluid.

    // Added for FCSys:
    import FCSys.Units;

    extends Modelica.Fluid.Interfaces.PartialLumpedVolume;
    // Port definitions
    parameter Integer nPorts=0 "Number of ports" annotation (Evaluate=true,
        Dialog(
        connectorSizing=true,
        tab="General",
        group="Ports"));
    Modelica.Fluid.Vessels.BaseClasses.VesselFluidPorts_b ports[nPorts](
        redeclare each package Medium = Medium) "Fluid inlets and outlets"
      annotation (Placement(transformation(extent={{-40,-10},{40,10}}, origin={
              0,-100})));
    // Port properties
    parameter Boolean use_portsData=true
      "= false to neglect pressure loss and kinetic energy"
      annotation (Evaluate=true,Dialog(tab="General",group="Ports"));
    parameter Modelica.Fluid.Vessels.BaseClasses.VesselPortsData[nPorts]
      portsData if use_portsData "Data of inlet/outlet ports" annotation (
        Dialog(
        tab="General",
        group="Ports",
        enable=use_portsData));
    parameter SI.MassFlowRate m_flow_small(min=0) = system.m_flow_small
      "Regularization range at zero mass flow rate" annotation (Dialog(
        tab="Advanced",
        group="Port properties",
        enable=stiffCharacteristicForEmptyPort));
    /*
  parameter Medium.AbsolutePressure dp_small = system.dp_small
    "Turbulent flow if |dp| >= dp_small (regularization of zero flow)"
    annotation(Dialog(tab="Advanced",group="Ports"));
*/
    Medium.EnthalpyFlowRate ports_H_flow[nPorts];
    Medium.MassFlowRate ports_mXi_flow[nPorts, Medium.nXi];
    Medium.MassFlowRate[Medium.nXi] sum_ports_mXi_flow
      "Substance mass flows through ports";
    Medium.ExtraPropertyFlowRate ports_mC_flow[nPorts, Medium.nC];
    Medium.ExtraPropertyFlowRate[Medium.nC] sum_ports_mC_flow
      "Trace substance mass flows through ports";
    // Heat transfer through boundary
    parameter Boolean use_HeatTransfer=false
      "= true to use the HeatTransfer model"
      annotation (Dialog(tab="Assumptions", group="Heat transfer"));
    replaceable model HeatTransfer =
        Modelica.Fluid.Vessels.BaseClasses.HeatTransfer.IdealHeatTransfer
      constrainedby
      Modelica.Fluid.Vessels.BaseClasses.HeatTransfer.PartialVesselHeatTransfer
      "Wall heat transfer" annotation (Dialog(
        tab="Assumptions",
        group="Heat transfer",
        enable=use_HeatTransfer), choicesAllMatching=true);
    HeatTransfer heatTransfer(
      redeclare final package Medium = Medium,
      final n=1,
      final states={medium.state},
      final use_k=use_HeatTransfer) annotation (Placement(transformation(
          extent={{-10,-10},{30,30}},
          rotation=90,
          origin={-50,-10})));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort if
      use_HeatTransfer
      annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
    // Conservation of kinetic energy
    Medium.Density[nPorts] portDensities
      "densities of the fluid at the device boundary";
    SI.Velocity[nPorts] portVelocities
      "velocities of fluid flow at device boundary";
    SI.EnergyFlowRate[nPorts] ports_E_flow
      "flow of kinetic and potential energy at device boundary";
    // Note:  should use fluidLevel_start - portsData.height
    Real[nPorts] s(each start=fluidLevel_max)
      "curve parameters for port flows vs. port pressures; for further details see, Modelica Tutorial: Ideal switching devices";
    Real[nPorts] ports_penetration
      "penetration of port with fluid, depending on fluid level and port diameter";
    // treatment of pressure losses at ports
    SI.Area[nPorts] portAreas={Modelica.Constants.pi/4*portsData_diameter[i]^2
        for i in 1:nPorts};
    Medium.AbsolutePressure[nPorts] vessel_ps_static
      "static pressures inside the vessel at the height of the corresponding ports, zero flow velocity";

    // Added for FCSys:
    parameter Axis axis=Axis.x "Axis normal to the face";
    FCSys.Connectors.Face face(
      final axis=axis,
      final isobaric=false,
      final inviscidX=true,
      final inviscidY=true,
      final inviscidZ=true)
      "Connection to a face of a FCSys.Subregions.Species model"
      annotation (Placement(transformation(extent={{90,-10},{110,10}})));

  protected
    input SI.Height fluidLevel=0
      "level of fluid in the vessel for treating heights of ports";
    parameter SI.Height fluidLevel_max=1 "maximum level of fluid in the vessel";
    parameter SI.Area vesselArea=Modelica.Constants.inf
      "Area of the vessel used to relate to cross flow area of ports";
    // Treatment of use_portsData=false to neglect portsData and to not require its specification either in this case.
    // Remove portsData conditionally if use_portsData=false. Simplify their use in model equations by always
    // providing portsData_diameter and portsData_height, independent of the use_portsData setting.
    // Note:  this moreover serves as work-around if a tool doesn't support a zero sized portsData record.
    Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_diameter_internal=
        portsData.diameter if use_portsData and nPorts > 0;
    Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_height_internal=
        portsData.height if use_portsData and nPorts > 0;
    Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_zeta_in_internal=
        portsData.zeta_in if use_portsData and nPorts > 0;
    Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_zeta_out_internal=
        portsData.zeta_out if use_portsData and nPorts > 0;
    Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_diameter;
    Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_height;
    Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_zeta_in;
    Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_zeta_out;
  equation
    // Added for FCSys:
    face.normal.rho = (medium.d/medium.MM)*Units.mol/Units.m^3;
    face.thermal.T = medium.T*Units.K;

    mb_flow = sum(ports.m_flow) + medium.MM*face.normal.Ndot/Units.kat
      "Added term for FCSys";
    mbXi_flow = sum_ports_mXi_flow;
    mbC_flow = sum_ports_mC_flow;
    Hb_flow = sum(ports_H_flow) + sum(ports_E_flow) + medium.h*medium.MM*face.normal.Ndot
      /Units.kat "Added term for FCSys";
    Qb_flow = heatTransfer.Q_flows[1] + face.thermal.Qdot/Units.W
      "Added term for FCSys";
    // Only one connection allowed to a port to avoid unwanted ideal mixing
    for i in 1:nPorts loop
      assert(cardinality(ports[i]) <= 1, "
each ports[i] of volume can at most be connected to one component.
If two or more connections are present, ideal mixing takes
place with these connections, which is usually not the intention
of the modeler. Increase nPorts to add an additional port.
");
    end for;
    // Check for correct solution
    assert(fluidLevel <= fluidLevel_max,
      "Vessel is overflowing (fluidLevel > fluidLevel_max = " + String(
      fluidLevel) + ")");
    assert(fluidLevel > -1e-6*fluidLevel_max, "Fluid level (= " + String(
      fluidLevel) + ") is below zero meaning that the solution failed.");
    // Boundary conditions
    // treatment of conditional portsData
    connect(portsData_diameter, portsData_diameter_internal);
    connect(portsData_height, portsData_height_internal);
    connect(portsData_zeta_in, portsData_zeta_in_internal);
    connect(portsData_zeta_out, portsData_zeta_out_internal);
    if not use_portsData then
      portsData_diameter = zeros(nPorts);
      portsData_height = zeros(nPorts);
      portsData_zeta_in = zeros(nPorts);
      portsData_zeta_out = zeros(nPorts);
    end if;
    // actual definition of port variables
    for i in 1:nPorts loop
      if use_portsData then
        // dp = 0.5*zeta*d*v*|v|
        // Note:  assume vessel_ps_static for portDensities to avoid algebraic loops for ports.p
        portDensities[i] = noEvent(Medium.density(Medium.setState_phX(
            vessel_ps_static[i],
            actualStream(ports[i].h_outflow),
            actualStream(ports[i].Xi_outflow))));
        portVelocities[i] = smooth(0, ports[i].m_flow/portAreas[i]/
          portDensities[i]);
        // Note:  the penetration should not go too close to zero as this would prevent a vessel from running empty
        ports_penetration[i] = Modelica.Fluid.Utilities.regStep(
            fluidLevel - portsData_height[i] - 0.1*portsData_diameter[i],
            1,
            1e-3,
            0.1*portsData_diameter[i]);
      else
        // an infinite port diameter is assumed
        portDensities[i] = medium.d;
        portVelocities[i] = 0;
        ports_penetration[i] = 1;
      end if;
      // fluid flow through ports
      if fluidLevel >= portsData_height[i] then
        // regular operation: fluidLevel is above ports[i]
        // Note:  >= covers default values of zero as well
        if use_portsData then
          /* Without regularization
        ports[i].p = vessel_ps_static[i] + 0.5*ports[i].m_flow^2/portAreas[i]^2
                      * noEvent(if ports[i].m_flow>0 then zeta_in[i]/portDensities[i] else -zeta_out[i]/medium.d);
        */
          ports[i].p = vessel_ps_static[i] + (0.5/portAreas[i]^2*
            Modelica.Fluid.Utilities.regSquare2(
              ports[i].m_flow,
              m_flow_small,
              (portsData_zeta_in[i] - 1 + portAreas[i]^2/vesselArea^2)/
              portDensities[i]*ports_penetration[i],
              (portsData_zeta_out[i] + 1 - portAreas[i]^2/vesselArea^2)/medium.d
              /ports_penetration[i]));
          /*
        // alternative formulation m_flow=f(dp); not allowing the ideal portsData_zeta_in[i]=1 though
        ports[i].m_flow = smooth(2, portAreas[i]*Utilities.regRoot2(ports[i].p - vessel_ps_static[i], dp_small,
                                     2*portDensities[i]/portsData_zeta_in[i],
                                     2*medium.d/portsData_zeta_out[i]));
        */
        else
          ports[i].p = vessel_ps_static[i];
        end if;
        s[i] = fluidLevel - portsData_height[i];
      elseif s[i] > 0 or portsData_height[i] >= fluidLevel_max then
        // ports[i] is above fluidLevel and has inflow
        ports[i].p = vessel_ps_static[i];
        s[i] = ports[i].m_flow;
      else
        // ports[i] is above fluidLevel, preventing outflow
        ports[i].m_flow = 0;
        s[i] = (ports[i].p - vessel_ps_static[i])/Medium.p_default*(
          portsData_height[i] - fluidLevel);
      end if;
      ports[i].h_outflow = medium.h;
      ports[i].Xi_outflow = medium.Xi;
      ports[i].C_outflow = C;
      ports_H_flow[i] = ports[i].m_flow*actualStream(ports[i].h_outflow)
        "Enthalpy flow";
      ports_E_flow[i] = ports[i].m_flow*(0.5*portVelocities[i]*portVelocities[i]
         + system.g*portsData_height[i]) "Flow of kinetic and potential energy";
      ports_mXi_flow[i, :] = ports[i].m_flow*actualStream(ports[i].Xi_outflow)
        "Component mass flow";
      ports_mC_flow[i, :] = ports[i].m_flow*actualStream(ports[i].C_outflow)
        "Trace substance mass flow";
    end for;
    for i in 1:Medium.nXi loop
      sum_ports_mXi_flow[i] = sum(ports_mXi_flow[:, i]);
    end for;
    for i in 1:Medium.nC loop
      sum_ports_mC_flow[i] = sum(ports_mC_flow[:, i]);
    end for;
    connect(heatPort, heatTransfer.heatPorts[1]) annotation (Line(
        points={{-100,5.55112e-16},{-87,5.55112e-16},{-87,2.22045e-15},{-74,
            2.22045e-15}},
        color={191,0,0},
        smooth=Smooth.None));

    annotation (
      Documentation(info="<html>
<p>
This base class extends PartialLumpedVolume with a vector of fluid ports and a replaceable wall HeatTransfer model.
<p>
The following modeling assumption are made:
<ul>
<li>homogeneous medium, i.e., phase separation is not taken into account,</li>
<li>no kinetic energy in the fluid, i.e., kinetic energy dissipates into the internal energy,</li>
<li>pressure loss definitions at vessel ports assume incompressible fluid,</li>
<li>outflow of ambient media is prevented at each port assuming check valve behavior.
    If <code> fluidlevel &lt; portsData_height[i] </code>and &nbsp; <code> ports[i].p &lt; vessel_ps_static[i]</code> massflow at the port is set to 0.</li>
</ul>
</p>
Each port has a (hydraulic) diameter and a height above the bottom of the vessel, which can be configured using the &nbsp;<b><code>portsData</code></b> record.
Alternatively the impact of port geometries can be neglected with <code>use_portsData=false</code>. This might be useful for early
design studies. Note that this means to assume an infinite port diameter at the bottom of the vessel.
Pressure drops and heights of the ports as well as kinetic and potential energy fluid entering or leaving the vessel are neglected then.
<p>
The following variables need to be defined by an extending model:
<ul>
<li><code>input fluidVolume</code>, the volume of the fluid in the vessel,</li>
<li><code>vessel_ps_static[nPorts]</code>, the static pressures inside the vessel at the height of the corresponding ports, at zero flow velocity, and</li>
<li><code>Wb_flow</code>, work term of the energy balance, e.g., p*der(V) if the volume is not constant or stirrer power.</li>
</ul>
An extending model should define:
<ul>
<li><code>parameter vesselArea</code> (default: Modelica.Constants.inf m2), the area of the vessel, to be related to cross flow areas of the ports for the consideration of dynamic pressure effects.</li>
</ul>
Optionally the fluid level may vary in the vessel, which effects the flow through the ports at configurable <code>portsData_height[nPorts]</code>.
This is why an extending model with varying fluid level needs to define:
<ul>
<li><code>input fluidLevel (default: 0m)</code>, the level the fluid in the vessel, and</li>
<li><code>parameter fluidLevel_max (default: 1m)</code>, the maximum level that must not be exceeded. Ports at or above fluidLevel_max can only receive inflow.</li>
</ul>
An extending model should not access the <code>portsData</code> record defined in the configuration dialog,
as an access to <code>portsData</code> may fail for <code>use_portsData=false</code> or <code>nPorts=0</code>.
Instead the predefined variables
<ul>
<li><code>portsData_diameter[nPorts]</code></li>,
<li><code>portsData_height[nPorts]</code></li>,
<li><code>portsData_zeta_in[nPorts]</code></li>, and
<li><code>portsData_zeta_out[nPorts]</code></li>
</ul>
should be used if these values are needed.
</p>
</html>", revisions="<html>
<ul>
<li><i>Jan. 2009</i> by R&uuml;diger Franke: extended with
   <ul><li>portsData record and threat configurable port heights,</li>
       <li>consideration of kinetic and potential energy of fluid entering or leaving in energy balance</li>
   </ul>
</li>
<li><i>Dec. 2008</i> by R&uuml;diger Franke: derived from OpenTank, in order to make general use of configurable port diameters</i>
</ul>
</html>"),
      Diagram(coordinateSystem(preserveAspectRatio=true,extent={{-100,-100},{
              100,100}}), graphics),
      Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
              100}}), graphics={Text(
              extent={{-150,110},{150,150}},
              textString="%name",
              lineColor={0,0,255})}));
  end BCsBaseClassesPartialLumpedVessel;

  model ClosedVolume
    "Volume of fixed size, closed to the ambient, with inlet/outlet ports"
    // Copied from Modelica.Fluid.

    import Modelica.Constants.pi;

    // Mass and energy balance, ports
    extends FCSys.WorkInProgress.BCsBaseClassesPartialLumpedVessel(
      final fluidVolume=V,
      vesselArea=pi*(3/4*V)^(2/3),
      heatTransfer(surfaceAreas={4*pi*(3/4*V/pi)^(2/3)}));
    parameter SI.Volume V "Volume";
  equation
    Wb_flow = 0;
    for i in 1:nPorts loop
      vessel_ps_static[i] = medium.p;
    end for;
    annotation (
      defaultComponentName="volume",
      Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
              100}}), graphics={Ellipse(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,0},
              fillPattern=FillPattern.Sphere,
              fillColor={170,213,255}),Text(
              extent={{-150,12},{150,-18}},
              lineColor={0,0,0},
              textString="V=%V")}),
      Documentation(info="<html>
<p>
Ideally mixed volume of constant size with two fluid ports and one medium model.
The flow properties are computed from the upstream quantities, pressures are equal in both nodes and the medium model if <code>use_portsData=false</code>.
Heat transfer through a thermal port is possible, it equals zero if the port remains unconnected.
A spherical shape is assumed for the heat transfer area, with V=4/3*pi*r^3, A=4*pi*r^2.
Ideal heat transfer is assumed per default; the thermal port temperature is equal to the medium temperature.
</p>
<p>
If <code>use_portsData=true</code>, the port pressures represent the pressures just after the outlet (or just before the inlet) in the attached pipe.
The hydraulic resistances <code>portsData.zeta_in</code> and <code>portsData.zeta_out</code> determine the dissipative pressure drop between volume and port depending on
the direction of mass flow. See <a href=\"modelica://Modelica.Fluid.Vessels.BaseClasses.VesselPortsData\">VesselPortsData</a> and <i>[Idelchik, Handbook of Hydraulic Resistance, 2004]</i>.
</p>
</html>"),
      Diagram(coordinateSystem(preserveAspectRatio=true,extent={{-100,-100},{
              100,100}}), graphics));
  end ClosedVolume;

  model SubRegionsSpeciesSpecies0Amount "Species with zero particle number"
    extends Subregions.Species.Species(
      final overrideEOS=true,
      final rho_IC=0,
      final derrho_IC=0,
      N(stateSelect=StateSelect.never),
      phi(each stateSelect=StateSelect.never),
      T(stateSelect=StateSelect.never));
    // Note:  StateSelect.never is necessary to avoid dynamic state selection
    // in Dymola 7.4.
    annotation (Documentation(info="<html>
  <p>For more information, see the <a href=\"modelica://FCSys.Subregions.Species.Species\">Species</a> model.</p></html>"));
  end SubRegionsSpeciesSpecies0Amount;

  model BCsAdaptersSpeciesFluid
    "<html>Adapter to connect a single fluid species between <a href=\"modelica://FCSys\">FCSys</a> and <a href=\"modelica://Modelica\">Modelica</a></html>"

    extends BCs.Adapters.Species.FluidNonionic;
    extends Modelica.Icons.UnderConstruction;

    Modelica.Electrical.Analog.Interfaces.NegativePin pin if Data.z <> 0
      "Modelica electrical pin" annotation (Placement(transformation(extent={{
              70,30},{90,50}}), iconTransformation(extent={{70,30},{90,50}})));

  equation
    // **Add electrical equations.
    annotation (
      Documentation(info="<html><p>The electrical connector (<code>pin</code>) is only included
    if the species is ionic.
    </p>
    <p>For additional information, see the
    <a href=\"modelica://FCSys.BCs.Adapters.Species.BaseClasses.PartialSpecies\">
    PartialSpecies</a> model.</p>
    </html>"),
      Icon(graphics={Line(
              points={{0,40},{80,40}},
              color={0,0,255},
              smooth=Smooth.None),Line(
              points={{0,60},{0,20}},
              color={0,0,0},
              smooth=Smooth.None,
              pattern=LinePattern.Dash)}),
      Diagram(graphics));
  end BCsAdaptersSpeciesFluid;

  model BCsAdaptersPhasesIonomer
    "<html>Adapter for ionomer between <a href=\"modelica://FCSys\">FCSys</a> and <a href=\"modelica://Modelica\">Modelica</a></html>"
    extends BCs.Adapters.Phases.BaseClasses.PartialPhase;
    extends Modelica.Icons.UnderConstruction;

    BCs.Adapters.Species.Solid C19HF37O5S(redeclare package Data =
          FCSys.Characteristics.C19HF37O5S.Ionomer)
      annotation (Placement(transformation(extent={{-10,10},{10,30}})));
    FCSys.WorkInProgress.BCsAdaptersSpeciesFluid 'H+'(redeclare package Data =
          FCSys.Characteristics.'H+'.Ionomer, redeclare package Medium =
          Modelica.Media.IdealGases.SingleGases.H2)
      annotation (Placement(transformation(extent={{-10,-30},{10,-10}})));
    // **Use model for H instead.

    BCs.Adapters.Species.FluidNonionic H2O(redeclare package Data =
          FCSys.Characteristics.H2O.Ionomer, redeclare package Medium =
          Modelica.Media.IdealGases.SingleGases.H2O)
      annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));
    BCs.Adapters.Junctions.Junction2 junction2
      annotation (Placement(transformation(extent={{60,-50},{40,-30}})));
    Modelica.Fluid.Interfaces.FluidPort_b fluidPort(redeclare final package
        Medium = Medium) "Modelica fluid port" annotation (Placement(
          transformation(extent={{70,-50},{90,-30}}), iconTransformation(extent
            ={{70,-50},{90,-30}})));
    Modelica.Electrical.Analog.Interfaces.NegativePin pin
      "Modelica electrical pin" annotation (Placement(transformation(extent={{
              70,30},{90,50}}), iconTransformation(extent={{70,30},{90,50}})));
  equation
    connect(C19HF37O5S.face.thermal, face.C19HF37O5S.thermal) annotation (Line(
        points={{-8,20},{-40,20},{-40,5.55112e-16},{-80,5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('H+'.face.normal, face.'H+'.normal) annotation (Line(
        points={{-8,-20},{-40,-20},{-40,5.55112e-16},{-80,5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('H+'.face.thermal, face.'H+'.thermal) annotation (Line(
        points={{-8,-20},{-40,-20},{-40,5.55112e-16},{-80,5.55112e-16}},
        color={127,127,127},
        smooth=Smooth.None));

    connect('H+'.pin, pin) annotation (Line(
        points={{8,-16},{60,-16},{60,40},{80,40}},
        color={0,0,255},
        smooth=Smooth.None));

    connect('H+'.heatPort, heatPort) annotation (Line(
        points={{8,-20},{40,-20},{40,5.55112e-16},{80,5.55112e-16}},
        color={191,0,0},
        smooth=Smooth.None));
    connect(C19HF37O5S.heatPort, heatPort) annotation (Line(
        points={{8,20},{40,20},{40,5.55112e-16},{80,5.55112e-16}},
        color={191,0,0},
        smooth=Smooth.None));
    connect(junction2.mixturePort, fluidPort) annotation (Line(
        points={{58,-40},{80,-40}},
        color={0,127,255},
        smooth=Smooth.None));
    annotation (Icon(graphics={Line(
              points={{0,60},{0,-60}},
              color={0,0,0},
              smooth=Smooth.None,
              pattern=LinePattern.Dash,
              thickness=0.5),Line(
              points={{0,0},{-80,0}},
              color={127,127,127},
              smooth=Smooth.None,
              thickness=0.5),Line(
              points={{0,40},{80,40}},
              color={0,0,255},
              smooth=Smooth.None),Line(
              points={{0,0},{80,0}},
              color={191,0,0},
              smooth=Smooth.None),Line(
              points={{0,-40},{80,-40}},
              color={0,127,255},
              smooth=Smooth.None)}));
  end BCsAdaptersPhasesIonomer;
end WorkInProgress;

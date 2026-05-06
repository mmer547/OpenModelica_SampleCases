package OilCircuit
  package Interfaces
    connector OilPort
      Real p "Pressure [Pa]";
      flow Real m_flow "Mass flow rate [kg/s]";
    end OilPort;

    connector PortIn
      extends OilPort;
      annotation(
        Icon(graphics = {Rectangle(fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}})}));
    end PortIn;

    connector PortOut
      extends OilPort;
      annotation(
        Icon(graphics = {Rectangle(fillColor = {255, 0, 0}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}})}));
    end PortOut;
  end Interfaces;

  package Components
    model PressureBoundary
      outer SystemParameters system;
      Interfaces.PortIn port annotation(
        Placement(transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {0, 100}, extent = {{-10, -10}, {10, 10}})));
    equation
      port.p = system.p_ref;
      annotation(
        Icon(graphics = {Rectangle(fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(extent = {{-80, 100}, {80, -100}}, textString = "Pressure
    Boundary")}));
    end PressureBoundary;

    model PositiveDisplacementPump
      outer SystemParameters system;
      parameter Real Vdisp = 1.0e-6 "Displacement volume per revolution [m3/rev]";
      parameter Real n = 50 "Rotational speed [rev/s]";
      parameter Real eta_v = 0.8 "Volumetric efficiency [-]";
      Real Q "Volume flow rate [m3/s]";
      Real m_flow "Mass flow rate [kg/s]";
      Interfaces.PortIn inlet annotation(
        Placement(transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}})));
      Interfaces.PortOut outlet annotation(
        Placement(transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}})));
    equation
      Q = Vdisp*n*eta_v;
      m_flow = system.rho*Q;
      inlet.m_flow = m_flow;
      outlet.m_flow = -m_flow;
      annotation(
        Icon(graphics = {Rectangle(fillColor = {255, 255, 255}, extent = {{-100, 100}, {100, -100}}), Text(extent = {{-80, 100}, {80, -100}}, textString = "textPositive
    Displacement
    Pump")}));
    end PositiveDisplacementPump;

    model LinearResistance
      parameter Real R = 1.0e10 "Hydraulic resistance [Pa.s/m3]";
      outer SystemParameters system;
      Interfaces.PortIn a annotation(
        Placement(transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}})));
      Interfaces.PortOut b annotation(
        Placement(transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}})));
      Real Q "Volume flow rate from a to b [m3/s]";
    equation
      Q = (a.p - b.p)/R;
      a.m_flow = system.rho*Q;
      b.m_flow = -system.rho*Q;
      annotation(
        Icon(graphics = {Rectangle(fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(extent = {{-80, 100}, {80, -100}}, textString = "Linear
    Resistance")}));
    end LinearResistance;

    model Orifice
      outer SystemParameters system;
      parameter Real Cd = 0.7 "Discharge coefficient [-]";
      parameter Real A = 1.0e-6 "Flow area [m2]";
      Interfaces.PortIn a annotation(
        Placement(transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}})));
      Interfaces.PortOut b annotation(
        Placement(transformation(origin = {100, -2}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {100, -2}, extent = {{-10, -10}, {10, 10}})));
      Real dp;
      Real Q;
    equation
      dp = a.p - b.p;
      Q = Cd*A*sqrt(2/system.rho)*dp/sqrt(sqrt(dp*dp + system.dp_eps*system.dp_eps));
      a.m_flow = system.rho*Q;
      b.m_flow = -system.rho*Q;
      annotation(
        Icon(graphics = {Rectangle(fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(extent = {{-80, 100}, {80, -100}}, textString = "Orifice")}));
    end Orifice;

    model OilVolume
      outer SystemParameters system;
      parameter Real V = 1.0e-5 "Volume [m3]";
      Real p(start = system.p_start, fixed = true) "Pressure [Pa]";
      Interfaces.PortIn port annotation(
        Placement(transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {0, 100}, extent = {{-10, -10}, {10, 10}})));
    equation
      port.p = p;
      der(p) = system.beta/(system.rho*V)*port.m_flow;
      annotation(
        Icon(graphics = {Rectangle(fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(extent = {{-80, 100}, {80, -100}}, textString = "OilVolume")}));
    end OilVolume;

    model DrainChamber
      outer SystemParameters system;
      parameter Real V = 5e-5 "Volume [m3]";
      parameter Real D_inner = 0.02;
      parameter Real L_chamber = 0.03;
      parameter Real e = 0.5e-3;
      parameter Real omega = 2*Modelica.Constants.pi*50;
      parameter Real k_motion = 0.3;
      Real p(start = system.p_start, fixed = true) "Pressure [Pa]";
      Real Q_motion "Equivalent volume flow by orbital motion [m3/s]";
      Real A_eff;
      Real motionRamp;
      Interfaces.PortIn inlet annotation(
        Placement(transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}})));
      Interfaces.PortOut outlet annotation(
        Placement(transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}})));
    equation
      inlet.p = p;
      outlet.p = p;
      A_eff = k_motion*D_inner*L_chamber;
      motionRamp = if time < system.t_ramp then time/system.t_ramp else 1.0;
      Q_motion = motionRamp*A_eff*e*omega*sin(omega*time);
      der(p) = system.beta/(system.rho*V)*(inlet.m_flow + outlet.m_flow + system.rho*Q_motion);
      annotation(
        Icon(graphics = {Rectangle(fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(extent = {{-80, 100}, {80, -100}}, textString = "Drain
Chamber")}));
    end DrainChamber;
  end Components;

  package Examples
  end Examples;

  package Tests
    model test1
      Components.PressureBoundary pressureBoundary annotation(
        Placement(transformation(origin = {-86, -76}, extent = {{-10, -10}, {10, 10}})));
      inner SystemParameters system annotation(
        Placement(transformation(origin = {-92, 92}, extent = {{-10, -10}, {10, 10}})));
      Components.PositiveDisplacementPump positiveDisplacementPump annotation(
        Placement(transformation(origin = {-66, 12}, extent = {{-10, -10}, {10, 10}})));
      Components.LinearResistance linearResistance annotation(
        Placement(transformation(origin = {-32, 12}, extent = {{-10, -10}, {10, 10}})));
      Components.Orifice orifice annotation(
        Placement(transformation(origin = {0, 56}, extent = {{-10, -10}, {10, 10}})));
      Components.Orifice orifice1 annotation(
        Placement(transformation(origin = {0, 12}, extent = {{-10, -10}, {10, 10}})));
      Components.Orifice orifice2 annotation(
        Placement(transformation(origin = {0, -36}, extent = {{-10, -10}, {10, 10}})));
  Components.OilVolume oilVolume annotation(
        Placement(transformation(origin = {44, -10}, extent = {{-10, -10}, {10, 10}})));
  Components.LinearResistance linearResistance1 annotation(
        Placement(transformation(origin = {80, 8}, extent = {{-10, -10}, {10, 10}})));
    equation
      connect(linearResistance.b, orifice2.a) annotation(
        Line(points = {{-22, 12}, {-16, 12}, {-16, -36}, {-10, -36}}));
  connect(linearResistance.b, orifice1.a) annotation(
        Line(points = {{-22, 12}, {-10, 12}}));
  connect(linearResistance.b, orifice.a) annotation(
        Line(points = {{-22, 12}, {-18, 12}, {-18, 56}, {-10, 56}}));
  connect(positiveDisplacementPump.outlet, linearResistance.a) annotation(
        Line(points = {{-56, 12}, {-42, 12}}));
  connect(positiveDisplacementPump.inlet, pressureBoundary.port) annotation(
        Line(points = {{-76, 12}, {-86, 12}, {-86, -66}}));
  connect(orifice2.b, oilVolume.port) annotation(
        Line(points = {{10, -36}, {22, -36}, {22, 8}, {44, 8}, {44, 0}}));
  connect(orifice1.b, oilVolume.port) annotation(
        Line(points = {{10, 12}, {44, 12}, {44, 0}}));
  connect(orifice.b, oilVolume.port) annotation(
        Line(points = {{10, 56}, {44, 56}, {44, 0}}));
  connect(linearResistance1.b, pressureBoundary.port) annotation(
        Line(points = {{90, 8}, {96, 8}, {96, -58}, {-86, -58}, {-86, -66}}));
  connect(linearResistance1.a, oilVolume.port) annotation(
        Line(points = {{70, 8}, {44, 8}, {44, 0}}));
    end test1;
  end Tests;

  model SystemParameters
    parameter Real rho = 850 "Density [kg/m3]";
    parameter Real mu = 0.03 "Dynamic viscosity [Pa.s]";
    parameter Real beta = 1.5e9 "Bulk modulus [Pa]";
    parameter Real p_ref = 4.0e6 "Reference pressure [Pa]";
    parameter Real p_start = p_ref "Initial pressure [Pa]";
    parameter Real dp_eps = 100.0 "Regularization pressure [Pa]";
    annotation(
      defaultComponentName = "system",
      defaultComponentPrefixes = "inner");
  end SystemParameters;
end OilCircuit;

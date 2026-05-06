package Compression_Chamber
  package Components
    model CompressionChamber
      outer SystemParameters system;
      replaceable package Medium = Media.R410A;
      
      Interfaces.port_in suctionPort annotation(
        Placement(transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {-100, 80}, extent = {{-10, -10}, {10, 10}})));
      Interfaces.port_leak leakPort annotation(
        Placement(transformation(origin = {-100, -60}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {-100, -80}, extent = {{-10, -10}, {10, 10}})));
      Interfaces.port_out dischargePort annotation(
        Placement(transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}})));
    
      output Real Vcalc "下限処理前の計算上の体積 [m3]";
      output Real V "圧縮室内体積 [m3]";
      output Real d "圧縮室内ガス密度 [kg/m3]";
      output Real p "圧縮室内圧力 [Pa]";
      output Real u "圧縮室内ガスの比内部エネルギー [J/kg]";
      output Real h "圧縮室内ガスの比エンタルピー [J/kg]";
      output Real cv "圧縮室内ガスの定容比熱 [J/kg/K]";
    
      Real T(start = system.Ts, fixed = true, stateSelect = StateSelect.always) "圧縮室内ガス温度 [K]";
      Real m(start = system.m0, fixed = true, stateSelect = StateSelect.always) "圧縮室内ガス質量 [kg]";
    
      output Real m_flow_net "圧縮室への正味質量流入量 [kg/s]";
      output Real angle "周期内角度 [deg]";
      output Real timeInCycle "周期内時刻 [s]";
    
    protected
      constant Real pi = Modelica.Constants.pi;
      Medium.ThermodynamicState state;
    
    equation
  /*
        圧縮室体積モデル
    
        圧縮室体積は、時間に対して周期的に変化する既知関数として与える。
        Vcalc は下限処理前の計算上の体積である。
        Vcalc が system.Vmin 以下になった場合は、数値的な破綻を避けるため、
        実際に用いる体積 V を system.Vmin に制限する。
      */
      angle = time/system.period*360;
      timeInCycle = mod(time, system.period);
      
      Vcalc = 
        system.Vcc/2 
        + 0.5*system.Vcc*cos(pi*system.f*timeInCycle);
      if Vcalc <= system.Vmin then
        V = system.Vmin;
      else
        V = Vcalc;
      end if;
// 圧縮室状態
      d = m/V;
      state = Medium.setState_dT(d, T);
      p = Medium.pressure(state);
      u = Medium.specificInternalEnergy(state);
      h = Medium.specificEnthalpy(state);
      cv = Medium.specificHeatCapacityCv(state);
  /*
        各ポートに与える状態量
    
        各ポートの圧力は圧縮室圧力 p と等しいとする。
    
        漏れポートおよび吐出ポートから流出するガスは、
        圧縮室内ガスと同じ密度・比エンタルピーを持つと仮定する。
    
        主吸入ポートの密度・比エンタルピーは吸入側モデルから与える。
        そのため、このモデルでは suctionPort.d, suctionPort.h は指定しない。
      */
      suctionPort.p = p;
      
      leakPort.p = p;
      leakPort.d = d;
      leakPort.h = h;
    
      dischargePort.p = p;
      dischargePort.d = d;
      dischargePort.h = h;

      /*
        質量保存式
    
        本モデルでは、各ポートの m_flow は
        「圧縮室へ流入する向きを正」と定義する。
    
          m_flow > 0 : 圧縮室への流入
          m_flow < 0 : 圧縮室からの流出
    
        したがって、各ポートの質量流量の和が
        圧縮室内質量 m の時間変化率になる。
      */
      m_flow_net = 
        suctionPort.m_flow 
        + leakPort.m_flow 
        + dischargePort.m_flow;
      der(m) = m_flow_net;

      /*
        エネルギー保存式
    
        圧縮室内ガスの内部エネルギー変化を、
        圧力仕事と各ポートを通じたエンタルピー流入出で表す。
    
        -p*der(V) は圧縮・膨張に伴う圧力仕事を表す。
        der(V) < 0 のとき圧縮であり、この項は正になる。
    
        主吸入では、吸入側から流入するガスの比エンタルピー suctionPort.h を用いる。
        漏れおよび吐出では、圧縮室から流出するガスが
        圧縮室内ガスと同じ比エンタルピー h を持つと仮定する。
    
        現在想定している流れの向きは以下の通り。
    
          suctionPort.m_flow > 0    : 吸入側から圧縮室へ流入
          leakPort.m_flow < 0       : 圧縮室から吸入側へ漏れ
          dischargePort.m_flow < 0  : 圧縮室から吐出側へ流出
      */
      m*cv*der(T) = 
        -p*der(V) 
        + suctionPort.m_flow*(suctionPort.h - u) 
        + leakPort.m_flow*(h - u) 
        + dischargePort.m_flow*(h - u);
        
      /*
        周期リセット
    
        各周期の開始時に、圧縮室温度 T を吸入温度 system.Ts に戻す。
        これは、周期運転を簡略的に表現するためのモデル上の仮定である。
    
        必要に応じて、圧縮室内質量 m も初期値 system.m0 に戻すことができる。
        その場合は reinit(m, system.m0) を有効にする。
      */
      when sample(system.period, system.period) then
        reinit(T, system.Ts);

      end when;
      annotation(
        Icon(graphics = {Rectangle(fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(extent = {{-80, 100}, {80, -100}}, textString = "Compression
Chamber")}));
    end CompressionChamber;

    model DischargeValve
      outer SystemParameters system;
      Interfaces.port_in inlet annotation(
        Placement(transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}})));
      Interfaces.port_out outlet annotation(
        Placement(transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}})));
      output Real m_flow_out "吐出質量流量 [kg/s]";
      output Integer m_flow_out_flg(start = 0) "吐出弁開閉フラグ";
    equation
      /*
        吐出弁モデル
    
        吐出弁は、入口側圧力 inlet.p が出口側圧力 outlet.p より高い場合に開く
        一方向弁として扱う。
    
        inlet.p > outlet.p のとき、入口から出口へ流れる。
        inlet.p <= outlet.p のとき、弁は閉じており、流量は 0 とする。
    
        ここでは弁体の運動、弁ばね、慣性、流路面積の時間変化は考慮せず、
        圧力差に基づく準定常のオリフィス流れとして簡略化している。
      */
      if inlet.p > outlet.p then
        m_flow_out = system.CdOut*system.Aout*sqrt(2*inlet.d*abs(inlet.p - outlet.p));
        m_flow_out_flg = 1;
      else
        m_flow_out = 0;
        m_flow_out_flg = 0;
      end if;
      /*
        質量流量の符号規約
    
        m_flow_out は、弁を通過して inlet から outlet へ流れる向きを正とする。
    
        一方、各ポートの m_flow は、
        そのポートが接続されているコンポーネントへ流入する向きを正とする。
    
        そのため、入口ポート inlet では、流れは弁モデルへ流入する向きなので正、
        出口ポート outlet では、流れは弁モデルから流出する向きなので負になる。
      */  
      inlet.m_flow = m_flow_out;
      outlet.m_flow = -m_flow_out;
      /*
        流出側へ渡す物性
    
        吐出弁を通過するガスは、入口側と同じ密度・比エンタルピーを持つと仮定する。
        したがって、outlet.d および outlet.h には inlet 側の値をそのまま与える。
      */
      outlet.d = inlet.d;
      outlet.h = inlet.h;
      annotation(
        Icon(graphics = {Rectangle(fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(extent = {{-80, 100}, {80, -100}}, textString = "Discharge
    Valve")}));
    end DischargeValve;

    model SuctionChamber
      outer SystemParameters system;
      replaceable package Medium = Media.R410A;
      Interfaces.port_out suctionPort annotation(
        Placement(transformation(origin = {100, 80}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {100, 80}, extent = {{-10, -10}, {10, 10}})));
      Interfaces.port_leak leakPort annotation(
        Placement(transformation(origin = {100, -80}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {100, -80}, extent = {{-10, -10}, {10, 10}})));
      output Real ps "吸入側圧力 [Pa]";
      output Real Ts "吸入側温度 [K]";
      output Real ds "吸入側密度 [kg/m3]";
      output Real hs "吸入側比エンタルピー [J/kg]";
      output Real m_flow_suction "主吸入質量流量 [kg/s]";
      output Real m_flow_leak_out "圧縮室から吸入側への漏れ質量流量 [kg/s]";
      output Integer m_flow_leak_out_flg(start = 0) "吸入側への漏れ発生フラグ";
    
      Real timeInCycle "周期内時刻 [s]";
      Real suctionDuration "主吸入が発生する時間幅 [s]";
    protected
      Medium.ThermodynamicState state_suction;
    equation
      /*
        吸入側境界条件
    
        吸入側は、一定圧力 ps、一定温度 Ts の境界として扱う。
        圧力と温度は SystemParameters で定義された値を使用する。
      */
      ps = system.ps;
      Ts = system.Ts;
  /*
        吸入側冷媒状態
    
        吸入側の圧力 ps と温度 Ts から冷媒物性を計算し、
        密度 ds と比エンタルピー hs を求める。
      */
      state_suction = Medium.setState_pT(ps, Ts);
      ds = Medium.density(state_suction);
      hs = Medium.specificEnthalpy(state_suction);
  /*
        主吸入流量モデル
    
        各周期の先頭から suctionAngle [deg] に相当する時間だけ、
        一定の主吸入質量流量 system.Qm が流れると仮定する。
    
        timeInCycle は現在時刻を 1 周期内に折り返した時刻である。
        suctionDuration は suctionAngle に対応する吸入継続時間である。
      */
      timeInCycle = mod(time, system.period);
      suctionDuration = system.suctionAngle/360*system.period;
      if timeInCycle > 0 and timeInCycle < suctionDuration then
        m_flow_suction = system.Qm;
      else
        m_flow_suction = 0;
      end if;
  /*
        主吸入ポートの質量流量
    
        suctionPort.m_flow は、SuctionChamber へ流入する向きを正とする。
        主吸入は吸入側から圧縮室へ流出するため、
        SuctionChamber から見た suctionPort.m_flow は負になる。
      */
      suctionPort.m_flow = -m_flow_suction;
  /*
        主吸入ポートに与える物性
    
        主吸入ポートから圧縮室へ流れるガスは、
        吸入側状態と同じ密度 ds、比エンタルピー hs を持つと仮定する。
    
        圧力 suctionPort.p は圧縮室側で与える。
        そのため、このモデルでは suctionPort.p は指定しない。
      */
      suctionPort.d = ds;
      suctionPort.h = hs;
  /*
        吸入側への漏れ流量モデル
    
        leakPort は CompressionChamber.leakPort と直接接続される。
        接続により、leakPort.p は圧縮室側で与えられた圧力と等しくなる。
    
        圧縮室圧力 leakPort.p が吸入側圧力 ps より高い場合、
        圧縮室から吸入側へ漏れが発生すると仮定する。
    
        漏れ流量は、圧力差に基づく準定常のオリフィス流れとして表す。
        逆向き、すなわち吸入側から圧縮室への漏れ戻りは考慮しない。
      */
      if leakPort.p > ps then
        m_flow_leak_out = system.CdLeak*system.Aleak*sqrt(2*leakPort.d*abs(leakPort.p - ps));
        m_flow_leak_out_flg = 1;
      else
        m_flow_leak_out = 0;
        m_flow_leak_out_flg = 0;
      end if;
  /*
        漏れポートの質量流量
    
        leakPort.m_flow は、SuctionChamber へ流入する向きを正とする。
        圧縮室から吸入側へ漏れる場合、SuctionChamber へ流入するため正になる。
      */
      leakPort.m_flow = m_flow_leak_out;
      annotation(
        Icon(graphics = {Rectangle(fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(extent = {{-80, 100}, {80, -100}}, textString = "Suction
Chamber")}));
    end SuctionChamber;

    model DischargeBoundary
      outer SystemParameters system;
      
      Interfaces.port_in port annotation(
        Placement(transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}})));
    equation
      /*
        吐出側境界条件
    
        吐出側は、一定圧力 system.pd の圧力境界として扱う。
        このモデルでは、吐出側の容積、温度変化、圧力脈動は考慮しない。
      */
      port.p = system.pd;
  /*
        質量流量の扱い
    
        port.m_flow は、この境界モデルでは直接指定しない。
        接続先のモデル、例えば DischargeValve 側で与えられた流量により、
        connect 方程式を通じて決まる。
    
        したがって、このモデルは流量を決める部品ではなく、
        圧力のみを与える境界条件として機能する。
      */
      annotation(
        Icon(graphics = {Rectangle(fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(extent = {{-80, 100}, {80, -100}}, textString = "Discharge
Boundary")}));
    end DischargeBoundary;

  end Components;

  package Media
    package CO2
      extends ExternalMedia.Media.CoolPropMedium(mediumName = "CO2", substanceNames = {"CO2"});
    end CO2;

    package R32
      extends ExternalMedia.Media.CoolPropMedium(mediumName = "R32", substanceNames = {"R32"});
    end R32;

    package R410A
      extends ExternalMedia.Media.CoolPropMedium(mediumName = "R410A", substanceNames = {"R410A"});
    end R410A;
  end Media;

  package Interfaces
    connector base_port
      flow Real m_flow "Mass flow rate into component [kg/s]";
      Real p "Pressure [Pa]";
      Real h "Specific enthalpy [J/kg]";
      Real d "Density [kg/m3]";
    end base_port;

    connector port_in
      extends base_port;
      annotation(
        Icon(graphics = {Rectangle(fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}})}));
    end port_in;

    connector port_out
      extends base_port;
      annotation(
        Icon(graphics = {Rectangle(fillColor = {255, 0, 0}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}})}));
    end port_out;
    
    connector port_leak
      extends base_port;
      annotation(
        Icon(graphics = {Rectangle(fillColor = {0, 255, 0}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}})}));
    end port_leak;
  end Interfaces;

  package Examples
    model SimpleChamberTest
  inner SystemParameters system annotation(
        Placement(transformation(origin = {-90, 90}, extent = {{-10, -10}, {10, 10}})));
  Components.SuctionChamber suctionChamber annotation(
        Placement(transformation(origin = {-90, 0}, extent = {{-10, -10}, {10, 10}})));
  Components.CompressionChamber compressionChamber annotation(
        Placement(transformation(origin = {-30, 0}, extent = {{-10, -10}, {10, 10}})));
  Components.DischargeValve dischargeValve annotation(
        Placement(transformation(origin = {30, 0}, extent = {{-10, -10}, {10, 10}})));
  Components.DischargeBoundary dischargeBoundary annotation(
        Placement(transformation(origin = {90, 0}, extent = {{-10, -10}, {10, 10}})));
    equation
  connect(compressionChamber.dischargePort, dischargeValve.inlet) annotation(
        Line(points = {{-20, 0}, {20, 0}}));
  connect(dischargeValve.outlet, dischargeBoundary.port) annotation(
        Line(points = {{40, 0}, {80, 0}}));
  connect(suctionChamber.leakPort, compressionChamber.leakPort) annotation(
        Line(points = {{-80, -8}, {-40, -8}}));
  connect(suctionChamber.suctionPort, compressionChamber.suctionPort) annotation(
        Line(points = {{-80, 8}, {-40, 8}}));
    end SimpleChamberTest;
  end Examples;

  model SystemParameters
    /*
      システム共通パラメータ
  
      本モデルは inner/outer 機構により、各コンポーネントから参照される。
      圧縮室、吸入側、吐出側、吐出弁、漏れ流路などで共通に使う
      運転条件および幾何パラメータをここで定義する。
    */
    
    // ------------------------------------------------------------
    // 運転条件
    // ------------------------------------------------------------
    parameter Real f = 60
      "回転周波数 [Hz]";
  
    parameter Real period = 1/f
      "1周期の時間 [s]";
  
    // ------------------------------------------------------------
    // 圧縮室パラメータ
    // ------------------------------------------------------------
    parameter Real Vmin = 0.000001
      "圧縮室の最小体積 [m3]";
  
    parameter Real Vcc = 0.0000565
      "圧縮室の押しのけ容積 [m3]";
  
    parameter Real m0 = 0.000648365
      "圧縮室内ガス質量の初期値 [kg]";
  
    // ------------------------------------------------------------
    // 吸入側・吐出側の境界条件
    // ------------------------------------------------------------
    parameter Real ps = 0.55e6
      "吸入側圧力 [Pa]";
  
    parameter Real Ts = 273.15 + 13.0
      "吸入側温度 [K]";
  
    parameter Real pd = 1.9e6
      "吐出側圧力 [Pa]";
  
    // ------------------------------------------------------------
    // 吐出弁パラメータ
    // ------------------------------------------------------------
    parameter Real CdOut = 0.9
      "吐出流量係数 [-]";
  
    parameter Real Aout = 1.0e-6
      "吐出口有効面積 [m2]";
  
    // ------------------------------------------------------------
    // 漏れ流路パラメータ
    // ------------------------------------------------------------
    parameter Real CdLeak = 0.1
      "漏れ流量係数 [-]";
  
    parameter Real Aleak = 1.0e-8
      "吸入側漏れ相当面積 [m2]";
  
    // ------------------------------------------------------------
    // 主吸入モデルパラメータ
    // ------------------------------------------------------------
    parameter Real Qm = 8.20604e-8
      "主吸入質量流量 [kg/s]";
  
    parameter Real suctionAngle = 10
      "主吸入が発生する角度範囲 [deg]";
  
  initial equation
    /*
      初期方程式
  
      現時点では、ここでは追加の初期条件を指定しない。
      各状態量の初期値は、各コンポーネント内の start 属性で与える。
    */

    annotation(
      defaultComponentName="system",
      defaultComponentPrefixes = "inner",
      Diagram,
      Icon(graphics = {Rectangle(fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}), Text(extent = {{-80, 100}, {80, -100}}, textString = "System
Parameters")}));
  end SystemParameters;

  annotation(
    uses(Modelica(version = "4.1.0")));
end Compression_Chamber;

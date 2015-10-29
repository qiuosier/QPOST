function [ mpc ] = QP_PST2MatPower( bus, line )
%% This function convert data file in PST format to MatPower format
% Default Base MVA = 100

%% MATPOWER Case Format : Version 2
mpc.version = '2';
% system MVA base
mpc.baseMVA = 100;

%% Convert bus data
% Define constants for named column indices of PST
% PST Bus Data
PST_BUS_I   = 1;
PST_VM      = 2;
PST_VA      = 3;
PST_PG      = 4;
PST_QG      = 5;
PST_PD      = 6;
PST_QD      = 7;
PST_GS      = 8;
PST_BS      = 9;
PST_BUS_TYPE= 10;
% The following fields may not exist in PST file
PST_QMAX    = 11;
PST_QMIN    = 12;
PST_BASE_KV = 13;
PST_VMAX    = 14;
PST_VMIN    = 15;
% Number of buses
N_bus = size(bus,1);
% Number of fields in pst bus data
N_field = size(bus,2);
% Use default values to complete missing data in pst bus data
if N_field < 11
    % No Max Q given, set to 3 pu
    bus = [bus, 300 * ones(N_bus,1)];
end
if N_field < 12
    % No Min Q given, set to -3 pu
    bus = [bus, -300 * ones(N_bus,1)];
end
if N_field < 13
    % No base KV given, set to 345 kV
    bus = [bus, 345 * ones(N_bus,1)];
end
if N_field < 14
    % No V MAX given, set to 1.1 pu
    bus = [bus, 1.1 * ones(N_bus,1)];
end
if N_field < 15
    % No V MIN given, set to 0.9 pu
    bus = [bus, 0.9 * ones(N_bus,1)];
end
% Adjust the power values. pu is used in pst while matpower uses MW/MVar
bus(:,[PST_PG,PST_QG,PST_PD,PST_QD]) = ...
    bus(:,[PST_PG,PST_QG,PST_PD,PST_QD]) * mpc.baseMVA;
% Convert bus type from pst to Matpower definition
% PST Def:
PST_REF = 1;
PST_PV  = 2; %#ok<NASGU>
PST_PQ  = 3;
% Matpower Def:
PQ      = 1;
PV      = 2;
REF     = 3;
NONE    = 4; %#ok<NASGU>
for i = 1:N_bus
    switch bus(i,PST_BUS_TYPE)
        case PST_REF
            bus(i,PST_BUS_TYPE) = REF;
        case PST_PQ
            bus(i,PST_BUS_TYPE) = PQ;
    end
end
% Matpower Definition
% Ref: Matpower 4.1 user manual Table B-1, appendix B or Matpower m-files
%   columns 1-13 must be included in input matrix (in case file)
%    1  BUS_I       bus number (positive integer)
%    2  BUS_TYPE    bus type (1 = PQ, 2 = PV, 3 = ref, 4 = isolated)
%    3  PD          Pd, real power demand (MW)
%    4  QD          Qd, reactive power demand (MVAr)
%    5  GS          Gs, shunt conductance (MW demanded at V = 1.0 p.u.)
%    6  BS          Bs, shunt susceptance (MVAr injected at V = 1.0 p.u.)
%    7  BUS_AREA    area number, (positive integer)
%    8  VM          Vm, voltage magnitude (p.u.)
%    9  VA          Va, voltage angle (degrees)
%    10 BASE_KV     baseKV, base voltage (kV)
%    11 ZONE        zone, loss zone (positive integer)
%    12 VMAX        maxVm, maximum voltage magnitude (p.u.)
%    13 VMIN        minVm, minimum voltage magnitude (p.u.)
mpc.bus = [
    bus(:,[PST_BUS_I,PST_BUS_TYPE,PST_PD,PST_QD,PST_GS,PST_BS]), ...% 1-6
    ones(N_bus,1), bus(:,[PST_VM,PST_VA,PST_BASE_KV]), ...          % 7-10
    ones(N_bus,1), bus(:,[PST_VMAX,PST_VMIN])];                     % 11-13

%% Convert generator buses to generator data
% Find Ref and PV buses
[Row_REF, ~] = find(bus(:,PST_BUS_TYPE) == REF);
[Row_PV, ~] = find(bus(:,PST_BUS_TYPE) == PV);
% Generator buses
GenBus = [Row_REF', Row_PV'];
% Number of generators
N_gen = length(GenBus);
% Generator data
% Matpower Definition
% Ref: Matpower 4.1 user manual Table B-2, appendix B or Matpower m-files
%   columns 1-21 must be included in input matrix (in case file)
%    1  GEN_BUS     bus number
%    2  PG          Pg, real power output (MW)
%    3  QG          Qg, reactive power output (MVAr)
%    4  QMAX        Qmax, maximum reactive power output (MVAr)
%    5  QMIN        Qmin, minimum reactive power output (MVAr)
%    6  VG          Vg, voltage magnitude setpoint (p.u.)
%    7  MBASE       mBase, total MVA base of machine, defaults to baseMVA
%    8  GEN_STATUS  status, > 0 - in service, <= 0 - out of service
%    9  PMAX        Pmax, maximum real power output (MW)
%    10 PMIN        Pmin, minimum real power output (MW)
%    11 PC1         Pc1, lower real power output of PQ capability curve (MW)
%    12 PC2         Pc2, upper real power output of PQ capability curve (MW)
%    13 QC1MIN      Qc1min, minimum reactive power output at Pc1 (MVAr)
%    14 QC1MAX      Qc1max, maximum reactive power output at Pc1 (MVAr)
%    15 QC2MIN      Qc2min, minimum reactive power output at Pc2 (MVAr)
%    16 QC2MAX      Qc2max, maximum reactive power output at Pc2 (MVAr)
%    17 RAMP_AGC    ramp rate for load following/AGC (MW/min)
%    18 RAMP_10     ramp rate for 10 minute reserves (MW)
%    19 RAMP_30     ramp rate for 30 minute reserves (MW)
%    20 RAMP_Q      ramp rate for reactive power (2 sec timescale) (MVAr/min)
%    21 APF         area participation factor
mpc.gen = [
    bus(GenBus,[PST_BUS_I, PST_PG, PST_QG, PST_QMAX, PST_QMIN, PST_VM]), ...   % 1-6
    mpc.baseMVA*ones(N_gen,1), ones(N_gen,1), 500*ones(N_gen,1), ...% 7-9
    10*ones(N_gen,1), zeros(N_gen,11)];

%% Convert line/branch data
% Line ratings are set to 250 MVA
% PST Line Data
PST_F_BUS   = 1;
PST_T_BUS   = 2;
PST_BR_R    = 3;
PST_BR_X    = 4;
PST_BR_B    = 5;
PST_TAB     = 6;
PST_SHIFT   = 7;
% Number of lines
N_line = size(line,1);
% Matpower Definition
% Ref: Matpower 4.1 user manual Table B-2, appendix B or Matpower m-files
%   columns 1-11 must be included in input matrix (in case file)
%    1  F_BUS       f, from bus number
%    2  T_BUS       t, to bus number
%    3  BR_R        r, resistance (p.u.)
%    4  BR_X        x, reactance (p.u.)
%    5  BR_B        b, total line charging susceptance (p.u.)
%    6  RATE_A      rateA, MVA rating A (long term rating)
%    7  RATE_B      rateB, MVA rating B (short term rating)
%    8  RATE_C      rateC, MVA rating C (emergency rating)
%    9  TAP         ratio, transformer off nominal turns ratio
%    10 SHIFT       angle, transformer phase shift angle (degrees)
%    11 BR_STATUS   initial branch status, 1 - in service, 0 - out of service
%    12 ANGMIN      minimum angle difference, angle(Vf) - angle(Vt) (degrees)
%    13 ANGMAX      maximum angle difference, angle(Vf) - angle(Vt) (degrees)
mpc.branch = [
    line(:, [PST_F_BUS, PST_T_BUS, PST_BR_R, PST_BR_X, PST_BR_B]) , ...
    250 * ones(N_line,3), line(:,[PST_TAB, PST_SHIFT]), ...
    ones(N_line,1), -360*ones(N_line,1), 360*ones(N_line,1)];

end


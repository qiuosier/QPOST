function [ bus ] = QP_SolvePowerFlow( bus, line, printResult )
% This function uses Matpower to solve power flow
% Input/output data are in PST format
% printResult = 1 : Print power flow results
% printResult = 0 : No output

% Solve Power Flow
% Convert to Matpower format
mpc = QP_PST2MatPower(bus, line);
% Solve power flow
mpopt = mpoption('OUT_ALL',-printResult,'VERBOSE',printResult);
PFSolution = runpf(mpc, mpopt);
% Update power system data from power flow results
% Define PST constants
PST_VM      = 2;
PST_VA      = 3;
PST_PG      = 4;
PST_QG      = 5;
% Define MP constants
MP_VM       = 8;
MP_VA       = 9;
MP_GEN_BUS  = 1;
MP_PG       = 2;
MP_QG       = 3;
% Update bus voltage and angle
bus(:,[PST_VM, PST_VA]) = PFSolution.bus(:,[MP_VM, MP_VA]);
% Update real and reactive power
bus(PFSolution.gen(:, MP_GEN_BUS)', [PST_PG, PST_QG]) = ... 
    PFSolution.gen(:, [MP_PG, MP_QG]) / mpc.baseMVA;

end


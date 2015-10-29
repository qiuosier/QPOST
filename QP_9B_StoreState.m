% Version $\Delta$. Qiu Qin, December 8, 2014. All Rights Reserved.
% This script is for QP_9B_RunModel.m
% This file should not be used independently

% Store final values of state variabls
deltaState = deltaValues;
omegaState = omegaValues.Data(end,:);
%RC_State0 = RC_State.Data(end,:);
% Store mechanical states (Voltages and Currents)
if size(MS_delta.Time,1)
    ddelta.Time = ddelta.Time + MS_delta.Time(end);
    omegaValues.Time = omegaValues.Time + MS_w.Time(end);
end
MS_delta = append(MS_delta,ddelta);
MS_w = append(MS_w, omegaValues);
% Store electrical states (Voltages and Currents) and outputs
if size(ES_output.Time,1)
    y_output.Time = y_output.Time + ES_output.Time(end);
    u_input.Time = u_input.Time + ES_input.Time(end);
    RC_State.Time = RC_State.Time + ES_state.Time(end);
end
% Concatenate the inputs, outputs and states
ES_output = append(ES_output,y_output);
ES_input = append(ES_input,u_input);
ES_state = append(ES_state,RC_State);
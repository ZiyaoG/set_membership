% Assuming Omega_t, Omega_t1 (Omega_{t+1}), M_t, and B are given and known:
% Omega_t = [...];
% Omega_t1 = [...];
% M_t = [...];
% B = [...];
clear all;
close all;
load("datas\quadrotors\sincos05.mat");

% Define the variables
J = sdpvar(3,3,'full');

% LMI
LMI = J * (Omega_t1 - Omega_t) - M_t - Omega_t * J * Omega_t - B <= 0;

% Solve the LMI
options = sdpsettings('verbose', 1, 'solver', 'lmilab'); % Choose a solver that is available and appropriate
sol = solvesdp(LMI, [], options);

% Check if the solution is found
if sol.problem == 0
    % Extract the solution
    J_value = double(J);
    disp('Solution found:');
    disp(J_value);
else
    disp('No solution found');
end
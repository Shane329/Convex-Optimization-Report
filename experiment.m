% Clear workspace and command window
clear;
clc;

% Define the ELD problem data
NG = 10;  % Number of generators

% Cost coefficients [c_i, b_i, a_i] for quadratic cost: c_i*P_i^2 + b_i*P_i + a_i
cost_coeff = [0.007, 7, 240;   % Generator 1
              0.0095, 10, 200;  % Generator 2
              0.009, 8.5, 220;  % Generator 3
              0.009, 11, 200;   % Generator 4
              0.008, 10.5, 220; % Generator 5
              0.0075, 12, 120;  % Generator 6
              0.0075, 14, 130;  % Generator 7
              0.0075, 14, 130;  % Generator 8
              0.0075, 14, 130;  % Generator 9
              0.0075, 14, 130]; % Generator 10

% Generator limits [P_min, P_max] in MW
P_limits = [0.66, 3.35;   % Generator 1
            0.9, 3.7;     % Generator 2
            0.8, 3.6;     % Generator 3
            0.66, 3.35;   % Generator 4
            0.72, 3.45;   % Generator 5
            0.66, 2.97;   % Generator 6
            0.88, 3.5;    % Generator 7
            0.754, 3.33;  % Generator 8
            0.9, 3.9;     % Generator 9
            0.56, 2.35];  % Generator 10

% Total power demand (P_load) in MW
P_load = sum(P_limits(:,1)) + 0.5*(sum(P_limits(:,2)) - sum(P_limits(:,1))); % Midpoint demand

% Loss matrix (B-coefficients) scaled by 1e-4
B = [0.14 0.17 0.15 0.19 0.26 0.22 0.22 0.19 0.26 0.15;
     0.17 0.6  0.13 0.16 0.15 0.2  0.2  0.7  0.15 0.13;
     0.15 0.13 0.65 0.17 0.24 0.19 0.19 0.13 0.24 0.65;
     0.19 0.16 0.17 0.71 0.3  0.25 0.25 0.18 0.3  0.17;
     0.26 0.15 0.24 0.3  0.69 0.32 0.32 0.16 0.69 0.24;
     0.22 0.2  0.19 0.25 0.32 0.85 0.85 0.21 0.32 0.19;
     0.34 0.23 0.25 0.43 0.18 0.97 0.67 0.28 0.18 0.25;
     0.38 0.56 0.38 0.56 0.37 0.55 0.38 0.56 0.37 0.38;
     0.43 0.23 0.43 0.23 0.42 0.27 0.43 0.23 0.42 0.43;
     0.45 0.51 0.45 0.51 0.48 0.58 0.45 0.51 0.48 0.45] * 1e-4;

B0 = zeros(10,1);  % B0i coefficients
B00 = 0;           % Constant loss term

% Noise levels (percentage deviation)
%noise_levels = [0, 0.1, 0.2, 0.3, 0.4, 0.5]; % Low, medium, high noise
num_trials = 50; % Number of trials for each noise level
noise_levels = linspace(0,0.5,21);

% Initialize results
results_qp = struct();
results_lambda = struct();

%%

% Run the experiment for each noise level
for i = 1:length(noise_levels)
    sigma = noise_levels(i);
    costs_qp = zeros(1, num_trials);
    costs_lambda = zeros(1, num_trials);
    
    % Replace '.' with '_' in the field name
    field_name = ['sigma_', strrep(num2str(sigma), '.', '_')];
    
    for trial = 1:num_trials
        % Generate lognormal noise (always positive, multiplicative)
        noise = exp(sigma * randn(size(cost_coeff)));
        
        % Apply noise to coefficients
        noisy_cost_coeff = cost_coeff .* noise;
        
        % Solve using QP
        [P_opt_qp, cost_qp] = solve_eld_qp(noisy_cost_coeff, P_limits, P_load, B, B0, B00);
        costs_qp(trial) = cost_qp;
        
        % Solve using Lambda Iteration
        [P_opt_lambda, cost_lambda] = solve_eld_lambda_with_losses(noisy_cost_coeff, P_limits, P_load, B, B0, B00);
        costs_lambda(trial) = cost_lambda;
    end
    
    % Store results
    results_qp.(field_name) = costs_qp;
    results_lambda.(field_name) = costs_lambda;
end


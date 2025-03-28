% QP Solver Function
function [P_opt, cost] = solve_eld_qp(cost_coeff, P_limits, P_load, B, B0, B00)
    % Ensure cost_coeff is properly structured: [a_i, b_i, c_i]
    a = cost_coeff(:, 1);
    b = cost_coeff(:, 2);
    
    % Construct QP matrices
    H = 2 * diag(a);
    f = b;
    lb = P_limits(:, 1);
    ub = P_limits(:, 2);
    
    % Check convexity
    if any(a <= 0)
        error('Non-convex cost function: Quadratic coefficients must be positive.');
    end
    
    % Initial guess (midpoint of bounds)
    P_init = mean(P_limits, 2);
    
    % Nonlinear constraint for power balance with losses
    power_balance = @(P) deal([], P'*B*P + B0'*P + B00 - P_load);
    
    % Solver options
    options = optimoptions('fmincon', ...
                          'Display', 'off', ...
                          'Algorithm', 'interior-point', ...
                          'MaxIterations', 1000);
    
    % Solve
    try
        P_opt = fmincon(@(P) 0.5*P'*H*P + f'*P, P_init, ...
                        [], [], [], [], lb, ub, power_balance, options);
        cost = sum(a .* P_opt.^2 + b .* P_opt + cost_coeff(:, 3));
    catch ME
        warning('QP failed: %s', ME.message);
        P_opt = P_init;  % Fallback to initial guess
        cost = NaN;
    end
end

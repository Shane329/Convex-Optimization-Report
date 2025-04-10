% QP Solver Function

function [P_opt, cost, iter] = solve_eld_qp_iterative(cost_coeff, P_limits, P_load, B, B0, B00, tol, max_iter)
    % Initialize
    P_prev = mean(P_limits, 2);
    iter = 0;
    converged = false;
    tol = 1e-4;                % Tolerance for power balance
    max_iter = 250;            % Maximum iterations
    
    while ~converged && iter < max_iter
        iter = iter + 1;
        
        % Solve with linearized constraint
        [P_opt, cost, exitflag] = solve_eld_qp_linearized(cost_coeff, P_limits, P_load, B, B0, B00, P_prev);
        
        % Check convergence
        if norm(P_opt - P_prev) < tol
            converged = true;
        end
        
        P_prev = P_opt;
        
        if exitflag <= 0
            break;  % Exit if QP failed
        end
    end
    
    if ~converged
        warning('Iterative QP did not fully converge after %d iterations', iter);
    end
end

function [P_opt, cost, exitflag] = solve_eld_qp_linearized(cost_coeff, P_limits, P_load, B, B0, B00, P_prev)
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
    
    % Initial guess (use previous solution if available, otherwise midpoint)
    if nargin < 7 || isempty(P_prev)
        P_prev = mean(P_limits, 2);
    end
    
    % Linearize the power balance constraint around P_prev
    Aeq = 2 * P_prev' * B + B0';
    beq = P_load - (B00 - P_prev' * B * P_prev);
    
    % Solver options
    options = optimoptions('quadprog', ...
                          'Display', 'off', ...
                          'MaxIterations', 1000);
    
    % Solve QP
    try
        [P_opt, ~, exitflag] = quadprog(H, f, [], [], Aeq, beq, lb, ub, P_prev, options);
        
        if exitflag <= 0
            warning('QP may not have converged properly (exitflag = %d)', exitflag);
            P_opt = P_prev;  % Fallback to previous solution
        end
        
        cost = sum(a .* P_opt.^2 + b .* P_opt + cost_coeff(:, 3));
    catch ME
        warning('QP failed: %s', ME.message);
        P_opt = P_prev;  % Fallback to previous solution
        cost = NaN;
        exitflag = -1;
    end
end

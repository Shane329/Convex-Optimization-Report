% Lambda Iteration Solver Function WITH LOSSES (Fixed version)
function [P_opt, cost_opt] = solve_eld_lambda_with_losses(cost_coeff, P_limits, P_load, B, B0, B00)
    NG = size(cost_coeff, 1);
    
    % Lambda iteration parameters
    lambda_initial = 20;       % Initial guess for λ (typical marginal cost range)
    delta_lambda = 1;          % Initial step size for λ updates
    tol = 1e-4;                % Tolerance for power balance
    max_iter = 250;            % Maximum iterations
    lambda_min = 0;            % Absolute minimum λ
    lambda_max = 100;          % Absolute maximum λ
    
    % Initialize variables
    lambda = lambda_initial;
    iter = 0;
    converged = false;
    P = zeros(NG, 1);
    
    while ~converged && iter < max_iter
        % Calculate generator outputs with penalty factors
        for i = 1:NG
            % Compute incremental loss (dP_L/dP_i)
            dPdP_i = 2*sum(B(i,:).*P') + B0(i);
            
            % Calculate penalty factor L_i
            L_i = 1/(1 - dPdP_i);
            
            % Update generator output
            P(i) = (lambda/L_i - cost_coeff(i,2)) / (2*cost_coeff(i,1));
            
            % Enforce generator limits
            P(i) = max(P(i), P_limits(i,1));
            P(i) = min(P(i), P_limits(i,2));
        end
        
        % Calculate total losses
        P_loss = P' * B * P + B0' * P + B00;
        
        % Check power balance (ΔP = total_gen - (load + losses))
        DeltaP = sum(P) - (P_load + P_loss);
        
        if abs(DeltaP) < tol
            converged = true;
        else
            % Update lambda using your specified logic
            if DeltaP > 0  % Over-generation: reduce λ
                lambda = lambda - delta_lambda;
            else           % Under-generation: increase λ
                lambda = lambda + delta_lambda;
            end
            
            % Optional: Adaptive step size refinement
            delta_lambda = delta_lambda * 0.95;  % Reduce step size each iteration
            
            % Enforce absolute bounds (safeguard)
            lambda = max(lambda_min, min(lambda_max, lambda));
        end
        
        iter = iter + 1;
    end
    
    % Calculate final outputs
    P_opt = P;
    cost_opt = sum(cost_coeff(:,1).*P.^2 + cost_coeff(:,2).*P + cost_coeff(:,3));
    
    if ~converged
        warning('Lambda iteration did not converge in %d iterations (Final ΔP = %.2e)', ...
                max_iter, DeltaP);
    end
end

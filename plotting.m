%%

% Extract and compute statistics for each noise level
sigma_values = linspace(0, 0.5, 21);
mean_qp = zeros(size(sigma_values));
std_qp = zeros(size(sigma_values));
mean_lambda = zeros(size(sigma_values));
std_lambda = zeros(size(sigma_values));

for i = 1:length(sigma_values)
    sigma = sigma_values(i);
    field_name = ['sigma_', strrep(num2str(sigma), '.', '_')];
    
    % QP statistics
    qp_costs = results_qp.(field_name);
    mean_qp(i) = mean(qp_costs);
    std_qp(i) = std(qp_costs);
    
    % Lambda statistics
    lambda_costs = results_lambda.(field_name);
    mean_lambda(i) = mean(lambda_costs);
    std_lambda(i) = std(lambda_costs);
end

% Create figure with shaded error regions
figure;
hold on;
grid on;

% QP plot with shaded region
qp_plot = plot(sigma_values, mean_qp, 'b-', 'LineWidth', 2, 'DisplayName', 'QP');
fill([sigma_values, fliplr(sigma_values)], ...
     [mean_qp - std_qp, fliplr(mean_qp + std_qp)], ...
     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Lambda plot with shaded region
lambda_plot = plot(sigma_values, mean_lambda, 'r--', 'LineWidth', 2, 'DisplayName', 'Lambda Iteration');
fill([sigma_values, fliplr(sigma_values)], ...
     [mean_lambda - std_lambda, fliplr(mean_lambda + std_lambda)], ...
     'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Formatting
xlabel('Noise Level (\sigma)');
ylabel('Optimal Cost ($)');
title('Mean Optimal Cost with \pm1 Standard Deviation');
legend([qp_plot, lambda_plot], 'Location', 'best');
set(gca, 'FontSize', 12);
xlim([0, 0.5]);
box on;

% Optional: Add annotation for interpretation
text(0.3, min(mean_lambda - std_lambda), ...
    'Shaded regions show \pm1\sigma variance', ...
    'FontSize', 10, 'BackgroundColor', 'white');

%%
% Calculate mean cost differences (ΔC = QP - LI)
sigma_values = linspace(0, 0.5, 21);
delta_C = zeros(size(sigma_values));

for i = 1:length(sigma_values)
    sigma = sigma_values(i);
    field_name = ['sigma_', strrep(num2str(sigma), '.', '_')];
    
    % Get mean costs (now QP - LI)
    mean_QP = mean(results_qp.(field_name));
    mean_LI = mean(results_lambda.(field_name));
    
    % Store difference
    delta_C(i) = mean_QP - mean_LI;
end

% Quadratic regression: ΔC = ασ² + βσ + γ
X = [sigma_values'.^2, sigma_values', ones(size(sigma_values'))];
coefficients = X\delta_C';
alpha = coefficients(1);
beta = coefficients(2);
gamma = coefficients(3);

% Calculate R-squared
y_pred = X*coefficients;
SS_res = sum((delta_C' - y_pred).^2);
SS_tot = sum((delta_C' - mean(delta_C)).^2);
R_squared = 1 - (SS_res/SS_tot);

% Create the plot
figure;
hold on;
grid on;

% Scatter plot of actual differences
scatter(sigma_values, delta_C, 40, 'filled', 'MarkerFaceColor', [0.2 0.4 0.8]);

% Regression line
x_fit = linspace(0, 0.5, 100);
y_fit = alpha*x_fit.^2 + beta*x_fit + gamma;
plot(x_fit, y_fit, 'r-', 'LineWidth', 2);

% Annotation with equation (fixed display)
eqn_text = sprintf('\\DeltaC = %.3fσ² + %.3fσ + %.3f\nR² = %.3f', ...
                  alpha, beta, gamma, R_squared);
text(0.05, max(delta_C)*0.9, eqn_text, ...
     'FontSize', 11, ...
     'BackgroundColor', 'white', ...
     'EdgeColor', 'k');

% Labels and title
xlabel('Noise Level (σ)', 'FontSize', 12);
ylabel('Cost Difference \DeltaC = C_{QP} - C_{LI} ($)', 'FontSize', 12);
title('QP vs. Lambda Iteration Performance Gap', 'FontSize', 14);
legend('Actual \DeltaC', 'Quadratic Fit', 'Location', 'best');

%%
% Calculate std deviation differences (QP - LI)
sigma_values = linspace(0, 0.5, 21);
std_diff = arrayfun(@(s) std(results_qp.(['sigma_', strrep(num2str(s), '.', '_')])) - ...
                   std(results_lambda.(['sigma_', strrep(num2str(s), '.', '_')])), sigma_values);

% Regular quadratic fit (Δσ = ασ² + βσ + γ)
X = [sigma_values'.^2, sigma_values', ones(21,1)];
coefficients = X \ std_diff';
alpha = coefficients(1);
beta = coefficients(2);
gamma = coefficients(3);

% Calculate fitted curve
x_fit = linspace(0, 0.5, 100);
y_fit = alpha*x_fit.^2 + beta*x_fit + gamma;

% Calculate R²
y_pred = X*coefficients;
SS_res = sum((std_diff' - y_pred).^2);
SS_tot = sum((std_diff' - mean(std_diff)).^2);
R_squared = 1 - (SS_res/SS_tot);

% Create clean plot
figure;
hold on;
grid on;

% Plot data points
scatter(sigma_values, std_diff, 40, 'filled', 'MarkerFaceColor', [0.2 0.4 0.8]);

% Plot quadratic fit
plot(x_fit, y_fit, 'r-', 'LineWidth', 2);

% Add equation and R²
text(0.05, max(std_diff)*0.8, ...
     sprintf('\\Deltaσ = %.3fσ² + %.3fσ + %.3f\nR² = %.3f', alpha, beta, gamma, R_squared), ...
     'FontSize', 11, 'BackgroundColor', 'white');

% Labels and title
xlabel('Noise Level (σ)');
ylabel('Standard Deviation Difference (σ_{QP} - σ_{LI}) ($)');
title('Quadratic Fit of Standard Deviation Differences');
legend('Data', 'Quadratic Fit', 'Location', 'best');

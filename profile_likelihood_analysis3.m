function profile_likelihood_analysis3(dataset, rnacol, gdmscol, T7conc)
    % Load the existing fitting function
    fit_params = fittinglsq2(dataset, rnacol, gdmscol, T7conc);

    % Load your data
    d = xlsread(dataset);
    t_data = d(:,1);
    mRNA_data = d(:,rnacol);
    GdmS_star_data = d(:,gdmscol);

    % Combine data into a single matrix for fitting
    data = [mRNA_data, GdmS_star_data];

    % Define the constants [T7, DNA, RNase]
    T7 = T7conc;
    DNA = 5.6;
    RNase = 270;
    constants = [T7, DNA, RNase];

    % Define initial conditions for the ODEs, with placeholders for TsR and TlR
    initial_conditions = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];

    % Define lower bounds for the parameters
    lb = [0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 100, 100];

    % Define upper bounds for the parameters
    ub = [10, 10, 10, 10, 10, 10, 10, 10, 10, 1000, 1000];

    % Define the objective function for fitting
    objective_function = @(params) fit_odes(params, t_data, initial_conditions, constants, data);

    % Number of parameters
    num_params = length(fit_params);

    % Range for parameter variation (±1000% of the fitted values)
    range_factor = 10;

    % Initialize storage for profile likelihood results
    num_points = 20;
    profile_likelihood = zeros(num_points, num_params);
    parameter_values = zeros(num_points, num_params);

    % Generate parameter ranges
    for i = 1:num_params
        param_range = linspace(fit_params(i) * (1 - range_factor), fit_params(i) * (1 + range_factor), num_points);
        parameter_values(:, i) = param_range';
    end

    % Perform profile likelihood analysis using proflik
    for i = 1:num_params
        param_range = parameter_values(:, i);
        other_params_initial = fit_params(setdiff(1:num_params, i));
        obj_fixed = @(fixed_value) obj_func_with_fixed_param(objective_function, fixed_value, i, other_params_initial, lb, ub, num_params);
        [~, prof_lik, ~] = proflik(obj_fixed, param_range);
        profile_likelihood(:, i) = prof_lik;
    end

    % Plot the profile likelihood results
    figure;
    for i = 1:num_params
        subplot(ceil(num_params / 2), 2, i);
        plot(parameter_values(:, i), profile_likelihood(:, i), '-o');
        xlabel(['Parameter ', num2str(i)]);
        ylabel('Profile Likelihood');
        title(['Profile Likelihood for Parameter ', num2str(i)]);
        grid on;
    end
end

% Function to compute the objective with one parameter fixed
function [nll, other_params] = obj_func_with_fixed_param(objective_function, fixed_value, fixed_index, other_params_initial, lb, ub, num_params)
    fixed_params = zeros(1, num_params);
    fixed_params(fixed_index) = fixed_value;
    other_params = lsqnonlin(@(other_params) objective_function(update_params(fixed_params, other_params, fixed_index)), other_params_initial, lb(setdiff(1:num_params, fixed_index)), ub(setdiff(1:num_params, fixed_index)), optimoptions('lsqnonlin', 'Display', 'off', 'MaxFunctionEvaluations', 100000, 'MaxIterations', 100000));
    nll = objective_function(update_params(fixed_params, other_params, fixed_index));
end

% Function to update parameters for optimization
function updated_params = update_params(fixed_params, other_params, fixed_index)
    updated_params = fixed_params;
    updated_params(setdiff(1:length(fixed_params), fixed_index)) = other_params;
end

function residuals = fit_odes(params, t_data, initial_conditions, constants, data)
    % Update initial conditions with fitted values of TsR and TlR
    initial_conditions(8) = params(10);
    initial_conditions(9) = params(11);
    % Solve the ODE
    [~, y] = ode45(@(t,y) coupled_odes(t, y, params(1:9), constants), t_data, initial_conditions);
    % Extract model predictions
    model_predictions = [y(:,2), y(:,7)];
    % Calculate residuals
    residuals = model_predictions - data;
end
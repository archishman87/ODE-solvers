function profile_likelihood_analysis(dataset, rnacol, gdmscol, T7conc)
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
    ub = [100, 100, 100, 100, 100, 100, 100, 100, 100, 1000, 1000];

    % Define the objective function for fitting
    objective_function = @(params) fit_odes(params, t_data, initial_conditions, constants, data);

    % Number of parameters
    num_params = length(fit_params);

    % Range for parameter variation (±1000% of the fitted values)
    range_factor = 0.1;

    % Initialize storage for profile likelihood results
    num_points = 10;
    profile_likelihood = zeros(num_points, num_params);
    parameter_values = zeros(num_points, num_params);

    for i = 1:num_params
        % Vary parameter i and optimize the others
        param_range = linspace(fit_params(i) * (1 - range_factor), fit_params(i) * (1 + range_factor), num_points);
        for j = 1:length(param_range)
            try
                % Fix parameter i
                fixed_params = fit_params;
                fixed_params(i) = param_range(j);

                % Optimize other parameters
                other_params_initial = fit_params(setdiff(1:num_params, i));
                obj_fixed = @(other_params) objective_function(update_params(fixed_params, other_params, i));
                options = optimoptions('lsqnonlin', 'Display', 'iter', 'MaxFunctionEvaluations', 100000, 'MaxIterations', 100000);
                [optimized_params, resnorm] = lsqnonlin(obj_fixed, other_params_initial, lb(setdiff(1:num_params, i)), ub(setdiff(1:num_params, i)), options);

                % Store profile likelihood and parameter values
                profile_likelihood(j, i) = resnorm;
                parameter_values(j, i) = param_range(j);
            catch
                % Handle ODE solver failures by setting high residuals
                profile_likelihood(j, i) = Inf;
                parameter_values(j, i) = param_range(j);
                disp(['ODE solver failed for parameter ', num2str(i), ' at value ', num2str(param_range(j))]);
            end
        end
    end

    % Plot the profile likelihood results
    figure;
    for i = 1:num_params
        subplot(ceil(num_params / 2), 2, i);
        plot(parameter_values(:, i), profile_likelihood(:, i), 'o');
        xlabel(['Parameter ', num2str(i)]);
        ylabel('Profile Likelihood');
        title(['Profile Likelihood for Parameter ', num2str(i)]);
        grid on;
    end
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
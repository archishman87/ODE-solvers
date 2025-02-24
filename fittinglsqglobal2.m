% Main script for fitting ODEs to data
function parameters = fittinglsqglobal2(dataset, rnacol_start, gdmscol_start, T7concs)
    % Number of datasets
    num_datasets = length(T7concs);
    
    % Load your data
    d = xlsread(dataset);
    t_data = d(:,1);
    mRNA_data = d(:,rnacol_start:rnacol_start + num_datasets - 1);
    GdmS_star_data = d(:,gdmscol_start:gdmscol_start + num_datasets - 1);

    % Initialize cell arrays to store data
    data = cell(num_datasets, 1);
    constants = cell(num_datasets, 1);
    
    % Combine data into cell arrays
    for i = 1:num_datasets
        data{i} = [mRNA_data(:,i), GdmS_star_data(:,i)];
        constants{i} = [T7concs(i), 5.6, 270]; % [T7, DNA, RNase]
    end

    % Define initial guesses for the parameters [kr, kb1, kb2, kp, k1, k2, k3, kmat, kd, TsR_initial, TlR_initial]
    initial_params = [0.01, 5, 1, 2, 0.1, 0.1, 0.1, 0.1, 0.01, 300.0, 300.0];

    % Define lower bounds for the parameters
    lb = [0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 10, 10];

    % Define upper bounds for the parameters
    ub = [10, 10, 10, 10, 100, 10, 10, 10, 10, 1000, 1000];

    % Open a parallel pool
    if isempty(gcp('nocreate'))
        parpool; % You can specify the number of workers as an argument, e.g., parpool('local', 4)
    end

    % Define the objective function for fitting
    objective_function = @(params) fit_odes(params, t_data, constants, data, num_datasets);

    % Perform the fitting using lsqnonlin
    options = optimoptions('lsqnonlin', 'Display', 'iter', 'MaxFunctionEvaluations', 100000, 'MaxIterations', 100000);
    [fit_params, resnorm] = lsqnonlin(objective_function, initial_params, lb, ub, options);

    format long g
    % 
    % % Display the fitted parameters
    % disp('Fitted Parameters:');
    % disp(fit_params(1:9));
    % 
    % % Display the initial values of TsR and TlR
    % disp('Fitted initial value of TsR at time = 0:');
    % disp(fit_params(10));
    % disp('Fitted initial value of TlR at time = 0:');
    % disp(fit_params(11));

    parameters = fit_params';

    % Update initial conditions with fitted values of TsR and TlR
    initial_conditions = zeros(9, 1);
    initial_conditions(8) = fit_params(10);
    initial_conditions(9) = fit_params(11);

    T73d = T7concs.*ones(length(t_data),1);
    % Plot the fit results
    figure;
    export_data = cell(num_datasets, 1);
    
    for i = 1:num_datasets
        [t_fit, y_fit] = ode45(@(t,y) coupled_odes(t, y, fit_params(1:9), constants{i}), t_data, initial_conditions);
        export_data{i} = [t_data, T73d(:,i), mRNA_data(:,i), GdmS_star_data(:,i), y_fit(:,2), y_fit(:,7)];
        plot3(t_data, T73d(:,i),mRNA_data(:,i), 'g.','MarkerSize',15, 'DisplayName', 'mRNA data');
        hold on;
        plot3(t_data, T73d(:,i),GdmS_star_data(:,i), 'm.','MarkerSize',15, 'DisplayName', 'GdmS* data');
        plot3(t_fit, T73d(:,i),y_fit(:,2), 'g-','LineWidth',2, 'DisplayName', 'mRNA fit');
        plot3(t_fit, T73d(:,i),y_fit(:,7), 'm-','LineWidth',2, 'DisplayName', 'GdmS* fit');

        xlabel('Time (hours)');
        ylabel('[T7] (nM)');
        zlabel('[mRNA] or [GdmS^*] (nM)')
        %title(['Fitting ODE Model to Data for T7 = ', num2str(T7concs(i))]);
        %legend('show');
        grid on;
        box on;
    end
    % Prepare data for export
    
    % Combine all data into one table
    combined_data = cell2mat(export_data);
    combined_table = array2table(combined_data, 'VariableNames', ...
        {'Time', 'T7', 'mRNA_data', 'GdmS_star_data', 'mRNA_fit', 'GdmS_star_fit'});

    % Write the table to an Excel file
    writetable(combined_table, 'globaloutput-avg.csv');
end

% Function to compute the model predictions and residuals
function residuals = fit_odes(params, t_data, constants, data, num_datasets)
    residuals = [];
    for i = 1:num_datasets
        % Update initial conditions with fitted values of TsR and TlR
        initial_conditions = zeros(9, 1);
        initial_conditions(8) = params(10);
        initial_conditions(9) = params(11);
        % Solve the ODE
        [~, y] = ode45(@(t,y) coupled_odes(t, y, params(1:9), constants{i}), t_data, initial_conditions);
        % Extract model predictions
        model_predictions = [y(:,2), y(:,7)];
        % Calculate residuals
        residuals = [residuals; model_predictions - data{i}];
    end

    
end

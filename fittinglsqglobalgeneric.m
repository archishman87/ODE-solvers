% Main script for fitting ODEs to data
function parameters = fittinglsqglobalgeneric(dataset, rnacol_start, gdmscol_start)

    T7concs = [470 1470 1470 470];
    RNaseconcs = [270 270 370 370];
    DNA = 5.6;
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
        constants{i} = [T7concs(i), DNA, RNaseconcs(i)]; % [T7, DNA, RNase]
    end

    % Define initial guesses for the parameters [kr, kb1, kb2, kp, k1, k2, k3, kmat, kd, TsR_initial, TlR_initial]
    initial_params = [0.001, 5, 1, 2, 0.1, 0.1, 0.1, 0.1, 0.01, 300.0, 300.0];

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

    TR3d = T7concs./RNaseconcs.*ones(length(t_data),1);
    % Plot the fit results
    figure;
    export_data = cell(num_datasets, 1);
    
    for i = 1:num_datasets
        [t_fit, y_fit] = ode45(@(t,y) coupled_odes(t, y, fit_params(1:9), constants{i}), t_data, initial_conditions);
        export_data{i} = [t_data, TR3d(:,i), mRNA_data(:,i), GdmS_star_data(:,i), y_fit(:,2), y_fit(:,7)];
        plot3(t_data, TR3d(:,i),mRNA_data(:,i)*10, 'g.','MarkerSize',15, 'DisplayName', 'mRNA data');
        hold on;
        plot3(t_data, TR3d(:,i),GdmS_star_data(:,i), 'm.','MarkerSize',15, 'DisplayName', 'GdmS* data');
        plot3(t_fit, TR3d(:,i),y_fit(:,2)*10, 'g-','LineWidth',2, 'DisplayName', 'mRNA fit');
        plot3(t_fit, TR3d(:,i),y_fit(:,7), 'm-','LineWidth',2, 'DisplayName', 'GdmS* fit');

        %xlabel('Time (hours)');
        %ylabel('[DNA] (nM)');
        %zlabel('[mRNA] or [GdmS^*] (nM)')
        %title(['Fitting ODE Model to Data for T7 = ', num2str(T7concs(i))]);
        %legend('show');
        grid on;
        box on;
    end
    % Prepare data for export

    % Time span for the simulation
    tspan = [0 6];
    num_time_points = 1000; % Number of points in time

    % Prepare figure
    %figure;
    hold on;

    % Initialize matrices for meshgrid
    TR_grid = [];
    t_grid = [];
    mRNA_grid = [];
    GdmS_star_grid = [];
    TRconcs = T7concs./RNaseconcs;
    
    for i = 1:length(T7concs)
        T7 = T7concs(i);
        RNase = RNaseconcs(i);
        constants = [T7, DNA, RNase];
        % Solve the ODE system
        [t, y] = ode45(@(t,y) coupled_odes(t, y, parameters, constants), linspace(tspan(1), tspan(2), num_time_points), initial_conditions);

        % Store results in grid matrices
        TR_grid = [TR_grid; TRconcs.*ones(size(t))'];
        t_grid = [t_grid; t'];
        mRNA_grid = [mRNA_grid; y(:,2)'];
        GdmS_star_grid = [GdmS_star_grid; y(:,7)'];
    end

    % Convert grids to mesh format
    [TR_mesh, t_mesh] = meshgrid(TRconcs, linspace(tspan(1), tspan(2), num_time_points));

    % Interpolate the data
    mRNA_mesh = griddata(t_grid(:), TR_grid(:), mRNA_grid(:), t_mesh, TR_mesh);
    GdmS_star_mesh = griddata(t_grid(:), TR_grid(:), GdmS_star_grid(:), t_mesh, TR_mesh);

    % Plot mRNA mesh
    surf(t_mesh, TR_mesh, mRNA_mesh*10,'FaceColor','g', 'EdgeColor', 'none', 'FaceAlpha', 0.2, 'DisplayName', 'mRNA');
    hold on;
    % Plot GdmS_star mesh
    surf(t_mesh, TR_mesh, GdmS_star_mesh,'FaceColor','m' ,'EdgeColor', 'none', 'FaceAlpha', 0.2, 'DisplayName', 'GdmS^*');

    % Customize the plot
    xlabel('Time (hours)');
    ylabel('[T7]/[RNase]');
    zlabel(' '); % Leave the zlabel empty for custom text placement
    % Get current axes
    ax = gca;
    zlim([0 500])

    % Get the position of the zlabel for custom text placement
    labelPos = ax.ZLabel.Position;
    % Add multi-colored text near the z-axis label
text(labelPos(1), 3, 75, ...
    ['\color{green}[mRNA] \color{black}or \color{magenta}[GdmS^*] \color{black}(nM)'], ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'bottom', ...
    'Rotation', 90, ...
    'FontSize', 15, ...
    'FontWeight','bold');

% Adjust plot limits to ensure the text is visible
% zlim([min(z) - 1, max(z) + 1]);
    % zlabel('[mRNA] or [GdmS^*] (nM)');
    %title('3D Plot of mRNA and GdmS^* Concentrations');
    %legend('show');
    grid on;
    
    view(3);
    box on
    ax.GridAlpha = 0.25;
    ax.BoxStyle = 'full';
    set(ax,'FontSize',15,'FontWeight','bold')
    rotate3d on
    pbaspect([1 1 1])
    view(-60,15)
    xticks([0 1 2 3 4 5 6])
    %yticks([0 500 1000 1500 2000])
    zticks([0:50:500])
    % Adjust positions of the labels
% Adjust positions of the labels
% Adjust rotation of the labels
ax.XLabel.Rotation = 30;  % Rotate x-label to match the x-axis rotation
ax.YLabel.Rotation = -5;   % Rotate y-label to match the y-axis rotation




    
    % Combine all data into one table
    % combined_data = cell2mat(export_data);
    % combined_table = array2table(combined_data, 'VariableNames', ...
    %     {'Time', 'TR', 'mRNA_data', 'GdmS_star_data', 'mRNA_fit', 'GdmS_star_fit'});
    % 
    % % Write the table to an Excel file
    % writetable(combined_table, 'globaloutput-dna.csv');
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

% Main script to run the simulation for multiple T7 conditions
function [] = runall3D(T7_concs)
    % Define initial conditions for the ODEs
    initial_mRNAim = 0.0;
    initial_mRNA = 0.0;
    initial_GdmS = 0.0;
    initial_GdmS1 = 0.0;
    initial_GdmS2 = 0.0;
    initial_GdmS3 = 0.0;
    initial_GdmS_star = 0.0;
    initial_TsR = 91.69390911;
    initial_TlR = 905.0532033;

    initial_conditions = [initial_mRNAim, initial_mRNA, initial_GdmS, initial_GdmS1, initial_GdmS2, initial_GdmS3, initial_GdmS_star, initial_TsR, initial_TlR];

    % Define constants [DNA, RNase]
    DNA = 5.6;
    RNase = 270;

    % Define parameters 
    kr = 0.016079516;
    kb1 = 1.22851158;
    kb2 = 0.000100154;
    kp = 0.107670335;
    k1 = 22.91161736;
    k2 = 1.958401957;
    k3 = 1.956331343;
    kmat = 1.957970722;
    kd = 0.004291751;

    params = [kr, kb1, kb2, kp, k1, k2, k3, kmat, kd];

    % Time span for the simulation
    tspan = [0 6];
    num_time_points = 1000; % Number of points in time
    time_points = linspace(tspan(1), tspan(2), num_time_points); % Define the time points array

    % Prepare figure
    figure;
    hold on;

    colors = lines(length(T7_concs)); % Color map for different T7 concentrations

    % Initialize matrices to store results
    T7_grid = [];
    t_grid = [];
    mRNA_grid = [];
    GdmS_star_grid = [];

    for i = 1:length(T7_concs)
        T7 = T7_concs(i);
        constants = [T7, DNA, RNase];

        % Solve the ODE system
        [t, y] = ode45(@(t,y) coupled_odes(t, y, params, constants), time_points, initial_conditions);

        % Store results in grid matrices
        T7_grid = [T7_grid; T7 * ones(size(t))'];
        t_grid = [t_grid; t'];
        mRNA_grid = [mRNA_grid; y(:,2)'];
        GdmS_star_grid = [GdmS_star_grid; y(:,7)'];

        % Plot points for mRNA and GdmS_star
        plot3(t, T7 * ones(size(t)), y(:,2), 'o', 'Color', colors(i,:), 'DisplayName', ['mRNA (T7 = ' num2str(T7) ')']);
        plot3(t, T7 * ones(size(t)), y(:,7), 'o', 'Color', colors(i,:), 'DisplayName', ['GdmS^* (T7 = ' num2str(T7) ')']);
    end

    % Plot lines parallel to the time axis
    for i = 1:length(T7_concs)
        plot3(t_grid(i,:), T7_concs(i) * ones(size(t_grid(i,:))), mRNA_grid(i,:), 'k');
        plot3(t_grid(i,:), T7_concs(i) * ones(size(t_grid(i,:))), GdmS_star_grid(i,:), 'k');
    end

    % Plot lines parallel to the T7 axis
    for j = 1:num_time_points
        for i = 1:length(T7_concs)-1
            % Extract concentrations at specific time points
            mRNA_start = mRNA_grid(i, j);
            GdmS_star_start = GdmS_star_grid(i, j);
            mRNA_end = mRNA_grid(i+1, j);
            GdmS_star_end = GdmS_star_grid(i+1, j);

            % Plot lines
            plot3(time_points(j) * ones(1,2), [T7_concs(i) T7_concs(i+1)], [mRNA_start mRNA_end], 'k--');
            plot3(time_points(j) * ones(1,2), [T7_concs(i) T7_concs(i+1)], [GdmS_star_start GdmS_star_end], 'k--');
        end
    end

    xlabel('Time');
    ylabel('T7 Concentration');
    zlabel('Concentration');
    %title('3D Grid Plot of mRNA and GdmS^* Concentrations');
    %legend('show');
    grid on;
    view(3);
end



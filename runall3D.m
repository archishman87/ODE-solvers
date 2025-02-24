% Main script to run the simulation for multiple T7 conditions
T7_concs = [470, 970, 1470, 1720, 2040];
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

    % Prepare figure
    %figure;
    hold on;

    % Initialize matrices for meshgrid
    T7_grid = [];
    t_grid = [];
    mRNA_grid = [];
    GdmS_star_grid = [];

    for i = 1:length(T7_concs)
        T7 = T7_concs(i);
        constants = [T7, DNA, RNase];

        % Solve the ODE system
        [t, y] = ode45(@(t,y) coupled_odes(t, y, params, constants), linspace(tspan(1), tspan(2), num_time_points), initial_conditions);

        % Store results in grid matrices
        T7_grid = [T7_grid; T7 * ones(size(t))'];
        t_grid = [t_grid; t'];
        mRNA_grid = [mRNA_grid; y(:,2)'];
        GdmS_star_grid = [GdmS_star_grid; y(:,7)'];
    end

    % Convert grids to mesh format
    [T7_mesh, t_mesh] = meshgrid(T7_concs, linspace(tspan(1), tspan(2), num_time_points));

    % Interpolate the data
    mRNA_mesh = griddata(t_grid(:), T7_grid(:), mRNA_grid(:), t_mesh, T7_mesh);
    GdmS_star_mesh = griddata(t_grid(:), T7_grid(:), GdmS_star_grid(:), t_mesh, T7_mesh);

    % Plot mRNA mesh
    surf(t_mesh, T7_mesh, mRNA_mesh,'FaceColor','g', 'EdgeColor', 'none', 'FaceAlpha', 0.2, 'DisplayName', 'mRNA');
    hold on;
    % Plot GdmS_star mesh
    surf(t_mesh, T7_mesh, GdmS_star_mesh,'FaceColor','m' ,'EdgeColor', 'none', 'FaceAlpha', 0.2, 'DisplayName', 'GdmS^*');

    % Customize the plot
    xlabel('Time');
    ylabel('T7 Concentration');
    zlabel('Concentration');
    %title('3D Plot of mRNA and GdmS^* Concentrations');
    %legend('show');
    grid on;
    view(3);

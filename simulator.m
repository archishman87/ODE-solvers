% Load the ODE function and parameters
run('coupled_odes.m');
run('coupled_ode_run.m');

% Define the parameter set
params = [0.016079516, 1.22851158, 0.000100154, 0.107670335, 22.91161736, 1.958401957, 1.956331343, 1.957970722, 0.004291751];

% Define constants [T7, DNA, RNase]
constants = [470, 5.6, 270]; % Example with T7 = 470

% Define time span for the simulation
tspan = [0 6];

% Define the range for initial conditions of TsR and TlR
TsR_range = linspace(50, 200, 10);
TlR_range = linspace(400, 1000, 10);

% Initialize matrices to store mRNA and GdmS_star concentrations
mRNA_conc = zeros(length(TsR_range), length(TlR_range));
GdmS_star_conc = zeros(length(TsR_range), length(TlR_range));

% Loop over each combination of TsR and TlR initial conditions
for i = 1:length(TsR_range)
    for j = 1:length(TlR_range)
        % Initial conditions
        initial_conditions = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, TsR_range(i), TlR_range(j)];
        
        % Solve the ODE
        [t, y] = ode45(@(t,y) coupled_odes(t, y, params, constants), tspan, initial_conditions);
        
        % Store the final concentrations of mRNA and GdmS_star
        mRNA_conc(i, j) = y(end, 2); % mRNA concentration
        GdmS_star_conc(i, j) = y(end, 7); % GdmS_star concentration
    end
end

% Plot the phase diagram
figure;
hold on;
for i = 1:length(TsR_range)
    for j = 1:length(TlR_range)
        plot(mRNA_conc(i, j), GdmS_star_conc(i, j), 'ko', 'MarkerFaceColor', 'k');
    end
end
xlabel('mRNA Concentration');
ylabel('GdmS^* Concentration');
title('Phase Diagram: mRNA vs GdmS^*');
grid on;
hold off;

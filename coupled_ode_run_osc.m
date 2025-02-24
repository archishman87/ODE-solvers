function [] = coupled_ode_run_osc(T7,f)
% Main script
% Define initial conditions for the ODEs
initial_mRNAim = 0.0;
initial_mRNA = 0.0;
initial_GdmS = 0.0;
initial_GdmS1 = 0.0;
initial_GdmS2 = 0.0;
initial_GdmS3 = 0.0;
initial_GdmS_star = 0.0;

% Define constants [T7, DNA, RNase]
DNA = 5.6;
RNase = 270;
constants = [T7, DNA, RNase];

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

params = [kr,kb1,kb2,kp,k1,k2,k3,kmat,kd];

% Time span for the simulation
tspan = [0 6];

% Solve the ODE system
[t, y] = ode45(@(t,y) coupled_odes_osc(t, y, params, constants,f), tspan, [initial_mRNAim, initial_mRNA, initial_GdmS, initial_GdmS1, initial_GdmS2, initial_GdmS3, initial_GdmS_star]);

% Calculate TsR and TlR for each time step
initial_TsR = 324.88736;
initial_TlR = 362.1430036;


decay_rate = 0.5;  % Adjust the decay rate for how fast the oscillation decays
    TsR = initial_TsR * exp(-decay_rate * t) .* (0.5 * sin(2 * pi * t / f) + 0.5);  % Decaying sine wave for TsR
    TlR = initial_TlR * exp(-decay_rate * t).* (0.5 * sin(2 * pi * t / f) + 0.5);  % Decaying sine wave for TlR  % Flat oscillation between 0% and 100% of initial_TlR

% Plot the results
figure;
plot(t, y(:,2)/5, 'g', 'LineWidth', 2, 'DisplayName', 'mRNA');
hold on;
plot(t, (y(:,7))/10e5, 'm', 'LineWidth', 2,'DisplayName', 'GdmS^*');
plot(t, TsR, 'b', 'LineWidth', 2, 'DisplayName', 'TsR');  % Plot oscillating TsR
plot(t, TlR, 'r', 'LineWidth', 2, 'DisplayName', 'TlR');  % Plot oscillating TlR
xlabel('Time');
ylabel('Concentration');
legend;
title('Time Evolution of Species Concentrations');
grid on;
end

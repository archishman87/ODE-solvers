% Define the rate constants
k1 = 1; % example value for k1
k2 = 0.5; % example value for k2
k3 = 0.3; % example value for k3
k4 = 0.2; % example value for k4

kel1 = 0.01; % example value for kel,1
kel2 = 0.05; % example value for kel,2
kel3 = 0.03; % example value for kel,3
kel4 = 0.02; % example value for kel,4

% Define the initial conditions
C0 = 0;     % Initial concentration of Choline (C)
PhC0 = 0;   % Initial concentration of Phosphocholine (PhC)
CDPC0 = 0;  % Initial concentration of CDP-Choline (CDPC)
PtC0 = 0;   % Initial concentration of Phosphatidylcholine (PtC)

initial_conditions = [C0, PhC0, CDPC0, PtC0];

% Define the time span for the simulation
tspan = [0 10]; % example time span from 0 to 100



% Solve the ODEs using ode45
[t, Y] = ode45(@(t, Y) odes(t, Y, k1, k2, k3, k4, kel1, kel2, kel3, kel4), tspan, initial_conditions);

% Plot the results
figure;
plot(t, Y(:,1), 'r', 'LineWidth', 2); hold on;
plot(t, Y(:,2), 'g', 'LineWidth', 2);
plot(t, Y(:,3), 'b', 'LineWidth', 2);
plot(t, Y(:,4), 'm', 'LineWidth', 2);
hold off;

xlabel('Time');
ylabel('Concentration');
legend('C (Choline)', 'PhC (Phosphocholine)', 'CDPC (CDP-Choline)', 'PtC (Phosphatidylcholine)');
title('Concentration of Components Over Time');
grid on;

% Define the system of ODEs
function dYdt = odes(t, Y, k1, k2, k3, k4, kel1, kel2, kel3, kel4)
    C = Y(1);
    PhC = Y(2);
    CDPC = Y(3);
    PtC = Y(4);
    
    dCdt = k1 - k2*C - kel1*C;
    dPhCdt = k2*C - k3*PhC - kel2*PhC;
    dCDPCdt = k3*PhC - k4*CDPC - kel3*CDPC;
    dPtCdt = k4*CDPC - kel4*PtC;
    
    dYdt = [dCdt; dPhCdt; dCDPCdt; dPtCdt];
end


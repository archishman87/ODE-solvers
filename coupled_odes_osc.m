% Define the system of differential equations
function dydt = coupled_odes_osc(t, y, params, constants,f)
    % Unpack parameters
    kr = params(1);
    kb1 = params(2);
    kb2 = params(3);
    kp = params(4);
    k1 = params(5);
    k2 = params(6);
    k3 = params(7);
    kmat = params(8);
    kd = params(9);
    kdp = 0.01;  % Degradation rate constant for GdmS_star

    % Unpack constants
    T7 = constants(1);
    DNA = constants(2);
    RNase = constants(3);
    
    % Unpack variables
    mRNAim = y(1);
    mRNA = y(2);
    GdmS = y(3);
    GdmS1 = y(4);
    GdmS2 = y(5);
    GdmS3 = y(6);
    GdmS_star = y(7);
    
    % Initial values for TsR and TlR
    initial_TsR = 324.88736;
    initial_TlR = 362.1430036;
    
    % Oscillating TsR and TlR between 0% and 100% of their initial values
    % The sine wave oscillates between 0 and 1, and is scaled to the initial value
    decay_rate = 0.5;  % Adjust the decay rate for how fast the oscillation decays
    TsR = initial_TsR * exp(-decay_rate * t).* (0.5 * sin(2 * pi * t / f) + 0.5);  % Decaying sine wave for TsR
    TlR = initial_TlR * exp(-decay_rate * t).* (0.5 * sin(2 * pi * t / f) + 0.5);  % Decaying sine wave for TlR  % Flat oscillation between 0% and 100% of initial_TlR
    
    % Differential equations
    dydt = zeros(7,1);  % Adjusted to match the number of variables
    dydt(1) = kr * DNA * TsR * T7 - kb1 * mRNAim + kb2 * mRNA;
    dydt(2) = -kd * mRNA * RNase + kb1 * mRNAim - kb2 * mRNA;
    dydt(3) = kp * mRNA * TlR * TsR / 2 - kmat * GdmS;
    dydt(4) = k1 * GdmS - k2 * GdmS1;
    dydt(5) = k2 * GdmS1 - k3 * GdmS2;
    dydt(6) = k3 * GdmS2 - kmat * GdmS3;
    dydt(7) = kmat * GdmS3 - kdp * GdmS_star;  % Degradation of GdmS_star
end

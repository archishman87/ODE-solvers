% Define the system of differential equations
function dydt = coupled_odes(t, y, params, constants)
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
    TsR = y(8);
    TlR = y(9);
    
    % Differential equations
    dydt = zeros(9,1);
    dydt(1) = kr * DNA * TsR * T7 - kb1 * mRNAim + kb2 * mRNA;
    dydt(2) = -kd * mRNA * RNase + kb1 * mRNAim - kb2 * mRNA;
    dydt(3) = kp * mRNA * TlR * TsR / 2 - kmat * GdmS;
    dydt(4) = k1 * GdmS - k2 * GdmS1;
    dydt(5) = k2 * GdmS1 - k3 * GdmS2;
    dydt(6) = k3 * GdmS2 - kmat * GdmS3;
    dydt(7) = kmat * GdmS3;
    dydt(8) = -kr * DNA * TsR * T7 - kp * mRNA * TlR * TsR / 2;
    dydt(9) = -kp * mRNA * TlR * TsR / 2;
end







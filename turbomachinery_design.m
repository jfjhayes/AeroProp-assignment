%% ENG5313: Aerospace Propulsion M Coursework 
% Setion 3 - Turbomachinery Design %
clc;
clf;
clear;

% Constants %
gamma = 1.4;
P_std = 101.325;                            % Pressure standard day (kN/m^2)
T_std = 288.15;                             % Temp standard day (K)
rho_SL = 1.225;                             % Density @ SL (kg/m^3)
g = 9.81;                                   % Standard acceleration of gravity (m/s^2)
R = 287;                             % Specific gas constant of dry air (kJ/kgK)
C_P = 1.005*1000;                                % Specific heat capacity of dry air (kJ/kgK)
Q_r = 43100;                                % Calorific value of fuel (kJ/kgK)

% Post-combustion constants %
gamma_g = 1.333;
C_Pg = 1.148*1000;                               % Specific heat capacity of HOT air (kJ/kgK)


% What do you want to design ? %
compressor = true;
turbine = false;

if compressor

    N_IP = 158.3;             % Change based on turbine results?
    N_HP = 175;

    % Constraints %
    M_rel1 = 1.21;
    DH_IP = 0.77;
    DH_HP = 0.8;

    % Get values from section 2 %
    % [T_o23, P_o23, T_o25, P_o25, T_o3, P_o3, mdot_core, eta_IPC, eta_HPC, Ma_2] = cycle_analyser3(B, pi_cH, pi_cI, T_o4, pi_fC)
    [T_oIP1, P_oIP1, T_oHP1, P_oHP1, T_oHP2, P_oHP2, mdot_core, eta_IPC, eta_HPC, ~] = cycle_analyser3(10.48, 10, 2.5, 2090, 2.27);

    % Initialise % 
    T_oIP2 = T_oHP1;
    P_oIP2 = P_oHP1;

    Delta_T_IP = T_oIP2 - T_oIP1;
    Delta_T_HP = T_oHP2 - T_oHP1;

    Ma_2 = 0.6;

    % No inlet guide vanes % 
    C_1 = Ma_2 * sqrt((gamma * R *T_oIP1)/(1+(gamma-1)/2 * Ma_2^2));
    C_a = C_1;


    % IPC static temps, pressures, and densities % 
    T_sIP1 = T_oIP1 - (C_a^2 / (2 * C_P));
    T_sIP2 = T_oIP2 - C_a^2 / (2 * C_P); 
    P_sIP1 = P_oIP1 * (T_sIP1 / T_oIP1)^(gamma/(gamma-1)); 
    P_sIP2 = P_oIP2 * (T_sIP2 / T_oIP2)^(gamma/(gamma-1)); 
    rho_IP1 = P_sIP1 / (R/1000 * T_sIP1);
    rho_IP2 = P_sIP2 / (R/1000 * T_sIP2); 

    a_IP1 = sqrt(gamma * R * T_sIP1);                           % local speed of sound
    V_IP1t = a_IP1 * M_rel1;                                    % relative tip velocity
    U_IP1t = (V_IP1t^2 - C_a^2)^0.5;                             % blade tip velocity

    r_IP1t = U_IP1t / (2 * pi * N_IP);                                          % inlet tip radius
    HTR_ratio_IP1 = sqrt(1 - (mdot_core / (r_IP1t^2 * rho_IP1 * pi * C_a)));    % HTR ratio (from continuity

    r_IP1h = r_IP1t * HTR_ratio_IP1;                                            % inlet hub height
    r_IPm = 0.5 * (r_IP1h + r_IP1t);                                            % mean radius
    h_IP1 = r_IP1t - r_IP1h;                                                    % inlet blade height
    A_IP1 = mdot_core / (rho_IP1 * C_a);                                        % inlet annulus area
    A_IP2 = mdot_core / (rho_IP2 * C_a);                                        % outlet annulus area
    h_IP2 = A_IP2 / (2 * pi * r_IPm);                                           % outlet blade height (constant mean r)
    r_IP2t = r_IPm + (h_IP2 / 2);                                               % outlet tip radius
    r_IP2h = r_IPm - (h_IP2 / 2);                                               % outlet hub radius

    U_m_IP = 2 * pi * r_IPm * N_IP;   
    beta_IP1m = atand(U_m_IP / C_a);
    V_IP1m = C_a / cosd(beta_IP1m);
    V_IP2m = DH_IP * V_IP1m;                                        % Use DH criterion              
    beta_IP2m = acosd(C_a / V_IP2m);
    Delta_T_IP_stage = (C_a * U_m_IP * (tand(beta_IP1m) - tand(beta_IP2m)))/C_P;
    N_stage_IP = Delta_T_IP / Delta_T_IP_stage + 1;

    % HPC static temps, pressures, and densities % 
    T_sHP1 = T_oHP1 - C_a^2 / (2 * C_P); 
    T_sHP2 = T_oHP2 - C_a^2 / (2 * C_P); 
    P_sHP1 = P_oHP1 * (T_sHP1 / T_oHP1)^(gamma/(gamma-1)); 
    P_sHP2 = P_oHP2 * (T_sHP2 / T_oHP2)^(gamma/(gamma-1)); 
    rho_HP1 = P_sHP1 / (R/1000 * T_sHP1); 
    rho_HP2 = P_sHP2 / (R/1000 * T_sHP2); 

    a_HP1 = sqrt(gamma * R * T_sHP1);                           % local speed of sound
    V_HP1t = a_HP1 * M_rel1;                                    % relative tip velocity
    U_HP1t = (V_HP1t^2 - C_a^2)^0.5;                             % blade tip velocity

    r_HP1t = U_HP1t / (2 * pi * N_HP);                                          % inlet tip radius
    HTR_ratio_HP1 = sqrt(1 - (mdot_core / (r_HP1t^2 * rho_HP1 * pi * C_a)));    % HTR ratio (from continuity)
    r_HP1h = r_HP1t * HTR_ratio_HP1;                                            % inlet hub height
    r_HPm = 0.5 * (r_HP1h + r_HP1t);                                            % mean radius
    h_HP1 = r_HP1t - r_HP1h;                                                    % inlet blade height
    A_HP1 = mdot_core / (rho_HP1 * C_a);                                        % inlet annulus area
    A_HP2 = mdot_core / (rho_HP2 * C_a);                                        % outlet annulus area
    h_HP2 = A_HP2 / (2 * pi * r_HPm);                                           % outlet blade height (constant mean r)
    r_HP2t = r_HPm + (h_HP2 / 2);                                               % outlet tip radius
    r_HP2h = r_HPm - (h_HP2 / 2);                                               % outlet hub radius

    U_m_HP = 2 * pi * r_HPm * N_HP;   
    beta_HP1m = atand(U_m_HP / C_a);
    V_HP1m = C_a / cosd(beta_HP1m);
    V_HP2m = DH_HP * V_HP1m;                                        % Use DH criterion              
    beta_HP2m = acosd(C_a / V_HP2m);
    Delta_T_HP_stage = (C_a * U_m_HP * (tand(beta_HP1m) - tand(beta_HP2m)))/C_P;
    N_stage_HP = Delta_T_HP / Delta_T_HP_stage + 1;

    fprintf('Hub-tip radius ratio at IP compressor inlet: %.4f\n', HTR_ratio_IP1);
    fprintf('Annulus hub radius at IP compressor inlet: %.4f m\n', r_IP1h);
    fprintf('Annulus tip radius at IP compressor inlet: %.4f m\n', r_IP1t);
    fprintf('Annulus hub radius at IP compressor exit: %.4f m\n', r_IP2h);
    fprintf('Annulus tip radius at IP compressor exit: %.4f m\n', r_IP2t);
    fprintf('Blade height at IP compressor exit: %.4f m\n', h_IP2);
    fprintf('Number of compressor stages required for IP compressor: %.0f\n', N_stage_IP);
    fprintf('IP Shaft speed: %.2f rpm\n', N_IP*60);

    fprintf('Hub-tip radius ratio at HP compressor inlet: %.4f\n', HTR_ratio_HP1);
    fprintf('Annulus hub radius at HP compressor inlet: %.4f m\n', r_HP1h);
    fprintf('Annulus tip radius at HP compressor inlet: %.4f m\n', r_HP1t);
    fprintf('Annulus hub radius at HP compressor exit: %.4f m\n', r_HP2h);
    fprintf('Annulus tip radius at HP compressor exit: %.4f m\n', r_HP2t);
    fprintf('Blade height at HP compressor exit: %.4f m\n', h_HP2);
    fprintf('Number of compressor stages required for HP compressor: %.0f\n', N_stage_HP);
    fprintf('HP Shaft speed: %.2f rpm\n', N_HP*60);


end

if turbine
    
    % Constraints % 
    psi_limit = 2.695;          % Stage loading limit (toleranced)

    N_IP = 158.3;             % Change based on turbine results?
    N_HP = 175;

    alpha_1_HPT = 0;

    % From section 2 %
    eta_HPT = 0.9367;
    eta_IPT = 0.9318;
    eta_LPT = 0.9396;

    % Critical Pressure Ratio % 
    P_crit = ((gamma_g + 1) / 2)^((gamma_g) / (gamma_g - 1));

    % Stuff to Vary %
    r_HPTm = 0.4;              % Allowed to be between 0.3 and 0.4 (prob tolerance later)
    r_IPTm = 0.425;
    U_HPT = r_HPTm * (2 * pi * N_HP);               % Mean blade speed?
    U_IPT = r_IPTm * (2 * pi * N_IP);               % Mean blade speed?


    % T_o4, P_o4, eta_HPT, T_o45, P_o45, T_o47, P_o47
    % 4 - HPT Entry
    % 45 - IPT Entry
    % 47 - LPT Entry

    [T_o4, P_o4, eta_HPT, T_o45, P_o45, T_o47, P_o47, mdot_core] = cycle_analyser4(10.48, 10, 2.5, 2090, 2.27);

    Delta_T = T_o4 - T_o47;             % temp drop over entire turbine
    pi_T = P_o4 / P_o47;                % pressure drop over entire turbine

    % HPT %
    psi_HPT = (C_Pg * (T_o4 - T_o45)) / ((2 * pi * r_HPTm * N_HP)^2);       % HPT stage loading (should be less than 2.695)
    alpha_3_HPT = deg2rad(30);
    phi_HPT = 0.6;

    beta_3_HPT = atan(1 / phi_HPT + tan(alpha_3_HPT));
    lambda_HPT = 0.5 * ((2 * phi_HPT * tan(beta_3_HPT)) - psi_HPT);         % reaction

    beta_2_HPT = atan(tan(beta_3_HPT) - (2 * lambda_HPT)/phi_HPT);
    alpha_2_HPT = atan((1 / phi_HPT) + tan(beta_2_HPT));

    C_a2_HPT = U_HPT * phi_HPT; 
    C_2_HPT = C_a2_HPT / cos(alpha_2_HPT);
    T_2_HPT = T_o4 - ((C_2_HPT)^2) / (2 * C_Pg);

    a_HPT = sqrt(gamma_g * R * T_2_HPT);         % local speed of sound
    P_ratio_HPT = P_o4 / P_o45;
    P_2_HPT = P_o4 / P_ratio_HPT;

    if P_ratio_HPT < P_crit
        P_exit_HPT = P_2_HPT;
    else
        P_exit_HPT = P_o4 / P_crit;
    end

    % Interstage annulus
    A_2_HPT = (R/1000 * T_2_HPT * mdot_core) / (P_exit_HPT * C_a2_HPT);
    C_a3_HPT = C_a2_HPT;
    C_a1_HPT = C_a3_HPT / (cos(alpha_3_HPT));

    % Turbine inlet static properties
    T_1_HPT = T_o4 - (((C_a1_HPT)^2) / (2 * C_Pg));
    P_1_HPT = P_o4 * (T_1_HPT / T_o4)^((gamma_g)/(gamma_g-1));
    A_1_HPT = (mdot_core * R/1000 * T_1_HPT) / (P_1_HPT * C_a1_HPT);

    % Outlet
    T_o3_HPT = T_o45; 
    C_3_HPT = C_a1_HPT;
    T_3_HPT = T_o3_HPT - ((C_3_HPT^2) / (2 * C_Pg));
    P_o3_HPT = P_o45;
    P_3_HPT = P_o3_HPT * (T_3_HPT / T_o3_HPT)^((gamma_g)/(gamma_g-1));
    A_3_HPT = (R/1000 * T_3_HPT * mdot_core) / (P_3_HPT * C_a3_HPT);

    % Radii % 
    h1_HPT = (N_HP / U_HPT) * A_1_HPT;
    rt1_HPT = r_HPTm + (h1_HPT / 2);
    rr1_HPT = r_HPTm - (h1_HPT / 2);
    HTR_ratio_HPT1 = rt1_HPT / rr1_HPT;

    h2_HPT = (N_HP / U_HPT) * A_2_HPT;
    rt2_HPT = r_HPTm + (h2_HPT / 2);
    rr2_HPT = r_HPTm - (h2_HPT / 2);
    HTR_ratio_HPT2 = rt2_HPT / rr2_HPT;

    h3_HPT = (N_HP / U_HPT) * A_3_HPT;
    rt3_HPT = r_HPTm + (h3_HPT / 2);
    rr3_HPT = r_HPTm - (h3_HPT / 2);
    HTR_ratio_HPT3 = rt3_HPT / rr3_HPT;

    % Angles %
    alpha_1t_HPT = atan((r_HPTm/rt1_HPT) * tan(alpha_1_HPT));
    alpha_1r_HPT = atan((r_HPTm/rr1_HPT) * tan(alpha_1_HPT));

    alpha_2t_HPT = atan((r_HPTm/rt2_HPT) * tan(alpha_2_HPT));
    alpha_2r_HPT = atan((r_HPTm/rr2_HPT) * tan(alpha_2_HPT));

    alpha_3t_HPT = atan((r_HPTm/rt3_HPT) * tan(alpha_3_HPT));
    alpha_3r_HPT = atan((r_HPTm/rr3_HPT) * tan(alpha_3_HPT));

    b2rt1_HPT = (r_HPTm / rt2_HPT) * tan(alpha_2_HPT);
    b2rt2_HPT = (rt2_HPT / r_HPTm) * (U_HPT / C_a2_HPT);
    beta_2t_HPT = atan(b2rt1_HPT - b2rt2_HPT);

    b2rr1_HPT = (r_HPTm / rr2_HPT) * tan(alpha_2_HPT);
    b2rr2_HPT = (rr2_HPT / r_HPTm) * (U_HPT / C_a2_HPT);
    beta_2r_HPT = atan(b2rr1_HPT - b2rr2_HPT);

    b3rt1_HPT = (r_HPTm / rt3_HPT) * tan(alpha_3_HPT);
    b3rt2_HPT = (rt3_HPT / r_HPTm) * (U_HPT / C_a3_HPT);
    beta_3t_HPT = atan(b3rt1_HPT + b3rt2_HPT);

    b3rr1_HPT = (r_HPTm / rr3_HPT) * tan(alpha_3_HPT);
    b3rr2_HPT = (rr3_HPT / r_HPTm) * (U_HPT / C_a3_HPT);
    beta_3r_HPT = atan(b3rr1_HPT + b3rr2_HPT);

    phi_t_HPT = 1 / (tan(beta_3t_HPT) - tan(alpha_3t_HPT));
    phi_r_HPT = 1 / (tan(beta_3r_HPT) - tan(alpha_3r_HPT));

    lambda_ts1_HPT = (phi_t_HPT/2)* (tan(beta_3t_HPT) - tan(beta_2t_HPT));
    lambda_rs1_HPT = (phi_r_HPT/2)* (tan(beta_3r_HPT) - tan(beta_2r_HPT));

    % Blade speeds %
    U_r2s1_HPT = rr2_HPT * (2 * pi * N_HP);
    U_t2s1_HPT = rt2_HPT * (2 * pi * N_HP);
    U_r3s1_HPT = rr3_HPT * (2 * pi * N_HP);
    U_t3s1_HPT = rt3_HPT * (2 * pi * N_HP);

    % Miscellaneous % 
    C_t1s1_HPT = C_a1_HPT / cos(alpha_1t_HPT);
    C_r1s1_HPT = C_a1_HPT / cos(alpha_1r_HPT);
    C_1s1_HPT = C_a1_HPT / cos(alpha_1_HPT);

    C_t2s1_HPT = C_a2_HPT / cos(alpha_2t_HPT);
    C_r2s1_HPT = C_a2_HPT / cos(alpha_2r_HPT);
    C_2s1_HPT = C_a2_HPT / cos(alpha_2_HPT);

    V_r2s1_HPT = C_a2_HPT / cos(beta_2r_HPT);
    V_t2s1_HPT = C_a2_HPT / cos(beta_2r_HPT);
    V_2s1_HPT = C_a2_HPT / cos(beta_2_HPT);

    C_t3s1_HPT = C_a3_HPT / cos(alpha_3t_HPT);
    C_r3s1_HPT = C_a3_HPT / cos(alpha_3r_HPT);
    C_3s1_HPT = C_a3_HPT / cos(alpha_3_HPT);

    V_r3s1_HPT = C_a3_HPT / cos(beta_3r_HPT);
    V_t3s1_HPT = C_a3_HPT / cos(beta_3r_HPT);
    V_3s1_HPT = C_a3_HPT / cos(beta_3_HPT);

    fprintf('\n*** HIGH-PRESSURE (HP) TURBINE ***\n\n');

    % Inlet/Exit Geometry and Flow Parameters
    fprintf('Hub-tip radius ratio at HP turbine inlet (-): %f\n', HTR_ratio_HPT1);
    fprintf('Annulus hub radius at HP turbine inlet (m): %f\n', rr1_HPT);
    fprintf('Annulus tip radius at HP turbine inlet (m): %f\n', rt1_HPT);
    fprintf('Annulus hub radius at HP turbine exit (m): %f\n', rr3_HPT);
    fprintf('Annulus tip radius at HP turbine exit (m): %f\n', rt3_HPT);
    fprintf('Blade height at HP turbine inlet (m): %f\n', h1_HPT);
    fprintf('Blade height at HP turbine exit (m): %f\n', h3_HPT);
    fprintf('Mid (mean) height radius for HP turbine (m): %f\n', r_HPTm);
    fprintf('Stage loading, psi: %f\n', psi_HPT);
    fprintf('Flow coefficient, phi: %f\n', phi_HPT);
    fprintf('Mid-height stage axial velocity (m/s): %f\n', C_a2_HPT);
    fprintf('\n*** HP TURBINE - REACTION & VELOCITIES (Hub, Mid, Tip) ***\n\n');
    fprintf('Reaction, lambda (deg): Hub = %f, Mid = %f, Tip = %f\n', ...
        rad2deg(lambda_rs1_HPT), rad2deg((lambda_rs1_HPT + lambda_ts1_HPT)/2), rad2deg(lambda_ts1_HPT));
    fprintf('Blade speed at turbine rotor inlet U (m/s): Hub = %f, Mid = %f, Tip = %f\n', ...
        U_r2s1_HPT, U_HPT, U_t2s1_HPT);
    fprintf('Blade speed at turbine rotor exit U (m/s): Hub = %f, Mid = %f, Tip = %f\n', ...
        U_r3s1_HPT, U_HPT, U_t3s1_HPT);
    fprintf('Inlet angle, alpha1 (deg): Hub = %f, Mid = %f, Tip = %f\n', ...
        rad2deg(alpha_1r_HPT), rad2deg(alpha_1_HPT), rad2deg(alpha_1t_HPT));
    fprintf('Nozzle exit angle, alpha2 (deg): Hub = %f, Mid = %f, Tip = %f\n', ...
        rad2deg(alpha_2r_HPT), rad2deg(alpha_2_HPT), rad2deg(alpha_2t_HPT));
    fprintf('Relative nozzle exit angle / blade inlet angle, beta2 (deg): Hub = %f, Mid = %f, Tip = %f\n', ...
        rad2deg(beta_2r_HPT), rad2deg(beta_2_HPT), rad2deg(beta_2t_HPT));
    fprintf('Relative blade exit angle, beta3 (deg): Hub = %f, Mid = %f, Tip = %f\n', ...
        rad2deg(beta_3r_HPT), rad2deg(beta_3_HPT), rad2deg(beta_3t_HPT));
    fprintf('Absolute blade exit angle, alpha3 (deg): Hub = %f, Mid = %f, Tip = %f\n', ...
        rad2deg(alpha_3r_HPT), rad2deg(alpha_3_HPT), rad2deg(alpha_3t_HPT));
    fprintf('Nozzle Inlet velocity, C1 (m/s): Hub = %f, Mid = %f, Tip = %f\n', ...
        C_r1s1_HPT, C_1s1_HPT, C_t1s1_HPT);
    fprintf('Nozzle exit velocity, C2 (m/s): Hub = %f, Mid = %f, Tip = %f\n', ...
        C_r2s1_HPT, C_2s1_HPT, C_t2s1_HPT);
    fprintf('Relative nozzle exit velocity / blade inlet velocity, V2r (m/s): Hub = %f, Mid = %f, Tip = %f\n', ...
        V_r2s1_HPT, V_2s1_HPT, V_t2s1_HPT);
    fprintf('Relative blade exit velocity, V3r (m/s): Hub = %f, Mid = %f, Tip = %f\n', ...
        V_r3s1_HPT, V_3s1_HPT, V_t3s1_HPT);
    fprintf('Absolute blade exit velocity, C3r (m/s): Hub = %f, Mid = %f, Tip = %f\n\n', ...
        C_r3s1_HPT, C_3s1_HPT, C_t3s1_HPT);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % IPT - NEED TO FIGURE OUT VALUES WHAT GET PASSED BETWEEN HPT AND IPT % 
    psi_IPT = (C_Pg * (T_o45 - T_o47)) / ((2 * pi * r_IPTm * N_IP)^2);      % IPT stage loading (shoudl be less than 2.695)
    phi_IPT = 0.6;

    alpha_1_IPT = alpha_3_HPT;
    alpha_3_IPT = deg2rad(17.6);            % Exit angle limit (toleranced)


    beta_3_IPT = atan(1 / phi_IPT + tan(alpha_3_IPT));
    lambda_IPT = 0.5 * ((2 * phi_IPT * tan(beta_3_IPT)) - psi_IPT);         % reaction

    beta_2_IPT = atan(tan(beta_3_IPT) - (2 * lambda_IPT)/phi_IPT);
    alpha_2_IPT = atan((1 / phi_IPT) + tan(beta_2_IPT));

    C_a2_IPT = U_IPT * phi_IPT; 
    C_2_IPT = C_a2_IPT / cos(alpha_2_IPT);
    T_2_IPT = T_o45 - ((C_2_IPT)^2) / (2 * C_Pg);

    a_IPT = sqrt(gamma_g * R * T_2_IPT);         % local speed of sound
    P_ratio_IPT = P_o45 / P_o47;
    P_2_IPT = P_o47 / P_ratio_IPT;

    if P_ratio_IPT < P_crit
        P_exit_IPT = P_2_IPT;
    else
        P_exit_IPT = P_o45 / P_crit;
    end

    % Interstage annulus
    A_2_IPT = (R/1000 * T_2_IPT * mdot_core) / (P_exit_IPT * C_a2_IPT);
    C_a3_IPT = C_a2_IPT;
    C_a1_IPT = C_a3_IPT / (cos(alpha_3_IPT));

    % Turbine inlet static properties
    T_1_IPT = T_o45 - (((C_a1_IPT)^2) / (2 * C_Pg));
    P_1_IPT = P_o45 * (T_1_IPT / T_o45)^((gamma_g)/(gamma_g-1));
    A_1_IPT = (mdot_core * R/1000 * T_1_IPT) / (P_1_IPT * C_a1_IPT);

    % Outlet
    T_o3_IPT = T_o47; 
    C_3_IPT = C_a1_IPT;
    T_3_IPT = T_o3_IPT - ((C_3_IPT^2) / (2 * C_Pg));
    P_o3_IPT = P_o47;
    P_3_IPT = P_o3_IPT * (T_3_IPT / T_o3_IPT)^((gamma_g)/(gamma_g-1));
    A_3_IPT = (R/1000 * T_3_IPT * mdot_core) / (P_3_IPT * C_a3_IPT);

    % Radii % 
    h1_IPT = (N_IP / U_IPT) * A_1_IPT;
    rt1_IPT = r_IPTm + (h1_IPT / 2);
    rr1_IPT = r_IPTm - (h1_IPT / 2);
    HTR_ratio_IPT1 = rt1_IPT / rr1_IPT;

    h2_IPT = (N_IP / U_IPT) * A_2_IPT;
    rt2_IPT = r_IPTm + (h2_IPT / 2);
    rr2_IPT = r_IPTm - (h2_IPT / 2);
    HTR_ratio_IPT2 = rt2_IPT / rr2_IPT;

    h3_IPT = (N_IP / U_IPT) * A_3_IPT;
    rt3_IPT = r_IPTm + (h3_IPT / 2);
    rr3_IPT = r_IPTm - (h3_IPT / 2);
    HTR_ratio_IPT3 = rt3_IPT / rr3_IPT;

    % Angles %
    alpha_1t_IPT = atan((r_IPTm/rt1_IPT) * tan(alpha_1_IPT));
    alpha_1r_IPT = atan((r_IPTm/rr1_IPT) * tan(alpha_1_IPT));

    alpha_2t_IPT = atan((r_IPTm/rt2_IPT) * tan(alpha_2_IPT));
    alpha_2r_IPT = atan((r_IPTm/rr2_IPT) * tan(alpha_2_IPT));

    alpha_3t_IPT = atan((r_IPTm/rt3_IPT) * tan(alpha_3_IPT));
    alpha_3r_IPT = atan((r_IPTm/rr3_IPT) * tan(alpha_3_IPT));

    b2rt1_IPT = (r_IPTm / rt2_IPT) * tan(alpha_2_IPT);
    b2rt2_IPT = (rt2_IPT / r_IPTm) * (U_IPT / C_a2_IPT);
    beta_2t_IPT = atan(b2rt1_IPT - b2rt2_IPT);

    b2rr1_IPT = (r_IPTm / rr2_IPT) * tan(alpha_2_IPT);
    b2rr2_IPT = (rr2_IPT / r_IPTm) * (U_IPT / C_a2_IPT);
    beta_2r_IPT = atan(b2rr1_IPT - b2rr2_IPT);

    b3rt1_IPT = (r_IPTm / rt3_IPT) * tan(alpha_3_IPT);
    b3rt2_IPT = (rt3_IPT / r_IPTm) * (U_IPT / C_a3_IPT);
    beta_3t_IPT = atan(b3rt1_IPT + b3rt2_IPT);

    b3rr1_IPT = (r_IPTm / rr3_IPT) * tan(alpha_3_IPT);
    b3rr2_IPT = (rr3_IPT / r_IPTm) * (U_IPT / C_a3_IPT);
    beta_3r_IPT = atan(b3rr1_IPT + b3rr2_IPT);

    phi_t_IPT = 1 / (tan(beta_3t_IPT) - tan(alpha_3t_IPT));
    phi_r_IPT = 1 / (tan(beta_3r_IPT) - tan(alpha_3r_IPT));

    lambda_ts1_IPT = (phi_t_IPT/2)* (tan(beta_3t_IPT) - tan(beta_2t_IPT));
    lambda_rs1_IPT = (phi_r_IPT/2)* (tan(beta_3r_IPT) - tan(beta_2r_IPT));

    % Blade speeds %
    U_r2s1_IPT = rr2_IPT * (2 * pi * N_IP);
    U_t2s1_IPT = rt2_IPT * (2 * pi * N_IP);
    U_r3s1_IPT = rr3_IPT * (2 * pi * N_IP);
    U_t3s1_IPT = rt3_IPT * (2 * pi * N_IP);

    % Miscellaneous % 
    C_t1s1_IPT = C_a1_IPT / cos(alpha_1t_IPT);
    C_r1s1_IPT = C_a1_IPT / cos(alpha_1r_IPT);
    C_1s1_IPT = C_a1_IPT / cos(alpha_1_IPT);

    C_t2s1_IPT = C_a2_IPT / cos(alpha_2t_IPT);
    C_r2s1_IPT = C_a2_IPT / cos(alpha_2r_IPT);
    C_2s1_IPT = C_a2_IPT / cos(alpha_2_IPT);

    V_r2s1_IPT = C_a2_IPT / cos(beta_2r_IPT);
    V_t2s1_IPT = C_a2_IPT / cos(beta_2r_IPT);
    V_2s1_IPT = C_a2_IPT / cos(beta_2_IPT);

    C_t3s1_IPT = C_a3_IPT / cos(alpha_3t_IPT);
    C_r3s1_IPT = C_a3_IPT / cos(alpha_3r_IPT);
    C_3s1_IPT = C_a3_IPT / cos(alpha_3_IPT);

    V_r3s1_IPT = C_a3_IPT / cos(beta_3r_IPT);
    V_t3s1_IPT = C_a3_IPT / cos(beta_3r_IPT);
    V_3s1_IPT = C_a3_IPT / cos(beta_3_IPT);

    fprintf('\n*** INTERMEDIATE-PRESSURE (IP) TURBINE ***\n\n');
    fprintf('Hub-tip radius ratio at IP turbine inlet (-): %f\n', HTR_ratio_IPT1);
    fprintf('Annulus hub radius at IP turbine inlet (m): %f\n', rr1_IPT);
    fprintf('Annulus tip radius at IP turbine inlet (m): %f\n', rt1_IPT);
    fprintf('Annulus hub radius at IP turbine exit (m): %f\n', rr3_IPT);
    fprintf('Annulus tip radius at IP turbine exit (m): %f\n', rt3_IPT);
    fprintf('Blade height at IP turbine inlet (m): %f\n', h1_IPT);
    fprintf('Blade height at IP turbine exit (m): %f\n', h3_IPT);
    fprintf('Mid (mean) height radius for IP turbine (m): %f\n', r_IPTm);
    fprintf('Stage loading, psi: %f\n', psi_IPT);
    fprintf('Flow coefficient, phi: %f\n', phi_IPT);
    fprintf('Mid-height axial velocity (m/s): %f\n', C_a2_IPT);
    fprintf('Shaft speed (rpm): %f\n', N_IP);
    fprintf('\n*** IP TURBINE - REACTION & VELOCITIES (Hub, Mid, Tip) ***\n\n');
    fprintf('Reaction, lambda (deg): Hub = %f, Mid = %f, Tip = %f\n', ...
        rad2deg(lambda_rs1_IPT), rad2deg((lambda_rs1_IPT + lambda_ts1_IPT)/2), rad2deg(lambda_ts1_IPT));
    fprintf('Blade speed at turbine rotor inlet U (m/s): Hub = %f, Mid = %f, Tip = %f\n', ...
        U_r2s1_IPT, U_IPT, U_t2s1_IPT);
    fprintf('Blade speed at turbine rotor exit U (m/s): Hub = %f, Mid = %f, Tip = %f\n', ...
        U_r3s1_IPT, U_IPT, U_t3s1_IPT);
    fprintf('Inlet angle, alpha1 (deg): Hub = %f, Mid = %f, Tip = %f\n', ...
        rad2deg(alpha_1r_IPT), rad2deg(alpha_1_IPT), rad2deg(alpha_1t_IPT));
    fprintf('Nozzle exit angle, alpha2 (deg): Hub = %f, Mid = %f, Tip = %f\n', ...
        rad2deg(alpha_2r_IPT), rad2deg(alpha_2_IPT), rad2deg(alpha_2t_IPT));
    fprintf('Relative nozzle exit angle / blade inlet angle, beta2 (deg): Hub = %f, Mid = %f, Tip = %f\n', ...
        rad2deg(beta_2r_IPT), rad2deg(beta_2_IPT), rad2deg(beta_2t_IPT));
    fprintf('Relative blade exit angle, beta3 (deg): Hub = %f, Mid = %f, Tip = %f\n', ...
        rad2deg(beta_3r_IPT), rad2deg(beta_3_IPT), rad2deg(beta_3t_IPT));
    fprintf('Absolute blade exit angle, alpha3 (deg): Hub = %f, Mid = %f, Tip = %f\n', ...
        rad2deg(alpha_3r_IPT), rad2deg(alpha_3_IPT), rad2deg(alpha_3t_IPT));
    fprintf('Nozzle Inlet velocity, C1 (m/s): Hub = %f, Mid = %f, Tip = %f\n', ...
        C_r1s1_IPT, C_1s1_IPT, C_t1s1_IPT);
    fprintf('Nozzle exit velocity, C2 (m/s): Hub = %f, Mid = %f, Tip = %f\n', ...
        C_r2s1_IPT, C_2s1_IPT, C_t2s1_IPT);
    fprintf('Relative nozzle exit velocity / blade inlet velocity, V2r (m/s): Hub = %f, Mid = %f, Tip = %f\n', ...
        V_r2s1_IPT, V_2s1_IPT, V_t2s1_IPT);
    fprintf('Relative blade exit velocity, V3r (m/s): Hub = %f, Mid = %f, Tip = %f\n', ...
        V_r3s1_IPT, V_3s1_IPT, V_t3s1_IPT);
    fprintf('Absolute blade exit velocity, C3r (m/s): Hub = %f, Mid = %f, Tip = %f\n', ...
        C_r3s1_IPT, C_3s1_IPT, C_t3s1_IPT);
end




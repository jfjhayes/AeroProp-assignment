function results = turbineFunction(N_HP, N_IP, r_HPTm, r_IPTm, alpha_3_HPT, phi_HPT, alpha_3_IPT, phi_IPT)
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

    % Constraints % 
    psi_limit = 2.695;          % Stage loading limit (toleranced)

    %N_IP = 158.3;             % Change based on turbine results?
    %N_HP = 175;

    alpha_1_HPT = 0;

    % From section 2 %
    eta_HPT = 0.9367;
    eta_IPT = 0.9318;
    eta_LPT = 0.9396;

    % Critical Pressure Ratio % 
    P_crit = ((gamma_g + 1) / 2)^((gamma_g) / (gamma_g - 1));

    % Stuff to Vary %
    %r_HPTm = 0.4;              % Allowed to be between 0.3 and 0.4 (prob tolerance later)
    %r_IPTm = 0.4;
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
    %alpha_3_HPT = deg2rad(40);
    %phi_HPT = 0.8;

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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%% INTERMEDIATE PRESSURE STUFF BELOW DONT GET CONFUSEEEEEED AGAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    psi_IPT = (C_Pg * (T_o45 - T_o47)) / ((2 * pi * r_IPTm * N_IP)^2);      % IPT stage loading (shoudl be less than 2.695)
    %phi_IPT = 0.8;

    alpha_1_IPT = alpha_3_HPT;
    %alpha_3_IPT = deg2rad(17.6);            % Exit angle limit (toleranced)


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

    results = struct();

    % Temperature Data
    results.T_o4   = T_o4;
    results.T_o45  = T_o45;
    results.T_o47  = T_o47;

    % Pressure Data
    results.P_o4   = P_o4;
    results.P_o45  = P_o45;
    results.P_o47  = P_o47;

    % Blade Speeds
    results.U_HPT  = U_HPT;
    results.U_IPT  = U_IPT;

    % Annulus Areas for HPT
    results.A1_HPT = A_1_HPT;
    results.A2_HPT = A_2_HPT;
    results.A3_HPT = A_3_HPT;

    % Annulus Areas for IPT
    results.A1_IPT = A_1_IPT;
    results.A2_IPT = A_2_IPT;
    results.A3_IPT = A_3_IPT;

    % Blade Angles for HPT
    results.alpha_1t_HPT = alpha_1t_HPT;
    results.alpha_1r_HPT = alpha_1r_HPT;
    results.alpha_2t_HPT = alpha_2t_HPT;
    results.alpha_2r_HPT = alpha_2r_HPT;
    results.alpha_3t_HPT = alpha_3t_HPT;
    results.alpha_3r_HPT = alpha_3r_HPT;

    % Blade Angles for IPT
    results.alpha_1t_IPT = alpha_1t_IPT;
    results.alpha_1r_IPT = alpha_1r_IPT;
    results.alpha_2t_IPT = alpha_2t_IPT;
    results.alpha_2r_IPT = alpha_2r_IPT;
    results.alpha_3t_IPT = alpha_3t_IPT;
    results.alpha_3r_IPT = alpha_3r_IPT;

    results.psi_HPT = psi_HPT;
    results.lambda_HPT = lambda_HPT;
    results.psi_IPT = psi_IPT;
    results.lambda_IPT = lambda_IPT;
    results.h1_HPT = h1_HPT;
    results.h2_HPT = h2_HPT;
    results.h3_HPT = h3_HPT;

    results.h1_IPT = h1_IPT;
    results.h2_IPT = h2_IPT;
    results.h3_IPT = h3_IPT;

    results.beta_2r_HPT = beta_2r_HPT;
    results.beta_2t_HPT = beta_2t_HPT;
    results.beta_3r_HPT = beta_3r_HPT;
    results.beta_3t_HPT = beta_3t_HPT;

    results.beta_2r_IPT = beta_2r_IPT;
    results.beta_2t_IPT = beta_2t_IPT;
    results.beta_3r_IPT = beta_3r_IPT;
    results.beta_3t_IPT = beta_3t_IPT;


    results.rt1_HPT = rt1_HPT;
    results.rr1_HPT = rr1_HPT;

    results.lambda_rs1_HPT = lambda_rs1_HPT;
    results.lambda_ts1_HPT = lambda_ts1_HPT;
    results.lambda_rs1_IPT = lambda_rs1_IPT;
    results.lambda_ts1_IPT = lambda_ts1_IPT;
    



end

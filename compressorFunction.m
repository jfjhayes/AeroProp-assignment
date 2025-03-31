function [N_stage_HP, N_stage_IP, r_IPm, r_HPm] =  compressorFunction(M_rel1, N_IP, N_HP)

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

end
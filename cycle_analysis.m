%% ENG5313: Aerospace Propulsion M Coursework 
% Setion 2 - Cycle Analysis %
clc;
clf;
clear;

%% Section 1 Parameters
Thrust_Loading = 0.43;                      % T_SL / MTOW
Wing_Loading = 5.8;                         % MTOW / S

% Constants %
gamma = 1.4;
P_std = 101.325;                            % Pressure standard day (kN/m^2)
T_std = 288.15;                             % Temp standard day (K)
T_std_h = 312.6;                            % Temp standard hot day (k)
rho_SL = 1.225;                             % Density @ SL (kg/m^3)
g = 9.81;                                   % Standard acceleration of gravity (m/s^2)
R = 287;                                    % Specific gas constant of dry air (J/kgK)

% Conversion factors %
ft2m = 0.3408;                              % feet to metres
sqft2m2 = 0.092903;                         % sqft to m^2
lbft2Pa = 47.88;                            % pound per foot to Pascal
lbf2N = 4.4482;                             % pound force to Newton
lb2kg = 0.4536;                             % pound to kg
nmi2m = 1852;                               % nautical miles to m

% Imperial parameters % 
MTOW_imp = 174200;                          % MTOW (lb)
MPW_imp = 44100;                            % Max payload mass (lb)
PW_imp = 0.8 * MPW_imp;                     % Actual payload mass (lb)
OEW_imp = 97700;                            % Operating empty mass (lb)
M_f_imp = MTOW_imp - PW_imp - OEW_imp;      % Fuel mass (lb)
S_imp = 1330;                               % Wing area (sqft)

Alt_cr_imp = 39800;                         % Cruise altitude (ft)
Rng_cr_imp = 3500;                          % Cruise range (nmi)

% Aircraft parameters (metric) % 
MTOW = MTOW_imp * lb2kg;                    % MTOW (kg)
MPW = MPW_imp * lb2kg;                      % Max payload mass (kg)
PW = PW_imp * lb2kg;                        % Actual payload mass (kg)
OEW = OEW_imp * lb2kg;                      % Operating empty mass (kg)
M_f = M_f_imp * lb2kg;                      % Fuel mass (kg)
S = S_imp * sqft2m2;                        % Wing area (m^2)

Ma_cr = 0.78;                               % Cruise speed (Mach no.)
Ma_cr_MAX = 0.82;                           % Max cruise speed (Mach no.)
Alt_cr = Alt_cr_imp * ft2m;                 % Cruise altitude (m)
Rng_cr = Rng_cr_imp * nmi2m;                % Cruise range (m)

C_D0 = 0.023;                               % Drag polar 1
k = 0.0334;                                 % Drag polar 2

% Fuel fractions %
M_f_climb = 0.04 * M_f;                                 % Mass fuel burned during climb (kg)
M_f_reserve = 0.15 * M_f;                               % Mass fuel reserve (kg)
M_f_descent = 0.04 * M_f;                               % Mass fuel burned during descent (kg)
M_cr = MTOW - M_f_climb;                                % Mass @ start of cruise (kg)
M_cr_end =  M_f_reserve + M_f_descent + OEW + PW;       % Mass @ end of cruise (kg)

% Atmospheric cruise conditions @ approx 13.5 km from ISA table % 
delta_cr = 0.1506;                          % P/P_std
theta_cr = 0.7519;                          % T/T_std for standard day
sigma_cr = delta_cr / theta_cr;                                  
theta_cr0 = theta_cr * (1 + ((gamma-1)/2) * Ma_cr^2);
delta_cr0 = delta_cr * (1 + ((gamma-1)/2) * Ma_cr^2)^(gamma/(gamma-1));
TR = 1.07;                                  % Throttle ratio
alpha_cr = delta_cr0 * (1 - (0.49 * sqrt(Ma_cr)));
beta_cr = M_cr / MTOW;

% Cruise parameters % 
T_cr = theta_cr * T_std;                    % Static temp @ approx 13.5 km on standard day (K)
rho_cr = sigma_cr * rho_SL;                 % Air density @ cruise (kg/m^3)
U_cr = Ma_cr * sqrt(gamma * R * T_cr);      % Flight speed (m/s)
q_cr = 0.5 * rho_cr * (U_cr^2);             % Dynamic pressure @ cruise (kg/m^3)

% Coefficients in cruise - LIFT = WEIGHT, THRUST = DRAG %
L_cr = M_cr * g;                            % Lift = weight
CL_cr = L_cr / (q_cr * S);                  % Coefficient of lift @ cruise
CD_cr = C_D0 + (k * CL_cr^2);               % Coefficient of drag @ cruise, k2=0
D_cr = CD_cr * q_cr * S;                    % Total drag force at flight speed

% Thrusts %
Thrust_SL = Thrust_Loading * MTOW * g;           % Sea-level thrust (N)
Thrust_cr = D_cr;                                % Thrust = Drag @ cruise (N)
SFC_cr = - (U_cr / g) * ((L_cr / D_cr) / Rng_cr) * log(M_cr_end / M_cr);    % (kg/Ns)


% REQUIREMENTS FROM SECTION 1 %
SFC_req = SFC_cr * 1000;     % SFC (kg/kNs)
Thrust_req = Thrust_cr / 2000;  % Single-engine cruise thrust requirement (kN)

%% Section 2 - Cycle Analysis
% Redefine pre-combustion constants - WORK EXLUSIVELY IN KILONEWTONS FROM NOW ON % 
gamma = 1.4;
P_std = 101.325;                            % Pressure standard day (kN/m^2)
T_std = 288.15;                             % Temp standard day (K)
rho_SL = 1.225;                             % Density @ SL (kg/m^3)
g = 9.81;                                   % Standard acceleration of gravity (m/s^2)
R = 0.287;                                  % Specific gas constant of dry air (kJ/kgK)
C_P = 1.005;                                % Specific heat capacity of dry air (kJ/kgK)
Q_r = 43100;                                % Calorific value of fuel (kJ/kgK)

% Post-combustion constants %
gamma_g = 1.333;
C_Pg = 1.148;                               % Specific heat capacity of HOT air (kJ/kgK)

% Efficiencies %
e_cH = 0.915;                               % HP compressor polytropic efficiency 
e_cI = 0.915;                               % IP compressor polytropic efficiency
e_tH = 0.93;                                % HP turbine polytropic efficiency
e_tI = 0.93;                                % IP turbine polytropic efficiency
e_tL = 0.92;                                % LP turbine polytropic efficiency

eta_i = 0.98;                               % Intake isentropic efficiency
eta_b = 0.99;                               % Combustion burner isentropic efficiency
eta_f = 0.91;                               % Fan isentropic efficiency
eta_n = 0.985;                              % Exit nozzle (core + bypass) isentropic efficiency
eta_mech = 0.985;                           % Mechanical efficiency

DeltaP_comb = 1 - 0.04;                     % Total pressure drop in combustion chamber
DeltaP_byp = 1 - 0.02;                      % Total pressure drop in bypass duct

% Limit Values - DO NOT EXCEED %
T_o4_MAX = 1900;                            % Max turbine entry temperature (K)
pi_o_MAX = 45;                              % Max overall pressure ratio (from fan inlet-compressor exit)
pi_f_MAX = 2.1;                             % Max fan pressure ratio
B_MAX = 11.5;                               % Max bypass ratio (not the car)
M_2_MIN = 0.4;                              % Min fan inlet Mach no. (Ma)
M_2_MAX = 0.62;                             % Max fan inlet Mach No. (Ma)
V_19_9_MIN = 0.65;                          % Min nozzle exit velocity (bypass/core) ratio
V_19_9_MAX = 0.85;                          % Max nozzle exit velocity (bypass/core) ratio

% Cruise Conditions @ Approx 13.5km from ISA Tables - COPY FROM ABOVE %
TR = 1.07;                                  % Throttle ratio
Alt_cr = Alt_cr_imp * ft2m;                 % Cruise altitude (m)
delta_cr = 0.1506;                          % P/P_std
theta_cr = 0.7519;                          % T/T_std for standard day
sigma_cr = delta_cr / theta_cr;                                  
theta_cr0 = theta_cr * (1 + ((gamma-1)/2) * Ma_cr^2);
delta_cr0 = delta_cr * (1 + ((gamma-1)/2) * Ma_cr^2)^(gamma/(gamma-1));
alpha_cr = delta_cr0 * (1 - (0.49 * sqrt(Ma_cr)));
beta_cr = M_cr / MTOW;

T_cr = theta_cr * T_std;                    % Static temp @ approx 13.5 km on standard day (K)
T_a = T_cr;                                 % Set ambient static temp (K)
P_a = delta_cr * P_std;                     % Set ambient static pressure (kN/m^2)
rho_cr = sigma_cr * rho_SL;                 % Air density @ cruise (kg/m^3)
U_cr = Ma_cr * sqrt(gamma * R * T_cr);      % Flight speed (m/s)
q_cr = 0.5 * rho_cr * (U_cr^2);             % Dynamic pressure @ cruise (kN/m^2)

% Stagnation characteristics INCLUDE KE %
T_oa = T_a * (1 + ((gamma-1)/2) * Ma_cr^2);                         % Ambient stagnation temp (K)
P_oa = P_a * (1 + ((gamma-1)/2) * Ma_cr^2)^(gamma/(gamma-1));       % Ambient stagnation pressure (kN/m^2)

%% Turbine Properties 
B = 11;                                     % Fixed bypass ratio (for now)

% Stuff to Vary % 
T_o4 = 1900;                                % Turbine entry temperature (K)
pi_cH = 10;                                  % HPC pressure ratio
pi_cI = 3;                                  % IPC pressure ratio
pi_fC = 1.5;                               % Fan core pressure ratio

pi_o = pi_fC * pi_cI * pi_cH;               % Overall pressure ratio

if pi_o > pi_o_MAX
    fprintf('Overall Pressure Ratio exceeds limits at %.2f\n', pi_o);
end

pi_fB = pi_fC - 0.05;                        % Fan bypass pressure ratio w/ correction factor

% Engine Inlet Conditions
V_f = Ma_cr * sqrt(gamma * (R*1000) * T_a);                             % Flight velocity using ambient static temperature (convert gamma into right units for this bit)
T_o2 = T_oa;                                                            % Inlet stagnation temp (K) - TEMP CONSTANT
P_o2 = P_a * (1 + (eta_i*((gamma-1)/2) * Ma_cr^2)^(gamma/(gamma-1)));   % Inlet stgnation pressure (kN/m^2) - 15 kN/m^2 = 15 kPa
pi_d = P_o2 / P_oa;                                                     % Total pressure recovery factor

% Fan exit / Bypass Entry - Station 21
P_o21B = P_o2 * pi_fB;                                      % Station 21 bypass stagnation pressure
T_o21sB = T_o2 * (P_o21B / P_o2)^((gamma -1)/gamma);        % S21 ideal bypass temp
T_o21B = ((T_o21sB - T_o2)/eta_f) + T_o2;                   % S21 actual bypass temp

% Fan exit / Core Entry - Station 21
P_o21C = P_o2 * pi_fC;                                      % Station 21 core stagnation pressure
T_o21sC = T_o2 * (P_o21C / P_o2)^((gamma -1)/gamma);        % S21 ideal core temp
T_o21C = ((T_o21sC - T_o2)/eta_f) + T_o2;                   % S21 actual core temp

% Bypass Entry - Station 13 %
P_o13 = P_o21B;            % Fan exit pressure w/ pressure drop 
T_o13 = T_o21B;            % Fan exit temperature, fan is isentropic compressor so T constant 

% Bypass nozzle exit - Station 19 WITH BYPASS DUCT PRESSURE DROP %
P_o19 = DeltaP_byp * P_o13;     % Bypass nozzle exit stagnation pressure
T_o19 = T_o13;                  % Bypass nozzle exit stagnation temp (pressure drop does not affect)

% Critical presure & TEMP - NOTE THAT EFFICIENCIES ARE INCLUDED %
P_crit_B = P_o19 / (1/(1 - ((1/eta_n)*((gamma-1)/(gamma+1))))^(gamma/(gamma-1)));   % Bypass nozzle exit critical pressure
T_crit_B = T_o19 / (1/(1 - ((1/eta_n)*((gamma-1)/(gamma+1)))));                     % Bypass nozzle exit critical temp
pi_choke_B = P_o19 / P_oa;                                                          % Choke pressure ratio

% Is exit nozzle choked? %
if pi_choke_B > P_crit_B
    P_19 = P_crit_B;                        % Actual exit static pressure (equals critical pressure bcos choke)
    T_19 = T_crit_B;                        % Actual exit static temp (equals critical temp bcos choke)
    V_19 = sqrt(gamma * R*1000 * T_19);     % Actual exit velocity (equal to sonic velocity bcos choke)
else
    P_19 = P_a;                                             % Bypass nozzle exit pressure (equal to ambient)
    T_o19s = T_o19 * (P_a / P_o19)^((gamma -1)/gamma);      % Ideal bypass nozzle exit temp
    T_19 = ((T_o19s - T_o19)/eta_n) + T_o19;                % Actual bypass nozzle exit temp
    Ma_19 = sqrt(2 * C_P * (T_o19s - T_19));                % Exit Mach number
    a_19 = sqrt(gamma * R*1000 * T_19);                     % Local speed of sound
    V_19 = Ma_19 * a_19;                                    % Exit velocity
end

% Intermediate Pressure Compressor (Entry @ Station 23) - Conditions at 23 same as 21C %
P_o23 = P_o21C;     % Pressure at IPC entry
T_o23 = T_o21C;     % Temp at IPC entry

eta_IPC = ((pi_cI ^ ((gamma-1)/gamma)) - 1) / ((pi_cI ^ ((gamma-1)/(gamma*e_cI))) - 1);     % IPC isentropic efficiency (using Farokhi Eq. 4.22)
P_o25 = pi_cI * P_o23;                                                                      % Pressure at core (HPC) inlet
T_o25 = T_o23 * (1 + ((1/eta_IPC) * (((P_o25/P_o23)^((gamma-1)/gamma)) - 1)));              % Temp at core (HPC) inlet

% High Pressure Compressor (Entry @ Station 25) % 
eta_HPC = ((pi_cH ^ ((gamma-1)/gamma)) - 1) / ((pi_cH ^ ((gamma-1)/(gamma*e_cH))) - 1);     % HPC isentropic efficiency (using Farokhi Eq. 4.22)
P_o3 = pi_cH * P_o25;                                                                       % Pressure at compressor exit
T_o3 = T_o25 * (1 + ((1/eta_HPC) * (((P_o3/P_o25)^((gamma-1)/gamma)) - 1)));                % Temp at compressor exit

% Turbine Entry - Station 4 %
P_o4 = DeltaP_comb * P_o3;                                      % Pressure loss at turbine entry
f = ((C_Pg * T_o4) - (C_P * T_o3))/ (Q_r - (C_Pg * T_o4));      % Fuel-air ratio
f_a = f / eta_b;                                                % Actual fuel-air ratio

% High Pressure Turbine (Entry @ Station 4) %
T_o45 = T_o4 - (C_P * (T_o3 - T_o25)) / (eta_mech * (1 + f_a) * C_Pg);
eta_HPT = (1 - (T_o45/T_o4)) / (1 - (T_o45/T_o4)^(1/e_tH));
P_o45 = P_o4 * (1 - (1 - (T_o45 / T_o4)) / (eta_HPT))^(gamma_g / (gamma_g - 1));

% Intermediate Pressure Turbine (Entry @ Station 45) % 
T_o47 = T_o45 - (C_P * (T_o25 - T_o23)) / (eta_mech * (1 + f_a) * C_Pg);
eta_IPT = (1 - (T_o47/T_o45)) / (1 - (T_o47/T_o45))^(1/e_tI);
P_o47 = P_o45 * (1 - (1 - (T_o47 / T_o45)) / (eta_IPT))^(gamma_g / (gamma_g - 1));

% Low Pressure Turbine (Entry @ Station 47) - ALSO DRIVES BYPASS FAN %
T_o5 = T_o47 - ((C_P * (T_o23 - T_o2)) + (B * C_P * (T_o13 - T_o2))) / (eta_mech * (1 + f_a) * C_Pg);
eta_LPT = (1 - (T_o5/T_o47)) / (1 - (T_o5/T_o47))^(1/e_tL);
P_o5 = P_o47 * (1 - (1 - (T_o5 / T_o47)) / (eta_LPT))^(gamma_g / (gamma_g - 1));

% Critical presure & TEMP - NOTE THAT EFFICIENCIES ARE INCLUDED %
P_crit_C_ratio = 1 / ((1 - ((1/eta_n)*((gamma_g-1)/(gamma_g+1))))^(gamma_g/(gamma_g-1)));       % Bypass nozzle exit critical pressure ratio
P_crit_C = P_o5 / P_crit_C_ratio;                                                               % Bypass nozzle exit critical pressure
pi_choke_C = P_o5 / P_oa;                                                                       % Choke pressure ratio

% Is exit nozzle choked? %
if pi_choke_C > P_crit_C_ratio
    fprintf('NOZZLE CHOKED \n')
    P_9 = P_crit_C;                                     % Actual exit static pressure (equals critical pressure bcos choke)
    T_9 = T_o5 * (P_9/P_o5)^((gamma_g-1)/gamma_g);
    V_9 = sqrt(gamma_g * R*1000 * T_9);                 % Actual exit velocity (equal to sonic velocity bcos choke)
else
    fprintf('NOZZLE CLEAR \n')
    P_9 = P_a;
    T_o9s = T_o5 * (P_9 / P_o5)^((gamma_g-1)/gamma_g);
    T_9 = eta_n * (T_o9s - T_o5) + T_o5;
    Ma_9 = sqrt(2 * C_P*1000 * (T_o9s - T_9));
    a_9 = sqrt(gamma_g * R*1000 * T_9);
    V_9 = Ma_9 * a_9;
end

% Velocity ratio % 
V_19_9 = V_19 / V_9;

% Engine Performance
F_S_core = 1 / (B + 1) * (((1 + f_a) * V_9) - V_f) + (((1 + f_a) * R*1000 * T_9) / (V_9 * P_9 * (B + 1)) * (P_9 - P_a));
F_S_bypass = ((B / (B + 1)) * (V_19 - V_f)) + (( (R*1000 * T_19 * B) / (V_19 * P_19 * (B + 1)) ) * (P_19 - P_a));
F_S = F_S_core + F_S_bypass;
sfc = f_a / ((B + 1) * F_S) * 1000;

% Efficiencies % 
eta_P = (F_S * V_f) / ((F_S * V_f) + 0.5 * (((1/(B+1)) * (V_9 - V_f)^2) + ((B/(B+1)) * (V_19 - V_f)^2)));                       % Propulsive efficiency
eta_th = ((F_S * V_f) + 0.5 * (((1/(B+1)) * (V_9 - V_f)^2) + ((B/(B+1)) * (V_19 - V_f)^2))) / ((f_a / (B + 1)) * Q_r*1000);     % Thermal efficiency
eta_o = (F_S * V_f) / ((f_a / (B + 1)) * Q_r*1000);                                                                             % Overall efficiency

% Engine mass flow % 
mdot_in = (Thrust_req*1000) / F_S;      % Inlet mass flow
mdot_core = (1/ (B+1)) * mdot_in;       % Core mass flow
mdot_bypass = (B/ (B+1)) * mdot_in;     % Bypass mass flow

% Nozzle exit areas %
A_9 = ((1+f_a) * mdot_core) / ((P_9 / (R*T_9)) * V_9);          % Core nozzle exit area
A_19 = ((1+f_a) * mdot_bypass) / ((P_19 / (R*T_19)) * V_19);    % Bypass nozzle exit area
A_exit = A_9 + A_19;                                            % Total nozzle exit area
d_exit = 2 * sqrt(A_exit / pi);                                 % Nozzle exit diameter

% Engine inlet area - Method of Choked Empty Annulus %
HTR_ratio = 0.35;                                                               % Ratio of Fan hub-tip radius (from NLA example)
T_crit_F_ideal = (T_o2) / (1 + ((gamma-1) / 2));                                % Ideal inlet critical temp
P_crit_F = P_o2 * (T_crit_F_ideal / T_o2)^(gamma/(gamma-1));                    % Inlet critical pressure
T_crit_F = T_o2 - (eta_i * (T_o2 - T_crit_F_ideal));                            % Inlet actual critical temp

V_crit_F = sqrt(gamma * R*1000 * T_crit_F);                                     % Velocity of choked flow (EQUAL TO SONIC OH YEAH)
mdot_emptyannulus_2_area = (P_crit_F * V_crit_F) / (R * T_crit_F);              % Mass flow per unit area for choked empty annulus
mdot_actualannulus_2_area = 0.88 * mdot_emptyannulus_2_area;                    % Value from NLA example
A_in = mdot_in / mdot_actualannulus_2_area;                                     % Fan area
r_t = sqrt(A_in / (pi * (1 - (HTR_ratio)^2)));                                  % Tip radius
d_in = 2 * r_t;                                                                 % Fan inlet diameter

% Fan inlet mach number % 
T_2 = T_o2 - eta_i * (T_o2 - T_crit_F_ideal);       % Actual static temperature at fan inlet
P_2 = P_o2 * (T_2 / T_o2)^(gamma/(gamma-1));        % Static pressure at fan inlet
rho_2 = P_2 / (R * T_2);                            % Density at fan inlet
V_2 = mdot_actualannulus_2_area / rho_2;            % Fan inlet velocity
a_2 = sqrt(gamma * R*1000 * T_2);                   % Speed of sound at fan inlet
Ma_2 = V_2 / a_2;                                   % Fan inlet Mach number

% Thrusts % 
F_T_core = mdot_in * F_S_core;
F_T_bypass = mdot_in * F_S_bypass;
F_T = F_T_core + F_T_bypass;

% Thrust Ratio % 
F_T_ratio = (F_T_core / mdot_core) / (F_T_bypass / mdot_bypass);

% ADD OVERALL TEXTUAL OUTPUT

% Textual Output %
printText = true;

if printText
    fprintf('Cruise SFC Requirement: %.6f kg/kNs\n', SFC_req);
    fprintf('Single-engine cruise Thrust Requirement: %.4f kN\n', Thrust_req);
    fprintf('Overall pressure ratio: %.4f \n', pi_o);
    fprintf('Turbine entry temperature: %.4f (K)\n', T_o4);
    fprintf('IPC PR: %.4f \n', pi_cI);
    fprintf('HPC PR: %.4f \n', pi_cH);
    fprintf('Core Fan PR: %.4f \n', pi_fC);
    fprintf('---------------INTAKE---------------\n');
    fprintf('Ambient static pressure %.4f (kN/m^2)\n', P_a);
    fprintf('Ambient static temp: %.4f (K)\n', T_a);
    fprintf('Flight velocity: %.4f (m/s)\n', V_f);
    fprintf('Inlet total temperature %.4f (K)\n', T_oa);
    fprintf('Inlet total pressure %.4f (kN/m^2)\n', P_oa);
    fprintf('Fan entry pressure %.4f (kN/m^2)\n', P_o2);
    fprintf('Fan entry temperature %.4f (K)\n', T_o2);
    fprintf('Fan inlet Mach No: %.4f \n', Ma_2);
    fprintf('Engine inlet mass flow: %.4f (kg/s)\n', mdot_in);
    fprintf('---------------BYPASS---------------\n');
    fprintf('Fan exit pressure: %.4f (kN/m^2)\n', P_o13);
    fprintf('Fan exit temp: %.4f (K)\n', T_o13);
    fprintf('Bypass nozzle exit total pressure: %.4f (kN/m^2)\n', P_o19);
    fprintf('Bypass nozzle exit actual pressure: %.4f (kN/m^2)\n', P_19);
    fprintf('Bypass nozzle exit temperature: %.4f (K)\n', T_19);
    fprintf('Bypass nozzle exit velocity: %.4f (m/s)\n', V_19);
    fprintf('Bypass mass flow: %.4f (kg/s)\n', mdot_bypass);
    fprintf('----------------CORE----------------\n');
    fprintf('IPC inlet pressure: %.4f (kN/m^2)\n', P_o23);
    fprintf('IPC inlet temperature: %.4f (K)\n', T_o23);
    fprintf('HPC inlet pressure: %.4f (kN/m^2)\n', P_o25);
    fprintf('HPC inlet temperature: %.4f (K)\n', T_o25);
    fprintf('Compressor exit pressure: %.4f (kN/m^2)\n', P_o3);
    fprintf('Compressor exit temperature: %.4f (K)\n', T_o3);
    fprintf('Actual fuel-air ratio: %.4f \n', f_a);
    fprintf('Turbine entry pressure: %.4f (kN/m^2) \n', P_o4);
    fprintf('Turbine entry temperature: %.4f (K)\n', T_o4);
    fprintf('HPT exit pressure: %.4f (kN/m^2) \n', P_o45);
    fprintf('HPT exit temperature: %.4f (K) \n', T_o45);
    fprintf('IPT exit pressure: %.4f (kN/m^2) \n', P_o47);
    fprintf('IPT exit temperature: %.4f (K) \n', T_o47);
    fprintf('LPT exit pressure: %.4f (kN/m^2) \n', P_o5);
    fprintf('LPT exit temperature: %.4f (K) \n', T_o5);
    fprintf('Core nozzle exit actual pressure: %.4f (kN/m^2)\n', P_9);
    fprintf('Core nozzle exit temperature: %.4f (K)\n', T_9);
    fprintf('Core nozzle exit velocity: %.4f (m/s)\n', V_9);
    fprintf('Core mass flow: %.4f (kg/s)\n', mdot_core);
    fprintf('--------------OVERALL--------------\n');
    fprintf('Specific Thrust: %.4f (Ns/kg) \n', F_S);
    fprintf('Specific Fuel Consumption: %.4f (kg/kNs) \n', sfc);
    fprintf('Propulsive efficiency: %.4f \n', eta_P);
    fprintf('Thermal efficiency: %.4f \n', eta_th);
    fprintf('Overall efficiency: %.4f \n', eta_o);
    fprintf('Fan hub-tip diameter ratio: %.4f \n', HTR_ratio);
    fprintf('Fan inlet area: %.4f (m^2) \n', A_in);
    fprintf('Fan inlet diameter: %.4f (m) \n', d_in);
    fprintf('Total nozzle exit area: %.4f (m^2) \n', A_exit);
    fprintf('Total nozzle exit diameter: %.4f (m) \n', d_exit);
    fprintf('Total thrust: %.4f (kN) \n', F_T/1000);
    fprintf('Thrust ratio: %.4f \n', F_T_ratio);
end

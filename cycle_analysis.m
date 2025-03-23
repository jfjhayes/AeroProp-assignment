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
fprintf('Cruise SFC Requirement: %.10f kg/kNs\n', SFC_req);
fprintf('Single-engine cruise Thrust Requirement: %.10f kN\n', Thrust_req);

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

DeltaP_comb = 0.04;                         % Total pressure drop in combustion chamber
DeltaP_byp = 0.02;                          % Total pressure drop in bypass duct

% Limit Values - DO NOT EXCEED %
T_o4_MAX = 1900;                            % Max turbine entry temperature (K)
pi_o3_o2_MAX = 45;                          % Max pressure ratio (from fan inlet-compressor exit)
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
P_a = delta_cr * P_std;                     % Set ambient static pressure (kg/m^3)
rho_cr = sigma_cr * rho_SL;                 % Air density @ cruise (kg/m^3)
U_cr = Ma_cr * sqrt(gamma * R * T_cr);      % Flight speed (m/s)
q_cr = 0.5 * rho_cr * (U_cr^2);             % Dynamic pressure @ cruise (kg/m^3)


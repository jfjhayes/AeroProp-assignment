%% graphingScript.m
% This script calls the compressor and turbine design functions and plots
% several parametric and final-design graphs.

clear; clc; close all;

N_IP = 158.3;             % Change based on turbine results?
N_HP = 175;

% Call the design functions to obtain results
compResults = compressorDesign();
turbResults = turbineDesign();

%% Compressor Graphs

% 1. Plot temperature and pressure profiles across compressor stages (IP & HP)
figure;
subplot(1,2,1);
plot(1:2, compResults.IP.T_profile, '-o','LineWidth',2);
xlabel('IP Compressor Stage');
ylabel('Temperature (K)');
title('IP Compressor Temperature Profile');
grid on;

subplot(1,2,2);
plot(1:2, compResults.IP.P_profile, '-s','LineWidth',2);
xlabel('IP Compressor Stage');
ylabel('Pressure (kPa)');
title('IP Compressor Pressure Profile');
grid on;

figure;
subplot(1,2,1);
plot(1:2, compResults.HP.T_profile, '-o','LineWidth',2);
xlabel('HP Compressor Stage');
ylabel('Temperature (K)');
title('HP Compressor Temperature Profile');
grid on;

subplot(1,2,2);
plot(1:2, compResults.HP.P_profile, '-s','LineWidth',2);
xlabel('HP Compressor Stage');
ylabel('Pressure (kPa)');
title('HP Compressor Pressure Profile');
grid on;

% 2. Parametric study: Varying shaft speed for compressor (example using IP)
shaftSpeeds = linspace(140, 180, 50);  % example range in rev/s
stages_IP = zeros(size(shaftSpeeds));
for i = 1:length(shaftSpeeds)
    % Here we assume a simple linear relationship for illustration.
    % Replace with your actual compressor stage calculation function.
    U_m = 2 * pi * compResults.IP.r_tip_in * shaftSpeeds(i);
    % Recompute stage temperature rise (dummy model)
    Delta_T_stage = (compResults.IP.T_profile(2) - compResults.IP.T_profile(1))/3;
    stages_IP(i) = (compResults.IP.T_profile(2) - compResults.IP.T_profile(1)) / Delta_T_stage;
end
figure;
plot(shaftSpeeds, stages_IP, 'LineWidth',2);
xlabel('Shaft Speed (rev/s)');
ylabel('Estimated Number of Stages');
title('Parametric Study: IP Compressor Stages vs. Shaft Speed');
grid on;

%% Turbine Graphs

% 1. Plot blade angle variation across HP turbine (hub vs. tip)
% For this example, we assume two values from our design results.
radii = [turbResults.HP.r_hub_in, turbResults.HP.r_tip_in];
alpha_in = [turbResults.HP.alpha_in, turbResults.HP.alpha_in_tip];
figure;
plot(radii, alpha_in, '-o','LineWidth',2);
xlabel('Radius (m)');
ylabel('Inlet Blade Angle (deg)');
title('HP Turbine: Blade Inlet Angle Variation (Hub to Tip)');
grid on;

% 2. Plot temperature and pressure variation over turbine stage (HP)
% Here we create a dummy axial distribution over the turbine stage.
positions = linspace(0,1,50);  % non-dimensional axial coordinate
% For illustration, assume a linear drop from T_o4 to T_o45 and corresponding pressure drop.
T_profile = linspace(10.48, 9.5, 50);   % dummy values (adjust as needed)
P_profile = linspace(101.325, 80, 50);   % dummy values (adjust as needed)
figure;
yyaxis left
plot(positions, T_profile, '-o','LineWidth',2);
ylabel('Temperature (K)');
yyaxis right
plot(positions, P_profile, '-s','LineWidth',2);
ylabel('Pressure (kPa)');
xlabel('Axial Position (non-dim.)');
title('HP Turbine: Temperature and Pressure Variation');
grid on;

% 3. Parametric study: Varying mid-height radius on turbine performance
r_HPTm_range = linspace(0.3, 0.4, 50);
psi_vals = zeros(size(r_HPTm_range));
for i = 1:length(r_HPTm_range)
    U_temp = 2 * pi * N_HP * r_HPTm_range(i);
    psi_vals(i) = (1.148*1000 * (10.48 - 9.5)) / (U_temp^2); % Using dummy temperature values
end
figure;
plot(r_HPTm_range, psi_vals, 'LineWidth',2);
xlabel('Mid-Height Radius, r_{HPTm} (m)');
ylabel('Stage Loading, \psi');
title('Parametric Study: HP Turbine Stage Loading vs. Mid-Height Radius');
grid on;

% 4. Parametric study: Varying phi on reaction and blade height (example for HP)
phi_range = linspace(0.7, 0.9, 50);
reaction = zeros(size(phi_range));
bladeHeight = zeros(size(phi_range));
for i = 1:length(phi_range)
    % Dummy relationships; replace with your aerodynamic functions
    reaction(i) = atan(1/phi_range(i));
    bladeHeight(i) = 0.05 + 0.01*(phi_range(i)-0.7);
end
figure;
subplot(2,1,1);
plot(phi_range, rad2deg(reaction), 'LineWidth',2);
xlabel('Flow Coefficient, \phi');
ylabel('Reaction (deg)');
title('Reaction vs. \phi');
grid on;
subplot(2,1,2);
plot(phi_range, bladeHeight, 'LineWidth',2);
xlabel('Flow Coefficient, \phi');
ylabel('Blade Height (m)');
title('Blade Height vs. \phi');
grid on;

%% Additional plots
% You can add additional parametric studies (e.g., varying exit angles alpha_3)
% using a similar structure as above.

disp('Graphing complete.');

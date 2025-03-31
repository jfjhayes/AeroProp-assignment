% Define a range for the mean radius
r_HPTm_range = linspace(0.3, 0.4, 10);

% Preallocate arrays for blade angles at tip and root (example for HPT exit angle)
alpha_3t_arr = zeros(size(r_HPTm_range));
alpha_3r_arr = zeros(size(r_HPTm_range));

% Fixed parameters
N_HP         = 175;
N_IP         = 158.3;
r_IPTm       = 0.4;
alpha_3_HPT  = deg2rad(40);
phi_HPT      = 0.8;
alpha_3_IPT  = deg2rad(17.6);
phi_IPT      = 0.8;

for i = 1:length(r_HPTm_range)
    r_HPTm = r_HPTm_range(i);
    
    % Call the turbine function with varying mean radius
    results = turbineFunction(N_HP, N_IP, r_HPTm, r_IPTm, alpha_3_HPT, phi_HPT, alpha_3_IPT, phi_IPT);
    
    % Extract the tip and root exit blade angles for HPT (alpha_3t and alpha_3r)
    alpha_3t_arr(i) = rad2deg(results.alpha_3t_HPT);  % converting back to degrees for plotting
    alpha_3r_arr(i) = rad2deg(results.alpha_3r_HPT);
end

% Plot blade angle variation with mean radius
figure;
plot(r_HPTm_range, alpha_3t_arr, 'r-o', r_HPTm_range, alpha_3r_arr, 'b-s');
xlabel('Mean Radius r_{HPTm} (m)');
ylabel('Exit Blade Angle (deg)');
legend('\alpha_{3t}', '\alpha_{3r}');
title('Radial Distribution of Exit Blade Angles vs. Mean Radius');

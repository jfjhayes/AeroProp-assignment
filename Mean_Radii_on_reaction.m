% Mean Radius on Reaction

r_HPTm_range = linspace(0.25, 0.45, 200);
r_IPTm_range = linspace(0.25, 0.45, 200);

N_HP_range = linspace(100,200,200);
N_IP_range = linspace(100,200,200);

%N_IP = 158.3;
N_HP = 175;
alpha_3_HPT  = deg2rad(40);
phi_HPT      = 0.8;
alpha_3_IPT  = deg2rad(17.6);
phi_IPT      = 0.8;

[r_HPTm_grid, N_HP_grid] = meshgrid(r_HPTm_range, N_HP_range);
lambda_HPT_grid = zeros(size(r_HPTm_grid));

for i = 1:numel(r_HPTm_grid)
    results = turbineFunction(N_HP_grid(i), N_IP, r_HPTm_grid(i), r_IPTm, alpha_3_HPT, phi_HPT, alpha_3_IPT, phi_IPT);
    lambda_HPT_grid(i) = rad2deg((results.lambda_rs1_HPT + results.lambda_ts1_HPT)/2);
end

figure;
contourf(r_HPTm_grid, lambda_HPT_grid, N_HP_grid, 20, 'LineColor', 'none');
shading interp;
h = colorbar;
ylabel(h, '$N_{HP}$', 'Interpreter', 'latex');
xlabel('$r_{HPT,m}$', 'Interpreter', 'latex');
ylabel('$\Lambda_{HPT}$', 'Interpreter', 'latex');
set(gca, "TickLabelInterpreter", 'latex');
saveas(gcf, 'HPMeanRadius_v_Reaction_for_Speed.eps', 'epsc');


N_IP = 158.3;
[r_IPTm_grid, N_IP_grid] = meshgrid(r_IPTm_range, N_IP_range);
lambda_IPT_grid = zeros(size(r_IPTm_grid));

for i = 1:numel(r_IPTm_grid)
    results = turbineFunction(N_HP, N_IP_grid(i), r_HPTm, r_IPTm_grid(i), alpha_3_HPT, phi_HPT, alpha_3_IPT, phi_IPT);
    lambda_IPT_grid(i) = rad2deg((results.lambda_rs1_IPT + results.lambda_ts1_IPT)/2);
end

figure;
contourf(r_IPTm_grid, lambda_IPT_grid, N_IP_grid, 20, 'LineColor', 'none');
shading interp;
h = colorbar;
ylabel(h, '$N_{IP}$', 'Interpreter', 'latex');
xlabel('$r_{IPT,m}$', 'Interpreter', 'latex');
ylabel('$\Lambda_{IPT}$', 'Interpreter', 'latex');
set(gca, "TickLabelInterpreter", 'latex');
saveas(gcf, 'IPMeanRadius_v_Reaction_for_Speed.eps', 'epsc');

%results.psi_HPT
%results.lambda_HPT
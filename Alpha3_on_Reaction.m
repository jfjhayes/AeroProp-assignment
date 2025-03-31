% Mean Radius on Stage Loading and Reaction

alpha_3_HPT_range = linspace(0, 65, 200);
alpha_3_IPT_range = linspace(0, 65, 200);

N_HP_range = linspace(100,200,200);
N_IP_range = linspace(100,200,200);

alpha_3_HPT  = deg2rad(40);
phi_HPT      = 0.8;
alpha_3_IPT  = deg2rad(17.6);
phi_IPT      = 0.8;

[alpha_3_HPT_grid, N_HP_grid] = meshgrid(alpha_3_HPT_range, N_HP_range);
lambda_HPT_grid = zeros(size(alpha_3_HPT_grid));
N_IP = 158.3;
r_HPTm = 0.4;
r_IPTm = 0.4;

for i = 1:numel(alpha_3_HPT_grid)
    results = turbineFunction(N_HP_grid(i), N_IP, r_HPTm, r_IPTm, deg2rad(alpha_3_HPT_grid(i)), phi_HPT, alpha_3_IPT, phi_IPT);
    lambda_HPT_grid(i) = rad2deg((results.lambda_rs1_HPT + results.lambda_ts1_HPT)/2);
end

figure;
contourf(alpha_3_HPT_grid, lambda_HPT_grid, N_HP_grid, 20, 'LineColor', 'none');
shading interp;
h = colorbar;
ylabel(h, '$N_{HP} (rps)$', 'Interpreter', 'latex');
xlabel('$\alpha_{3,HPT}$ (deg)', 'Interpreter', 'latex');
ylabel('$\Lambda_{HPT}$ (deg)', 'Interpreter', 'latex');
%yline(0.04, 'k--','Label', '$h_{2,limit}$', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');
set(gca, "TickLabelInterpreter", 'latex');
saveas(gcf, 'HPT_Alpha3_v_Reaction_for_Speed.eps', 'epsc');

[alpha_3_IPT_grid, N_IP_grid] = meshgrid(alpha_3_IPT_range, N_IP_range);
lambda_IPT_grid = zeros(size(alpha_3_IPT_grid));
N_IP = 158.3;
r_HPTm = 0.4;
r_IPTm = 0.4;

for i = 1:numel(alpha_3_IPT_grid)
    results = turbineFunction(N_HP, N_IP_grid(i), r_HPTm, r_IPTm, alpha_3_HPT, phi_HPT, deg2rad(alpha_3_IPT_grid(i)), phi_IPT);
    lambda_IPT_grid(i) = rad2deg((results.lambda_rs1_IPT + results.lambda_ts1_IPT)/2);
end

figure;
contourf(alpha_3_IPT_grid, lambda_IPT_grid, N_IP_grid, 20, 'LineColor', 'none');
shading interp;
h = colorbar;
ylabel(h, '$N_{IP} (rps)$', 'Interpreter', 'latex');
xlabel('$\alpha_{3,IPT}$ (deg)', 'Interpreter', 'latex');
ylabel('$\Lambda_{IPT}$ (deg)', 'Interpreter', 'latex');
%yline(0.04, 'k--','Label', '$h_{2,limit}$', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');
set(gca, "TickLabelInterpreter", 'latex');
saveas(gcf, 'IPT_Alpha3_v_Reaction_for_Speed.eps', 'epsc');

%results.psi_HPT
%results.lambda_HPT
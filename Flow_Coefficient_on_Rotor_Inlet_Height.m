% Mean Radius on Stage Loading and Reaction

phi_HPT_range = linspace(0.4, 1.6, 200);
phi_IPT_range = linspace(0.4, 1.6, 200);

N_HP_range = linspace(100,200,200);
N_IP_range = linspace(100,200,200);

alpha_3_HPT  = deg2rad(40);
phi_HPT      = 0.8;
alpha_3_IPT  = deg2rad(17.6);
phi_IPT      = 0.8;

[phi_HPT_grid, N_HP_grid] = meshgrid(phi_HPT_range, N_HP_range);
h2_HPT_grid = zeros(size(phi_HPT_grid));
N_IP = 158.3;
r_HPTm = 0.4;
r_IPTm = 0.4;

for i = 1:numel(phi_HPT_grid)
    results = turbineFunction(N_HP_grid(i), N_IP, r_HPTm, r_IPTm, alpha_3_HPT, phi_HPT_grid(i), alpha_3_IPT, phi_IPT);
    h2_HPT_grid(i) = results.h2_HPT;
end

figure;
contourf(phi_HPT_grid, h2_HPT_grid, N_HP_grid, 20, 'LineColor', 'none');
shading interp;
h = colorbar;
ylabel(h, '$N_{HP}$', 'Interpreter', 'latex');
xlabel('$\phi_{HPT}$', 'Interpreter', 'latex');
ylabel('$h_{2,HPT}$ (m)', 'Interpreter', 'latex');
%yline(0.04, 'k--','Label', '$h_{2,limit}$', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');
set(gca, "TickLabelInterpreter", 'latex');
saveas(gcf, 'HPT_Flow_Coefficient_v_Blade_Height_for_Speed.eps', 'epsc');

N_HP = 175;
[phi_IPT_grid, N_IP_grid] = meshgrid(phi_IPT_range, N_IP_range);
h2_IPT_grid = zeros(size(phi_IPT_grid));
r_HPTm = 0.4;
r_IPTm = 0.4;

for i = 1:numel(phi_IPT_grid)
    results = turbineFunction(N_HP, N_IP_grid(i), r_HPTm, r_IPTm, alpha_3_HPT, phi_HPT, alpha_3_IPT, phi_IPT_grid(i));
    h2_IPT_grid(i) = results.h2_IPT;
end

figure;
contourf(phi_IPT_grid, h2_IPT_grid, N_IP_grid, 20, 'LineColor', 'none');
shading interp;
h = colorbar;
ylabel(h, '$N_{IP}$', 'Interpreter', 'latex');
xlabel('$\phi_{IPT}$', 'Interpreter', 'latex');
ylabel('$h_{2,IPT}$ (m)', 'Interpreter', 'latex');
set(gca, "TickLabelInterpreter", 'latex');
saveas(gcf, 'IPT_Flow_Coefficient_v_Blade_Height_for_Speed.eps', 'epsc');

%results.psi_HPT
%results.lambda_HPT
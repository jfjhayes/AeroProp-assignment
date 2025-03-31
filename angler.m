N_IP         = 158.3;
r_HPTm       = 0.4;
r_IPTm       = 0.425;
alpha_3_HPT  = deg2rad(30);
phi_HPT      = 0.6;
alpha_3_IPT  = deg2rad(17.6);
phi_IPT      = 0.6;

results = turbineFunction(N_HP, N_IP, r_HPTm, r_IPTm, alpha_3_HPT, phi_HPT, alpha_3_IPT, phi_IPT);

HP_r = [0 0.5 1];

HP_a2 = [ rad2deg(results.alpha_2r_HPT), ...
          rad2deg((results.alpha_2r_HPT + results.alpha_2t_HPT)/2), ...
          rad2deg(results.alpha_2t_HPT) ];

HP_a3 = [ rad2deg(results.alpha_3r_HPT), ...
          rad2deg((results.alpha_3r_HPT + results.alpha_3t_HPT)/2), ...
          rad2deg(results.alpha_3t_HPT) ];

HP_b2 = [ rad2deg(results.beta_2r_HPT), ...
          rad2deg((results.beta_2r_HPT + results.beta_2t_HPT)/2), ...
          rad2deg(results.beta_2t_HPT) ];

HP_b3 = [ rad2deg(results.beta_3r_HPT), ...
          rad2deg((results.beta_3r_HPT + results.beta_3t_HPT)/2), ...
          rad2deg(results.beta_3t_HPT) ];

IP_offset = 2; 
IP_r = [0 0.5 1] + IP_offset;

IP_a1 = [ rad2deg(results.alpha_1r_IPT), ...
          rad2deg((results.alpha_1r_IPT + results.alpha_1t_IPT)/2), ...
          rad2deg(results.alpha_1t_IPT) ];

IP_a2 = [ rad2deg(results.alpha_2r_IPT), ...
          rad2deg((results.alpha_2r_IPT + results.alpha_2t_IPT)/2), ...
          rad2deg(results.alpha_2t_IPT) ];

IP_a3 = [ rad2deg(results.alpha_3r_IPT), ...
          rad2deg((results.alpha_3r_IPT + results.alpha_3t_IPT)/2), ...
          rad2deg(results.alpha_3t_IPT) ];

IP_b2 = [ rad2deg(results.beta_2r_IPT), ...
          rad2deg((results.beta_2r_IPT + results.beta_2t_IPT)/2), ...
          rad2deg(results.beta_2t_IPT) ];

IP_b3 = [ rad2deg(results.beta_3r_IPT), ...
          rad2deg((results.beta_3r_IPT + results.beta_3t_IPT)/2), ...
          rad2deg(results.beta_3t_IPT) ];


figure;
hold on;
grid on;
xlabel('Blade Section', 'Interpreter','latex');
ylabel('Blade Angle (degrees)', 'Interpreter', 'latex');
%title('Combined Blade Angle Diagram for HPT and IPT');

plot(HP_r, HP_a2, 'r-o', 'DisplayName','HPT $\alpha_{2}$');
plot(HP_r, HP_b2, 'r--s', 'DisplayName','HPT $\beta_{2}$');
plot(HP_r, HP_b3, 'r-.^', 'DisplayName','HPT $\beta_{3}$');
plot(HP_r, HP_a3, 'r:', 'DisplayName','HPT $\alpha_{3}$');

% IPT data - shifted %
plot(IP_r, IP_a1, 'b-o', 'DisplayName','IPT $\alpha_{1}$');
plot(IP_r, IP_a2, 'b--s', 'DisplayName','IPT $\alpha_{2}$');
plot(IP_r, IP_b2, 'b-.^', 'DisplayName','IPT $\beta_{2}$');
plot(IP_r, IP_b3, 'b:', 'DisplayName','IPT $\beta_{3}$');
plot(IP_r, IP_a3, 'k-o', 'DisplayName','IPT $\alpha_{3}$');

xline(0.5, '-', 'HandleVisibility','off');           % HPT mean
xline(IP_offset+0.5, '-', 'HandleVisibility','off');   % IPT mean

xticks([0 0.5 1, IP_offset, IP_offset+0.5, IP_offset+1]);
xticklabels({'Hub (HPT)','Mean (HPT)','Tip (HPT)', 'Hub (IPT)','Mean (IPT)','Tip (IPT)'});
legend('Location','best', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'tex');
hold off;
saveas(gcf, 'Angle_Diagram.eps', 'epsc');

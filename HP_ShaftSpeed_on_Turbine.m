N_HP_range = linspace(100, 200, 100); 
N_IP_range = linspace(100, 200, 100) 

T_o4_arr    = zeros(size(N_HP_range));
T_o45_arr   = zeros(size(N_HP_range));
T_o47_arr   = zeros(size(N_HP_range));
U_HPT_arr   = zeros(size(N_HP_range));
A1_HPT_arr  = zeros(size(N_HP_range));
A2_HPT_arr  = zeros(size(N_HP_range));
A3_HPT_arr  = zeros(size(N_HP_range));

U_IPT_arr   = zeros(size(N_IP_range));
A1_IPT_arr  = zeros(size(N_IP_range));
A2_IPT_arr  = zeros(size(N_IP_range));
A3_IPT_arr  = zeros(size(N_IP_range));

% Constants and fixed parameters
N_IP         = 158.3;
r_HPTm       = 0.4;
r_IPTm       = 0.4;
alpha_3_HPT  = deg2rad(40);
phi_HPT      = 0.8;
alpha_3_IPT  = deg2rad(17.6);
phi_IPT      = 0.8;

for i = 1:length(N_HP_range)
    N_HP = N_HP_range(i);
    
    results = turbineFunction(N_HP, N_IP, r_HPTm, r_IPTm, alpha_3_HPT, phi_HPT, alpha_3_IPT, phi_IPT);
    
    % Store
    T_o4_arr(i)    = results.T_o4;
    T_o45_arr(i)   = results.T_o45;
    T_o47_arr(i)   = results.T_o47;
    U_HPT_arr(i)   = results.U_HPT;
    A1_HPT_arr(i)  = results.A1_HPT;
    A2_HPT_arr(i)  = results.A2_HPT;
    A3_HPT_arr(i)  = results.A3_HPT;
end

N_HP = 175;
for i = 1:length(N_IP_range)
    N_IP = N_IP_range(i);
    
    results = turbineFunction(N_HP, N_IP, r_HPTm, r_IPTm, alpha_3_HPT, phi_HPT, alpha_3_IPT, phi_IPT);
    
    % Store
    T_o4_arr(i)    = results.T_o4;
    T_o45_arr(i)   = results.T_o45;
    T_o47_arr(i)   = results.T_o47;
    U_IPT_arr(i)   = results.U_IPT;
    A1_IPT_arr(i)  = results.A1_IPT;
    A2_IPT_arr(i)  = results.A2_IPT;
    A3_IPT_arr(i)  = results.A3_IPT;
end

% Plot Annulus Areas vs. N_HP
figure;
plot(N_HP_range, A1_HPT_arr, 'm-', N_HP_range, A2_HPT_arr, 'c-', N_HP_range, A3_HPT_arr, 'b-');
xlabel('High-pressure Shaft Speed (rps)', 'Interpreter', 'latex');
ylabel('Annulus Area ($m^2$)', 'Interpreter', 'latex');
%title('Annulus Area Variation with Shaft Speed');
xline(175, 'k--','Label', '$N_{HP,limit}$', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');
xline(158.3, 'k--','Label', '$N_{IP,limit}$', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');
legend('$A_{1,HPT}$', '$A_{2,HPT}$', '$A_{3,HPT}$', 'Interpreter', 'latex');
set(gca, "TickLabelInterpreter", 'latex');
saveas(gcf, 'ShaftSpeeds_V_HPT_Area.eps', 'epsc');

% Plot Annulus Areas vs. N_HP
figure;
plot(N_IP_range, A1_IPT_arr, 'm', N_IP_range, A2_HPT_arr, 'c-', N_IP_range, A3_HPT_arr, 'b-');
xlabel('Intermediate-pressure Shaft Speed (rps)', 'Interpreter', 'latex');
ylabel('Annulus Area ($m^2$)', 'Interpreter', 'latex');
%title('Annulus Area Variation with Shaft Speed');
xline(175, 'k--','Label', '$N_{HP,limit}$', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');
xline(158.3, 'k--','Label', '$N_{IP,limit}$', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');
legend('$A_{1,IPT}$', '$A_{2,IPT}$', '$A_{3,IPT}$', 'Interpreter', 'latex');
set(gca, "TickLabelInterpreter", 'latex');
saveas(gcf, 'ShaftSpeeds_V_IPT_Area.eps', 'epsc');

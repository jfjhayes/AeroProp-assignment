clc; clear; close all;

% Define range for M_rel1
M_rel1 = 1.21;

N_IP_values = linspace(100, 200, 100);
N_HP = 175; 

N_stage_IP = zeros(size(N_IP_values));
N_stage_HP = zeros(size(N_IP_values));

% VARY IP SHAFT SPEED
for i = 1:length(N_IP_values)
    [N_stage_HP(i), N_stage_IP(i),~,~] = compressorFunction(M_rel1, N_IP_values(i), N_HP);
end

N_HP_values = linspace(100, 200, 100);
N_IP = 158.3;

N_stage_IP2 = zeros(size(N_HP_values));
N_stage_HP2 = zeros(size(N_HP_values));

% VARY HP SHAFT SPEED
for i = 1:length(N_HP_values)
    [N_stage_HP2(i), N_stage_IP2(i),~,~] = compressorFunction(M_rel1, N_IP, N_HP_values(i));
end

figure;
yyaxis left
% Effect of IP Shaft Speed on Number of IP stages
plot(N_IP_values, N_stage_IP, 'b', 'LineWidth', 1.5); hold on;
% Effect of HP Shaft Speed on Number of IP stages
plot(N_HP_values, N_stage_IP2, 'b--', 'LineWidth', 1.5); 
ylabel('No. of IPC Stages', 'Interpreter', 'latex');

yyaxis right
% Effect of IP Shaft Speed on Number of HP stages
plot(N_IP_values, N_stage_HP, 'r--', 'LineWidth', 1.5); hold on;
% Effect of HP Shaft Speed on Number of HP stages
plot(N_HP_values, N_stage_HP2, 'r-', 'LineWidth', 1.5); 
ylabel('No. of HPC Stages', 'Interpreter', 'latex');


xlabel('Shaft Speed (rps)', 'Interpreter', 'latex');
xline(175, 'k--','Label', '$N_{HP,limit}$', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');
xline(158.3, 'k--','Label', '$N_{IP,limit}$', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');
legend('$N_{IP}$', '$N_{HP}$', '$N_{IP}$', '$N_{HP}$', 'Interpreter', 'latex');
set(gca, "TickLabelInterpreter", 'latex');
saveas(gcf, 'ShaftSpeeds_V_No_Stages_Compressor.eps', 'epsc');

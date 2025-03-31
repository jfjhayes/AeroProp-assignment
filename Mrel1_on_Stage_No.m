clc; clear; close all;

% Define range for M_rel1
M_rel1_values = linspace(0.5, 1.5, 50);
N_IP = 158.3; 
N_HP = 175; 

N_stage_IP = zeros(size(M_rel1_values));
N_stage_HP = zeros(size(M_rel1_values));

for i = 1:length(M_rel1_values)
    [N_stage_HP(i), N_stage_IP(i), ~,~] = compressorFunction(M_rel1_values(i), N_IP, N_HP);
end

figure;
yyaxis left
plot(M_rel1_values, N_stage_IP, 'b', 'LineWidth', 1.5);
ylabel('No. of IPC Stages', 'Interpreter', 'latex');
yyaxis right
plot(M_rel1_values, N_stage_HP, 'r', 'LineWidth', 1.5);
ylabel('No. of HPC Stages', 'Interpreter', 'latex');
xlabel('Relative Mach No. (Ma)', 'Interpreter', 'latex');
xline(1.21, 'k--','Label', '$Ma_{t,rel,limit}$', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');
legend('IPC', 'HPC', 'Interpreter', 'latex');
set(gca, "TickLabelInterpreter", 'latex');
saveas(gcf, 'M_rel1_V_Stage_No_Compressor.eps', 'epsc');


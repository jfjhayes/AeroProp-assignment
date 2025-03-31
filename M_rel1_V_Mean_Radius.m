clc; clear; close all;

% Define range for M_rel1
M_rel1_values = linspace(0.5, 1.5, 50);
N_IP = 158.3; 
N_HP = 175; 

r_IPm = zeros(size(M_rel1_values));
r_HPm = zeros(size(M_rel1_values));

for i = 1:length(M_rel1_values)
    [~,~, r_IPm(i), r_HPm(i)] = compressorFunction(M_rel1_values(i), N_IP, N_HP);
end

figure;
yyaxis left
plot(M_rel1_values, r_IPm, 'b', 'LineWidth', 1.5);
ylabel('IP Mean Radius (m)', 'Interpreter', 'latex');
yyaxis right
plot(M_rel1_values, r_HPm, 'r', 'LineWidth', 1.5);
ylabel('HP Mean Radius (m)', 'Interpreter', 'latex');
xlabel('Relative Mach No. (Ma)', 'Interpreter', 'latex');
xline(1.21, 'k--','Label', '$Ma_{t,rel,limit}$', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');
legend('IPC', 'HPC', 'Interpreter', 'latex');
set(gca, "TickLabelInterpreter", 'latex');
saveas(gcf, 'M_rel1_V_Mean_Radius_Compressor.eps', 'epsc');



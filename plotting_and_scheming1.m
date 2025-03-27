%% ENG5313: Aerospace Propulsion M Coursework 
% Setion 2 - Cycle Analysis %
clc;
clear;
close;

exportMode = true;
graph = 6;


% 1 - pi_fC vs SFC for B
% 2 - pi_fC vs Thrust for B
% 3 - pi_cH vs SFC for B
% 4 - pi_cH vs Thrust for B
% 5 - pi_cI vs SFC for B
% 6 - pi_cI vs Thrust for B
% 7 - pi_o vs SFC for B
% 8 - pi_o vs Thrust for B
% 9 - T_o4 vs SFC for B
% 10 - T_o4 vs Thrust for B
% 11 - pi_o vs Engine Inlet Diameter for B
% 12 - T_o4 vs Engine Inlet Diameter for B
% 13 - pi_o vs Engine Exit Diameter for B
% 14 - T_o4 vs Engine Exit Diameter for B
% 15 - pi_o vs V_19_9 for B
% 16 - T_o4 v V_19_9 for B

if graph == 1
 
    % CORE FAN PRESSURE RATIO vs SFC for B
    HPC = 10;
    IPC = 2.5;
    T_o4 = 2090;
 
    % Define parameter ranges
    pi_fC = linspace(1.2, 2.5, 200); % Fan core pressure ratio range
    B = linspace(4, 11.5, 200);          % Bypass ratio range
 
    % Meshgrid for evaluation
    [pi_fC_grid, B_grid] = meshgrid(pi_fC, B);
    sfc_grid = nan(size(pi_fC_grid)); % Preallocate with NaNs
 
    % Compute SFC for each (pi_fC, B)
    % BE AWARE T_o4 toleranced
    for i = 1:numel(pi_fC_grid)
        [~, sfc, ~, ~, ~, ~] = cycle_analyser(B_grid(i), HPC, IPC, T_o4, pi_fC_grid(i));
       
        % Ensure only real, positive SFC values are stored
        if isreal(sfc) && sfc > 0
            sfc_grid(i) = sfc;
        end
    end
 
    % Actual pi_o_limit: 49.5
    % Ideal pi_o_limit: 45
 
    pi_fC_limit = 60 / (HPC * IPC);
 
    figure;
    contourf(pi_fC_grid, sfc_grid, B_grid, 20, 'LineColor', 'none');
    shading interp;
    h = colorbar;
    ylabel(h, '$B$', 'Interpreter', 'latex');
    xlabel('$\pi_{fC}$', 'Interpreter', 'latex');
    ylabel('$SFC$ (kgkN$^{-1}$s$^{-1}$)', 'Interpreter', 'latex');
    set(gca, "TickLabelInterpreter", 'latex');
    ylim([0.005, 0.057]);
 
    % Actual y-limit: 0.0144529
    % Ideal y-limit:0.013139
    yline(0.0144529, 'k--','Label', '$SFC_{req}$', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');             % SFC requirement
    xline(pi_fC_limit, 'k--','Label', '$\pi_{fC_{lim}}$', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');         % pi_fC limit from pi_o  
   
    if exportMode
        saveas(gcf, 'core_fan_PR_v_SFC_forB.eps', 'epsc');
    end
 
    % Good values are T_o4 2090, HPC 10, IPC 2.5, pi_fC 2.27, B 10.49 - SFC MET


elseif graph == 2

    % CORE FAN PRESSURE RATIO VS SPECIFIC THRUST FOR B
    HPC = 10;
    IPC = 2.5;
    T_o4 = 2090;
 
    % Define parameter ranges
    pi_fC = linspace(1.2, 2.5, 200); % Fan core pressure ratio range
    B = linspace(4, 11.5, 200);          % Bypass ratio range
 
    % Meshgrid for evaluation
    [pi_fC_grid, B_grid] = meshgrid(pi_fC, B);
    F_S_grid = nan(size(pi_fC_grid)); % Preallocate with NaNs
 
    % Compute SFC for each (pi_fC, B)
    % BE AWARE T_o4 toleranced
    for i = 1:numel(pi_fC_grid)
        [F_S, ~, ~, ~, ~, ~] = cycle_analyser(B_grid(i), HPC, IPC, T_o4, pi_fC_grid(i));
       
        % Ensure only real, positive SFC values are stored
        if isreal(F_S) && F_S > 0
            F_S_grid(i) = F_S;
        end
    end
 
    % Actual pi_o_limit: 49.5
    % Ideal pi_o_limit: 45
 
    pi_fC_limit = 60 / (HPC * IPC);
 
    figure;
    contourf(pi_fC_grid, F_S_grid, B_grid, 20, 'LineColor', 'none');
    shading interp;
    h = colorbar;
    set(h, 'TickLabelInterpreter', 'latex');
    ylabel(h, '$B$', 'Interpreter', 'latex');
    xlabel('$\pi_{fC}$', 'Interpreter', 'latex');
    ylabel('$F_S$ (Nskg$^{-1}$)', 'Interpreter', 'latex');
    set(gca, "TickLabelInterpreter", 'latex');
    %ylim([0., 0.057]);
 
    xline(pi_fC_limit, 'k--','Label', '$\pi_{fC_{lim}}$', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');         % pi_fC limit from pi_o  
   
    if exportMode
        saveas(gcf, 'core_fan_PR_v_FS_forB.eps', 'epsc');
    end



elseif graph == 3

    % HPC vs SFC for B
    pi_fC = 2.27;
    IPC = 2.5;
    T_o4 = 2090;
 
    % Define parameter ranges
    pi_cH = linspace(8, 12.5, 200); % Fan core pressure ratio range
    B = linspace(4, 11.5, 200);          % Bypass ratio range
 
    % Meshgrid for evaluation
    [pi_cH_grid, B_grid] = meshgrid(pi_cH, B);
    sfc_grid = nan(size(pi_cH_grid)); % Preallocate with NaNs
 
    % Compute SFC for each (pi_fC, B)
    % BE AWARE T_o4 toleranced
    for i = 1:numel(pi_cH_grid)
        [~, sfc, ~, ~, ~, ~] = cycle_analyser(B_grid(i), pi_cH_grid(i), IPC, T_o4, pi_fC);
       
        % Ensure only real, positive SFC values are stored
        if isreal(sfc) && sfc > 0
            sfc_grid(i) = sfc;
        end
    end
 
    % Actual pi_o_limit: 49.5
    % Ideal pi_o_limit: 45
 
    pi_cH_limit = 60 / (pi_fC * IPC);
 
    figure;
    contourf(pi_cH_grid, sfc_grid, B_grid, 20, 'LineColor', 'none');
    shading interp;
    h = colorbar;
    set(h, 'TickLabelInterpreter', 'latex');
    ylabel(h, '$B$', 'Interpreter', 'latex');
    xlabel('$\pi_{fC}$', 'Interpreter', 'latex');
    ylabel('$SFC$ (kgkN$^{-1}$s$^{-1}$)', 'Interpreter', 'latex');
    set(gca, "TickLabelInterpreter", 'latex');
    ylim([0.01, 0.027]);
 
    % Actual y-limit: 0.0144529
    % Ideal y-limit:0.013139
    yline(0.0144529, 'k--','Label', '$SFC_{req}$', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');             % SFC requirement
    xline(pi_cH_limit, 'k--','Label', '$\pi_{fC_{lim}}$', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');         % pi_fC limit from pi_o  
   
    if exportMode
        saveas(gcf, 'HPC_v_SFC_forB.eps', 'epsc');
    end

elseif graph == 4

    % HPC vs F_S for B
    pi_fC = 2.27;
    IPC = 2.5;
    T_o4 = 2090;
 
    % Define parameter ranges
    pi_cH = linspace(6, 12.5, 200);     % HPC range
    B = linspace(4, 11.5, 200);          % Bypass ratio range
 
    % Meshgrid for evaluation
    [pi_cH_grid, B_grid] = meshgrid(pi_cH, B);
    F_S_grid = nan(size(pi_cH_grid)); % Preallocate with NaNs
 
    % Compute SFC for each (pi_fC, B)
    % BE AWARE T_o4 toleranced
    for i = 1:numel(pi_cH_grid)
        [F_S, ~, ~, ~, ~, ~] = cycle_analyser(B_grid(i), pi_cH_grid(i), IPC, T_o4, pi_fC);
       
        % Ensure only real, positive SFC values are stored
        if isreal(F_S) && F_S > 0
            F_S_grid(i) = F_S;
        end
    end
 
    % Actual pi_o_limit: 49.5
    % Ideal pi_o_limit: 45
 
    pi_cH_limit = 60 / (pi_fC * IPC);
 
    figure;
    contourf(pi_cH_grid, F_S_grid, B_grid, 20, 'LineColor', 'none');
    shading interp;
    h = colorbar;
    set(h, 'TickLabelInterpreter', 'latex');
    ylabel(h, '$B$', 'Interpreter', 'latex');
    xlabel('$\pi_{cH}$', 'Interpreter', 'latex');
    ylabel('$F_S$ (Nskg$^{-1}$)', 'Interpreter', 'latex');
    set(gca, "TickLabelInterpreter", 'latex');
    ylim([100, 375]);
 
    % Actual y-limit: 23.6456
    % Ideal y-limit: 21.496 kN
    xline(pi_cH_limit, 'k--','Label', '$\pi_{fC_{lim}}$', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');         % pi_fC limit from pi_o  
   
    if exportMode
        saveas(gcf, 'HPC_v_FS_forB.eps', 'epsc');
    end



elseif graph == 5
    

    % IPC vs SFC for B
    pi_fC = 2.27;
    HPC = 10;
    T_o4 = 2090;
 
    % Define parameter ranges
    pi_cI = linspace(2.1, 4, 200); % Fan core pressure ratio range
    B = linspace(4, 11.5, 200);          % Bypass ratio range
 
    % Meshgrid for evaluation
    [pi_cI_grid, B_grid] = meshgrid(pi_cI, B);
    sfc_grid = nan(size(pi_cI_grid)); % Preallocate with NaNs
 
    % Compute SFC for each (pi_fC, B)
    % BE AWARE T_o4 toleranced
    for i = 1:numel(pi_cI_grid)
        [~, sfc, ~, ~, ~, ~] = cycle_analyser(B_grid(i), HPC, pi_cI_grid(i), T_o4, pi_fC);
       
        % Ensure only real, positive SFC values are stored
        if isreal(sfc) && sfc > 0
            sfc_grid(i) = sfc;
        end
    end
 
    % Actual pi_o_limit: 49.5
    % Ideal pi_o_limit: 45
 
    pi_cI_limit = 60 / (pi_fC * HPC);
 
    figure;
    contourf(pi_cI_grid, sfc_grid, B_grid, 20, 'LineColor', 'none');
    shading interp;
    h = colorbar;
    set(h, 'TickLabelInterpreter', 'latex');
    ylabel(h, '$B$', 'Interpreter', 'latex');
    xlabel('$\pi_{fI}$', 'Interpreter', 'latex');
    ylabel('$SFC$ (kgkN$^{-1}$s$^{-1}$)', 'Interpreter', 'latex');
    set(gca, "TickLabelInterpreter", 'latex');
    ylim([0.01, 0.027]);
 
    % Actual y-limit: 0.0144529
    % Ideal y-limit:0.013139
    yline(0.0144529, 'k--','Label', '$SFC_{req}$', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');             % SFC requirement
    xline(pi_cI_limit, 'k--','Label', '$\pi_{fC_{lim}}$', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');         % pi_fC limit from pi_o  
   
    if exportMode
        saveas(gcf, 'IPC_v_SFC_forB.eps', 'epsc');
    end

elseif graph == 6

    % IPC vs F_S for B
    pi_fC = 2.27;
    HPC = 10;
    T_o4 = 2090;
 
    % Define parameter ranges
    pi_cI = linspace(2.1, 4, 200); % Fan core pressure ratio range
    B = linspace(4, 11.5, 200);          % Bypass ratio range
 
    % Meshgrid for evaluation
    [pi_cI_grid, B_grid] = meshgrid(pi_cI, B);
    F_S_grid = nan(size(pi_cI_grid)); % Preallocate with NaNs
 
    % Compute SFC for each (pi_fC, B)
    % BE AWARE T_o4 toleranced
    for i = 1:numel(pi_cI_grid)
        [F_S, ~, ~, ~, ~, ~] = cycle_analyser(B_grid(i), HPC, pi_cI_grid(i), T_o4, pi_fC);
       
        % Ensure only real, positive SFC values are stored
        if isreal(F_S) && F_S > 0
            F_S_grid(i) = F_S;
        end
    end
 
    % Actual pi_o_limit: 49.5
    % Ideal pi_o_limit: 45
 
    pi_cI_limit = 60 / (pi_fC * HPC);
 
    figure;
    contourf(pi_cI_grid, F_S_grid, B_grid, 20, 'LineColor', 'none');
    shading interp;
    h = colorbar;
    set(h, 'TickLabelInterpreter', 'latex');
    ylabel(h, '$B$', 'Interpreter', 'latex');
    xlabel('$\pi_{cI}$', 'Interpreter', 'latex');
    ylabel('$F_S$ (Nskg$^{-1}$)', 'Interpreter', 'latex');
    set(gca, "TickLabelInterpreter", 'latex');
    ylim([175, 375]);
 
    % Actual y-limit: 23.6456
    % Ideal y-limit: 21.496 kN
    xline(pi_cI_limit, 'k--','Label', '$\pi_{fC_{lim}}$', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');         % pi_fC limit from pi_o  
   
    if exportMode
        saveas(gcf, 'IPC_v_FS_forB.eps', 'epsc');
    end

elseif graph == 7
    % pi_o vs SFC for B

    T_o4 = 2090;
    pi_o = linspace(33, 60, 200);
    B = linspace(4, 11.5, 200);          % Bypass ratio range
 
    % Meshgrid for evaluation
    [pi_o_grid, B_grid] = meshgrid(pi_o, B);
    sfc_grid = nan(size(pi_o_grid)); % Preallocate with NaNs
 
    % Compute SFC for each (pi_fC, B)
    % BE AWARE T_o4 toleranced
    for i = 1:numel(pi_o_grid)

        % Good values are T_o4 2090, HPC 10, IPC 2.5, pi_fC 2.27 - pi_o = 

        HPC(i) = (10 / 56.75) * pi_o_grid(i);
        IPC(i) = (2.5 / 56.75) * pi_o_grid(i);
        pi_fC(i) = (2.27 / 56.75) * pi_o_grid(i);

        [~, sfc, ~, ~, ~, ~] = cycle_analyser(B_grid(i), HPC(i), IPC(i), T_o4, pi_fC(i));
       
        % Ensure only real, positive SFC values are stored
        if isreal(sfc) && sfc > 0
            sfc_grid(i) = sfc;
        end
    end
 
    % Actual pi_o_limit: 49.5
    % Ideal pi_o_limit: 45
 
    figure;
    contourf(pi_o_grid, sfc_grid, B_grid, 20, 'LineColor', 'none');
    shading interp;
    h = colorbar;
    set(h, 'TickLabelInterpreter', 'latex');
    ylabel(h, '$B$', 'Interpreter', 'latex');
    xlabel('$\pi_{o}$', 'Interpreter', 'latex');
    ylabel('$SFC$ (kgkN$^{-1}$s$^{-1}$)', 'Interpreter', 'latex');
    set(gca, "TickLabelInterpreter", 'latex');
    ylim([0.01, 0.07]);
 
    % Actual y-limit: 0.0144529
    % Ideal y-limit:0.013139
    yline(0.0144529, 'k--','Label', '$SFC_{req}$', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');             % SFC requirement
    xline(49.5, 'k--','Label', '$\pi_{o_{lim}}$', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');         % pi_fC limit from pi_o  
   
    if exportMode
        saveas(gcf, 'OVERALL_PRESSURE_RATIO_v_SFC_forB.eps', 'epsc');
    end

elseif graph == 8

    % pi_o vs F_S for B

    T_o4 = 2090;
    pi_o = linspace(30, 60, 200);
    B = linspace(4, 11.5, 200);          % Bypass ratio range
    
    % Meshgrid for evaluation
    [pi_o_grid, B_grid] = meshgrid(pi_o, B);
    F_S_grid = nan(size(pi_o_grid)); % Preallocate with NaNs
    
    % Compute SFC for each (pi_fC, B)
    % BE AWARE T_o4 toleranced
    for i = 1:numel(pi_o_grid)

        % Good values are T_o4 2090, HPC 10, IPC 2.5, pi_fC 2.27 - pi_o = 

        HPC(i) = (10 / 56.75) * pi_o_grid(i);
        IPC(i) = (2.5 / 56.75) * pi_o_grid(i);
        pi_fC(i) = (2.27 / 56.75) * pi_o_grid(i);

        [F_S, ~, ~, ~, ~, ~] = cycle_analyser(B_grid(i), HPC(i), IPC(i), T_o4, pi_fC(i));
       
        % Ensure only real, positive SFC values are stored
        if isreal(F_S) && F_S > 0
            F_S_grid(i) = F_S;
        end
    end
 
    figure;
    contourf(pi_o_grid, F_S_grid, B_grid, 20, 'LineColor', 'none');
    shading interp;
    h = colorbar;
    set(h, 'TickLabelInterpreter', 'latex');
    ylabel(h, '$B$', 'Interpreter', 'latex');
    xlabel('$\pi_{o}$', 'Interpreter', 'latex');
    ylabel('$F_S$ (Nskg$^{-1}$)', 'Interpreter', 'latex');
    set(gca, "TickLabelInterpreter", 'latex');
    %ylim([175, 375]);
 
    xline(49.5, 'k--','Label', '$\pi_{o_{lim}}$', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex'); 
   
    if exportMode
        saveas(gcf, 'OVERALL_PR_v_FS_forB.eps', 'epsc');
    end




elseif graph == 9
    

    % T_o4 vs SFC for B
    pi_fC = 2.27;
    HPC = 10;
    IPC = 2.5;
 
    % Define parameter ranges
    T_o4 = linspace(1700, 2150, 100); % Fan core pressure ratio range
    B = linspace(4, 11.5, 100);          % Bypass ratio range
 
    % Meshgrid for evaluation
    [T_o4_grid, B_grid] = meshgrid(T_o4, B);
    sfc_grid = nan(size(T_o4_grid)); % Preallocate with NaNs
 
    % Compute SFC for each (pi_fC, B)
    % BE AWARE T_o4 toleranced
    for i = 1:numel(T_o4_grid)
        [~, sfc, ~, ~, ~, ~] = cycle_analyser(B_grid(i), HPC, IPC, T_o4_grid(i), pi_fC);
       
        % Ensure only real, positive SFC values are stored
        if isreal(sfc) && sfc > 0
            sfc_grid(i) = sfc;
        end
    end
 
    % Actual pi_o_limit: 49.5
    % Ideal pi_o_limit: 45
 
    T_o4_limit = 2090;
 
    figure;
    contourf(T_o4_grid, sfc_grid, B_grid, 20, 'LineColor', 'none');
    shading interp;
    h = colorbar;
    set(h, 'TickLabelInterpreter', 'latex');
    ylabel(h, '$B$', 'Interpreter', 'latex');
    xlabel('$T_{o4}$ (K)', 'Interpreter', 'latex');
    ylabel('$SFC$ (kgkN$^{-1}$s$^{-1}$)', 'Interpreter', 'latex');
    set(gca, "TickLabelInterpreter", 'latex');
    ylim([0.012, 0.027]);
 
    % Actual y-limit: 0.0144529
    % Ideal y-limit:0.013139
    yline(0.0144529, 'k--','Label', '$SFC_{req}$', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');             % SFC requirement
    xline(T_o4_limit, 'k--','Label', '$\pi_{fC_{lim}}$', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');         % pi_fC limit from pi_o  
   
    if exportMode
        saveas(gcf, 'TET_v_SFC_forB.eps', 'epsc');
    end

elseif graph == 10

    % T_o4 vs F_S for B
    pi_fC = 2.27;
    HPC = 10;
    IPC = 2.5;
 
    % Define parameter ranges
    T_o4 = linspace(1700, 2150, 100); % Fan core pressure ratio range
    B = linspace(4, 11.5, 100);          % Bypass ratio range
 
    % Meshgrid for evaluation
 
    % Meshgrid for evaluation
    [T_o4_grid, B_grid] = meshgrid(T_o4, B);
    F_S_grid = nan(size(T_o4_grid)); % Preallocate with NaNs
 
    % Compute SFC for each (pi_fC, B)
    % BE AWARE T_o4 toleranced
    for i = 1:numel(T_o4_grid)
        [F_S, ~, ~, ~, ~, ~] = cycle_analyser(B_grid(i), HPC, IPC, T_o4_grid(i), pi_fC);
       
        % Ensure only real, positive SFC values are stored
        if isreal(F_S) && F_S > 10 && F_S < 350
            F_S_grid(i) = F_S;
        end
    end
 
    % Actual pi_o_limit: 49.5
    % Ideal pi_o_limit: 45
 
    T_o4_limit = 2090;
 
    figure;
    contourf(T_o4_grid, F_S_grid, B_grid, 20, 'LineColor', 'none');
    shading interp;
    h = colorbar;
    set(h, 'TickLabelInterpreter', 'latex');
    ylabel(h, '$B$', 'Interpreter', 'latex');
    xlabel('$T_o4$ (K)', 'Interpreter', 'latex');
    ylabel('$F_S$ (Nskg$^{-1}$)', 'Interpreter', 'latex');
    set(gca, "TickLabelInterpreter", 'latex');
    ylim([215, 375]);
 
    % Actual y-limit: 23.6456
    % Ideal y-limit: 21.496 kN
    xline(T_o4_limit, 'k--','Label', '$\pi_{fC_{lim}}$', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');         % pi_fC limit from pi_o  
   
    if exportMode
        saveas(gcf, 'TET_v_FS_forB.eps', 'epsc');
    end



elseif graph == 11


    % pi_o vs d_in for B

    T_o4 = 2090;
    pi_o = linspace(33, 60, 200);
    B = linspace(4, 11.5, 200);          % Bypass ratio range
 
    % Meshgrid for evaluation
    [pi_o_grid, B_grid] = meshgrid(pi_o, B);
    d_in_grid = nan(size(pi_o_grid)); % Preallocate with NaNs
 
    % Compute SFC for each (pi_fC, B)
    % BE AWARE T_o4 toleranced
    for i = 1:numel(pi_o_grid)

        % Good values are T_o4 2090, HPC 10, IPC 2.5, pi_fC 2.27 - pi_o = 

        HPC(i) = (10 / 56.75) * pi_o_grid(i);
        IPC(i) = (2.5 / 56.75) * pi_o_grid(i);
        pi_fC(i) = (2.27 / 56.75) * pi_o_grid(i);

        [~, ~, ~, ~, d_in, ~] = cycle_analyser(B_grid(i), HPC(i), IPC(i), T_o4, pi_fC(i));
       
        % Ensure only real, positive SFC values are stored
        if isreal(d_in) && d_in > 0 
            d_in_grid(i) = d_in;
        end
    end
 
    % Actual pi_o_limit: 49.5
    % Ideal pi_o_limit: 45
 
    figure;
    contourf(pi_o_grid, d_in_grid, B_grid, 20, 'LineColor', 'none');
    shading interp;
    h = colorbar;
    set(h, 'TickLabelInterpreter', 'latex');
    ylabel(h, '$B$', 'Interpreter', 'latex');
    xlabel('$\pi_o$', 'Interpreter', 'latex');
    ylabel('$d_{in}$ (m)', 'Interpreter', 'latex');
    set(gca, "TickLabelInterpreter", 'latex');
    %ylim([0.01, 0.07]);
 
    % Actual y-limit: 0.0144529
    % Ideal y-limit:0.013139
    yline((2.06*1.1), 'k--','Label', '$d_{in_{PW1100G-JM}}$', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', 'Interpreter', 'latex');             % SFC requirement
    xline(56.75, 'k--','Label', '$\pi_{o_{lim}}$', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');         % pi_fC limit from pi_o  
   
    if exportMode
        saveas(gcf, 'OVERALL_PRESSURE_RATIO_v_D_in_forB.eps', 'epsc');
    end

elseif graph == 12

    % T_o4 vs d_in for B

    HPC = 10;
    IPC = 2.5;
    pi_fC = 2.27;

    T_o4 = linspace(1700, 2150, 200);
    B = linspace(4, 11.5, 200);          % Bypass ratio range
    
    % Meshgrid for evaluation
    [T_o4_grid, B_grid] = meshgrid(T_o4, B);
    d_in_grid = nan(size(T_o4_grid)); % Preallocate with NaNs
    
    % Compute SFC for each (pi_fC, B)
    % BE AWARE T_o4 toleranced
    for i = 1:numel(T_o4_grid)

        % Good values are T_o4 2090, HPC 10, IPC 2.5, pi_fC 2.27 - pi_o = 

        [~, ~, ~, ~, d_in, ~] = cycle_analyser(B_grid(i), HPC, IPC, T_o4_grid(i), pi_fC);
        
        % Ensure only real, positive SFC values are stored
        if isreal(d_in) && d_in > 1.5
            d_in_grid(i) = d_in;
        end
    end
    
    % Actual pi_o_limit: 49.5
    % Ideal pi_o_limit: 45
    
    figure;
    contourf(T_o4_grid, d_in_grid, B_grid, 20, 'LineColor', 'none');
    shading interp;
    h = colorbar;
    set(h, 'TickLabelInterpreter', 'latex');
    ylabel(h, '$B$', 'Interpreter', 'latex');
    xlabel('$T_{04}$ (K)', 'Interpreter', 'latex');
    ylabel('$d_{in}$ (m)', 'Interpreter', 'latex');
    set(gca, "TickLabelInterpreter", 'latex');
    ylim([1.5, 2.3]);
    
    % Actual y-limit: 0.0144529
    % Ideal y-limit:0.013139
    yline((2.06*1.1), 'k--','Label', '$d_{in_{PW1100G-JM}}$', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');             % SFC requirement
    xline(2090, 'k--','Label', '$T_{o4_{lim}}$', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');         % pi_fC limit from pi_o  
    
    if exportMode
        saveas(gcf, 'TET_v_D_in_forB.eps', 'epsc');
    end


elseif graph == 13


    % pi_o vs d_exit for B

    T_o4 = 2090;
    pi_o = linspace(33, 60, 200);
    B = linspace(4, 11.5, 200);          % Bypass ratio range
 
    % Meshgrid for evaluation
    [pi_o_grid, B_grid] = meshgrid(pi_o, B);
    d_exit_grid = nan(size(pi_o_grid)); % Preallocate with NaNs
 
    % Compute SFC for each (pi_fC, B)
    % BE AWARE T_o4 toleranced
    for i = 1:numel(pi_o_grid)

        % Good values are T_o4 2090, HPC 10, IPC 2.5, pi_fC 2.27 - pi_o = 

        HPC(i) = (10 / 56.75) * pi_o_grid(i);
        IPC(i) = (2.5 / 56.75) * pi_o_grid(i);
        pi_fC(i) = (2.27 / 56.75) * pi_o_grid(i);

        [~, ~, ~, ~, ~, d_exit] = cycle_analyser(B_grid(i), HPC(i), IPC(i), T_o4, pi_fC(i));
       
        % Ensure only real, positive SFC values are stored
        if isreal(d_exit) && d_exit > 0 
            d_exit_grid(i) = d_exit;
        end
    end
 
    % Actual pi_o_limit: 49.5
    % Ideal pi_o_limit: 45
 
    figure;
    contourf(pi_o_grid, d_exit_grid, B_grid, 20, 'LineColor', 'none');
    shading interp;
    h = colorbar;
    set(h, 'TickLabelInterpreter', 'latex');
    ylabel(h, '$B$', 'Interpreter', 'latex');
    xlabel('$\pi_o$', 'Interpreter', 'latex');
    ylabel('$d_{exit}$ (m)', 'Interpreter', 'latex');
    set(gca, "TickLabelInterpreter", 'latex');
    %ylim([0.01, 0.07]);
 
    % Actual y-limit: 0.0144529
    % Ideal y-limit:0.013139
    %yline((2.06*1.1), 'k--','Label', '$d_{in_{PW1100G-JM}}$', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'right', 'Interpreter', 'latex');             % SFC requirement
    xline(56.75, 'k--','Label', '$\pi_{o_{lim}}$', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');         % pi_fC limit from pi_o  
   
    if exportMode
        saveas(gcf, 'OVERALL_PRESSURE_RATIO_v_D_exit_forB.eps', 'epsc');
    end

elseif graph == 14

    % T_o4 vs d_in for B

    HPC = 10;
    IPC = 2.5;
    pi_fC = 2.27;

    T_o4 = linspace(1700, 2150, 200);
    B = linspace(4, 11.5, 200);          % Bypass ratio range
    
    % Meshgrid for evaluation
    [T_o4_grid, B_grid] = meshgrid(T_o4, B);
    d_exit_grid = nan(size(T_o4_grid)); % Preallocate with NaNs
    
    % Compute SFC for each (pi_fC, B)
    % BE AWARE T_o4 toleranced
    for i = 1:numel(T_o4_grid)

        % Good values are T_o4 2090, HPC 10, IPC 2.5, pi_fC 2.27 - pi_o = 

        [~, ~, ~, ~, ~, d_exit] = cycle_analyser(B_grid(i), HPC, IPC, T_o4_grid(i), pi_fC);
        
        % Ensure only real, positive SFC values are stored
        if isreal(d_exit) && d_exit > 0.89
            d_exit_grid(i) = d_exit;
        end
    end
    
    % Actual pi_o_limit: 49.5
    % Ideal pi_o_limit: 45
    
    figure;
    contourf(T_o4_grid, d_exit_grid, B_grid, 20, 'LineColor', 'none');
    shading interp;
    h = colorbar;
    set(h, 'TickLabelInterpreter', 'latex');
    ylabel(h, '$B$', 'Interpreter', 'latex');
    xlabel('$T_{04}$ (K)', 'Interpreter', 'latex');
    ylabel('$d_{exit}$ (m)', 'Interpreter', 'latex');
    set(gca, "TickLabelInterpreter", 'latex');
    %ylim([1.5, 2.3]);
    
    % Actual y-limit: 0.0144529
    % Ideal y-limit:0.013139
    %yline((2.06*1.1), 'k--','Label', '$d_{in_{PW1100G-JM}}$', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');             % SFC requirement
    xline(2090, 'k--','Label', '$T_{o4_{lim}}$', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');         % pi_fC limit from pi_o  
    
    if exportMode
        saveas(gcf, 'TET_v_D_exit_forB.eps', 'epsc');
    end

elseif graph == 15


    % pi_o vs V_19_9 for B

    T_o4 = 2090;
    pi_o = linspace(33, 60, 200);
    B = linspace(4, 11.5, 200);          % Bypass ratio range
 
    % Meshgrid for evaluation
    [pi_o_grid, B_grid] = meshgrid(pi_o, B);
    V_19_9_grid = nan(size(pi_o_grid)); % Preallocate with NaNs
 
    % Compute SFC for each (pi_fC, B)
    % BE AWARE T_o4 toleranced
    for i = 1:numel(pi_o_grid)

        % Good values are T_o4 2090, HPC 10, IPC 2.5, pi_fC 2.27 - pi_o = 

        HPC(i) = (10 / 56.75) * pi_o_grid(i);
        IPC(i) = (2.5 / 56.75) * pi_o_grid(i);
        pi_fC(i) = (2.27 / 56.75) * pi_o_grid(i);

        [~, ~, ~, ~, ~, ~, V_19_9,~] = cycle_analyser2(B_grid(i), HPC(i), IPC(i), T_o4, pi_fC(i));
       
        % Ensure only real, positive SFC values are stored
        if isreal(V_19_9) && V_19_9 > 0 
            V_19_9_grid(i) = V_19_9;
        end
    end
 
    % Actual pi_o_limit: 49.5
    % Ideal pi_o_limit: 45
 
    figure;
    contourf(pi_o_grid, V_19_9_grid, B_grid, 20, 'LineColor', 'none');
    shading interp;
    h = colorbar;
    set(h, 'TickLabelInterpreter', 'latex');
    ylabel(h, '$B$', 'Interpreter', 'latex');
    xlabel('$\pi_o$', 'Interpreter', 'latex');
    ylabel('$V_{19}/V_9$', 'Interpreter', 'latex');
    set(gca, "TickLabelInterpreter", 'latex');
    ylim([0.29, 0.9]);
 
    yline(0.65, 'k--','Label', '$V_{19}/V_9$ low', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');             % SFC requirement
    yline(0.85, 'k--','Label', '$V_{19}/V_9 $ high', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');             % SFC requirement
    xline(56.75, 'k--','Label', '$\pi_{o_{lim}}$', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');         % pi_fC limit from pi_o  
   
    if exportMode
        saveas(gcf, 'OVERALL_PRESSURE_RATIO_v_V199_forB.eps', 'epsc');
    end

elseif graph == 16

    % T_o4 vs d_in for B

    HPC = 10;
    IPC = 2.5;
    pi_fC = 2.27;

    T_o4 = linspace(1700, 2150, 200);
    B = linspace(4, 11.5, 200);          % Bypass ratio range
    
    % Meshgrid for evaluation
    [T_o4_grid, B_grid] = meshgrid(T_o4, B);
    V_19_9_grid = nan(size(T_o4_grid)); % Preallocate with NaNs
    
    % Compute SFC for each (pi_fC, B)
    % BE AWARE T_o4 toleranced
    for i = 1:numel(T_o4_grid)

        % Good values are T_o4 2090, HPC 10, IPC 2.5, pi_fC 2.27 - pi_o = 

        [~, ~, ~, ~, ~, ~, V_19_9,~] = cycle_analyser2(B_grid(i), HPC, IPC, T_o4_grid(i), pi_fC);
        
        % Ensure only real, positive SFC values are stored
        if isreal(V_19_9) && V_19_9 > 0.6
            V_19_9_grid(i) = V_19_9;
        end
    end
    
    % Actual pi_o_limit: 49.5
    % Ideal pi_o_limit: 45
    
    figure;
    contourf(T_o4_grid, V_19_9_grid, B_grid, 20, 'LineColor', 'none');
    shading interp;
    h = colorbar;
    set(h, 'TickLabelInterpreter', 'latex');
    ylabel(h, '$B$', 'Interpreter', 'latex');
    xlabel('$T_{04}$ (K)', 'Interpreter', 'latex');
    ylabel('$V_{19}/V_9$', 'Interpreter', 'latex');
    set(gca, "TickLabelInterpreter", 'latex');
    ylim([0.6, 0.9]);
    
    % Actual y-limit: 0.0144529
    % Ideal y-limit:0.013139
    yline(0.65, 'k--','Label', '$V_{19}/V_9$ low', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');             % SFC requirement
    yline(0.85, 'k--','Label', '$V_{19}/V_9 $ high', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'center', 'Interpreter', 'latex');             % SFC requirement
    xline(2090, 'k--','Label', '$T_{o4_{lim}}$', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');         % pi_fC limit from pi_o  
    
    if exportMode
        saveas(gcf, 'TET_v_V199_forB.eps', 'epsc');
    end

end
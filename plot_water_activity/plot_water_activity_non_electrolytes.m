function plot_water_activity_non_electrolytes()
% Plot water activity for NON-ELECTROLYTES (Feb 2026)
% Organic compounds, sugars, and alcohols
    close all 
    clear
    clc 
    
    % --- 1. SETUP & PATHS ---
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(filepath, '..', 'calculate_mf'));
    addpath(fullfile(filepath, '..', 'util'));
    addpath(fullfile(filepath, '..', 'data'));
    addpath(fullfile(filepath, '..', 'plotting')); % For PlotDefaults_Slides
    
    fig_out_dir = fullfile(filepath, '..', 'figures', 'non_electrolytes');
    if ~exist(fig_out_dir, 'dir')
        mkdir(fig_out_dir);
    end
    
    T = 25; 
    MWw = 18.015; % g/mol for water
    
    % --- 2. COMPOUND DEFINITIONS ---
    salt_data = load_salt_data();
    
    % Filter to only the NON-ELECTROLYTES
    non_electrolytes = {'Glycerol', 'Sucrose', 'Glucose', 'Urea', 'Fructose', ...
                        'Ethanol', 'Methanol', 'Xylose', 'Arabinose', 'Xylitol'};
    
    keep = cellfun(@(r) any(strcmp(r{1}, non_electrolytes)), salt_data);
    salt_data = salt_data(keep);
    
    % --- 3. DATA PROCESSING LOOP ---
    num_points = 100;
    all_compounds_data = struct();
    
    fprintf('Processing %d non-electrolytes...\n', length(salt_data));
    
    for s = 1:length(salt_data)
        compound_name = salt_data{s}{1};
        MW = salt_data{s}{2};
        RH_min = salt_data{s}{3};
        RH_max = salt_data{s}{4};
        func_name = salt_data{s}{5};
        func_args = salt_data{s}{6};
        
        RH_vec = linspace(RH_min + 0.001, RH_max - 0.001, num_points);
        mf_compound = zeros(size(RH_vec));
        mf_water = zeros(size(RH_vec));
        molality = zeros(size(RH_vec));
        x_water = zeros(size(RH_vec));
        
        success = true;
        for i = 1:length(RH_vec)
            try
                if func_args == 1
                    mf_compound(i) = feval(func_name, RH_vec(i), T);
                else
                    mf_compound(i) = feval(func_name, RH_vec(i));
                end
                mf_water(i) = 1 - mf_compound(i);
                
                % Molality: moles compound per kg water
                molality(i) = (mf_compound(i) / MW) / (mf_water(i) / 1000);
                
                % Mole fraction (no dissociation)
                n_w = mf_water(i) / MWw;
                n_c = mf_compound(i) / MW;
                x_water(i) = n_w / (n_w + n_c);
                
            catch ME
                warning('Error processing %s at RH=%.4f: %s', compound_name, RH_vec(i), ME.message);
                success = false;
                break;
            end
        end
        
        if success
            % Calculate Activity Coefficients
            gamma_w = RH_vec ./ x_water;
            
            % Store data
            all_compounds_data.(compound_name) = struct();
            all_compounds_data.(compound_name).RH = RH_vec;
            all_compounds_data.(compound_name).mf_compound = mf_compound;
            all_compounds_data.(compound_name).molality = molality;
            all_compounds_data.(compound_name).x_water = x_water;
            all_compounds_data.(compound_name).gamma_w = gamma_w;
            all_compounds_data.(compound_name).display_name = compound_name;
            fprintf('  Processed: %s\n', compound_name);
        end
    end
    
    compound_names = fieldnames(all_compounds_data);
    num_compounds = length(compound_names);
    fprintf('Successfully processed %d non-electrolytes\n', num_compounds);
    
    % --- Color Generation ---
    colors = generate_colors(num_compounds);
    
    %% --- PLOT 1: Water Activity vs Mass Fraction ---
    figure('Position', [100, 100, 1200, 800]);
    hold on; grid on; box on;
    
    for s = 1:num_compounds
        compound_name = compound_names{s};
        data = all_compounds_data.(compound_name);
        
        plot(data.mf_compound, data.RH, 'LineWidth', 3, ...
             'color', colors(s,:), 'DisplayName', data.display_name);
    end
    
    ax = gca;
    ax.FontSize = 20;
    xlabel('Mass Fraction of Compound', 'FontSize', 26, 'FontWeight', 'bold');
    ylabel('Water Activity (RH)', 'FontSize', 26, 'FontWeight', 'bold');
    title('Non-Electrolytes: Water Activity vs Mass Fraction', ...
          'FontSize', 28, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 14);
    
    set(gcf, 'color', 'w');
    print(fullfile(fig_out_dir, 'water_activity_vs_mf_non_electrolytes'), '-dtiff', '-r600');
    fprintf('Saved: water_activity_vs_mf_non_electrolytes.tif\n');
    
    %% --- PLOT 2: Water Activity vs Molality ---
    figure('Position', [100, 100, 1200, 800]);
    hold on; grid on; box on;
    
    for s = 1:num_compounds
        compound_name = compound_names{s};
        data = all_compounds_data.(compound_name);
        
        plot(data.molality, data.RH, 'LineWidth', 3, ...
             'color', colors(s,:), 'DisplayName', data.display_name);
    end
    
    ax = gca;
    ax.FontSize = 20;
    xlabel('Molality (mol/kg water)', 'FontSize', 26, 'FontWeight', 'bold');
    ylabel('Water Activity (RH)', 'FontSize', 26, 'FontWeight', 'bold');
    title('Non-Electrolytes: Water Activity vs Molality', ...
          'FontSize', 28, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 14);
    
    set(gcf, 'color', 'w');
    print(fullfile(fig_out_dir, 'water_activity_vs_molality_non_electrolytes'), '-dtiff', '-r600');
    fprintf('Saved: water_activity_vs_molality_non_electrolytes.tif\n');
    
    %% --- PLOT 3: Activity Coefficient vs Molality ---
    figure('Position', [100, 100, 1200, 800]);
    hold on; grid on; box on;
    
    for s = 1:num_compounds
        compound_name = compound_names{s};
        data = all_compounds_data.(compound_name);
        
        plot(data.molality, data.gamma_w, 'LineWidth', 3, ...
             'color', colors(s,:), 'DisplayName', data.display_name);
    end
    
    % Add ideal reference line (gamma = 1)
    xlim_vals = xlim;
    plot(xlim_vals, [1 1], 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal (\gamma_w = 1)');
    
    ax = gca;
    ax.FontSize = 20;
    xlabel('Molality (mol/kg water)', 'FontSize', 26, 'FontWeight', 'bold');
    ylabel('\gamma_w (Molecular Basis)', 'FontSize', 26, 'FontWeight', 'bold');
    title('Non-Electrolytes: Water Activity Coefficient vs Molality', ...
          'FontSize', 28, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 14);
    
    set(gcf, 'color', 'w');
    print(fullfile(fig_out_dir, 'gamma_vs_molality_non_electrolytes'), '-dtiff', '-r600');
    fprintf('Saved: gamma_vs_molality_non_electrolytes.tif\n');
    
end

%% Helper function for color generation
function colors = generate_colors(num_compounds)
    colors = zeros(num_compounds, 3);
    base_colors = [
        0.0000, 0.4470, 0.7410;  % Blue
        0.8500, 0.3250, 0.0980;  % Orange
        0.9290, 0.6940, 0.1250;  % Yellow
        0.4940, 0.1840, 0.5560;  % Purple
        0.4660, 0.6740, 0.1880;  % Green
        0.3010, 0.7450, 0.9330;  % Cyan
        0.6350, 0.0780, 0.1840;  % Dark Red
        0.0000, 0.5000, 0.0000;  % Dark Green
        0.7500, 0.0000, 0.7500;  % Magenta
    ];
    
    for i = 1:num_compounds
        if i <= size(base_colors, 1)
            colors(i,:) = base_colors(i,:);
        else
            hue = mod((i-1) * 0.38196601125, 1.0);
            colors(i,:) = hsv2rgb([hue, 0.85, 0.85]);
        end
    end
end

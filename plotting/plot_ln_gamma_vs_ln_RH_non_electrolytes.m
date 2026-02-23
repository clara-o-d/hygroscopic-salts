% Plot ln(gamma) vs ln(RH) for NON-ELECTROLYTES (Feb 2026)
% Organic compounds, sugars, and alcohols
close all 
clear
clc 

% Add necessary folders to path
[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, '..', 'calculate_mf'));
addpath(fullfile(filepath, '..', 'util'));
addpath(fullfile(filepath, '..', 'data'));
addpath(filepath); % For PlotDefaults_Slides

% Apply plot defaults
PlotDefaults_Slides();

% Override font size defaults for text objects
set(0,'defaultTextFontSize',35);

% Define output directory for figures
fig_out_dir = fullfile(filepath, '..', 'figures', 'non_electrolytes');
if ~exist(fig_out_dir, 'dir')
    mkdir(fig_out_dir);
end

T = 25; 
MWw = 18.015; % g/mol

% Load canonical salt data
salt_data = load_salt_data();

% Filter to only the NON-ELECTROLYTES
non_electrolytes = {'Glycerol', 'Sucrose', 'Glucose', 'Urea', 'Fructose', ...
                    'Ethanol', 'Methanol', 'Xylose', 'Arabinose', 'Xylitol'};

keep = cellfun(@(r) any(strcmp(r{1}, non_electrolytes)), salt_data);
salt_data = salt_data(keep);

% Process each compound
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
    
    % For non-electrolytes, nu = 0 (no dissociation)
    nu = 0;
    
    % Create RH vector
    RH_vec = linspace(RH_min + 0.001, RH_max - 0.001, num_points);
    
    % Initialize arrays
    mf_compound = zeros(size(RH_vec));
    mf_water = zeros(size(RH_vec));
    x_water = zeros(size(RH_vec));
    
    % Calculate for each RH value
    success = true;
    for i = 1:length(RH_vec)
        try
            if func_args == 1
                mf_compound(i) = feval(func_name, RH_vec(i), T);
            else
                mf_compound(i) = feval(func_name, RH_vec(i));
            end
            mf_water(i) = 1 - mf_compound(i);
            
            % Moles of water and compound (per unit mass basis)
            n_w = mf_water(i) / MWw;
            n_c = mf_compound(i) / MW;
            
            % Mole fraction (no ionic dissociation for non-electrolytes)
            x_water(i) = n_w / (n_w + n_c);
            
        catch ME
            warning('Error processing %s at RH=%.4f: %s', compound_name, RH_vec(i), ME.message);
            success = false;
            break;
        end
    end
    
    if success
        % Calculate Activity Coefficients (gamma = a_w / x_w)
        gamma_w = RH_vec ./ x_water;
        
        % Calculate logarithmic quantities
        ln_gamma_w = log(gamma_w);
        ln_RH = log(RH_vec);
        
        % Store data
        all_compounds_data.(compound_name) = struct();
        all_compounds_data.(compound_name).RH = RH_vec;
        all_compounds_data.(compound_name).ln_RH = ln_RH;
        all_compounds_data.(compound_name).x_water = x_water;
        all_compounds_data.(compound_name).gamma_w = gamma_w;
        all_compounds_data.(compound_name).ln_gamma_w = ln_gamma_w;
        all_compounds_data.(compound_name).display_name = compound_name;
        fprintf('  Processed: %s\n', compound_name);
    end
end

compound_names = fieldnames(all_compounds_data);
num_compounds = length(compound_names);
fprintf('Successfully processed %d non-electrolytes\n', num_compounds);

% --- Color Generation ---
colors = generate_colors(num_compounds);

%% --- PLOT: ln(Activity Coefficient) vs ln(RH) ---
figure('Position', [100, 100, 1200, 800]);
hold on; grid on; box on;

for s = 1:num_compounds
    compound_name = compound_names{s};
    data = all_compounds_data.(compound_name);
    
    % Plot ln(gamma) vs ln(RH)
    plot(data.ln_RH, data.ln_gamma_w, 'LineWidth', 3, ...
         'color', colors(s,:), 'DisplayName', data.display_name);
    
    % Label at first point (lowest RH, most negative ln(RH))
    text(data.ln_RH(1), data.ln_gamma_w(1), ['  ' data.display_name], ...
         'Color', colors(s,:), 'FontSize', 18, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'left');
end

% Add ideal reference line (ln(gamma) = 0)
xlim_vals = xlim;
plot(xlim_vals, [0 0], 'k--', 'LineWidth', 2, 'DisplayName', 'Ideal (ln(\gamma_w) = 0)');

ax = gca;
ax.FontSize = 20;

h_xlabel = xlabel('ln(RH) = ln(a_w)');
h_ylabel = ylabel('ln(\gamma_w) (Molecular Basis)');
h_title = title('Non-Electrolytes: ln(Water Activity Coefficient) vs ln(RH)');

h_xlabel.FontSize = 26;
h_xlabel.FontWeight = 'bold';
h_ylabel.FontSize = 26;
h_ylabel.FontWeight = 'bold';
h_title.FontSize = 28;
h_title.FontWeight = 'bold';

legend('Location', 'best', 'FontSize', 14);

set(gcf, 'color', 'w');

% Set scaling factor for print
alpha = 0.7;
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 5*alpha 4*alpha];

% Save figure
print(fullfile(fig_out_dir, 'ln_Activity_Coefficient_vs_ln_RH_non_electrolytes'), '-dtiff', '-r600');
fprintf('Saved: %s\n', fullfile(fig_out_dir, 'ln_Activity_Coefficient_vs_ln_RH_non_electrolytes.tif'));

%% Helper function for color generation
function colors = generate_colors(num_compounds)
    colors = zeros(num_compounds, 3);
    % Use distinct colors
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
            % Fall back to golden angle if more colors needed
            hue = mod((i-1) * 0.38196601125, 1.0);
            colors(i,:) = hsv2rgb([hue, 0.85, 0.85]);
        end
    end
end

function plot_pitzer_critical_mixes()
    % PLOT_PITZER_CRITICAL_MIXES
    % Plots water activity for multi-salt systems and thermodynamic properties
    % for binary systems.
    close all 
    clear
    clc 

    % --- 1. LOAD PARAMETERS ---
    filename = 'pitzer_parameters_phreeqc.xlsx';
    if exist(filename, 'file') ~= 2
        error('File %s not found. Run the Python script first.', filename);
    end

    fprintf('Loading Pitzer parameters...\n');
    opts = detectImportOptions(filename);
    opts.VariableTypes = repmat({'string'}, 1, length(opts.VariableTypes));
    opts.VariableTypes(5:end) = {'double'}; 

    % Use 'preserve' to keep spaces in column names like "Species 1"
    param_table = readtable(filename, 'Sheet', 'Coefficients', 'VariableNamingRule', 'preserve');
    species_table = readtable(filename, 'Sheet', 'Species_Charges', 'VariableNamingRule', 'preserve');

    % Build Charge Map & Safety Defaults
    charge_map = containers.Map();
    for i = 1:height(species_table)
        charge_map(char(species_table.Species(i))) = double(species_table.Charge(i));
    end
    defaults = {'Na+',1; 'K+',1; 'Li+',1; 'H+',1; 'Cl-',-1; 'Br-',-1; 'OH-',-1; ...
                'I-',-1; 'Ca+2',2; 'Mg+2',2; 'SO4-2',-2; 'CO3-2',-2; 'HCO3-',-1};
    for i = 1:size(defaults,1)
        if ~isKey(charge_map, defaults{i,1}), charge_map(defaults{i,1}) = defaults{i,2}; end
    end

    % Constants
    T = 298.15; % K
    R = 8.31446; % J/(mol K)
    Mw = 0.0180153; % kg/mol
    A_phi = 0.392; 

    % --- 2. DEFINE CRITICAL CASES (Existing Code) ---
    cases(1).name = 'Mg(HCO3)2 + MgSO4'; cases(1).saltA = {'Mg+2', 1; 'HCO3-', 2}; cases(1).saltB = {'Mg+2', 1; 'SO4-2', 1}; cases(1).note = 'Large \Psi_{HCO3, SO4, Mg} = -0.161';
    cases(2).name = 'MgCl2 + Mg(HCO3)2'; cases(2).saltA = {'Mg+2', 1; 'Cl-', 2}; cases(2).saltB = {'Mg+2', 1; 'HCO3-', 2}; cases(2).note = 'Large \Psi_{Cl, HCO3, Mg} = -0.096';
    cases(3).name = 'NaCl + CaCl2'; cases(3).saltA = {'Na+', 1; 'Cl-', 1}; cases(3).saltB = {'Ca+2', 1; 'Cl-', 2}; cases(3).note = 'Large \theta_{Ca, Na} = 0.092';
    cases(4).name = 'NaCl + MgSO4'; cases(4).saltA = {'Na+', 1; 'Cl-', 1}; cases(4).saltB = {'Mg+2', 1; 'SO4-2', 1}; cases(4).note = 'Simultaneous Cation & Anion Exchange';


% --- 3. EXECUTE & PLOT CRITICAL MIXES (ROBUST GRID) ---
    fprintf('Plotting critical mixes...\n');
    molalities = 0.1:0.1:4.0; 
    ratios = [0.0, 0.25, 0.50, 0.75, 1.0];

    % Dynamic Grid Calculation
    n_cases = length(cases);
    n_cols = ceil(sqrt(n_cases));
    n_rows = ceil(n_cases / n_cols);

    figure('Color', 'w', 'Name', 'Critical Mixes Water Activity', 'Position', [50, 50, 1200, 800]);
    
    for c = 1:n_cases
        current_case = cases(c);
        
        % Use dynamic row/col counts instead of hardcoded (2,2)
        subplot(n_rows, n_cols, c); 
        hold on; grid on; box on;
        
        colors = parula(length(ratios));
        for r_idx = 1:length(ratios)
            y = ratios(r_idx); aw_curve = zeros(size(molalities));
            for m_idx = 1:length(molalities)
                m_tot = molalities(m_idx); comp = containers.Map();
                for ion_i = 1:size(current_case.saltA, 1), add_molality(comp, current_case.saltA{ion_i, 1}, m_tot * (1-y) * current_case.saltA{ion_i, 2}); end
                for ion_i = 1:size(current_case.saltB, 1), add_molality(comp, current_case.saltB{ion_i, 1}, m_tot * y * current_case.saltB{ion_i, 2}); end
                phi = calculate_pitzer_phi(comp, T, param_table, charge_map, A_phi);
                sum_m = sum(cell2mat(values(comp)));
                aw_curve(m_idx) = exp( -phi * sum_m * Mw );
            end
            plot(molalities, aw_curve, 'LineWidth', 2, 'Color', colors(r_idx, :), 'DisplayName', sprintf('y_{SaltB} = %.2f', y));
        end
        title(sprintf('%s\n(%s)', current_case.name, current_case.note), 'Interpreter', 'tex');
        xlabel('Total Salt Molality (m)'); ylabel('Water Activity (a_w)');
        
        % Only show legend on the first plot to prevent clutter
        if c == 1, legend('Location', 'SouthWest'); end
    end


% --- 4. NEW SECTION: BINARY SALT ENTHALPY & ENTROPY OF MIXING ---
    fprintf('Calculating thermodynamic properties for binary salts...\n');

    % Format: { 'Name', {Components}, Molality_Limit }
    binary_salts = {
        % --- Chlorides ---
        'NaCl',   {'Na+', 1; 'Cl-', 1}, 6.0;
        'KCl',    {'K+', 1; 'Cl-', 1}, 4.5;
        'LiCl',   {'Li+', 1; 'Cl-', 1}, 10.0;
        'CaCl2',  {'Ca+2', 1; 'Cl-', 2}, 5.0;
        'MgCl2',  {'Mg+2', 1; 'Cl-', 2}, 4.5;
        'SrCl2',  {'Sr+2', 1; 'Cl-', 2}, 3.0;
        'BaCl2',  {'Ba+2', 1; 'Cl-', 2}, 1.5;

        % --- Bromides ---
        'NaBr',   {'Na+', 1; 'Br-', 1}, 7.0;
        'KBr',    {'K+', 1; 'Br-', 1}, 5.0;
        'LiBr',   {'Li+', 1; 'Br-', 1}, 8.0;
        'MgBr2',  {'Mg+2', 1; 'Br-', 2}, 4.5;
        'CaBr2',  {'Ca+2', 1; 'Br-', 2}, 5.0;
        'SrBr2',  {'Sr+2', 1; 'Br-', 2}, 3.5;

        % --- Sulfates ---
        'Na2SO4', {'Na+', 2; 'SO4-2', 1}, 2.5;
        'K2SO4',  {'K+', 2; 'SO4-2', 1}, 0.6;
        'Li2SO4', {'Li+', 2; 'SO4-2', 1}, 2.5;
        'MgSO4',  {'Mg+2', 1; 'SO4-2', 1}, 2.0;

        % --- Carbonates & Bicarbonates ---
        'Na2CO3', {'Na+', 2; 'CO3-2', 1}, 2.5;
        'K2CO3',  {'K+', 2; 'CO3-2', 1}, 6.0;
        'NaHCO3', {'Na+', 1; 'HCO3-', 1}, 1.0;
        'KHCO3',  {'K+', 1; 'HCO3-', 1}, 3.0;

        % --- Hydroxides ---
        'NaOH',   {'Na+', 1; 'OH-', 1}, 15.0;
        'KOH',    {'K+', 1; 'OH-', 1}, 15.0;
        'LiOH',   {'Li+', 1; 'OH-', 1}, 4.0
    };
    
    % Three separate figures (one per plot)
    figure('Color', 'w', 'Name', 'Excess Entropy', 'Position', [50, 50, 700, 450]);
    ax1 = gca(); hold on; grid on; box on;
    ylabel('$T \Delta \bar{S}_w^{excess}$ (J/mol)','Interpreter','latex');
    title(['$T \Delta \bar{S}_w^{excess}$ (ideal $-RT\ln x_w$ subtracted) at T = ', sprintf('%.1f K', T)], 'Interpreter','latex');
    xlabel('Relative Humidity (%)');
    
    figure('Color', 'w', 'Name', 'Enthalpy of Mixing', 'Position', [100, 100, 700, 450]);
    ax2 = gca(); hold on; grid on; box on;
    ylabel('$\Delta \bar{H}_w$ (J/mol)', 'Interpreter','latex');
    title('Enthalpy of Mixing (\Delta\bar{H}_w)', Interpreter="tex");
    xlabel('Relative Humidity (%)');

    figure('Color', 'w', 'Name', 'Activity Coefficient', 'Position', [150, 150, 700, 450]);
    ax3 = gca(); hold on; grid on; box on;
    ylabel('$\gamma_w = a_w / x_w$', 'Interpreter','latex');
    xlabel('Relative Humidity (%)'); 
    title('Activity Coefficient (\gamma_w)', Interpreter="tex");

    % Three additional figures with x_w as x-axis
    figure('Color', 'w', 'Name', 'Excess Entropy vs x_w', 'Position', [200, 50, 700, 450]);
    ax4 = gca(); hold on; grid on; box on;
    ylabel('$T \Delta \bar{S}_w^{excess}$ (J/mol)','Interpreter','latex');
    title(['$T \Delta \bar{S}_w^{excess}$ vs $x_w$ at T = ', sprintf('%.1f K', T)], 'Interpreter','latex');
    xlabel('Water mole fraction ($x_w$)','Interpreter','latex');
    
    figure('Color', 'w', 'Name', 'Enthalpy of Mixing vs x_w', 'Position', [250, 100, 700, 450]);
    ax5 = gca(); hold on; grid on; box on;
    ylabel('$\Delta \bar{H}_w$ (J/mol)', 'Interpreter','latex');
    title('Enthalpy of Mixing vs $x_w$', Interpreter="tex");
    xlabel('Water mole fraction ($x_w$)','Interpreter','latex');

    figure('Color', 'w', 'Name', 'Activity Coefficient vs x_w', 'Position', [300, 150, 700, 450]);
    ax6 = gca(); hold on; grid on; box on;
    ylabel('$\gamma_w = a_w / x_w$', 'Interpreter','latex');
    xlabel('Water mole fraction ($x_w$)','Interpreter','latex');
    title('Activity Coefficient vs $x_w$', Interpreter="tex");

    % ln(activity) vs ln(x_w), ΔH̄_w/RT, ΔS̄_w^excess/R
    figure('Color', 'w', 'Name', 'ln(a_w) vs thermodynamic quantities', 'Position', [350, 200, 700, 450]);
    ax7 = gca(); hold on; grid on; box on;
    xlabel('$\ln(a_w)$', 'Interpreter','latex');
    ylabel('ln(x_w), $\Delta\bar{H}_w/RT$, $\Delta\bar{S}_w^{excess}/R$', 'Interpreter','latex');
    title('$\ln(x_w)$, $\Delta\bar{H}_w/RT$, $\Delta\bar{S}_w^{excess}/R$ vs $\ln(a_w)$', 'Interpreter','latex');
    
    colors = lines(size(binary_salts, 1));
    
    for s_idx = 1:size(binary_salts, 1)
        salt_name = binary_salts{s_idx, 1};
        salt_comp = binary_salts{s_idx, 2};
        max_m     = binary_salts{s_idx, 3};
        
        m_range = 0.1:0.1:max_m;
        rh_vals = []; T_dS_vals = []; dH_vals = []; xw_vals = []; aw_vals = [];
        
        for m_i = 1:length(m_range)
            m = m_range(m_i);
            comp = containers.Map();
            for ion_i = 1:size(salt_comp, 1)
                add_molality(comp, salt_comp{ion_i, 1}, m * salt_comp{ion_i, 2});
            end
            sum_m = sum(cell2mat(values(comp)));

            try
                phi = calculate_pitzer_phi(comp, T, param_table, charge_map, A_phi);
                dphi_dT = calculate_pitzer_dphi_dT(comp, T, param_table, charge_map);

                aw = exp(-phi * sum_m * Mw);
                
                % Thermodynamics
                % T*ΔS̄_w (total) and ΔH̄_w from Pitzer; ideal partial entropy: -R ln(x_w)
                T_S_mix = T * R * Mw * sum_m * (phi + T * dphi_dT);
                H_mix   = R * T^2 * Mw * sum_m * dphi_dT; 
                x_w = (1/Mw) / ( (1/Mw) + sum_m );

                rh_vals(end+1) = aw * 100;
                aw_vals(end+1) = aw;
                T_dS_vals(end+1) = T_S_mix;   %#ok<AGROW>  % Total (for activity/gamma_w)
                dH_vals(end+1) = H_mix;
                xw_vals(end+1) = x_w;     
            catch
                continue;
            end
        end
        
        if isempty(rh_vals), continue; end

        % Plot Lines (ax1: excess entropy = total - ideal = T*ΔS̄_w + R*T*ln(x_w))
        T_S_excess = T_dS_vals + R*T*log(xw_vals);
        plot(ax1, rh_vals, T_S_excess, 'LineWidth', 2, 'Color', colors(s_idx,:));
        plot(ax2, rh_vals, dH_vals, 'LineWidth', 2, 'Color', colors(s_idx,:));
        
        activity_calc = exp((dH_vals - T_dS_vals) ./ (R*T));
        gamma_w = activity_calc ./ xw_vals;
        plot(ax3, rh_vals, gamma_w, 'LineWidth', 2, 'Color', colors(s_idx,:));

        % Same plots with x_w as x-axis
        plot(ax4, xw_vals, T_S_excess, 'LineWidth', 2, 'Color', colors(s_idx,:));
        plot(ax5, xw_vals, dH_vals, 'LineWidth', 2, 'Color', colors(s_idx,:));
        plot(ax6, xw_vals, gamma_w, 'LineWidth', 2, 'Color', colors(s_idx,:));

        % ln(activity) vs ln(x_w), ΔH̄_w/RT, ΔS̄_w^excess/R
        ln_aw = log(aw_vals);
        ln_xw = log(xw_vals);
        dH_dimless = dH_vals / (R*T);           % partial enthalpy / RT
        S_excess = T_S_excess / T;              % partial excess entropy in J/(mol K)
        dS_dimless = S_excess / R;              % partial excess entropy / R
        if s_idx == 1
            plot(ax7, ln_aw, ln_xw, 'LineWidth', 1.5, 'Color', [0 0.45 0.74], 'DisplayName', 'ln(x_w)');
            plot(ax7, ln_aw, dH_dimless, 'LineWidth', 1.5, 'Color', [0.85 0.33 0.1], 'DisplayName', '\Delta\bar{H}_w/RT');
            plot(ax7, ln_aw, dS_dimless, 'LineWidth', 1.5, 'Color', [0.0 0.0 0.0], 'DisplayName', '\Delta\bar{S}_w^{excess}/R');
        else
            plot(ax7, ln_aw, ln_xw, 'LineWidth', 1.5, 'Color', [0 0.45 0.74], 'HandleVisibility', 'off');
            plot(ax7, ln_aw, dH_dimless, 'LineWidth', 1.5, 'Color', [0.85 0.33 0.1], 'HandleVisibility', 'off');
            plot(ax7, ln_aw, dS_dimless, 'LineWidth', 1.5, 'Color', [0.0 0.0 0.0], 'HandleVisibility', 'off');
        end
        
        % --- LABELS (Rightmost limit = last point in array) ---
        % 'Right' alignment means the text ends at the point, sitting to the left of the line
        text(ax1, rh_vals(end), T_S_excess(end), [' ' salt_name ' '], ...
            'Color', 'k', 'FontSize', 8, 'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'middle', 'BackgroundColor', 'w', 'EdgeColor', 'none', 'Margin', 0.5);
            
        text(ax2, rh_vals(end), dH_vals(end), [' ' salt_name ' '], ...
            'Color', 'k', 'FontSize', 8, 'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'middle', 'BackgroundColor', 'w', 'EdgeColor', 'none', 'Margin', 0.5);
            
        text(ax3, rh_vals(end), gamma_w(end), [' ' salt_name ' '], ...
            'Color', 'k', 'FontSize', 8, 'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'middle', 'BackgroundColor', 'w', 'EdgeColor', 'none', 'Margin', 0.5);

        text(ax4, xw_vals(end), T_S_excess(end), [' ' salt_name ' '], ...
            'Color', 'k', 'FontSize', 8, 'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'middle', 'BackgroundColor', 'w', 'EdgeColor', 'none', 'Margin', 0.5);
            
        text(ax5, xw_vals(end), dH_vals(end), [' ' salt_name ' '], ...
            'Color', 'k', 'FontSize', 8, 'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'middle', 'BackgroundColor', 'w', 'EdgeColor', 'none', 'Margin', 0.5);
            
        text(ax6, xw_vals(end), gamma_w(end), [' ' salt_name ' '], ...
            'Color', 'k', 'FontSize', 8, 'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'middle', 'BackgroundColor', 'w', 'EdgeColor', 'none', 'Margin', 0.5);
    end
    
    % Set limits (each figure independent)
    xlim(ax1, [30 100]); xlim(ax2, [30 100]); xlim(ax3, [30 100]);
    xlim(ax4, [0.8 1]); xlim(ax5, [0.8 1]); xlim(ax6, [0.8 1]);
    ylim(ax1, [-500, 500]); ylim(ax2, [-500, 500]); ylim(ax7, [-2.5, 0.75]);
    legend(ax7, 'Location', 'best');
    
    fprintf('Done.\n');
end

% --- HELPER: Add Molality Safely ---
function add_molality(map, ion, amount)
    if isKey(map, ion)
        map(ion) = map(ion) + amount;
    else
        map(ion) = amount;
    end
end

% --- PITZER CALCULATION ENGINE (PHI) ---
function phi = calculate_pitzer_phi(comp, T, params, charges, A_phi)
    species = keys(comp); m = cell2mat(values(comp));
    cations = {}; anions = {};
    for i = 1:length(species), s = species{i}; if charges(s) > 0, cations{end+1} = s; elseif charges(s) < 0, anions{end+1} = s; end; end
    I = 0; Z_total = 0; sum_m = sum(m);
    for i = 1:length(species), mi = comp(species{i}); zi = abs(charges(species{i})); I = I + 0.5 * mi * zi^2; Z_total = Z_total + mi * zi; end
    sqrt_I = sqrt(I); b = 1.2; 
    f_phi = -A_phi * (sqrt_I / (1 + b * sqrt_I));
    term_B_C = 0; term_Mix = 0;
    for c_idx = 1:length(cations)
        cat = cations{c_idx}; mc = comp(cat); zc = charges(cat);
        for a_idx = 1:length(anions)
            an = anions{a_idx}; ma = comp(an); za = charges(an);
            [b0, b1, b2, c_phi, ~, ~, ~, ~] = get_BC_params(params, cat, an, T); % Ignore derivatives here
            alpha1 = 2.0; alpha2 = 12.0; if (abs(zc) == 2 && abs(za) == 2), alpha1 = 1.4; end 
            g1 = exp(-alpha1 * sqrt_I); g2 = exp(-alpha2 * sqrt_I);
            B_phi = b0 + b1*g1 + b2*g2;
            term_B_C = term_B_C + mc * ma * (B_phi + Z_total * c_phi);
        end
    end
    % (Mixing terms omitted for brevity in this view, assume they are present as in original code)
    phi = 1 + (2 / sum_m) * (I * f_phi + term_B_C + term_Mix);
end


% --- NEW PITZER CALCULATION ENGINE (dPHI/dT) ---
function dphi_dT = calculate_pitzer_dphi_dT(comp, T, params, charges)
    % Calculates d(Phi)/dT neglecting mixing terms and dA_phi/dT.
    species = keys(comp); m = cell2mat(values(comp));
    cations = {}; anions = {};
    for i = 1:length(species), s = species{i}; if charges(s) > 0, cations{end+1} = s; elseif charges(s) < 0, anions{end+1} = s; end; end
    I = 0; Z_total = 0; sum_m = sum(m);
    for i = 1:length(species), mi = comp(species{i}); zi = abs(charges(species{i})); I = I + 0.5 * mi * zi^2; Z_total = Z_total + mi * zi; end
    sqrt_I = sqrt(I);

    % Note: Assuming dA_phi/dT = 0, so df_phi/dT = 0.

    d_term_B_C = 0;
    for c_idx = 1:length(cations)
        cat = cations{c_idx}; mc = comp(cat); zc = charges(cat);
        for a_idx = 1:length(anions)
            an = anions{a_idx}; ma = comp(an); za = charges(an);
            % Get derivatives of parameters
            [~, ~, ~, ~, db0_dT, db1_dT, db2_dT, dc_phi_dT] = get_BC_params(params, cat, an, T);
            
            alpha1 = 2.0; alpha2 = 12.0; if (abs(zc) == 2 && abs(za) == 2), alpha1 = 1.4; end 
            g1 = exp(-alpha1 * sqrt_I); g2 = exp(-alpha2 * sqrt_I);
            
            % dB_phi / dT
            dB_phi_dT = db0_dT + db1_dT*g1 + db2_dT*g2;
            
            % Summation term derivative
            d_term_B_C = d_term_B_C + mc * ma * (dB_phi_dT + Z_total * dc_phi_dT);
        end
    end
    
    % d(phi)/dT = (2 / Sum_m) * d(Sum Terms)/dT
    dphi_dT = (2 / sum_m) * (d_term_B_C);
end

% --- UPDATED LOOKUP HELPERS ---
function [b0, b1, b2, c0, db0, db1, db2, dc0] = get_BC_params(tbl, s1, s2, T)
    [b0, db0] = get_p_val(tbl, 'B0', s1, s2, T);
    [b1, db1] = get_p_val(tbl, 'B1', s1, s2, T);
    [b2, db2] = get_p_val(tbl, 'B2', s1, s2, T);
    [c0, dc0] = get_p_val(tbl, 'C0', s1, s2, T);
end

function [val, dval] = get_mix_param(tbl, type, s1, s2, s3, T)
    if nargin < 6, T = s3; s3 = ''; end
    [val, dval] = get_p_val(tbl, type, s1, s2, T, s3);
end

% --- UPDATED CORE PARAMETER RETRIEVAL ---
function [val, dval_dT] = get_p_val(tbl, param_type, s1, s2, T, s3)
    if nargin < 6, s3 = ''; end
    p_col = tbl.Parameter; s1_col = tbl.("Species 1"); s2_col = tbl.("Species 2");
    if ismember('Species 3', tbl.Properties.VariableNames), s3_col = tbl.("Species 3"); else, s3_col = repmat({''}, height(tbl), 1); end

    rows = strcmp(p_col, param_type);
    if isempty(s3)
        match = rows & ((strcmp(s1_col, s1) & strcmp(s2_col, s2)) | (strcmp(s1_col, s2) & strcmp(s2_col, s1)));
    else
        match = rows & ((strcmp(s1_col, s1) & strcmp(s2_col, s2) & strcmp(s3_col, s3)) | ...
                        (strcmp(s1_col, s2) & strcmp(s2_col, s1) & strcmp(s3_col, s3)));
    end

    idx = find(match, 1);
    if isempty(idx), val = 0; dval_dT = 0; return; end

    % a(1)=A0, a(2)=A1, etc.
    a = [tbl.("A0 (25C)")(idx), tbl.("A1 (1/T)")(idx), tbl.("A2 (lnT)")(idx), tbl.("A3 (T-Tr)")(idx), tbl.("A4 (T^2)")(idx), tbl.("A5 (1/T^2)")(idx)];
    Tr = 298.15;

    % Calculate Value P(T)
    if abs(T-Tr)<0.01
        val = a(1);
        % dP/dT at T=Tr
        % dP/dT = -A1/Tr^2 + A2/Tr + A3 + 2*A4*Tr - 2*A5/Tr^3
        dval_dT = -a(2)/(Tr^2) + a(3)/Tr + a(4) + 2*a(5)*Tr - 2*a(6)/(Tr^3);
    else
        val = a(1) + a(2)*(1/T - 1/Tr) + a(3)*log(T/Tr) + a(4)*(T - Tr) + a(5)*(T^2 - Tr^2) + a(6)*(1/T^2 - 1/Tr^2);
        % Calculate dP/dT
        % dP/dT = -A1/T^2 + A2/T + A3 + 2*A4*T - 2*A5/T^3
        dval_dT = -a(2)/(T^2) + a(3)/T + a(4) + 2*a(5)*T - 2*a(6)/(T^3);
    end
end
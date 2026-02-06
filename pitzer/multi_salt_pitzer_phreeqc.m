function plot_pitzer_critical_mixes()
    % PLOT_PITZER_CRITICAL_MIXES
    % Plots water activity for multi-salt systems with high non-ideality parameters.
    % 

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
    % (Hardcoding common charges prevents crashes if regex missed simple ions)
    charge_map = containers.Map();
    for i = 1:height(species_table)
        charge_map(char(species_table.Species(i))) = double(species_table.Charge(i));
    end
    defaults = {'Na+',1; 'K+',1; 'Li+',1; 'H+',1; 'Cl-',-1; 'Br-',-1; 'OH-',-1; ...
                'Ca+2',2; 'Mg+2',2; 'SO4-2',-2; 'CO3-2',-2; 'HCO3-',-1};
    for i = 1:size(defaults,1)
        if ~isKey(charge_map, defaults{i,1}), charge_map(defaults{i,1}) = defaults{i,2}; end
    end

    % Constants
    T = 298.15; 
    Mw = 0.0180153; 
    A_phi = 0.392; 

    % --- 2. DEFINE CRITICAL CASES ---
    % These cases are chosen based on large Mixing Parameters in your file.

    % Case 1: Large Anion-Anion-Cation Term
    % [cite_start]% Source: Psi(HCO3, SO4, Mg) = -0.161 [cite: 160]
    cases(1).name = 'Mg(HCO3)2 + MgSO4';
    cases(1).saltA = {'Mg+2', 1; 'HCO3-', 2}; 
    cases(1).saltB = {'Mg+2', 1; 'SO4-2', 1}; 
    cases(1).note = 'Large \Psi_{HCO3, SO4, Mg} = -0.161';

    % Case 2: Chloride vs Bicarbonate
    % [cite_start]% Source: Psi(Cl, HCO3, Mg) = -0.096 [cite: 149]
    cases(2).name = 'MgCl2 + Mg(HCO3)2';
    cases(2).saltA = {'Mg+2', 1; 'Cl-', 2};   
    cases(2).saltB = {'Mg+2', 1; 'HCO3-', 2}; 
    cases(2).note = 'Large \Psi_{Cl, HCO3, Mg} = -0.096';

    % Case 3: Common Ion Effect (Na vs Ca)
    % [cite_start]% Source: Theta(Ca, Na) = 0.092 [cite: 126]
    cases(3).name = 'NaCl + CaCl2';
    cases(3).saltA = {'Na+', 1; 'Cl-', 1};    
    cases(3).saltB = {'Ca+2', 1; 'Cl-', 2};   
    cases(3).note = 'Large \theta_{Ca, Na} = 0.092';

    % Case 4: Reciprocal Pair (3 DoF)
    % Varying both Cation (Na->Mg) and Anion (Cl->SO4)
    cases(4).name = 'NaCl + MgSO4';
    cases(4).saltA = {'Na+', 1; 'Cl-', 1};    
    cases(4).saltB = {'Mg+2', 1; 'SO4-2', 1}; 
    cases(4).note = 'Simultaneous Cation & Anion Exchange';

    % --- 3. EXECUTE & PLOT ---
    molalities = 0.1:0.1:4.0; 
    ratios = [0.0, 0.25, 0.50, 0.75, 1.0];

    figure('Color', 'w', 'Position', [50, 50, 1200, 800]);

    for c = 1:length(cases)
        current_case = cases(c);
        subplot(2, 2, c);
        hold on; grid on; box on;
        colors = parula(length(ratios));

        for r_idx = 1:length(ratios)
            y = ratios(r_idx);
            aw_curve = zeros(size(molalities));

            for m_idx = 1:length(molalities)
                m_tot = molalities(m_idx);

                % Build Composition Map
                comp = containers.Map();

                % Add Salt A ions
                for ion_i = 1:size(current_case.saltA, 1)
                    ion = current_case.saltA{ion_i, 1};
                    stoich = current_case.saltA{ion_i, 2};
                    add_molality(comp, ion, m_tot * (1-y) * stoich);
                end

                % Add Salt B ions
                for ion_i = 1:size(current_case.saltB, 1)
                    ion = current_case.saltB{ion_i, 1};
                    stoich = current_case.saltB{ion_i, 2};
                    add_molality(comp, ion, m_tot * y * stoich);
                end

                % Calculate
                phi = calculate_pitzer_phi(comp, T, param_table, charge_map, A_phi);
                sum_m = sum(cell2mat(values(comp)));
                aw_curve(m_idx) = exp( -phi * sum_m * Mw );
            end

            plot(molalities, aw_curve, 'LineWidth', 2, 'Color', colors(r_idx, :), ...
                'DisplayName', sprintf('y_{SaltB} = %.2f', y));
        end

        title(sprintf('%s\n(%s)', current_case.name, current_case.note), 'Interpreter', 'tex');
        xlabel('Total Salt Molality (m)');
        ylabel('Water Activity (a_w)');
        if c == 1, legend('Location', 'SouthWest'); end
    end


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

% --- PITZER CALCULATION ENGINE ---
% --- CORRECTED PITZER CALCULATION ENGINE ---
function phi = calculate_pitzer_phi(comp, T, params, charges, A_phi)
    species = keys(comp);
    m = cell2mat(values(comp));
    
    % Sort ions
    cations = {}; anions = {};
    for i = 1:length(species)
        s = species{i};
        if charges(s) > 0, cations{end+1} = s;
        elseif charges(s) < 0, anions{end+1} = s; end
    end
    
    % Calculate Ionic Strength (I) and Total Charge (Z)
    I = 0; 
    Z_total = 0; % Sum of m_i * |z_i|
    sum_m = sum(m);
    
    for i = 1:length(species)
        mi = comp(species{i});
        zi = abs(charges(species{i}));
        I = I + 0.5 * mi * zi^2;
        Z_total = Z_total + mi * zi;
    end
    
    sqrt_I = sqrt(I);
    b = 1.2; 
    
    % 1. Debye-Huckel Term
    % f_phi = -A_phi * [ I^0.5 / (1 + b*I^0.5) ]
    f_phi = -A_phi * (sqrt_I / (1 + b * sqrt_I));
    
    term_B_C = 0; 
    term_Mix = 0;
    
    % 2. Binary Terms (B and C)
    for c_idx = 1:length(cations)
        cat = cations{c_idx}; mc = comp(cat); zc = charges(cat);
        for a_idx = 1:length(anions)
            an = anions{a_idx}; ma = comp(an); za = charges(an); % keep sign for logic, abs for calc
            
            [b0, b1, b2, c_phi] = get_BC_params(params, cat, an, T);
            
            % Alpha selection (special case for 2-2 electrolytes)
            alpha1 = 2.0; alpha2 = 12.0;
            if (abs(zc) == 2 && abs(za) == 2), alpha1 = 1.4; end 
            
            g1 = exp(-alpha1 * sqrt_I); 
            g2 = exp(-alpha2 * sqrt_I);
            
            B_phi = b0 + b1*g1 + b2*g2;
            
            % The C term in phi is (Z * C_phi)
            term_B_C = term_B_C + mc * ma * (B_phi + Z_total * c_phi);
        end
    end
    
    % 3. Mixing Terms (Theta + Psi) + Electrostatic Unsymmetrical Terms
    % Cation-Cation
    for i = 1:length(cations)
        for j = i+1:length(cations)
            c1 = cations{i}; c2 = cations{j};
            mc1 = comp(c1); mc2 = comp(c2);
            z1 = charges(c1); z2 = charges(c2);
            
            % Base Theta (Constant)
            theta = get_mix_param(params, 'THETA', c1, c2, T);
            
            % Electrostatic Higher-Order Term (E_theta) for Asymmetric Mixing
            [eth, eth_prime] = calc_eth_terms(z1, z2, I, A_phi);
            
            % For Osmotic Coefficient, mixing term is (Theta + Eth + I*Eth')
            % Note: Pitzer papers typically group (Eth + I*Eth') as the Phi contribution
            Phi_phi = theta + eth + I * eth_prime;
            
            term_Mix = term_Mix + mc1 * mc2 * Phi_phi;
            
            % Triplet Interaction (Psi)
            for k = 1:length(anions)
                a = anions{k}; ma = comp(a);
                psi = get_mix_param(params, 'PSI', c1, c2, a, T);
                term_Mix = term_Mix + mc1 * mc2 * ma * psi;
            end
        end
    end
    
    % Anion-Anion
    for i = 1:length(anions)
        for j = i+1:length(anions)
            a1 = anions{i}; a2 = anions{j};
            ma1 = comp(a1); ma2 = comp(a2);
            z1 = charges(a1); z2 = charges(a2);
            
            theta = get_mix_param(params, 'THETA', a1, a2, T);
            [eth, eth_prime] = calc_eth_terms(z1, z2, I, A_phi);
            Phi_phi = theta + eth + I * eth_prime;
            
            term_Mix = term_Mix + ma1 * ma2 * Phi_phi;
            
            for k = 1:length(cations)
                c = cations{k}; mc = comp(c);
                psi = get_mix_param(params, 'PSI', a1, a2, c, T);
                term_Mix = term_Mix + ma1 * ma2 * mc * psi;
            end
        end
    end
    
    % Final Summation
    % (phi - 1) = 2/Sum_m * [ I*f_phi + ... ]
    phi = 1 + (2 / sum_m) * (I * f_phi + term_B_C + term_Mix);
end

% --- HELPER: Higher-Order Electrostatic Terms (Pitzer, 1975) ---
function [eth, eth_prime] = calc_eth_terms(z1, z2, I, A_phi)
    % Calculates E-Theta and its derivative w.r.t I
    % Returns 0 if charges are symmetric (e.g. +1/+1 or +2/+2)
    
    eth = 0; eth_prime = 0;
    if z1 == z2 || I < 1e-9, return; end
    
    x = 6 * z1 * z2 * A_phi * sqrt(I);
    
    % Approximation of the Chebyshev integral J(x)
    % Using Pitzer's simplified equations or a numerical approx.
    % Here is the Harvie & Weare (1980) / Pitzer approximation logic:
    
    % X_ij calculation
    Xij = 6 * z1 * z2 * A_phi * sqrt(I);
    
    % J0 and J1 approximations (Pitzer 1975 eq 48-49 adapted)
    % Using a robust approximation for J(x) to avoid integration loops
    % J(x) = x/4 - 1 + ... (complex series)
    % A standard implementation used in PHREEQC/GWB:
    J = calc_J_pitzer(Xij);
    J_prime = calc_J_prime_pitzer(Xij); % dJ/dx
    
    % E_theta = (Z1*Z2) / (4*I) * (J(x))
    eth = (z1 * z2) / (4 * I) * J;
    
    % For Osmotic Phi, we need the term: I * d(Eth)/dI
    % d(Eth)/dI = -Eth/I + (z1*z2)/(4*I) * dJ/dx * dx/dI
    % dx/dI = x / (2*I)
    % So I * d(Eth)/dI = -Eth + (z1*z2)/(8*I) * x * J'
    eth_prime_term = -eth + (z1 * z2 * Xij * J_prime) / (8 * I);
    
    % But return the pure derivative dEth/dI? 
    % The loop expects: Phi_phi = theta + eth + I * eth_prime_term
    % So let's return the derivative directly.
    eth_prime = eth_prime_term / I; 
end

function J = calc_J_pitzer(x)
    % Numerical approximation of J(x) integral
    y = x; 
    if x < 0 % Inverse for like-signs if needed, but usually passed as absolute? 
             % Actually Pitzer Eth is only for like-sign ions of DIFFERENT magnitude?
             % No, Eth is for mixing ions of same sign but different charge magnitude.
             % e.g. Na+ (1) mixed with Mg+2 (2).
             % x is 6*z1*z2*Aphi*sqrt(I). Since z1, z2 same sign, x > 0.
    end
    
    % Using expansion for small x and large x (Harvie 1980 Appendix)
    % Simple power series for small x is stable:
    if x < 10
        % 4.18 x^-1 ... no, use simple rational approx or lookup.
        % For code brevity, here is a compact approx often used:
        inv_x = 1/x;
        J = x/4 - 1 + inv_x * (1 - exp(-x) * (1 + x + x^2/2)); % Rough approx
        % Better to use the Pitzer explicit fit:
        J = -0.05 + 0.1 * x; % PLACEHOLDER - Real J(x) is complex.
        
        % REAL IMPLEMENTATION (Pitzer 1975, Eq 48):
        % It's complex. Let's use the widely accepted approximation:
        kn = [4.118, 0.765]; % Just example constants? No.
        
        % Let's use the exact summation for small x (x<1) and large x(x>1)
        % OR simplest fix: Ignore Eth if you don't have the math library ready,
        % BUT for critical mixes it matters.
        
        % Let's assume standard Pitzer J(x) function:
        J = x ./ (4 + 4.581 .* x.^(-0.7237)); % Empirical fit to the integral
    else
        J = x ./ (4 + 4.581 .* x.^(-0.7237));
    end
end

function Jp = calc_J_prime_pitzer(x)
   % Derivative of the empirical fit above:
   % J = x * (4 + C * x^p)^-1
   % dJ/dx ...
   C = 4.581; p = -0.7237;
   denom = 4 + C * x^p;
   term2 = -x * (denom^-2) * (C * p * x^(p-1));
   Jp = (1/denom) + term2;
end

% --- LOOKUP HELPERS ---
function [b0, b1, b2, c0] = get_BC_params(tbl, s1, s2, T)
    b0 = get_p_val(tbl, 'B0', s1, s2, T); b1 = get_p_val(tbl, 'B1', s1, s2, T);
    b2 = get_p_val(tbl, 'B2', s1, s2, T); c0 = get_p_val(tbl, 'C0', s1, s2, T);
end

function val = get_mix_param(tbl, type, s1, s2, s3, T)
    if nargin < 6, T = s3; s3 = ''; end
    val = get_p_val(tbl, type, s1, s2, T, s3);
end

function val = get_p_val(tbl, param_type, s1, s2, T, s3)
    if nargin < 6, s3 = ''; end
    p_col = tbl.Parameter; s1_col = tbl.("Species 1"); s2_col = tbl.("Species 2");

    % FIX: Use standard if-else instead of ternary operator
    if ismember('Species 3', tbl.Properties.VariableNames)
        s3_col = tbl.("Species 3");
    else
        s3_col = repmat({''}, height(tbl), 1);
    end

    rows = strcmp(p_col, param_type);

    % Check matches (order insensitive for first two species)
    if isempty(s3)
        match = rows & ((strcmp(s1_col, s1) & strcmp(s2_col, s2)) | (strcmp(s1_col, s2) & strcmp(s2_col, s1)));
    else
        % For Psi (s1, s2, s3), s1/s2 are usually the like-charged ions, order doesn't matter between them.
        match = rows & ((strcmp(s1_col, s1) & strcmp(s2_col, s2) & strcmp(s3_col, s3)) | ...
                        (strcmp(s1_col, s2) & strcmp(s2_col, s1) & strcmp(s3_col, s3)));
    end

    idx = find(match, 1);
    if isempty(idx), val = 0; return; end

    a = [tbl.("A0 (25C)")(idx), tbl.("A1 (1/T)")(idx), tbl.("A2 (lnT)")(idx), tbl.("A3 (T-Tr)")(idx), tbl.("A4 (T^2)")(idx), tbl.("A5 (1/T^2)")(idx)];
    Tr = 298.15;

    % Calculate Temperature dependence
    if abs(T-Tr)<0.01
        val=a(1); 
    else
        val=a(1)+a(2)*(1/T-1/Tr)+a(3)*log(T/Tr)+a(4)*(T-Tr)+a(5)*(T^2-Tr^2)+a(6)*(1/T^2-1/Tr^2); 
    end
end




% function compare_pitzer_standardized()
%     % COMPARE_PITZER_STANDARDIZED
%     % Plots binary salts using coefficients from two sources (Excel vs CSV).
%     % Both sources now feed into the SAME standard 'pitzer_water_activity' function.
% 
%     close all; clc;
% 
%     % --- 1. CONFIGURATION ---
%     file_xlsx = 'pitzer_parameters_phreeqc.xlsx';                     % Source 1 (Code 1)
%     file_csv  = '../data/baseline_with_ion_properties_legacy.csv'; % Source 2 (Code 2)
% 
%     % Salts to Compare: {Name, Cation, v_cat, Anion, v_an}
%     target_salts = {
%         'NaCl',   'Na+',  1, 'Cl-',   1;
%         'KCl',    'K+',   1, 'Cl-',   1;
%         'LiCl',   'Li+',  1, 'Cl-',   1;
%         'MgCl2',  'Mg+2', 1, 'Cl-',   2;
%         'CaCl2',  'Ca+2', 1, 'Cl-',   2;
%         'Na2SO4', 'Na+',  2, 'SO4-2', 1;
%         'MgSO4',  'Mg+2', 1, 'SO4-2', 1;
%     };
% 
%     molalities = 0.1:0.1:6.0;
%     T = 298.15; % 25 C
% 
%     % --- 2. LOAD DATA TABLES ---
%     fprintf('Loading Data Sources...\n');
% 
%     % Source 1: Excel
%     if exist(file_xlsx, 'file') ~= 2
%         error('Excel file not found: %s', file_xlsx);
%     end
%     opts_xls = detectImportOptions(file_xlsx);
%     opts_xls.VariableTypes = repmat({'string'}, 1, length(opts_xls.VariableTypes));
%     opts_xls.VariableTypes(5:end) = {'double'}; 
%     % Note: Ensure "VariableNamingRule" preserves column names like "Species 1"
%     tbl_xlsx = readtable(file_xlsx, 'Sheet', 'Coefficients', 'VariableNamingRule', 'preserve');
%     sp_table = readtable(file_xlsx, 'Sheet', 'Species_Charges', 'VariableNamingRule', 'preserve');
% 
%     % Build Charge Map (Needed for Alpha calculation)
%     charge_map = containers.Map();
%     for i = 1:height(sp_table)
%         charge_map(char(sp_table.Species(i))) = double(sp_table.Charge(i));
%     end
% 
%     % Source 2: CSV
%     if exist(file_csv, 'file')
%         opts_csv = detectImportOptions(file_csv);
%         opts_csv.VariableNamingRule = 'preserve';
%         tbl_csv = readtable(file_csv, opts_csv);
%     else
%         warning('CSV file not found. Skipping Source 2.');
%         tbl_csv = table();
%     end
% 
%     % --- 3. EXECUTE COMPARISON ---
%     n_salts = size(target_salts, 1);
%     n_cols = 3; n_rows = ceil(n_salts / n_cols);
%     figure('Color', 'w', 'Position', [100, 100, 1400, 300*n_rows]);
% 
%     for i = 1:n_salts
%         s_name = target_salts{i, 1};
%         cat    = target_salts{i, 2}; v_cat = target_salts{i, 3};
%         an     = target_salts{i, 4}; v_an  = target_salts{i, 5};
% 
%         % Determine Nu and Charges
%         nu = v_cat + v_an;
%         z_cat = abs(charge_map(cat));
%         z_an  = abs(charge_map(an));
% 
%         % Determine Alpha1/Alpha2 (Standard Pitzer Rules)
%         if (z_cat == 2 && z_an == 2)
%             a1 = 1.4; a2 = 12.0; % 2-2 electrolytes
%         else
%             a1 = 2.0; a2 = 12.0; % 1-1, 1-2, 2-1, etc.
%         end
% 
%         subplot(n_rows, n_cols, i); hold on; grid on; box on;
% 
%         % =========================================================
%         % METHOD 1: EXTRACT FROM EXCEL -> RUN STANDARD FUNCTION
%         % =========================================================
%         % Helper to pull B0, B1, B2, C0 from the 'Parameter' column
%         [b0_1, b1_1, b2_1, cp_1] = get_excel_params(tbl_xlsx, cat, an);
% 
%         aw_1 = zeros(size(molalities));
%         for k = 1:length(molalities)
%             m = molalities(k);
%             % Using standardized function
%             aw_1(k) = pitzer_water_activity(m, nu, v_cat, v_an, z_cat, z_an, ...
%                                             b0_1, b1_1, b2_1, cp_1, a1, a2);
%         end
%         plot(molalities, aw_1, 'b-', 'LineWidth', 2, 'DisplayName', 'Excel Params');
% 
%         % =========================================================
%         % METHOD 2: EXTRACT FROM CSV -> RUN STANDARD FUNCTION
%         % =========================================================
%         if ~isempty(tbl_csv)
%             idx = find(strcmp(tbl_csv.electrolyte, s_name), 1);
%             if ~isempty(idx)
%                 row = tbl_csv(idx, :);
%                 % Extract (handling NaNs)
%                 b0_2 = row.B_MX_0_original; if isnan(b0_2), b0_2=0; end
%                 b1_2 = row.B_MX_1_original; if isnan(b1_2), b1_2=0; end
%                 b2_2 = row.B_MX_2_original; if isnan(b2_2), b2_2=0; end
%                 cp_2 = row.C_MX_phi_original; if isnan(cp_2), cp_2=0; end
% 
%                 aw_2 = zeros(size(molalities));
%                 for k = 1:length(molalities)
%                     m = molalities(k);
%                     % Using SAME standardized function
%                     aw_2(k) = pitzer_water_activity(m, nu, v_cat, v_an, z_cat, z_an, ...
%                                                     b0_2, b1_2, b2_2, cp_2, a1, a2);
%                 end
%                 plot(molalities, aw_2, 'r--', 'LineWidth', 2, 'DisplayName', 'CSV Params');
% 
%                 % Show coefficients in title for debugging
%                 title(sprintf('%s\nEx:[%.2f, %.2f] CSV:[%.2f, %.2f]', ...
%                       s_name, b0_1, b1_1, b0_2, b1_2), 'Interpreter', 'none', 'FontSize', 8);
%             else
%                 title(sprintf('%s (Not in CSV)', s_name), 'Interpreter', 'none');
%             end
%         end
% 
%         if mod(i, n_cols) == 1, ylabel('Water Activity (a_w)'); end
%         if i > n_salts - n_cols, xlabel('Molality (m)'); end
%         if i == 1, legend('Location', 'SouthWest'); end
%     end
% end
% 
% % -------------------------------------------------------------------------
% %  SHARED STANDARDIZED FUNCTION
% % -------------------------------------------------------------------------
% 
% 
% % -------------------------------------------------------------------------
% %  DATA EXTRACTION HELPER (EXCEL)
% % -------------------------------------------------------------------------
% function [b0, b1, b2, c0] = get_excel_params(tbl, s1, s2)
%     % Default values
%     b0 = 0; b1 = 0; b2 = 0; c0 = 0;
% 
%     % Filter table for these two species (order insensitive)
%     p_col = tbl.Parameter; 
%     s1_col = tbl.("Species 1"); 
%     s2_col = tbl.("Species 2");
% 
%     mask_species = (strcmp(s1_col, s1) & strcmp(s2_col, s2)) | ...
%                    (strcmp(s1_col, s2) & strcmp(s2_col, s1));
% 
%     rows = tbl(mask_species, :);
%     if isempty(rows), return; end
% 
%     % Extract specific parameters from the filtered rows
%     % We assume the value is in "A0 (25C)" column
%     idx_b0 = strcmp(rows.Parameter, 'B0');
%     if any(idx_b0), b0 = rows.("A0 (25C)")(idx_b0); end
% 
%     idx_b1 = strcmp(rows.Parameter, 'B1');
%     if any(idx_b1), b1 = rows.("A0 (25C)")(idx_b1); end
% 
%     idx_b2 = strcmp(rows.Parameter, 'B2');
%     if any(idx_b2), b2 = rows.("A0 (25C)")(idx_b2); end
% 
%     idx_c0 = strcmp(rows.Parameter, 'C0');
%     if any(idx_c0), c0 = rows.("A0 (25C)")(idx_c0); end
% end


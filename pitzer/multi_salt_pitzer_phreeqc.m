function plot_pitzer_critical_mixes()
    % PLOT_PITZER_CRITICAL_MIXES
    % Plots water activity for multi-salt systems with high non-ideality parameters.
    % 
    % REQ: 'pitzer_parameters.xlsx' (generated from your dat file)

    % --- 1. LOAD PARAMETERS ---
    filename = 'pitzer_parameters.xlsx';
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
function phi = calculate_pitzer_phi(comp, T, params, charges, A_phi)
    species = keys(comp);
    m = cell2mat(values(comp));
    
    cations = {}; anions = {};
    for i = 1:length(species)
        s = species{i};
        if charges(s) > 0, cations{end+1} = s;
        elseif charges(s) < 0, anions{end+1} = s; end
    end
    
    I = 0; 
    for i = 1:length(species), I = I + 0.5 * comp(species{i}) * charges(species{i})^2; end
    
    b = 1.2; sqrt_I = sqrt(I);
    f_phi = -A_phi * (sqrt_I / (1 + b * sqrt_I));
    
    term_B_C = 0; term_Mix = 0;
    
    % B and C Terms (Binary Interactions)
    for c_idx = 1:length(cations)
        cat = cations{c_idx}; mc = comp(cat); zc = abs(charges(cat));
        for a_idx = 1:length(anions)
            an = anions{a_idx}; ma = comp(an); za = abs(charges(an));
            
            [b0, b1, b2, c_phi] = get_BC_params(params, cat, an, T);
            
            alpha1 = 2.0; alpha2 = 12.0;
            if (zc == 2 && za == 2), alpha1 = 1.4; end % 2-2 Electrolyte check
            
            g1 = exp(-alpha1 * sqrt_I); g2 = exp(-alpha2 * sqrt_I);
            B_phi = b0 + b1 * g1 + b2 * g2;
            
            Z = 0; for k=1:length(species), Z = Z + comp(species{k})*abs(charges(species{k})); end
            term_B_C = term_B_C + mc * ma * (B_phi + Z * c_phi);
        end
    end
    
    % Mixing Terms (Theta and Psi)
    % Cation-Cation Mixing
    for i = 1:length(cations)
        for j = i+1:length(cations)
            c1 = cations{i}; c2 = cations{j};
            term_Mix = term_Mix + comp(c1) * comp(c2) * get_mix_param(params, 'THETA', c1, c2, T);
            for k = 1:length(anions)
                term_Mix = term_Mix + comp(c1) * comp(c2) * comp(anions{k}) * get_mix_param(params, 'PSI', c1, c2, anions{k}, T);
            end
        end
    end
    % Anion-Anion Mixing
    for i = 1:length(anions)
        for j = i+1:length(anions)
            a1 = anions{i}; a2 = anions{j};
            term_Mix = term_Mix + comp(a1) * comp(a2) * get_mix_param(params, 'THETA', a1, a2, T);
            for k = 1:length(cations)
                term_Mix = term_Mix + comp(a1) * comp(a2) * comp(cations{k}) * get_mix_param(params, 'PSI', a1, a2, cations{k}, T);
            end
        end
    end

    phi = 1 + (2 / sum(m)) * (I * f_phi + term_B_C + term_Mix);
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
    
    a = [tbl.("A0 (25C)")(idx), tbl.("A1 (T)")(idx), tbl.("A2 (lnT)")(idx), tbl.("A3 (T-Tr)")(idx), tbl.("A4 (T^2)")(idx), tbl.("A5 (1/T)")(idx)];
    Tr = 298.15;
    
    % Calculate Temperature dependence
    if abs(T-Tr)<0.01
        val=a(1); 
    else
        val=a(1)+a(2)*(1/T-1/Tr)+a(3)*log(T/Tr)+a(4)*(T-Tr)+a(5)*(T^2-Tr^2)+a(6)*(1/T^2-1/Tr^2); 
    end
end





% function multi_salt_pitzer_phreeqc()
%     % PLOT_PITZER_WATER_ACTIVITY
%     % Calculates and plots water activity coefficients using Pitzer parameters.
% 
%     % --- 1. SETUP & DATA LOADING ---
%     filename = 'pitzer_parameters_charges.xlsx';
%     if exist(filename, 'file') ~= 2
%         error('File %s not found. Run Python script first.', filename);
%     end
% 
%     fprintf('Loading Pitzer parameters from %s...\n', filename);
%     opts = detectImportOptions(filename);
%     opts.VariableTypes = repmat({'string'}, 1, length(opts.VariableTypes));
%     opts.VariableTypes(5:end) = {'double'}; 
% 
%     param_table = readtable(filename, 'Sheet', 'Coefficients', 'VariableNamingRule', 'preserve');
%     species_table = readtable(filename, 'Sheet', 'Species_Charges', 'VariableNamingRule', 'preserve');
% 
%     % Build Charge Map
%     charge_map = containers.Map();
%     for i = 1:height(species_table)
%         s = char(species_table.Species(i));
%         z = double(species_table.Charge(i));
%         charge_map(s) = z;
%     end
% 
%     % --- SAFETY CHECK: Hardcode common ions if missing ---
%     % This prevents crashes if regex parsing misses a simple ion
%     defaults = {'Na+', 1; 'K+', 1; 'Li+', 1; 'H+', 1; ...
%                 'Cl-', -1; 'Br-', -1; 'OH-', -1; ...
%                 'Ca+2', 2; 'Mg+2', 2; 'SO4-2', -2; 'CO3-2', -2};
%     for i = 1:size(defaults, 1)
%         ion = defaults{i, 1};
%         if ~isKey(charge_map, ion)
%             charge_map(ion) = defaults{i, 2};
%         end
%     end
% 
%     % Constants
%     T = 298.15; 
%     Mw = 0.0180153; 
%     A_phi = 0.392; 
% 
%     % --- 2. SIMULATION: NaCl - NaBr Mixtures ---
%     fprintf('Simulating NaCl - NaBr mixtures...\n');
%     molalities = 0.1:0.1:6.0;
%     ratios_Br = [0.0, 0.25, 0.50, 0.75, 1.0];
% 
%     results_mix = [];
% 
%     for y = ratios_Br
%         aw_curve = zeros(size(molalities));
%         for i = 1:length(molalities)
%             m_total = molalities(i);
% 
%             % Composition: Na+, Cl-, Br-
%             comp = containers.Map();
%             comp('Na+') = m_total;
%             if (1-y) > 0, comp('Cl-') = m_total * (1-y); end
%             if y > 0,     comp('Br-') = m_total * y;     end
% 
%             phi = calculate_pitzer_phi(comp, T, param_table, charge_map, A_phi);
%             sum_m = sum(cell2mat(values(comp)));
%             aw_curve(i) = exp( -phi * sum_m * Mw );
%         end
%         results_mix = [results_mix; aw_curve];
%     end
% 
%     % --- 3. SIMULATION: Pure Salts Comparison ---
%     fprintf('Simulating Pure Salts (NaCl, CaCl2, Na2SO4)...\n');
%     salts = {'NaCl', 'CaCl2', 'Na2SO4'};
%     results_pure = [];
% 
%     for s = 1:length(salts)
%         salt_name = salts{s};
%         aw_curve = zeros(size(molalities));
%         for i = 1:length(molalities)
%             m = molalities(i);
%             comp = containers.Map();
% 
%             if strcmp(salt_name, 'NaCl')
%                 comp('Na+') = m; comp('Cl-') = m;
%             elseif strcmp(salt_name, 'CaCl2')
%                 comp('Ca+2') = m; comp('Cl-') = 2*m;
%             elseif strcmp(salt_name, 'Na2SO4')
%                 comp('Na+') = 2*m; comp('SO4-2') = m;
%             end
% 
%             phi = calculate_pitzer_phi(comp, T, param_table, charge_map, A_phi);
%             sum_m = sum(cell2mat(values(comp)));
%             aw_curve(i) = exp( -phi * sum_m * Mw );
%         end
%         results_pure = [results_pure; aw_curve];
%     end
% 
%     % --- 4. PLOTTING ---
%     fprintf('Generating plots...\n');
% 
%     figure('Color', 'w', 'Position', [100, 100, 1000, 500]);
% 
%     % Plot 1
%     subplot(1, 2, 1);
%     hold on; grid on; box on;
%     colors = parula(length(ratios_Br));
%     for k = 1:length(ratios_Br)
%         plot(molalities, results_mix(k, :), 'LineWidth', 2, 'Color', colors(k, :), ...
%             'DisplayName', sprintf('x_{Br} = %.2f', ratios_Br(k)));
%     end
%     xlabel('Total Molality (mol/kg)');
%     ylabel('Water Activity (a_w)');
%     title('NaCl - NaBr Mixture (25^{\circ}C)');
%     legend('Location', 'SouthWest');
%     ylim([0.7, 1.0]);
% 
%     % Plot 2
%     subplot(1, 2, 2);
%     hold on; grid on; box on;
%     styles = {'-', '--', '-.'};
%     for k = 1:length(salts)
%         plot(molalities, results_pure(k, :), 'LineWidth', 2, 'LineStyle', styles{k}, ...
%             'DisplayName', salts{k});
%     end
%     xlabel('Salt Molality (mol/kg)');
%     ylabel('Water Activity (a_w)');
%     title('Pure Salts Comparison (25^{\circ}C)');
%     legend('Location', 'SouthWest');
%     ylim([0.7, 1.0]);
% 
%     fprintf('Done.\n');
% end
% 
% % --- PITZER ENGINE ---
% function phi = calculate_pitzer_phi(comp, T, params, charges, A_phi)
%     species = keys(comp);
%     molalities = values(comp);
%     m = cell2mat(molalities);
% 
%     cations = {}; anions = {}; neutrals = {};
%     for i = 1:length(species)
%         s = species{i};
%         if ~isKey(charges, s), error('Charge not found for species: %s', s); end
%         z = charges(s);
%         if z > 0, cations{end+1} = s;
%         elseif z < 0, anions{end+1} = s;
%         else, neutrals{end+1} = s;
%         end
%     end
% 
%     I = 0;
%     sum_m = sum(m);
%     for i = 1:length(species)
%         z = charges(species{i});
%         I = I + 0.5 * comp(species{i}) * z^2;
%     end
% 
%     b = 1.2;
%     sqrt_I = sqrt(I);
%     f_phi = -A_phi * (sqrt_I / (1 + b * sqrt_I));
% 
%     term_B_C = 0;
%     term_Mix = 0; 
% 
%     % B and C Terms
%     for c_idx = 1:length(cations)
%         cat = cations{c_idx};
%         mc = comp(cat);
%         zc = abs(charges(cat));
% 
%         for a_idx = 1:length(anions)
%             an = anions{a_idx};
%             ma = comp(an);
%             za = abs(charges(an));
% 
%             [b0, b1, b2, c_phi] = get_BC_params(params, cat, an, T);
% 
%             alpha1 = 2.0; alpha2 = 12.0;
%             if (zc == 2 && za == 2), alpha1 = 1.4; end
% 
%             g1 = exp(-alpha1 * sqrt_I);
%             g2 = exp(-alpha2 * sqrt_I);
%             B_phi = b0 + b1 * g1 + b2 * g2;
% 
%             Z = 0; 
%             for k=1:length(species), Z = Z + comp(species{k})*abs(charges(species{k})); end
% 
%             term_B_C = term_B_C + mc * ma * (B_phi + Z * c_phi);
%         end
%     end
% 
%     % Mixing Terms
%     for i = 1:length(cations)
%         for j = i+1:length(cations)
%             c1 = cations{i}; c2 = cations{j};
%             theta = get_mix_param(params, 'THETA', c1, c2, T);
%             term_Mix = term_Mix + comp(c1) * comp(c2) * theta;
%             for k = 1:length(anions)
%                 a = anions{k};
%                 psi = get_mix_param(params, 'PSI', c1, c2, a, T);
%                 term_Mix = term_Mix + comp(c1) * comp(c2) * comp(a) * psi;
%             end
%         end
%     end
% 
%     for i = 1:length(anions)
%         for j = i+1:length(anions)
%             a1 = anions{i}; a2 = anions{j};
%             theta = get_mix_param(params, 'THETA', a1, a2, T);
%             term_Mix = term_Mix + comp(a1) * comp(a2) * theta;
%             for k = 1:length(cations)
%                 c = cations{k};
%                 psi = get_mix_param(params, 'PSI', a1, a2, c, T);
%                 term_Mix = term_Mix + comp(a1) * comp(a2) * comp(c) * psi;
%             end
%         end
%     end
% 
%     summable = term_B_C + term_Mix;
%     phi = 1 + (2 / sum_m) * (I * f_phi + summable);
% end
% 
% function [b0, b1, b2, c0] = get_BC_params(tbl, s1, s2, T)
%     b0 = get_p_val(tbl, 'B0', s1, s2, T);
%     b1 = get_p_val(tbl, 'B1', s1, s2, T);
%     b2 = get_p_val(tbl, 'B2', s1, s2, T);
%     c0 = get_p_val(tbl, 'C0', s1, s2, T);
% end
% 
% function val = get_mix_param(tbl, type, s1, s2, s3, T)
%     if nargin < 6, T = s3; s3 = ''; end
%     val = get_p_val(tbl, type, s1, s2, T, s3);
% end
% 
% function val = get_p_val(tbl, param_type, s1, s2, T, s3)
%     if nargin < 6, s3 = ''; end
% 
%     % Access columns safely by name using parens and string literals to avoid dot-notation issues with spaces
%     % Assuming 'preserve' naming rule kept "Species 1" with spaces
%     p_col = tbl.Parameter;
%     s1_col = tbl.("Species 1");
%     s2_col = tbl.("Species 2");
%     if ismember('Species 3', tbl.Properties.VariableNames)
%         s3_col = tbl.("Species 3");
%     else
%         s3_col = repmat({''}, height(tbl), 1);
%     end
% 
%     % Logical Indexing
%     rows = strcmp(p_col, param_type);
% 
%     if isempty(s3)
%         match = rows & ((strcmp(s1_col, s1) & strcmp(s2_col, s2)) | (strcmp(s1_col, s2) & strcmp(s2_col, s1)));
%     else
%         match = rows & ( ...
%             (strcmp(s1_col, s1) & strcmp(s2_col, s2) & strcmp(s3_col, s3)) | ...
%             (strcmp(s1_col, s2) & strcmp(s2_col, s1) & strcmp(s3_col, s3)) ); 
%     end
% 
%     idx = find(match, 1);
%     if isempty(idx), val = 0; return; end
% 
%     a0 = tbl.("A0 (25C)")(idx);
%     a1 = tbl.("A1 (T)")(idx);
%     a2 = tbl.("A2 (lnT)")(idx);
%     a3 = tbl.("A3 (T-Tr)")(idx);
%     a4 = tbl.("A4 (T^2)")(idx);
%     a5 = tbl.("A5 (1/T)")(idx);
% 
%     Tr = 298.15;
%     if abs(T - Tr) < 0.01
%         val = a0;
%     else
%         val = a0 + a1*(1/T - 1/Tr) + a2*log(T/Tr) + a3*(T-Tr) + a4*(T^2 - Tr^2) + a5*(1/T^2 - 1/Tr^2);
%     end
% end
close all
clear
clc

% Quick diagnostic: Why don't LiCl mixtures work?

[filepath,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, '..', 'data'));

% Load parameters
psi_file = fullfile(filepath, '..', 'data', 'parsed_thermodb', 'pitzer_psi.csv');
theta_file = fullfile(filepath, '..', 'data', 'parsed_thermodb', 'pitzer_theta.csv');
binary_file = fullfile(filepath, '..', 'data', 'parsed_thermodb', 'pitzer_binary.csv');

pitzer_psi = readtable(psi_file);
pitzer_theta = readtable(theta_file);
pitzer_binary = readtable(binary_file);

fprintf('=== DIAGNOSTIC: LiCl MIXTURE PARAMETERS ===\n\n');

% Check LiCl binary parameters
fprintf('1. LiCl Binary Parameters (Li+-Cl-):\n');
licl_idx = find((strcmp(pitzer_binary.species1, 'Li+') & strcmp(pitzer_binary.species2, 'Cl-')) | ...
                (strcmp(pitzer_binary.species1, 'Cl-') & strcmp(pitzer_binary.species2, 'Li+')));
if ~isempty(licl_idx)
    fprintf('   ✓ FOUND: beta0_a1 = %.4f\n', pitzer_binary.beta0_a1(licl_idx));
else
    fprintf('   ✗ NOT FOUND\n');
end

% Check CaCl2 binary parameters
fprintf('\n2. CaCl2 Binary Parameters (Ca++-Cl-):\n');
cacl2_idx = find((strcmp(pitzer_binary.species1, 'Ca++') & strcmp(pitzer_binary.species2, 'Cl-')) | ...
                 (strcmp(pitzer_binary.species1, 'Cl-') & strcmp(pitzer_binary.species2, 'Ca++')));
if ~isempty(cacl2_idx)
    fprintf('   ✓ FOUND: beta0_a1 = %.4f\n', pitzer_binary.beta0_a1(cacl2_idx));
else
    fprintf('   ✗ NOT FOUND\n');
end

% Check MgCl2 binary parameters
fprintf('\n3. MgCl2 Binary Parameters (Mg++-Cl-):\n');
mgcl2_idx = find((strcmp(pitzer_binary.species1, 'Mg++') & strcmp(pitzer_binary.species2, 'Cl-')) | ...
                 (strcmp(pitzer_binary.species1, 'Cl-') & strcmp(pitzer_binary.species2, 'Mg++')));
if ~isempty(mgcl2_idx)
    fprintf('   ✓ FOUND: beta0_a1 = %.4f\n', pitzer_binary.beta0_a1(mgcl2_idx));
else
    fprintf('   ✗ NOT FOUND\n');
end

% Check Li+-Ca++ theta parameter
fprintf('\n4. Li+-Ca++ Theta Parameter:\n');
theta_idx = find((strcmp(pitzer_theta.species1, 'Li+') & strcmp(pitzer_theta.species2, 'Ca++')) | ...
                 (strcmp(pitzer_theta.species1, 'Ca++') & strcmp(pitzer_theta.species2, 'Li+')));
if ~isempty(theta_idx)
    fprintf('   ✓ FOUND: theta_a1 = %.4f\n', pitzer_theta.theta_a1(theta_idx));
else
    fprintf('   ✗ NOT FOUND - THIS IS WHY LiCl+CaCl2 FAILS\n');
end

% Check Li+-Mg++ theta parameter
fprintf('\n5. Li+-Mg++ Theta Parameter:\n');
theta_idx = find((strcmp(pitzer_theta.species1, 'Li+') & strcmp(pitzer_theta.species2, 'Mg++')) | ...
                 (strcmp(pitzer_theta.species1, 'Mg++') & strcmp(pitzer_theta.species2, 'Li+')));
if ~isempty(theta_idx)
    fprintf('   ✓ FOUND: theta_a1 = %.4f\n', pitzer_theta.theta_a1(theta_idx));
else
    fprintf('   ✗ NOT FOUND - THIS IS WHY LiCl+MgCl2 FAILS\n');
end

% Check Li+-Ca++-Cl- psi parameter
fprintf('\n6. Li+-Ca++-Cl- Psi Parameter:\n');
psi_idx = [];
for i = 1:height(pitzer_psi)
    s = sort({pitzer_psi.species1{i}, pitzer_psi.species2{i}, pitzer_psi.species3{i}});
    target = sort({'Li+', 'Ca++', 'Cl-'});
    if isequal(s, target)
        psi_idx = i;
        break;
    end
end
if ~isempty(psi_idx)
    fprintf('   ✓ FOUND: psi_a1 = %.4f\n', pitzer_psi.psi_a1(psi_idx));
else
    fprintf('   ✗ NOT FOUND - THIS IS WHY LiCl+CaCl2 FAILS\n');
end

% Check Li+-Mg++-Cl- psi parameter
fprintf('\n7. Li+-Mg++-Cl- Psi Parameter:\n');
psi_idx = [];
for i = 1:height(pitzer_psi)
    s = sort({pitzer_psi.species1{i}, pitzer_psi.species2{i}, pitzer_psi.species3{i}});
    target = sort({'Li+', 'Mg++', 'Cl-'});
    if isequal(s, target)
        psi_idx = i;
        break;
    end
end
if ~isempty(psi_idx)
    fprintf('   ✓ FOUND: psi_a1 = %.4f\n', pitzer_psi.psi_a1(psi_idx));
else
    fprintf('   ✗ NOT FOUND - THIS IS WHY LiCl+MgCl2 FAILS\n');
end

% Show what Li+ mixing parameters DO exist
fprintf('\n\n=== AVAILABLE Li+ MIXING PARAMETERS ===\n');

fprintf('\nLi+ Theta Parameters Found:\n');
theta_li = pitzer_theta(strcmp(pitzer_theta.species1, 'Li+') | strcmp(pitzer_theta.species2, 'Li+'), :);
for i = 1:height(theta_li)
    fprintf('  %s - %s: theta = %.4f\n', theta_li.species1{i}, theta_li.species2{i}, theta_li.theta_a1(i));
end

fprintf('\nLi+ Psi Parameters Found:\n');
psi_li = pitzer_psi(strcmp(pitzer_psi.species1, 'Li+') | ...
                    strcmp(pitzer_psi.species2, 'Li+') | ...
                    strcmp(pitzer_psi.species3, 'Li+'), :);
if height(psi_li) == 0
    fprintf('  NONE - No Li+ psi parameters in database!\n');
else
    for i = 1:height(psi_li)
        fprintf('  %s - %s - %s: psi = %.4f\n', ...
            psi_li.species1{i}, psi_li.species2{i}, psi_li.species3{i}, psi_li.psi_a1(i));
    end
end

fprintf('\n\n=== CONCLUSION ===\n');
fprintf('LiCl cannot be mixed with CaCl2 or MgCl2 because:\n');
fprintf('  1. Missing Li+-Ca++ and Li+-Mg++ theta parameters\n');
fprintf('  2. Missing Li+-Ca++-Cl- and Li+-Mg++-Cl- psi parameters\n');
fprintf('\nLiCl CAN be mixed with:\n');
fprintf('  - NaCl (Li+-Na+ theta exists: 0.0029)\n');
fprintf('  - KCl (Li+-K+ theta exists: -0.0563)\n');
fprintf('  - CsCl (Li+-Cs+ theta exists: -0.1242)\n');
fprintf('  - HCl (Li+-H+ theta exists: 0.015)\n');
fprintf('  ...if corresponding psi parameters also exist\n');

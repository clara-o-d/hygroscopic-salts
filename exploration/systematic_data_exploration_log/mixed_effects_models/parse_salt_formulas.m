function [cations, anions] = parse_salt_formulas(salt_names)
    % Parse salt formulas to extract cation and anion names
    % Uses pattern matching for common ion patterns
    
    n = length(salt_names);
    cations = cell(n, 1);
    anions = cell(n, 1);
    
    % Define cation patterns (element + charge indicator)
    cation_patterns = {'Li', 'Na', 'K', 'Rb', 'Cs', 'NH4', 'Ag', ...
                       'Mg', 'Ca', 'Sr', 'Ba', 'Zn', 'Ni', 'Cu', 'Mn', 'H'};
    
    % Define anion patterns
    anion_patterns = {'Cl', 'Br', 'I', 'F', 'NO3', 'NO2', 'SO4', 'ClO4', 'ClO3', 'OH', 'BrO3'};
    
    for i = 1:n
        salt = salt_names{i};
        
        % Try to match cation
        cation_found = '';
        for j = 1:length(cation_patterns)
            % Use strncmp for Octave compatibility instead of startsWith
            if strncmp(salt, cation_patterns{j}, length(cation_patterns{j}))
                cation_found = cation_patterns{j};
                break;
            end
        end
        
        % Try to match anion (from the end or middle)
        anion_found = '';
        for j = 1:length(anion_patterns)
            % Use strfind for Octave compatibility instead of contains
            anion_pos = strfind(salt, anion_patterns{j});
            if ~isempty(anion_pos)
                anion_found = anion_patterns{j};
                break;
            end
        end
        
        % Handle special cases
        if isempty(cation_found)
            cation_found = 'Unknown';
        end
        if isempty(anion_found)
            anion_found = 'Unknown';
        end
        
        % Store with charge notation
        if strcmp(cation_found, 'NH4')
            cations{i} = 'NH4+';
        elseif strcmp(cation_found, 'H')
            cations{i} = 'H+';
        elseif ismember(cation_found, {'Li', 'Na', 'K', 'Rb', 'Cs', 'Ag'})
            cations{i} = [cation_found '+'];
        elseif ismember(cation_found, {'Mg', 'Ca', 'Sr', 'Ba', 'Zn', 'Ni', 'Cu', 'Mn'})
            cations{i} = [cation_found '2+'];
        else
            cations{i} = cation_found;
        end
        
        if ismember(anion_found, {'Cl', 'Br', 'I', 'F', 'OH'})
            anions{i} = [anion_found '-'];
        elseif ismember(anion_found, {'NO3', 'NO2', 'ClO4', 'ClO3', 'BrO3'})
            anions{i} = [anion_found '-'];
        elseif strcmp(anion_found, 'SO4')
            anions{i} = 'SO4^2-';
        else
            anions{i} = anion_found;
        end
    end
end

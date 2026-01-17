function polynomial_fitter_tabs()
    % 1. GET DATA
    solutes = load_all_data(); 
    
    % 2. PREPARE FIGURE
    f = figure('Name', 'Solute Fits', 'NumberTitle', 'off', 'Color', 'w');
    tabgp = uitabgroup(f); 
    
    colors = lines(length(solutes)); % Generate distinct colors for the summary
    
    fprintf('--- POLYNOMIAL FITTING RESULTS (4th Order) ---\n');

    % ---------------------------------------------------------
    % 3. GENERATE INDIVIDUAL TABS & CALCULATE FITS
    % ---------------------------------------------------------
    for i = 1:length(solutes)
        % Extract Data
        name = solutes(i).name;
        XY = solutes(i).data;
        x = XY(:, 1);
        y = XY(:, 2);
        
        % Perform Fit (4th Order)
        [p, S] = polyfit(x, y, 4);
        
        % Generate Smooth Fit Line (Strictly within data range - NO EXTRAPOLATION)
        x_fit = linspace(min(x), max(x), 200);
        y_fit = polyval(p, x_fit);
        
        % Store fit data for the summary tab later
        solutes(i).p = p;
        solutes(i).x_fit = x_fit;
        solutes(i).y_fit = y_fit;
        solutes(i).rsq = 1 - S.normr^2 / ((length(y)-1) * var(y));
        
        % --- PRINT STATS ---
        fprintf('\n>>> %s <<<\n', name);
        fprintf('  Input Range:  [%.4f, %.4f]\n', min(x), max(x));
        fprintf('  Output Range: [%.4f, %.4f]\n', min(y), max(y));
        fprintf('  Coeffs:       %s\n', mat2str(p, 4));
        fprintf('  R-squared:    %.5f\n', solutes(i).rsq);
        
        % --- CREATE INDIVIDUAL TAB ---
        t = uitab(tabgp, 'Title', name);
        ax = axes('Parent', t);
        hold(ax, 'on'); grid(ax, 'on');
        
        % Plot Data & Fit
        plot(ax, x, y, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Data');
        plot(ax, x_fit, y_fit, 'r-', 'LineWidth', 1.5, 'DisplayName', '4th Order Fit');
        
        legend(ax, 'Location', 'best');
        title(ax, [name ' Fit (R^2 = ' num2str(solutes(i).rsq, '%.4f') ')']);
        xlabel(ax, 'Mass Fraction (Salt)'); ylabel(ax, 'RH');
    end
    
    % ---------------------------------------------------------
    % 4. CREATE SUMMARY TAB (LAST)
    % ---------------------------------------------------------
    t_summary = uitab(tabgp, 'Title', 'All Solutes Summary');
    ax_summary = axes('Parent', t_summary);
    hold(ax_summary, 'on'); grid(ax_summary, 'on');
    
    for i = 1:length(solutes)
        c = colors(i, :);
        x = solutes(i).data(:, 1);
        y = solutes(i).data(:, 2);
        
        % Plot Data Points (Circles)
        plot(ax_summary, x, y, 'o', 'Color', c, 'MarkerFaceColor', c, ...
            'HandleVisibility', 'off'); % Hide points from legend to keep it clean
        
        % Plot Fitted Lines (Solid lines)
        plot(ax_summary, solutes(i).x_fit, solutes(i).y_fit, '-', ...
            'Color', c, 'LineWidth', 2, 'DisplayName', solutes(i).name);
    end
    
    title(ax_summary, 'Comparison of Fitted Polynomials');
    xlabel(ax_summary, 'Mass Fraction (Salt)'); ylabel(ax_summary, 'RH');
    legend(ax_summary, 'Location', 'bestoutside'); % Move legend outside if crowded
end
% ==========================================
%           DATA SECTION
% ==========================================
function solutes = load_all_data()
    % Just copy-paste this block for each new solute and increment the index (i)
    i = 1;
    solutes(i).name = 'NaCl';
    solutes(i).data = [
        0.011553555	0.9934
        0.017230794	0.99
        0.028391848	0.9835
        0.055216011	0.967
        0.080598843	0.95
        0.104653474	0.932
        0.127481497	0.913
        0.149174401	0.894
        0.169814798	0.872
        0.189477472	0.852
        0.208230288	0.83
        0.22613497	0.808
        0.243247784	0.784
        0.259620126	0.762
    ];

    i = i + 1;
    solutes(i).name = 'CsCl';
    solutes(i).data = [
        0.239366969	0.993
        0.320671385	0.99
        0.440320241	0.984
        0.611419917	0.97
        0.702399388	0.955
        0.758858582	0.94
        0.79731155	0.924
        0.825187547	0.91
        0.846322946	0.894
        0.862898916	0.879
        0.876247187	0.863
        0.887226869	0.847
        0.896417029	0.833
        0.904222197	0.817
    ];
    
    % --- PASTE NEXT DATA HERE ---
    % i = i + 1;
    % solutes(i).name = 'Solute C';
    % solutes(i).data = [ ... ];

    i = i + 1;
    solutes(i).name = 'NH4Cl';
    solutes(i).data = [
        0.03078762	0.993
        0.045481299	0.99
        0.073571417	0.984
        0.137059194	0.968
        0.192403459	0.952
        0.241076621	0.937
        0.2842163	0.92
        0.322715366	0.905
        0.357284399	0.889
        0.38849595	0.873
        0.416816504	0.855
        0.442629952	0.842
        0.4662551	0.824
        0.487958898	0.812
    ];

    i = i + 1;
    solutes(i).name = 'KCl';
    solutes(i).data = [
        0.058116552	0.9935
        0.084713208	0.9903
        0.133641229	0.9839
        0.23577341	0.9682
        0.316364899	0.9525
        0.381580325	0.9363
        0.435436871	0.9201
        0.480664441	0.904
        0.519183043	0.887
        0.552382396	0.869
        0.581293181	0.852
    ];

    i = i + 1;
    solutes(i).name = 'NaNO3';
    solutes(i).data = [
        0.004073519	0.9996
        0.019694343	0.9983
        0.038664995	0.9967
        0.074279543	0.9935
        0.107317582	0.9905
        0.138287177	0.9875
        0.167068715	0.9845
        0.193911605	0.9816
        0.219225078	0.9787
        0.242928181	0.9758
        0.26521286	0.973
        0.286631502	0.9701
    ];

    i = i + 1;
    solutes(i).name = 'AgNO3';
    solutes(i).data = [
        0.444222535	0.986
        0.615408929	0.974
        0.705972934	0.965
        0.761832573	0.955
        0.762104836	0.956
        0.7999844	0.946
        0.800240395	0.947
        0.827457925	0.941
        0.827696025	0.94
        0.848479669	0.933
        0.864847618	0.928
        0.878095495	0.922
        0.888720807	0.916
        0.888819893	0.917
        0.898010266	0.911
        0.905882465	0.907
        0.912368947	0.902
        0.917954349	0.896
        0.9231079	0.892
        0.927171229	0.888
        0.927585158	0.886
        0.93504415	0.879
        0.941054429	0.87
        0.946352652	0.862
    ];

    i = i + 1;
    solutes(i).name = 'KI';
    solutes(i).data = [
        0.001527241	1
        0.003049825	0.9999
        0.007589841	0.9998
        0.015065338	0.9997
        0.029683483	0.9993
        0.071045396	0.9983
        0.132665518	0.9966
        0.234253654	0.9934
        0.314539471	0.9901
        0.379587539	0.9868
        0.433359918	0.9836
        0.478554623	0.9803
        0.517072489	0.977
        0.550291342	0.9737
    ];


end
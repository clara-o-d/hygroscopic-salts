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
        degree = solutes(i).degree;
        x = XY(:, 1);
        y = XY(:, 2);
        
        % Perform Fit (4th Order, unless )
        if check_constraints(i)
            % --- CONSTRAINED FIT (Intercept forced to 1) ---  # NOT USED
            % 1. Create Design Matrix (x^4 down to x^1, NO column of ones)
            M = [x.^4, x.^3, x.^2, x];
            % 2. Solve for coefficients using matrix division (\)
            coeffs = M \ (y - 1);
            % 3. Construct the full polynomial vector p
            p = [coeffs', 1];
            % 4. Manually calculate S.normr for your R^2 calc later
            S.normr = norm(y - polyval(p, x));
        else
            [p, S] = polyfit(x, y, degree);
        end
        
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
function solute_fixed_intercept = check_constraints(i)
    % Just copy-paste this block for each new solute and increment the index (i)
    solutes_with_fixed_intercept =  false(1,20);  %[false,false,false,false,false,false,false,false];
    solute_fixed_intercept = solutes_with_fixed_intercept(i);
end
function solutes = load_all_data()
    standard_degree = 4;

    i = 1;
    solutes(i).name = 'NaCl';
    solutes(i).degree = standard_degree;
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
    solutes(i).degree = standard_degree;
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

    i = i + 1;
    solutes(i).name = 'NH4Cl';
    solutes(i).degree = standard_degree;
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
    solutes(i).degree = standard_degree;
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
    solutes(i).degree = standard_degree;
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
    solutes(i).degree = standard_degree;
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
    solutes(i).degree = standard_degree;
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



    i = i + 1;
    solutes(i).name = 'Na2SO4';
    solutes(i).degree = standard_degree;
    solutes(i).data = [
0.100719344	0.9957
0.183006403	0.992
0.251496838	0.9882
0.309392075	0.9848
0.358974154	0.9815
0.401913661	0.9781
0.439461664	0.9749
0.472573618	0.9718
0.528301665	0.966
0.626865463	0.95
0.691357835	0.935
0.736841933	0.918
0.770642044	0.899
    ];


    i = i + 1;
    solutes(i).name = 'K2SO4';
    solutes(i).degree = standard_degree;
    solutes(i).data = [
0.144246377	0.9958
0.252124683	0.992
0.3358491	0.9884
0.402714979	0.985
0.457348455	0.9815
0.502824907	0.9783
0.54126851	0.975
0.574193596	0.972
    ];

    i = i + 1;
    solutes(i).name = '(NH4)2SO4';
    solutes(i).degree = standard_degree;
    solutes(i).data = [
0.088359762	0.9959
0.162372343	0.9921
0.225269727	0.9886
0.279380946	0.9853
0.326426824	0.982
0.367706346	0.9788
0.404218455	0.9756
0.436743954	0.9724
0.492189721	0.9662
0.592478634	0.95
0.659687859	0.935
0.707867006	0.919
0.74409618	0.902
0.772330809	0.885
0.79495413	0.867
0.813487718	0.85
0.828948628	0.831
    ];

    i = i + 1;
    solutes(i).name = 'MgSO4';
    solutes(i).degree = standard_degree;
    solutes(i).data = [
0.138560575	0.996
0.194374538	0.9941
0.243396052	0.9924
0.286793937	0.9906
0.325483392	0.9888
0.360191246	0.987
0.391502051	0.9851
0.419891271	0.9832
0.445749594	0.9812
0.491116515	0.9768
0.529618533	0.9717
0.562704239	0.966
0.591441443	0.9598
0.616634576	0.9529
0.667839923	0.932
0.706978248	0.905
    ];


    i = i + 1;
    solutes(i).name = 'MnSO4';
    solutes(i).degree = standard_degree;
    solutes(i).data = [
0.202001706	0.9961
0.275206471	0.9945
0.336108851	0.9929
0.387569665	0.9913
0.431626528	0.9898
0.469770078	0.9882
0.503115971	0.9866
0.532515816	0.985
0.558630928	0.9832
0.602987609	0.9794
0.639242947	0.975
0.669430677	0.9703
0.694956372	0.9649
0.716822589	0.959
0.759857484	0.941
0.791537797	0.919
0.815833577	0.892
0.835057266	0.862
    ];


    i = i + 1;
    solutes(i).name = 'Li2SO4';
    solutes(i).degree = standard_degree;
    solutes(i).data = [
0.062879486	0.9956
0.118319127	0.9915
0.167565583	0.9875
0.211601723	0.9833
0.251212883	0.9793
0.287034126	0.9753
0.319584536	0.971
0.349292542	0.9667
0.401550986	0.9585
0.501614567	0.935
0.573009459	0.911
0.626512436	0.883
0.668100294	0.853
    ];

    i = i + 1;
    solutes(i).name = 'NiSO4';
    solutes(i).degree = standard_degree;
    solutes(i).data = [
0.210037636	0.9962
0.285114108	0.9946
0.347158848	0.9931
0.399293953	0.9916
0.443717964	0.9901
0.482023851	0.9887
0.515394081	0.9871
0.544724846	0.9855
0.570707752	0.9838
0.614687875	0.98
0.650494053	0.9758
0.68021129	0.971
0.705271036	0.965
0.726688655	0.9588
0.768707977	0.939
    ];


    i = i + 1;
    solutes(i).name = 'CuSO4';
    solutes(i).degree = standard_degree;
    solutes(i).data = [
0.220466691	0.9963
0.297865343	0.9947
0.361282602	0.9932
0.414193065	0.9917
0.459008085	0.9901
0.497453564	0.9885
0.530797354	0.9868
0.559991756	0.9851
0.58576594	0.9833
0.629205677	0.9795
0.664399319	0.975
    ];

    i = i + 1;
    solutes(i).name = 'ZnSO4';
    solutes(i).degree = standard_degree;
    solutes(i).data = [
0.224418821	0.9962
0.302666232	0.9946
0.366571989	0.993
0.419748019	0.9914
0.464687308	0.9898
0.503166047	0.9882
0.536483979	0.9865
0.565614112	0.9847
0.59129932	0.9828
0.634520837	0.9788
0.669475003	0.9743
0.698326811	0.9692
0.722546007	0.9634
0.743165428	0.957
0.783406608	0.938
0.812745838	0.913
    ];

    i = i + 1;
    solutes(i).name = 'LiNO3';
    solutes(i).degree = standard_degree;
    solutes(i).data = [
0.025718446	0.9967
0.050147184	0.9934
0.073380854	0.99
0.095505058	0.9866
0.116597412	0.9831
0.208844138	0.9643
0.283647181	0.945
0.34552699	0.9239
0.397566268	0.9033
0.441939476	0.8818
0.480224371	0.8567
0.513593547	0.8329
0.542936646	0.8087
0.568940847	0.8723
0.592145363	0.7591
0.612979232	0.7353
    ];


    i = i + 1;
    solutes(i).name = 'KNO3';
    solutes(i).degree = standard_degree;
    solutes(i).data = [
0.05369094412962968	0.9967
0.1019102317026735	0.9938
0.14545373560523572	0.9909
0.18497011602333824	0.9881
0.22099335347118487	0.9854
0.3619894454673626	0.9726
0.45976849663153735	0.9623
0.5315598394276059	0.9528
0.5865086838340733	0.944
0.6299197409623073	0.9376
0.66508168409534	0.9315
    ];

    i = i + 1;
    solutes(i).name = 'Ba(NO3)2';
    solutes(i).degree = standard_degree;
    solutes(i).data = [
        0.274900309	0.9958
        0.431249889	0.9923
        0.532133564	0.989
        0.602619979	0.9859
    ];

    i = i + 1;
    solutes(i).name = 'Ca(NO3)2';
    solutes(i).degree = standard_degree;
    solutes(i).data = [
0.13004135961569413	0.9955
0.23015327449593306	0.9911
0.30960195937378304	0.9867
0.37418633802400036	0.9825
0.42772107668489917	0.9781
0.5991661588103011	0.9566
0.6915673591463117	0.9296
0.7493482218959041	0.9043
0.7888958514783003	0.8741
0.817664582385093	0.8458
0.8395326592457981	0.8119
0.8567170475456024	0.7791
0.8705769176680611	0.7462
0.8819919290733174	0.7123
0.8915565535282179	0.6785
0.899686983296087	0.6464
    ];


end
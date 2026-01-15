clear all;
%relative humidity = 100 * water activity
%water activity is what's given in the paper

SL = [1,2,4,8];

wa = [.118, .292, .392, .519, .620, .629, .645, .654, .665, .676, .686, .697, .706, .715, .732]; %water activity data from paper
RH = wa .* 100; %assume ideal gas to convert to relative humidity. I'm not sure how good this assumption is but it's what we use.

WaterUptake = [0, .5, .8, 1.5, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3, 3.2, 3.4, 3.7]; %moles water per moles ammonium nitrate from paper


x = linspace(0,100,1000);

%mole ratio is given and I want to conver to gram ratio.
%the internet tells me the conversion factor going from mole h20 / mole
%nitrate = .22507
WUgrams = WaterUptake .* .22507;


x1 = linspace(0,100,1000);
basePolynomial = 1.8117e-04 .* x.^2 -0.0027  .* x; %curve found using curveFitter toolbox (you should download it if you haven't already)

hydrogel = zeros(length(SL),length(basePolynomial));

LHS = zeros(length(SL),1);
for i = 1:length(SL)
    LHS = (1+(1/SL(1,i)));
    hydrogel(i,:) = basePolynomial ./ LHS;
end


figure()
plot(x1,hydrogel(1,:))
hold on
plot(x1,hydrogel(2,:))
plot(x1,hydrogel(3,:))
plot(x1,hydrogel(4,:))
legend('SL = 1', 'SL = 2', 'SL = 4', 'SL = 8')
xlabel('Relative Humidity')
ylabel('Predicted Grams H2O / Grams NHO4')
title('Second Order Model Prediction')


%for polynomial degree three
%fitted curve is 1.6791e-05 * x^3 + -0.0013 * x^2 + 0.0587 * x - 0.5504
LHS1 = zeros(length(SL),1);
thirdOrder = 1.9173e-06 .* x .^3 - 2.4989e-05 .* x.^2 + 0.0025 .* x;
hydrogelThirdOrder = zeros(length(SL),length(basePolynomial));
for i = 1:length(SL)
    LHS = (1 + (1/SL(1,i))); %this is the formula we use to predict hydrogel uptake as a function of salt concentration in the gel
    hydrogelThirdOrder(i,:) = thirdOrder ./ LHS;
end

%hydrogelThirdOrder = thirdOrder / (1+(1/SL1));
figure()
plot(x1,hydrogelThirdOrder)
hold on
plot(x1,thirdOrder)
plot(RH,WUgrams,'rx')
xlabel('Relative Humidity')
ylabel('Third Order Prediction Grams H2O / Grams NHO4')
title('Third Order Model Prediction')
legend('SL = 1', 'SL = 2', 'SL = 4', 'SL = 8','Just Salt','Experimental Data')



%enfore that the intercept is zero
%plot the salt on this graph as well
%plot experimental measurements from the paper as well


%% --- ADDED SECTION: Figure 3 (Activity Coefficient Plots) ---
% Using the Third Order Model data calculated above

figure('Position', [100, 100, 1200, 500]); % Make figure wider for subplots
conv_factor = 0.22507; % Conversion factor used in original code

% --- Prepare Subplot 1 (Gamma vs Mole Fraction) ---
subplot(1, 2, 1);
hold on; box on;
xlabel('Mole Fraction of Water (x_w)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Activity Coefficient (\gamma_w)', 'FontSize', 12, 'FontWeight', 'bold');
title('Gamma vs Mole Fraction', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
xlim([0 1]); 
ylim([0 2]); % Adjust ylim as needed

% --- Prepare Subplot 2 (Gamma vs Relative Humidity) ---
subplot(1, 2, 2);
hold on; box on;
xlabel('Relative Humidity (%)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Activity Coefficient (\gamma_w)', 'FontSize', 12, 'FontWeight', 'bold');
title('Gamma vs Relative Humidity', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
xlim([0 100]); 
ylim([0 2]); % Adjust ylim as needed

% --- Loop through SL and Plot on Both ---
for i = 1:length(SL)
    % Get Uptake (Grams Water / Grams Dry Gel)
    uptake_g_g = hydrogelThirdOrder(i, :);
    
    % Convert to Moles Water / Moles Basis
    % Ratio = (Grams Water / Grams Basis) / ConversionFactor
    uptake_m_m = uptake_g_g ./ conv_factor;
    
    % Calculate Mole Fraction of Water (xw)
    % xw = moles_water / (moles_water + moles_basis)
    xw = uptake_m_m ./ (uptake_m_m + 1);
    
    % Calculate Activity Coefficient (gamma_w)
    % Definition: RH/100 = gamma_w * xw  => gamma_w = (RH/100) ./ xw
    aw = x1 ./ 100; % Convert RH % to Water Activity (0-1)
    gamma_w = aw ./ xw;
    
    % Plot on Subplot 1
    subplot(1, 2, 1);
    plot(xw, gamma_w, 'LineWidth', 2, 'DisplayName', ['SL = ' num2str(SL(i))]);
    
    % Plot on Subplot 2
    subplot(1, 2, 2);
    plot(x1, gamma_w, 'LineWidth', 2, 'DisplayName', ['SL = ' num2str(SL(i))]);
end

% --- Calculate "Just Salt" Data ---
% The variable 'thirdOrder' contains the pure salt uptake (Grams/Grams)
uptake_g_g_salt = thirdOrder; 
uptake_m_m_salt = uptake_g_g_salt ./ conv_factor;
xw_salt = uptake_m_m_salt ./ (uptake_m_m_salt + 1);
aw_salt = x1 ./ 100;
gamma_salt = aw_salt ./ xw_salt;

% --- Plot "Just Salt" on Both ---
subplot(1, 2, 1);
plot(xw_salt, gamma_salt, 'k--', 'LineWidth', 2, 'DisplayName', 'Just Salt');
% Add Ideal Reference
plot([0 1], [1 1], 'k:', 'LineWidth', 1.5, 'DisplayName', 'Ideal');
legend('show', 'Location', 'best');

subplot(1, 2, 2);
plot(x1, gamma_salt, 'k--', 'LineWidth', 2, 'DisplayName', 'Just Salt');
% Add Ideal Reference
plot([0 100], [1 1], 'k:', 'LineWidth', 1.5, 'DisplayName', 'Ideal');
legend('show', 'Location', 'best');
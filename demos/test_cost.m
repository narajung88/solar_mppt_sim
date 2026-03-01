% =========================================================================
% TEST SCRIPT: MPPT Cost Function Optimization
% =========================================================================
% Tests all helper functions and the main fmincon optimization.
% Uses fake but physically realistic data for a 24-panel array.
% Run section by section (Ctrl+Enter) or all at once.
% =========================================================================

clear; clc; close all;

%% -------------------------------------------------------------------------
%  SECTION 1: FAKE SYSTEM PARAMETERS
%  Realistic values for a generic 300W commercial panel
% --------------------------------------------------------------------------

% Panel electrical parameters (single-diode model, per panel)
params.Iph0             = 9.0;      % [A]     photocurrent at STC (1000 W/m², 25°C)
params.I0               = 1e-10;    % [A]     dark saturation current
params.n                = 1.2;      % [-]     diode ideality factor
params.Rs               = 0.3;      % [Ohm]   series resistance per panel
params.Rsh              = 200;      % [Ohm]   shunt resistance per panel
params.Voc_cell         = 40;       % [V]     open circuit voltage per panel
params.alpha            = 0.003;    % [-]     mismatch loss coefficient

% Array configuration
N = 24;                             % total number of panels (FIXED)
s = 4;                              % panels in series per string (fixed by voltage spec)
T = 298.15;                         % [K]  operating temperature (25°C)

% Cost parameters
params.cost_per_mppt    = 150;      % [$]  cost of one MPPT unit
params.cost_per_string  = 20;       % [$]  wiring cost per string

% Acceptable power loss threshold
epsilon = 0.1;                     % tolerate up to 10% power loss

fprintf('=== SYSTEM SETUP ===\n');
fprintf('Total panels N = %d\n', N);
fprintf('Series per string s = %d\n', s);
fprintf('Temperature T = %.1f K (%.1f C)\n', T, T-273.15);
fprintf('Max acceptable power loss: %.0f%%\n\n', epsilon*100);


%% -------------------------------------------------------------------------
%  SECTION 2: FAKE IRRADIANCE DATA
%  Simulate partial shading — some panels get less sun
% --------------------------------------------------------------------------

% Uniform irradiance (clear day, no shading)
G_uniform = 1000 * ones(1, N);

% Partially shaded irradiance (e.g. shadow across one corner of array)
G_shaded = 1000 * ones(1, N);
G_shaded(1:6) = 600;    % first 6 panels shaded to 600 W/m²
G_shaded(7:8) = 800;    % next 2 partially shaded

fprintf('=== IRRADIANCE PROFILES ===\n');
fprintf('Uniform:  min=%.0f  max=%.0f  std=%.1f W/m²\n', ...
    min(G_uniform), max(G_uniform), std(G_uniform));
fprintf('Shaded:   min=%.0f  max=%.0f  std=%.1f W/m²\n\n', ...
    min(G_shaded), max(G_shaded), std(G_shaded));


%% -------------------------------------------------------------------------
%  SECTION 3: TEST cell_current()
%  Check that the implicit diode equation solves correctly
% --------------------------------------------------------------------------

fprintf('=== TEST: cell_current() ===\n');
V_test = [0, 10, 20, 30, 35, 38, 40];
fprintf('%-10s %-10s %-10s\n', 'V (V)', 'I (A)', 'P (W)');
for i = 1:length(V_test)
    I = cell_current(V_test(i), params.Iph0, params.I0, params.n, ...
                     26e-3*(T/298.15), params.Rs, params.Rsh);
    fprintf('%-10.1f %-10.4f %-10.4f\n', V_test(i), I, V_test(i)*I);
end
fprintf('\n');

% Sanity checks
I_sc = cell_current(0, params.Iph0, params.I0, params.n, 26e-3, params.Rs, params.Rsh);
assert(I_sc > 0, 'FAIL: Short circuit current should be positive');
assert(I_sc < params.Iph0 * 1.1, 'FAIL: Isc should be close to Iph0');
fprintf('PASS: cell_current() sanity checks\n\n');


%% -------------------------------------------------------------------------
%  SECTION 4: TEST array_power_mpp()
%  Check that MPP power scales correctly with s, p, G
% --------------------------------------------------------------------------

fprintf('=== TEST: array_power_mpp() ===\n');

P_1s1p = array_power_mpp(1, 1, 1000, T, params);
P_2s1p = array_power_mpp(2, 1, 1000, T, params);
P_1s2p = array_power_mpp(1, 2, 1000, T, params);

fprintf('1 series, 1 parallel:  P = %.2f W\n', P_1s1p);
fprintf('2 series, 1 parallel:  P = %.2f W  (should be ~2x above)\n', P_2s1p);
fprintf('1 series, 2 parallel:  P = %.2f W  (should be ~2x above)\n', P_1s2p);

assert(P_2s1p > P_1s1p * 1.8, 'FAIL: Doubling series should roughly double power');
assert(P_1s2p > P_1s1p * 1.8, 'FAIL: Doubling parallel should roughly double power');
fprintf('PASS: array_power_mpp() scaling checks\n\n');


%% -------------------------------------------------------------------------
%  SECTION 5: TEST total_power() — uniform vs shaded, varying M
%  More MPPTs should help more under shading than under uniform irradiance
% --------------------------------------------------------------------------

fprintf('=== TEST: total_power() vs M ===\n');
fprintf('%-6s %-18s %-18s\n', 'M', 'P_uniform (W)', 'P_shaded (W)');

M_range = 1:N/s;
P_uniform_arr = zeros(size(M_range));
P_shaded_arr  = zeros(size(M_range));

for M = M_range
    P_uniform_arr(M) = total_power(M, s, N, G_uniform, T, params);
    P_shaded_arr(M)  = total_power(M, s, N, G_shaded,  T, params);
    fprintf('%-6d %-18.2f %-18.2f\n', M, P_uniform_arr(M), P_shaded_arr(M));
end

% Under shading, adding MPPTs should help more than under uniform irradiance
gain_uniform = P_uniform_arr(end) - P_uniform_arr(1);
gain_shaded  = P_shaded_arr(end)  - P_shaded_arr(1);
fprintf('\nPower gain from M=1 to M=%d:\n', M_range(end));
fprintf('  Uniform irradiance: +%.2f W\n', gain_uniform);
fprintf('  Shaded irradiance:  +%.2f W  (should be larger)\n', gain_shaded);
assert(gain_shaded >= gain_uniform, ...
    'FAIL: Shading should benefit more from additional MPPTs');
fprintf('PASS: total_power() shading sensitivity check\n\n');


%% -------------------------------------------------------------------------
%  SECTION 6: TEST system_cost()
%  Cost should increase monotonically with M
% --------------------------------------------------------------------------

fprintf('=== TEST: system_cost() ===\n');
fprintf('%-6s %-12s\n', 'M', 'Cost ($)');
costs = zeros(size(M_range));
for M = M_range
    costs(M) = system_cost(M, s, params);
    fprintf('%-6d %-12.2f\n', M, costs(M));
end
assert(all(diff(costs) > 0), 'FAIL: Cost should increase with M');
fprintf('PASS: system_cost() monotonically increases with M\n\n');


%% -------------------------------------------------------------------------
%  SECTION 7: MAIN OPTIMIZATION — fmincon
%  Find minimum cost M subject to power constraint
% --------------------------------------------------------------------------

fprintf('=== OPTIMIZATION: fmincon ===\n');

for scenario = 1:2
    if scenario == 1
        G_vector = G_uniform;
        label = 'Uniform irradiance';
    else
        G_vector = G_shaded;
        label = 'Partial shading';
    end

    % Step 1: find P_max (best possible power with maximum M)
    M_max = N / s;
    P_max = total_power(M_max, s, N, G_vector, T, params);
    P_threshold = (1 - epsilon) * P_max;

    fprintf('\nScenario: %s\n', label);
    fprintf('P_max (M=%d) = %.2f W\n', M_max, P_max);
    fprintf('P_threshold (%.0f%% of max) = %.2f W\n', (1-epsilon)*100, P_threshold);

    % Step 2: nonlinear power constraint — P(M) >= P_threshold
    nonlcon = @(x) deal( ...
        P_threshold - total_power(round(x(1)), s, N, G_vector, T, params), ...
        [] );

    % Step 3: cost objective
    J = @(x) system_cost(round(x(1)), s, params);

    % Step 4: run fmincon
    x0  = [2];
    lb  = [1];
    ub  = [M_max];
    options = optimoptions('fmincon', 'Display', 'off', ...
                           'StepTolerance', 1e-3);

    [x_opt, J_opt] = fmincon(J, x0, [], [], [], [], lb, ub, nonlcon, options);
    M_opt = round(x_opt(1));
    P_opt = total_power(M_opt, s, N, G_vector, T, params);

    fprintf('Optimal M = %d MPPTs\n', M_opt);
    fprintf('Power at optimal M = %.2f W (%.1f%% of max)\n', ...
            P_opt, 100*P_opt/P_max);
    fprintf('System cost = $%.2f\n', J_opt);
end


%% -------------------------------------------------------------------------
%  SECTION 8: PLOT — Power & Cost vs M
% --------------------------------------------------------------------------

figure('Name', 'MPPT Tradeoff Analysis', 'Position', [100 100 1000 420]);

% --- Left plot: Power vs M ---
subplot(1,2,1);
plot(M_range, P_uniform_arr, 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor','b');
hold on;
plot(M_range, P_shaded_arr,  'r-s', 'LineWidth', 1.5, 'MarkerFaceColor','r');
yline(P_threshold, 'k--', 'Power threshold (95% of P_{max})', ...
      'LabelHorizontalAlignment','left');
xlabel('Number of MPPTs (M)');
ylabel('Array Power (W)');
title('Power vs MPPT Count');
legend('Uniform irradiance', 'Partial shading', 'Location', 'southeast');
grid on;

% --- Right plot: Marginal power gain vs cost ---
subplot(1,2,2);
marginal_P_shaded = diff(P_shaded_arr);
marginal_cost     = diff(costs);
bar_data = [marginal_P_shaded', marginal_cost'];
b = bar(M_range(2:end), bar_data);
b(1).FaceColor = [0.2 0.6 0.9];
b(2).FaceColor = [0.9 0.3 0.3];
xlabel('Adding M-th MPPT');
ylabel('Marginal Power (W)  /  Cost ($)');
title('Marginal Gain vs Marginal Cost');
legend('ΔPower (W)', 'ΔCost ($)', 'Location', 'northeast');
grid on;

sgtitle('Solar Array MPPT Optimization — Test Results');


%% =========================================================================
%  HELPER FUNCTIONS
% ==========================================================================

function I = cell_current(V, Iph, I0, n, Vt, Rs, Rsh)
    % Clamp exponent to prevent overflow before fzero evaluates
    f = @(I) I - Iph + I0*(exp(min((V + I*Rs)/(n*Vt), 500)) - 1) ...
             + (V + I*Rs)/Rsh;
    
    % Use a bracket [0, Iph] instead of a single starting point
    % Current must physically lie between 0 and Iph
    try
        I = fzero(f, [0, Iph]);
    catch
        I = 0;   % if it still fails, panel is effectively off
    end
    I = max(I, 0);
end

function P_mpp = array_power_mpp(s, p, G, T, params)
    % Estimate MPP power for a sub-array of s series, p parallel panels
    Vt  = 26e-3 * (T / 298.15);
    Iph = params.Iph0 * (G / 1000);   % scale with irradiance

    % Sweep voltage to find MPP (reliable, works for any nonlinearity)
    V_range = linspace(0, s * params.Voc_cell * 0.99, 200);
    P_best  = 0;
    for i = 1:length(V_range)
        V_cell = V_range(i) / s;
        I_cell = cell_current(V_cell, Iph, params.I0, params.n, Vt, ...
                              params.Rs, params.Rsh);
        P = V_range(i) * (p * I_cell);
        if P > P_best, P_best = P; end
    end
    P_mpp = P_best;
end

function P_total = total_power(M, s, N, G_vector, T, params)
    % Total power from M MPPTs, each managing a group of panels
    panels_per_mppt = floor(N / M);
    p               = floor(panels_per_mppt / s);
    if p < 1, P_total = 0; return; end

    P_total = 0;
    for m = 1:M
        idx     = (m-1)*panels_per_mppt + 1 : m*panels_per_mppt;
        G_group = G_vector(idx);
        G_mean  = mean(G_group);
        sigma_G = std(G_group);

        P_group = array_power_mpp(s, p, G_mean, T, params);

        % Mismatch loss: penalize intra-group irradiance variance
        if G_mean > 0
            mismatch = params.alpha * (sigma_G / G_mean)^2 * P_group;
        else
            mismatch = 0;
        end

        P_total = P_total + P_group - mismatch;
    end
end

function C = system_cost(M, s, params)
    % Cost of MPPT hardware and wiring only (panels are fixed cost, excluded)
    C = M * params.cost_per_mppt ...
      + M * s * params.cost_per_string;
end
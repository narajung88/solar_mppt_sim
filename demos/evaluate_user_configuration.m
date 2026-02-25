function evaluate_user_configuration()
% EVALUATE_USER_CONFIGURATION  Sweeps sun angles across a user-defined
%   panel grouping configuration and reports total + per-string power.
%
%   Prerequisites:
%       1. Run build_solar_array.m once to generate solar_array_sim.slx
%       2. Ensure solar_array_sim.slx is on your MATLAB path
%       3. Simscape & Simscape Electrical toolboxes installed
%
%   Output:
%       - Plot of total delivered power vs sun angle
%       - Console report of peak power per string and overall

clc; clear; close all;

%% =========================================================================
%  CONFIGURATION — edit this section to define your panel groupings
%  Each cell entry = one string (panels connected in series to one MPPT)
%  Panel numbers correspond to indices in panelCellCounts and panelTilt
% ==========================================================================
config = { [1 2 3], [4 5], [6 7], [8], [9 10 11] };

%% =========================================================================
%  SYSTEM PARAMETERS
% ==========================================================================

% Tilt angle of each panel (degrees) — one value per panel index
panelTilt = [15 15 15 5 5 0 0 20 20 20 20];

% Sun angles to sweep (degrees from horizontal)
sunAngles = 0:5:80;

% Beam irradiance (W/m^2)
Sun.Gbeam = 1000;

% Maxeon solar cell electrical parameters
pvCell.Voc  = 0.730;   % V
pvCell.Isc  = 6.18;    % A
pvCell.Vmpp = 0.632;   % V
pvCell.Impp = 5.89;    % A

% Number of cells per panel (one per panel index)
panelCellCounts = [36 36 36 20 20 26 26 46 46 30 30];

% Load resistance (ohms) — must match R_load in model workspace
R_load = 10;  %#ok<NASGU>  assigned to base workspace below

% MPPT operating voltage window (V)
Vmin_MPPT = 15;
Vmax_MPPT = 60;

%% =========================================================================
%  SETUP
% ==========================================================================
model = 'solar_array_sim';

if ~bdIsLoaded(model)
    load_system(model);
end

nStrings  = numel(config);
nAngles   = numel(sunAngles);

% Output matrices
stringPowerMatrix = zeros(nStrings, nAngles);
totalPower        = zeros(1, nAngles);

% Push R_load to base workspace so Simscape can resolve it
assignin('base', 'R_load', R_load);

%% =========================================================================
%  MAIN SWEEP
% ==========================================================================
fprintf('Starting sweep over %d sun angles...\n', nAngles);

for a = 1:nAngles
    Sun.theta  = sunAngles(a);
    systemPower = 0;

    for s = 1:nStrings
        stringPanels = config{s};
        nPanels      = numel(stringPanels);

        %% --- MPPT window pre-check ---
        Vmpp_est = sum(panelCellCounts(stringPanels)) * pvCell.Vmpp;
        if Vmpp_est < Vmin_MPPT || Vmpp_est > Vmax_MPPT
            fprintf('  [SKIP] String %d at %d°: Vmpp_est=%.1fV outside [%.0f, %.0f]V\n', ...
                s, sunAngles(a), Vmpp_est, Vmin_MPPT, Vmax_MPPT);
            continue
        end

        %% --- Compute irradiance vector for this string ---
        G = zeros(1, nPanels);
        for k = 1:nPanels
            p           = stringPanels(k);
            theta_inc   = Sun.theta - panelTilt(p);
            G(k)        = Sun.Gbeam * max(0, cosd(theta_inc));
        end

        %% --- Build simulation input ---
        in = Simulink.SimulationInput(model);

        % String geometry
        in = in.setVariable('Nseries', nPanels);
        in = in.setVariable('Gvec',    G);
        in = in.setVariable('NCells',  panelCellCounts(stringPanels));

        % Cell electrical parameters
        in = in.setVariable('Isc',  pvCell.Isc);
        in = in.setVariable('Voc',  pvCell.Voc);
        in = in.setVariable('Vmpp', pvCell.Vmpp);
        in = in.setVariable('Impp', pvCell.Impp);

        % Load
        in = in.setVariable('R_load', R_load);

        %% --- Run simulation ---
        try
            out = sim(in);
        catch ME
            warning('String %d at %d° failed to simulate: %s', ...
                s, sunAngles(a), ME.message);
            continue
        end

        %% --- Extract steady-state power ---
        % Model must log signal named "P_out" (To Workspace or signal logging)
        try
            logs    = out.logsout;
            Pts     = logs.get('P_out').Values;
            nPts    = numel(Pts.Data);
            avgWin  = max(1, nPts - 50);           % last 50 samples or all
            Pstring = mean(Pts.Data(avgWin:end));
        catch
            % Fallback: try yout directly if logsout not configured
            warning('Could not read P_out from logsout for string %d. Check signal logging.', s);
            Pstring = 0;
        end

        stringPowerMatrix(s, a) = Pstring;
        systemPower = systemPower + Pstring;
    end % string loop

    totalPower(a) = systemPower;
    fprintf('  Sun angle %3d°: Total power = %.2f W\n', sunAngles(a), systemPower);

end % angle loop

fprintf('\nSweep complete.\n');

%% =========================================================================
%  RESULTS — Plot
% ==========================================================================
figure('Name', 'Total Delivered Power vs Sun Angle', 'NumberTitle', 'off');
plot(sunAngles, totalPower, 'b-o', 'LineWidth', 2, 'MarkerSize', 5);
xlabel('Sun Angle (degrees)');
ylabel('Total Delivered Power (W)');
title('User-Defined MPPT Configuration — Power vs Sun Angle');
grid on;

% Overlay per-string curves
hold on;
colors = lines(nStrings);
for s = 1:nStrings
    plot(sunAngles, stringPowerMatrix(s,:), '--', ...
        'Color', colors(s,:), 'LineWidth', 1.2, ...
        'DisplayName', sprintf('String %d: panels %s', s, mat2str(config{s})));
end
legend('Total', arrayfun(@(s) sprintf('String %d', s), 1:nStrings, ...
    'UniformOutput', false){:}, 'Location', 'best');
hold off;

%% =========================================================================
%  RESULTS — Console report
% ==========================================================================
[maxPower, idx] = max(totalPower);
fprintf('\n========================================\n');
fprintf('  PEAK TOTAL POWER: %.2f W at %d° sun angle\n', maxPower, sunAngles(idx));
fprintf('========================================\n');

fprintf('\n--- Per-String Peak Power ---\n');
for s = 1:nStrings
    [maxStrPower, idxStr] = max(stringPowerMatrix(s,:));
    fprintf('  String %d  (Panels %-10s):  %.2f W  at %d°\n', ...
        s, mat2str(config{s}), maxStrPower, sunAngles(idxStr));
end
fprintf('\n');

end % end evaluate_user_configuration

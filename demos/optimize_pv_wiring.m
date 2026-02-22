function evaluate_user_configuration()

clc; clear; close all;

% Get user input to define configuration
% Each cell entry = one string (panels in series)
% Input must be entered in in this format: [1 2 3]; [4 5]; [6 7]; [8]; [9 10 11]

config = { [1 2 3], [4 5], [6 7], [8], [9 10 11] };
%%userInput = input('Enter panel groupings (e.g. [1 2 3]; [4 5]; [6]): ', 's');

% Convert string into cell array
%%groups = strsplit(userInput, ';');
%%config = cell(1, length(groups));

%%for i = 1:length(groups)
%%    config{i} = str2num(groups{i}); 
%%end

%%Define system settings (sun angles)
panelTilt = [15 15 15 5 5 0 0 20 20 20 20]; % degrees, arbitrary values right now

sunAngles = 0:5:80;
Sun.Gbeam = 1000;

%%Factor in the tilt of the panels

%%Maxeon solar cell data
pvCell.Voc  = 0.730;
pvCell.Isc  = 6.18;
pvCell.Vmpp = 0.632;
pvCell.Impp = 5.89;

%%Cell counts per panel: Each panel corresponds to a certain one on the
%%aeroshell
panelCellCounts = [36 36 36 20 20 26 26 46 46 30 30];

%%PowerMR MPPT limits (from the datasheet)
Vmin_MPPT = 15;
Vmax_MPPT = 60;

%Initialize a matrix to store the power values of individual array
%groupings
stringPowerMatrix = zeros(length(config), length(sunAngles));

%%Sweeping sun angles

totalPower = zeros(size(sunAngles));

for a = 1:length(sunAngles)
    
    model = "solar_array_sim";  % your Simulink/Simscape model
load_system(model);

for a = 1:length(sunAngles)
    Sun.theta = sunAngles(a);
    systemPower = 0;

    for s = 1:length(config)
        stringPanels = config{s};

        % Optional: keep your MPPT window pre-check using estimated Vmpp
        Vmpp_est = 0;
        for p = stringPanels
            Vmpp_est = Vmpp_est + panelCellCounts(p) * pvCell.Vmpp;
        end
        if Vmpp_est < Vmin_MPPT || Vmpp_est > Vmax_MPPT
            continue
        end

        % Build irradiance vector for ONLY panels in this string
        G = zeros(1, numel(stringPanels));
        for k = 1:numel(stringPanels)
            p = stringPanels(k);
            theta_inc = Sun.theta - panelTilt(p);
            thetaFactor = max(0, cosd(theta_inc));
            G(k) = Sun.Gbeam * thetaFactor;  % W/m^2
        end

        % Create simulation input and set variables
        in = Simulink.SimulationInput(model);

        % Example model parameters you create in the model workspace:
        %   Nseries = number of panels in the string
        %   Gvec    = irradiance vector for each panel in the string
        in = in.setVariable("Nseries", numel(stringPanels));
        in = in.setVariable("Gvec", G);

        % (Optional) also pass temperature, PV params, etc.
        % in = in.setVariable("Tcell", 25);

        out = sim(in);

        % Extract delivered power (you must log it as logsout signal "P_out")
        logs = out.logsout;
        Pts = logs.get("P_out").Values;     % timeseries
        Pstring = mean(Pts.Data(end-50:end)); % steady average near end

        stringPowerMatrix(s,a) = Pstring;
        systemPower = systemPower + Pstring;
    end

    totalPower(a) = systemPower;
end
end

%%Plot our results

figure;
plot(sunAngles, totalPower, 'LineWidth',2);
xlabel('Sun Angle (degrees)');
ylabel('Total Delivered Power (W)');
title('User-Defined MPPT Configuration');
grid on;

%%Print our peak power based on sun angle

[maxPower, idx] = max(totalPower);
fprintf('\nMaximum delivered power = %.2f W at %d° sun angle\n', ...
        maxPower, sunAngles(idx));

%%Print power from each individual string
fprintf('\n--- Per-String Peak Power ---\n');

for s = 1:length(config)
    [maxStrPower, idxStr] = max(stringPowerMatrix(s,:));
    fprintf('String %d (Panels %s): %.2f W at %d°\n', ...
        s, mat2str(config{s}), maxStrPower, sunAngles(idxStr));
end

end

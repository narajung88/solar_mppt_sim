%%script to run simulation
clear;
% Set workspace variables
R_load = 10;        % ohms - adjust as needed, using a constant power load need to change
Voc    = 0.6;       % volts per cell
Isc    = 8.0;       % amps
NCells = [36, 40, 48, 56, 60, 100, 100, 100, 100, 100, 100];  % cells per panel, length >= max panel index
Gvec   = 1000 * ones(1, 11); % irradiance in W/m^2

% Run simulation
model = 'actual_solar_array_sim';  % or whatever your .slx is named
load_system(model);
simOut = sim(model);

% Results are in simOut or logged to workspace depending on your To Workspace blocks
% If you used To Workspace blocks, results appear directly as variables:
% P_out_str1, P_out_str2, etc.

% Plot results
figure;
hold on;

for s = 1:5
    varName = sprintf('P_out_str%d', s);

    if isprop(simOut, varName)
        data = simOut.(varName);   % Access inside simOut

        plot(data.time, data.signals.values);
    end
end

xlabel('Time (s)');
ylabel('Power (W)');
legend(arrayfun(@(s) sprintf('String %d', s), 1:5, 'UniformOutput', false));
title('String Power Output');
grid on;
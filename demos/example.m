%% Load Sweep Script
clear;

% Workspace parameters
Voc    = 0.6;
Isc    = 8.0;
NCells = [36, 40, 48, 56, 60, 100, 100, 100, 100, 100, 100];
Gvec   = 1000 * ones(1, 11);

model = 'single_solar_array_demo';
load_system(model);

% Sweep range
R_vals = linspace(0.1, 20, 60);   % Adjust range if needed

P_max = zeros(size(R_vals));

for k = 1:length(R_vals)
    
    R_load = R_vals(k);   % Update load
    
    simOut = sim(model, 'StopTime', '1');
    
    % Assuming you are measuring total array power as P_total
    data = simOut.P_total;  
    
    % Take final steady-state value
    P_max(k) = max(data.Data);
end

% Find best load
[bestPower, idx] = max(P_max);
bestR = R_vals(idx);

fprintf('Best R = %.3f ohms\n', bestR);
fprintf('Max Power = %.3f W\n', bestPower);

% Plot result
figure;
plot(R_vals, P_max, 'LineWidth', 2);
xlabel('Load Resistance (Ohms)');
ylabel('Maximum Power (W)');
title('Load Sweep - Maximum Power vs Resistance');
grid on;
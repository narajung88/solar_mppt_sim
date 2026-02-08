%% PV MPPT Wiring Demo
% 10 panels, arbitrary MPPT wiring, heterogeneous shading
% Author: demo framework
% ---------------------------------------------------------

clear; clc; close all;

%% -------------------------------
% 1. SYSTEM PARAMETERS
% -------------------------------

Npanels = 10;

% Reference panel parameters (can individualize later)
for i = 1:Npanels
    panels(i).Voc_ref = 40;     % Open-circuit voltage [V]
    panels(i).Isc_ref = 9;      % Short-circuit current [A]
    panels(i).Vmp_ref = 33;     % MPP voltage [V]
    panels(i).Imp_ref = 8.5;    % MPP current [A]
    panels(i).alpha_I = 0.0005; % Current temp coefficient
    panels(i).beta_V  = -0.002; % Voltage temp coefficient
end

%% -------------------------------
% 2. SHADING / TEMPERATURE
% -------------------------------

% Irradiance per panel (W/m^2)
G = [1000 900 750 600 500 1000 950 800 650 400];

% Temperature per panel (°C)
T = 25 * ones(1, Npanels);

%% -------------------------------
% 3. MPPT WIRING CONFIGURATION
% -------------------------------

% Each MPPT has parallel strings
% Each string is a list of panels in series

wiring.MPPT(1).strings = { [1 2], [3] };
wiring.MPPT(2).strings = { [4 5 6] };
wiring.MPPT(3).strings = { [7], [8 9 10] };

% MPPT voltage limits
Vmin = 20;
Vmax = 150;

%% -------------------------------
% 4. SIMULATION
% -------------------------------

Ptotal = 0;

figure('Name','MPPT P–V Curves'); hold on;
colors = lines(length(wiring.MPPT));

for k = 1:length(wiring.MPPT)
    mppt = buildMPPT(panels, wiring.MPPT(k), G, T);
    [Pk, Vk] = findMPP(mppt.V, mppt.I, Vmin, Vmax);
    Ptotal = Ptotal + Pk;

    plot(mppt.V, mppt.P, 'LineWidth', 1.5, 'Color', colors(k,:));
    plot(Vk, Pk, 'o', 'Color', colors(k,:), 'MarkerSize', 8, 'LineWidth', 2);

    fprintf('MPPT %d: Pmax = %.1f W at V = %.1f V\n', k, Pk, Vk);
end

xlabel('Voltage (V)');
ylabel('Power (W)');
title('MPPT Power–Voltage Curves');
grid on;
legend(arrayfun(@(k) sprintf('MPPT %d',k),1:length(wiring.MPPT),'UniformOutput',false));

fprintf('\nTOTAL SYSTEM POWER: %.1f W\n', Ptotal);

%% =========================================================
% LOCAL FUNCTIONS
% =========================================================

function [V, I] = panelIV(panel, G, T)
    % Generates IV curve for a single panel
    Ns = 200;

    Voc = panel.Voc_ref * (1 + panel.beta_V*(T-25));
    Isc = panel.Isc_ref * (G / 1000) * (1 + panel.alpha_I*(T-25));

    V = linspace(0, Voc, Ns);

    % Empirical IV shape (fast + stable)
    I = Isc * (1 - (V / Voc).^1.3);
    I(I < 0) = 0;
end

function [V, I] = seriesStringIV(panels, panelIdx, G, T)
    % Series connection: voltages add, current limited
    Ns = 200;
    V = zeros(1, Ns);
    I = inf(1, Ns);

    for idx = panelIdx
        [Vp, Ip] = panelIV(panels(idx), G(idx), T(idx));
        V = V + Vp;
        I = min(I, Ip);
    end
end

function [V, I] = parallelStringsIV(strings)
    % Parallel connection: currents add, voltage equal
    V = strings{1}.V;
    I = zeros(size(V));

    for s = 1:length(strings)
        I = I + interp1(strings{s}.V, strings{s}.I, V, 'linear', 0);
    end
end

function mppt = buildMPPT(panels, wiringMPPT, G, T)
    % Build aggregate IV curve for one MPPT

    strings = cell(1, length(wiringMPPT.strings));

    for s = 1:length(wiringMPPT.strings)
        idx = wiringMPPT.strings{s};
        [V, I] = seriesStringIV(panels, idx, G, T);
        strings{s}.V = V;
        strings{s}.I = I;
    end

    [mppt.V, mppt.I] = parallelStringsIV(strings);
    mppt.P = mppt.V .* mppt.I;
end

function [Pmax, Vmp] = findMPP(V, I, Vmin, Vmax)
    % Static MPPT (global max power point)
    valid = (V >= Vmin) & (V <= Vmax);
    P = V(valid) .* I(valid);

    [Pmax, idx] = max(P);
    Vvalid = V(valid);
    Vmp = Vvalid(idx);
end

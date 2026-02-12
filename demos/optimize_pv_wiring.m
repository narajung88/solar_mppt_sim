function optimize_pv_wiring()

clc; clear;

%% ===============================
% USER PARAMETERS
% ===============================

Npanels = 10;
MaxMPPTs = 3;

G = [1000 900 750 600 500 1000 950 800 650 400];
T = 25 * ones(1,Npanels);

Vmin = 20;
Vmax = 200;

%% ===============================
% PANEL MODEL
% ===============================

panel.Voc = 40;
panel.Isc = 9;
panel.alpha_I = 0.0005;
panel.beta_V = -0.002;

%% ===============================
% SEARCH
% ===============================

bestPower = 0;
bestConfig = [];

panels = 1:Npanels;

for K = 1:MaxMPPTs
    
    % Generate all ways to assign panels to K groups (ordered split only)
    splits = nchoosek(1:Npanels-1, K-1);
    
    if K == 1
        configs = {{panels}};
    else
        configs = {};
        for i = 1:size(splits,1)
            idx = [0 splits(i,:) Npanels];
            config = cell(1,K);
            for k = 1:K
                config{k} = panels(idx(k)+1:idx(k+1));
            end
            configs{end+1} = config;
        end
    end
    
    % Evaluate each config
    for c = 1:length(configs)
        
        config = configs{c};
        totalPower = 0;
        
        for k = 1:length(config)
            
            panels_k = config{k};
            strings = split_into_series_strings(panels_k);
            
            mppt = buildMPPT(panel, strings, G, T);
            Pk = computeMPP(mppt.V, mppt.I, Vmin, Vmax);
            
            totalPower = totalPower + Pk;
        end
        
        if totalPower > bestPower
            bestPower = totalPower;
            bestConfig = config;
        end
    end
end

%% ===============================
% RESULT
% ===============================

fprintf("\n===== OPTIMAL WIRING =====\n");
fprintf("Total Power: %.2f W\n\n", bestPower);

for k = 1:length(bestConfig)
    fprintf("MPPT %d panels: %s\n", k, mat2str(bestConfig{k}));
end

end

%% =============================================================
% Helper functions
%% =============================================================

function strings = split_into_series_strings(panels)
% Simple strategy: put all panels in one series string
strings = {panels};
end

function mppt = buildMPPT(panel, strings, G, T)

for s = 1:length(strings)
    idx = strings{s};
    [V,I] = seriesIV(panel, idx, G, T);
    str{s}.V = V;
    str{s}.I = I;
end

[mppt.V, mppt.I] = parallelIV(str);
mppt.P = mppt.V .* mppt.I;

end

function [V,I] = seriesIV(panel, idx, G, T)

Ns = 200;
V = linspace(0, panel.Voc*length(idx), Ns);
I = inf(1,Ns);

for i = idx
    Voc_i = panel.Voc * (1 + panel.beta_V*(T(i)-25));
    Isc_i = panel.Isc * (G(i)/1000);
    
    Vp = linspace(0, Voc_i, Ns);
    Ip = Isc_i * (1 - (Vp/Voc_i).^1.3);
    Ip(Ip<0)=0;
    
    I = min(I, Ip);
end
end

function [V,I] = parallelIV(strings)

V = strings{1}.V;
I = zeros(size(V));

for s = 1:length(strings)
    I = I + interp1(strings{s}.V, strings{s}.I, V, 'linear',0);
end
end

function Pmax = computeMPP(V,I,Vmin,Vmax)

valid = (V>=Vmin & V<=Vmax);
P = V(valid).*I(valid);
Pmax = max(P);

end

function optimize_pv_wiring()

clc; clear; close all;

%% =====================================================
% BASIC SETTINGS
%% =====================================================

Npanels = 10;

% Sun conditions
Sun.Gbeam = 1000;     % W/m^2
Sun.theta = 30;       % Sun angle (degrees)

% Cell electrical parameters
pvCell.Voc = 0.6;     % Volts per cell
pvCell.Isc = 9;       % Amps at 1000 W/m^2
pvCell.beta_V = -0.002;

Vmin = 5;
Vmax = 600;

%% =====================================================
% PANEL CELL COUNTS
%% =====================================================

panelCellCounts = [ ...
    36, ...  % 3x12
    36, ...  % 3x12
    36, ...  % 3x12
    20, ...  % 4x5
    20, ...  % 4x5
    26, ...
    26, ...
    46, ... 
    46, ...
    30];     % 3x10

%% =====================================================
% DEFINE GROUPS + ANGLES PER PANEL
%% =====================================================

panelGroups = cell(1, Npanels);

for p = 1:Npanels
    
    Ncells = panelCellCounts(p);
    groupSize = floor(Ncells/3);
    
    g1 = 1:groupSize;
    g2 = groupSize+1 : 2*groupSize;
    g3 = 2*groupSize+1 : Ncells;
    
    panelGroups{p}.cells = {g1, g2, g3};
    panelGroups{p}.tilt  = [0, 20, 45];   % degrees
    
end

%% =====================================================
% COMPUTE IRRADIANCE PER CELL
%% =====================================================

G = cell(1, Npanels);

for p = 1:Npanels
    Ncells = panelCellCounts(p);
    G{p} = computeGroupIrradiance(Ncells, panelGroups{p}, Sun);
end

%% =====================================================
% CONNECT ALL PANELS IN SERIES
%% =====================================================

[V,I] = seriesIV(pvCell, 1:Npanels, G);

Pmax = computeMPP(V, I, Vmin, Vmax);

fprintf("\nTotal Maximum Power: %.2f W\n", Pmax);

end

%% =====================================================
% IRRADIANCE FROM ANGLES
%% =====================================================

function G_cells = computeGroupIrradiance(Ncells, groups, Sun)

G_cells = zeros(1, Ncells);

for g = 1:length(groups.tilt)
    
    theta = abs(Sun.theta - groups.tilt(g));
    Geff = Sun.Gbeam * cosd(theta);
    
    if Geff < 0
        Geff = 0;
    end
    
    G_cells(groups.cells{g}) = Geff;
end

end

%% =====================================================
% SERIES CONNECTION (cells → panels → string)
%% =====================================================

function [V,I] = seriesIV(pvCell, panelIdx, G)

Ns = 300;

V_total = zeros(1, Ns);
I_total = inf(1, Ns);

for p = panelIdx
    
    G_cells = G{p};
    Ncells  = length(G_cells);
    
    V_panel = linspace(0, pvCell.Voc * Ncells, Ns);
    I_panel = inf(1, Ns);
    
    for c = 1:Ncells
        
        Voc_cell = pvCell.Voc;
        Isc_cell = pvCell.Isc * (G_cells(c)/1000);
        
        Vc = linspace(0, Voc_cell, Ns);
        Ic = Isc_cell * (1 - (Vc/Voc_cell).^1.3);
        Ic(Ic < 0) = 0;
        
        I_panel = min(I_panel, Ic);
    end
    
    V_total = V_total + V_panel;
    I_total = min(I_total, I_panel);
end

V = V_total;
I = I_total;

end

%% =====================================================
% FIND MAXIMUM POWER
%% =====================================================

function Pmax = computeMPP(V, I, Vmin, Vmax)

valid = (V >= Vmin & V <= Vmax);
P = V(valid) .* I(valid);
Pmax = max(P);

end

function evaluate_user_configuration()

clc; clear; close all;

%% =====================================================
% USER INPUT: DEFINE STRING CONFIGURATION
% Each cell entry = one string (panels in series)
%% =====================================================

config = { [1 2 3], [4 5], [6 7], [8], [9 10 11] };

%% =====================================================
% SYSTEM SETTINGS
%% =====================================================

sunAngles = 0:5:80;
Sun.Gbeam = 1000;

%% =====================================================
% MAXEON GEN III Me1 CELL DATA
%% =====================================================

pvCell.Voc  = 0.730;
pvCell.Isc  = 6.18;
pvCell.Vmpp = 0.632;
pvCell.Impp = 5.89;

panelCellCounts = [36 36 36 20 20 26 26 46 46 30 30];

%% =====================================================
% MPPT LIMITS (PowerMR 10A)
%% =====================================================

Vmin_MPPT = 15;
Vmax_MPPT = 60;

%% =====================================================
% SWEEP SUN ANGLE
%% =====================================================

totalPower = zeros(size(sunAngles));

for a = 1:length(sunAngles)
    
    Sun.theta = sunAngles(a);
    thetaFactor = max(0, cosd(Sun.theta));
    
    systemPower = 0;
    
    for s = 1:length(config)
        
        stringPanels = config{s};
        
        % ---- Compute string Vmpp ----
        Vmpp_est = 0;
        for p = stringPanels
            Vmpp_est = Vmpp_est + panelCellCounts(p) * pvCell.Vmpp;
        end
        
        % ---- Check MPPT voltage window ----
        if Vmpp_est < Vmin_MPPT || Vmpp_est > Vmax_MPPT
            continue
        end
        
        % ---- Estimate raw string power ----
        Praw = 0;
        for p = stringPanels
            Ncells = panelCellCounts(p);
            Praw = Praw + Ncells * pvCell.Vmpp * pvCell.Impp * thetaFactor;
        end

        
        % ---- Apply MPPT power caps ----
        if Vmpp_est < 25
            Pstring = min(Praw,150);
        elseif Vmpp_est < 48
            Pstring = min(Praw,250);
        else
            Pstring = min(Praw,400);
        end
        
        systemPower = systemPower + Pstring;
        
    end
    
    totalPower(a) = systemPower;
    
end

%% =====================================================
% PLOT RESULT
%% =====================================================

figure;
plot(sunAngles, totalPower, 'LineWidth',2);
xlabel('Sun Angle (degrees)');
ylabel('Total Delivered Power (W)');
title('User-Defined MPPT Configuration');
grid on;

%% =====================================================
% PRINT PEAK POWER
%% =====================================================

[maxPower, idx] = max(totalPower);
fprintf('\nMaximum delivered power = %.2f W at %d° sun angle\n', ...
        maxPower, sunAngles(idx));

end

function build_solar_array_top_level()
clc;
model = 'solar_array_sim';

if bdIsLoaded(model), close_system(model, 0); end
if exist([model '.slx'], 'file'), delete([model '.slx']); end
new_system(model); open_system(model);

set_param(model, 'Solver',   'ode23t');
set_param(model, 'StopTime', '0.1');

config   = { [1 2 3], [4 5], [6 7], [8], [9 10 11] };
nStrings = numel(config);

%% --- Infrastructure ---
add_block('nesl_utility/Solver Configuration', [model '/SolverConfig']);
set_param([model '/SolverConfig'], 'Position', posXY(60, 60));

add_block('fl_lib/Electrical/Electrical Elements/Electrical Reference', [model '/Ground']);
set_param([model '/Ground'], 'Position', posXY(900, 500));

add_block('fl_lib/Electrical/Electrical Elements/Resistor', [model '/R_Load']);
set_param([model '/R_Load'], 'R', 'R_load');
set_param([model '/R_Load'], 'Position', posXY(820, 280));

add_block('fl_lib/Thermal/Thermal Elements/Thermal Reference', [model '/ThermalRef']);
set_param([model '/ThermalRef'], 'Position', posXY(500, 600));

% Resistor: LConn1=left, RConn1=right
add_line(model, 'R_Load/RConn1',      'Ground/LConn1', 'autorouting','on');
add_line(model, 'SolverConfig/RConn1','R_Load/LConn1', 'autorouting','on');

%% --- One Solar Cell block per string ---
for s = 1:nStrings
    panelIDs = config{s};
    yBase    = 80 + (s-1)*150;
    xPV      = 220;

    pvName = sprintf('S%d_PV', s);
    pvPath = [model '/' pvName];

    add_block('ee_lib/Sources/Solar Cell', pvPath);
    set_param(pvPath, 'Position', posXY(xPV, yBase));
    set_param(pvPath, 'Isc',      sprintf('Isc * mean(Gvec(%s)) / 1000', mat2str(panelIDs)));
    set_param(pvPath, 'Voc',      'Voc');
    set_param(pvPath, 'N_series', sprintf('sum(NCells(%s))',              mat2str(panelIDs)));

    % LConn1 = thermal
    add_line(model, sprintf('%s/LConn1', pvName), 'ThermalRef/LConn1', 'autorouting','on');

    % LConn2 = electrical+, RConn1 = electrical-
    % Negative terminal -> Ground
    add_line(model, sprintf('%s/RConn1', pvName), 'Ground/LConn1', 'autorouting','on');

    xSens = xPV + 170;

    % Current sensor in series on positive terminal
    % LConn1=in, RConn1=out, RConn2=PS signal
    iSensName = sprintf('S%d_I_sens', s);
    add_block('fl_lib/Electrical/Electrical Sensors/Current Sensor', [model '/' iSensName]);
    set_param([model '/' iSensName], 'Position', posXY(xSens, yBase));

    add_line(model, sprintf('%s/LConn2', pvName),    sprintf('%s/LConn1', iSensName), 'autorouting','on');
    add_line(model, sprintf('%s/RConn1', iSensName), 'R_Load/LConn1',                'autorouting','on');

    % Voltage sensor across +bus and Ground
    % LConn1=+, RConn1=-, RConn2=PS signal
    vSensName = sprintf('S%d_V_sens', s);
    add_block('fl_lib/Electrical/Electrical Sensors/Voltage Sensor', [model '/' vSensName]);
    set_param([model '/' vSensName], 'Position', posXY(xSens + 160, yBase - 80));

    add_line(model, sprintf('%s/LConn1', vSensName), 'R_Load/LConn1', 'autorouting','on');
    add_line(model, sprintf('%s/RConn1', vSensName), 'Ground/LConn1', 'autorouting','on');

    % PS-Simulink converters
    psVName = sprintf('S%d_PS_V', s);
    psIName = sprintf('S%d_PS_I', s);
    add_block('nesl_utility/PS-Simulink Converter', [model '/' psVName]);
    add_block('nesl_utility/PS-Simulink Converter', [model '/' psIName]);
    set_param([model '/' psVName], 'Position', posXY(xSens + 340, yBase - 80));
    set_param([model '/' psIName], 'Position', posXY(xSens + 340, yBase));

    add_line(model, sprintf('%s/RConn2', vSensName), sprintf('%s/RConn1', psVName), 'autorouting','on');
    add_line(model, sprintf('%s/RConn2', iSensName), sprintf('%s/RConn1', psIName), 'autorouting','on');

    % Product V*I
    pCalcName = sprintf('S%d_P_calc', s);
    add_block('simulink/Math Operations/Product', [model '/' pCalcName]);
    set_param([model '/' pCalcName], 'Position', posXY(xSens + 500, yBase - 40));

    add_line(model, [psVName '/1'], [pCalcName '/1'], 'autorouting','on');
    add_line(model, [psIName '/1'], [pCalcName '/2'], 'autorouting','on');

    % To Workspace
    pLogName = sprintf('S%d_P_log', s);
    add_block('simulink/Sinks/To Workspace', [model '/' pLogName]);
    set_param([model '/' pLogName], ...
        'VariableName', sprintf('P_out_str%d', s), ...
        'SaveFormat',   'Structure With Time');
    set_param([model '/' pLogName], 'Position', posXY(xSens + 650, yBase - 40));

    add_line(model, [pCalcName '/1'], [pLogName '/1'], 'autorouting','on');
end

save_system(model);
fprintf('\nBuilt and saved "%s.slx"\n', model);
fprintf('Workspace vars needed: R_load, Voc, Isc, NCells, Gvec\n');
end

function pos = posXY(x, y)
    w = 60; h = 40;
    pos = [x, y, x+w, y+h];
end
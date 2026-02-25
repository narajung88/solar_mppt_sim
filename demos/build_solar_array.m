function build_solar_array()
% BUILD_SOLAR_ARRAY  Programmatically constructs the solar_array_sim.slx
%   Simscape Electrical model. Run this ONCE to generate the model, then
%   use evaluate_user_configuration.m to run sweeps against it.
%
%   Model structure:
%       - 5 String subsystems, each with PV panels in series + MPPT placeholder
%       - All strings connected in parallel to a common resistive load
%       - P_out signal logged per string and as total
%
%   Requirements:
%       - Simscape & Simscape Electrical toolboxes
%       - MATLAB R2021b or later recommended

clc;

model = 'solar_array_sim';

%% --- Create or overwrite model ---
if bdIsLoaded(model)
    close_system(model, 0);
end
if exist([model '.slx'], 'file')
    delete([model '.slx']);
end
new_system(model);
open_system(model);

% Set solver for Simscape (variable-step, ode23t works well for electrical)
set_param(model, 'Solver', 'ode23t');
set_param(model, 'StopTime', '0.1');  % short sim; adjust as needed

%% --- Top-level blocks ---

% Solver Configuration (required for every Simscape model)
add_block('nesl_utility/Solver Configuration', [model '/SolverConfig']);
set_param([model '/SolverConfig'], 'Position', posXY(50, 50));

% Electrical reference (ground)
add_block('fl_lib/Electrical/Electrical Elements/Electrical Reference', ...
    [model '/Ground']);
set_param([model '/Ground'], 'Position', posXY(150, 400));

% Load resistor (shared across all parallel strings)
add_block('fl_lib/Electrical/Electrical Elements/Resistor', ...
    [model '/R_Load']);
set_param([model '/R_Load'], 'R', 'R_load');  % set R_load in workspace
set_param([model '/R_Load'], 'Position', posXY(600, 200));

%% --- String subsystem config ---
% Must match config in evaluate_user_configuration.m
config = { [1 2 3], [4 5], [6 7], [8], [9 10 11] };
nStrings = numel(config);

stringHandles = zeros(1, nStrings);

for s = 1:nStrings
    nPanels = numel(config{s});
    yPos = 100 + (s-1) * 120;
    
    subsysName = sprintf('String_%d', s);
    subsysPath = [model '/' subsysName];
    
    % Add subsystem block
    add_block('built-in/Subsystem', subsysPath);
    set_param(subsysPath, 'Position', posXY(300, yPos));
    
    % Build internals of this string subsystem
    build_string_subsystem(subsysPath, nPanels, s);
    
    stringHandles(s) = get_param(subsysPath, 'Handle');
end

%% --- Connect strings in parallel to load ---
% All string positive outputs connect to load top terminal
% All string negative outputs connect to ground
for s = 1:nStrings
    subsysName = sprintf('String_%d', s);
    
    % These line connections assume your subsystem has ports named
    % 'p' (positive out) and 'n' (negative out)
    % Adjust port names if yours differ
    add_line(model, [subsysName '/p'], 'R_Load/p', 'autorouting', 'on');
    add_line(model, [subsysName '/n'], 'Ground/LConn1', 'autorouting', 'on');
end

% Connect load bottom to ground
add_line(model, 'R_Load/n', 'Ground/LConn1', 'autorouting', 'on');

%% --- Save model ---
save_system(model);
fprintf('Model "%s.slx" built and saved successfully.\n', model);
fprintf('Open it in Simulink to verify connections before running evaluate_user_configuration.\n');

end % end build_solar_array


%% =========================================================================
function build_string_subsystem(subsysPath, nPanels, stringIndex)
% Builds the internals of one string subsystem:
%   nPanels Solar Cell blocks in series -> PS math for P=V*I -> logging

    % Delete default content in new subsystem
    delete_block(find_system(subsysPath, 'SearchDepth', 1, 'type', 'block'));

    %% Add Solar Cell blocks
    for k = 1:nPanels
        cellName  = sprintf('PV_%d', k);
        cellPath  = [subsysPath '/' cellName];

        % NOTE: Verify this library path in YOUR MATLAB version by dragging
        % a Solar Cell block into a blank model and running:
        %   get_param('tempModel/Solar Cell', 'ReferenceBlock')
        add_block('ee_lib/Sources/Solar Cell', cellPath);
        set_param(cellPath, 'Position', posXY(100 + (k-1)*150, 150));

        % Parameterize — these variable names are set by evaluate_user_configuration
        set_param(cellPath, 'Isc',    'Isc');
        set_param(cellPath, 'Voc',    'Voc');
        set_param(cellPath, 'Vmpp',   'Vmpp');
        set_param(cellPath, 'Impp',   'Impp');
        set_param(cellPath, 'Ncells', sprintf('NCells(%d)', k));
        set_param(cellPath, 'Irr',    sprintf('Gvec(%d)', k));
    end

    %% Wire Solar Cells in series: PV_1.p -> PV_2.n -> PV_3.n ...
    for k = 1:nPanels-1
        add_line(subsysPath, ...
            sprintf('PV_%d/p', k), ...
            sprintf('PV_%d/n', k+1), ...
            'autorouting', 'on');
    end

    %% Voltage sensor across the whole string (for P = V*I logging)
    add_block('fl_lib/Electrical/Electrical Sensors/Voltage Sensor', ...
        [subsysPath '/V_sens']);
    set_param([subsysPath '/V_sens'], 'Position', posXY(100 + nPanels*150, 100));

    add_block('fl_lib/Electrical/Electrical Sensors/Current Sensor', ...
        [subsysPath '/I_sens']);
    set_param([subsysPath '/I_sens'], 'Position', posXY(100 + nPanels*150, 200));

    %% PS Simulink converters (to bring physical signals into Simulink math)
    add_block('nesl_utility/PS-Simulink Converter', [subsysPath '/PS_V']);
    add_block('nesl_utility/PS-Simulink Converter', [subsysPath '/PS_I']);
    set_param([subsysPath '/PS_V'], 'Position', posXY(200 + nPanels*150, 100));
    set_param([subsysPath '/PS_I'], 'Position', posXY(200 + nPanels*150, 200));

    %% Multiply V*I to get power
    add_block('simulink/Math Operations/Product', [subsysPath '/P_calc']);
    set_param([subsysPath '/P_calc'], 'Position', posXY(300 + nPanels*150, 150));

    %% To Workspace block for logging P_out
    add_block('simulink/Sinks/To Workspace', [subsysPath '/P_log']);
    set_param([subsysPath '/P_log'], ...
        'VariableName', sprintf('P_out_str%d', stringIndex), ...
        'SaveFormat',   'Structure With Time');
    set_param([subsysPath '/P_log'], 'Position', posXY(400 + nPanels*150, 150));

    %% Wire sensors -> converters -> product -> log
    add_line(subsysPath, 'V_sens/V', 'PS_V/I',   'autorouting', 'on');
    add_line(subsysPath, 'I_sens/I', 'PS_I/I',   'autorouting', 'on');
    add_line(subsysPath, 'PS_V/O',   'P_calc/1', 'autorouting', 'on');
    add_line(subsysPath, 'PS_I/O',   'P_calc/2', 'autorouting', 'on');
    add_line(subsysPath, 'P_calc/1', 'P_log/1',  'autorouting', 'on');

    %% Add subsystem output ports (p = positive rail, n = negative rail)
    add_block('built-in/Out1', [subsysPath '/p']);
    add_block('built-in/Out1', [subsysPath '/n']);

    % Connect first panel negative and last panel positive to output ports
    add_line(subsysPath, 'PV_1/n',              'n/1', 'autorouting', 'on');
    add_line(subsysPath, sprintf('PV_%d/p', nPanels), 'p/1', 'autorouting', 'on');

end % end build_string_subsystem


%% =========================================================================
function pos = posXY(x, y)
% Returns a [left top right bottom] position vector for set_param
    w = 60; h = 40;
    pos = [x, y, x+w, y+h];
end

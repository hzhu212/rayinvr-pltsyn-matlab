function gui_pltsyn(working_dir)
% rayinvr-pltsyn 的图形界面，实现交互式绘图。
% 本脚本依赖以下 Packages，请事先安装到 MATLAB 中：
%     1. GUI Layout Toolbox: https://www.mathworks.com/matlabcentral/fileexchange/47982-gui-layout-toolbox
%     2. Widgets Toolbox: https://www.mathworks.com/matlabcentral/fileexchange/66235-widgets-toolbox
%
% Usage:
% gui_pltsyn(working_dir)
% @param: working_dir: the directory that contains "plotdata.pltsyn.mat".(i.e. the *.in folder)

    working_dir = 'D:\Archive\Research\rayinvr\rayinvr-data\examples\e3';
    % if nargin < 1
    %     fprintf('Argument "working_dir" is required.\n');
    %     return;
    % end

    data_name = 'plotdata.pltsyn.mat';
    data_path = fullfile(working_dir, data_name);
    if ~exist(data_path, 'file')
        fprintf([
            'The working_dir "%s" does not contain the plot-data file: "%s". '...
            'You should run "start_pltsyn" first to generate plot-data.\n'], ...
            working_dir, data_name);
    end

    % obj properties: labels, xtraces, xinc, data, xlabel, ylabel, xlim, ylim, vred
    obj = load_plotdata(data_path);

    settings.colors = {'r', 'g', 'b', 'c', 'm', 'y', [1,0.65,0], [0.5,0.2,0.9], [0.6,0.8,0.2], [0.4,0.2,0.2], [0.4,0.4,1]};
    settings.enable_color = false;
    settings.xlim = obj.xlim;       % x data limit
    settings.ylim = obj.ylim;       % y data limit
    settings.xtick = sort(obj.xtraces{1});
    vred = obj.vred;
    if isempty(vred) || vred == 0, vred = inf; end
    settings.vred_orig = vred;      % original v-reduce
    settings.vred = vred;           % v-reduce
    settings.xinc = obj.xinc;       % trace interval in x
    settings.ntpi = 5;              % traces per inch
    settings.mspi = 200;            % ms per inch
    settings.amp_lim = 0;           % set upper limit for amplitude when plotting

    % fill waveform:
    %   1. NO fill;
    %   2. fill with Black-White;
    %   3. fill with Red-Blue.
    %   4. custom: set peak_fill and trough_fill.
    settings.fill_wave_options = {'NO', 'Black-White', 'Red-Blue', 'Custom'};
    settings.fill_wave = 1;
    settings.peak_fill = [1, 1, 1];
    settings.trough_fill = [1, 1, 1];

    % gui is a struct with many properties attached to it
    gui = struct();

    % graphic handles of main window
    [xshots, raygroups] = get_raygroups_of_each_shot(obj.labels);
    gui.h = create_interface(xshots, raygroups, 'shot-', 'raygroup-');

    % handle window resize
    set(gui.h.window, 'ResizeFcn', @on_window_resize);

    % graphic handles of settings window
    gui.h.st = struct();
    % gui properties
    gui.pp = struct();

    % expand all the tree branches by default
    on_expand_all();

    % init view limit
    fun_calc_view_limit();

    % init axes
    % set_axes();


%% -----------------------------------------------------------------------------
function [obj] = load_plotdata(data_path)
% load_plotdata: load plotdata of pltsyn from .mat file generated by
% `pltsyn-matlab/fun_store_pltsyn_plt.m`
    load(data_path, 'obj');
end


%% -----------------------------------------------------------------------------
function [xshots, raygroups] = get_raygroups_of_each_shot(labels)
% get_raygroups_of_each_shot: split each label to xshot and ray-group code,
% and then group by xshot. e.g. labels={'2.00-3.2', '2.00-4.2', '3.00-3.2'},
% splited xshots={'2.00', '2.00', '3.00'}, splited ray-groups={'3.2', '4.2', '3.2'}.
% Then group by xshot: xshots={'2.00', '3.00'}, raygroups={{'3.2', '4.2'}, {'3.2'}}.

    labels = sort(labels);
    xshots_split = {};
    raygroups_split = {};
    for ii = 1:length(labels)
        tmp = strsplit(labels{ii}, '-');
        xshots_split{end+1} = tmp{1};
        raygroups_split{end+1} = tmp{2};
    end
    g = findgroups(xshots_split);
    xshots = unique(xshots_split);
    raygroups = splitapply(@(c) {c}, raygroups_split, g);
end


%% -----------------------------------------------------------------------------
function [h] = create_interface(tree_branches, tree_nodes, branch_prefix, node_prefix)
% initialize interface
% tree_branches, tree_nodes, branch_prefix, node_prefix are used to create CheckboxTree widget.

    h = struct();
    h.window = figure(...
        'Name', 'GUI of rayinvr-pltsyn', 'NumberTitle', 'off', 'MenuBar', 'none', ...
        'ToolBar', 'figure');

    % + File menu
    fileMenu = uimenu(h.window, 'Label', 'File');
    uimenu(fileMenu, 'Label', 'Export Plot', 'Callback', @on_export_plot);
    % uimenu(fileMenu, 'Label', 'Change Working Directory', 'Callback', @on_change_working_dir);

    % + Help menu
    helpMenu = uimenu(h.window, 'Label', 'Help');
    uimenu(helpMenu, 'Label', 'README', 'Callback', @on_readme);

    % Arrange the main interface
    mainLayout = uix.HBoxFlex('Parent', h.window, 'Spacing', 3);

    % + Create the panels
    controlPanel = uix.BoxPanel('Parent', mainLayout, 'Title', 'Select Ray Groups', 'HelpFcn', @on_control_panel_help);
    viewPanel = uix.BoxPanel('Parent', mainLayout, 'Title', 'Seismic View', 'Padding', 0);

    % + Adjust the main layout
    set(mainLayout, 'Widths', [200, -1]);

    % + Create the view
    h.viewAxes = axes('Parent', uicontainer('Parent', viewPanel));
    % Minimize the margin of axes
    set(h.viewAxes, 'LooseInset', get(h.viewAxes, 'TightInset'));
    % Add context menu to plot area
    ctxmenu1 = uicontextmenu(h.window);
    h.viewAxes.UIContextMenu = ctxmenu1;
    uimenu(ctxmenu1, 'Label', 'Enable Color', 'Callback', @on_enable_color);
    uimenu(ctxmenu1, 'Label', 'Disable Color', 'Callback', @on_disable_color);

    % + Create the controls
    controlLayout = uix.VBox('Parent', controlPanel, 'Padding', 3, 'Spacing', 3);
    h.checkboxTree = uiw.widget.CheckboxTree(...
        'Parent', controlLayout, ...
        'MouseClickedCallback', @on_mouse_clicked, ...
        'SelectionChangeFcn', @on_select_changed ...
        );
    h.checkboxTree.Root.Name = 'Select All';
    for ii = 1:length(tree_branches)
        branch_id = tree_branches{ii};
        branch = uiw.widget.CheckboxTreeNode('Name', strcat(branch_prefix, branch_id), 'Parent', h.checkboxTree.Root);
        branch.UserData = branch_id;
        branch_nodes = tree_nodes{ii};
        for jj = 1:length(branch_nodes)
            node_id = branch_nodes{jj};
            node = uiw.widget.CheckboxTreeNode('Name', strcat(node_prefix, node_id), 'Parent', branch);
            node.TooltipString = sprintf('Tick raygroups to plot.\n Select a raygroup to scale plots according to.');
            node_path = strcat(branch_id, '-', node_id);
            node.UserData = node_path;
        end
    end

    % context menu for checkboxTree
    ctxmenu2 = uicontextmenu(h.window);
    h.checkboxTree.UIContextMenu = ctxmenu2;
    uimenu(ctxmenu2, 'Label', 'Expand All', 'Callback', @on_expand_all);
    uimenu(ctxmenu2, 'Label', 'Collapse All', 'Callback', @on_collapse_all);

    buttonBox = uix.VBox('Parent', controlLayout, 'Spacing', 3);
    % Settings button
    settingButton = uicontrol(...
        'Style', 'PushButton', 'Parent', buttonBox, 'String', 'Settings', ...
        'Callback', @on_settings);
    % Apply button
    applyButton = uicontrol(...
        'Style', 'PushButton', 'Parent', buttonBox, 'String', 'Apply', ...
        'Callback', @on_apply);

    set(controlLayout, 'Heights', [-1, 65]); % Make the list fill the space
end

function [] = on_window_resize(~, ~)
% when window resize, recalculate view limit and redraw
    fun_calc_view_limit();
    redraw();
end

%% -----------------------------------------------------------------------------
function [] = set_axes()
% initialize axes
    ax = gui.h.viewAxes;
    set(ax, 'XAxisLocation', 'top', 'YDir', 'reverse', 'FontName', 'Consolas');
    xlabel_ = 'Distance (km)';
    if isinf(settings.vred)
        ylabel_ = 'Time (s)';
    else
        ylabel_ = sprintf('Time-Distance/%.2f (s)', settings.vred);
    end
    xlabel(ax, xlabel_, 'FontName', 'Consolas', 'FontSize', 11);
    ylabel(ax, ylabel_, 'FontName', 'Consolas', 'FontSize', 11);

    % apply view limit
    xlim(ax, gui.pp.vxlim);
    ylim(ax, gui.pp.vylim);
    % yticks(ax, settings.ylim(1):0.1:settings.ylim(2));

    % set x and y sliders
    if isfield(gui.h, 'xslider')
        delete(gui.h.xslider);
    end
    if abs(diff(gui.pp.vxlim)) < abs(diff(settings.xlim))
        gui.h.xslider = set_xslider(ax);
    end

    if isfield(gui.h, 'yslider')
        delete(gui.h.yslider);
    end
    if abs(diff(gui.pp.vylim)) < abs(diff(settings.ylim))
        gui.h.yslider = set_yslider(ax);
    end

    box on;
    ax.YAxis.TickDirection = 'out';

end

function [h] = set_xslider(ax)
    % get ax position in pixels
    old = get(ax, 'Units');
    set(ax, 'Units', 'pixels');
    axpos = ax.Position;
    set(ax, 'Units', old);

    xsliderpos = [axpos(1), axpos(2)-20-1, axpos(3), 20];
    min_val = settings.xlim(1);
    max_val = settings.xlim(2) - (gui.pp.vxlim(2) - gui.pp.vxlim(1));
    % how many views that the data limit can split into
    nview = ceil(diff(settings.xlim)/diff(gui.pp.vxlim));
    h = uicontrol(...
        'Style', 'slider', 'Parent', ax.Parent, 'Units', 'pixels', ...
        'Position', xsliderpos, 'SliderStep', [1/nview/10, 1/nview], ...
        'BackgroundColor', [220,220,220]/256, ...
        'Min', min_val, 'Max', max_val, 'Value', gui.pp.vxlim(1), ...
        'Callback', @scrollx);

    function [] = scrollx(src, ~)
        gui.pp.vxlim = gui.pp.vxlim + (src.Value - gui.pp.vxlim(1));
        xlim(ax, gui.pp.vxlim);
    end
end

function [h] = set_yslider(ax)
    % get ax position in pixels
    old = get(ax, 'Units');
    set(ax, 'Units', 'pixels');
    axpos = ax.Position;
    set(ax, 'Units', old);

    ysliderpos = [axpos(1)+axpos(3)+1, axpos(2), 20, axpos(4)];

    % y slider is should be reversed because the y axis is reversed
    min_val = settings.ylim(1);
    max_val = settings.ylim(2) - (gui.pp.vylim(2) - gui.pp.vylim(1));
    % how many views that the data limit can split into
    nview = ceil(diff(settings.ylim)/diff(gui.pp.vylim));
    h = uicontrol(...
        'Style', 'slider', 'Parent', ax.Parent, 'Units', 'pixels', ...
        'Position', ysliderpos, 'SliderStep', [1/nview/10, 1/nview], ...
        'BackgroundColor', [220,220,220]/256, ...
        'Min', min_val, 'Max', max_val, 'Value', (max_val - gui.pp.vylim(1)), ...
        'Callback', @scrolly);

    function [] = scrolly(src, ~)
        gui.pp.vylim = gui.pp.vylim + ((src.Max - src.Value) - gui.pp.vylim(1));
        ylim(ax, gui.pp.vylim);
        % yticks(ax, settings.ylim(1):0.1:settings.ylim(2));
    end
end


%% -----------------------------------------------------------------------------
function [] = redraw()
% replot according to selected data

    % get all selected raygroups
    selected = {};
    checked = gui.h.checkboxTree.CheckedNodes;
    for ii = 1:length(checked)
        selected = [selected, fun_get_end_nodes(checked(ii))];
    end

    % plot selected raygroups
    fun_plot_raygroups(selected);
    set_axes();
end

%% -----------------------------------------------------------------------------
function [res] = fun_get_end_nodes(node)
% get tree end nodes recursively.
    if isempty(node.Children)
        res = {node.UserData};
        return;
    end
    res = {};
    for ii = 1:length(node.Children)
        res = [res, fun_get_end_nodes(node.Children(ii))];
    end
end

%% -----------------------------------------------------------------------------
function [] = fun_plot_raygroups(raygroups)
% fun_plot_raygroups: plot raygroups
% @param: raygroups: a cell array containing ray group names. e.g. {'2.00-3.2', '2.00-4.2'}
    function [] = plot_raygroup(ax, data, xtraces, max_amp, name, color)
    % plot_raygroup: plot a single raygroup
        % get amplitude scale rate
        if isempty(settings.amp_lim) || settings.amp_lim == 0
            xinc = min(abs(diff(sort(xtraces))));
            settings.amp_lim = xinc * 0.5;
            settings.amp_lim = round(settings.amp_lim, 4);
        end
        scale_rate = settings.amp_lim / abs(max_amp);

        % max shift for vred
        max_xtrace = max(xtraces);
        max_shift = max_xtrace/settings.vred_orig - max_xtrace/settings.vred;

        for ii = 1:length(data)
            % plot a single trace
            xydata = data{ii};
            x = xydata(1, :);
            y = xydata(2, :);
            % shift for vred
            shift = xtraces(ii)/settings.vred_orig - xtraces(ii)/settings.vred;
            y = y + shift;

            % align v-reduce shift with zero amplitude wave
            x = [xtraces(ii), x, xtraces(ii)];
            if shift >= 0
                y = [settings.ylim(1), y, settings.ylim(2)+max_shift];
            else
                y = [settings.ylim(1)+max_shift, y, settings.ylim(2)];
            end

            % scale amplitude
            amp = x - xtraces(ii);
            amp = amp * scale_rate;
            x = xtraces(ii) + amp;
            curve = plot(ax, x, y, '-', 'Color', color_, 'LineWidth', 0.5, 'DisplayName', name);
            % only show legend for the first trace of ray group
            if ii ~= 1
                set(get(get(curve,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            end

            % fill peak and trough with colors
            if settings.fill_wave ~= 1
                % fill peak
                pos = find(amp >= 0);
                h = fill(ax, x(pos), y(pos), settings.peak_fill, 'EdgeColor', 'none');
                set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

                % fill trough
                neg = find(amp <= 0);
                h = fill(ax, x(neg), y(neg), settings.trough_fill, 'EdgeColor', 'none');
                set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            end
        end
    end

    % clear axes
    cla(gui.h.viewAxes);
    hold(gui.h.viewAxes, 'on');

    % plot raygroups
    % max_amp = fun_get_max_amplitude(gui.h.scaleEdit.String);
    max_amp = max(cellfun(@fun_get_max_amplitude, raygroups));
    for ii = 1:length(raygroups)
        color_ = 'k';
        if settings.enable_color
            idx = mod(ii, length(settings.colors));
            if idx == 0, idx = len(settings.colors); end
            color_ = settings.colors{idx};
        end
        idx = find(cellfun(@(s)strcmp(s, raygroups{ii}), obj.labels), 1);
        plot_raygroup(gui.h.viewAxes, obj.data{idx}, obj.xtraces{idx}, max_amp, raygroups{ii}, color_);
    end
    legend(gui.h.viewAxes, 'Location', 'northeast');
end

%% -----------------------------------------------------------------------------
function [] = on_expand_all(~, ~)
% on_expand_all: expand all nodes of checkboxTree
    nodes = gui.h.checkboxTree.Root.Children;
    for ii = 1:length(nodes)
        gui.h.checkboxTree.expandNode(nodes(ii));
    end
end

%% -----------------------------------------------------------------------------
function [] = on_collapse_all(src, event)
% on_collapse_all: collapse all nodes of checkboxTree
    nodes = gui.h.checkboxTree.Root.Children;
    for ii = 1:length(nodes)
        gui.h.checkboxTree.collapseNode(nodes(ii));
    end
end

%% -----------------------------------------------------------------------------
function [] = on_mouse_clicked(src, event)
% on_mouse_clicked: when mouse clicked. check if its single-click or double-click
    persistent chk
    if isempty(chk)
        chk = 1;
        pause(0.2); %Add a delay to distinguish single click from a double click
        if chk == 1
            chk = [];
            gui.pp.clickType = 'single';
        end
    else
        chk = [];
        gui.pp.clickType = 'double';
    end
end

%% -----------------------------------------------------------------------------
function [] = on_select_changed(~, event)
% on_node_select: when any tree node is selected
    if ~isempty(event.Nodes)
        node = event.Nodes(1);
        % gui.h.scaleEdit.String = node.UserData;
        pause(0.2);
        if strcmp(gui.pp.clickType, 'single'), return; end
        % when double click, fast plot the selected raygroup
        raygroups = fun_get_end_nodes(node);
        fun_plot_raygroups(raygroups);
    end
end

%% on_control_panel_help: help control panel
function [] = on_control_panel_help(~, ~)
    fprintf([
        'Help for control panel:\n\t1. Tick one or multiple ray groups to plot.\n\t', ...
        '2. Select one ray group to scale seismic amplitude to.\n\n']);
    input('Press <Enter> to continue ...');
end

%% -----------------------------------------------------------------------------
function [] = on_enable_color(~, ~)
% on_enable_color: enable color
    settings.enable_color = true;
end

%% -----------------------------------------------------------------------------
function [] = on_disable_color(~, ~)
% on_enable_color: disable color
    settings.enable_color = false;
end

%% -----------------------------------------------------------------------------
function [] = on_apply(~, ~)
% on_apply: on apply button clicked
    redraw();
end

function [] = on_readme(~, ~)
% on_readme: show readme
    !start gui_pltsyn_readme.txt
end

%% -----------------------------------------------------------------------------
function [] = on_export_plot(~, ~)
% on_export_plot: export axes to a new figure
    if ~isfield(gui.h, 'exportWindow') || ~ishandle(gui.h.exportWindow)
        gui.h.exportWindow = figure();
    else
        figure(gui.h.exportWindow);
    end
    newAxes = copyobj(gui.h.viewAxes, gui.h.exportWindow);
    % The original position is copied too, so adjust it.
    set(newAxes, 'Units', 'normalized', 'Position', get(groot, 'DefaultAxesPosition'));
end

%% -----------------------------------------------------------------------------
function [scaleby] = fun_get_max_amplitude(raygroup)
% fun_get_max_amplitude: get the max amplitude of specific raygroup(s)
    scaleby = 0;
    idx = find(cellfun(@(lb)startsWith(lb, raygroup), obj.labels));
    xtraces = obj.xtraces(idx);
    data = obj.data(idx);
    for ii = 1:length(xtraces)
        for jj = 1:length(xtraces{ii})
            xydata = data{ii}{jj};
            xtrace = xtraces{ii}(jj);
            max_amp = max(abs(xydata(1, :) - xtrace));
            scaleby = max(scaleby, max_amp);
        end
    end
end

%% -----------------------------------------------------------------------------


%% Settings window
%% -----------------------------------------------------------------------------
function [] = on_settings(src, event)
% on_settings: open Settings window
    if isfield(gui.h, 'settingsWindow') && ishandle(gui.h.settingsWindow)
        figure(gui.h.settingsWindow);
        return;
    end

    gui.h.settingsWindow = figure(...
        'Name', 'Settings', 'NumberTitle', 'off', 'MenuBar', 'none', ...
        'ToolBar', 'none');

    layout = uix.VBox('Parent', gui.h.settingsWindow, 'Padding', 10, 'Spacing', 10);
    mainLayout = uix.HBox('Parent', layout, 'Padding', 0, 'Spacing', 30);
    buttonArea = uix.HBox('Parent', layout, 'Padding', 0, 'Spacing', 8);
    set(layout, 'Heights', [-1, 28]);

    leftLayout = uix.VBox('Parent', mainLayout, 'Padding', 0, 'Spacing', 10);
    rightLayout = uix.VBox('Parent', mainLayout, 'Padding', 0, 'Spacing', 10);
    set(mainLayout, 'Widths', [300, 200]);

    gui.h.st.xlimText = uiw.widget.EditableText(...
        'Parent', leftLayout, 'Value', mat2str(settings.xlim), ...
        'Label', 'X Axis Range', 'LabelLocation', 'left', 'LabelWidth', 100,...
        'Callback', @on_num_pair_edited);

    gui.h.st.ylimText = uiw.widget.EditableText(...
        'Parent', leftLayout, 'Value', mat2str(settings.ylim), ...
        'Label', 'Y Axis Range', 'LabelLocation', 'left', 'LabelWidth', 100,...
        'Callback', @on_num_pair_edited);

    gui.h.st.vredText = uiw.widget.EditableText(...
        'Parent', leftLayout, 'Value', mat2str(settings.vred), ...
        'Label', 'v-reduce', 'LabelLocation', 'left', 'LabelWidth', 100,...
        'Callback', @on_vred_edited);

    gui.h.st.ntpiText = uiw.widget.EditableText(...
        'Parent', leftLayout, 'Value', mat2str(settings.ntpi), ...
        'Label', 'Traces per Inch', 'LabelLocation', 'left', 'LabelWidth', 100,...
        'Callback', @on_num_edited);

    gui.h.st.mspiText = uiw.widget.EditableText(...
        'Parent', leftLayout, 'Value', mat2str(settings.mspi), ...
        'Label', 'T(ms) per Inch', 'LabelLocation', 'left', 'LabelWidth', 100,...
        'Callback', @on_num_edited);

    gui.h.st.amplimText = uiw.widget.EditableText(...
        'Parent', leftLayout, 'Value', mat2str(settings.amp_lim), ...
        'Label', 'Max Amplitude', 'LabelLocation', 'left', 'LabelWidth', 100,...
        'Callback', @on_num_edited);

    gui.h.st.fillWavePopup = uiw.widget.Popup(...
        'Parent', leftLayout, 'Items', settings.fill_wave_options, ...
        'Label', 'Fill Waveform', 'LabelLocation', 'left', 'LabelWidth', 100,...
        'SelectedIndex', settings.fill_wave, 'Callback', @on_fill_wave_selected);

    gui.h.st.fillPeakColorSelector = uiw.widget.ColorSelector(...
        'Parent', leftLayout, 'Label', 'peak-fill', 'LabelWidth', 100, ...
        'LabelHorizontalAlignment', 'right', ...
        'Value', settings.peak_fill, 'Enable', 'off');
    gui.h.st.fillTroughColorSelector = uiw.widget.ColorSelector(...
        'Parent', leftLayout, 'Label', 'trough-fill', 'LabelWidth', 100, ...
        'LabelHorizontalAlignment', 'right', ...
        'Value', settings.trough_fill, 'Enable', 'off');

    uix.Empty('Parent', leftLayout);
    set(leftLayout, 'Heights', [ones(1, 9) * 25, -1]);

    gui.h.st.colorList = uiw.widget.ListWithButtons(...
        'Parent', rightLayout, ...
        'Items', cellfun(@mat2str, settings.colors, 'UniformOutput', false), ...
        'AllowAdd', true, ... %Requires callback implementation
        'AllowCopy', false, ... %Requires callback implementation
        'AllowDelete', true, ... %Requires callback implementation
        'AllowEdit', false, ... %Requires callback implementation
        'AllowMove', true, ... %Callback is optional
        'AllowPlot', false, ... %Requires callback implementation
        'AllowReverse', false, ... %Requires callback implementation
        'AllowRun', false, ... %Requires callback implementation
        'Callback', @color_list_callback, ...
        'ButtonLocation', 'right', ...
        'Label', 'Colors for plotting', ...
        'LabelLocation', 'top', ...
        'LabelHeight', 18);

    gui.h.st.colorSelector = uiw.widget.ColorSelector(...
        'Parent', rightLayout, ...
        'Value', settings.colors{1}, ...
        'Callback', @on_color_edited, ...
        'LabelVisible', 'off', ...
        'Label', '', ...
        'LabelLocation', 'right', ...
        'LabelWidth', 30);

    uix.Empty('Parent', rightLayout);
    set(rightLayout, 'Heights', [250, 25, -1]);

    uix.Empty('Parent', buttonArea);
    uicontrol('Style', 'PushButton', 'Parent', buttonArea, 'String', 'OK', 'Callback', @on_settings_ok);
    uicontrol('Style', 'PushButton', 'Parent', buttonArea, 'String', 'Cancel', 'Callback', @on_settings_cancel);
    uicontrol('Style', 'PushButton', 'Parent', buttonArea, 'String', 'Apply', 'Callback', @on_settings_apply);
    set(buttonArea, 'Widths', [-1, 80, 80, 80]);
end

%% -----------------------------------------------------------------------------
function [] = on_num_edited(src, event)
% when EditableText change, validate the value
    if strcmp(event.Interaction, 'Edit')
        is_valid = validate_num(event.NewValue);
        if ~is_valid
            fun_on_invalid_value(src, event);
            src.Value = event.OldValue;
        end
    end
end

%% -----------------------------------------------------------------------------
function [] = on_num_pair_edited(src, event)
% when EditableText change, validate the value
    if strcmp(event.Interaction, 'Edit')
        is_valid = validate_num_pair(event.NewValue);
        if ~is_valid
            fun_on_invalid_value(src, event);
            src.Value = event.OldValue;
        end
    end
end

%% -----------------------------------------------------------------------------
function [] = on_vred_edited(src, event)
% when EditableText change, validate the value
    on_num_edited(src, event);
    if strcmp(event.Interaction, 'Edit')
        if str2num(event.NewValue) == 0
            msg = sprintf('Vred can NOT be zero, maybe you want "Inf"');
            fun_on_invalid_value(src, event, msg);
            src.Value = 'Inf';
        end
    end
end

%% -----------------------------------------------------------------------------
%% fun_on_invalid_value: when EditableText get invalid value
function [] = fun_on_invalid_value(src, event, msg)
    if nargin < 3
        msg = sprintf('Invalid value for "%s" = %s', strip(strip(src.Label),'right',':'), src.Value);
    end
    h = errordlg(msg, 'Value Error');
    % h.Position(3:4) = [250, 70];
    htext = findobj(h, 'Type', 'Text');
    htext.FontSize = 9;
    set(h, 'Resize', 'on');
end

%% -----------------------------------------------------------------------------
function [] = on_fill_wave_selected(src, event)
% when fillWavePopup selection
    enable = 'off';
    if strcmp(lower(src.Value), lower('NO'))
        set(gui.h.st.fillPeakColorSelector, 'Value', [1, 1, 1]);
        set(gui.h.st.fillTroughColorSelector, 'Value', [1, 1, 1]);
    elseif strcmp(lower(src.Value), lower('Black-White'))
        set(gui.h.st.fillPeakColorSelector, 'Value', [0, 0, 0]);
        set(gui.h.st.fillTroughColorSelector, 'Value', [1, 1, 1]);
    elseif strcmp(lower(src.Value), lower('Red-Blue'))
        set(gui.h.st.fillPeakColorSelector, 'Value', [1, 0, 0]);
        set(gui.h.st.fillTroughColorSelector, 'Value', [0, 0, 1]);
    elseif strcmp(lower(src.Value), lower('Custom'))
        enable = 'on';
    end

    set(gui.h.st.fillPeakColorSelector, 'Enable', enable);
    set(gui.h.st.fillTroughColorSelector, 'Enable', enable);
end

%% -----------------------------------------------------------------------------
function [] = color_list_callback(src, event)
% handle color list event
    % color selected
    if strcmp(event.Interaction, 'Selection')
        gui.h.st.colorSelector.Value = src.SelectedItems{1};

    % delete color
    elseif strcmp(event.Interaction, 'Delete')
        idx = src.SelectedIndex;
        src.Items(src.SelectedIndex) = [];
        if idx > length(src.Items), idx = idx - 1; end
        src.SelectedIndex = idx;

    % add color to the end of list
    elseif strcmp(event.Interaction, 'Add')
        src.Items{end+1} = fun_ensure_str(gui.h.st.colorSelector.Value);
        src.SelectedIndex = length(src.Items);
    end
end

%% -----------------------------------------------------------------------------
function [] = on_color_edited(src, event)
% when color is edited in ColorSelector
    if strcmp(event.Interaction, 'Edit')
        src.Value = [0, 0, 0];
        src.Value = event.NewValue;
    end
    idx = gui.h.st.colorList.SelectedIndex;
    gui.h.st.colorList.Items{idx} = fun_ensure_str(src.Value);
    % avoid auto-resetting SelectedIndex
    gui.h.st.colorList.SelectedIndex = idx;
end

%% -----------------------------------------------------------------------------
function [] = on_settings_ok(~, ~)
% when settings OK, save values and close settings window
    on_settings_apply();
    close(gui.h.settingsWindow);
end

%% -----------------------------------------------------------------------------
function [] = on_settings_cancel(~, ~)
% when canceling settings, restore values and close settings window
    close(gui.h.settingsWindow);
end

%% -----------------------------------------------------------------------------
function [] = on_settings_apply(~, ~)
% when applying settings, save values and keep settings window alive
    settings.colors = cellfun(@eval, gui.h.st.colorList.Items, 'UniformOutput', false);
    settings.xlim = str2num(gui.h.st.xlimText.Value);
    settings.ylim = str2num(gui.h.st.ylimText.Value);
    settings.vred = str2num(gui.h.st.vredText.Value);
    settings.ntpi = str2num(gui.h.st.ntpiText.Value);
    settings.mspi = str2num(gui.h.st.mspiText.Value);
    settings.amp_lim = str2num(gui.h.st.amplimText.Value);
    settings.fill_wave = gui.h.st.fillWavePopup.SelectedIndex;
    settings.peak_fill = gui.h.st.fillPeakColorSelector.Value;
    settings.trough_fill = gui.h.st.fillTroughColorSelector.Value;
    % update view limit
    fun_calc_view_limit();
    redraw();
end

function [] = fun_calc_view_limit()
% fun_calc_view_limit: calculate view limit of x and y directory
    % get ax position in inches
    ax = gui.h.viewAxes;
    old = get(ax, 'Units');
    set(ax, 'Units', 'inches');
    axpos = ax.Position;
    set(ax, 'Units', old);

    % calculate view length
    width_inch = axpos(3);
    height_inch = axpos(4);
    xlen = width_inch * settings.ntpi * settings.xinc;
    ylen = height_inch * settings.mspi / 1e3;

    % if view limit is not set yet, use data limit as default
    if isfield(gui.pp, 'vxlim') && ~isempty(gui.pp.vxlim)
        old_vxlim = gui.pp.vxlim;
        old_vylim = gui.pp.vylim;
    else
        old_vxlim = settings.xlim;
        old_vylim = settings.ylim;
    end

    % calculate x and y view limit by ntpi and mspi respectively
    % apply left-alignment on new view and old view
    gui.pp.vxlim = old_vxlim(1) + [0, xlen];
    gui.pp.vylim = old_vylim(1) + [0, ylen];

    % % when view limit out of range, adjust it by upper limit
    % if gui.pp.vxlim(2) > settings.xlim(2)
    %     shift = gui.pp.vxlim(2) - settings.xlim(2);
    %     gui.pp.vxlim = gui.pp.vxlim - shift;
    %     gui.pp.vxlim(1) = max(gui.pp.vxlim(1), settings.xlim(1));
    % end
    % if gui.pp.vylim(2) > settings.ylim(2)
    %     shift = gui.pp.vylim(2) - settings.ylim(2);
    %     gui.pp.vylim = gui.pp.vylim - shift;
    %     gui.pp.vylim(1) = max(gui.pp.vylim(1), settings.ylim(1));
    % end
end

%% -----------------------------------------------------------------------------
function [value] = fun_ensure_str(value)
% ensure a value is char array or string type
    if ~ischar(value) && ~isstring(value)
        value = mat2str(value);
    end
end

%% -----------------------------------------------------------------------------
function [is_valid] = validate_num(value)
% check if a string/char can be parsed into a num
    value = str2num(value);
    is_valid = ~isempty(value);
end

%% -----------------------------------------------------------------------------
function [is_valid] = validate_num_pair(value)
% check if a string/char can be parsed into a pair of nums(2 nums)
    value = str2num(value);
    is_valid = false;
    if ~isempty(value) && numel(value) == 2
        is_valid = true;
    end
end

end

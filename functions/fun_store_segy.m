function fun_store_segy(figs, working_dir)
% 将每个 shot 的地震数据保存到一个 segy 文件
% 依赖 Seislab 包：https://www.mathworks.com/matlabcentral/fileexchange/53109-seislab-3-02

    % Seislab 包初始化
    presets;

    path_rin = fun_trans_rin2m(fullfile(working_dir, 'r.in'));
    path_sin = fun_trans_rin2m(fullfile(working_dir, 's.in'));
    clear(path_rin);
    clear(path_sin);
    s_rin = run2struct(path_rin);
    s_sin = run2struct(path_sin);

    % 地震数据头信息
    base = struct();
    base.type = 'seismic';
    base.tag = 'unspecified';

    % 时间轴参数
    base.first = s_sin.tmin * 1000;
    base.last = s_sin.tmax * 1000;
    base.step = 2;
    base.units = 'ms';
    time_axis = (base.first:base.step:base.last)';

    % header 信息
    base.header_info = { ...
        'ds_seqno', 'n/a', 'Trace sequence number within line'; ...
        'cdp', 'n/a', 'CDP sequence number'; ...
        'cdp_x', 'm', 'X-coordinate of CDP'; ...
        'cdp_y', 'm', 'Y-coordinate of CDP'; ...
        'iline_no', 'n/a', 'In-line number'; ...
        'xline_no', 'n/a', 'Cross-line number'; ...
        'offset', 'm', 'Offset'; ...
    };

    % 接收点 X 坐标转化为以米为单位
    xreceivers = int32(s_rin.xmins * 1000):int32(s_rin.xincs * 1000):int32(s_rin.xmaxs * 1000);
    ntraces = length(xreceivers);
    base.headers = [ ...
        1:ntraces; ...
        1:ntraces; ...
        zeros(1, ntraces); ...
        xreceivers; ...
        ones(1, ntraces); ...
        1:ntraces; ...
        xreceivers; ...
    ];

    % write_segy_file 的额外参数，用于指定如何保存 headers
    segy_header_params = {'headers',{'cdp_x',73,4},{'cdp_y',77,4},{'iline_no',181,4},{'xline_no',185,4}};

    % 从 figure 对象中提取数据
    all_curves = [];
    for ii = 1:numel(figs)
        fig = figs(ii);
        all_curves = [all_curves; get(get(fig, 'Children'), 'Children')];
    end

    traces = {};
    labels = {};
    xtraces = [];
    for ii = 1:numel(all_curves)
        curve = all_curves(ii);
        if isempty(curve.UserData), continue; end

        % 需要记录所有 trace 的 x 坐标，用于发现和补齐某些没有震动的 trace。
        xtraces(end+1) = curve.UserData.xtrace;
        labels{end+1} = curve.UserData.tag;

        % 按照地震数据纵轴的要求，对 curve 进行线性插值。同时去掉偏移量。
        traces{end+1} = interp1(curve.YData * 1000, curve.XData - curve.UserData.xtrace, time_axis);
    end
    xtraces = int32(xtraces * 1000);

    % 按照炮点以及事件分组
    [labels, idx] = sort(labels);
    traces = traces(idx);
    xtraces = xtraces(idx);
    g = findgroups(labels);
    labels = unique(labels);
    traces = splitapply(@(c) {c}, traces, g);
    xtraces = splitapply(@(c) {c}, xtraces, g);

    % 将 traces 按照其 x 坐标排序。
    % 同时，将缺失的（没有震动） trace，以 0 振幅补齐。
    default_trace = zeros(size(time_axis));
    for ii = 1:numel(traces)
        missed = setdiff(xreceivers, xtraces{ii});
        [~, idx] = sort([xtraces{ii}, missed]);
        ts = [traces{ii}, repmat({default_trace}, 1, length(missed))];
        traces{ii} = ts(idx);
    end

    % 保存数据
    % 为了在导入 HRS 时不用再重新命名数据集的名称，使用路径的一部分作为 segy 文件的前缀，以示区分
    splited_path = strsplit(working_dir, filesep());
    file_prefix = [strjoin(splited_path(end-1:end), '-'), '-'];
    outdir = fullfile(working_dir, 'segy');
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end
    for ii = 1:numel(labels)
        label = labels{ii};

        % % 只保存叠加后的数据，忽略单个事件的数据
        % if ~endsWith(label, 'stacked')
        %     continue;
        % end

        % 根据炮点坐标计算接收点的 offset
        tmp = strsplit(label, '-');
        shot_position = int32(str2double(tmp{1}) * 1000);
        offset = xreceivers - shot_position;

        seismic = struct(base);
        seismic.name = label;
        seismic.traces = cell2mat(traces{ii});
        seismic.headers(end, :) = offset;

        file = fullfile(outdir, [file_prefix, label, '.sgy']);
        write_segy_file(seismic, file, segy_header_params);
    end
end

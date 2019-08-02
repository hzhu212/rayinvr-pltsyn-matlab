function fun_store_segy(figs, working_dir)
% 将每个 shot 的地震数据保存到一个 segy 文件
% 依赖 Seislab 包：https://www.mathworks.com/matlabcentral/fileexchange/53109-seislab-3-02

    % Seislab 包初始化
    presets;

    path_rin = fullfile(working_dir, 'r.in');
    path_sin = fullfile(working_dir, 's.in');
    s_rin = run2struct(fun_trans_rin2m(path_rin));
    s_sin = run2struct(fun_trans_rin2m(path_sin));

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
        'ds_seqno', NaN, 'Trace sequence number within line'; ...
        'depth', 'm', 'Source depth below surface'; ...
        'offset', 'm', 'offset'; ...
    };
    offset = s_rin.xmins:s_rin.xincs:s_rin.xmaxs;
    base.headers = [ ...
        1:length(offset); ...
        ones(1, length(offset)) * 1201; ...
        offset * 1000; ...
    ];

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
        traces{end+1} = interp1(curve.YData * 1000, curve.XData, time_axis) - curve.UserData.xtrace;
    end

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
        missed = setdiff(offset, xtraces{ii});
        [~, idx] = sort([xtraces{ii}, missed]);
        ts = [traces{ii}, repmat({default_trace}, 1, length(missed))];
        traces{ii} = ts(idx);
    end

    % 保存数据
    outdir = fullfile(working_dir, 'segy');
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end
    for ii = 1:numel(labels)
        % % 只保存叠加后的数据，忽略单个事件的数据
        % if ~endsWith(labels{ii}, 'stacked')
        %     continue;
        % end

        file = fullfile(outdir, [labels{ii}, '.sgy']);
        seismic = struct(base);
        seismic.name = labels{ii};
        seismic.traces = cell2mat(traces{ii});
        write_segy_file(seismic, file, {'print', 0});
    end
end

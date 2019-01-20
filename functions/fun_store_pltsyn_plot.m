function fun_store_pltsyn_plot(figs, working_dir)
% store plot data to a .mat file

    file = fullfile(working_dir, 'plotdata.pltsyn.mat');
    all_curves = [];
    for ii = 1:length(figs)
        fig = figs(ii);
        all_curves = [all_curves; get(get(fig, 'Children'), 'Children')];
    end

    fig1 = get(figs(1));
    % vred = fig1.UserData.vred;

    data = {};
    labels = {};
    xtraces = [];
    for ii = 1:length(all_curves)
        curve = all_curves(ii);
        if isempty(curve.UserData), continue; end
        xdata = curve.XData;
        ydata = curve.YData;
        xtrace = curve.UserData.xtrace;
        % % We should store original un-reduced time, although the plot used reduced data.
        % if ~isempty(vred) && vred ~= 0
        %     ydata = ydata + xtrace/vred;
        % end
        data{end+1} = [xdata; ydata];
        xtraces(end+1) = xtrace;
        labels{end+1} = curve.UserData.tag;
    end

    % group data by labels
    [labels, idx] = sort(labels);
    data = data(idx);
    xtraces = xtraces(idx);
    g = findgroups(labels);
    labels = unique(labels);
    data = splitapply(@(c) {c}, data, g);
    xtraces = splitapply(@(c) {c}, xtraces, g);

    % save all data and information to a struct
    obj = struct();
    obj.labels = labels;
    obj.xtraces = xtraces;
    % trace interval
    obj.xinc = min(abs(diff(sort(xtraces{1}))));
    obj.data = data;
    obj.xlabel = fig1.CurrentAxes.XLabel.String;
    obj.ylabel = fig1.CurrentAxes.YLabel.String;
    obj.xlim = fig1.CurrentAxes.XAxis.Limits;
    obj.ylim = fig1.CurrentAxes.YAxis.Limits;
    obj.vred = fig1.UserData.vred;

    save(file, 'obj');
end

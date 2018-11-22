function fun_store_plotdata(figs, working_dir)
% store plot data to a .mat file

    file = fullfile(working_dir, 'plotdata.pltsyn.mat');
    all_curves = [];
    for ii = 1:length(figs)
        fig = figs(ii);
        all_curves = [all_curves, get(get(fig, 'Children'), 'Children')];
    end
    all_curves = all_curves(:);

    data = {};
    labels = {};
    xtraces = [];
    for ii = 1:length(all_curves)
        curve = all_curves(ii);
        if isempty(curve.UserData), continue; end
        data{end+1} = [curve.XData; curve.YData];
        xtraces(end+1) = curve.UserData.xtrace;
        labels{end+1} = curve.UserData.tag;
    end
    data{end+1} = xtraces;
    data{end+1} = strjoin(labels, '|');
    data = data';
    save(file, 'data');
end

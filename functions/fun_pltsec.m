function [figs] = fun_pltsyn(vred, xomit, nskip, ipol, tol, nptsw, itrev, idump, iamp, twin, imeth, iroute)
% pltsec routine for PLTSYN
% plot synthetic seismic sections using the time, amplitude,
% and phase of each arrival read in from unit 11

    global path_pltsyn_par path_pltsyn_com;
    global path_txout;
    global fid_11 fid_12 fid_18 fid_20;
    run(path_pltsyn_par);
    run(path_pltsyn_com);

    % declare variables
    figs = [];
    namp = 0;

    % real
    dist_ = zeros(1, pseis);
    s = zeros(1, ppseis);
    r = zeros(1, ppseis+121);
    time = zeros(1, ppseis);
    sect = zeros(pseis, parriv, 3);
    wavlet = zeros(1, pwavlt);
    xamp = zeros(1, parriv*pseis);
    tamp = zeros(1, parriv*pseis);
    hil = zeros(1, 121);
    sp = zeros(1, ppseis);
    tp = zeros(1, ppseis);

    % integer
    na = int32(zeros(1, pseis));
    ipamp = int32(zeros(1, parriv*pseis));

    wavlet = [...
        .01,.072,.121,.12,.044,-.113,-.32,-.512,-.607,-.538,-.283,.112,.542,...
        .875,1.0,.875,.542,.112,-.283,-.538,-.607,-.512,-.32,-.113,.044,.12,...
        .121,.072,.01, repelem(0., 71)];

    nseg = 0;
    p121 = ppseis + 121;
    pi2 = pi * 2;
    dens = dens / 2;
    fifi = double(ifill);

    % NEW
    % sample rate(s), 0.001 -- 1ms
    sample_rate = 0.003;
    npts = floor((tmax - tmin) / sample_rate) + 1;

    nptsp = npts + 121;
    nomit = 0;
    if vred ~= 0, rvred = 1 / vred;
    else, rvred = 0; end
    % nwin = round(twin * sps);

    % 40
    nomit = find(xomit(1:pseis) < -99998, 1);
    if isempty(nomit), nomit = pseis;
    else, nomit = nomit - 1; end

    % 50
    % 80
    time(1:npts) = tmin:sample_rate:tmax;

    % 90
    for ii = 2:2:60
        fudge = ii * 0.00018604 - 0.00037207;
        hil(61+ii-1) = 0.6366198 / (ii-1) - fudge;
        hil(61-ii+1) = -hil(61+ii-1);
    end

    if iconv == 1
        if iwavlt < 2
            if iwavlt == 0
                nwave = 29;
                nwave1 = 30;
            else
                nwave = nptsw;
                nwave1 = nwave + 1;
                % 430
                wavlet = sin(2*pi/(nptsw-1) .* ((1:nptsw)-1));
            end
        else
            nwave1 = 1;
            % 450
            wavlet = fscanf(fid_20, '%f');
            nwave = length(wavlet);
            % 440
            nwave1 = nwave + 1;
            if nwave == 0
                fprintf('***  w.in is empty  ***\n');
                return;
            end
        end

        % 445
        wsum = sum(wavlet(1:nwave));
        wmean = wsum / nwave;
        % 455
        wavlet(1:nwave) = wavlet(1:nwave) - wmean;
        % 465
        wmax = max(abs(wavlet(1:nwave)));
        if wmax > 0.000001
            % 475
            wavlet(1:nwave) = wavlet(1:nwave) / wmax;
        else
            fprintf('***  source wavelet is zero  ***\n');
            return;
        end
        if abs(wavlet(1) < 0.001)
            wavlet(1) = 0.001 * sign(wavlet(1));
        end
        if abs(wavlet(nwave) < 0.001)
            wavlet(nwave) = 0.001 * sign(wavlet(nwave));
        end
        if ipol == -1
            % 130
            wavlet(1:nwave) = -wavlet(1:nwave);
        end
        twave = (nwave - 1) / sps;
    end

    % 1000: for each xshot
    tline = fgetl(fid_11);
    while ischar(tline)
        tmp = sscanf(tline, '%10f%10d');
        assert(length(tmp) == 2 && abs(tmp(2)) == 1, '*** unexpected error in reading sect_ext.out ***');
        if isempty(tmp)
            fprintf('***  insufficient data in sect_ext.out  ***\n');
            return;
        end
        [xshot, idum] = deal(tmp(1), tmp(2));
        if iplot >= -1
            fig_current = figure('numbertitle', 'off', 'name', sprintf('xshot=%.3f km', xshot));
            figs(end+1) = fig_current;
            set(fig_current, 'Position', [300, 200, 800, 500]);
            axes('FontName', 'Consolas');
            set(gca(), 'XAxisLocation', 'top', 'YDir', 'reverse');
            xlabel('DISTANCE (km)', 'FontName', 'Consolas', 'FontSize', 11);
            if vred == 0
                nchart = 8;
                tlab = 'TIME (s)';
            else
                tlab = ['T-D/', sprintf('%.2f', vred), ' (s)'];
                nchart = length(tlab);
            end
            ylabel(tlab, 'FontName', 'Consolas', 'FontSize', 11);
            xlim([xmin, xmax]);
            ylim([tmin, tmax]);
            title(sprintf('xshot=%.3f km', xshot));
            box on;
            iplots = 1;
        end

        maxamp = 0;
        ntplt = 0;
        nsexsp = 0;
        nhpts = npts;
        if iplot >= -1, cla(); end

        nseis = 0;

        % 991
        xshotn = -999999;

        % 180: for each seismogram
        tline = fgetl(fid_11);
        while ischar(tline)
            tmp = sscanf(tline, '%10f%10d');
            [dists, nas] = deal(tmp(1), tmp(2));
            if nas < 0
                xshotn = dists;
                % go to 99
                break;
            end
            nseis = nseis + 1;
            dist_(nseis) = dists;
            na(nseis) = nas;
            if na(nseis) > parriv
                % 555
                fprintf('***  max # of arrivals on seismogram at%10.3f km exceeded ***\n', dist_(nseis));
                return;
            end
            if na(nseis) > 0
                % 190: for each arrival
                for jj = 1:na(nseis)
                    tline = fgetl(fid_11);
                    tmp = sscanf(tline, '%12f%12f%12f%4f');
                    assert(length(tmp) == 4, '***  insufficient data in sect_ext.out  ***\n');
                    sect(nseis, jj, 1:4) = tmp;
                    if vred > 0
                        sect(nseis, jj, 1) = sect(nseis, jj, 1) - abs(dist_(nseis) - xshot) / vred;
                    end
                end
            end
            % go to 180
            tline = fgetl(fid_11);
        end

        if iamp > 0
            [xrs, trs, urs, irs] = fun_load_txin(path_txout);
            idx_shots = find(irs <= 0);
            shots = xrs(idx_shots);
            idx_cur_shot = idx_shots(find(abs(shots - xshot) < 1e-4, 1));
            idx_next_shot = idx_shots(find(idx_shots == idx_cur_shot, 1) + 1);
            namp = idx_next_shot - idx_cur_shot - 1;
            if namp <= 0
                fprintf('***  error in file tx.out  ***\n');
                return;
            else
                fprintf(fid_18, '%10.3f%10.3f%10.3f%10d\n', xshot, 0, 0, 0);
            end
            xamp = xrs((idx_cur_shot+1):(idx_next_shot-1));
            tamp = trs((idx_cur_shot+1):(idx_next_shot-1));
            ipamp = irs((idx_cur_shot+1):(idx_next_shot-1));
        end
        idir = 1;

        % length of these arrays or cells will be equal to `nseis`
        traces = {};
        xreceivers = [];
        raycodes = [];

        % 70
        % for each seismogram
        for ii = 1:nseis
            npts = nhpts;
            if mod(ii, nskip) ~= 0, continue; end
            if dist_(ii) < xmin || dist_(ii) > xmax, continue; end
            if nomit > 0
                % 75
                if any(xomit(1:nomit) - dist_(ii) < 0.001), continue; end
            end
            nsexsp = nsexsp + 1;

            % 270
            % for each event
            % each seismogram will receive several events, like: 3.2, 4.2, ...
            for jj = 1:na(ii)
                % 170
                r(1:nptsp) = 0;

                if iconv == 1 && (sect(ii,jj,1)+twave < tmin || sect(ii,jj,1) > tmax)
                    continue;
                end
                cur_raycode = sect(ii,jj,4);
                % nstart = round((sect(ii,jj,1) - tmin) * sps);
                nstart = round((sect(ii,jj,1) - tmin) / sample_rate);
                sn = sin(sect(ii,jj,3)) * sect(ii,jj,2);
                kpos = nstart + 61;
                if kpos > 0 && kpos <= nptsp
                    r(kpos) = r(kpos) + cos(sect(ii,jj,3)) * sect(ii,jj,2);
                end
                % TODO
                % 370
                for k = 2:2:120
                    kpos = nstart + k;
                    if kpos > 0 && kpos <= nptsp
                        r(kpos) = r(kpos) - hil(k) * sn;
                    end
                end

                if nsmth > 0
                    % 280
                    for jj = 1:nsmth
                        [r, ~] = fun_smooth(r, nptsp);
                    end
                end
                if iconv == 1
                    % TODO
                    % 470
                    for jj = 1:npts
                        if jj <= nwave - 60
                            kmax = jj - 59;
                        else
                            kmax = nwave;
                        end
                        % 570
                        k = 1:kmax;
                        s(jj) = sum(wavlet(k) .* r(jj-k+60));
                    end
                else
                    % 460
                    s(1:npts) = r((1:npts) + 60);
                end

                if inmo == 1
                    s = fun_nmo(s, npts, sps, tmin, dist_(ii), vrms);
                end

                if idump == 1
                    % FIXME
                    % 35
                    fprintf(fid_12, '%4d', s(1:npts));
                end

                if iplot < -1, continue; end
                if tol > 0 && npts > 2
                    % find the points that need to be plot: the amplitude should reach the tolerance
                    idx_amp = find(abs(s(1:npts)) >= tol);
                    % find all the boundary of the points that need to be plot. set their amplitude to zero,
                    % in case the tol is so large as to make a tilt straight line between the ploting-point
                    % and unploting-point.
                    idx_zero = setdiff([idx_amp-1, idx_amp+1], [idx_amp, 0, npts+1]);
                    idx_plot = sort([idx_amp, idx_zero]);
                    % add start point and end point
                    idx_plot = unique([1, idx_plot, npts]);
                    nplt = length(idx_plot);
                    sp(1:nplt) = s(idx_plot);
                    tp(1:nplt) = time(idx_plot);
                else
                    nplt = npts;
                    % 710
                    sp(1:nplt) = s(1:nplt);
                    tp(1:nplt) = time(1:nplt);
                end
                sp = sp + dist_(ii);

                ntplt = ntplt + nplt;
                if idir == 1
                    nl1 = 1;
                    nl2 = 2;
                    nl3 = nplt;
                    istep = 1;
                else
                    nl1 = nplt;
                    nl2 = nplt - 1;
                    nl3 = 1;
                    istep = -1;
                end
                % 620
                traces{end+1} = [reshape(tp(nl1:istep:nl3), [], 1), reshape(sp(nl1:istep:nl3), [], 1)];
                xreceivers(end+1) = dist_(ii);
                raycodes(end+1) = cur_raycode;
            end
            % end
            idir = -idir;
        end

        % plot seismic records for each shot
        plot_per_shot(fig_current, xshot, traces, xreceivers, raycodes);

        if iplot >= -1 && nsexsp > 0 && tol > 0
            per = ntplt / (nsexsp * npts) * 100;
            % 105
            fprintf('***  %6.2f%% of points plotted using spmin=%9.3e  ***\n', per, tol);
        end

        % 665
        fprintf('>>>  shot at %10.3f km completed  <<<\n', xshot);
        if xshotn > -999998
            xshot = xshotn;
        end
        % go to 1000
    end

    if namp > 0
        fprintf(fid_18, '%10.3f%10.3f%10.3f%10d\n', 0, 0, 0, -1);
    end
end


function plot_per_shot(fig, xshot, traces, xreceivers, raycodes)
% plot seismic traces for each shot
% plot each event separatly with different colors

    % TODO: different shots and events have different number of receivers

    colors = ['r', 'g', 'b', 'c', 'm', 'y'];
    all_codes = unique(raycodes);
    all_xreceivers = sort(unique(xreceivers));

    per_event_data = {};

    for ii = 1:length(all_codes)
        code = all_codes(ii);
        color_ = colors(mod(ii,length(colors)) + 1);
        tag = sprintf('%.4f-%.1f', xshot, code);

        % find traces for this event
        idx = raycodes == code;
        use_xreceivers = xreceivers(idx);
        use_traces = traces(idx);

        % sort traces by x position of receivers
        [use_xreceivers, sorted_idx] = sort(use_xreceivers);
        use_traces = use_traces(sorted_idx);

        plot_per_shot_event(fig, tag, color_, use_traces, use_xreceivers);

        per_event_data{end+1} = use_traces;
    end

    %% calculate stacked
    % collect traces at each x receiver
    traces_stacked = per_event_data{1};
    for ii = 2:length(per_event_data)
        traces = per_event_data{ii};
        for jj = 1:length(traces)
            traces_stacked{jj} = [traces_stacked{jj}; traces{jj}];
        end
    end

    % stack traces at each x receiver
    for ii = 1:length(all_xreceivers)
        xtrace = all_xreceivers(ii);
        tr = traces_stacked{ii};
        tb = array2table(tr, 'VariableNames', {'time', 'amp'});
        tb.amp = tb.amp - xtrace;
        tb = grpstats(tb, 'time', 'sum');
        tb = sortrows(tb, 'time', 'asc');
        new_tr = [tb.time, tb.sum_amp + xtrace];
        traces_stacked{ii} = new_tr;
    end

    tag = sprintf('%.4f-all', xshot);
    color_ = 'k';
    plot_per_shot_event(fig, tag, color_, traces_stacked, all_xreceivers);
end

function plot_per_shot_event(fig, tag, color_, traces, xreceivers)
% plot seismic traces for each shot's specific event

    global my_xscale my_xclip;
    figure(fig);
    hold on;
    % plot for per event
    for ii = 1:length(traces)
        xtrace = xreceivers(ii);
        tr = traces{ii};
        [time, amp] = deal(tr(:, 1), tr(:, 2));
        [x, y] = deal(amp, time);

        % scale and clip
        amp = x - xtrace;
        amp = amp * my_xscale;
        if my_xclip > 0
            idx = abs(amp) > my_xclip;
            amp(idx) = sign(amp(idx)) * my_xclip;
        end
        x = xtrace + amp;

        curve = plot(x, y, [color_, ':']);
        % alpha(curve, 0.3);
        curve.UserData.tag = tag;
        curve.UserData.xtrace = xtrace;
    end
    hold off;
end

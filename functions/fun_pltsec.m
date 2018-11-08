function fun_pltsyn(vred, xomit, nskip, ipol, tol, nptsw, itrev, idump, iamp, twin, imeth, iroute)
% pltsec routine for PLTSYN
% plot synthetic seismic sections using the time, amplitude,
% and phase of each arrival read in from unit 11

    global path_pltsyn_par path_pltsyn_com;
    global path_txout;
    global fid_11 fid_12 fid_18 fid_20;
    run(path_pltsyn_par);
    run(path_pltsyn_com);

    % declare variables
    figs = {};
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
    % tol = tol / xscale;
    % tol = 1e-5;
    npts = round((tmax - tmin) * sps) + 1;

    % NEW
    % sample rate(s), 0.001 -- 1ms
    sample_rate = 0.003;
    npts = floor((tmax - tmin) / sample_rate) + 1;

    nptsp = npts + 121;
    nomit = 0;
    if vred ~= 0, rvred = 1 / vred;
    else, rvred = 0; end
    nwin = round(twin * sps);

    % 40
    nomit = find(xomit(1:pseis) < -99998, 1);
    if isempty(nomit), nomit = pseis;
    else, nomit = nomit - 1; end
    % for ii = 1:pseis
    %     if xomit(ii) < -99998, break; end
    %     nomit = nomit + 1;
    % end

    % 50
    % 80
    time(1:npts) = tmin:sample_rate:tmax;
    % for ii = 1:npts
    %     if itrev ~= 1
    %         time(ii) = (ii - 1) / sps / tscale + orig;
    %     else
    %         time(ii) = (tmin + (ii-1)/sps - tmax) / tscale + orig;
    %     end
    % end

    if ishade == 1
        ptmm = abs(tscale) * sps;
        if ptmm < dens
            nvaplt = round(dens * npts / ptmm);
        else
            nvaplt = npts;
        end
        if nvaplt > ppseis
            fprintf('***  interpolated trace for shading too long  ***\n');
            return;
        end
        if nvaplt > npts
            % fun_intrtr(time, npts, nvaplt);
            time(1:nvaplt) = interp1(1:npts, time(1:npts), 1:nvaplt);
        end
    end

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

    if itrev ~= 1
        yorig = orig;
        iside = 1;
    else
        yorig = orig + tmm;
        iside = -1;
    end

    % 1000: for each xshot
    tline = fgetl(fid_11);
    while ischar(tline)
        tmp = sscanf(tline, '%10f%10d');
        assert(length(tmp) == 2 && abs(tmp(2)) == 1, '*** unexpected error in reading sect.out ***');
        if isempty(tmp)
            fprintf('***  insufficient data in sect.out  ***\n');
            return;
        end
        [xshot, idum] = deal(tmp(1), tmp(2));
        if iplot >= -1
            fig_current = figure('numbertitle', 'off', 'name', sprintf('xshot=%.3f km', xshot));
            figs{end+1} = fig_current;
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
                if iscale == 2
                    range_ = abs(dist_(nseis) - xshot);
                    if range_ == 0, range_ = 0.001; end
                    rr = range_ ^ rcor;
                end

                % 190: for each arrival
                for jj = 1:na(nseis)
                    tline = fgetl(fid_11);
                    tmp = sscanf(tline, '%12f');
                    assert(length(tmp) == 3, '***  insufficient data in sect.out  ***\n');
                    sect(nseis, jj, 1:3) = tmp;
                    if iscale == 2
                        samp = abs(sect(nseis, jj, 2) * rr);
                        if samp > maxamp, maxamp = samp; end
                    end
                    if vred > 0
                        sect(nseis, jj, 1) = sect(nseis, jj, 1) - abs(dist_(nseis) - xshot) / vred;
                    end
                end
            end
            % go to 180
            tline = fgetl(fid_11);
        end

        % 99
        if iscale == 2
            scalef = amp / maxamp;
            fprintf('***  scalef=%15.5f  ***\n', scalef * xnorm * rcor);
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

        traces = {};
        % 70
        for ii = 1:nseis
            npts = nhpts;
            if mod(ii, nskip) ~= 0, continue; end
            if dist_(ii) < xmin || dist_(ii) > xmax, continue; end
            if nomit > 0
                % 75
                if any(xomit(1:nomit) - dist_(ii) < 0.001), continue; end
            end
            nsexsp = nsexsp + 1;
            % 170
            r(1:nptsp) = 0;
            if na(ii) > 0
                % 270
                for jj = 1:na(ii)
                    if iconv == 1 && (sect(ii,jj,1)+twave < tmin || sect(ii,jj,1) > tmax)
                        continue;
                    end
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

            % 109
            for jj = 1:namp
                if abs(xamp(jj) - dist_(ii)) < 0.001
                    n0 = round((tamp(jj) - tmin) * sps) + 1;
                    n1 = n0 + nwin;
                    if n0 < 1, n0 = 1; end
                    if n1 < 1, n1 = 1; end
                    if n0 > npts, n0 = npts; end
                    if n1 > npts, n1 = npts; end
                    nwina = n1 - n0 + 1;
                    minamp = 0;
                    maxamp = 0;

                    if n1 > n0
                        if imeth ~= 1
                            % 360
                            minamp = min(s(n0:n1));
                            maxamp = max(s(n0:n1));
                            maxamp = (abs(minamp) + abs(maxamp)) / 2;
                        else
                            maxamp = sum(abs(s(n0:n1))) / nwina;
                        end
                        if maxamp > 0, maxamp = log10(maxamp); end
                        fprintf(fid_18, '%10.3f%10.3f%10.3f%10d\n', xamp(jj), maxamp, 0, ipamp(jj));
                    end
                end
            end

            if iplot < -1, continue; end
            % if iscale == 0
            %     % 480
            %     max_ = max(abs(s(1:npts)));
            %     if max_ > 0
            %         f = amp / max_;
            %         % 490
            %         s(1:npts) = f * s(1:npts);
            %     end
            % end
            % if iscale > 0
            %     range_ = abs(dist_(ii) - xshot);
            %     if range_ == 0, range_ = 0.001; end
            % end
            % if iscale == 1
            %     srn = scalef * (range_ / xnorm) ^ rcor;
            %     % 580
            %     s(1:npts) = srn * s(1:npts);
            % end
            % if iscale == 2
            %     sr = scalef * range_ ^ rcor;
            %     % 610
            %     s(1:npts) = sr * s(1:npts);
            % end
            % if clip > 0 && iscale == 1
            %     % 590
            %     idx = abs(s(1:npts)) > clip;
            %     s(idx) = sign(s(idx)) * clip;
            % end
            % tmean = (dist_(ii) - xmin) / xscale + orig;
            % 600
            % s(1:npts) = s(1:npts) / xscale + tmean;
            if ishade == 1
                if nvaplt > npts
                    s(1:nvaplt) = interp1(1:npts, s(1:npts), 1:nvaplt);
                end
                npts = nvaplt;
            end
            if tol > 0 && npts > 2
                nplt = 1;
                sp(nplt) = s(1);
                tp(nplt) = time(1);

                % 700
                % idxp = abs(s(2:npts-1) - tmean) >= tol;
                idxp = find(abs(s(2:npts-1)) >= tol);
                idxp = unique([idxp-1, idxp, idxp+1]);
                countp = length(idxp);
                if countp > 0
                    nplt = nplt + countp;
                    sp(2:2+countp-1) = s(idxp);
                    tp(2:2+countp-1) = time(idxp);
                end
                % for jj = 2:npts-1
                %     if abs(s(jj) - tmean) >= tol
                %         nplt = nplt + 1;
                %         sp(nplt) = s(jj);
                %         tp(nplt) = time(jj);
                %     end
                % end

                nplt = nplt + 1;
                sp(nplt) = s(npts);
                tp(nplt) = time(npts);
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
            idir = -idir;
        end
        plot_seis(traces, fig_current);
        % return;

        if iplot >= -1 && nsexsp > 0 && tol > 0
            per = ntplt / (nsexsp * npts) * 100;
            % 105
            fprintf('***  %6.2f%% of points plotted using spmin=%9.3e  ***\n', per, (tol*xscale));
        end

        % 665
        fprintf('>>>  shot at %10.3f km completed  <<<\n', xshot);
        if xshotn > -999998
            xshot = xshotn;
        end

        % go to 1000
        % tline = fgetl(fid_11);
    end

    if namp > 0
        fprintf(fid_18, '%10.3f%10.3f%10.3f%10d\n', 0, 0, 0, -1);
    end
end


function plot_seis(traces, fig)
% plot seismic traces

    figure(fig);
    hold on;
    for ii = 1:length(traces)
        tr = traces{ii};
        [time, amp] = deal(tr(:, 1), tr(:, 2));
        [x, y] = deal(amp, time);
        % x = x(1) + (x - x(1)) * 3;
        plot(x, y, 'k-');
    end
    hold off;
end

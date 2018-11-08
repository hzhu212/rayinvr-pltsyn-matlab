function [x] = fun_nmo(x, npts, sps, tmin, dist_, vrms)
% apply constant vrms velocity nmo to trace x

    global path_pltsyn_par;
    run(path_pltsyn_par);

    y = zeros(1, ppseis);

    term = (2 * dist_ / vrms) ^ 2;
    for ii = 1:npts
        t0 = tmin + (ii - 1) / sps;
        t = sqrt(t0 ^ 2 + term);
        n1 = ifix((t - tmin) * sps + 1);
        n2 = n1 + 1;
        if(n1 < 1 || n2 > npts)
            y(ii) = 0;
        else
            t1 = tmin + (n1 - 1) / sps;
            y(ii) = (x(n2)-x(n1)) * (t-t1) * sps + x(n1);
        end
    end

    x(1:npts) = y(1:npts);
end

function main(options)
% pltsyn main function

    % globals
    global fid_11 fid_12 fid_17 fid_18 fid_19 fid_20;
    global path_txout;
    global path_pltsyn_par path_pltsyn_com;

    % option variables
    path_in = options.path_in;
    path_sin_m = options.path_sin_m;
    path_win = fullfile(path_in, 'w.in');
    path_synout = fullfile(path_in, 'syn.out');
    path_sectout = fullfile(path_in, 'sect.out');
    path_txout = fullfile(path_in, 'tx.out');
    path_ampout = fullfile(path_in, 'amp.out');
    path_pout = fullfile(path_in, 'p.out');

    % init .par and .com parameters
    path_pltsyn_par = 'pltsyn_par.m';
    path_pltsyn_com = 'pltsyn_com.m';
    run(path_pltsyn_par);
    run(path_pltsyn_com);

    % init variables
    xomit = ones(1, pseis) * -99999;

    % default parameter values
    iroute = 1;
    twin = 0.25;
    imeth = 1;
    iamp = 0;
    itrev = 0;
    idump = 0;
    nptsw = 19;
    spmin = 1e-5;
    ipol = 0;
    nskip = 1;
    vred = 8;

    % read s.in variables
    run(path_sin_m);

    % IO units
    fid_11 = fopen(path_sectout, 'r');
    if idump == 1, fid_12 = fopen(path_synout, 'w'); end
    if iwavlt == 2, fid_20 = fopen(path_win, 'r'); end
    if iplot == -1, iplot=-2;
    elseif(iplot == 0), iplot=-1;
    elseif(iplot == 2), iplot=0;
    else, ; end
    if iplot == -1 || iplot == 0, fid_19 = fopen(path_pout, 'w'); end

    xscale = (xmax - xmin) / xmm;
    tscale = (tmax - tmin) / tmm;
    if itrev == 1, tscale = -tscale; end

    ipol = -2 * ipol + 1;
    if amp < 0., amp = (xmax - xmin) / 100; end
    if spmin < 0., spmin = (xmax - xmin) / 100000; end

    if iamp > 0
        fid_17 = fopen(path_txout, 'r');
        fid_18 = fopen(path_ampout, 'w');
    end

    if iroute ~= 1, ibcol = 0; end

    fun_pltsec(vred,xomit,nskip,ipol,spmin,nptsw,itrev,idump,iamp,twin,imeth,iroute);

    % if iplots == 1, fun_plotnd(1); end

    fclose('all');
end

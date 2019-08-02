function main(options)
% pltsyn main function

    % globals
    global fid_11 fid_12 fid_17 fid_18 fid_19 fid_20;
    global path_txout;
    global path_pltsyn_par path_pltsyn_com;

    % option variables
    working_dir = options.working_dir;
    path_sin_m = options.path_sin_m;
    path_win = fullfile(working_dir, 'w.in');
    path_synout = fullfile(working_dir, 'syn.out');
    path_sectout = fullfile(working_dir, 'sect_ext.out');
    path_txout = fullfile(working_dir, 'tx.out');
    path_ampout = fullfile(working_dir, 'amp.out');
    path_pout = fullfile(working_dir, 'p.out');

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
    spmin = 1e-6;
    ipol = 0;
    nskip = 1;
    vred = 8;

    % default values for custom parameters (not defined by Zelt rayinvr).
    % all prefixed by 'my_'.
    % these are only default values, please set them in `s.in`.
    global my_xscale my_xclip;
    % scale rate for x axis
    my_xscale = 1;
    % clip for amptitude. if <= 0, do not clip, otherwise, clip to `my_xclip`
    my_xclip = 0;

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

    % xscale = (xmax - xmin) / xmm;
    % tscale = (tmax - tmin) / tmm;
    % if itrev == 1, tscale = -tscale; end

    ipol = -2 * ipol + 1;
    if amp < 0., amp = (xmax - xmin) / 100; end
    if spmin < 0., spmin = (xmax - xmin) / 100000; end

    if iamp > 0
        fid_17 = fopen(path_txout, 'r');
        fid_18 = fopen(path_ampout, 'w');
    end

    % if iroute ~= 1, ibcol = 0; end

    [figs] = fun_pltsec(vred, xomit, nskip, ipol, spmin, nptsw, itrev, idump, iamp, twin, imeth, iroute);

    fclose('all');

    fun_store_pltsyn_plot(figs, working_dir);
    fun_store_segy(figs, working_dir);
end

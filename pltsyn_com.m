% common blocks for PLTSYN

% common /blk1/
global xmin xmax xtmin xtmax xmm ntickx ndecix xscale tmin tmax ttmin ttmax ...
    tmm ntickt ndecit tscale albht iplots orig;

% common /blk2/
global rcor xnorm scalef amp iscale clip nsmth iconv sps iwavlt;

% common /blk3/
global inmo vrms;

% common /blk4/
global ishade ifill dens;

% common /cplot/
global iplot isep iseg nseg xwndow ywndow ibcol ifcol sf;


%% blkdat.f
% Block data for PLTSYN

% assign default and initial values to common block parameters
global is_inited;

if isempty(is_inited)
    xmin = 0;
    xmax = 300;
    xmm = 250;
    ntickx = -1;
    ndecix = -2;

    tmin = 0;
    tmax = 10;
    tmm = 125;
    ntickt = -1;
    ndecit = -2;
    albht = 2.5;

    xtmin = -999999;
    xtmax = -999999;
    ttmin = -999999;
    ttmax = -999999;

    rcor = 1;
    xnorm = 100;
    scalef = 100;
    amp = -1;
    iscale = 0;
    clip = 0;

    nsmth = 0;
    iconv = 1;
    sps = 60;
    iwavlt = 0;

    ishade = 0;
    ifill = 1;
    dens = 5;
    inmo = 0;
    vrms = 2.5;

    iplot = 1;
    iplots = 0;
    orig = 12.5;
    iseg = 0;
    nseg = 0;

    xwndow = 0;
    ywndow = 0;
    ibcol = 0;
    ifcol = 1;
    sf = 1.2;

    is_inited = true;
end

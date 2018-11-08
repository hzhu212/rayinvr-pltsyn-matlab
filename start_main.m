function start_main(path_in)
% starting script for main.m

    if nargin < 1
        path_in = fullfile(pwd, 'data', 'examples', 'e3');
    end

    % add functions to path.
    % `genpath` will add all the sub-directories too.
    addpath(genpath('./functions'));

    path_rin = fullfile(path_in, 'r.in');
    path_sin = fullfile(path_in, 's.in');
    options.path_in = path_in;
    options.path_sin_m = fun_trans_rin2m(path_sin);

    % clear cached r_in.m, or the modification for r.in won't work
    clear(options.path_sin_m);

    main(options);
end

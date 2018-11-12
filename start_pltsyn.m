function start_pltsyn(working_dir)
% start_main(working_dir)
%
% Starting script for rayinvr pltsyn.
% @param working_dir: the directory contains *.in files.

    if nargin < 1
        fprintf('Argument "working_dir" is required.\n');
        return;
    end

    % add functions to path.
    % `genpath` will add all the sub-directories too.
    addpath(genpath('./functions'));

    path_rin = fullfile(working_dir, 'r.in');
    path_sin = fullfile(working_dir, 's.in');
    options.working_dir = working_dir;
    options.path_sin_m = fun_trans_rin2m(path_sin);

    % clear cached r_in.m, or the modification for r.in won't work
    clear(options.path_sin_m);

    main(options);
end

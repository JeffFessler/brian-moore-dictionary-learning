function [vars, opts] = road_dls_par2b()
% Syntax:   [vars, opts] = road_dls_par2b();

% Knobs
vars.inpath   = 'data/road396.mat';
vars.rawpath  = 'road_dls_par2b/data.mat';
vars.outpath  = 'road_dls_par2b.mat';
vars.idim     = [128, nan, 200];
vars.nIters   = 20;
vars.np       = 0;
vars.algo     = 'paint2';
vars.method   = 1;
vars.SNR      = inf;
% sweep
vars.p        = [0.5, 0.6, 0.7, 0.8, 0.9];
vars.seed     = [1];
vars.lambda   = [1e-4];
vars.mu       = [0.09];
vars.dr       = [nan]; % FIXED D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dls opts
opts.A        = 1;
%opts.M
%opts.xdim
opts.pdim     = [8, 8, 5];
opts.pgap     = [2, 2, 1];
opts.type     = 'hard';
%opts.dr
opts.ddim     = [prod(opts.pdim(1:(end - 1))), opts.pdim(end)];
opts.fixedD   = true; % FIXED D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%opts.nIters
opts.nItersDB = 1;
opts.nItersX  = -1;
opts.D0       = dctmtx(prod(opts.pdim))';
%opts.X0
%opts.Xtrue
opts.flag     = 1;

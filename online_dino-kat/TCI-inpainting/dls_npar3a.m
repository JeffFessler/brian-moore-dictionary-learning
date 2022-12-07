function [vars, opts, vars0] = dls_npar3a()
% Syntax:   [vars, opts, vars0] = dls_npar3a();

% Knobs
vars.inpath   = 'data/coastguard2.mat';
vars.rawpath  = 'dls_npar3a/data.mat';
vars.outpath  = 'dls_npar3a.mat';
vars.idim     = [nan, nan, 150]; % 89,150
vars.nIters   = 20;
vars.np       = 0;
vars.SNR      = 19; %% NOISY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sweep
vars.p        = [0.7];
vars.seed     = [1];
vars.lambda   = 0.01 * logspace(-0.7,0.7,5);
vars.mu       = 0.25 * logspace(-0.7,0.7,5);
vars.dr       = [1, 5];

% dls opts
opts.A        = 1;
%opts.M
%opts.xdim
opts.pdim     = [8, 8, 5];
opts.pgap     = [4, 4, 1];
opts.type     = 'hard';
%opts.dr
opts.ddim     = [prod(opts.pdim(1:(end - 1))), opts.pdim(end)];
opts.fixedD   = false;
%opts.nIters
opts.nItersDB = 1;
opts.nItersX  = -1;
opts.D0       = dctmtx(prod(opts.pdim))';
%opts.X0
%opts.Xtrue
opts.flag     = 1;

% interpVideo opts
vars0.algo     = 'grid2';
vars0.method   = 'cubic';

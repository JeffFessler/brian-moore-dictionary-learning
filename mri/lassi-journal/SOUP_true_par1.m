function [vars, opts] = SOUP_true_par1()
% Syntax:   [vars, opts] = SOUP_true_par1();

% Knobs
vars.inpath   = 'otazo_full.mat';
vars.rawpath  = 'SOUP_true_par1/data.mat';
vars.outpath  = 'SOUP_true_par1.mat';
vars.lambdaB  = 0.03 * logspace(-1,1,9);
vars.dr       = [1, 2, 3, 4, 5];
vars.nIters   = 50;
vars.pdim     = [8, 8, 5];
vars.pgap     = [2, 2, 2];

% SOUP-DIL opts
opts.type    = 'hard';
%opts.dr     = dr;
opts.ddim    = [prod(vars.pdim(1:(end - 1))), vars.pdim(end)];
opts.fixedD  = false;
%opts.nIters = vars.nIters;
opts.D0      = dctmtx(prod(vars.pdim));
%opts.B0     = zeros(size(opts.D0,2),size(Y,2));
opts.flag    = 1;

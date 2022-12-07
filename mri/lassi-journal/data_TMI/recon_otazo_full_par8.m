function [vars, opts, vars0, opts0] = recon_otazo_full_par8()
% Syntax:   [vars, opts, vars0, opts0] = recon_otazo_full_par8();

% Knobs
vars.inpath    = 'otazo_full.mat';
vars.outpath   = 'recon_otazo_full_par8.mat';
vars.ps        = 1 ./ 8; % [4, 8, 12, 16, 20, 24]
vars.seed      = 3;
vars.r         = 1;
vars.lambdaL   = nan;
vars.lambdaS   = 0.01;
vars.lambdaB   = 0.02;
vars.lambdaB0  = 0.02;
vars.dr        = 1;
vars.nIters    = 50;

% DRPCA opts
%opts.A       = A;
%opts.sdim    = [ny, nx, nt];
opts.pdim     = [8, 8, 5];
opts.pgap     = [2, 2, 2];
opts.type     = 'soft';
%opts.dr      = dr;
opts.ddim     = [prod(opts.pdim(1:(end - 1))), opts.pdim(end)];
%opts.nIters  = vars.nIters;
opts.nItersDB = 1;
opts.nItersLS = 5;
%opts.L0      = reshape(init.Lhat,[],nt);
%opts.S0      = reshape(init.Shat,[],nt);
opts.D0       = dctmtx(prod(opts.pdim));
%opts.Xtrue   = reshape(data.Xtrue,[],nt);
opts.accel    = false;
opts.tau      = 1;
opts.flag     = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RPCA initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Knobs
vars0.r        = nan;
vars0.lambdaL  = 1.1955; % visual
%vars0.lambdaL = 0.53;   % mmse
vars0.lambdaS  = 0.01;
vars0.nIters   = 250;

% RPCA opts
%opts0.A       = A;
%opts0.T       = T;
%opts0.r       = vars0.r;
%opts0.nIters  = vars0.nIters;
%opts0.L0      = reshape(data.Xfft,[],nt);
%opts0.S0      = zeros(ny * nx,nt);
%opts0.Xtrue   = reshape(data.Xtrue,[],nt);
opts0.accel    = false;
opts0.tau      = 1;
opts0.flag     = 1;
function [vars, opts, vars0, opts0] = LASSI_full_par1c()
% Syntax:   [vars, opts, vars0, opts0] = LASSI_full_par1c();

% Knobs
vars.inpath    = 'otazo_full.mat';
vars.rawpath   = 'LASSI_full_par1c/data.mat';
vars.outpath   = 'LASSI_full_par1c.mat';
vars.ps        = 1 ./ [4, 8, 12, 16, 20, 24];
vars.seed      = 1:5;
vars.r         = nan;
vars.lambdaL   = 0.5;
vars.lambdaS   = 0.02;
vars.lambdaB   = 0.03;
vars.lambdaB0  = 0.03;
vars.dr        = 1;
vars.nIters    = 50;

% DRPCA opts
%opts.A       = A;
%opts.sdim    = [ny, nx, nt];
opts.pdim     = [8, 8, 5];
opts.pgap     = [2, 2, 2];
opts.type     = 'hard';
%opts.dr      = dr;
opts.ddim     = [prod(opts.pdim(1:(end - 1))), opts.pdim(end)];
opts.fixedD   = false;                              % Optimize dictionary
%opts.nIters  = vars.nIters;
opts.nItersDB = 1;
opts.nItersLS = 5;
%opts.L0      = reshape(init.Lhat,[],nt);
%opts.S0      = reshape(init.Shat,[],nt);
opts.D0       = kron(kron(dctmtx(opts.pdim(3)), ... % Separable DCT
                          dctmtx(opts.pdim(2))), ...
                          dctmtx(opts.pdim(1)));
%opts.Xtrue   = reshape(data.Xtrue,[],nt);
opts.accel    = false;
opts.tau      = 1;
opts.flag     = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Knobs
vars0.r        = nan;
%vars0.lambdaL = 1.1955; % visual
vars0.lambdaL  = 0.53;   % mmse
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

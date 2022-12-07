function [vars, opts, vars0, opts0] = invivo_LASSI_full_Sonly_par1()
% Syntax:   [vars, opts, vars0, opts0] = invivo_LASSI_full_Sonly_par1();

% Knobs
vars.inpath   = 'invivo_full.mat';
vars.rawpath  = 'invivo_LASSI_full_Sonly_par1/data.mat';
vars.outpath  = 'invivo_LASSI_full_Sonly_par1.mat';
vars.nLines   = [5, 10, 15, 20, 25, 30];
vars.seed     = 1:5;
vars.r        = nan;
vars.lambdaL  = 1e100;
vars.lambdaS  = 0.001;
vars.lambdaB  = 0.01;
vars.lambdaB0 = 0.01;
vars.dr       = 5;
vars.nIters   = 50;

% DRPCA opts
%opts.A       = A;
%opts.sdim    = [ny, nx, nt];
opts.pdim     = [8, 8, 5];
opts.pgap     = [2, 2, 2];
opts.type     = 'hard';
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
opts.tau      = 2; % Since lambdaL = inf
opts.flag     = 1;

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RPCA initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Knobs
vars0.r       = nan;
vars0.lambdaL = 0.158;
vars0.lambdaS = 0.0025;
vars0.nIters  = 250;

% RPCA opts
%opts0.A      = A;
%opts0.T      = T;
%opts0.r      = vars0.r;
%opts0.nIters = vars0.nIters;
%opts0.L0     = reshape(Xfft,[],nt);
%opts0.S0     = zeros(ny * nx,nt);
opts0.accel   = false;
opts0.tau     = 1;
opts0.flag    = 1;
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ktSLR initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Knobs
vars0.p        = 1;
vars0.mu1      = 0;
vars0.mu2      = 0.001;
vars0.nItersO  = 5;
vars0.nItersI  = 10;
vars0.nItersCG = 5;

% ktSLR opts
%opts0.nItersO  = vars0.nItersO;
%opts0.nItersI  = vars0.nItersI;
%opts0.nItersCG = vars0.nItersCG;
opts0.beta1     = 1;
opts0.beta2     = 1;
opts0.beta1rate = 50;
opts0.beta2rate = 25;
opts0.stepSize  = [1, 1, 0.303];

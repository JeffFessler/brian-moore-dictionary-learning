function [vars, opts, vars0, opts0] = recon_otazo_full_par11()
% Syntax:   [vars, opts, vars0, opts0] = recon_otazo_full_par11();

% Knobs
vars.inpath    = 'otazo_full.mat';
vars.outpath   = 'recon_otazo_full_par11.mat';
vars.ps        = 1 ./ 20; % [4, 8, 12, 16, 20, 24]
vars.seed      = 3;

% NO LASSI
opts = struct();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RPCA initialization
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

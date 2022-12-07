function otazo_R8_lassi_Sonly_run(idx)
% Syntax: otazo_R8_lassi_Sonly_run();
%         otazo_R8_lassi_Sonly_run(idx);

% Get parameters
basename = mfilename();
basename = basename(1:(end - 4));
[vars, opts] = eval(sprintf('%s_par%d();',basename,idx));

% Add dependencies to path
addpath('./deps_lassi');

% Load undersampled data (Y, mask, samp, Xfft, Xtrue)
data = load(vars.inpath);
[ny, nx, nt] = size(data.Xtrue);
A = Emat_xyt(data.mask,data.samp,[ny, nx, nt]);

% Load L + S initialization (Lhat, Shat, lambdaL, lambdaS, nIters)
init = load(vars.initpath);

% Run LASSI
opts.A       = A;
opts.r       = vars.r;
opts.sdim    = [ny, nx, nt];
opts.dr      = vars.dr;
opts.nIters  = vars.nIters;
opts.L0      = zeros(ny * nx,nt); % 0
opts.S0      = reshape(init.Lhat + init.Shat,[],nt); % L + S
opts.Xtrue   = reshape(data.Xtrue,[],nt);
lB = logspace(log10(vars.lambdaB0),log10(vars.lambdaB),vars.nIters);
[Lhat, Shat, Dhat, ~, stats] = drpca(data.Y,vars.lambdaL,vars.lambdaS,lB,opts); %#ok

% Remove large params
opts.A     = nan;
opts.L0    = nan;
opts.S0    = nan;
opts.D0    = nan;
opts.Xtrue = nan;

% Save results
save(vars.outpath,'Lhat','Shat','Dhat','vars','opts','stats');
fprintf('DONE\n');

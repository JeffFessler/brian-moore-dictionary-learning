function recon_otazo_full_run(idx)
% Syntax: recon_otazo_full_run(idx);

% Get parameters
basename = mfilename();
basename = basename(1:(end - 4));
[vars, opts, vars0, opts0] = eval(sprintf('%s_par%d();',basename,idx));

% Determine initialization scheme
runRPCA  = isfield(vars0,'lambdaL');
runKTSLR = isfield(vars0,'mu1');
runLASSI = isfield(vars,'dr');

% Generate undersampled data
[Y, A, ~, Xtrue, Xfft] = generateCardiacPerfData(vars.ps,inf,vars.seed,vars.inpath);
[ny, nx, nt] = size(Xtrue);

% Initialization
if runRPCA
    % Add dependencies to path
    addpath('./deps_rpca');
    
    % Run RPCA
    opts0.A      = A;
    opts0.T      = TempFFT(3,[ny, nx, nt]);
    opts0.r      = vars0.r;
    opts0.nIters = vars0.nIters;
    opts0.L0     = reshape(Xfft,[],nt);
    opts0.S0     = zeros(ny * nx,nt);
    opts0.Xtrue  = reshape(Xtrue,[],nt);
    [L0, S0, stats0] = robustPCA(Y,vars0.lambdaL,vars0.lambdaS,opts0); %#ok
    if isfield(vars0,'L0SX') && (vars0.L0SX == true)
        % L = 0, S = X
        Linit = zeros(ny,nx,nt);
        Sinit = L0 + S0;
    else
        % L = L, S = S
        Linit = L0;
        Sinit = S0;
    end
elseif runKTSLR
    % Add dependencies to path
    addpath('./deps_ktslr');
    
    % Run k-t SLR
    opts0.nItersO  = vars0.nItersO;
    opts0.nItersI  = vars0.nItersI;
    opts0.nItersCG = vars0.nItersCG;
    [X0, stats0] = ktSLR(Y,A,vars0.p,vars0.mu1,vars0.mu2,Xfft,opts0); %#ok
    
    % Construct LASSI initialization
    Linit = zeros(ny,nx,nt);
    Sinit = X0;
else
    % Unrecognized algorithm
    error('Unrecognized algorithm');
end

% LASSI
if runLASSI
    % Add dependencies to path
    addpath('./deps_lassi');

    % Run LASSI
    opts.A       = A;
    opts.r       = vars.r;
    opts.sdim    = [ny, nx, nt];
    opts.dr      = vars.dr;
    opts.nIters  = vars.nIters;
    opts.L0      = reshape(Linit,[],nt);
    opts.S0      = reshape(Sinit,[],nt);
    opts.Xtrue   = reshape(Xtrue,[],nt);
    lB = logspace(log10(vars.lambdaB0),log10(vars.lambdaB),vars.nIters);
    [Lhat, Shat, Dhat, ~, stats] = drpca(Y,vars.lambdaL,vars.lambdaS,lB,opts);%#ok
    %[Lhat, Shat, Dhat, Bhat, stats] = drpca(Y,vars.lambdaL,vars.lambdaS,lB,opts);%#ok
end

% Save results
varsToSave = {'vars','vars0'};
if runRPCA
    % RPCA outputs
    varsToSave = [varsToSave, {'L0','S0','stats0'}];
end
if runKTSLR
    % KTSLR outputs
    varsToSave = [varsToSave, {'X0','stats0'}];
end
if runLASSI
    % LASSI outputs
    varsToSave = [varsToSave, {'Lhat','Shat','Dhat','stats'}];
    %varsToSave = [varsToSave, {'Lhat','Shat','Dhat','Bhat','stats'}];
end
save(vars.outpath,varsToSave{:});
fprintf('DONE\n');

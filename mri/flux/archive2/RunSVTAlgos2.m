function RunSVTAlgos2(idx)
% Syntax:   RunSVTAlgos2(idx);

%--------------------------------------------------------------------------
% Simulation setup
%--------------------------------------------------------------------------
% runpbs parameters
MAX_IDX = 250; % Max # indices
offset = 0; % data#.mat offset

% Data parameters
nc_list = [8]; % # coils
%nt_list = [2 4 6 8 10 12 16 20 24 30 32 40]; % # time points
%SNR_list = linspace(10,60,10); % SNR values
nt_list = [2 4 6 8 12 14 21 24 28 42]; % # time points
SNR_list = 40:5:60; % SNR values
seed_list = 1; % rng seeds

% Algo parameters
% L + S
algoStr = 'L+S PGM'; % {'L+S PGM','L&S ADMM'}
proxS = 'mixed21'; % {'soft','mixed21'}
T = 1; % {1,TempFFT(2)}
lambdaL_list = logspace(log10(0.001),log10(10),11);
lambdaS_list = logspace(log10(0.001),log10(10),11);
%{
% L & S
algoStr = 'L&S ADMM'; % {'L+S PGM','L&S ADMM'}
proxS = 'soft'; % {'soft','mixed21'}
T = 1; % {1,TempFFT(2)}
lambdaL_list = logspace(log10(0.001),log10(10),11);
lambdaS_list = logspace(log10(0.001),log10(10),11);
%}

% Paths
irtpath = '/home/brimoor/Research/irt/'; % IRT toolbox location
proxpath = '/home/brimoor/Research/optshrink4mri/prox/'; % Proximal functions folder
%inpath = '/home/brimoor/Research/optshrink4mri/dynobj.mat'; % Ground truth data
%outpath = '/home/brimoor/Research/optshrink4mri/svt_LpS_soft/data.mat'; % Output data path
inpath = '/home/brimoor/Research/optshrink4mri/dynobj2.mat'; % Ground truth data
outpath = '/home/brimoor/Research/optshrink4mri/svt_LpS_mixed21_2/data.mat'; % Output data path
%--------------------------------------------------------------------------

% Get indices for this iteration
Nc = length(nc_list);
Nt = length(nt_list);
NSNR = length(SNR_list);
Nseed = length(seed_list);
NlambdaL = length(lambdaL_list);
NlambdaS = length(lambdaS_list);
N = Nc * Nt * NSNR * Nseed * NlambdaL * NlambdaS;
inds = getInds(idx,N,MAX_IDX);
Ninds = length(inds);
if (Ninds == 0)
    fprintf('Nothing to do for index %i/%i\n',idx,MAX_IDX);
    return;
end
[nc_inds nt_inds SNR_inds seed_inds lambdaL_inds lambdaS_inds] = ind2sub([Nc Nt NSNR Nseed NlambdaL NlambdaS],inds);

% Get parameters for this iteration
nc = nc_list(nc_inds(:));
nt = nt_list(nt_inds(:));
SNR = SNR_list(SNR_inds(:));
seed = seed_list(seed_inds(:));
lambdaL = lambdaL_list(lambdaL_inds(:));
lambdaS = lambdaS_list(lambdaS_inds(:));

% Set paths
irtsetup(irtpath); % Add irt toolbox to path
addpath(regexprep(proxpath,'\','/')); % Add prox/ directory to path

% Load ground truth data
gtData = load(inpath);
Xtrue_full = gtData.dyn_obj;
dce = gtData.dce;
clear gtData;

% Perform reconstructions
labels = {'nc','nt','SNR','seed','lambdaL','lambdaS','NRMSE','NRMSE_ROI1','NRMSE_ROI2','NRMSE_ROI3','NRMSE_ROI'}; %#ok
data = [nc(:) nt(:) SNR(:) seed(:) lambdaL(:) lambdaS(:) nan(Ninds,5)];
NRMSEfcn = @(Xhat,Xtrue) norm(Xhat(:) - Xtrue(:)) / norm(Xtrue(:));
for i = 1:Ninds
    % Start simulation timer
    stimer = tic;
    
    % Perform reconstruction
    [~,recon,~,~,mask,ROIs,~,splineInterp] = run_svt_algo(nc(i),nt(i),SNR(i),seed(i),lambdaL(i),proxS,lambdaS(i),algoStr,T,Xtrue_full,dce);
    
    % Interpolate reconstruction
    Xhat = splineInterp(abs(embed(recon.X,mask)));
    
    % Compute NRMSEs
    nf = size(Xtrue_full,3);
    data(i,end - 4) = NRMSEfcn(Xhat,Xtrue_full);
    M1 = repmat(ROIs{1},[1 1 nf]); % Lesion #1
    data(i,end - 3) = NRMSEfcn(Xhat(M1),Xtrue_full(M1));
    M2 = repmat(ROIs{2},[1 1 nf]); % Lesion #2
    data(i,end - 2) = NRMSEfcn(Xhat(M2),Xtrue_full(M2));
    M3 = repmat(ROIs{3},[1 1 nf]); % Lesion #3
    data(i,end - 1) = NRMSEfcn(Xhat(M3),Xtrue_full(M3));
    M = (M1 | M2 | M3); % All lesions
    data(i,end) = NRMSEfcn(Xhat(M),Xtrue_full(M));
    
    % Display progress
    fprintf('*** Simulation %i/%i complete [Time = %.2fs]\n',i,Ninds,toc(stimer));
end

% Make output directory, if necessary
[outdir outfile outext] = fileparts(outpath);
if ~exist(outdir,'dir')
    mkdir(outdir);
end

% Save data to file
wtimer = tic;
fileIdx = offset + idx;
outpathn = sprintf('%s/%s%i%s',outdir,outfile,fileIdx,outext); % Append #
save(outpathn,'labels','data');
fprintf('*** File "%s" written [Time = %.2fs]\n',outpathn,toc(wtimer));

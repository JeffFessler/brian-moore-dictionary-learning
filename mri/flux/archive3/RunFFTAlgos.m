function RunFFTAlgos(idx)
% Syntax:   RunFFTAlgos(idx);

%--------------------------------------------------------------------------
% Simulation setup
%--------------------------------------------------------------------------
% runpbs parameters
MAX_IDX = 250; % Max # indices
offset = 0; % data#.mat offset

% Data parameters
nc_list = [8]; % # coils
%nt_list = [2 4 6 8 10 12 16 20 24 30 32 40]; % #1
%SNR_list = linspace(10,60,10);               % #1
%nt_list = [2 4 6 8 12 14 21 24 28 42];       % #2
%SNR_list = 40:5:60;                          % #2
nt_list = [2 4 6 8 12 16 21 24 28 42];        % #3
SNR_list = 40:5:60;                           % #3
seed_list = 1:200; % rng seeds

% Paths
irtpath = '/home/brimoor/Research/irt/'; % IRT toolbox location
proxpath = '/home/brimoor/Research/optshrink4mri/prox/'; % Proximal functions folder
%inpath = '/home/brimoor/Research/optshrink4mri/dynobj.mat'; % Ground truth data
%outpath = '/home/brimoor/Research/optshrink4mri/fft/data.mat'; % Output data path
%inpath = '/home/brimoor/Research/optshrink4mri/dynobj2.mat'; % Ground truth data
%outpath = '/home/brimoor/Research/optshrink4mri/fft_2/data.mat'; % Output data path
inpath = '/home/brimoor/Research/optshrink4mri/dynobj3.mat'; % Ground truth data
outpath = '/home/brimoor/Research/optshrink4mri/fft_3/data.mat'; % Output data path
%--------------------------------------------------------------------------

% Get indices for this iteration
Nc = length(nc_list);
Nt = length(nt_list);
NSNR = length(SNR_list);
Nseed = length(seed_list);
N = Nc * Nt * NSNR * Nseed;
inds = getInds(idx,N,MAX_IDX);
Ninds = length(inds);
if (Ninds == 0)
    fprintf('Nothing to do for index %i/%i\n',idx,MAX_IDX);
    return;
end
[nc_inds nt_inds SNR_inds seed_inds] = ind2sub([Nc Nt NSNR Nseed],inds);

% Get parameters for this iteration
nc = nc_list(nc_inds(:));
nt = nt_list(nt_inds(:));
SNR = SNR_list(SNR_inds(:));
seed = seed_list(seed_inds(:));

% Set paths
irtsetup(irtpath); % Add irt toolbox to path
addpath(regexprep(proxpath,'\','/')); % Add prox/ directory to path

% Load ground truth data
gtData = load(inpath);
Xtrue_full = gtData.dyn_obj;
dce = gtData.dce;
clear gtData;

% Perform simulations
labels = {'nc','nt','SNR','seed','NRMSE','NRMSE_ROI1','NRMSE_ROI2','NRMSE_ROI3','NRMSE_ROI'}; %#ok
data = [nc(:) nt(:) SNR(:) seed(:) nan(Ninds,5)];
NRMSEfcn = @(Xhat,Xtrue) norm(Xhat(:) - Xtrue(:)) / norm(Xtrue(:)); % NRMSE function handle
for i = 1:Ninds
    % Start simulation timer
    stimer = tic;
    
    % Perform reconstruction
    [~,Xfft,~,ROIs,~,splineInterp] = run_fft_algo(nc(i),nt(i),SNR(i),seed(i),Xtrue_full,dce);
    
    % Interpolate reconstruction
    Xhat = splineInterp(abs(Xfft));
    
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

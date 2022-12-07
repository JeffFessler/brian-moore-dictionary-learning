function RunPCAAlgos(idx)
% Syntax:   RunPCAAlgos(idx);

%--------------------------------------------------------------------------
% Simulation setup
%--------------------------------------------------------------------------
% runpbs parameters
MAX_IDX = 250; % Max # indices
offset = 0; % data#.mat offset

%{
% Data parameters #1
p_list = [8 12 16 32 48 64 96 128] / 128;
SNR_list = 40:5:60;
seed_list = 1;
%}

% Data parameters #2
p_list = [4 8 12 16 20 24 28 32] / 128;
SNR_list = 10:10:60;
seed_list = 1;

% Algo parameters
% L + S
algoStr = 'L+S PGM'; % {'L+S PGM','L&S ADMM'}
proxS = 'soft'; % {'soft','mixed21'}
T = TempFFT(2); % {1,TempFFT(2)}
r_list = 1:10;
lambdaS_list = logspace(log10(0.001),log10(10),11);
%{
% L & S
algoStr = 'L&S ADMM'; % {'L+S PGM','L&S ADMM'}
proxS = 'mixed21'; % {'soft','mixed21'}
T = 1; % {1,TempFFT(2)}
r_list = 1:10;
lambdaS_list = logspace(log10(0.001),log10(10),11);
%}

% Paths
irtpath = '/home/brimoor/Research/irt/'; % IRT toolbox location
proxpath = '/home/brimoor/Research/otazo/prox/'; % Proximal functions folder
inpath = '/home/brimoor/Research/otazo/otazo_true.mat'; % Ground truth data
outpath = '/home/brimoor/Research/otazo/pca_LpS_soft_T/data.mat'; % Output data path
%--------------------------------------------------------------------------

% Get indices for this iteration
Np = length(p_list);
NSNR = length(SNR_list);
Nseed = length(seed_list);
Nr = length(r_list);
NlambdaS = length(lambdaS_list);
N = Np * NSNR * Nseed * Nr * NlambdaS;
inds = getInds(idx,N,MAX_IDX);
Ninds = length(inds);
if (Ninds == 0)
    fprintf('Nothing to do for index %i/%i\n',idx,MAX_IDX);
    return;
end
[p_inds SNR_inds seed_inds r_inds lambdaS_inds] = ind2sub([Np NSNR Nseed Nr NlambdaS],inds);

% Get parameters for this iteration
p = p_list(p_inds(:));
SNR = SNR_list(SNR_inds(:));
seed = seed_list(seed_inds(:));
r = r_list(r_inds(:));
lambdaS = lambdaS_list(lambdaS_inds(:));

% Set paths
irtsetup(irtpath); % Add irt toolbox to path
addpath(regexprep(proxpath,'\','/')); % Add prox/ directory to path

% Perform reconstructions
labels = {'p','SNR','seed','r','lambdaS','NRMSE','NRMSE_ROI1','NRMSE_ROI2'}; %#ok
data = [p(:) SNR(:) seed(:) r(:) lambdaS(:) nan(Ninds,3)];
for i = 1:Ninds
    % Start simulation timer
    stimer = tic;
    
    % Perform reconstruction
    [~,recon,Xtrue,ROIs] = run_pca_algo(p(i),SNR(i),seed(i),inpath,r(i),proxS,lambdaS(i),algoStr,T);
    [ny nx] = size(ROIs{1});
    nt = size(Xtrue,2);
    Xhat = reshape(recon.X,[ny nx nt]);
    Xtrue = reshape(Xtrue,[ny nx nt]);
    
    % Compute NRMSEs
    data(i,end - 2) = computeNRMSE(Xhat,Xtrue,0,ROIs); % Whole image
    data(i,end - 1) = computeNRMSE(Xhat,Xtrue,1,ROIs); % Heart-only
    data(i,end - 0) = computeNRMSE(Xhat,Xtrue,2,ROIs); % Body-only
    
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

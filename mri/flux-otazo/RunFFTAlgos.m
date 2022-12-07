function RunFFTAlgos(idx)
% Syntax:   RunFFTAlgos(idx);

%--------------------------------------------------------------------------
% Simulation setup
%--------------------------------------------------------------------------
% runpbs parameters
MAX_IDX = 250; % Max # indices
offset = 0; % data#.mat offset

% Data parameters #2
p_list = [4 8 12 16 20 24 28 32] / 128;
SNR_list = 10:10:60;
seed_list = 1:100;

% Paths
irtpath = '/home/brimoor/Research/irt/'; % IRT toolbox location
proxpath = '/home/brimoor/Research/otazo/prox/'; % Proximal functions folder
inpath = '/home/brimoor/Research/otazo/otazo_true.mat'; % Ground truth data
outpath = '/home/brimoor/Research/otazo/fft/data.mat'; % Output data path
%--------------------------------------------------------------------------

% Get indices for this iteration
Np = length(p_list);
NSNR = length(SNR_list);
Nseed = length(seed_list);
N = Np * NSNR * Nseed;
inds = getInds(idx,N,MAX_IDX);
Ninds = length(inds);
if (Ninds == 0)
    fprintf('Nothing to do for index %i/%i\n',idx,MAX_IDX);
    return;
end
[p_inds SNR_inds seed_inds] = ind2sub([Np NSNR Nseed],inds);

% Get parameters for this iteration
p = p_list(p_inds(:));
SNR = SNR_list(SNR_inds(:));
seed = seed_list(seed_inds(:));

% Set paths
irtsetup(irtpath); % Add irt toolbox to path
addpath(regexprep(proxpath,'\','/')); % Add prox/ directory to path

% Perform reconstructions
labels = {'p','SNR','seed','NRMSE','NRMSE_ROI1','NRMSE_ROI2'}; %#ok
data = [p(:) SNR(:) seed(:) nan(Ninds,3)];
for i = 1:Ninds
    % Start simulation timer
    stimer = tic;
    
    % Perform reconstruction
    [~,recon,Xtrue,ROIs] = run_fft_algo(p(i),SNR(i),seed(i),inpath);
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

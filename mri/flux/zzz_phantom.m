%% Generate dynobj

% Knobs
dceFcn = @ir_mri_dce_obj3;
%inpath = 'C:/Users/brimoor/Documents/Google Drive/MATLAB/Misc/irt/mri/phantom/';
inpath = '/Users/Brian/Google Drive/MATLAB/Misc/irt/mri/phantom/';
outpath = 'dynobj3.mat';

% Generate dynamic object
n_tr_merge = 100;
n_tr_round = 24;
[dyn_obj dce] = dceFcn('chat',0, ...
                       'brainweb_dir',inpath, ...
                       'n_tr_merge',n_tr_merge, ...
                       'n_tr_round',n_tr_round, ...
                       'phase','');

%{
% Playback movie
movie.video = permute(dyn_obj,[2 1 3]);
movie.Fv = 60;
opts.mag = 3;
PlayMovie(movie,opts);
%}

% Save object
save(outpath,'dyn_obj','dce');

%% Compute norms

% Knobs
inpath = 'dynobj3.mat';
nc = 8; % # coils

% Load data
load(inpath);

% Compute norms
%ntl = [2 4 6 8 10 12 16 20 24 30 32 40]; % #1
%ntl = [2 4 6 8 12 14 21 24 28 42]; % #2
ntl = [2 4 6 8 12 16 21 24 28 42]; % #3
normAl = zeros(size(ntl));
for i = 1:length(ntl)
    % Generate phantom
    nt = ntl(i);
    SNR = 40; % Irrelevant
    seed = 1; % Irrelevant
    [~,~,~,~,~,~,~,A] = GenerateFesslerPhantom(nc,nt,SNR,seed,dyn_obj,dce);
    
    % Compute norm
    AAtfcn = @(x) A' * (A * x);
    opts = struct('issym',true,'isreal',false,'disp',1);
    %opts = struct('issym',true,'isreal',false,'maxit',5,'disp',1);
    normAli = sqrt(eigs(AAtfcn,size(A,2),1,'LM',opts)) %#ok
    normAl(i) = normAli;
end

%% Compute scaling constants

% Knobs
inpath = 'dynobj3.mat';
outpath = 'constants3.mat';
%nt = [2 4 6 8 10 12 16 20 24 30 32 40]; % #1
%nt = [2 4 6 8 12 14 21 24 28 42]; % #2
nt = [2 4 6 8 12 16 21 24 28 42]; % #3

% Load dynobj
load(inpath);
ny = size(dyn_obj,1);
nx = size(dyn_obj,2);

% Compute parameters
C = struct('CL',nan,'soft',nan,'mixed21',nan,'n',nan);
constants = repmat(C,max(nt),1);
m = length(nt);
for i = 1:m
    % Form Xtrue_mat
    Xtrue = reshape(dyn_obj,[ny nx (size(dyn_obj,3) / nt(i)) nt(i)]);
    Xtrue = squeeze(mean(Xtrue,3));
    mask = (conv2(dce.labels,ones(3),'same') > 0);
    Xtrue_mat = masker(Xtrue,mask);
    
    % Compute parameters
    Ci = struct();
    CL = svd(Xtrue_mat);
    Ci.CL = CL(2);
    Ci.soft = max(abs(Xtrue_mat(:)));
    Ci.mixed21 = max(sqrt(sum(abs(Xtrue_mat).^2,2)));
    Ci.n = max(size(Xtrue_mat));
    constants(nt(i)) = Ci;
end

% Save values
save(outpath,'constants');

%% Plot contrast curves

% Knobs
inpath = 'dynobj3.mat';
nc = 8; % Irrelevant
nt = 24; % Irrelevant
SNR = 40; % Irrelevant
seed = 1; % Irrelevant
cm = linspecer(3);

% Load data
load(inpath);
[Xtrue Xfft mask ROIs nd splineInterp Y A] = GenerateFesslerPhantom(nc,nt,SNR,seed,dyn_obj,dce);
Xtrue_full = permute(dyn_obj,[2 1 3]);
M1 = permute(ROIs{1},[2 1 3]);
M2 = permute(ROIs{2},[2 1 3]);
M3 = permute(ROIs{3},[2 1 3]);
M = M1 | M2 | M3;
[ny nx nt] = size(Xtrue_full);

% Scale to [0 1]
Xtrue_full = Xtrue_full / max(Xtrue_full(:));

% Mean image
Xbar = mean(Xtrue_full,3);
XX = repmat(im2uint8(Xbar),[1 1 3]);

% Get static vals
X1 = Xtrue_full(:,:,1);
%{
imshow(X1);
[x y] = ginput(); delete(gcf);
inds = sub2ind(size(X1),round(y),round(x));
%}
inds = [14118; 10641; 21408];
static_vals = X1(inds(:)) %#ok
%static_vals = [0.17029; 0.11321; 0.037934]; % Otazo
%static_vals = [0.25687; 0.15830; 0.052760]; % Phantom #1
%static_vals = [0.17199; 0.10599; 0.035326]; % Phantom #2
%static_vals = [0.1440 0.0876 0.0288]; % Phantom #3

% Colored mask image
CM = uint8(floor(255 * cm));
%MM = 255 * ones(ny,nx,3,'uint8'); % White background
%MM = zeros(ny,nx,3,'uint8'); % Black background
MM = XX; % Mean background
for i = 1:3
    Mi = eval(sprintf('M%i;',i));
    for j = 1:3
        MMj = MM(:,:,j);
        MMj(Mi) = CM(i,j);
        MM(:,:,j) = MMj;
    end
end

% Constrast curves
Tbar = zeros(3,nt);
Tmin = zeros(3,nt);
Tmax = zeros(3,nt);
for i = 1:3
    Mi = eval(sprintf('M%i;',i));
    for j = 1:nt
        Xj = Xtrue_full(:,:,j);
        Xij = Xj(Mi);
        Tbar(i,j) = mean(Xij(:));
        Tmin(i,j) = min(Xij(:));
        Tmax(i,j) = max(Xij(:));
    end
end

%{
% Smooth max curves
w = 5; % Filter width
h = ones(1,w) / w;
for i = 1:5
    Tmin(i,:) = imfilter(Tmin(i,:),h,'replicate','same');
    Tmax(i,:) = imfilter(Tmax(i,:),h,'replicate','same');
end
%}

%{
%--------------------------------------------------------------------------
% Plot results #1
%--------------------------------------------------------------------------
figure('Color','w');

% Mean image
subplot(2,2,1);
imshow(XX);

% Colored mask image
subplot(2,2,2);
imshow(MM);

% Mean trajectories
subplot(2,2,3);
for i = 1:3
    plot(1:nt,Tbar(i,:),'-','Color',cm(i,:));
    hold on;
end
for j = 1:length(static_vals)
    plot(1:nt,static_vals(j) * ones(1,nt),'k:');
end
xlabel('Frame');
title('Mean Trajectories');
set(gca,'XLim',[1 nt]);
set(gca,'YLim',[0 max(Xtrue_full(:))]);
box on;

% Max trajectories
subplot(2,2,4);
ph = zeros(1,3);
for i = 1:3
    ph(i) = plot(1:nt,Tmax(i,:),'-','Color',cm(i,:));
    hold on;
end
for j = 1:length(static_vals)
    plot(1:nt,static_vals(j) * ones(1,nt),'k:');
end
xlabel('Frame');
title('Max Trajectories');
set(gca,'XLim',[1 nt]);
set(gca,'YLim',[0 max(Xtrue_full(:))]);
box on;

% Save figure
set(gcf,'Position',[579 343 764 517]);
%export_fig -pdf -transparent phantom_traj1
%--------------------------------------------------------------------------
%}

%--------------------------------------------------------------------------
% Plot results #2
%--------------------------------------------------------------------------
figure('Color','w');

% Images
subplot(3,3,1:6);
imshow(cell2mov({XX MM},round(size(XX,2) / 2),255));

% Trajectories
Nbars = 5; inds = round(linspace(2,nt - 1,Nbars));
for i = 1:3
    subplot(3,3,6 + i);
    hold on;
    xvec = 1:nt;
    plot(xvec,Tbar(i,:),'-','Color',cm(i,:));
    plot(xvec,Tmax(i,:),':','Color',cm(i,:));
    plot(xvec,Tmin(i,:),':','Color',cm(i,:));
    Li = Tbar(i,inds) - Tmin(i,inds);
    Ui = Tmax(i,inds) - Tbar(i,inds);
    errorbar(xvec(inds),Tbar(i,inds),Li,Ui,'.','Color',cm(i,:));
    for j = 1:length(static_vals)
        plot(1:nt,static_vals(j) * ones(1,nt),'k:');
    end
    xlabel('Frame');
    title(sprintf('Region %i',i));
    set(gca,'XLim',[1 nt]);
    set(gca,'YLim',[0 max(Xtrue_full(:))]);
    box on;
end

% Save figure
set(gcf,'Position',[579 343 764 517]);
%export_fig -pdf -transparent phantom_traj1
%export_fig -pdf -transparent phantom_traj2
%export_fig -pdf -transparent phantom_traj3
%--------------------------------------------------------------------------

%% Plot masks

% Knobs
maskpath = './masks.mat'; % Mask path
dim = [1000 500]; % Figure dimensions

% Load masks
load(maskpath);

% Create figure
scrsz = get(0,'ScreenSize');
xyc = 0.5 * scrsz(3:4);
pos = [(xyc - 0.5 * dim) dim];
figure('Position',pos);

% Display masks
[ni nj] = bestSubplotShape(5);
subplot(ni,nj,1); imshow(mask,[]); title('X');
subplot(ni,nj,2); imshow(ROI,[]); title('ROI');
subplot(ni,nj,3); imshow(ROI1,[]); title('ROI1');
subplot(ni,nj,4); imshow(ROI2,[]); title('ROI2');
subplot(ni,nj,5); imshow(ROI3,[]); title('ROI3');
drawnow;

% Save figure
%export_fig -pdf -transparent masks

%% Generate NRMSE curves

%{
% Knobs
fftpath = '/home/brimoor/Research/optshrink4mri/fft_3/data.mat';
svtpath = '/home/brimoor/Research/optshrink4mri/svt_LpS_soft_3/data.mat';
optpath = '/home/brimoor/Research/optshrink4mri/opt_LpS_soft_3/data.mat';
outpath = '/home/brimoor/Research/optshrink4mri/nrmse3/NRMSE_LpS_soft_3.mat';
%outpath = '/home/brimoor/Research/optshrink4mri/nrmse3/NRMSE_LpS_soft_full_3.mat';
%}

%{
% Knobs
fftpath = '/home/brimoor/Research/optshrink4mri/fft_3/data.mat';
svtpath = '/home/brimoor/Research/optshrink4mri/svt_LpS_soft_T_3/data.mat';
optpath = '/home/brimoor/Research/optshrink4mri/opt_LpS_soft_T_3/data.mat';
outpath = '/home/brimoor/Research/optshrink4mri/nrmse3/NRMSE_LpS_soft_T_3.mat';
%outpath = '/home/brimoor/Research/optshrink4mri/nrmse3/NRMSE_LpS_soft_full_T_3.mat';
%}

%{
% Knobs
fftpath = '/home/brimoor/Research/optshrink4mri/fft_3/data.mat';
svtpath = '/home/brimoor/Research/optshrink4mri/svt_LpS_mixed21_3/data.mat';
optpath = '/home/brimoor/Research/optshrink4mri/opt_LpS_mixed21_3/data.mat';
outpath = '/home/brimoor/Research/optshrink4mri/nrmse3/NRMSE_LpS_mixed21_3.mat';
%outpath = '/home/brimoor/Research/optshrink4mri/nrmse3/NRMSE_LpS_mixed21_full_3.mat';
%}

%{
% Knobs
fftpath = '/home/brimoor/Research/optshrink4mri/fft_3/data.mat';
svtpath = '/home/brimoor/Research/optshrink4mri/svt_LaS_soft_3/data.mat';
optpath = '/home/brimoor/Research/optshrink4mri/opt_LaS_soft_3/data.mat';
outpath = '/home/brimoor/Research/optshrink4mri/nrmse3/NRMSE_LaS_soft_3.mat';
%outpath = '/home/brimoor/Research/optshrink4mri/nrmse3/NRMSE_LaS_soft_full_3.mat';
%}

%{
% Knobs
fftpath = '/home/brimoor/Research/optshrink4mri/fft_3/data.mat';
svtpath = '/home/brimoor/Research/optshrink4mri/svt_LaS_soft_T_3/data.mat';
optpath = '/home/brimoor/Research/optshrink4mri/opt_LaS_soft_T_3/data.mat';
outpath = '/home/brimoor/Research/optshrink4mri/nrmse3/NRMSE_LaS_soft_T_3.mat';
%outpath = '/home/brimoor/Research/optshrink4mri/nrmse3/NRMSE_LaS_soft_full_T_3.mat';
%}

%{
% Knobs
fftpath = '/home/brimoor/Research/optshrink4mri/fft_3/data.mat';
svtpath = '/home/brimoor/Research/optshrink4mri/svt_LaS_mixed21_3/data.mat';
optpath = '/home/brimoor/Research/optshrink4mri/opt_LaS_mixed21_3/data.mat';
outpath = '/home/brimoor/Research/optshrink4mri/nrmse3/NRMSE_LaS_mixed21_3.mat';
%outpath = '/home/brimoor/Research/optshrink4mri/nrmse3/NRMSE_LaS_mixed21_full_3.mat';
%}

% Generate NRMSE curves
GenerateNRMSEcurves(fftpath,svtpath,optpath,outpath);

%% Plot fixed nc heatmaps

% Knobs
nc = 8; % # coils
types = {'X','ROI'}; % Types of heatmaps to plot
SNRs = [40 60]; % SNR range
NRMSEpath = './nrmse3/NRMSE_LpS_soft_3.mat';
%NRMSEpath = './nrmse3/NRMSE_LpS_mixed21_3.mat';
%NRMSEpath = './nrmse3/NRMSE_LaS_soft_3.mat';
%NRMSEpath = './nrmse3/NRMSE_LaS_mixed21_3.mat';
%dim = [1500 300]; % Figure dimensions
dim = [1000 200]; % Figure dimensions

% Load NRMSE data
load(NRMSEpath);

% Compute figure position
scrsz = get(0,'ScreenSize');
xyc = 0.5 * scrsz(3:4);
pos = [(xyc - 0.5 * dim) dim];

% Locate desired nc
[~,idx_nc_fft] = min(abs(labels_fft.nc - nc));
[~,idx_nc_svt] = min(abs(labels_svt.nc - nc));
[~,idx_nc_opt] = min(abs(labels_opt.nc - nc));

% Locate valid SNRs
idx_SNR1_fft = find(labels_fft.SNR >= SNRs(1),1,'first');
idx_SNR2_fft = find(labels_fft.SNR <= SNRs(2),1,'last');
inds_SNR_fft = idx_SNR1_fft:idx_SNR2_fft;
idx_SNR1_svt = find(labels_svt.SNR >= SNRs(1),1,'first');
idx_SNR2_svt = find(labels_svt.SNR <= SNRs(2),1,'last');
inds_SNR_svt = idx_SNR1_svt:idx_SNR2_svt;
idx_SNR1_opt = find(labels_opt.SNR >= SNRs(1),1,'first');
idx_SNR2_opt = find(labels_opt.SNR <= SNRs(2),1,'last');
inds_SNR_opt = idx_SNR1_opt:idx_SNR2_opt;

% Generate heatmaps
typeinds = struct('X',1,'ROI1',2,'ROI2',3,'ROI3',4,'ROI',5);
for i = 1:length(types)
    % Get plotting index
    ind = typeinds.(types{i});
    
    % Create figure
    figure('name',sprintf('nc = %i',nc),'Position',pos);
    [ni nj] = bestSubplotShape(3);
    
    % Compute master axis limits
    fcn = @(clim,X) [min(clim(1),min(X(:))) max(clim(2),max(X(:)))];
    clim = [inf -inf];
    clim = fcn(clim,NRMSE_fft{ind}(idx_nc_fft,:,inds_SNR_fft));
    clim = fcn(clim,NRMSE_svt{ind}(idx_nc_svt,:,inds_SNR_svt));
    clim = fcn(clim,NRMSE_opt{ind}(idx_nc_opt,:,inds_SNR_opt));
    
    % FFT
    subplot(ni,nj,1);
    pcolor(labels_fft.nt,labels_fft.SNR(inds_SNR_fft),squeeze(NRMSE_fft{ind}(idx_nc_fft,:,inds_SNR_fft))');
    colorbar;
    caxis(clim);
    shading interp;
    box on;
    xlabel('# frames');
    ylabel('SNR (dB)');
    title(sprintf('FFT - NRMSE [%s]',types{i}));
    
    % SVT
    subplot(ni,nj,2);
    pcolor(labels_svt.nt,labels_svt.SNR(inds_SNR_svt),squeeze(NRMSE_svt{ind}(idx_nc_svt,:,inds_SNR_svt))');
    colorbar;
    caxis(clim);
    shading interp;
    box on;
    xlabel('# frames');
    ylabel('SNR (dB)');
    title(sprintf('SVT - NRMSE [%s]',types{i}));
    
    % OptShrink
    subplot(ni,nj,3);
    pcolor(labels_opt.nt,labels_opt.SNR(inds_SNR_opt),squeeze(NRMSE_opt{ind}(idx_nc_opt,:,inds_SNR_opt))');
    colorbar;
    caxis(clim);
    shading interp;
    box on;
    xlabel('# frames');
    ylabel('SNR (dB)');
    title(sprintf('OptShrink - NRMSE [%s]',types{i}));
    
    % Save figure
    %export_fig -pdf -transparent NRMSE_LaS_nc_8_X
    %export_fig -pdf -transparent NRMSE_LaS_nc_8_ROI
end

%% Plot fixed (nc,SNR) curves

% Knobs
nc = 8; % # coils
%SNR = {10 20 30 40 50 60}; % SNR
%SNR = {40 45 50 55 60}; % SNR
SNR = {55}; % SNR
types = {'X','ROI','ROI1','ROI2','ROI3'}; % Types of curves to plot
NRMSEpath = './nrmse3/NRMSE_LpS_soft_3.mat';
%NRMSEpath = './nrmse3/NRMSE_LpS_soft_T_3.mat';
%NRMSEpath = './nrmse3/NRMSE_LpS_mixed21_3.mat';
%NRMSEpath = './nrmse3/NRMSE_LaS_soft_3.mat';
%NRMSEpath = './nrmse3/NRMSE_LaS_soft_T_3.mat';
%NRMSEpath = './nrmse3/NRMSE_LaS_mixed21_3.mat';
%dim = [1100 550]; % Figure dimensions
dim = [1000 500]; % Figure dimensions

% Load NRMSE data
load(NRMSEpath);

% Compute figure position
scrsz = get(0,'ScreenSize');
xyc = 0.5 * scrsz(3:4);
pos = [(xyc - 0.5 * dim) dim];

% Locate desired nc
[~,idx_nc_fft] = min(abs(labels_fft.nc - nc));
[~,idx_nc_svt] = min(abs(labels_svt.nc - nc));
[~,idx_nc_opt] = min(abs(labels_opt.nc - nc));

% Generate NRMSE curves
for j = 1:length(SNR)
    % Locate desired SNR
    [~,idx_SNR_fft] = min(abs(labels_fft.SNR - SNR{j}));
    [~,idx_SNR_svt] = min(abs(labels_svt.SNR - SNR{j}));
    [~,idx_SNR_opt] = min(abs(labels_opt.SNR - SNR{j}));
    
    % Create figure
    figure('name',sprintf('nc = %i, SNR = %i',nc,SNR{j}),'Position',pos);
    Ntypes = length(types);
    [ni nj] = bestSubplotShape(Ntypes);
    
    % Generate curves
    typeinds = struct('X',1,'ROI1',2,'ROI2',3,'ROI3',4,'ROI',5);
    for i = 1:Ntypes
        % Initialize subplot
        subplot(ni,nj,i);
        hold on;
        
        % Get plotting index
        ind = typeinds.(types{i});
        
        % Plot curves
        p(1) = plot(labels_fft.nt,squeeze(NRMSE_fft{ind}(idx_nc_fft,:,idx_SNR_fft))','g-o'); % FFT
        p(2) = plot(labels_svt.nt,squeeze(NRMSE_svt{ind}(idx_nc_svt,:,idx_SNR_svt))','b-o'); % SVT
        p(3) = plot(labels_opt.nt,squeeze(NRMSE_opt{ind}(idx_nc_opt,:,idx_SNR_opt))','r-o'); % OptShrink
        
        % Format axes
        box on;
        xlabel('# frames');
        ylabel('NRMSE');
        title(sprintf('NRMSE [%s]',types{i}));
        if (i == ceil(0.5 * Ntypes))
            % Add legend (center plot only)
            legend(p,'FFT','SVT','OptShrink','Location','Best');
        end
    end
end

% Save figure
%export_fig -pdf -transparent NRMSE_LpS_soft_SNR_55
%export_fig -pdf -transparent NRMSE_LpS_soft_T_SNR_55
%export_fig -pdf -transparent NRMSE_LpS_mixed21_SNR_55
%export_fig -pdf -transparent NRMSE_LaS_soft_SNR_55
%export_fig -pdf -transparent NRMSE_LaS_soft_T_SNR_55
%export_fig -pdf -transparent NRMSE_LaS_mixed21_SNR_55

%% Robustness to rank parameter

erase;

%--------------------------------------------------------------------------
% Knobs
%--------------------------------------------------------------------------
% Data generation
nc = 8;
nt = 21;
SNR = [40 45 50 55 60];

% Algo info
lambdaL_max = inf; % {2,inf}
%lambdaS = 0.01; % Fixed sparsity parameter
region = 'X'; % {'X','ROI1','ROI2','ROI3','ROI'}

% NRMSE info
NRMSEpath = './nrmse3/NRMSE_LpS_soft_full_3.mat';
%NRMSEpath = './nrmse3/NRMSE_LpS_soft_T_full_3.mat';
%NRMSEpath = './nrmse3/NRMSE_LpS_mixed21_full_3.mat';
%NRMSEpath = './nrmse3/NRMSE_LaS_soft_full_3.mat';
%NRMSEpath = './nrmse3/NRMSE_LaS_soft_T_full_3.mat';
%NRMSEpath = './nrmse3/NRMSE_LaS_mixed21_full_3.mat';

% Asymptote info
%constPath = 'constants.mat'; % Constants path
%constPath = 'constants2.mat'; % Constants path
%constPath = 'constants3.mat'; % Constants path

% Plot info
gapx = 0.5;
gapy = 0.025;
%--------------------------------------------------------------------------

% Get info
kvps = struct('X',1,'ROI1',2,'ROI2',3,'ROI3',4,'ROI',5);
metric = kvps.(region);
if ~isempty(strfind(NRMSEpath,'LpS'))
    algoStr = 'LpS';
else
    algoStr = 'LaS';
end
if ~isempty(strfind(NRMSEpath,'soft'))
    proxS = 'soft';
else
    proxS = 'mixed21';
end

%--------------------------------------------------------------------------
% Get lambdaS values
%--------------------------------------------------------------------------
nSNR = length(SNR);
if exist('lambdaS','var')
    % Fixed lambdaS value
    lambdaS_opt = repmat(lambdaS,1,nSNR);
    lambdaS_svt = repmat(lambdaS,1,nSNR);
    
    % Title strings
    lStr = 'fixed';
    titlestr_opt = sprintf('OptShrink [nt = %i, \\lambda_S = %01.03f]',nt,lambdaS);
    titlestr_svt = sprintf('SVT [nt = %i, \\lambda_S = %01.03f]',nt,lambdaS);
else
    % mNRMSE lambdaS values
    data = load(regexprep(NRMSEpath,'_full','')); % NRMSE data
    lambdaS_opt = zeros(1,nSNR);
    lambdaS_svt = zeros(1,nSNR);
    for i = 1:nSNR
        % OptShrink
        [~,nc_mmse_opt_idx] = min(abs(data.labels_opt.nc - nc));
        [~,nt_mmse_opt_idx] = min(abs(data.labels_opt.nt - nt));
        [~,SNR_mmse_opt_idx] = min(abs(data.labels_opt.SNR - SNR(i)));
        inds_mmse_opt = data.idxopt_opt{metric}{nc_mmse_opt_idx,nt_mmse_opt_idx,SNR_mmse_opt_idx};
        lambdaS_opt(i) = data.labels_opt.lambdaS(inds_mmse_opt(2));
        
        % SVT
        [~,nc_svt_idx] = min(abs(data.labels_svt.nc - nc));
        [~,nt_svt_idx] = min(abs(data.labels_svt.nt - nt));
        [~,SNR_svt_idx] = min(abs(data.labels_svt.SNR - SNR(i)));
        inds_mmse_svt = data.idxopt_svt{metric}{nc_svt_idx,nt_svt_idx,SNR_svt_idx};
        lambdaS_svt(i) = data.labels_svt.lambdaS(inds_mmse_svt(2));
    end
    
    % Title strings
    lStr = 'mmse';
    titlestr_opt = sprintf('OptShrink [nt = %i, \\lambda_S = OPT]',nt);
    titlestr_svt = sprintf('SVT [nt = %i, \\lambda_S = OPT]',nt);
end
%--------------------------------------------------------------------------

% Load data
load(NRMSEpath);

% Initialize figure
cm = hsv(nSNR + 1);
ylim = [inf -inf];
%figure('Position',[750 500 560 230]);
figure('Position',[750 500 685 280]);
%figure;

% Loop over SNRs
phndl = zeros(1,nSNR);
pstr = cell(1,nSNR);
for i = 1:nSNR
    % Get OptShrink NRMSEs
    [~,nc_opt_idx] = min(abs(labels_opt.nc - nc));
    [~,nt_opt_idx] = min(abs(labels_opt.nt - nt));
    [~,SNR_opt_idx] = min(abs(labels_opt.SNR - SNR(i)));
    [~,lambdaS_opt_idx] = min(abs(labels_opt.lambdaS - lambdaS_opt(i)));
    NRMSEopt = squeeze(NRMSE_opt{metric}(nc_opt_idx,nt_opt_idx,SNR_opt_idx,1,:,lambdaS_opt_idx));
    r = labels_opt.r;
    
    % Get SVT NRMSEs
    [~,nc_svt_idx] = min(abs(labels_svt.nc - nc));
    [~,nt_svt_idx] = min(abs(labels_svt.nt - nt));
    [~,SNR_svt_idx] = min(abs(labels_svt.SNR - SNR(i)));
    [~,lambdaS_svt_idx] = min(abs(labels_svt.lambdaS - lambdaS_svt(i)));
    NRMSEsvt = squeeze(NRMSE_svt{metric}(nc_svt_idx,nt_svt_idx,SNR_svt_idx,1,:,lambdaS_svt_idx));
    lambdaL = labels_svt.lambdaL;
    
    % Clip to lambdaL range
    indsL = (lambdaL < lambdaL_max);
    NRMSEsvt = NRMSEsvt(indsL);
    lambdaL = lambdaL(indsL);
    
    % Compute common y-limits
    ylim = [min([ylim(1); NRMSEopt(:); NRMSEsvt(:)]) ...
            max([ylim(2); NRMSEopt(:); NRMSEsvt(:)])];
    
    % OptShrink
    subplot(1,2,1);
    phndl(i) = plot(r,NRMSEopt,'-o','Color',cm(i,:));
    pstr{i} = sprintf('SNR = %i dB',SNR(i));
    hold on;
    
    % SVT
    subplot(1,2,2);
    semilogx(lambdaL,NRMSEsvt,'-o','Color',cm(i,:));
    hold on;
end
xlim_svt = [min(lambdaL(:)) max(lambdaL(:))];

% Add optimal ratio asymptotes
if exist('constPath','var')
    cdata = load(constPath);
    subplot(1,2,2);
    %for i = 1:nSNR
    i = round(nSNR / 2);
    
    % Compute asymptote
    n = cdata.constants(nt).n;
    CL = cdata.constants(nt).CL;
    CS = cdata.constants(nt).(proxS);
    lLa = sqrt(n) * CS * lambdaS_svt(i) / CL;
    
    % Update xlim
    if lLa < xlim_svt(1)
        xlim_svt(1) = (1 - gapx) * lLa;
    elseif lLa > xlim_svt(2)
        xlim_svt(2) = (1 + gapx) * lLa;
    end
    
    % Add asymptote to graph
    %plot([lLa lLa],ylim + gap * [-1 1],'--','Color',cm(i,:));
    plot([lLa lLa],ylim + gapy * [-1 1],'k--');
    %end
end

% Finalize OptShrink subplot
subplot(1,2,1);
xlabel('r');
ylabel(sprintf('NRMSE(%s)',region));
title(titlestr_opt);
legend(phndl,pstr{:});
set(gca,'XLim',[0 max(r(:))]);
set(gca,'YLim',ylim + gapy * [-1 1]);

% Finalize SVT subplot
subplot(1,2,2);
xlabel('\lambda_L');
ylabel(sprintf('NRMSE(%s)',region));
title(titlestr_svt);
set(gca,'XLim',xlim_svt);
set(gca,'YLim',ylim + gapy * [-1 1]);

% Save figure
drawnow; pause(0.1);
figStr = sprintf('sensitivity_%s_%s_%s_%s',algoStr,proxS,region,lStr);
if exist('constPath','var')
    figStr = [figStr '_asymp'];
end
export_fig('-pdf','-transparent',figStr);

%% Fixed (nc,nt,SNR) tests

% Knobs
nc = 8;
nt = 20;
SNR = 50;
seed = 1;
algoStr = 'L+S PGM';
T = TempFFT(2); % {1,TempFFT(2)}
proxS = 'soft';
%inpath = 'dynobj.mat';
%inpath = 'dynobj2.mat';
inpath = 'dynobj3.mat';

% Parameter values
r = 1:5;
lambdaL = logspace(log10(0.001),log10(10),5);
lambdaS = logspace(log10(0.001),log10(10),5);

% Load data
gtData = load(inpath);
Xtrue_full = gtData.dyn_obj;
dce = gtData.dce;
clear gtData;
nf = size(Xtrue_full,3);

% Loop over parameter values
M = length(r); % Assumed same as length(lambdaL)
N = length(lambdaS);
NRMSEopt = nan(M,N,2);
NRMSEsvt = nan(M,N,2);
for i = M:-1:1
    for j = 1:N
        % OptShrink reconstruction
        %[~,recon_opt,~,~,mask,ROIs,~,splineInterp] = run_optshrink_algo(nc,nt,SNR,seed,r(i),proxS,lambdaS(j),algoStr,T,Xtrue_full,dce);
        
        % SVT reconstruction
        [~,recon_svt,~,~,mask,ROIs,~,splineInterp] = run_svt_algo(nc,nt,SNR,seed,lambdaL(i),proxS,lambdaS(j),algoStr,T,Xtrue_full,dce);
        
        % Interpolate results
        %Xhat_opt = splineInterp(abs(embed(recon_opt.X,mask)));
        Xhat_svt = splineInterp(abs(embed(recon_svt.X,mask)));
        
        % Compute NRMSEs
        NRMSEfcn = @(Xhat,Xtrue) norm(Xhat(:) - Xtrue(:)) / norm(Xtrue(:));
        %NRMSEopt(i,j,1) = NRMSEfcn(Xhat_opt,Xtrue_full);
        NRMSEsvt(i,j,1) = NRMSEfcn(Xhat_svt,Xtrue_full);
        M3 = repmat(ROIs{3},[1 1 nf]);
        %NRMSEopt(i,j,2) = NRMSEfcn(Xhat_opt(M3),Xtrue_full(M3));
        NRMSEsvt(i,j,2) = NRMSEfcn(Xhat_svt(M3),Xtrue_full(M3));
    end
end

% Find mNRMSE parameters
vec = @(X)X(:);
%NRMSEopt(isnan(NRMSEopt)) = inf;
NRMSEsvt(isnan(NRMSEsvt)) = inf;
%[nrmseOpt1 Opt1idx] = min(vec(NRMSEopt(:,:,1)));
%[nrmseOpt4 Opt4idx] = min(vec(NRMSEopt(:,:,2)));
[nrmseSVT1 SVT1idx] = min(vec(NRMSEsvt(:,:,1)));
[nrmseSVT4 SVT4idx] = min(vec(NRMSEsvt(:,:,2)));
%[Opt1i Opt1j] = ind2sub([M N],Opt1idx);
%[Opt4i Opt4j] = ind2sub([M N],Opt4idx);
[SVT1i SVT1j] = ind2sub([M N],SVT1idx);
[SVT4i SVT4j] = ind2sub([M N],SVT4idx);
%Opt1_results = [r(Opt1i) lambdaS(Opt1j) nrmseOpt1] %#ok
%Opt4_results = [r(Opt4i) lambdaS(Opt4j) nrmseOpt4] %#ok
SVT1_results = [lambdaL(SVT1i) lambdaS(SVT1j) nrmseSVT1] %#ok
SVT4_results = [lambdaL(SVT4i) lambdaS(SVT4j) nrmseSVT4] %#ok

%% Spot check mNRMSE reconstructions

%--------------------------------------------------------------------------
% Knobs
%--------------------------------------------------------------------------
% Data generation
nc = 8;
nt = 21;
SNR = 55;
seed = 1;
proxS = 'mixed21';
%inpath = 'dynobj.mat';
%inpath = 'dynobj2.mat';
inpath = 'dynobj3.mat';

%{
% L + S
algoStr = 'L+S PGM';
NRMSEpath = './nrmse3/NRMSE_LpS_mixed21_3.mat';
T = 1; % {1,TempFFT(2)}
%}

% L & S
algoStr = 'L&S ADMM';
NRMSEpath = './nrmse3/NRMSE_LaS_mixed21_3.mat';
T = 1; % {1,TempFFT(2)}

% NRMSE metric
metric = 1; % {1 = X, 2 = lesion #1, ..., 5 = ROI}
%--------------------------------------------------------------------------

% Load data
load(NRMSEpath);
gtData = load(inpath);
Xtrue_full = gtData.dyn_obj;
dce = gtData.dce;
clear gtData;
nf = size(Xtrue_full,3);

% Set parameters manually
%r = 1;
%lambdaL = 0.37276;
%lambdaS1 = 1;
%lambdaS2 = 1;

%--------------------------------------------------------------------------
% OptShrink parameters
%--------------------------------------------------------------------------
% Locate optimal indices
[~,nc_opt_idx] = min(abs(labels_opt.nc - nc));
[~,nt_opt_idx] = min(abs(labels_opt.nt - nt));
[~,SNR_opt_idx] = min(abs(labels_opt.SNR - SNR));
inds = idxopt_opt{metric}{nc_opt_idx,nt_opt_idx,SNR_opt_idx};

% Get optimal parameters
r = labels_opt.r(inds(1));
lambdaS1 = labels_opt.lambdaS(inds(2));
OptShrink_params = [r lambdaS1];
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% SVT parameters
%--------------------------------------------------------------------------
% Locate optimal indices
[~,nc_svt_idx] = min(abs(labels_svt.nc - nc));
[~,nt_svt_idx] = min(abs(labels_svt.nt - nt));
[~,SNR_svt_idx] = min(abs(labels_svt.SNR - SNR));
inds = idxopt_svt{metric}{nc_svt_idx,nt_svt_idx,SNR_svt_idx};

% Get optimal parameters
lambdaL = labels_svt.lambdaL(inds(1));
lambdaS2 = labels_svt.lambdaS(inds(2));
SVT_params = [lambdaL lambdaS2];
%--------------------------------------------------------------------------

%%

% OptShrink reconstruction
[~,recon_opt,~,Xfft,mask,ROIs,~,splineInterp] = run_optshrink_algo(nc,nt,SNR,seed,r,proxS,lambdaS1,algoStr,T,Xtrue_full,dce);

% SVT reconstruction
[~,recon_svt,~,~,~,~,~,~] = run_svt_algo(nc,nt,SNR,seed,lambdaL,proxS,lambdaS2,algoStr,T,Xtrue_full,dce);

% Interpolate results
Xhat_fft = splineInterp(abs(Xfft));
Xhat_opt = splineInterp(abs(embed(recon_opt.X,mask)));
Xhat_svt = splineInterp(abs(embed(recon_svt.X,mask)));
if strcmp(algoStr,'L+S PGM')
    % L + S only
    Lhat_opt = splineInterp(abs(embed(recon_opt.L,mask)));
    Shat_opt = splineInterp(abs(embed(recon_opt.S,mask)));
    Lhat_svt = splineInterp(abs(embed(recon_svt.L,mask)));
    Shat_svt = splineInterp(abs(embed(recon_svt.S,mask)));
end

% Compute NRMSEs
NRMSEfcn = @(Xhat,Xtrue) norm(Xhat(:) - Xtrue(:)) / norm(Xtrue(:));
switch metric
    case 1
        % X
        nrmseopt = NRMSEfcn(Xhat_opt,Xtrue_full);
        nrmsesvt = NRMSEfcn(Xhat_svt,Xtrue_full);
    case 2
        % Lesion #1
        M1 = repmat(ROIs{1},[1 1 nf]);
        nrmseopt = NRMSEfcn(Xhat_opt(M1),Xtrue_full(M1));
        nrmsesvt = NRMSEfcn(Xhat_svt(M1),Xtrue_full(M1));
    case 3
        % Lesion #2
        M2 = repmat(ROIs{2},[1 1 nf]);
        nrmseopt = NRMSEfcn(Xhat_opt(M2),Xtrue_full(M2));
        nrmsesvt = NRMSEfcn(Xhat_svt(M2),Xtrue_full(M2));
    case 4
        % Lesion #3
        M3 = repmat(ROIs{3},[1 1 nf]);
        nrmseopt = NRMSEfcn(Xhat_opt(M3),Xtrue_full(M3));
        nrmsesvt = NRMSEfcn(Xhat_svt(M3),Xtrue_full(M3));
    case 5
        % All lesions
        M1 = repmat(ROIs{1},[1 1 nf]);
        M2 = repmat(ROIs{2},[1 1 nf]);
        M3 = repmat(ROIs{3},[1 1 nf]);
        M = (M1 | M2 | M3);
        nrmseopt = NRMSEfcn(Xhat_opt(M),Xtrue_full(M));
        nrmsesvt = NRMSEfcn(Xhat_svt(M),Xtrue_full(M));
end

% Display NRMSEs
OptShrink_NRMSE = [nrmseopt NRMSE_opt{metric}(nc_opt_idx,nt_opt_idx,SNR_opt_idx)] %#ok
SVT_NRMSE = [nrmsesvt NRMSE_svt{metric}(nc_svt_idx,nt_svt_idx,SNR_svt_idx)] %#ok

% Display parameters
OptShrink_params = OptShrink_params %#ok
SVT_params = SVT_params %#ok

%% Playback reconstructions as movie

% Knobs
mag = 1.5; % Magnification
fps = 30; % Frames per second

% Concatenate data
if strcmp(algoStr,'L+S PGM')
    % L + S
    M = {Xtrue_full Xhat_opt Lhat_opt Shat_opt;
         Xhat_fft Xhat_svt Lhat_svt Shat_svt};
    xlabels = {'','L + S','L','S'};
    ylabels = {'OptShrink','SVT'};
else
    % L & S
    M = {Xtrue_full Xhat_opt;
         Xhat_fft Xhat_svt};
    xlabels = {'','X'};
    ylabels = {'OptShrink','SVT'};
end

% Format data
for i = 1:numel(M)
    % Permute rows/columns
    M{i} = permute(M{i},[2 1 3]);
    
    % Display each on full-scale
    %M{i} = M{i} - min(M{i}(:));
    %M{i} = M{i} / max(M{i}(:));
    
    % Display on Xtrue scale
    M{i} = M{i} - min(Xtrue_full(:));
    M{i} = M{i} / max(Xtrue_full(:));
    M{i} = min(max(M{i},0),1);
end
M = cat(1,cat(2,M{1,:}),cat(2,M{2,:}));

% Play movie
movie = struct();
movie.video = M;
movie.Fv = fps;
opts = struct();
opts.mag = mag;
opts.xlabels = xlabels;
opts.ylabels = ylabels;
PlayMovie(movie,opts);

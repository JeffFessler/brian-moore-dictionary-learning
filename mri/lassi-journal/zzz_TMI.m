%% [VIEW] otazo

% Knobs
ps       = 1/8;
SNR      = inf;
seed     = 1;
inpath   = 'otazo_full.mat';
roipath  = 'otazo_roi.mat';
outpath  = 'tmi_otazo_roi.png';
roiFrame = 20;

% Generate invivo data
[Y, ~, ~, Xtrue, Xfft] = generateCardiacPerfData(ps,SNR,seed,inpath);
mask = (abs(Y(:,:,:,1)) ~= 0);
lim  = getLim(abs(Xtrue));

% Visualize data
PlayMovie(cat(2,mask,Xtrue,Xfft));

%{
% Pick ROI
figure;
imshow(abs(Xtrue(:,:,roiFrame)),lim);
axis on;
%}

% Save ROI
ROI.rows = 36:89;
ROI.cols = 32:86;
save(roipath,'ROI');

% Visualize ROI
%PlayMovie(Xtrue(ROI.rows,ROI.cols,:));

% Save ROI image
Xroi = addROI(abs(Xtrue(:,:,roiFrame)),ROI,lim);
imshow(Xroi);
imwrite(Xroi,outpath,'png');

fprintf('DONE\n');

%% [VIEW] invivo

% Knobs
nLines   = 15;
SNR      = inf;
seed     = 1;
inpath   = 'invivo_full.mat';
roipath  = 'invivo_roi.mat';
outpath  = 'tmi_invivo_roi.png';
roiFrame = 25;

% Generate invivo data
[Y, ~, ~, Xtrue, Xfft] = generateInvivoData(nLines,SNR,seed,inpath);
mask = (abs(Y) ~= 0);
lim  = getLim(abs(Xtrue));

% Undersampling factor
ps = nnz(mask) / numel(mask) %#ok

% Visualize data
PlayMovie(cat(2,mask,Xtrue,Xfft));

%{
% Pick ROI
figure;
imshow(abs(Xtrue(:,:,roiFrame)));
axis on;
%}

% Save ROI
ROI.rows = 78:124;
ROI.cols = 28:72;
save(roipath,'ROI');

% Visualize ROI
%PlayMovie(Xtrue(ROI.rows,ROI.cols,:));

% Save ROI image
Xroi = addROI(abs(Xtrue(:,:,roiFrame)),ROI,lim);
Xroi = Xroi(31:160,1:90,:);         % Truncate image
Xroi = flipdim(flipdim(Xroi,1),2);  % Flip image
imshow(Xroi);
imwrite(Xroi,outpath,'png');

fprintf('DONE\n');

%% [VIEW] pincat

% Knobs
nLines   = 18;
SNR      = inf;
seed     = 1;
inpath   = 'pincat_full.mat';
roipath  = 'pincat_roi.mat';
outpath  = 'tmi_pincat_roi.png';
roiFrame = 16;

% Generate invivo data
[Y, ~, ~, Xtrue, Xfft] = generateInvivoData(nLines,SNR,seed,inpath);
mask = (abs(Y) ~= 0);
lim  = getLim(abs(Xtrue));

% Undersampling factor
ps = nnz(mask) / numel(mask) %#ok

% Visualize data
%PlayMovie(cat(2,mask,Xtrue,Xfft));

%{
% Pick ROI
figure;
imshow(abs(Xtrue(:,:,roiFrame)));
axis on;
%}

% Save ROI
ROI.rows = 35:89;
ROI.cols = 50:101;
save(roipath,'ROI');

% Visualize ROI
%PlayMovie(Xtrue(ROI.rows,ROI.cols,:));

% Save ROI image
Xroi = addROI(abs(Xtrue(:,:,roiFrame)),ROI,lim);
imshow(Xroi);
imwrite(Xroi,outpath,'png');

fprintf('DONE\n');

%% [SPOT] robustPCA

% Knobs
%{
inpath  = 'otazo_R8.mat';
outpath = 'otazo_R8_lps_mse.mat';
r       = nan;
lambdaL = 1.1955;
lambdaS = 0.01;
nIters  = 250;
%}
%{
inpath  = 'invivo_R8.mat';
outpath = 'invivo_R8_lps_nrmse.mat';
r       = nan;
lambdaL = 0.178;
lambdaS = 0.003;
nIters  = 250;
%}
%{
inpath  = 'pincat_R8.mat';
outpath = 'pincat_R8_lps_nrmse.mat';
r       = nan;
lambdaL = 0.2069;
lambdaS = 0.0030;
nIters  = 250;
%}
inpath  = 'pincat_R8n.mat';
outpath = 'pincat_R8n_lps_nrmse.mat';
r       = nan;
lambdaL = 3.793;
lambdaS = 0.15; % or 1e100
nIters  = 250;

% Load dependencies
addpath('./deps_rpca');

% Load undersampled data
load(inpath);
[ny, nx, nt] = size(Xtrue);
%A = Emat_xyt(mask,samp,[ny, nx, nt]);
A  = Afft(mask,[ny, nx, nt]);
T  = TempFFT(3,[ny, nx, nt]);

% Standard back-propagation
%Xfft = reshape(A' * Y,[ny, nx, nt]);

% Robust PCA
opts.A      = A;
opts.T      = T;
opts.r      = r;
opts.nIters = nIters;
opts.L0     = reshape(Xfft,[],nt);
opts.S0     = zeros(ny * nx,nt);
opts.Xtrue  = reshape(Xtrue,[],nt);
opts.accel  = false;
opts.tau    = 1;
opts.flag   = 1;
[Lhat, Shat, stats] = robustPCA(Y,lambdaL,lambdaS,opts);
Lhat = reshape(Lhat,ny,nx,nt);
Shat = reshape(Shat,ny,nx,nt);

%{
% Play movie
optsPM.xlabels = {'Xtrue','Xfft','Lhat + Shat','Lhat','Shat'};
PlayMovie(cat(2,Xtrue,Xfft,Lhat + Shat,Lhat,Shat),optsPM);
%}

% Save results, if requested
if exist('outpath','var') && ~isempty(outpath)
    time  = stats.time;
    cost  = stats.cost;
    nrmse = stats.nrmse;
    delta = stats.delta;
    save(outpath,'lambdaL','lambdaS','r','nIters', ...
                 'time','cost','nrmse','delta','Lhat','Shat');
end

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
cfigure();

% Cost
subplot(1,2,1);
plot(1:nIters,stats.cost,'b-o');
xlabel('Iteration');
title('Cost');
axis tight; padAxis();

% NRMSE
subplot(1,2,2);
plot(1:nIters,stats.nrmse,'b-o');
xlabel('Iteration');
title('NRMSE');
axis tight; padAxis();
%--------------------------------------------------------------------------

%% [SPOT] drpca

% Knobs
%inpath   = 'otazo_R8.mat';
%initpath = 'otazo_R8_lps_visual.mat';
%lambdaL  = 0.5;
%lambdaS  = 0.02;
%lambdaB  = logspace(log10(0.06),log10(0.03),nIters);
%dr       = 1;
inpath    = 'invivo_R8.mat';
initpath  = 'invivo_R8_lps_nrmse.mat';
lambdaL   = 0.5;
lambdaS   = 0.02;
lambdaB   = 0.03;
dr        = 1;
nIters    = 3;

% Load dependencies
addpath('./deps_lassi');

% Load undersampled data
load(inpath);
[ny, nx, nt] = size(Xtrue);

% Load L+S initialization
init = load(initpath);
L0   = init.Lhat;
S0   = init.Shat;

% DRPCA
%opts.A       = Emat_xyt(mask,samp,[ny, nx, nt]);
opts.A        = Afft(mask,[ny, nx, nt]);
opts.sdim     = [ny, nx, nt];
opts.pdim     = [8, 8, 5];
opts.pgap     = [2, 2, 2];
opts.type     = 'hard';
opts.dr       = dr;
opts.ddim     = [64, 5];
opts.nIters   = nIters;
opts.nItersDB = 1;
opts.nItersLS = 5;
opts.L0       = reshape(L0,[],nt);
opts.S0       = reshape(S0,[],nt);
opts.D0       = dctmtx(prod(opts.pdim));
opts.Xtrue    = reshape(Xtrue,[],nt);
opts.accel    = false;
opts.tau      = 1;
opts.flag     = 1;
[Lhat, Shat, ~, ~,stats] = drpca(Y,lambdaL,lambdaS,lambdaB,opts);
Lhat = reshape(Lhat,ny,nx,nt);
Shat = reshape(Shat,ny,nx,nt);

% Play movie
optsPM.xlabels = {'Xtrue','Xfft','Lhat + Shat','Lhat','Shat'};
PlayMovie(cat(2,Xtrue,Xfft,Lhat + Shat,Lhat,Shat),optsPM);

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
cfigure([950, 250]);

% Cost
subplot(1,3,1);
plot(stats.cost,'b-o');
xlabel('Iteration');
title('Cost');
axis tight; padAxis();

% NRMSE
subplot(1,3,2);
plot(stats.nrmse,'b-o');
xlabel('Iteration');
title('NRMSE');
axis tight; padAxis();

% Sparsity
subplot(1,3,3);
plot(stats.sparsity,'b-o');
xlabel('Iteration');
title('Sparsity');
axis tight; padAxis();
%--------------------------------------------------------------------------

%% [SPOT] ktslr

%{
% Knobs
inpath    = 'otazo_R8.mat';
outpath   = 'otazo_R8_ktslr_nrmse.mat';
p         = 1;
mu1       = 0.2;
mu2       = 0.003;
%}
%{
inpath    = 'invivo_R8.mat';
outpath   = 'invivo_R8_ktslr_nrmse.mat';
p         = 1;
mu1       = 0.002;
mu2       = 0.001;
%}
%{
inpath    = 'pincat_R8.mat';
%outpath  = 'pincat_R8_ktslr_nrmse.mat';
p         = nan;
mu1       = nan
mu2       = nan;
%}
%{
inpath    = 'pincat_R8n.mat';
outpath   = 'pincat_R8n_ktslr_nrmse.mat';
p         = 1;
mu1       = 2;
mu2       = 0.01;
%}
nItersO   = 5;
nItersI   = 10;
nItersCG  = 5;
beta1     = 1;
beta2     = 1;
beta1rate = 50;
beta2rate = 25;
stepSize  = [1, 1, 0.303];

% Load dependencies
addpath('./deps_ktslr');

% Load undersampled data
load(inpath);
[ny, nx, nt] = size(Xtrue);
%A = Emat_xyt(mask,samp,[ny, nx, nt]);
A  = Afft(mask,[ny, nx, nt]);

% Standard back-propagation
%Xfft = reshape(A' * Y,[ny, nx, nt]);

% Run k-t SLR
opts.nItersO   = nItersO;
opts.nItersI   = nItersI;
opts.nItersCG  = nItersCG;
opts.beta1     = beta1;
opts.beta2     = beta2;
opts.beta1rate = beta1rate;
opts.beta2rate = beta2rate;
opts.stepSize  = stepSize;
opts.Xtrue     = Xtrue;
[Xhat, stats]  = ktSLR(Y,A,p,mu1,mu2,Xfft,opts);

%{
% Play movie
optsPM.xlabels = {'Xtrue','Xfft','Xhat'};
PlayMovie(cat(2,Xtrue,Xfft,Xhat),optsPM);
%}

% Save results, if requested
if exist('outpath','var') && ~isempty(outpath)
    time  = stats.time;
    cost  = stats.cost;
    nrmse = stats.nrmse;
    delta = stats.delta;
    save(outpath,'p','mu1','mu2','nItersO','nItersI','nItersCG', ...
                 'beta1','beta2','beta1rate','beta2rate','stepSize', ...
                 'time','cost','nrmse','delta','Xhat');
end

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
cfigure();

% Cost
subplot(1,2,1);
plot(1:stats.nIters,stats.cost,'b-o');
xlabel('Iteration');
title('Cost');
axis tight; padAxis();

% NRMSE
subplot(1,2,2);
plot(1:stats.nIters,stats.nrmse,'b-o');
xlabel('Iteration');
title('NRMSE');
axis tight; padAxis();
%--------------------------------------------------------------------------

%% [PLOT] LASSI R8

% Knobs
%inpath  = 'LASSI_R8_par1.mat';
%inpath  = 'LASSI_R8_par2.mat';
%inpath  = 'LASSI_R8_par3.mat';
%inpath  = 'LASSI_R8_par4.mat';
%inpath  = 'LASSI_R8_par1b.mat';
%inpath  = 'LASSI_R8_par2b.mat';
%inpath  = 'LASSI_R8_par2c.mat';
%inpath  = 'LASSI_R8_par3b.mat';
%inpath  = 'LASSI_R8_par3c.mat';
%inpath  = 'LASSI_R8_par4b.mat';
%inpath  = 'invivo_LASSI_R8_par1.mat';
%inpath  = 'invivo_LASSI_R8_par1b.mat';
%inpath  = 'invivo_LASSI_R8_par2.mat';
%inpath  = 'invivo_LASSI_R8_par3.mat';
%inpath  = 'invivo_LASSI_R8_par4.mat';
%inpath  = 'pincat_LASSI_R8_par1.mat';
%inpath  = 'pincat_LASSI_R8_par1b.mat';
%inpath  = 'pincat_LASSI_R8_par2.mat';
%inpath  = 'pincat_LASSI_R8_par3.mat';
%inpath  = 'pincat_LASSI_R8_par4.mat';
%inpath  = 'pincat_LASSI_R8_par1n.mat';
%inpath  = 'pincat_LASSI_R8_par2n.mat';
%inpath  = 'pincat_LASSI_R8_par3n.mat';
%inpath  = 'pincat_LASSI_R8_par4n.mat';
%inpath  = 'pincat_LASSI_R8_par1nb.mat';
%inpath  = 'pincat_LASSI_R8_par2nb.mat';
%inpath  = 'pincat_LASSI_R8_par3nb.mat';
%inpath  = 'pincat_LASSI_R8_par4nb.mat';
%inpath  = 'pincat_LASSI_R8_par1nc.mat';
%inpath  = 'pincat_LASSI_R8_par1nd.mat';
%inpath  = 'pincat_LASSI_R8_par1ne.mat';
%inpath  = 'pincat_LASSI_R8_par1nf.mat';
%inpath  = 'pincat_LASSI_R8_par1ne2.mat';
%inpath  = 'pincat_LASSI_R8_par1nf2.mat';
%inpath  = 'pincat_LASSI_R8_par1ng.mat';
%inpath  = 'pincat_LASSI_R8_par1nh.mat';
%inpath  = 'invivo_LASSI_R8_par1c.mat';
%inpath  = 'invivo_LASSI_R8_par1d.mat';
%inpath  = 'invivo_LASSI_R8_Sonly_par2.mat';
%inpath  = 'pincat_LASSI_R8_Sonly_par1n.mat';
%inpath  = 'pincat_LASSI_R8_Sonly_par2n.mat';
%inpath  = 'invivo_LASSI_R8_par1e.mat';
%inpath  = 'invivo_LASSI_R8_par1f.mat';
%inpath  = 'invivo_LASSI_R8_Sonly_par3.mat';
%inpath  = 'pincat_LASSI_R8_par1ng2.mat';
%inpath  = 'pincat_LASSI_R8_par1nh2.mat';
%inpath  = 'invivo_LASSI_R8_par1c2.mat';
%inpath  = 'invivo_LASSI_R8_par1d2.mat';
%inpath  = 'invivo_LASSI_R8_par1c3.mat';
%inpath  = 'LASSI_R8_par1c.mat';
%inpath  = 'LASSI_R8_par1d.mat';
%inpath  = 'LASSI_R8_par1e.mat';
%inpath  = 'LASSI_R8_par1f.mat';
%inpath  = 'LASSI_R8_par1f2.mat';
%inpath  = 'LASSI_R8_par1g.mat';
%inpath  = 'LASSI_R8_par1h.mat';
%inpath  = 'LASSI_R8_par1h2.mat';
%inpath  = 'LASSI_R8_par1i.mat';
%inpath  = 'LASSI_R8_par5.mat';
%inpath  = 'LASSI_R8_par5b.mat';
%inpath  = 'LASSI_R8_par6.mat';
%inpath  = 'LASSI_R8_par6b.mat';
%inpath  = 'LASSI_R8_par7.mat';
%inpath  = 'LASSI_R8_par7b.mat';
%inpath  = 'LASSI_R8_par7c.mat';
%inpath  = 'LASSI_R8_par8.mat';
%inpath  = 'LASSI_R8_par8b.mat';
%inpath  = 'LASSI_R8_par9.mat';
SAVE_FIGURE = true;

% Load data (Lhat, Shat, vars, stats)
load(inpath);

% Compute optimal parameters
nrmse = stats.nrmse(:,:,:,:,:,end);
[mnrmse, idx] = min(nrmse(:));
[ir, il, is, ib, id] = ind2sub(size(nrmse),idx);

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
%cfigure([1523, 640]); % Windows
cfigure([1141, 500]); % Mac
ax = [];

% r/lambdaL
if isnan(vars.lambdaL(1))
    % OptShrink
    subplot(2,4,1);
    if (numel(vars.r) > 1), ax(end + 1) = gca(); end
    plot(vars.r,squeeze(nrmse(:,il,is,ib,id)),'b-o');
    xlabel('r');
    ylabel('NRMSE');
    title(sprintf('r = %d',vars.r(ir)));
    axis tight; padAxis();
elseif all(vars.lambdaL < 1e10)
    % SVT
    subplot(2,4,1);
    if (numel(vars.lambdaL) > 1), ax(end + 1) = gca(); end
    plot(vars.lambdaL,squeeze(nrmse(ir,:,is,ib,id)),'b-o');
    xlabel('\lambda_L');
    ylabel('NRMSE');
    title(sprintf('\\lambda_L = %.3f',vars.lambdaL(il)));
    axis tight; padAxis();
end

% lambdaS
subplot(2,4,2);
if (numel(vars.lambdaS) > 1), ax(end + 1) = gca(); end
plot(vars.lambdaS,squeeze(nrmse(ir,il,:,ib,id)),'b-o');
xlabel('\lambda_S');
ylabel('NRMSE');
title(sprintf('\\lambda_S = %.3f',vars.lambdaS(is)));
axis tight; padAxis();

% lambdaB
subplot(2,4,3);
if (numel(vars.lambdaB) > 1), ax(end + 1) = gca(); end
plot(vars.lambdaB,squeeze(nrmse(ir,il,is,:,id)),'b-o');
xlabel('\lambda_B');
ylabel('NRMSE');
title(sprintf('\\lambda_B = %.3f',vars.lambdaB(ib)));
axis tight; padAxis();

% dr
subplot(2,4,4);
if (numel(vars.dr) > 1), ax(end + 1) = gca(); end
plot(vars.dr,squeeze(nrmse(ir,il,is,ib,:)),'b-o');
xlabel('r_d');
ylabel('NRMSE');
title(sprintf('r_d = %d',vars.dr(id)));
axis tight; padAxis();

% Optimal cost trajectory
subplot(2,4,5);
plot(1:vars.nIters,squeeze(stats.cost(ir,il,is,ib,id,:)),'b-');
xlabel('Iteration');
ylabel('cost');
title('Optimal cost');
axis tight; padAxis();

% Optimal MSE trajectory
subplot(2,4,6);
plot(1:vars.nIters,squeeze(stats.nrmse(ir,il,is,ib,id,:)),'b-');
xlabel('Iteration');
ylabel('NRMSE');
title(sprintf('Optimal NRMSE (%.2f%%)',100 * mnrmse));
axis tight; padAxis();

% Optimal delta
subplot(2,4,7);
semilogy(1:vars.nIters,squeeze(stats.delta(ir,il,is,ib,id,:)),'b-');
xlabel('Iteration');
ylabel('delta');
title('Optimal delta');
axis tight; padAxis();

% Optimal sparsity trajectory
subplot(2,4,8);
plot(1:vars.nIters,squeeze(stats.sparsity(ir,il,is,ib,id,:)),'b-');
xlabel('Iteration');
ylabel('Sparsity %');
title('Optimal sparsity');
axis tight; padAxis();

% Match axes
if ~isempty(ax), matchAxes(ax,[],'y'); end

% Save figure
if SAVE_FIGURE
    [~, name] = fileparts(inpath);
    export_fig('-pdf','-transparent',[name, '.pdf']);
    fprintf('DONE\n');
end
%--------------------------------------------------------------------------

%% [PLOT] LASSI R8 (overcompleteness)

% Knobs
inpath = 'LASSI_R8_Dover2_par1.mat';
SAVE_FIGURE = true;

% Load data
load(inpath);
nrmse = squeeze(stats.nrmse);
nrmse = nrmse(:,:,end);
nd    = size(nrmse,2);

% Compute optimal parameters
[~, idxOpt] = min(nrmse,[],1);
lbOpt       = vars.lambdaB(idxOpt);

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
cfigure([817, 450]);

ax       = zeros(1,nd);
[nh, nw] = bestSubplotShape(nd);
for i = 1:nd
    ax(i) = subplot(nh,nw,i);
    plot(vars.lambdaB,nrmse(:,i),'b-o');
    xlabel('\lambda_B');
    ylabel('NRMSE');
    title(sprintf('NRMSE = %.2f%%, \\lambda_B = %.3f, np = %d',100 * nrmse(idxOpt(i),i),lbOpt(i),vars.np(i)));
    axis tight;
end
matchAxes(ax,[],'y');

% Save figure
if SAVE_FIGURE
    [~, name] = fileparts(inpath);
    export_fig('-pdf','-transparent',[name, '.pdf']);
    fprintf('DONE\n');
end
%--------------------------------------------------------------------------

%% [PLOT] RPCA R8

% Knobs
%inpath = 'RPCA_R8_par1.mat';
%inpath = 'invivo_RPCA_R8_par1.mat';
%inpath = 'invivo_RPCA_R8_par2.mat';
%inpath = 'invivo_RPCA_R8_par3.mat';
%inpath = 'invivo_RPCA_R8_par4.mat';
%inpath = 'pincat_RPCA_R8_par1.mat';
%inpath = 'pincat_RPCA_R8_par2.mat';
%inpath = 'pincat_RPCA_R8_par1n.mat';
%inpath = 'pincat_RPCA_R8_par2n.mat';
%inpath = 'invivo_RPCA_R8_par1b.mat';
SAVE_FIGURE = true;

% Load data (Lhat, Shat, vars, stats)
load(inpath);

% Compute optimal parameters
nrmse = stats.nrmse(:,:,:,end);
[mnrmse, idx] = min(nrmse(:));
[ir, il, is] = ind2sub(size(nrmse),idx);

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
cfigure();
ax = [];

% r/lambdaL
if isnan(vars.lambdaL(1))
    % OptShrink
    subplot(2,2,1);
    if (numel(vars.r) > 1), ax(end + 1) = gca(); end
    plot(vars.r,squeeze(nrmse(:,il,is)),'b-o');
    xlabel('r');
    ylabel('NRMSE');
    title(sprintf('r = %d',vars.r(ir)));
    axis tight; padAxis();
else
    % SVT
    subplot(2,2,1);
    if (numel(vars.lambdaL) > 1), ax(end + 1) = gca(); end
    semilogx(vars.lambdaL,squeeze(nrmse(ir,:,is)),'b-o');
    xlabel('\lambda_L');
    ylabel('NRMSE');
    title(sprintf('\\lambda_L = %.3f',vars.lambdaL(il)));
    axis tight; padAxis();
end

% lambdaS
subplot(2,2,2);
if (numel(vars.lambdaS) > 1), ax(end + 1) = gca(); end
semilogx(vars.lambdaS,squeeze(nrmse(ir,il,:)),'b-o');
xlabel('\lambda_S');
ylabel('NRMSE');
title(sprintf('\\lambda_S = %.4f',vars.lambdaS(is)));
axis tight; padAxis();

% Optimal cost trajectory
subplot(2,2,3);
plot(1:vars.nIters,squeeze(stats.cost(ir,il,is,:)),'b-');
xlabel('Iteration');
ylabel('cost');
title('Optimal cost');
axis tight; padAxis();

% Optimal MSE trajectory
subplot(2,2,4);
plot(1:vars.nIters,squeeze(stats.nrmse(ir,il,is,:)),'b-');
xlabel('Iteration');
ylabel('NRMSE');
title(sprintf('Optimal NRMSE (%.2f%%)',100 * mnrmse));
axis tight; padAxis();

% Match axes
if ~isempty(ax), matchAxes(ax,[],'y'); end

% Save figure
if SAVE_FIGURE
    [~, name] = fileparts(inpath);
    export_fig('-pdf','-transparent',[name, '.pdf']);
    fprintf('DONE\n');
end
%--------------------------------------------------------------------------

%% [PLOT] KTSLR R8

% Knobs
%inpath = 'KTSLR_R8_par1.mat';
%inpath = 'invivo_KTSLR_R8_par1.mat';
%inpath = 'pincat_KTSLR_R8_par1.mat';
%inpath = 'pincat_KTSLR_R8_par1n.mat';
%inpath = 'pincat_KTSLR_R8_par2n.mat';
inpath  = 'pincat_KTSLR_R8_par3n.mat';
SAVE_FIGURE = true;

% Load data (Xhat, vars, stats)
load(inpath);
nIters = vars.nItersO * vars.nItersI;

% Compute optimal parameters
nrmse = stats.nrmse(:,:,:,end);
[mnrmse, idx] = min(nrmse(:));
[ip, i1, i2] = ind2sub(size(nrmse),idx);

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
cfigure([772, 425]);
ax = [];

% p
subplot(2,3,1);
if (numel(vars.p) > 1), ax(end + 1) = gca(); end
plot(vars.p,squeeze(nrmse(:,i1,i2)),'b-o');
xlabel('p');
ylabel('NRMSE');
title(sprintf('p = %.2f',vars.p(ip)));
axis tight; padAxis();

% mu1
subplot(2,3,2);
if (numel(vars.mu1) > 1), ax(end + 1) = gca(); end
semilogx(vars.mu1,squeeze(nrmse(ip,:,i2)),'b-o');
xlabel('\mu_1');
ylabel('NRMSE');
title(sprintf('\\mu_1 = %.4f',vars.mu1(i1)));
axis tight; padAxis();

% mu2
subplot(2,3,3);
if (numel(vars.mu2) > 1), ax(end + 1) = gca(); end
semilogx(vars.mu2,squeeze(nrmse(ip,i1,:)),'b-o');
xlabel('\mu_2');
ylabel('NRMSE');
title(sprintf('\\mu_2 = %.4f',vars.mu2(i2)));
axis tight; padAxis();

% Optimal cost trajectory
subplot(2,3,4);
plot(1:nIters,squeeze(stats.cost(ip,i1,i2,:)),'b-');
xlabel('Iteration');
ylabel('cost');
title('Optimal cost');
axis tight; padAxis();

% Optimal MSE trajectory
subplot(2,3,5);
plot(1:nIters,squeeze(stats.nrmse(ip,i1,i2,:)),'b-');
xlabel('Iteration');
ylabel('NRMSE');
title(sprintf('Optimal NRMSE (%.2f%%)',100 * mnrmse));
axis tight; padAxis();

% Match axes
if ~isempty(ax), matchAxes(ax,[],'y'); end

% Save figure
if SAVE_FIGURE
    [~, name] = fileparts(inpath);
    export_fig('-pdf','-transparent',[name, '.pdf']);
    fprintf('DONE\n');
end
%--------------------------------------------------------------------------

%% [TABLE] sweep otazo

%{
% Knobs
data(1).inpath = './LASSI_full_par1.mat';
data(2).inpath = './LASSI_full_par2.mat';
data(3).inpath = './LASSI_full_par3.mat';
data(4).inpath = './LASSI_full_par4.mat';
data(5).inpath = './RPCA_full_par1.mat';
data(6).inpath = './RPCA_full_par2.mat';
data(7).inpath = './KTSLR_full_par1.mat';
data(1).label = 'LASSI - SVT - hard';
data(2).label = 'LASSI - SVT - soft';
data(3).label = 'LASSI - OPT - hard';
data(4).label = 'LASSI - OPT - soft';
data(5).label = 'RPCA - SVT';
data(6).label = 'RPCA - OPT';
data(7).label = 'KTSLR - TV only';
outpath = 'sweep_otazo.pdf';
sortIdx = 2;
trials  = [1:5];
%}

% Knobs
%data(1).inpath = './LASSI_full_par1.mat'; % SVT, L+S visual
%data(1).inpath = './LASSI_full_par1b.mat'; % SVT, L+S mmse
%data(1).inpath = './LASSI_full_par1e.mat'; % SVT, L+S mmse, L = 0, S = X
data(1).inpath = './LASSI_full_par3.mat'; % OPT, L+S visual
%data(1).inpath = './LASSI_full_par3b.mat'; % OPT, L+S mmse
%data(1).inpath = './LASSI_full_par3c.mat'; % OPT, L+S mmse, L = 0, S = X
data(2).inpath = './LASSI_full_Sonly_par1.mat'; % L+S visual
%data(2).inpath = './LASSI_full_Sonly_par2.mat'; % L+S mmse
%data(3).inpath = './RPCA_full_par1.mat'; % L+S visual
data(3).inpath = './RPCA_full_par1b.mat'; % L+S mmse
data(4).inpath = './KTSLR_full_par2.mat'; % k-t SLR mmse
data(1).label = 'LASSI';
data(2).label = 'DINO-KAT';
data(3).label = 'L+S';
data(4).label = 'k-t SLR';
outpath1 = 'tmi_sweep_otazo1.pdf';
outpath2 = 'tmi_sweep_otazo2.pdf';
sortIdx  = 4;
trials   = [3];

% Load data
nData = numel(data);
cm    = linspecer(nData);
for i = 1:nData
    % Load data
    datai  = load(data(i).inpath);
    
    % Record results (otazo)
    data(i).ps    = datai.vars.ps; %#ok
    data(i).nrmse = nanmean(datai.stats.nrmse(:,trials,end),2); %#ok
end

% Compute effective accelerations
accel = round(1 ./ data(1).ps);
astr  = arrayfun(@(a)sprintf('%dx',a),accel,'UniformOutput',false);

%--------------------------------------------------------------------------
% Print table
%--------------------------------------------------------------------------
% Headings
labels = struct();
labels.align = 'center';
labelFcn     = @(d) sprintf('%dx',d);
labels.strs  = [{'Acceleration'}, arrayfun(labelFcn,accel,'UniformOutput',false)];

% Row data
rows = struct();
for i = 1:nData
    rows(i).name = sprintf('NRMSE (%s) %%',data(i).label);
    rows(i).data = 100 * data(i).nrmse;
end
rows(nData).sep = true;
rows(nData + 1).name = 'Improvement over DINO-KAT (dB)';
rows(nData + 1).data = 20 * log10(data(2).nrmse ./ data(1).nrmse);
rows(nData + 2).name = 'Improvement over L+S (dB)';
rows(nData + 2).data = 20 * log10(data(3).nrmse ./ data(1).nrmse);
rows(nData + 3).name = 'Improvement over k-t SLR (dB)';
rows(nData + 3).data = 20 * log10(data(4).nrmse ./ data(1).nrmse);
[rows.fmt]    = deal('%.1f');
[rows.nAlign] = deal('center');
[rows.dAlign] = deal('right');

% Print table
printTable(labels,rows);
%--------------------------------------------------------------------------

%{
% Sort data, if requested
if exist('sortIdx','var') && ~isempty(sortIdx)
    [~, idx] = sort(cellfun(@(n)n(sortIdx),{data.nrmse}),'descend');
end

%--------------------------------------------------------------------------
% Plot results #1
%--------------------------------------------------------------------------
cfigure([282, 251]);

% Knobs
mk = 'oxd*';
markerSize = 8;
fontSize   = 12;
lineWidth  = 1;

% Plot NRMSE
phndl = zeros(1,nData);
for i = 1:nData
    ii = idx(i);
    phndl(i) = semilogx(data(ii).ps,data(ii).nrmse,'-','Marker',mk(ii),'MarkerSize',markerSize,'LineWidth',lineWidth,'Color',cm(ii,:));
    hold on;
end
xlabel('Column sampling probability');
ylabel('NRMSE');
title('Cardiac perfusion dataset');
legend(phndl,data(idx).label);
axis tight; padAxis();
SetFigFontSize(fontSize);
export_fig('-pdf','-transparent',outpath1);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Plot results #2
%--------------------------------------------------------------------------
cfigure([282, 251]);

% Knobs
mk = 'oxd*';
markerSize = 8;
fontSize   = 12;
lineWidth  = 1;

% Plot NRMSE
phndl = zeros(1,nData);
for i = 1:nData
    ii = idx(i);
    phndl(i) = plot(accel,data(ii).nrmse,'-','Marker',mk(ii),'MarkerSize',markerSize,'LineWidth',lineWidth,'Color',cm(ii,:));
    hold on;
end
xlabel('Acceleration');
ylabel('NRMSE');
title('Cardiac perfusion dataset');
legend(phndl,data(idx).label,'Location','NW');
axis tight; padAxis();
set(gca,'XTick',accel); set(gca,'XTickLabel',astr);
SetFigFontSize(fontSize);
export_fig('-pdf','-transparent',outpath2);
%--------------------------------------------------------------------------
%}

%% [TABLE] sweep invivo

%{
% Knobs
data(1).inpath = './invivo_LASSI_full_par1.mat';
data(2).inpath = './invivo_LASSI_full_par2.mat';
data(3).inpath = './invivo_LASSI_full_par3.mat';
data(4).inpath = './invivo_LASSI_full_par4.mat';
data(5).inpath = './invivo_RPCA_full_par1.mat';
data(6).inpath = './invivo_RPCA_full_par2.mat';
data(7).inpath = './invivo_KTSLR_full_par1.mat';
data(1).label = 'LASSI - SVT - hard';
data(2).label = 'LASSI - SVT - soft';
data(3).label = 'LASSI - OPT - hard';
data(4).label = 'LASSI - OPT - soft';
data(5).label = 'RPCA - SVT';
data(6).label = 'RPCA - OPT';
data(7).label = 'KTSLR - TV only';
outpath = 'sweep_invivo.pdf';
sortIdx = 3;
trials  = [1:5];
%}

% Knobs
data(1).inpath = './invivo_LASSI_full_par1b.mat'; % r = 5, L = 0, S = ktSLR
%data(1).inpath = './invivo_LASSI_full_par1c.mat'; % r = 1, L = 0, S = ktSLR
%data(2).inpath = './invivo_LASSI_full_Sonly_par1.mat'; % r = 5, S = ktSLR
data(2).inpath = './invivo_LASSI_full_Sonly_par1b.mat'; % r = 1, S = ktSLR
data(3).inpath = './invivo_RPCA_full_par1b.mat';
data(4).inpath = './invivo_KTSLR_full_par1.mat';
data(1).label = 'LASSI';
data(2).label = 'DINO-KAT';
data(3).label = 'L+S';
data(4).label = 'k-t SLR';
outpath1 = 'tmi_sweep_invivo1.pdf';
outpath2 = 'tmi_sweep_invivo2.pdf';
sortIdx  = 3;
trials   = [5];

% Load data
nData = numel(data);
cm    = linspecer(nData);
for i = 1:nData
    % Load data
    datai  = load(data(i).inpath);
    
    % Record results (otazo)
    data(i).nLines = datai.vars.nLines; %#ok
    data(i).nrmse  = nanmean(datai.stats.nrmse(:,trials,end),2); %#ok
end

% Compute effective accelerations
%{
nCols = numel(data(1).nLines);
accel = zeros(1,nCols);
for i = 1:nCols
    Y = generateInvivoData(data(1).nLines(i),inf,trials(1),'invivo_full.mat');
    accel(i) = round(numel(Y) / nnz(Y));
end
accel %#ok
%}
accel = [23, 12, 8, 6, 5, 4];

% Flip the order of the data
accel = fliplr(accel);
for i = 1:nData
    data(i).nLines = fliplr(data(i).nLines(:)'); %#ok
    data(i).nrmse  = fliplr(data(i).nrmse(:)'); %#ok
end
astr = arrayfun(@(a)sprintf('%dx',a),accel,'UniformOutput',false);

%--------------------------------------------------------------------------
% Print table
%--------------------------------------------------------------------------
% Headings
labels = struct();
labels.align = 'center';
%labelFcn    = @(d) sprintf('%d',d);
%labels.strs = [{'# Lines'}, arrayfun(labelFcn,data(1).nLines,'UniformOutput',false)];
labelFcn     = @(d) sprintf('%dx',d);
labels.strs  = [{'Acceleration'}, arrayfun(labelFcn,accel,'UniformOutput',false)];

% Row data
rows = struct();
for i = 1:nData
    rows(i).name = sprintf('NRMSE (%s) %%',data(i).label);
    rows(i).data = 100 * data(i).nrmse;
end
rows(nData).sep = true;
rows(nData + 1).name = 'Improvement over DINO-KAT (dB)';
rows(nData + 1).data = 20 * log10(data(2).nrmse ./ data(1).nrmse);
rows(nData + 2).name = 'Improvement over L+S (dB)';
rows(nData + 2).data = 20 * log10(data(3).nrmse ./ data(1).nrmse);
rows(nData + 3).name = 'Improvement over k-t SLR (dB)';
rows(nData + 3).data = 20 * log10(data(4).nrmse ./ data(1).nrmse);
[rows.fmt]    = deal('%.1f');
[rows.nAlign] = deal('center');
[rows.dAlign] = deal('right');

% Print table
printTable(labels,rows);
%--------------------------------------------------------------------------

%{
% Sort data, if requested
if exist('sortIdx','var') && ~isempty(sortIdx)
    [~, idx] = sort(cellfun(@(n)n(sortIdx),{data.nrmse}),'descend');
end

%--------------------------------------------------------------------------
% Plot results #1
%--------------------------------------------------------------------------
cfigure([282, 251]);

% Knobs
mk = 'oxd*';
markerSize = 8;
fontSize   = 12;
lineWidth  = 1;

% Plot NRMSE
phndl = zeros(1,nData);
for i = 1:nData
    ii = idx(i);
    phndl(i) = plot(data(ii).nLines,data(ii).nrmse,'-','Marker',mk(ii),'MarkerSize',markerSize,'LineWidth',lineWidth,'Color',cm(ii,:));
    hold on;
end
xlabel('Number of radial lines');
ylabel('NRMSE');
title('Invivo dataset');
legend(phndl,data(idx).label);
axis tight; padAxis();
set(gca,'XTick',fliplr(data(1).nLines));
SetFigFontSize(fontSize);
export_fig('-pdf','-transparent',outpath1);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Plot results #2
%--------------------------------------------------------------------------
cfigure([282, 251]);

% Knobs
mk = 'oxd*';
markerSize = 8;
fontSize   = 12;
lineWidth  = 1;

% Plot NRMSE
phndl = zeros(1,nData);
for i = 1:nData
    ii = idx(i);
    phndl(i) = semilogx(accel,data(ii).nrmse,'-','Marker',mk(ii),'MarkerSize',markerSize,'LineWidth',lineWidth,'Color',cm(ii,:));
    hold on;
end
xlabel('Acceleration');
ylabel('NRMSE');
title('Invivo dataset');
legend(phndl,data(idx).label,'Location','NW');
axis tight; padAxis();
set(gca,'XTick',accel); set(gca,'XTickLabel',astr);
SetFigFontSize(fontSize);
export_fig('-pdf','-transparent',outpath2);
%--------------------------------------------------------------------------
%}

%% [TABLE] sweep pincat

%{
% Knobs
data(1).inpath = './pincat_LASSI_full_par1.mat';
data(2).inpath = './pincat_LASSI_full_par2.mat';
data(3).inpath = './pincat_LASSI_full_par3.mat';
data(4).inpath = './pincat_LASSI_full_par4.mat';
data(5).inpath = './pincat_RPCA_full_par1.mat';
data(6).inpath = './pincat_RPCA_full_par2.mat';
data(7).inpath = './pincat_KTSLR_full_par1.mat';
data(1).label = 'LASSI - SVT - hard';
data(2).label = 'LASSI - SVT - soft';
data(3).label = 'LASSI - OPT - hard';
data(4).label = 'LASSI - OPT - soft';
data(5).label = 'RPCA - SVT';
data(6).label = 'RPCA - OPT';
data(7).label = 'KTSLR (TV only)';
outpath = 'sweep_pincat.pdf';
sortIdx = 3;
trials  = [1:5];
%}

% Knobs
%data(1).inpath = './pincat_LASSI_full_par1n.mat'; % r = 5, L = 0, S = ktSLR
data(1).inpath = './pincat_LASSI_full_par1nb.mat'; % r = 1, L = 0, S = ktSLR
%data(2).inpath = './pincat_LASSI_full_Sonly_par1n.mat'; % r = 5, S = ktSLR
data(2).inpath = './pincat_LASSI_full_Sonly_par1nb.mat'; % r = 1, S = ktSLR
data(3).inpath = './pincat_RPCA_full_par1nb.mat';
data(4).inpath = './pincat_KTSLR_full_par1nb.mat';
data(1).label = 'LASSI';
data(2).label = 'DINO-KAT';
data(3).label = 'L+S';
data(4).label = 'k-t SLR';
outpath1 = 'tmi_sweep_pincatn1.pdf';
outpath2 = 'tmi_sweep_pincatn2.pdf';
sortIdx  = 3;
trials   = [4];

% Load data
nData = numel(data);
cm    = linspecer(nData);
for i = 1:nData
    % Load data
    datai  = load(data(i).inpath);
    
    % Record results (otazo)
    data(i).nLines = datai.vars.nLines; %#ok
    data(i).nrmse  = nanmean(datai.stats.nrmse(:,trials,end),2); %#ok
end

% Compute effective accelerations
%{
nCols = numel(data(1).nLines);
accel = zeros(1,nCols);
for i = 1:nCols
    Y = generateInvivoData(data(1).nLines(i),45,trials(1),'pincat_full.mat');
    accel(i) = round(numel(Y) / nnz(Y));
end
accel %#ok
%}
accel = [27, 14, 9, 7, 6, 5];

% Flip the order of the data
accel = fliplr(accel);
for i = 1:nData
    data(i).nLines = fliplr(data(i).nLines(:)'); %#ok
    data(i).nrmse  = fliplr(data(i).nrmse(:)'); %#ok
end
astr = arrayfun(@(a)sprintf('%dx',a),accel,'UniformOutput',false);

%--------------------------------------------------------------------------
% Print table
%--------------------------------------------------------------------------
% Headings
labels = struct();
labels.align = 'center';
%labelFcn    = @(d) sprintf('%d',d);
%labels.strs = [{'# Lines'}, arrayfun(labelFcn,data(1).nLines,'UniformOutput',false)];
labelFcn     = @(d) sprintf('%dx',d);
labels.strs  = [{'Acceleration'}, arrayfun(labelFcn,accel,'UniformOutput',false)];

% Row data
rows = struct();
for i = 1:nData
    rows(i).name = sprintf('NRMSE (%s) %%',data(i).label);
    rows(i).data = 100 * data(i).nrmse;
end
rows(nData).sep = true;
rows(nData + 1).name = 'Improvement over DINO-KAT (dB)';
rows(nData + 1).data = 20 * log10(data(2).nrmse ./ data(1).nrmse);
rows(nData + 2).name = 'Improvement over L+S (dB)';
rows(nData + 2).data = 20 * log10(data(3).nrmse ./ data(1).nrmse);
rows(nData + 3).name = 'Improvement over k-t SLR (dB)';
rows(nData + 3).data = 20 * log10(data(4).nrmse ./ data(1).nrmse);
[rows.fmt]    = deal('%.1f');
[rows.nAlign] = deal('center');
[rows.dAlign] = deal('right');

% Print table
printTable(labels,rows);
%--------------------------------------------------------------------------

%{
% Sort data, if requested
if exist('sortIdx','var') && ~isempty(sortIdx)
    [~, idx] = sort(cellfun(@(n)n(sortIdx),{data.nrmse}),'descend');
end

%--------------------------------------------------------------------------
% Plot results #1
%--------------------------------------------------------------------------
cfigure([282, 251]);

% Knobs
mk = 'oxd*';
markerSize = 8;
fontSize   = 12;
lineWidth  = 1;

% Plot NRMSE
phndl = zeros(1,nData);
for i = 1:nData
    ii = idx(i);
    phndl(i) = plot(data(ii).nLines,data(ii).nrmse,'-','Marker',mk(ii),'MarkerSize',markerSize,'LineWidth',lineWidth,'Color',cm(ii,:));
    hold on;
end
xlabel('Number of radial lines');
ylabel('NRMSE');
title('PINCAT dataset');
legend(phndl,data(idx).label);
axis tight; padAxis();
set(gca,'XTick',fliplr(data(1).nLines));
SetFigFontSize(fontSize);
export_fig('-pdf','-transparent',outpath1);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Plot results #2
%--------------------------------------------------------------------------
cfigure([282, 251]);

% Knobs
mk = 'oxd*';
markerSize = 8;
fontSize   = 12;
lineWidth  = 1;

% Plot NRMSE
phndl = zeros(1,nData);
for i = 1:nData
    ii = idx(i);
    phndl(i) = semilogx(accel,data(ii).nrmse,'-','Marker',mk(ii),'MarkerSize',markerSize,'LineWidth',lineWidth,'Color',cm(ii,:));
    hold on;
end
xlabel('Acceleration');
ylabel('NRMSE');
title('PINCAT dataset');
legend(phndl,data(idx).label,'Location','NW');
axis tight; padAxis();
set(gca,'XTick',accel); set(gca,'XTickLabel',astr);
SetFigFontSize(fontSize);
export_fig('-pdf','-transparent',outpath2);
%--------------------------------------------------------------------------
%}

%% [TABLE] sweep otazo (ROI)

% Knobs
%data(1).inpath = 'LASSI_full_par1'; % SVT, L+S visual
%data(1).inpath = 'LASSI_full_par1b'; % SVT, L+S mmse
%data(1).inpath = 'LASSI_full_par1e'; % SVT, L+S mmse, L = 0, S = X
data(1).inpath = 'LASSI_full_par3'; % OPT, L+S visual
%data(1).inpath = 'LASSI_full_par3b'; % OPT, L+S mmse
%data(1).inpath = 'LASSI_full_par3c'; % OPT, L+S mmse, L = 0, S = X
data(2).inpath = 'LASSI_full_Sonly_par1'; % L+S visual
%data(2).inpath = 'LASSI_full_Sonly_par2'; % L+S mmse
%data(3).inpath = 'RPCA_full_par1.mat'; % L+S visual
data(3).inpath = 'RPCA_full_par1b'; % L+S mmse
data(4).inpath = 'KTSLR_full_par2'; % k-t SLR mmse
data(1).label = 'LASSI';
data(2).label = 'DINO-KAT';
data(3).label = 'L+S';
data(4).label = 'k-t SLR';
basedir  = '/Users/Brian/Desktop/lassi_journal';
inds     = 13:18;
truepath = 'otazo_full.mat';
roipath  = 'otazo_roi.mat';

% Load data
datat = load(truepath);
Xtrue = datat.Xtrue;
[ny, nx, nt] = deal(datat.ny,datat.nx,datat.nt);
datar = load(roipath,'ROI');
ROI   = datar.ROI;
nData = numel(data);
nP    = numel(inds);
for i = 1:nData
    for j = 1:nP
        % Load data
        path  = sprintf('%s/%s/data%d.mat',basedir,data(i).inpath,inds(j));
        datai = load(path);
        
        % Compute NRMSE
        if isfield(datai,'Lhat')
            Xhat = datai.Lhat + datai.Shat;
        else
            Xhat = datai.Xhat;
        end
        Xhat = reshape(Xhat,ny,nx,nt);
        nrmse = computeNRMSE(Xhat,Xtrue,ROI);
        
        % Record results
        data(i).ps(j)    = datai.params.ps; %#ok
        data(i).nrmse(j) = nrmse; %#ok
    end
end

% Compute effective accelerations
accel = round(1 ./ data(1).ps);

%--------------------------------------------------------------------------
% Print table
%--------------------------------------------------------------------------
% Headings
labels = struct();
labels.align = 'center';
labelFcn     = @(d) sprintf('%dx',d);
labels.strs  = [{'Acceleration'}, arrayfun(labelFcn,accel,'UniformOutput',false)];

% Row data
rows = struct();
for i = 1:nData
    rows(i).name = sprintf('NRMSE (%s) %%',data(i).label);
    rows(i).data = 100 * data(i).nrmse;
end
rows(nData).sep = true;
rows(nData + 1).name = 'Improvement over DINO-KAT (dB)';
rows(nData + 1).data = 20 * log10(data(2).nrmse ./ data(1).nrmse);
rows(nData + 2).name = 'Improvement over L+S (dB)';
rows(nData + 2).data = 20 * log10(data(3).nrmse ./ data(1).nrmse);
rows(nData + 3).name = 'Improvement over k-t SLR (dB)';
rows(nData + 3).data = 20 * log10(data(4).nrmse ./ data(1).nrmse);
[rows.fmt]    = deal('%.1f');
[rows.nAlign] = deal('center');
[rows.dAlign] = deal('right');

% Print table
printTable(labels,rows);
%--------------------------------------------------------------------------

%% [TABLE] sweep invivo (ROI)

% Knobs
data(1).inpath = 'invivo_LASSI_full_par1b'; % r = 5, L = 0, S = ktSLR
%data(1).inpath = 'invivo_LASSI_full_par1c'; % r = 1, L = 0, S = ktSLR
%data(2).inpath = 'invivo_LASSI_full_Sonly_par1'; % r = 5, S = ktSLR
data(2).inpath = 'invivo_LASSI_full_Sonly_par1b'; % r = 1, S = ktSLR
data(3).inpath = 'invivo_RPCA_full_par1b';
data(4).inpath = 'invivo_KTSLR_full_par1';
data(1).label = 'LASSI';
data(2).label = 'DINO-KAT';
data(3).label = 'L+S';
data(4).label = 'k-t SLR';
basedir  = '/Users/Brian/Desktop/lassi_journal';
inds     = 25:30;
truepath = 'invivo_full.mat';
roipath  = 'invivo_roi.mat';

% Load data
datat = load(truepath);
Xtrue = datat.Xtrue;
[ny, nx, nt] = deal(datat.ny,datat.nx,datat.nt);
datar = load(roipath,'ROI');
ROI   = datar.ROI;
nData = numel(data);
nP    = numel(inds);
for i = 1:nData
    for j = 1:nP
        % Load data
        path  = sprintf('%s/%s/data%d.mat',basedir,data(i).inpath,inds(j));
        datai = load(path);
        
        % Compute NRMSE
        if isfield(datai,'Lhat')
            Xhat = datai.Lhat + datai.Shat;
        else
            Xhat = datai.Xhat;
        end
        Xhat = reshape(Xhat,ny,nx,nt);
        nrmse = computeNRMSE(Xhat,Xtrue,ROI);
        
        % Record results
        data(i).nLines(j) = datai.params.nLines; %#ok
        data(i).nrmse(j)  = nrmse; %#ok
    end
end

% Effective accelerations
accel = [23, 12, 8, 6, 5, 4];

% Flip the order of the data
accel = fliplr(accel);
for i = 1:nData
    data(i).nLines = fliplr(data(i).nLines(:)'); %#ok
    data(i).nrmse  = fliplr(data(i).nrmse(:)'); %#ok
end

%--------------------------------------------------------------------------
% Print table
%--------------------------------------------------------------------------
% Headings
labels = struct();
labels.align = 'center';
labelFcn     = @(d) sprintf('%dx',d);
labels.strs  = [{'Acceleration'}, arrayfun(labelFcn,accel,'UniformOutput',false)];

% Row data
rows = struct();
for i = 1:nData
    rows(i).name = sprintf('NRMSE (%s) %%',data(i).label);
    rows(i).data = 100 * data(i).nrmse;
end
rows(nData).sep = true;
rows(nData + 1).name = 'Improvement over DINO-KAT (dB)';
rows(nData + 1).data = 20 * log10(data(2).nrmse ./ data(1).nrmse);
rows(nData + 2).name = 'Improvement over L+S (dB)';
rows(nData + 2).data = 20 * log10(data(3).nrmse ./ data(1).nrmse);
rows(nData + 3).name = 'Improvement over k-t SLR (dB)';
rows(nData + 3).data = 20 * log10(data(4).nrmse ./ data(1).nrmse);
[rows.fmt]    = deal('%.1f');
[rows.nAlign] = deal('center');
[rows.dAlign] = deal('right');

% Print table
printTable(labels,rows);
%--------------------------------------------------------------------------

%% [TABLE] sweep pincat (ROI)

% Knobs
%data(1).inpath = 'pincat_LASSI_full_par1n'; % r = 5, L = 0, S = ktSLR
data(1).inpath = 'pincat_LASSI_full_par1nb'; % r = 1, L = 0, S = ktSLR
%data(2).inpath = 'pincat_LASSI_full_Sonly_par1n'; % r = 5, S = ktSLR
data(2).inpath = 'pincat_LASSI_full_Sonly_par1nb'; % r = 1, S = ktSLR
data(3).inpath = 'pincat_RPCA_full_par1nb';
data(4).inpath = 'pincat_KTSLR_full_par1nb';
data(1).label = 'LASSI';
data(2).label = 'DINO-KAT';
data(3).label = 'L+S';
data(4).label = 'k-t SLR';
basedir  = '/Users/Brian/Desktop/lassi_journal';
inds     = 19:24;
truepath = 'pincat_full.mat';
roipath  = 'pincat_roi.mat';

% Load data
datat = load(truepath);
Xtrue = datat.Xtrue;
[ny, nx, nt] = deal(datat.ny,datat.nx,datat.nt);
datar = load(roipath,'ROI');
ROI   = datar.ROI;
nData = numel(data);
nP    = numel(inds);
for i = 1:nData
    for j = 1:nP
        % Load data
        path  = sprintf('%s/%s/data%d.mat',basedir,data(i).inpath,inds(j));
        datai = load(path);
        
        % Compute NRMSE
        if isfield(datai,'Lhat')
            Xhat = datai.Lhat + datai.Shat;
        else
            Xhat = datai.Xhat;
        end
        Xhat = reshape(Xhat,ny,nx,nt);
        nrmse = computeNRMSE(Xhat,Xtrue,ROI);
        
        % Record results
        data(i).nLines(j) = datai.params.nLines; %#ok
        data(i).nrmse(j)  = nrmse; %#ok
    end
end

% Effective accelerations
accel = [27, 14, 9, 7, 6, 5];

% Flip the order of the data
accel = fliplr(accel);
for i = 1:nData
    data(i).nLines = fliplr(data(i).nLines(:)'); %#ok
    data(i).nrmse  = fliplr(data(i).nrmse(:)'); %#ok
end

%--------------------------------------------------------------------------
% Print table
%--------------------------------------------------------------------------
% Headings
labels = struct();
labels.align = 'center';
labelFcn     = @(d) sprintf('%dx',d);
labels.strs  = [{'Acceleration'}, arrayfun(labelFcn,accel,'UniformOutput',false)];

% Row data
rows = struct();
for i = 1:nData
    rows(i).name = sprintf('NRMSE (%s) %%',data(i).label);
    rows(i).data = 100 * data(i).nrmse;
end
rows(nData).sep = true;
rows(nData + 1).name = 'Improvement over DINO-KAT (dB)';
rows(nData + 1).data = 20 * log10(data(2).nrmse ./ data(1).nrmse);
rows(nData + 2).name = 'Improvement over L+S (dB)';
rows(nData + 2).data = 20 * log10(data(3).nrmse ./ data(1).nrmse);
rows(nData + 3).name = 'Improvement over k-t SLR (dB)';
rows(nData + 3).data = 20 * log10(data(4).nrmse ./ data(1).nrmse);
[rows.fmt]    = deal('%.1f');
[rows.nAlign] = deal('center');
[rows.dAlign] = deal('right');

% Print table
printTable(labels,rows);
%--------------------------------------------------------------------------

%% [PLOT] Dinit experiments

%{
% Knobs
% L0 = L, S0 = S from L+S mmse
data(1).inpath = './LASSI_full_par1b.mat';
data(2).inpath = './LASSI_full_par1c.mat';
data(3).inpath = './LASSI_full_par1d.mat';
data(4).inpath = './LASSI_full_par1b2.mat';
data(5).inpath = './LASSI_full_par1c2.mat';
data(6).inpath = './LASSI_full_par1d2.mat';
data(1).label = 'Full DCT';
data(2).label = 'Separable DCT';
data(3).label = 'Random';
data(4).label = 'Full DCT - FIXED';
data(5).label = 'Separable DCT - FIXED';
data(6).label = 'Random - FIXED';
outpath = 'tmi_Dinit.pdf';
sortIdx = 2;
trials  = [1:4];
%}

%{
% Knobs
% L0 = L, S0 = S from L+S visual
data(1).inpath = './LASSI_full_par1b3.mat';
data(2).inpath = './LASSI_full_par1c3.mat';
data(3).inpath = './LASSI_full_par1d3.mat';
data(4).inpath = './LASSI_full_par1b4.mat';
data(5).inpath = './LASSI_full_par1c4.mat';
data(6).inpath = './LASSI_full_par1d4.mat';
data(1).label = 'Full DCT';
data(2).label = 'Separable DCT';
data(3).label = 'Random';
data(4).label = 'Full DCT - FIXED';
data(5).label = 'Separable DCT - FIXED';
data(6).label = 'Random - FIXED';
outpath = 'tmi_Dinit.pdf';
sortIdx = 2;
trials  = [1:5];
%}

% Knobs
% L0 = 0, S0 = L+S from L+S mmse
data(1).inpath = './LASSI_full_par1d5.mat';
%data(1).inpath = './LASSI_full_Drand_par1.mat'; % 5 D0's, trial = 3 only
data(2).inpath = './LASSI_full_par1b5.mat';
data(3).inpath = './LASSI_full_par1c5.mat';
data(4).inpath = './LASSI_full_par1c6.mat';
data(5).inpath = './LASSI_full_par1b6.mat';
data(6).inpath = './LASSI_full_par1d6.mat';
%data(6).inpath = './LASSI_full_Drand_par2.mat'; % 5 D0's, trial = 3 only
data(1).label = 'Random (init.)';
data(2).label = '1-D DCT (init.)';
data(3).label = 'Separable DCT (init.)';
data(4).label = 'Separable DCT (fixed)';
data(5).label = '1-D DCT (fixed)';
data(6).label = 'Random (fixed)';
%outpath = 'tmi_Dinit.pdf';
sortIdx = 6;
trials  = [3]; % Must be 3 if using Drand sims
%D0s    = [5];

% Load data
nData = numel(data);
for i = 1:nData
    % Load data
    datai  = load(data(i).inpath);
    
    % Record results (otazo)
    data(i).ps    = datai.vars.ps; %#ok
    if ndims(datai.stats.nrmse) == 4
        % Drand data
        data(i).nrmse = nanmean(datai.stats.nrmse(:,1,D0s,end),3); %#ok
    else
        % Regular data
        data(i).nrmse = nanmean(datai.stats.nrmse(:,trials,end),2); %#ok
    end
end
accel = round(1 ./ data(1).ps);
astr  = arrayfun(@(a)sprintf('%dx',a),accel,'UniformOutput',false);

% Sort data, if requested
if exist('sortIdx','var') && ~isempty(sortIdx)
    [~, idx] = sort(cellfun(@(n)n(sortIdx),{data.nrmse}),'descend');
end

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
%cfigure([419, 349]);
cfigure([332, 272]);

% Knobs
if nData > 4
    % Remove yellow
    cm = linspecer(nData + 1);
    cm = cm([1:4,6:end],:);
else
    cm = linspecer(nData);
end
if nData > 5
    % Switch blue and red
    cm = [cm(2,:); cm(1,:); cm(3:end,:)];
end
mk = 'o*v^d<>';
markerSize = 8;
fontSize   = 12;
lineWidth  = 1;

% Plot NRMSE
phndl = zeros(1,nData);
for i = 1:nData
    ii = idx(i);
    phndl(i) = plot(accel,data(ii).nrmse,'-','Marker',mk(ii),'MarkerSize',markerSize,'LineWidth',lineWidth,'Color',cm(ii,:));
    hold on;
end
xlabel('Acceleration');
ylabel('NRMSE');
title('LASSI dictionary initialization');
legend(phndl,data(idx).label,'Location','NW');
axis tight; padAxis();
set(gca,'XTick',accel); set(gca,'XTickLabel',astr);
ylim = get(gca,'YLim'); set(gca,'YLim',[ylim(1), (ylim(2) + 0.05)]); % HACK
SetFigFontSize(fontSize);
export_fig -pdf -transparent tmi_dict_init
%--------------------------------------------------------------------------

%% [PLOT] LSinit experiments

% Knobs
data(1).inpath = './LASSI_full_par1e.mat'; % L+S init
data(2).inpath = './LASSI_full_par1f.mat'; % ktSLR init
%data(3).inpath = './LASSI_full_par1g.mat'; % FFT init
data(3).inpath = './LASSI_full_par1g2.mat'; % FFT init (larger lambdaB)
data(4).inpath = './RPCA_full_par1b.mat'; % L+S mmse
data(5).inpath = './KTSLR_full_par2.mat'; % ktSLR mmse
data(6).inpath = './FFT_full_par.mat'; % FFT
data(1).label = 'LASSI (L+S init.)';
data(2).label = 'LASSI (k-t SLR init.)';
data(3).label = 'LASSI (baseline init.)';
data(4).label = 'L+S';
data(5).label = 'k-t SLR';
data(6).label = 'Baseline recon.';
trials = [3];

% Load data
nData = numel(data);
for i = 1:nData
    % Load data
    datai  = load(data(i).inpath);
    
    % Record results (otazo)
    data(i).ps    = datai.vars.ps; %#ok
    data(i).nrmse = nanmean(datai.stats.nrmse(:,trials,end),2); %#ok
end
accel = round(1 ./ data(1).ps);
astr  = arrayfun(@(a)sprintf('%dx',a),accel,'UniformOutput',false);

% Display NRMSEs
nrmse_lps_ktslr_fft = [data(1).nrmse(:), ...
                       data(2).nrmse(:), ...
                       data(3).nrmse(:)] %#ok

%{
%--------------------------------------------------------------------------
% Plot results #1
%--------------------------------------------------------------------------
% Knobs
%cm = linspecer(nData);
cm = linspecer(9); cm = cm([1:4, 6, 9],:);
mk = 'o*v^d<>';
markerSize = 8;
fontSize   = 12;
lineWidth  = 1;

% L+S initialization
cfigure([332, 272]);
phndl = zeros(1,2);
phndl(1) = plot(accel,data(4).nrmse,'-','Marker',mk(3),'MarkerSize',markerSize,'LineWidth',lineWidth,'Color',cm(3,:)); hold on;
phndl(2) = plot(accel,data(1).nrmse,'-','Marker',mk(1),'MarkerSize',markerSize,'LineWidth',lineWidth,'Color',cm(2,:)); hold on;
xlabel('Acceleration');
ylabel('NRMSE');
title('L+S initialization');
legend(phndl,data([4, 1]).label,'Location','NW');
axis tight; padAxis();
set(gca,'XTick',accel); set(gca,'XTickLabel',astr);
set(gca,'YLim',[0.0940, 0.2785]); % HACK
SetFigFontSize(fontSize);
export_fig -pdf -transparent tmi_LSinit1

% k-t SLR initialization
cfigure([332, 272]);
phndl = zeros(1,2);
phndl(1) = plot(accel,data(5).nrmse,'-','Marker',mk(4),'MarkerSize',markerSize,'LineWidth',lineWidth,'Color',cm(4,:)); hold on;
phndl(2) = plot(accel,data(2).nrmse,'-','Marker',mk(2),'MarkerSize',markerSize,'LineWidth',lineWidth,'Color',cm(1,:)); hold on;
xlabel('Acceleration');
ylabel('NRMSE');
title('k-t SLR initialization');
legend(phndl,data([5, 2]).label,'Location','NW');
axis tight; padAxis();
set(gca,'XTick',accel); set(gca,'XTickLabel',astr);
set(gca,'YLim',[0.0940, 0.2785]); % HACK
SetFigFontSize(fontSize);
export_fig -pdf -transparent tmi_LSinit2

% FFT initialization
cfigure([332, 272]);
phndl = zeros(1,2);
phndl(1) = plot(accel,data(6).nrmse,'-','Marker',mk(5),'MarkerSize',markerSize,'LineWidth',lineWidth,'Color',cm(5,:)); hold on;
phndl(2) = plot(accel,data(3).nrmse,'-','Marker',mk(6),'MarkerSize',markerSize,'LineWidth',lineWidth,'Color',cm(6,:)); hold on;
xlabel('Acceleration');
ylabel('NRMSE');
title('Baseline initialization');
legend(phndl,data([6, 3]).label,'Location','Best');
axis tight; padAxis();
set(gca,'XTick',accel); set(gca,'XTickLabel',astr);
%set(gca,'YLim',[0.0940, 0.2785]); % HACK
SetFigFontSize(fontSize);
export_fig -pdf -transparent tmi_LSinit3
%--------------------------------------------------------------------------
%}

%--------------------------------------------------------------------------
% Plot results #2
%--------------------------------------------------------------------------
% Knobs
%cm = linspecer(nData);
cm = linspecer(9); cm = cm([1:4, 6, 9],:);
mk = 'o*v^d<>';
markerSize = 8;
fontSize   = 12;
lineWidth  = 1;

% L+S initialization
cfigure([332, 272]);
phndl = zeros(1,2);
phndl(1) = plot(accel,data(4).nrmse,'-','Marker',mk(3),'MarkerSize',markerSize,'LineWidth',lineWidth,'Color',cm(3,:)); hold on;
phndl(2) = plot(accel,data(1).nrmse,'-','Marker',mk(1),'MarkerSize',markerSize,'LineWidth',lineWidth,'Color',cm(2,:)); hold on;
xlabel('Acceleration');
ylabel('NRMSE');
title('L+S initialization');
legend(phndl,data([4, 1]).label,'Location','NW');
axis tight; padAxis();
set(gca,'XTick',accel); set(gca,'XTickLabel',astr);
set(gca,'YLim',[0.0940, 0.2785]); % HACK
SetFigFontSize(fontSize);
export_fig -pdf -transparent tmi_LSinit1

% k-t SLR initialization
cfigure([332, 272]);
phndl = zeros(1,3);
phndl(1) = plot(accel,data(5).nrmse,'-','Marker',mk(4),'MarkerSize',markerSize,'LineWidth',lineWidth,'Color',cm(4,:)); hold on;
phndl(2) = plot(accel,data(2).nrmse,'-','Marker',mk(2),'MarkerSize',markerSize,'LineWidth',lineWidth,'Color',cm(1,:)); hold on;
phndl(3) = plot(accel,data(3).nrmse,'-','Marker',mk(6),'MarkerSize',markerSize,'LineWidth',lineWidth,'Color',cm(6,:)); hold on;
xlabel('Acceleration');
ylabel('NRMSE');
title('k-t SLR initialization');
legend(phndl,data([5, 2, 3]).label,'Location','NW');
axis tight; padAxis();
set(gca,'XTick',accel); set(gca,'XTickLabel',astr);
set(gca,'YLim',[0.0940, 0.2785]); % HACK
SetFigFontSize(fontSize);
export_fig -pdf -transparent tmi_LSinit2
%--------------------------------------------------------------------------

%% [PLOT] LASSI variations experiments

% Knobs
code = [1, 1, 1, 1, 1, 1, 1, 1]; % mmse
outpath1 = 'tmi_lassi_vars1_mmse.pdf';
outpath2 = 'tmi_lassi_vars2_mmse.pdf';
%{
code = [2, 2, 2, 2, 2, 2, 2, 2]; % visual
outpath1 = 'tmi_lassi_vars1_visual.pdf';
outpath2 = 'tmi_lassi_vars2_visual.pdf';
%}
%{
code = [1, 1, -1, -1, 1, 1, 1, 1]; % best
outpath1 = 'tmi_lassi_vars1_best.pdf';
outpath2 = 'tmi_lassi_vars2_best.pdf';
%}

% L+S mmse
data1(1).inpath = './LASSI_full_par1e.mat';  % hard, L = 0, S = L+S mmse
data1(2).inpath = './LASSI_full_par2c.mat';  % soft, L = 0, S = L+S mmse
data1(3).inpath = './LASSI_full_par3c.mat';  % hard, L = 0, S = L+S mmse
data1(4).inpath = './LASSI_full_par4c.mat';  % soft, L = 0, S = L+S mmse
data1(5).inpath = './LASSI_full_par1h2.mat'; % hard, p = 0  , L = 0, S = L+S mmse
data1(6).inpath = './LASSI_full_par2d2.mat'; % soft, p = 0  , L = 0, S = L+S mmse
data1(7).inpath = './LASSI_full_par1i2.mat'; % hard, p = 0.5, L = 0, S = L+S mmse
data1(8).inpath = './LASSI_full_par2e2.mat'; % soft, p = 0.5, L = 0, S = L+S mmse

% L+S visual
data2(1).inpath = './LASSI_full_par1.mat';   % hard, L+S visual
data2(2).inpath = './LASSI_full_par2.mat';   % soft, L+S visual
data2(3).inpath = './LASSI_full_par3.mat';   % hard, L+S visual
data2(4).inpath = './LASSI_full_par4.mat';   % soft, L+S visual
data2(5).inpath = './LASSI_full_par1h3.mat'; % hard, p = 0  , L+S visual
data2(6).inpath = './LASSI_full_par2d3.mat'; % soft, p = 0  , L+S visual
data2(7).inpath = './LASSI_full_par1i3.mat'; % hard, p = 0.5, L+S visual
data2(8).inpath = './LASSI_full_par2e3.mat'; % soft, p = 0.5, L+S visual

% Knobs (universal)
data(1).label = 'LASSI (p = 1)'; % L0
data(2).label = 'LASSI (p = 1)'; % L1
data(3).label = 'LASSI (OPT)'; % L0
data(4).label = 'LASSI (OPT)'; % L1
data(5).label = 'LASSI (p = 0)'; % L0
data(6).label = 'LASSI (p = 0)'; % L1
data(7).label = 'LASSI (p = 0.5)'; % L0
data(8).label = 'LASSI (p = 0.5)'; % L1
trials  = [3];

% Load data
nData = numel(data);
for i = 1:nData
    if code(i) == 1
        % mmse
        datai = load(data1(i).inpath);
        
    elseif code(i) == 2
        % visual
        datai = load(data2(i).inpath);
    else
        % pointwise best
        data1i = load(data1(i).inpath);
        data2i = load(data2(i).inpath);
        datai.vars.ps = data1i.vars.ps;
        datai.stats.nrmse = min(data1i.stats.nrmse,data2i.stats.nrmse);
    end
    data(i).ps    = datai.vars.ps; %#ok
    data(i).nrmse = nanmean(datai.stats.nrmse(:,trials,end),2); %#ok
end

accel = round(1 ./ data(1).ps);
astr  = arrayfun(@(a)sprintf('%dx',a),accel,'UniformOutput',false);

% Print NRMSEs
nrmse_L0_p1_opt_p0_p0p5 = [data(1).nrmse(:), data(3).nrmse(:), data(5).nrmse(:), data(7).nrmse(:)] %#ok
nrmse_L1_p1_opt_p0_p0p5 = [data(2).nrmse(:), data(4).nrmse(:), data(6).nrmse(:), data(8).nrmse(:)] %#ok

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
% Knobs
%cm = linspecer(nData);
cm = linspecer(9); cm(5,:) = [];
mk = 'o*v^d<>s+x';
markerSize = 8;
fontSize   = 12;
lineWidth  = 1;

p1  = [5, 1, 7, 3]; np1 = numel(p1);
p2  = [6, 2, 8, 4]; np2 = numel(p2);

% L0 regularization
cfigure([332, 272]);
phndl = zeros(1,np1);
for i = 1:np1
    ii = p1(i);
    phndl(i) = plot(accel,data(ii).nrmse,'-','Marker',mk(ii),'MarkerSize',markerSize,'LineWidth',lineWidth,'Color',cm(ii,:));
    hold on;
end
xlabel('Acceleration');
ylabel('NRMSE');
title('L0 sparsity');
legend(phndl,data(p1).label,'Location','NW');
axis tight; ax1 = gca();
set(ax1,'XTick',accel); set(ax1,'XTickLabel',astr);
SetFigFontSize(fontSize);

% L1 regularization
cfigure([332, 272]);
phndl = zeros(1,np2);
for i = 1:np2
    ii = p2(i);
    phndl(i) = plot(accel,data(ii).nrmse,'-','Marker',mk(ii),'MarkerSize',markerSize,'LineWidth',lineWidth,'Color',cm(ii,:));
    hold on;
end
xlabel('Acceleration');
ylabel('NRMSE');
title('L1 sparsity');
legend(phndl,data(p2).label,'Location','NW');
axis tight; ax2 = gca();
set(ax2,'XTick',accel); set(ax2,'XTickLabel',astr);
SetFigFontSize(fontSize);

matchAxes([ax1, ax2]);
export_fig('-pdf','-transparent',outpath1,ax1);
export_fig('-pdf','-transparent',outpath2,ax2);
%--------------------------------------------------------------------------

%% [PLOT] LASSI atom rank experiments

% Knobs
data(1).inpath = './LASSI_full_par1e.mat';  % r = 1, L+S mmse, L = 0, S = X
data(2).inpath = './LASSI_full_par1e2.mat'; % r = 2, L+S mmse, L = 0, S = X
data(3).inpath = './LASSI_full_par1e3.mat'; % r = 3, L+S mmse, L = 0, S = X
data(4).inpath = './LASSI_full_par1e4.mat'; % r = 4, L+S mmse, L = 0, S = X
data(5).inpath = './LASSI_full_par1e5.mat'; % r = 5, L+S mmse, L = 0, S = X
trials  = [3];
outpath = './tmi_lowrank_atoms.pdf';

%{
% Knobs
data(1).inpath = './LASSI_full_Sonly_par2.mat';  % r = 1, L+S mmse
data(2).inpath = './LASSI_full_Sonly_par2b.mat'; % r = 2, L+S mmse
data(3).inpath = './LASSI_full_Sonly_par2c.mat'; % r = 3, L+S mmse
data(4).inpath = './LASSI_full_Sonly_par2d.mat'; % r = 4, L+S mmse
data(5).inpath = './LASSI_full_Sonly_par2e.mat'; % r = 5, L+S mmse
trials  = [3];
outpath = './tmi_lowrank_atoms_dinokat.pdf';
%}

% Load data
nData = numel(data);
for i = 1:nData
    % Load data
    datai  = load(data(i).inpath);
    
    % Record results (otazo)
    data(i).ps    = datai.vars.ps; %#ok
    data(i).nrmse = nanmean(datai.stats.nrmse(:,trials,end),2); %#ok
end

% Compute gains
gainDB = 20 * log10(data(5).nrmse ./ data(1).nrmse);
accel  = round(1 ./ data(1).ps);
astr   = arrayfun(@(a)sprintf('%dx',a),accel,'UniformOutput',false);

%--------------------------------------------------------------------------
% Plot gains
%--------------------------------------------------------------------------
%cfigure([287, 241]);
cfigure([332, 272]);

% Knobs
cm       = linspecer(1);
fontSize = 12;

bar(accel,gainDB,'FaceColor',cm);
xlabel('Acceleration');
ylabel('Gain over r = 5 (dB)');
title('LASSI: rank-1 atoms');
%title('DINO-KAT: rank-1 atoms');
axis tight;
set(gca,'YLim',[0, 0.25]);
padAxis([],[],'x');
set(gca,'XTick',accel); set(gca,'XTickLabel',astr);
SetFigFontSize(fontSize);
export_fig('-pdf','-transparent',outpath); 
%--------------------------------------------------------------------------

%% [PLOT] LASSI overcompleteness experiments

% Knobs
data(1).inpath = './LASSI_full_par1eb.mat'; % np = -200, r = 1, L+S mmse, L = 0, S = X
data(2).inpath = './LASSI_full_par1ec.mat'; % np = -100, r = 1, L+S mmse, L = 0, S = X
data(3).inpath = './LASSI_full_par1e.mat';  % np =    0, r = 1, L+S mmse, L = 0, S = X
data(4).inpath = './LASSI_full_par1ed.mat'; % np =  100, r = 1, L+S mmse, L = 0, S = X
data(5).inpath = './LASSI_full_par1ee.mat'; % np =  200, r = 1, L+S mmse, L = 0, S = X
data(1).label = '120 atoms';
data(2).label = '220 atoms';
data(3).label = '320 atoms';
data(4).label = '420 atoms';
data(5).label = '520 atoms';
sortIdx = 3;
trials  = [3]; % Must be trials = 3
outpath = 'tmi_LASSIover.pdf';

% Load data
nData = numel(data);
for i = 1:nData
    % Load data
    datai  = load(data(i).inpath);
    
    % Record results (otazo)
    data(i).ps    = datai.vars.ps; %#ok
    data(i).nrmse = nanmean(datai.stats.nrmse(:,trials,end),2); %#ok
end

% Knobs
inds   = [2:5];

% Extract interesting data
accel  = round(1 ./ data(1).ps(inds));
atoms  = [120, 220, 320, 420, 520];
nAccel = numel(accel);
nAtoms = numel(atoms);
nrmse  = zeros(nAccel,nData);
for i = 1:nData
    nrmse(:,i) = data(i).nrmse(inds);
end
[~, idx] = sort(nrmse(:,1),'descend');

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
%cfigure([326, 269]);
cfigure([332, 272]);

% Knobs
if nAccel > 4
    % Remove yellow
    cm = linspecer(nAccel + 1);
    cm = cm([1:4,6:end],:);
else
    cm = linspecer(nAccel);
end
if nData > 5
    % Switch blue and red
    cm = [cm(2,:); cm(1,:); cm(3:end,:)];
end
mk = 'o*v^d<>';
markerSize = 8;
fontSize   = 12;
lineWidth  = 1;

% Plot results
phndl = zeros(1,nAccel);
pstr  = cell(1,nAccel);
for i = 1:nAccel
    ii = idx(i);
    phndl(i) = plot(atoms,nrmse(ii,:),'-','Marker',mk(ii),'MarkerSize',markerSize,'LineWidth',lineWidth,'Color',cm(ii,:));
    pstr{i}  = sprintf('%dx accel.',accel(ii));
    hold on;
end
xlabel('Number of atoms');
ylabel('NRMSE');
title('LASSI dictionary size');
legend(phndl,pstr{:});
axis tight; padAxis();
ylim = get(gca,'YLim');
set(gca,'YLim',[ylim(1), (ylim(2) + 0.03)]);
SetFigFontSize(fontSize);
export_fig -pdf -transparent tmi_dict_size
%--------------------------------------------------------------------------

%{
%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
cfigure([437, 358]);

% Sort data, if requested
if exist('sortIdx','var') && ~isempty(sortIdx)
    [~, idx] = sort(cellfun(@(n)n(sortIdx),{data.nrmse}),'descend');
end

% Knobs
cm = linspecer(nData);
mk = 'o*v^d<>';
markerSize = 8;
fontSize   = 12;
lineWidth  = 1;

% Plot NRMSE
phndl = zeros(1,nData);
for i = 1:nData
    ii = idx(i);
    phndl(i) = semilogx(data(ii).ps,data(ii).nrmse,'-','Marker',mk(ii),'MarkerSize',markerSize,'LineWidth',lineWidth,'Color',cm(ii,:));
    hold on;
end
xlabel('ps');
ylabel('NRMSE');
title('LASSI dictionary size');
legend(phndl,data(idx).label);
axis tight; padAxis();
SetFigFontSize(fontSize);

% Save figure, if requested
if exist('outpath','var') && ~isempty(outpath)
    export_fig('-pdf','-transparent',outpath);
end
%--------------------------------------------------------------------------
%}

%% [PLOT] LASSI empirical convergence

%{
% LASSI variations (8x accel - 50 iters)
data(1).inpath = 'data_TMI/recon_otazo_full_par6.mat'; % SVT, hard, L+S visual
data(2).inpath = 'data_TMI/recon_otazo_full_par7.mat'; % SVT, soft, L+S visual
data(3).inpath = 'data_TMI/recon_otazo_full_par1.mat'; % OPT, hard, L+S visual
data(4).inpath = 'data_TMI/recon_otazo_full_par8.mat'; % OPT, soft, L+S visual
data(1).label = 'LASSI (SVT, L0)';
data(2).label = 'LASSI (SVT, L1)';
data(3).label = 'LASSI (OPT, L0)';
data(4).label = 'LASSI (OPT, L1)';
its = 1:50;
%}

% LASSI variations (8x accel - 200 iters)
data(1).inpath = 'data_TMI/recon_otazo_full_par206.mat'; % SVT, hard, L+S visual
data(2).inpath = 'data_TMI/recon_otazo_full_par207.mat'; % SVT, soft, L+S visual
data(3).inpath = 'data_TMI/recon_otazo_full_par201.mat'; % OPT, hard, L+S visual
data(4).inpath = 'data_TMI/recon_otazo_full_par208.mat'; % OPT, soft, L+S visual
data(1).label = 'LASSI (SVT, L0)';
data(2).label = 'LASSI (SVT, L1)';
data(3).label = 'LASSI (OPT, L0)';
data(4).label = 'LASSI (OPT, L1)';
its = 1:100;

%{
% LASSI variations (20x accel)
data(1).inpath = 'data_TMI/recon_otazo_full_par14.mat'; % SVT, hard, L+S visual
data(2).inpath = 'data_TMI/recon_otazo_full_par15.mat'; % SVT, soft, L+S visual
data(3).inpath = 'data_TMI/recon_otazo_full_par9.mat';  % OPT, hard, L+S visual
data(4).inpath = 'data_TMI/recon_otazo_full_par16.mat'; % OPT, soft, L+S visual
data(1).label = 'LASSI (SVT, L0)';
data(2).label = 'LASSI (SVT, L1)';
data(3).label = 'LASSI (OPT, L0)';
data(4).label = 'LASSI (OPT, L1)';
its = 1:50;
%}

%{
% Different datasets
data(1).inpath = 'data_TMI/recon_otazo_full_par1.mat';	% Otazo 8x
data(2).inpath = 'data_TMI/recon_invivo_full_par1.mat'; % invivo 8x
data(3).inpath = 'data_TMI/recon_pincat_full_par1.mat'; % PINCAT 9x
data(1).label = 'LASSI (otazo)';
data(2).label = 'LASSI (invivo)';
data(3).label = 'LASSI (PINCAT)';
its = 1:50;
%}

%{
% Different datasets
data(1).inpath = 'data_TMI/recon_otazo_full_par9.mat';	% Otazo 20x
data(2).inpath = 'data_TMI/recon_invivo_full_par5.mat'; % invivo 5x
data(3).inpath = 'data_TMI/recon_pincat_full_par5.mat'; % PINCAT 6x
data(1).label = 'LASSI (otazo)';
data(2).label = 'LASSI (invivo)';
data(3).label = 'LASSI (PINCAT)';
its = 1:50;
%}

% Load data
nData = numel(data);
clear stats;
for i = 1:nData
    % Load data
    datai = load(data(i).inpath,'stats');
    
    % Record stats
    stats(i) = datai.stats; %#ok
end

% Sort by metrics
[~, idxobj] = sort(cellfun(@(n) n(its(end)),{stats(1:2).cost}),'descend');
[~, idxnrm] = sort(cellfun(@(n) n(its(end)),{stats.nrmse}),'descend');
[~, idxcon] = sort(cellfun(@(n) n(its(end)),{stats.delta}),'descend');
[~, idxspa] = sort(cellfun(@(n) n(its(end)),{stats.sparsity}),'descend');

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
% Knobs
cm = linspecer(nData);
mk = 'o*v^<>d*';
markerSize = 8;
nMarkers   = 10;
fontSize   = 12;
lineWidth  = 1;
spacing    = 'x';

% Cost
cfigure([348, 296]);
nObj  = numel(idxobj);
phndl = zeros(1,nObj);
for i = 1:nObj
    ii = idxobj(i);
    phndl(i) = line_fewer_markers(its,stats(ii).cost(its),nMarkers,'-','Marker',mk(ii),'MarkerSize',markerSize,'LineWidth',lineWidth,'Spacing',spacing,'Color',cm(ii,:));
    %set(gca,'YScale','log');
    hold on; box on;
end
%axis tight; padAxis();
padAxis([],[],'x');
xTick = get(gca,'XTick'); xTick(1) = 1; set(gca,'XTick',xTick);
legend(phndl,data(idxobj).label);
title('Objective');
xlabel('Iteration');
SetFigFontSize(fontSize);
export_fig -pdf -transparent tmi_empconv_objective

% NRMSE
cfigure([348, 296]);
nNrm  = numel(idxnrm);
phndl = zeros(1,nNrm);
for i = 1:nNrm
    ii = idxnrm(i);
    phndl(i) = line_fewer_markers(its,stats(ii).nrmse(its),nMarkers,'-','Marker',mk(ii),'MarkerSize',markerSize,'LineWidth',lineWidth,'Spacing',spacing,'Color',cm(ii,:));
    hold on; box on;
end
%axis tight; padAxis();
padAxis([],[],'x');
xTick = get(gca,'XTick'); xTick(1) = 1; set(gca,'XTick',xTick);
legend(phndl,data(idxnrm).label);
title('NRMSE');
xlabel('Iteration');
SetFigFontSize(fontSize);
export_fig -pdf -transparent tmi_empconv_nrmse

% Iterate convergence
cfigure([348, 296]);
nCon  = numel(idxcon);
phndl = zeros(1,nCon);
for i = 1:nCon
    ii = idxcon(i);
    phndl(i) = line_fewer_markers(its,stats(ii).delta(its),nMarkers,'-','Marker',mk(ii),'MarkerSize',markerSize,'LineWidth',lineWidth,'Spacing',spacing,'Color',cm(ii,:));
    set(gca,'YScale','log');
    hold on; box on;
end
%axis tight; padAxis();
padAxis([],[],'x');
xTick = get(gca,'XTick'); xTick(1) = 1; set(gca,'XTick',xTick);
legend(phndl,data(idxcon).label);
title('Iterate convergence');
xlabel('Iteration');
SetFigFontSize(fontSize);
export_fig -pdf -transparent tmi_empconv_itconv

% Sparsity
cfigure([348, 296]);
nSpa  = numel(idxspa);
phndl = zeros(1,nSpa);
for i = 1:nSpa
    ii = idxspa(i);
    phndl(i) = line_fewer_markers(its,stats(ii).sparsity(its),nMarkers,'-','Marker',mk(ii),'MarkerSize',markerSize,'LineWidth',lineWidth,'Spacing',spacing,'Color',cm(ii,:));
    %set(gca,'YScale','log');
    hold on; box on;
end
%axis tight; padAxis();
padAxis([],[],'x');
xTick = get(gca,'XTick'); xTick(1) = 1; set(gca,'XTick',xTick);
legend(phndl,data(idxspa).label);
title('Nonzero coefficients (%)');
xlabel('Iteration');
SetFigFontSize(fontSize);
export_fig -pdf -transparent tmi_empconv_sparsity
%--------------------------------------------------------------------------

%% [PLOT] LASSI sparsity experiments (R8)

% Knobs
data(1).inpath = './LASSI_R8_par1i.mat';
data(1).label  = '8x accel.';
outpath = 'tmi_lassi_sparsity_R8.pdf';

% Load data
nData = numel(data);
for i = 1:nData
    % Load data
    datai  = load(data(i).inpath);
    
    % Record results
    data(i).lambdaB  = datai.vars.lambdaB; %#ok
    data(i).sparsity = squeeze(datai.stats.sparsity(1,1,1,:,1,end)); %#ok
    data(i).nrmse    = squeeze(datai.stats.nrmse(1,1,1,:,1,end)); %#ok
end

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
%cfigure([326, 269]);
cfigure([332, 272]);

% Knobs
cm = linspecer(nData);
mk = 'o*v^d<>';
markerSize = 8;
fontSize   = 12;
lineWidth  = 1;

% Plot results
phndl = zeros(1,nData);
pstr  = cell(1,nData);
for i = 1:nData
    phndl(i) = plot(data(i).sparsity,data(i).nrmse,'-','Marker',mk(i),'MarkerSize',markerSize,'LineWidth',lineWidth,'Color',cm(i,:));
    pstr{i}  = data(i).label;
    hold on;
end
xlabel('Nonzero coefficients (%)');
ylabel('NRMSE');
title('Dictionary sparsity');
legend(phndl,pstr{:},'Location','NW');
axis tight; padAxis();
SetFigFontSize(fontSize);
export_fig('-pdf','-transparent',outpath);
%--------------------------------------------------------------------------

%% [PLOT] LASSI sparsity experiments (full)

% Knobs
inpath  = 'LASSI_full_sparsity_par1.mat';
outpath = 'tmi_lassi_sparsity.pdf';

% Load data
datai = load(inpath);
nData = numel(datai.vars.ps);
for i = 1:nData
    data(i).lambdaB  = datai.vars.lambdaB; %#ok
    data(i).sparsity = squeeze(datai.stats.sparsity(i,1,:,end)); %#ok
    data(i).nrmse    = squeeze(datai.stats.nrmse(i,1,:,end)); %#ok
    data(i).accel    = round(1 ./ datai.vars.ps(i)); %#ok
end

% Filter data
accelWhitelist = [8, 12, 16, 20];
idx   = ismember([data.accel],accelWhitelist);
data  = data(idx);
nData = numel(data);

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
%cfigure([326, 269]);
cfigure([332, 272]);

GAP = 0.04;

% Knobs
cm = linspecer(nData);
mk = 'o*v^d<>';
markerSize = 8;
fontSize   = 12;
lineWidth  = 1;

% Plot results
phndl = zeros(1,nData);
pstr  = cell(1,nData);
for i = 1:nData
    ii = nData + 1 - i;
    phndl(ii) = plot(data(i).sparsity,data(i).nrmse,'-','Marker',mk(i),'MarkerSize',markerSize,'LineWidth',lineWidth,'Color',cm(i,:));
    pstr{ii}  = sprintf('%dx accel.',data(i).accel);
    hold on;
end
xlabel('Nonzero coefficients (%)');
ylabel('NRMSE');
title('LASSI: Effect of sparsity');
legend(phndl,pstr{:},'Location','Best');
SetFigFontSize(fontSize);

% Fix axis
axis tight;
set(gca,'XLim',[0, 100]);
ylim = get(gca,'YLim'); set(gca,'YLim',ylim + [0, GAP]);
padAxis();

export_fig('-pdf','-transparent',outpath);
%--------------------------------------------------------------------------

%% [PLOT] per frame errors

% Knobs (otazo 8x) [YES]
truthpath = 'otazo_full.mat';
lassipath = 'data_TMI/recon_otazo_full_par1.mat';
dinopath  = 'data_TMI/recon_otazo_full_par2.mat';
rpcapath  = 'data_TMI/recon_otazo_full_par3.mat';
ktslrpath = 'data_TMI/recon_otazo_full_par4.mat';
%fftspec   = @recon_otazo_full_par1;

%{
% Knobs (otazo 20x) [NO]
truthpath = 'otazo_full.mat';
lassipath = 'data_TMI/recon_otazo_full_par9.mat';
dinopath  = 'data_TMI/recon_otazo_full_par10.mat';
rpcapath  = 'data_TMI/recon_otazo_full_par11.mat';
ktslrpath = 'data_TMI/recon_otazo_full_par12.mat';
%fftspec   = @recon_otazo_full_par9;
%}

%{
% Knobs (invivo 8x) [YES]
truthpath = 'invivo_full.mat';
%lassipath = 'data_TMI/recon_invivo_full_par1.mat'; % r = 1 atoms
lassipath = 'data_TMI/recon_invivo_full_par5.mat'; % r = 5 atoms
dinopath  = 'data_TMI/recon_invivo_full_par2.mat';
rpcapath  = 'data_TMI/recon_invivo_full_par3.mat';
ktslrpath = 'data_TMI/recon_invivo_full_par4.mat';
%fftspec   = recon_invivo_full_par5;
%}

%{
% Knobs (invivo 5x) [NO]
truthpath = 'invivo_full.mat';
%lassipath = 'data_TMI/recon_invivo_full_par6.mat'; % r = 1 atoms
lassipath = 'data_TMI/recon_invivo_full_par10.mat'; % r = 5 atoms
dinopath  = 'data_TMI/recon_invivo_full_par7.mat';
rpcapath  = 'data_TMI/recon_invivo_full_par8.mat';
ktslrpath = 'data_TMI/recon_invivo_full_par9.mat';
%fftspec   = recon_invivo_full_par10;
%}

%{
% Knobs (pincat 9x) [YES]
truthpath = 'pincat_full.mat';
lassipath = 'data_TMI/recon_pincat_full_par1.mat';
dinopath  = 'data_TMI/recon_pincat_full_par2.mat';
rpcapath  = 'data_TMI/recon_pincat_full_par3.mat';
ktslrpath = 'data_TMI/recon_pincat_full_par4.mat';
%fftspec   = recon_pincat_full_par1;
%}

%{
% Knobs (pincat 14x) [NO]
truthpath = 'pincat_full.mat';
lassipath = 'data_TMI/recon_pincat_full_par5.mat';
dinopath  = 'data_TMI/recon_pincat_full_par6.mat';
rpcapath  = 'data_TMI/recon_pincat_full_par7.mat';
ktslrpath = 'data_TMI/recon_pincat_full_par8.mat';
%fftspec   = recon_pincat_full_par5;
%}

% Load data
truthdata = load(truthpath,'Xtrue');
lassidata = load(lassipath,'Lhat','Shat');
dinodata  = load(dinopath,'Lhat','Shat');
rpcadata  = load(rpcapath,'L0','S0');
ktslrdata = load(ktslrpath,'X0');
Xref = truthdata.Xtrue;
Xlas = reshape(lassidata.Lhat + lassidata.Shat,size(Xref));
Xdin = reshape(dinodata.Lhat + dinodata.Shat,size(Xref));
Xlps = reshape(rpcadata.L0 + rpcadata.S0,size(Xref));
Xkts = reshape(ktslrdata.X0,size(Xref));

%{
% FFT recon
addpath('data_TMI');
vars = fftspec();
[~, ~, ~, ~, Xfft] = generateCardiacPerfData(vars.ps,inf,vars.seed,vars.inpath);
%}

% Compute per-frame errors
nt = size(Xref,3);
NRMSE = @(Xhat,Xtrue) norm(Xhat(:) - Xtrue(:)) / norm(Xtrue(:));
errlas  = zeros(1,nt);
errlps  = zeros(1,nt);
errkts  = zeros(1,nt);
errdin  = zeros(1,nt);
%errfft  = zeros(1,nt);
for i = 1:nt
    errlas(i) = NRMSE(Xlas(:,:,i),Xref(:,:,i));
    errlps(i) = NRMSE(Xlps(:,:,i),Xref(:,:,i));
    errkts(i) = NRMSE(Xkts(:,:,i),Xref(:,:,i));
    errdin(i) = NRMSE(Xdin(:,:,i),Xref(:,:,i));
    %errfft(i) = NRMSE(Xfft(:,:,i),Xref(:,:,i));
end

%--------------------------------------------------------------------------
% Per-frame errors
%--------------------------------------------------------------------------
cm = linspecer(5);
fontSize   = 12;
lineWidth  = 1;

% Without DINO-KAT
cfigure([332, 272]);
phndl = zeros(1,3);
phndl(2) = plot(1:nt,errlps,'.-','MarkerSize',11,'LineWidth',lineWidth,'Color',cm(3,:)); hold on;
phndl(1) = plot(1:nt,errkts,'*-','MarkerSize',6,'LineWidth',lineWidth,'Color',cm(2,:));
phndl(3) = plot(1:nt,errlas,'-','LineWidth',lineWidth,'Color',cm(1,:));
%phndl(4) = plot(1:nt,errfft,'Marker',mk(5),'MarkerSize',markerSize,'LineWidth',lineWidth,'Color',cm(5,:));
xlabel('Frame');
ylabel('NRMSE');
title('Per-frame error');
legend(phndl,'k-t SLR','L+S','LASSI','Location','NE'); % otazo
%legend(phndl,'L+S','k-t SLR','LASSI','Location','NE'); % invivo, pincat
axis tight; padAxis();
xTick = get(gca,'XTick'); xTick(1) = 1; set(gca,'XTick',xTick);
SetFigFontSize(fontSize);
isDataset = @(t,r) strncmpi(t,r,numel(r));
if isDataset(truthpath,'otazo')
    export_fig -pdf -transparent tmi_frame_errors_otazo
elseif isDataset(truthpath,'invivo')
    export_fig -pdf -transparent tmi_frame_errors_invivo
elseif isDataset(truthpath,'pincat')
    export_fig -pdf -transparent tmi_frame_errors_pincat
end  

% With DINO-KAT
cfigure([332, 272]);
phndl = zeros(1,4);
phndl(2) = plot(1:nt,errlps,'v-','MarkerSize',6,'LineWidth',lineWidth,'Color',cm(4,:)); hold on;
phndl(1) = plot(1:nt,errkts,'*-','MarkerSize',6,'LineWidth',lineWidth,'Color',cm(3,:));
phndl(3) = plot(1:nt,errdin,'.-','MarkerSize',11,'LineWidth',lineWidth,'Color',cm(2,:));
phndl(4) = plot(1:nt,errlas,'-','LineWidth',lineWidth,'Color',cm(1,:));
%phndl(5) = plot(1:nt,errfft,'Marker',mk(5),'MarkerSize',markerSize,'LineWidth',lineWidth,'Color',cm(5,:));
xlabel('Frame');
ylabel('NRMSE');
title('Per-frame error');
legend(phndl,'k-t SLR','L+S','DINO-KAT','LASSI','Location','NE'); % otazo
%legend(phndl,'L+S','k-t SLR','DINO-KAT','LASSI','Location','NE'); % invivo, pincat
axis tight; padAxis();
xTick = get(gca,'XTick'); xTick(1) = 1; set(gca,'XTick',xTick);
SetFigFontSize(fontSize);
isDataset = @(t,r) strncmpi(t,r,numel(r));
if isDataset(truthpath,'otazo')
    export_fig -pdf -transparent tmi_frame_errors_otazo_withdk
elseif isDataset(truthpath,'invivo')
    export_fig -pdf -transparent tmi_frame_errors_invivo_withdk
elseif isDataset(truthpath,'pincat')
    export_fig -pdf -transparent tmi_frame_errors_pincat_withdk
end
%--------------------------------------------------------------------------

%% [PLOT] low-rank regularization stability

% OPT
data1(1).inpath = 'LASSI_R8_par9.mat';   % OPT
data1(1).label = 'LASSI (p = 1)';
outpath1 = 'tmi_lassi_stability1.pdf';

% p
data2(1).inpath = 'LASSI_R8_par7c.mat';  % p = 0.5
data2(2).inpath = 'LASSI_R8_par5c.mat';  % p = 0
data2(3).inpath = 'LASSI_R8_par10.mat';  % p = 1
data2(1).label = 'LASSI (p = 0.5)';
data2(2).label = 'LASSI (p = 0)';
data2(3).label = 'LASSI (p = 1)';
outpath2 = 'tmi_lassi_stability2.pdf';

% Load OPT data
nData1 = numel(data1);
for i = 1:nData1
    data1i = load(data1(i).inpath);
    data1(i).r     = data1i.vars.r; %#ok
    data1(i).nrmse = squeeze(data1i.stats.nrmse(:,1,1,1,1,end)); %#ok
end

% Load p data
nData2 = numel(data2);
for i = 1:nData2
    data2i = load(data2(i).inpath);
    data2(i).lambdaL = data2i.vars.lambdaL; %#ok
    data2(i).nrmse   = squeeze(data2i.stats.nrmse(1,:,1,1,1,end)); %#ok
end

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
% Knobs
cm = linspecer(nData1 + nData2);
mk = 'o*v^d<>';
markerSize = 8;
fontSize   = 12;
lineWidth  = 1;

% Plot OPT results
cfigure([332, 272]);
phndl1 = zeros(1,nData1);
for i = 1:nData1
    ii = i;
    phndl1(i) = plot(data1(i).r,data1(i).nrmse,'-','Marker',mk(ii),'MarkerSize',markerSize,'LineWidth',lineWidth,'Color',cm(ii,:));
    hold on;
end
xlabel('r');
ylabel('NRMSE');
legend(phndl1,data1.label,'Location','NW');
SetFigFontSize(fontSize);
axis tight; padAxis(); ax1 = gca();

% Plot p = 0.5 results
cfigure([332, 272]);
phndl2 = zeros(1,nData2);
for i = 1:nData2
    ii = nData1 + i;
    phndl2(i) = plot(data2(i).lambdaL,data2(i).nrmse,'-','Marker',mk(ii),'MarkerSize',markerSize,'LineWidth',lineWidth,'Color',cm(ii,:));
    hold on;
end
xlabel('\lambda_L');
ylabel('NRMSE');
legend(phndl2,data2.label,'Location','NE');
SetFigFontSize(fontSize);
axis tight; padAxis(); ax2 = gca();

% Save results
matchAxes([ax1, ax2],0,'y');
export_fig('-pdf','-transparent',outpath1,ax1);
export_fig('-pdf','-transparent',outpath2,ax2);
%--------------------------------------------------------------------------

%% [PLOT] SOUP-DIL representation error

% Knobs
data(1).inpath  = 'SOUP_true_par1.mat';
data(1).title   = 'Training DINO-KATs';
data(1).outpath = 'tmi_soup_repError.pdf';
sparsityRange = [0, 50];

% Load data
nData = numel(data);
for i = 1:nData
    % Load data
    datai  = load(data(i).inpath);
    
    % Extract data
    data(i).lambdaB  = datai.vars.lambdaB;
    data(i).dr       = datai.vars.dr;
    data(i).sparsity = squeeze(datai.stats.sparsity(:,:,end));
    data(i).repError = squeeze(datai.stats.repError(:,:,end));    
end

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
ax = zeros(1,nData);
for i = 1:nData
    cfigure([332, 272]);
    nCurves = numel(data(i).dr);
    
    % Knobs
    cm = linspecer(nCurves);
    mk = 'o*v^d<>';
    markerSize = 8;
    fontSize   = 12;
    lineWidth  = 1;
    
    % Plot results
    phndl = zeros(1,nCurves);
    pstr  = cell(1,nCurves);
    for j = 1:nCurves
        % Extract relevant data
        sparsity = data(i).sparsity(:,j);
        repError = data(i).repError(:,j);
        idx = (sparsity >= sparsityRange(1)) & ...
              (sparsity <= sparsityRange(2));        
        phndl(j) = plot(sparsity(idx),repError(idx),'-','Marker',mk(j),'MarkerSize',markerSize,'LineWidth',lineWidth,'Color',cm(j,:));
        pstr{j}  = sprintf('r = %d',data(i).dr(j));
        hold on;
    end
    xlabel('Nonzero coefficients (%)');
    ylabel('Representation error');
    title(data(i).title);
    legend(phndl,pstr{:},'Location','NE');
    axis tight; padAxis();
    SetFigFontSize(fontSize);
    export_fig('-pdf','-transparent',data(i).outpath);
end
%--------------------------------------------------------------------------

%% [IMAGE] otazo dictionary

% Knobs
lassipath = 'data_TMI/recon_otazo_full_par1.mat';
fontSize  = 12;
sep       = [0, 0, 0];

% Load data
lassidata = load(lassipath);

% Initial dictionary
D0 = reshape(dctmtx(320),8,8,5,320);
D0 = reshape(squeeze(D0(:,:,1,:)),64,320);

% Learned dictioanry
addpath('deps_lassi');
D  = reshape(lassidata.Dhat,8,8,5,320);
D  = reshape(squeeze(D(:,:,1,:)),64,320);

% Initial atoms
T0(:,:,1) = tilePatches(D0,[8, 8],[10, 32],1,sep(1),true);
T0(:,:,2) = tilePatches(D0,[8, 8],[10, 32],1,sep(2),true);
T0(:,:,3) = tilePatches(D0,[8, 8],[10, 32],1,sep(3),true);
T0 = im2uint8(T0);

% Atom magnitudes
Tm(:,:,1) = tilePatches(abs(D),[8, 8],[8, 40],1,sep(1),true);
Tm(:,:,2) = tilePatches(abs(D),[8, 8],[8, 40],1,sep(2),true);
Tm(:,:,3) = tilePatches(abs(D),[8, 8],[8, 40],1,sep(3),true);
Tm = im2uint8(Tm);

% Atom real-parts
Tr(:,:,1) = tilePatches(real(D),[8, 8],[10, 32],1,sep(1),true);
Tr(:,:,2) = tilePatches(real(D),[8, 8],[10, 32],1,sep(2),true);
Tr(:,:,3) = tilePatches(real(D),[8, 8],[10, 32],1,sep(3),true);
Tr = im2uint8(Tr);

% Atom imaginary-parts
Ti(:,:,1) = tilePatches(imag(D),[8, 8],[10, 32],1,sep(1),true);
Ti(:,:,2) = tilePatches(imag(D),[8, 8],[10, 32],1,sep(2),true);
Ti(:,:,3) = tilePatches(imag(D),[8, 8],[10, 32],1,sep(3),true);
Ti = im2uint8(Ti);

% Display initial dictionary
cfigure();
imshow(T0);
title('Initial atoms');
SetFigFontSize(fontSize);
export_fig -pdf -transparent tmi_atom_init
imwrite(T0,'tmi_atom_init.png','png');

% Display atom magnitudes
cfigure();
imshow(Tm);
title('Atom magnitudes');
SetFigFontSize(fontSize);
export_fig -pdf -transparent tmi_atom_mag
imwrite(Tm,'tmi_atom_mag.png','png');

% Display atom real-parts
cfigure();
imshow(Tr);
title('Atom real-parts');
SetFigFontSize(fontSize);
export_fig -pdf -transparent tmi_atom_real
imwrite(Tr,'tmi_atom_real.png','png');

% Display atom imaginary-parts
figure;
imshow(Ti);
title('Atom imaginary-parts');
SetFigFontSize(fontSize);
export_fig -pdf -transparent tmi_atom_imag
imwrite(Ti,'tmi_atom_imag.png','png');

%% [MOVIE] recons (LASSI vs. L+S)

% Knobs (8x)
truthpath = 'otazo_full.mat';
lassipath = 'data_TMI/recon_otazo_full_par1.mat';
rpcapath  = 'data_TMI/recon_otazo_full_par3.mat';
alpha = 1.00;  % Alpha < 1 to clip large values
gamma = 0.90;  % Gamma < 1 to decrease contrast
mag      = 1;
Fv       = 10;
fontSize = 14;
gap      = 20;
%outpath = './tmi_otazo_8x_lassi_vs_lps.avi';  % Compressed AVI
%outpath = './tmi_otazo_8x_lassi_vs_lps.aviu'; % Uncompressed AVI
outpath = './tmi_otazo_8x_lassi_vs_lps.gif';   % GIF

%{
% Knobs (20x)
truthpath = 'otazo_full.mat';
lassipath = 'data_TMI/recon_otazo_full_par9.mat';
rpcapath  = 'data_TMI/recon_otazo_full_par11.mat';
alpha = 1.00;  % Alpha < 1 to clip large values
gamma = 0.90;  % Gamma < 1 to decrease contrast
mag      = 1;
Fv       = 10;
fontSize = 14;
gap      = 20;
%outpath = './tmi_otazo_20x_lassi_vs_lps.avi';  % Compressed AVI
%outpath = './tmi_otazo_20x_lassi_vs_lps.aviu'; % Uncompressed AVI
outpath = './tmi_otazo_20x_lassi_vs_lps.gif';   % GIF
%}

% Load data
truthdata = load(truthpath,'Xtrue');
lassidata = load(lassipath,'Lhat','Shat');
rpcapath  = load(rpcapath,'L0','S0');
Xref = truthdata.Xtrue;
Llas = reshape(lassidata.Lhat,size(Xref));
Slas = reshape(lassidata.Shat,size(Xref));
Xlas = Llas + Slas;
Llps = reshape(rpcapath.L0,size(Xref));
Slps = reshape(rpcapath.S0,size(Xref));
Xlps = Llps + Slps;

% Check NRMSEs
NRMSElas = 100 * norm(Xlas(:) - Xref(:)) / norm(Xref(:)) %#ok
NRMSElps = 100 * norm(Xlps(:) - Xref(:)) / norm(Xref(:)) %#ok

% Correct images
M     = max(abs(Xref(:)));
Xref = min(max(alpha * abs(Xref) / M,0),1).^gamma;
Xlas = min(max(alpha * abs(Xlas) / M,0),1).^gamma;
Llas = min(max(alpha * abs(Llas) / M,0),1).^gamma;
Slas = min(max(alpha * abs(Slas) / M,0),1).^gamma;
Xlps = min(max(alpha * abs(Xlps) / M,0),1).^gamma;
Llps = min(max(alpha * abs(Llps) / M,0),1).^gamma;
Slps = min(max(alpha * abs(Slps) / M,0),1).^gamma;

% Generate movie #1
C = {Xref,Xlas,Llas,Slas;
     Xref,Xlps,Llps,Slps};
M = cell2mov(C,1,1);

%{
% Generate movie #2
[ny, nx, nt] = size(Xref);
XREF = cat(1,ones(floor(0.5 * ny),nx,nt), ...
             Xref, ...
             ones(ceil(0.5 * ny) + 1,nx,nt));
C = {Xlas,Llas,Slas;
     Xlps,Llps,Slps};
M = cell2mov({XREF, cell2mov(C,1,1)},1,1);
%}

% Play movie
movie.video = M;
movie.Fv = Fv;
opts.xlabels = {'Reference','L+S','L','S'};
opts.ylabels = {sprintf('L+S (%.1f%%)',NRMSElps), ...
                sprintf('LASSI (%.1f%%)',NRMSElas)};
opts.mag      = mag;
opts.fontSize = fontSize;
opts.gap      = gap;
P = PlayMovie(movie,opts);
if strcmpi(outpath((end - 2):end),'gif')
    % Save GIF
    P.SaveGIF(outpath);
else
    % Save movie
    P.SaveMovie(outpath);
end
P.Close();

%% [MOVIE] recons (by dataset)

flipFcn    = @(X) flipdim(flipdim(X,1),2);
flipInvivo = @(C) cellfun(flipFcn,C,'UniformOutput',false);

% Knobs (otazo 8x)
truthpath = 'otazo_full.mat';
lassipath = 'data_TMI/recon_otazo_full_par1.mat';
dinopath  = 'data_TMI/recon_otazo_full_par2.mat';
rpcapath  = 'data_TMI/recon_otazo_full_par3.mat';
ktslrpath = 'data_TMI/recon_otazo_full_par4.mat';
alpha = 1.00;  % Alpha > 1 to clip large values
gamma = 0.90;  % Gamma < 1 to increase contrast
mag      = 1;
Fv       = 10;
fontSize = 14;
gap      = 20;
rows     = 1:128;
cols     = 1:128;
%outpath1 = './tmi_otazo_8x_recon.avi';         % Compressed AVI
%outpath2 = './tmi_otazo_8x_recon_nodk.avi';    % Compressed AVI
%outpath1 = './tmi_otazo_8x_recon.aviu';        % Uncompressed AVI
%outpath2 = './tmi_otazo_8x_recon_nodk.aviu';   % Uncompressed AVI
outpath1 = './tmi_otazo_8x_recon.gif';          % GIF
outpath2 = './tmi_otazo_8x_recon_nodk.gif';     % GIF

%{
% Knobs (otazo 20x)
truthpath = 'otazo_full.mat';
lassipath = 'data_TMI/recon_otazo_full_par9.mat';
dinopath  = 'data_TMI/recon_otazo_full_par10.mat';
rpcapath  = 'data_TMI/recon_otazo_full_par11.mat';
ktslrpath = 'data_TMI/recon_otazo_full_par12.mat';
alpha = 1.00;  % Alpha > 1 to clip large values
gamma = 0.90;  % Gamma < 1 to increase contrast
mag      = 1;
Fv       = 10;
fontSize = 14;
gap      = 20;
rows     = 1:128;
cols     = 1:128;
%outpath1 = './tmi_otazo_20x_recon.avi';        % Compressed AVI
%outpath2 = './tmi_otazo_20x_recon_nodk.avi';   % Compressed AVI
%outpath1 = './tmi_otazo_20x_recon.aviu';       % Uncompressed AVI
%outpath2 = './tmi_otazo_20x_recon_nodk.aviu';  % Uncompressed AVI
outpath1 = './tmi_otazo_20x_recon.gif';         % GIF
outpath2 = './tmi_otazo_20x_recon_nodk.gif';    % GIF
%}

%{
% Knobs (invivo 8x)
truthpath = 'invivo_full.mat';
%lassipath = 'data_TMI/recon_invivo_full_par1.mat'; % r = 1 atoms
lassipath = 'data_TMI/recon_invivo_full_par5.mat';  % r = 5 atoms
dinopath  = 'data_TMI/recon_invivo_full_par2.mat';
rpcapath  = 'data_TMI/recon_invivo_full_par3.mat';
ktslrpath = 'data_TMI/recon_invivo_full_par4.mat';
alpha = 1.00;  % Alpha > 1 to clip large values
gamma = 0.75;  % Gamma < 1 to increase contrast
mag      = 1.5;
Fv       = 10;
fontSize = 14;
gap      = 20;
rows     = 21:170;
cols     = 1:90;
%outpath1 = './tmi_invivo_8x_recon.avi';        % Compressed AVI
%outpath2 = './tmi_invivo_8x_recon_nodk.avi';   % Compressed AVI
%outpath1 = './tmi_invivo_8x_recon.aviu';       % Uncompressed AVI
%outpath2 = './tmi_invivo_8x_recon_nodk.aviu';  % Uncompressed AVI
outpath1 = './tmi_invivo_8x_recon.gif';         % GIF
outpath2 = './tmi_invivo_8x_recon_nodk.gif';    % GIF
%}

%{
% Knobs (invivo 5x)
truthpath = 'invivo_full.mat';
%lassipath = 'data_TMI/recon_invivo_full_par6.mat'; % r = 1 atoms
lassipath = 'data_TMI/recon_invivo_full_par10.mat'; % r = 5 atoms
dinopath  = 'data_TMI/recon_invivo_full_par7.mat';
rpcapath  = 'data_TMI/recon_invivo_full_par8.mat';
ktslrpath = 'data_TMI/recon_invivo_full_par9.mat';
alpha = 1.00;  % Alpha > 1 to clip large values
gamma = 0.75;  % Gamma < 1 to increase contrast
mag      = 1.5;
Fv       = 10;
fontSize = 14;
gap      = 20;
rows     = 21:170;
cols     = 1:90;
%outpath1 = './tmi_invivo_5x_recon.avi';        % Compressed AVI
%outpath2 = './tmi_invivo_5x_recon_nodk.avi';   % Compressed AVI
%outpath1 = './tmi_invivo_5x_recon.aviu';       % Uncompressed AVI
%outpath2 = './tmi_invivo_5x_recon_nodk.aviu';  % Uncompressed AVI
outpath1 = './tmi_invivo_5x_recon.gif';         % GIF
outpath2 = './tmi_invivo_5x_recon_nodk.gif';    % GIF
%}

%{
% Knobs (pincat 9x)
truthpath = 'pincat_full.mat';
lassipath = 'data_TMI/recon_pincat_full_par1.mat';
dinopath  = 'data_TMI/recon_pincat_full_par2.mat';
rpcapath  = 'data_TMI/recon_pincat_full_par3.mat';
ktslrpath = 'data_TMI/recon_pincat_full_par4.mat';
alpha = 1.00;  % Alpha > 1 to clip large values
gamma = 1.00;  % Gamma < 1 to increase contrast
mag      = 1;
Fv       = 10;
fontSize = 14;
gap      = 20;
rows     = 1:128;
cols     = 1:128;
%outpath1 = './tmi_pincat_9x_recon.avi';        % Compressed AVI
%outpath2 = './tmi_pincat_9x_recon_nodk.avi';   % Compressed AVI
%outpath1 = './tmi_pincat_9x_recon.aviu';       % Uncompressed AVI
%outpath2 = './tmi_pincat_9x_recon_nodk.aviu';  % Uncompressed AVI
outpath1 = './tmi_pincat_9x_recon.gif';         % GIF
outpath2 = './tmi_pincat_9x_recon_nodk.gif';    % GIF
%}

%{
% Knobs (pincat 14x)
truthpath = 'pincat_full.mat';
lassipath = 'data_TMI/recon_pincat_full_par5.mat';
dinopath  = 'data_TMI/recon_pincat_full_par6.mat';
rpcapath  = 'data_TMI/recon_pincat_full_par7.mat';
ktslrpath = 'data_TMI/recon_pincat_full_par8.mat';
alpha = 1.00;  % Alpha > 1 to clip large values
gamma = 1.00;  % Gamma < 1 to increase contrast
mag      = 1;
Fv       = 10;
fontSize = 14;
gap      = 20;
rows     = 1:128;
cols     = 1:128;
%outpath1 = './tmi_pincat_14x_recon.avi';       % Compressed AVI
%outpath2 = './tmi_pincat_14x_recon_nodk.avi';  % Compressed AVI
%outpath1 = './tmi_pincat_14x_recon.aviu';      % Uncompressed AVI
%outpath2 = './tmi_pincat_14x_recon_nodk.aviu'; % Uncompressed AVI
outpath1 = './tmi_pincat_14x_recon.gif';        % GIF
outpath2 = './tmi_pincat_14x_recon_nodk.gif';   % GIF
%}

% Load data
truthdata = load(truthpath,'Xtrue');
lassidata = load(lassipath,'Lhat','Shat');
dinodata  = load(dinopath,'Lhat','Shat');
rpcapath  = load(rpcapath,'L0','S0');
ktslrpath = load(ktslrpath,'X0');
Xref = truthdata.Xtrue;
Xlas = reshape(lassidata.Lhat + lassidata.Shat,size(Xref));
Xdin = reshape(dinodata.Lhat + dinodata.Shat,size(Xref));
Xlps = reshape(rpcapath.L0 + rpcapath.S0,size(Xref));
Xkts = reshape(ktslrpath.X0,size(Xref));

% Check NRMSEs
NRMSElas = 100 * norm(Xlas(:) - Xref(:)) / norm(Xref(:));
NRMSEdin = 100 * norm(Xdin(:) - Xref(:)) / norm(Xref(:));
NRMSElps = 100 * norm(Xlps(:) - Xref(:)) / norm(Xref(:));
NRMSEkts = 100 * norm(Xkts(:) - Xref(:)) / norm(Xref(:));

% Correct images
LX    = max(abs(Xref(:)));
Xrefd = min(max(alpha * Xref / LX,0),1).^gamma;
Xlasd = min(max(alpha * Xlas / LX,0),1).^gamma;
Xdind = min(max(alpha * Xdin / LX,0),1).^gamma;
Xlpsd = min(max(alpha * Xlps / LX,0),1).^gamma;
Xktsd = min(max(alpha * Xkts / LX,0),1).^gamma;

% Clip images
Xrefd = Xrefd(rows,cols,:);
Xlasd = Xlasd(rows,cols,:);
Xdind = Xdind(rows,cols,:);
Xlpsd = Xlpsd(rows,cols,:);
Xktsd = Xktsd(rows,cols,:);

% Generate movie #1
C1 = {Xrefd,Xlasd,Xdind,Xlpsd,Xktsd};
if ~isempty(strfind(truthpath,'invivo'))
    % Flip invivo images
    C1 = flipInvivo(C1);
end
M1 = cell2mov(C1,1,1);

% Play movie #1
lbl = @(s,n) sprintf('%s (%.1f%%)',s,n);
movie1.video = M1;
movie1.Fv    = Fv;
opts1.xlabels = {'Reference', ...
                lbl('LASSI',NRMSElas), ...
                lbl('DINO-KAT',NRMSEdin), ...
                lbl('L+S',NRMSElps), ...
                lbl('k-t SLR',NRMSEkts)};
opts1.ylabels = {''};
opts1.mag      = mag;
opts1.fontSize = fontSize;
opts1.gap      = gap;
P1 = PlayMovie(movie1,opts1);
if strcmpi(outpath1((end - 2):end),'gif')
    % Save GIF
    P1.SaveGIF(outpath1);
else
    % Save movie
    P1.SaveMovie(outpath1);
end
P1.Close();

% Generate movie #2
C2 = {Xrefd,Xlasd,Xlpsd,Xktsd};
if ~isempty(strfind(truthpath,'invivo'))
    % Flip invivo images
    C2 = flipInvivo(C2);
end
M2 = cell2mov(C2,1,1);

% Play movie #2
lbl = @(s,n) sprintf('%s (%.1f%%)',s,n);
movie2.video = M2;
movie2.Fv    = Fv;
opts2.xlabels = {'Reference', ...
                lbl('LASSI',NRMSElas), ...
                lbl('L+S',NRMSElps), ...
                lbl('k-t SLR',NRMSEkts)};
opts2.ylabels = {''};
opts2.mag      = mag;
opts2.fontSize = fontSize;
opts2.gap      = gap;
P2 = PlayMovie(movie2,opts2);
if strcmpi(outpath2((end - 2):end),'gif')
    % Save GIF
    P2.SaveGIF(outpath2);
else
    % Save movie
    P2.SaveMovie(outpath2);
end
P2.Close();

%% [MOVIE] recons + error maps (by dataset)

flipFcn    = @(X) flipdim(flipdim(X,1),2);
flipInvivo = @(C) cellfun(flipFcn,C,'UniformOutput',false);

% Knobs (otazo 8x) [YES]
truthpath = 'otazo_full.mat';
lassipath = 'data_TMI/recon_otazo_full_par1.mat';
dinopath  = 'data_TMI/recon_otazo_full_par2.mat';
rpcapath  = 'data_TMI/recon_otazo_full_par3.mat';
ktslrpath = 'data_TMI/recon_otazo_full_par4.mat';
alphaX = 1.00;  % Alpha < 1 to clip large values
gammaX = 0.90;  % Gamma < 1 to decrease contrast
alphaE = 0.40;  % Alpha < 1 to clip large values
gammaE = 1.00;  % Gamma < 1 to decrease contrast
mag      = 1;
Fv       = 10;
fontSize = 14;
gap      = 20;
rows = 1:128;
cols = 1:128;
%outpath1 = './tmi_errors_otazo_8x.avi';       % Compressed AVI
%outpath2 = './tmi_errors_otazo_8x_nodk.avi';  % Compressed AVI
%outpath1 = './tmi_errors_otazo_8x.aviu';      % Uncompressed AVI
%outpath2 = './tmi_errors_otazo_8x_nodk.aviu'; % Uncompressed AVI
outpath1 = './tmi_errors_otazo_8x.gif';        % GIF
outpath2 = './tmi_errors_otazo_8x_nodk.gif';   % GIF

%{
% Knobs (otazo 20x) [NO]
truthpath = 'otazo_full.mat';
lassipath = 'data_TMI/recon_otazo_full_par9.mat';
dinopath  = 'data_TMI/recon_otazo_full_par10.mat';
rpcapath  = 'data_TMI/recon_otazo_full_par11.mat';
ktslrpath = 'data_TMI/recon_otazo_full_par12.mat';
alphaX = 1.00;  % Alpha < 1 to clip large values
gammaX = 0.90;  % Gamma < 1 to decrease contrast
alphaE = 0.75;  % Alpha < 1 to clip large values
gammaE = 1.00;  % Gamma < 1 to decrease contrast
mag      = 1;
Fv       = 10;
fontSize = 14;
gap      = 20;
rows = 1:128;
cols = 1:128;
%outpath1 = './tmi_errors_otazo_20x.avi';       % Compressed AVI
%outpath2 = './tmi_errors_otazo_20x_nodk.avi';  % Compressed AVI
%outpath1 = './tmi_errors_otazo_20x.aviu';      % Uncompressed AVI
%outpath2 = './tmi_errors_otazo_20x_nodk.aviu'; % Uncompressed AVI
outpath1 = './tmi_errors_otazo_20x.gif';        % GIF
outpath2 = './tmi_errors_otazo_20x_nodk.gif';   % GIF
%}

%{
% Knobs (invivo 8x) [YES]
truthpath = 'invivo_full.mat';
%lassipath = 'data_TMI/recon_invivo_full_par1.mat'; % r = 1 atoms
lassipath = 'data_TMI/recon_invivo_full_par5.mat'; % r = 5 atoms
dinopath  = 'data_TMI/recon_invivo_full_par2.mat';
rpcapath  = 'data_TMI/recon_invivo_full_par3.mat';
ktslrpath = 'data_TMI/recon_invivo_full_par4.mat';
alphaX = 1.00;  % Alpha < 1 to clip large values
gammaX = 0.75;  % Gamma < 1 to decrease contrast
alphaE = 0.48;  % Alpha < 1 to clip large values
gammaE = 1.00;  % Gamma < 1 to decrease contrast
mag      = 1.5;
Fv       = 10;
fontSize = 14;
gap      = 20;
rows = 31:160;
cols = 1:90;
%outpath1 = './tmi_errors_invivo_8x.avi';       % Compressed AVI
%outpath2 = './tmi_errors_invivo_8x_nodk.avi';  % Compressed AVI
%outpath1 = './tmi_errors_invivo_8x.aviu';      % Uncompressed AVI
%outpath2 = './tmi_errors_invivo_8x_nodk.aviu'; % Uncompressed AVI
outpath1 = './tmi_errors_invivo_8x.gif';        % GIF
outpath2 = './tmi_errors_invivo_8x_nodk.gif';   % GIF
%}

%{
% Knobs (invivo 5x) [NO]
truthpath = 'invivo_full.mat';
%lassipath = 'data_TMI/recon_invivo_full_par6.mat'; % r = 1 atoms
lassipath = 'data_TMI/recon_invivo_full_par10.mat'; % r = 5 atoms
dinopath  = 'data_TMI/recon_invivo_full_par7.mat';
rpcapath  = 'data_TMI/recon_invivo_full_par8.mat';
ktslrpath = 'data_TMI/recon_invivo_full_par9.mat';
alphaX = 1.00;  % Alpha < 1 to clip large values
gammaX = 0.75;  % Gamma < 1 to decrease contrast
alphaE = 0.50;  % Alpha < 1 to clip large values
gammaE = 1.00;  % Gamma < 1 to decrease contrast
mag      = 1.5;
Fv       = 10;
fontSize = 14;
gap      = 20;
rows = 31:160;
cols = 1:90;
%outpath1 = './tmi_errors_invivo_5x.avi';       % Compressed AVI
%outpath2 = './tmi_errors_invivo_5x_nodk.avi';  % Compressed AVI
%outpath1 = './tmi_errors_invivo_5x.aviu';      % Uncompressed AVI
%outpath2 = './tmi_errors_invivo_5x_nodk.aviu'; % Uncompressed AVI
outpath1 = './tmi_errors_invivo_5x.gif';        % GIF
outpath2 = './tmi_errors_invivo_5x_nodk.gif';   % GIF
%}

%{
% Knobs (pincat 9x) [YES]
truthpath = 'pincat_full.mat';
lassipath = 'data_TMI/recon_pincat_full_par1.mat';
dinopath  = 'data_TMI/recon_pincat_full_par2.mat';
rpcapath  = 'data_TMI/recon_pincat_full_par3.mat';
ktslrpath = 'data_TMI/recon_pincat_full_par4.mat';
alphaX = 1.00;  % Alpha < 1 to clip large values
gammaX = 1.00;  % Gamma < 1 to decrease contrast
alphaE = 0.80;  % Alpha < 1 to clip large values
gammaE = 1.00;  % Gamma < 1 to decrease contrast
mag      = 1;
Fv       = 10;
fontSize = 14;
gap      = 20;
rows = 1:128;
cols = 1:128;
%outpath1 = './tmi_errors_pincat_9x.avi';       % Compressed AVI
%outpath2 = './tmi_errors_pincat_9x_nodk.avi';  % Compressed AVI
%outpath1 = './tmi_errors_pincat_9x.aviu';      % Uncompressed AVI
%outpath2 = './tmi_errors_pincat_9x_nodk.aviu'; % Uncompressed AVI
outpath1 = './tmi_errors_pincat_9x.gif';        % GIF
outpath2 = './tmi_errors_pincat_9x_nodk.gif';   % GIF
%}

%{
% Knobs (pincat 14x) [NO]
truthpath = 'pincat_full.mat';
lassipath = 'data_TMI/recon_pincat_full_par5.mat';
dinopath  = 'data_TMI/recon_pincat_full_par6.mat';
rpcapath  = 'data_TMI/recon_pincat_full_par7.mat';
ktslrpath = 'data_TMI/recon_pincat_full_par8.mat';
alphaX = 1.00;  % Alpha < 1 to clip large values
gammaX = 1.00;  % Gamma < 1 to decrease contrast
alphaE = 0.80;  % Alpha < 1 to clip large values
gammaE = 1.00;  % Gamma < 1 to decrease contrast
mag      = 1;
Fv       = 10;
fontSize = 14;
gap      = 20;
rows = 1:128;
cols = 1:128;
%outpath1 = './tmi_errors_pincat_14x.avi';       % Compressed AVI
%outpath2 = './tmi_errors_pincat_14x_nodk.avi';  % Compressed AVI
%outpath1 = './tmi_errors_pincat_14x.aviu';      % Uncompressed AVI
%outpath2 = './tmi_errors_pincat_14x_nodk.aviu'; % Uncompressed AVI
outpath1 = './tmi_errors_pincat_14x.gif';        % GIF
outpath2 = './tmi_errors_pincat_14x_nodk.gif';   % GIF
%}

% Load data
truthdata = load(truthpath,'Xtrue');
lassidata = load(lassipath,'Lhat','Shat');
dinodata  = load(dinopath,'Lhat','Shat');
rpcapath  = load(rpcapath,'L0','S0');
ktslrpath = load(ktslrpath,'X0');
Xref = truthdata.Xtrue;
Xlas = reshape(lassidata.Lhat + lassidata.Shat,size(Xref));
Xdin = reshape(dinodata.Lhat + dinodata.Shat,size(Xref));
Xlps = reshape(rpcapath.L0 + rpcapath.S0,size(Xref));
Xkts = reshape(ktslrpath.X0,size(Xref));
[ny, nx, nt] = size(Xref);

% Check NRMSEs
NRMSElas = 100 * norm(Xlas(:) - Xref(:)) / norm(Xref(:));
NRMSEdin = 100 * norm(Xdin(:) - Xref(:)) / norm(Xref(:));
NRMSElps = 100 * norm(Xlps(:) - Xref(:)) / norm(Xref(:));
NRMSEkts = 100 * norm(Xkts(:) - Xref(:)) / norm(Xref(:));

% Format images for display
M      = max(abs(Xref(:)));
Xrefd = min(max(alphaX * abs(Xref) / M,0),1).^gammaX;
Xlasd = min(max(alphaX * abs(Xlas) / M,0),1).^gammaX;

% Error maps
Elas = abs(Xref - Xlas);
Edin = abs(Xref - Xdin);
Elps = abs(Xref - Xlps);
Ekts = abs(Xref - Xkts);

% Format error maps for display
E      = max([Elas(:);
              Edin(:);
              Elps(:);
              Ekts(:);]) * alphaE;
Elasd = E * (min(max(Elas / E,0),1).^gammaE);
Edind = E * (min(max(Edin / E,0),1).^gammaE);
Elpsd = E * (min(max(Elps / E,0),1).^gammaE);
Ektsd = E * (min(max(Ekts / E,0),1).^gammaE);
clim   = [0, E];

% Clip data
Xrefd = Xrefd(rows,cols,:);
Xlasd = Xlasd(rows,cols,:);
Elasd = Elasd(rows,cols,:);
Edind = Edind(rows,cols,:);
Elpsd = Elpsd(rows,cols,:);
Ektsd = Ektsd(rows,cols,:);

% Knobs
%cm  = jet(64);
cm  = parula(64);
%cm = hot(64);
%cm = spectral(64); cm = cm(end/2:end,:);
val  = [255, 255, 255];
gap0 = 1;
gap1 = 5;

% Construct color images
Xrefi = colors2im(Xrefd,gray(256),[0, 1]);
Xlasi = colors2im(Xlasd,gray(256),[0, 1]);
Elasi = colors2im(Elasd,cm,clim);
Edini = colors2im(Edind,cm,clim);
Elpsi = colors2im(Elpsd,cm,clim);
Ektsi = colors2im(Ektsd,cm,clim);

% Generate movie #1
CX1 = {Xrefi,Xlasi};
CE1 = {Elasi,Edini,Elpsi,Ektsi};
if ~isempty(strfind(truthpath,'invivo'))
    % Flip invivo images
    CX1 = flipInvivo(CX1);
    CE1 = flipInvivo(CE1);
end
MX1 = cell2mov(CX1,gap0,val);
ME1 = cell2mov(CE1,gap0,val);
M1  = cell2mov({MX1,ME1},gap1,val);

% Play movie #1
movie1.video = M1;
movie1.Fv    = Fv;
lbl = @(s,n) sprintf('%s (%.1f%%)',s,n);
opts1.xlabels  = {'Reference','LASSI', ...
                 lbl('LASSI',NRMSElas), ...
                 lbl('DINO-KAT',NRMSEdin), ...
                 lbl('L+S',NRMSElps), ...
                 lbl('k-t SLR',NRMSEkts)};
%{
opts1.xlabels = {'Reference', ...
                 'LASSI', ...
                 'LASSI error', ...
                 'DINO-KAT error', ...
                 'L+S error', ...
                 'k-t SLR error'};
%}
opts1.ylabels  = {''};
opts1.mag      = mag;
opts1.fontSize = fontSize;
opts1.gap      = gap;
P1 = PlayMovie(movie1,opts1);
if strcmpi(outpath1((end - 2):end),'gif')
    % Save GIF
    P1.SaveGIF(outpath1);
else
    % Save movie
    P1.SaveMovie(outpath1);
end
P1.Close();

% Generate movie #2
CX2 = {Xrefi,Xlasi};
CE2 = {Elasi,Elpsi,Ektsi};
if ~isempty(strfind(truthpath,'invivo'))
    % Flip invivo images
    CX2 = flipInvivo(CX2);
    CE2 = flipInvivo(CE2);
end
MX2 = cell2mov(CX2,gap0,val);
ME2 = cell2mov(CE2,gap0,val);
M2  = cell2mov({MX2,ME2},gap1,val);

% Play movie #2
movie2.video = M2;
movie2.Fv    = Fv;
lbl = @(s,n) sprintf('%s (%.1f%%)',s,n);
opts2.xlabels  = {'Reference','LASSI', ...
                 lbl('LASSI',NRMSElas), ...
                 lbl('L+S',NRMSElps), ...
                 lbl('k-t SLR',NRMSEkts)};
%{
opts2.xlabels = {'Reference', ...
                 'LASSI', ...
                 'LASSI error', ...
                 'L+S error', ...
                 'k-t SLR error'};
%}
opts2.ylabels  = {''};
opts2.mag      = mag;
opts2.fontSize = fontSize;
opts2.gap      = gap;
P2 = PlayMovie(movie2,opts2);
if strcmpi(outpath2((end - 2):end),'gif')
    % Save GIF
    P2.SaveGIF(outpath2);
else
    % Save movie
    P2.SaveMovie(outpath2);
end
P2.Close();

%% [IMAGE] recons (LASSI vs. L+S)

% Knobs (8x)
truthpath = 'otazo_full.mat';
lassipath = 'data_TMI/recon_otazo_full_par1.mat';
rpcapath  = 'data_TMI/recon_otazo_full_par3.mat';
fontSize  = 14;
f1 = 7;
f2 = 13;
alpha = 1.00;  % Alpha < 1 to clip large values
gamma = 0.90;  % Gamma < 1 to decrease contrast

%{
% Knobs (20x)
truthpath = 'otazo_full.mat';
lassipath = 'data_TMI/recon_otazo_full_par9.mat';
rpcapath  = 'data_TMI/recon_otazo_full_par11.mat';
fontSize  = 12;
f1 = 7;
f2 = 13;
alpha = 1.00;  % Alpha < 1 to clip large values
gamma = 0.90;  % Gamma < 1 to decrease contrast
%}

% Load data
truthdata = load(truthpath,'Xtrue');
lassidata = load(lassipath,'Lhat','Shat');
rpcapath  = load(rpcapath,'L0','S0');
Xref = truthdata.Xtrue;
Llas = reshape(lassidata.Lhat,size(Xref));
Slas = reshape(lassidata.Shat,size(Xref));
Xlas = Llas + Slas;
Llps = reshape(rpcapath.L0,size(Xref));
Slps = reshape(rpcapath.S0,size(Xref));
Xlps = Llps + Slps;

% Extract images
M     = max(abs(Xref(:)));
Xref1 = min(max(alpha * abs(Xref(:,:,f1)) / M,0),1).^gamma;
Xlas1 = min(max(alpha * abs(Xlas(:,:,f1)) / M,0),1).^gamma;
Llas1 = min(max(alpha * abs(Llas(:,:,f1)) / M,0),1).^gamma;
Slas1 = min(max(alpha * abs(Slas(:,:,f1)) / M,0),1).^gamma;
Xlps1 = min(max(alpha * abs(Xlps(:,:,f1)) / M,0),1).^gamma;
Llps1 = min(max(alpha * abs(Llps(:,:,f1)) / M,0),1).^gamma;
Slps1 = min(max(alpha * abs(Slps(:,:,f1)) / M,0),1).^gamma;
Xref2 = min(max(alpha * abs(Xref(:,:,f2)) / M,0),1).^gamma;
Xlas2 = min(max(alpha * abs(Xlas(:,:,f2)) / M,0),1).^gamma;
Llas2 = min(max(alpha * abs(Llas(:,:,f2)) / M,0),1).^gamma;
Slas2 = min(max(alpha * abs(Slas(:,:,f2)) / M,0),1).^gamma;
Xlps2 = min(max(alpha * abs(Xlps(:,:,f2)) / M,0),1).^gamma;
Llps2 = min(max(alpha * abs(Llps(:,:,f2)) / M,0),1).^gamma;
Slps2 = min(max(alpha * abs(Slps(:,:,f2)) / M,0),1).^gamma;

% Display LASSI results
C1 = {Xref1,Xlas1,Llas1,Slas1;
      Xref2,Xlas2,Llas2,Slas2};
M1 = cell2mov(C1,1,1);
opts1.xlabels = {'Reference','LASSI','L','S'};
opts1.ylabels = {sprintf('F%d',f2),sprintf('F%d',f1)};
PlayMovie(M1,opts1);
SetFigFontSize(fontSize);
export_fig -pdf -transparent tmi_recon_lassi

% Display L+S results
C2 = {Xref1,Xlps1,Llps1,Slps1;
      Xref2,Xlps2,Llps2,Slps2};
M2 = cell2mov(C2,1,1);
opts2.xlabels = {'Reference','L+S','L','S'};
opts2.ylabels = {sprintf('F%d',f2),sprintf('F%d',f1)};
PlayMovie(M2,opts2);
SetFigFontSize(fontSize);
export_fig -pdf -transparent tmi_recon_lps2

%% [IMAGE] recons + error maps (by dataset)

flipFcn    = @(X) flipdim(flipdim(X,1),2);
flipInvivo = @(C) cellfun(flipFcn,C,'UniformOutput',false);

% Knobs (otazo 8x) [YES]
truthpath = 'otazo_full.mat';
lassipath = 'data_TMI/recon_otazo_full_par1.mat';
dinopath  = 'data_TMI/recon_otazo_full_par2.mat';
rpcapath  = 'data_TMI/recon_otazo_full_par3.mat';
ktslrpath = 'data_TMI/recon_otazo_full_par4.mat';
f1 = 7;
f2 = 13;
alphaX = 1.00;  % Alpha < 1 to clip large values
gammaX = 0.90;  % Gamma < 1 to decrease contrast
alphaE = 0.40;  % Alpha < 1 to clip large values
gammaE = 1.00;  % Gamma < 1 to decrease contrast
gammaC = 1.75;  % Colormap gamma "correction"
fpos     = [923, 272];
fontSize = 13;
rows = 1:128;
cols = 1:128;
colobarFcn = @() cobar([],4,2,2);
Emanual = 0.12;

%{
% Knobs (otazo 20x) [NO]
truthpath = 'otazo_full.mat';
lassipath = 'data_TMI/recon_otazo_full_par9.mat';
dinopath  = 'data_TMI/recon_otazo_full_par10.mat';
rpcapath  = 'data_TMI/recon_otazo_full_par11.mat';
ktslrpath = 'data_TMI/recon_otazo_full_par12.mat';
f1 = 7;
f2 = 13;
alphaX = 1.00;  % Alpha < 1 to clip large values
gammaX = 0.90;  % Gamma < 1 to decrease contrast
alphaE = 0.75;  % Alpha < 1 to clip large values
gammaE = 1.00;  % Gamma < 1 to decrease contrast
gammaC = 1.75;  % Colormap gamma "correction"
fpos     = [923, 272];
fontSize = 13;
rows = 1:128;
cols = 1:128;
colobarFcn = @() cobar([],4,2,2);
Emanual = 0.12;
%}

%{
% Knobs (invivo 8x) [YES]
truthpath = 'invivo_full.mat';
%lassipath = 'data_TMI/recon_invivo_full_par1.mat'; % r = 1 atoms
lassipath = 'data_TMI/recon_invivo_full_par5.mat'; % r = 5 atoms
dinopath  = 'data_TMI/recon_invivo_full_par2.mat';
rpcapath  = 'data_TMI/recon_invivo_full_par3.mat';
ktslrpath = 'data_TMI/recon_invivo_full_par4.mat';
f1 = 13;
f2 = 45;
alphaX = 1.00;  % Alpha < 1 to clip large values
gammaX = 0.75;  % Gamma < 1 to decrease contrast
alphaE = 0.48;  % Alpha < 1 to clip large values
gammaE = 1.00;  % Gamma < 1 to decrease contrast
gammaC = 1.50;  % Colormap gamma "correction"
fpos = [923, 396];
fontSize = 14;
rows = 31:160;
cols = 1:90;
colobarFcn = @() cobar([],4,2,2);
Emanual = 0.12;
%}

%{
% Knobs (invivo 5x) [NO]
truthpath = 'invivo_full.mat';
%lassipath = 'data_TMI/recon_invivo_full_par6.mat'; % r = 1 atoms
lassipath = 'data_TMI/recon_invivo_full_par10.mat'; % r = 5 atoms
dinopath  = 'data_TMI/recon_invivo_full_par7.mat';
rpcapath  = 'data_TMI/recon_invivo_full_par8.mat';
ktslrpath = 'data_TMI/recon_invivo_full_par9.mat';
f1 = 13;
f2 = 45;
alphaX = 1.00;  % Alpha < 1 to clip large values
gammaX = 0.75;  % Gamma < 1 to decrease contrast
alphaE = 0.50;  % Alpha < 1 to clip large values
gammaE = 1.00;  % Gamma < 1 to decrease contrast
gammaC = 1.50;  % Colormap gamma "correction"
fpos = [923, 396];
fontSize = 14;
rows = 31:160;
cols = 1:90;
colobarFcn = @() cobar([],4,2,2);
Emanual = 0.12;
%}

%{
% Knobs (pincat 9x) [YES]
truthpath = 'pincat_full.mat';
lassipath = 'data_TMI/recon_pincat_full_par1.mat';
dinopath  = 'data_TMI/recon_pincat_full_par2.mat';
rpcapath  = 'data_TMI/recon_pincat_full_par3.mat';
ktslrpath = 'data_TMI/recon_pincat_full_par4.mat';
f1 = 16;
f2 = 25;
alphaX = 1.00;  % Alpha < 1 to clip large values
gammaX = 1.00;  % Gamma < 1 to decrease contrast
alphaE = 0.80;  % Alpha < 1 to clip large values
gammaE = 1.00;  % Gamma < 1 to decrease contrast
gammaC = 1.20;  % Colormap gamma "correction"
fpos     = [923, 272];
fontSize = 13;
rows = 1:128;
cols = 1:128;
colobarFcn = @() cobar([],4,2,2);
Emanual = 0.30;
%}

%{
% Knobs (pincat 14x) [NO]
truthpath = 'pincat_full.mat';
lassipath = 'data_TMI/recon_pincat_full_par5.mat';
dinopath  = 'data_TMI/recon_pincat_full_par6.mat';
rpcapath  = 'data_TMI/recon_pincat_full_par7.mat';
ktslrpath = 'data_TMI/recon_pincat_full_par8.mat';
f1 = 16;
f2 = 25;
alphaX = 1.00;  % Alpha < 1 to clip large values
gammaX = 1.00;  % Gamma < 1 to decrease contrast
alphaE = 0.80;  % Alpha < 1 to clip large values
gammaE = 1.00;  % Gamma < 1 to decrease contrast
gammaC = 1.20;  % Colormap gamma "correction"
fpos     = [923, 272];
fontSize = 12;
rows = 1:128;
cols = 1:128;
colobarFcn = @() cobar([],4,2,2);
Emanual = 0.30;
%}

% Load data
truthdata = load(truthpath,'Xtrue');
lassidata = load(lassipath,'Lhat','Shat');
dinodata  = load(dinopath,'Lhat','Shat');
rpcapath  = load(rpcapath,'L0','S0');
ktslrpath = load(ktslrpath,'X0');
Xref = truthdata.Xtrue;
Xlas = reshape(lassidata.Lhat + lassidata.Shat,size(Xref));
Xdin = reshape(dinodata.Lhat + dinodata.Shat,size(Xref));
Xlps = reshape(rpcapath.L0 + rpcapath.S0,size(Xref));
Xkts = reshape(ktslrpath.X0,size(Xref));
[ny, nx, nt] = size(Xref);

% Extract images
Xref1 = abs(Xref(:,:,f1));
Xlas1 = abs(Xlas(:,:,f1));
Xdin1 = abs(Xdin(:,:,f1));
Xlps1 = abs(Xlps(:,:,f1));
Xkts1 = abs(Xkts(:,:,f1));
Xref2 = abs(Xref(:,:,f2));
Xlas2 = abs(Xlas(:,:,f2));
Xdin2 = abs(Xdin(:,:,f2));
Xlps2 = abs(Xlps(:,:,f2));
Xkts2 = abs(Xkts(:,:,f2));

% Format images for display
M      = max(abs(Xref(:)));
Xref1d = min(max(alphaX * Xref1 / M,0),1).^gammaX;
Xref2d = min(max(alphaX * Xref2 / M,0),1).^gammaX;
Xlas1d = min(max(alphaX * Xlas1 / M,0),1).^gammaX;
Xlas2d = min(max(alphaX * Xlas2 / M,0),1).^gammaX;

% Error maps
Elas1 = abs(Xref1 - Xlas1);
Edin1 = abs(Xref1 - Xdin1);
Elps1 = abs(Xref1 - Xlps1);
Ekts1 = abs(Xref1 - Xkts1);
Elas2 = abs(Xref2 - Xlas2);
Edin2 = abs(Xref2 - Xdin2);
Elps2 = abs(Xref2 - Xlps2);
Ekts2 = abs(Xref2 - Xkts2);

% Format error maps for display
E      = max([Elas1(:); Elas2(:); ...
              Edin1(:); Edin2(:); ...
              Elps1(:); Elps2(:); ...
              Ekts1(:); Ekts2(:)]) * alphaE;
if exist('Emanual','var') && ~isempty(Emanual)
    E = Emanual;
end
Elas1d = E * (min(max(Elas1 / E,0),1).^gammaE);
Edin1d = E * (min(max(Edin1 / E,0),1).^gammaE);
Elps1d = E * (min(max(Elps1 / E,0),1).^gammaE);
Ekts1d = E * (min(max(Ekts1 / E,0),1).^gammaE);
Elas2d = E * (min(max(Elas2 / E,0),1).^gammaE);
Edin2d = E * (min(max(Edin2 / E,0),1).^gammaE);
Elps2d = E * (min(max(Elps2 / E,0),1).^gammaE);
Ekts2d = E * (min(max(Ekts2 / E,0),1).^gammaE);
clim = [0, E];

%--------------------------------------------------------------------------
% Display reconstructions
%--------------------------------------------------------------------------
% Knobs
%cm  = jet(64);
cm  = parula(64);
%cm = hot(64);
%cm = spectral(64); cm = cm(end/2:end,:);
val  = [255, 255, 255];
gap0 = 1;
gap1 = 5;

% Apply gamma "correction" to colormap
n = size(cm,1);
x = linspace(0,1,n);
cm0 = cm;
for i = 1:3
    cm(:,i) = pchip(x,cm0(:,i),x.^gammaC);
end

% Generate images
Xref1i = uint8(repmat(255 * Xref1d,[1 1 3]));
Xref2i = uint8(repmat(255 * Xref2d,[1 1 3]));
Xlas1i = uint8(repmat(255 * Xlas1d,[1 1 3]));
Xlas2i = uint8(repmat(255 * Xlas2d,[1 1 3]));
Elas1i = colors2im(Elas1d,cm,clim);
Elas2i = colors2im(Elas2d,cm,clim);
Edin1i = colors2im(Edin1d,cm,clim);
Edin2i = colors2im(Edin2d,cm,clim);
Elps1i = colors2im(Elps1d,cm,clim);
Elps2i = colors2im(Elps2d,cm,clim);
Ekts1i = colors2im(Ekts1d,cm,clim);
Ekts2i = colors2im(Ekts2d,cm,clim);

% Clip images
Xref1i = Xref1i(rows,cols,:);
Xlas1i = Xlas1i(rows,cols,:);
Elas1i = Elas1i(rows,cols,:);
Edin1i = Edin1i(rows,cols,:);
Elps1i = Elps1i(rows,cols,:);
Ekts1i = Ekts1i(rows,cols,:);
Xref2i = Xref2i(rows,cols,:);
Xlas2i = Xlas2i(rows,cols,:);
Elas2i = Elas2i(rows,cols,:);
Edin2i = Edin2i(rows,cols,:);
Elps2i = Elps2i(rows,cols,:);
Ekts2i = Ekts2i(rows,cols,:);

for i = 1:3
    % Construct images #1
    xlabels1 = {'Reference','LASSI','LASSI error','DINO-KAT error','L+S error','k-t SLR error'};
    if i == 1
        CX1      = {Xref1i,Xlas1i; Xref2i,Xlas2i}; % Both
        CE1      = {Elas1i,Edin1i,Elps1i,Ekts1i; Elas2i,Edin2i,Elps2i,Ekts2i}; % Both
        ylabels1 = {sprintf('F%d',f1), sprintf('F%d',f2)}; % Both
    elseif i == 2
        CX1      = {Xref1i,Xlas1i}; % First
        CE1      = {Elas1i,Edin1i,Elps1i,Ekts1i}; % First
        ylabels1 = {sprintf('F%d',f1)}; % First
    elseif i == 3
        CX1      = {Xref2i,Xlas2i}; % Second
        CE1      = {Elas2i,Edin2i,Elps2i,Ekts2i}; % Second
        ylabels1 = {sprintf('F%d',f2)}; % Second
    end
    if ~isempty(strfind(truthpath,'invivo'))
        % Flip invivo images
        CX1 = flipInvivo(CX1);
        CE1 = flipInvivo(CE1);
    end
    MX1 = cell2mov(CX1,gap0,val);
    ME1 = cell2mov(CE1,gap0,val);
    M1  = cell2mov({MX1,ME1},gap1,val);
    
    % Construct images #2
    xlabels2 = {'Reference','LASSI','LASSI error','L+S error','k-t SLR error'};
    if i == 1
        CX2      = {Xref1i,Xlas1i; Xref2i,Xlas2i}; % Both
        CE2      = {Elas1i,Elps1i,Ekts1i; Elas2i,Elps2i,Ekts2i}; % Both
        ylabels2 = {sprintf('F%d',f1), sprintf('F%d',f2)}; % Both
    elseif i == 2
        CX2      = {Xref1i,Xlas1i}; % First
        CE2      = {Elas1i,Elps1i,Ekts1i}; % First
        ylabels2 = {sprintf('F%d',f1)}; % First
    elseif i == 3
        CX2      = {Xref2i,Xlas2i}; % Second
        CE2      = {Elas2i,Elps2i,Ekts2i}; % Second
        ylabels2 = {sprintf('F%d',f2)}; % Second
    end
    if ~isempty(strfind(truthpath,'invivo'))
        % Flip invivo images
        CX2 = flipInvivo(CX2);
        CE2 = flipInvivo(CE2);
    end
    MX2 = cell2mov(CX2,gap0,val);
    ME2 = cell2mov(CE2,gap0,val);
    M2  = cell2mov({MX2,ME2},gap1,val);
    
    % Display results #1
    cfigure(fpos);
    imagesc(M1);
    axis equal; axis tight; axis off;
    addXLabels(xlabels1,[],0.01,'top' ,'FontSize',fontSize);
    addYLabels(ylabels1,[],0.002,'left','FontSize',fontSize);
    caxis(clim);
    colormap(cm);
    colobarFcn();
    SetFigFontSize(fontSize);
    isDataset = @(t,r) strncmpi(t,r,numel(r));
    if i == 1
        if isDataset(truthpath,'otazo')
            export_fig -pdf -transparent tmi_errors_otazo_8x_F7_F13
        elseif isDataset(truthpath,'invivo')
            export_fig -pdf -transparent tmi_errors_invivo_8x_F13_F45
        elseif isDataset(truthpath,'pincat')
            export_fig -pdf -transparent tmi_errors_pincat_9x_F16_F25
        end
    elseif i == 2
        if isDataset(truthpath,'otazo')
            export_fig -pdf -transparent tmi_errors_otazo_8x_F7
        elseif isDataset(truthpath,'invivo')
            export_fig -pdf -transparent tmi_errors_invivo_8x_F13
        elseif isDataset(truthpath,'pincat')
            export_fig -pdf -transparent tmi_errors_pincat_9x_F16
        end
    elseif i == 3
        if isDataset(truthpath,'otazo')
            export_fig -pdf -transparent tmi_errors_otazo_8x_F13
        elseif isDataset(truthpath,'invivo')
            export_fig -pdf -transparent tmi_errors_invivo_8x_F45
        elseif isDataset(truthpath,'pincat')
            export_fig -pdf -transparent tmi_errors_pincat_9x_F25
        end
    end
    close();
    
    % Display results #2
    cfigure(fpos);
    imagesc(M2);
    axis equal; axis tight; axis off;
    addXLabels(xlabels2,[],0,'top' ,'FontSize',fontSize);
    addYLabels(ylabels2,[],0,'left','FontSize',fontSize);
    caxis(clim);
    colormap(cm);
    colobarFcn();
    SetFigFontSize(fontSize);
    isDataset = @(t,r) strncmpi(t,r,numel(r));
    if i == 1
        if isDataset(truthpath,'otazo')
            export_fig -pdf -transparent tmi_errors_otazo_8x_F7_F13_nodk
        elseif isDataset(truthpath,'invivo')
            export_fig -pdf -transparent tmi_errors_invivo_8x_F13_F45_nodk
        elseif isDataset(truthpath,'pincat')
            export_fig -pdf -transparent tmi_errors_pincat_9x_F16_F25_nodk
        end
    elseif i == 2
        if isDataset(truthpath,'otazo')
            export_fig -pdf -transparent tmi_errors_otazo_8x_F7_nodk
        elseif isDataset(truthpath,'invivo')
            export_fig -pdf -transparent tmi_errors_invivo_8x_F13_nodk
        elseif isDataset(truthpath,'pincat')
            export_fig -pdf -transparent tmi_errors_pincat_9x_F16_nodk
        end
    elseif i == 3
        if isDataset(truthpath,'otazo')
            export_fig -pdf -transparent tmi_errors_otazo_8x_F13_nodk
        elseif isDataset(truthpath,'invivo')
            export_fig -pdf -transparent tmi_errors_invivo_8x_F45_nodk
        elseif isDataset(truthpath,'pincat')
            export_fig -pdf -transparent tmi_errors_pincat_9x_F25_nodk
        end
    end
    close();
end
%--------------------------------------------------------------------------

%% [IMAGE] x-t slices (by dataset)

flipFcn    = @(X) flipdim(flipdim(X,1),2);
flipInvivo = @(C) cellfun(flipFcn,C,'UniformOutput',false);

%{
% Knobs (otazo 8x) [YES]
truthpath = 'otazo_full.mat';
lassipath = 'data_TMI/recon_otazo_full_par1.mat';
dinopath  = 'data_TMI/recon_otazo_full_par2.mat';
rpcapath  = 'data_TMI/recon_otazo_full_par3.mat';
ktslrpath = 'data_TMI/recon_otazo_full_par4.mat';
alpha = 1.00;  % Alpha < 1 to clip large values
gamma = 0.90;  % Gamma < 1 to decrease contrast
fontSize = 13;
f    = 7;
row  = 64;
%rows = 25:102;
%cols = 27:91;
%fpos = [770, 257];
rows = 1:128;
cols = 1:128;
fpos = [770, 206];
%}

%{
% Knobs (otazo 20x) [NO]
truthpath = 'otazo_full.mat';
lassipath = 'data_TMI/recon_otazo_full_par9.mat';
dinopath  = 'data_TMI/recon_otazo_full_par10.mat';
rpcapath  = 'data_TMI/recon_otazo_full_par11.mat';
ktslrpath = 'data_TMI/recon_otazo_full_par12.mat';
alpha = 1.00;  % Alpha < 1 to clip large values
gamma = 0.90;  % Gamma < 1 to decrease contrast
fontSize = 13;
f    = 7;
row  = 64;
rows = 1:128;
cols = 1:128;
fpos = [770, 320];
%}

%{
% Knobs (invivo 8x) [YES]
truthpath = 'invivo_full.mat';
%lassipath = 'data_TMI/recon_invivo_full_par1.mat'; % r = 1 atoms
lassipath = 'data_TMI/recon_invivo_full_par5.mat'; % r = 5 atoms
dinopath  = 'data_TMI/recon_invivo_full_par2.mat';
rpcapath  = 'data_TMI/recon_invivo_full_par3.mat';
ktslrpath = 'data_TMI/recon_invivo_full_par4.mat';
alpha = 1.00;  % Alpha < 1 to clip large values
gamma = 0.75;  % Gamma < 1 to decrease contrast
fontSize = 13;
f    = 13;
row  = 100;
rows = 60:135;
cols = 20:85;
fpos = [770, 232];
%rows = 31:160;
%cols = 1:90;
%fpos = [770, 290];
%}

%{
% Knobs (invivo 5x) [NO]
truthpath = 'invivo_full.mat';
%lassipath = 'data_TMI/recon_invivo_full_par6.mat'; % r = 1 atoms
lassipath = 'data_TMI/recon_invivo_full_par10.mat'; % r = 5 atoms
dinopath  = 'data_TMI/recon_invivo_full_par7.mat';
rpcapath  = 'data_TMI/recon_invivo_full_par8.mat';
ktslrpath = 'data_TMI/recon_invivo_full_par9.mat';
alpha = 1.00;  % Alpha < 1 to clip large values
gamma = 0.75;  % Gamma < 1 to decrease contrast
fontSize = 13;
f    = 13;
row  = 95;
rows = 31:160;
cols = 1:90;
fpos = [770, 193];
%}

% Knobs (pincat 9x) [YES]
truthpath = 'pincat_full.mat';
lassipath = 'data_TMI/recon_pincat_full_par1.mat';
dinopath  = 'data_TMI/recon_pincat_full_par2.mat';
rpcapath  = 'data_TMI/recon_pincat_full_par3.mat';
ktslrpath = 'data_TMI/recon_pincat_full_par4.mat';
alpha = 1.00;  % Alpha < 1 to clip large values
gamma = 1.00;  % Gamma < 1 to decrease contrast
fontSize = 13;
f    = 16;
row  = 66;
%rows = 15:113;
%cols = 15:113;
%fpos = [770, 207];
rows = 1:128;
cols = 1:128;
fpos = [770, 207];

%{
% Knobs (pincat 14x) [NO]
truthpath = 'pincat_full.mat';
lassipath = 'data_TMI/recon_pincat_full_par5.mat';
dinopath  = 'data_TMI/recon_pincat_full_par6.mat';
rpcapath  = 'data_TMI/recon_pincat_full_par7.mat';
ktslrpath = 'data_TMI/recon_pincat_full_par8.mat';
alpha = 1.00;  % Alpha < 1 to clip large values
gamma = 1.00;  % Gamma < 1 to decrease contrast
fontSize = 13;
f    = 16;
row  = 64;
rows = 1:128;
cols = 1:128;
fpos = [770, 320];
%}

% Load data
truthdata = load(truthpath,'Xtrue');
lassidata = load(lassipath,'Lhat','Shat');
dinodata  = load(dinopath,'Lhat','Shat');
rpcapath  = load(rpcapath,'L0','S0');
ktslrpath = load(ktslrpath,'X0');
Xref = truthdata.Xtrue;
Xlas = reshape(lassidata.Lhat + lassidata.Shat,size(Xref));
Xdin = reshape(dinodata.Lhat + dinodata.Shat,size(Xref));
Xlps = reshape(rpcapath.L0 + rpcapath.S0,size(Xref));
Xkts = reshape(ktslrpath.X0,size(Xref));
[ny, nx, nt] = size(Xref);

% Extract x-t data
Xref1 = squeeze(Xref(row,cols,:));
Xlas1 = squeeze(Xlas(row,cols,:));
Xdin1 = squeeze(Xdin(row,cols,:));
Xlps1 = squeeze(Xlps(row,cols,:));
Xkts1 = squeeze(Xkts(row,cols,:));

% Reference image
Xref2 = Xref(rows,cols,f);

% Compute NRMSEs
nrmseFcn = @(X,Xtrue) norm(X(:) - Xtrue(:)) / norm(Xtrue(:));
nrmse_las = 100 * nrmseFcn(Xlas1,Xref1) %#ok
nrmse_din = 100 * nrmseFcn(Xdin1,Xref1) %#ok
nrmse_lps = 100 * nrmseFcn(Xlps1,Xref1) %#ok
nrmse_kts = 100 * nrmseFcn(Xkts1,Xref1) %#ok

% Format images for display
M      = max(abs(Xref(:)));
Xref1d = min(max(alpha * abs(Xref1) / M,0),1).^gamma;
Xref2d = min(max(alpha * abs(Xref2) / M,0),1).^gamma;
Xlas1d = min(max(alpha * abs(Xlas1) / M,0),1).^gamma;
Xdin1d = min(max(alpha * abs(Xdin1) / M,0),1).^gamma;
Xlps1d = min(max(alpha * abs(Xlps1) / M,0),1).^gamma;
Xkts1d = min(max(alpha * abs(Xkts1) / M,0),1).^gamma;

%--------------------------------------------------------------------------
% Display reconstructions
%--------------------------------------------------------------------------
% Knobs
val = [255, 255, 255];
gap = 1;
line = [0, 224, 0];

% Generate x-t images
Xref1i = uint8(repmat(255 * Xref1d,[1 1 3]));
Xlas1i = uint8(repmat(255 * Xlas1d,[1 1 3]));
Xdin1i = uint8(repmat(255 * Xdin1d,[1 1 3]));
Xlps1i = uint8(repmat(255 * Xlps1d,[1 1 3]));
Xkts1i = uint8(repmat(255 * Xkts1d,[1 1 3]));

% Generate reference image
r = row + 1 - rows(1);
Xref2i = uint8(repmat(255 * Xref2d,[1 1 3]));
Xref2i(r,:,:) = repmat(permute(line(:),[2, 3, 1]),[1, size(Xref2i,2), 1]);
if ~isempty(strfind(truthpath,'invivo'))
    % Flip invivo images
    Xref2i = flipFcn(Xref2i);
end

% Construct image
xlabels = {'Reference','LASSI','DINO-KAT','L+S','k-t SLR'};
C = {Xref1i, Xlas1i, Xdin1i, Xlps1i, Xkts1i};
M = cell2mov(C,gap,val);
if ~isempty(strfind(truthpath,'invivo'))
    % Flip invivo images
    M = flipdim(M,1);
end

% Display results
cfigure(fpos);
subplot(1,3,1);
imshow(Xref2i); title('Reference','FontSize',fontSize); axis normal;
subplot(1,3,2:3); imshow(M);
addXLabels(xlabels,[],0.01,'top' ,'FontSize',fontSize); axis normal;

% Save figure
isDataset = @(t,r) strncmpi(t,r,numel(r));
if isDataset(truthpath,'otazo')
    export_fig -pdf -transparent tmi_xt_otazo
elseif isDataset(truthpath,'invivo')
    export_fig -pdf -transparent tmi_xt_invivo
elseif isDataset(truthpath,'pincat')
    export_fig -pdf -transparent tmi_xt_pincat
end
%--------------------------------------------------------------------------

%% [IMAGE] recons + error maps (accel sweep)

% Knobs
truthpath = 'otazo_full.mat';
lassipath = {'data_LASSI_full_par3/data13.mat', ... %  4x
             'data_LASSI_full_par3/data15.mat', ... % 12x
             'data_LASSI_full_par3/data17.mat'};    % 20x
Xlabels = {'LASSI 4x'      ,'LASSI 12x'      ,'LASSI 20x'};
Elabels = {'LASSI 4x error','LASSI 12x error','LASSI 20x error'};
outpath = 'tmi_errors_otazo_accel_041220_F7_F13.pdf';
Emanual = 0.15;
f1 = 7;
f2 = 13;
alphaX = 1.00;  % Alpha < 1 to clip large values
gammaX = 0.90;  % Gamma < 1 to decrease contrast
alphaE = 1.00;  % Alpha < 1 to clip large values
gammaE = 1.00;  % Gamma < 1 to decrease contrast
gammaC = 1.75;  % Colormap gamma "correction"
fpos     = [923, 272];
fontSize = 12;
rows = 1:128;
cols = 1:128;
colobarFcn = @() cobar([],4,2,2);

% Load data
nPaths = numel(lassipath);
truthdata = load(truthpath,'Xtrue');
lassidata = cell(1,nPaths);
for i = 1:nPaths
    lassidata{i} = load(lassipath{i},'Lhat','Shat');
end
Xref = truthdata.Xtrue;
Xlas = cell(1,nPaths);
for i = 1:nPaths
    Xlas{i} = reshape(lassidata{i}.Lhat + lassidata{i}.Shat,size(Xref));
end
[ny, nx, nt] = size(Xref);

% Extract images
Xref1 = abs(Xref(:,:,f1));
Xref2 = abs(Xref(:,:,f2));
Xlas1 = cell(1,nPaths);
Xlas2 = cell(1,nPaths);
for i = 1:nPaths
    Xlas1{i} = abs(Xlas{i}(:,:,f1));
    Xlas2{i} = abs(Xlas{i}(:,:,f2));
end

% Format images for display
M      = max(abs(Xref(:)));
Xref1d = min(max(alphaX * Xref1 / M,0),1).^gammaX;
Xref2d = min(max(alphaX * Xref2 / M,0),1).^gammaX;
Xlas1d = cell(1,nPaths);
Xlas2d = cell(1,nPaths);
for i = 1:nPaths
    Xlas1d{i} = min(max(alphaX * Xlas1{i} / M,0),1).^gammaX;
    Xlas2d{i} = min(max(alphaX * Xlas2{i} / M,0),1).^gammaX;
end

% Error maps
Elas1 = cell(1,nPaths);
Elas2 = cell(1,nPaths);
for i = 1:nPaths
    Elas1{i} = abs(Xref1 - Xlas1{i});
    Elas2{i} = abs(Xref2 - Xlas2{i});
end

% Format error maps for display
E = max([vec([Elas1{:}]); vec([Elas2{:}])]) * alphaE;
if exist('Emanual','var') && ~isempty(Emanual)
    E = Emanual;
end
Elas1d = cell(1,nPaths);
Elas2d = cell(1,nPaths);
for i = 1:nPaths
    Elas1d{i} = E * (min(max(Elas1{i} / E,0),1).^gammaE);
    Elas2d{i} = E * (min(max(Elas2{i} / E,0),1).^gammaE);
end
clim = [0, E];

%--------------------------------------------------------------------------
% Display reconstructions
%--------------------------------------------------------------------------
% Knobs
%cm  = jet(64);
cm  = parula(64);
%cm = hot(64);
%cm = spectral(64); cm = cm(end/2:end,:);
val  = [255, 255, 255];
gap0 = 1;
gap1 = 5;

% Apply gamma "correction" to colormap
n = size(cm,1);
x = linspace(0,1,n);
cm0 = cm;
for i = 1:3
    cm(:,i) = pchip(x,cm0(:,i),x.^gammaC);
end

% Generate images
Xref1i = uint8(repmat(255 * Xref1d,[1 1 3]));
Xref2i = uint8(repmat(255 * Xref2d,[1 1 3]));
Xlas1i = cell(1,nPaths);
Xlas2i = cell(1,nPaths);
Elas1i = cell(1,nPaths);
Elas2i = cell(1,nPaths);
for i = 1:nPaths
    Xlas1i{i} = uint8(repmat(255 * Xlas1d{i},[1 1 3]));
    Xlas2i{i} = uint8(repmat(255 * Xlas2d{i},[1 1 3]));
    Elas1i{i} = colors2im(Elas1d{i},cm,clim);
    Elas2i{i} = colors2im(Elas2d{i},cm,clim);
    
    % Clip images
    Xref1i    = Xref1i(rows,cols,:);
    Xref2i    = Xref2i(rows,cols,:);
    Xlas1i{i} = Xlas1i{i}(rows,cols,:);
    Xlas2i{i} = Xlas2i{i}(rows,cols,:);
    Elas1i{i} = Elas1i{i}(rows,cols,:);
    Elas2i{i} = Elas2i{i}(rows,cols,:);
end

% Construct images
xlabels = [{'Reference'}, Xlabels, Elabels];
CX = cat(1,[{Xref1i}, Xlas1i], ...
           [{Xref2i}, Xlas2i]);
CE = cat(1,Elas1i, ...
           Elas2i);
ylabels = {sprintf('F%d',f1), sprintf('F%d',f2)};
MX = cell2mov(CX,gap0,val);
ME = cell2mov(CE,gap0,val);
M  = cell2mov({MX,ME},gap1,val);

% Display results
cfigure(fpos);
imagesc(M);
axis equal; axis tight; axis off;
addXLabels(xlabels,[],0.01 ,'top' ,'FontSize',fontSize);
addYLabels(ylabels,[],0.002,'left','FontSize',fontSize);
caxis(clim);
colormap(cm);
colobarFcn();
SetFigFontSize(fontSize);
if exist('outpath','var') && ~isempty(outpath)
    export_fig('-pdf','-transparent',outpath);
end
close();
%--------------------------------------------------------------------------

%% [IMAGE] x-t slices (accel sweep)

%{
% Knobs (otazo)
truthpath = 'otazo_full.mat';
lassipath = {'data_LASSI_full_par3/data13.mat', ... %  4x
             'data_LASSI_full_par3/data15.mat', ... % 12x
             'data_LASSI_full_par3/data17.mat'};    % 20x
Xlabels = {'LASSI 4x','LASSI 12x','LASSI 20x'};
outpath = 'tmi_xt_otazo_accel_041220.pdf';
alpha = 1.00;  % Alpha < 1 to clip large values
gamma = 0.90;  % Gamma < 1 to decrease contrast
fontSize = 12;
f    = 7;
row  = 64;
%rows = 25:102;
%cols = 27:91;
%fpos = [770, 257];
rows = 1:128;
cols = 1:128;
fpos = [770, 206];
%}

% Knobs (PINCAT)
truthpath = 'pincat_full.mat';
lassipath = {'data_pincat_LASSI_full_par1nb/data24.mat', ... % 5x
             'data_pincat_LASSI_full_par1nb/data21.mat', ... % 9x
             'data_pincat_LASSI_full_par1nb/data19.mat'};    % 27x
Xlabels = {'LASSI 5x','LASSI 9x','LASSI 27x'};
outpath = 'tmi_xt_pincat_accel_050927.pdf';
alpha = 1.00;  % Alpha < 1 to clip large values
gamma = 1.00;  % Gamma < 1 to decrease contrast
fontSize = 12;
f    = 16;
row  = 66;
%rows = 15:113;
%cols = 15:113;
%fpos = [770, 207];
rows = 1:128;
cols = 1:128;
fpos = [770, 207];

% Load data
nPaths = numel(lassipath);
truthdata = load(truthpath,'Xtrue');
lassidata = cell(1,nPaths);
for i = 1:nPaths
    lassidata{i} = load(lassipath{i},'Lhat','Shat');
end
Xref = truthdata.Xtrue;
Xlas = cell(1,nPaths);
for i = 1:nPaths
    Xlas{i} = reshape(lassidata{i}.Lhat + lassidata{i}.Shat,size(Xref));
end
[ny, nx, nt] = size(Xref);

% Extract x-t data
Xref1 = squeeze(Xref(row,cols,:));
Xlas1 = cell(1,nPaths);
for i = 1:nPaths
    Xlas1{i} = squeeze(Xlas{i}(row,cols,:));
end

% Reference image
Xref2 = Xref(rows,cols,f);

% Format images for display
M      = max(abs(Xref(:)));
Xref1d = min(max(alpha * abs(Xref1) / M,0),1).^gamma;
Xref2d = min(max(alpha * abs(Xref2) / M,0),1).^gamma;
Xlas1d = cell(1,nPaths);
for i = 1:nPaths
    Xlas1d{i} = min(max(alpha * abs(Xlas1{i}) / M,0),1).^gamma;
end

%--------------------------------------------------------------------------
% Display reconstructions
%--------------------------------------------------------------------------
% Knobs
val = [255, 255, 255];
gap = 1;
line = [0, 224, 0];

% Generate x-t images
Xref1i = uint8(repmat(255 * Xref1d,[1 1 3]));
Xlas1i = cell(1,nPaths);
for i = 1:nPaths
    Xlas1i{i} = uint8(repmat(255 * Xlas1d{i},[1 1 3]));
end

% Generate reference image
r = row + 1 - rows(1);
Xref2i = uint8(repmat(255 * Xref2d,[1 1 3]));
Xref2i(r,:,:) = repmat(permute(line(:),[2, 3, 1]),[1, size(Xref2i,2), 1]);

% Construct image
xlabels = [{'Reference'}, Xlabels];
C = [{Xref1i}, Xlas1i];
M = cell2mov(C,gap,val);

% Display results
cfigure(fpos);
subplot(1,3,1);
imshow(Xref2i); title('Reference','FontSize',fontSize); axis normal;
subplot(1,3,2:3); imshow(M);
addXLabels(xlabels,[],0.01,'top' ,'FontSize',fontSize); axis normal;

% Save figure
if exist('outpath','var') && ~isempty(outpath)
    export_fig('-pdf','-transparent',outpath);
end
%--------------------------------------------------------------------------

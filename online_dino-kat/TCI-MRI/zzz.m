%% otazo data

% Knobs
ps = 1 / 8; % 1 ./ [4, 8, 12, 16, 20, 24]
SNR = inf;
seed = 1;
inpath = '/Users/Brian/Google Drive/MATLAB/Datasets/MRI/otazo_full.mat';

% Generete undersampled data
rng(seed);
[Y, Xtrue, Xfft, b1] = generateCardiacPerfData2(ps,SNR,inpath);
M = any(Y,4);

% Baseline NRMSE
NRMSE = norm(Xfft(:) - Xtrue(:)) / norm(Xtrue(:)) %#ok

%{
% Display baseline
PlayMovie(cat(2,Xtrue,Xfft));
%}

% Compare operators
frames = 1:5;
nf = numel(frames);
[ny, nx, nt] = size(Xtrue);
A  = Emat_xyt(M,b1,[ny, nx, nt]);
Af = Emat_xyt(M(:,:,frames),b1,[ny, nx, nf]);
y  = A * Xtrue;
yf = Af * Xtrue(:,:,frames);
x  = reshape(A' * y,[ny, nx, nt]);
xf = reshape(Af' * yf,[ny, nx, nf]);
erry = norm(yf(:) - vec(y(:,:,frames,:))) %#ok
errx = norm(vec(x(:,:,frames)) - xf(:)) %#ok

%% invivo data

% Knobs
nLines = 15; % [5, 10, 15, 20, 25, 30];
SNR = inf;
seed = 1;
inpath = '/Users/Brian/Google Drive/MATLAB/Datasets/MRI/invivo_full.mat';

% Generete undersampled data
rng(seed);
[Y, Xtrue, Xfft] = generateInvivoData2(nLines,SNR,inpath);
M = (abs(Y) ~= 0);

% Baseline NRMSE
NRMSE = norm(Xfft(:) - Xtrue(:)) / norm(Xtrue(:)) %#ok

%{
% Display baseline
PlayMovie(cat(2,Xtrue,Xfft));
%}

% Compare operators
frames = 1:5;
nf = numel(frames);
[ny, nx, nt] = size(Xtrue);
A  = Afft(M,[ny, nx, nt]);
Af = Afft(M(:,:,frames),[ny, nx, nf]);
y  = A * Xtrue;
yf = Af * Xtrue(:,:,frames);
x  = reshape(A' * y,[ny, nx, nt]);
xf = reshape(Af' * yf,[ny, nx, nf]);
erry = norm(yf(:) - vec(y(:,:,frames,:))) %#ok
errx = norm(vec(x(:,:,frames)) - xf(:)) %#ok

%% PINCAT data

% Knobs
nLines = 18; % [5, 10, 15, 20, 25, 30];
SNR = 45;
seed = 1;
inpath = '/Users/Brian/Google Drive/MATLAB/Datasets/MRI/pincat_full.mat';

% Generete undersampled data
rng(seed);
[Y, Xtrue, Xfft] = generateInvivoData2(nLines,SNR,inpath);
M = (abs(Y) ~= 0);

% Baseline NRMSE
NRMSE = norm(Xfft(:) - Xtrue(:)) / norm(Xtrue(:)) %#ok

%{
% Display baseline
PlayMovie(cat(2,Xtrue,Xfft));
%}

% Compare operators
frames = 1:5;
nf = numel(frames);
[ny, nx, nt] = size(Xtrue);
A  = Afft(M,[ny, nx, nt]);
Af = Afft(M(:,:,frames),[ny, nx, nf]);
y  = A * Xtrue;
yf = Af * Xtrue(:,:,frames);
x  = reshape(A' * y,[ny, nx, nt]);
xf = reshape(Af' * yf,[ny, nx, nf]);
erry = norm(yf(:) - vec(y(:,:,frames,:))) %#ok
errx = norm(vec(x(:,:,frames)) - xf(:)) %#ok

%% effective acceleration
%
% nLines = [5, 10, 15, 20, 25, 30]
%
% invivo
% accel = [23, 12, 8, 6, 5, 4]
%
% PINCAT
% accel = [27, 14, 9, 7, 6, 5]
%

% Knobs
%inpath = '/Users/Brian/Google Drive/MATLAB/Datasets/MRI/invivo_full.mat';
%inpath = '/Users/Brian/Google Drive/MATLAB/Datasets/MRI/pincat_full.mat';

% Effective acclerations
nLines = [5, 10, 15, 20, 25, 30];
ps = zeros(size(nLines));
for i = 1:numel(nLines)
    Y = generateInvivoData2(nLines(i),inf,inpath);
    M = (abs(Y) ~= 0);
    ps(i) = nnz(M) / numel(M);
end
accel = round(1 ./ ps) %#ok

%% [PLOT] *_online*Dls_par*

%inpath = 'otazo_onlineUnitaryDls_par1.mat';
%inpath = 'otazo_onlineUnitaryDls_par2a.mat';
%inpath = 'otazo_onlineUnitaryDls_par2b.mat';
%inpath = 'otazo_onlineUnitaryDls_par2c.mat';
%inpath = 'otazo_onlineUnitaryDls_par3a.mat';
%inpath = 'otazo_onlineUnitaryDls_par3a2.mat';
%inpath = 'otazo_onlineUnitaryDls_par3b.mat';
%inpath = 'otazo_onlineUnitaryDls_par3c.mat';
%inpath = 'otazo_onlineUnitaryDls_par5a.mat';
%inpath = 'otazo_onlineUnitaryDls_par5b.mat';
%inpath = 'otazo_onlineUnitaryDls_par5c.mat';
%inpath = 'otazo_onlineUnitaryDls_par5e.mat';
%inpath = 'otazo_onlineUnitaryDls_par5f.mat';
%inpath = 'otazo_onlineDls_par1.mat';
%inpath = 'otazo_onlineDls_par1b.mat';
%inpath = 'otazo_onlineDls_par2a.mat';
%inpath = 'otazo_onlineDls_par2b.mat';
%inpath = 'otazo_onlineDls_par2c.mat';
%inpath = 'otazo_onlineDls_par2d.mat';
%inpath = 'otazo_onlineDls_par2e.mat';
%inpath = 'otazo_onlineDls_par2f.mat';
%inpath = 'otazo_onlineDls_par2f2.mat';
%inpath = 'otazo_onlineDls_par2g.mat';
%inpath = 'otazo_onlineDls_par2h.mat';
%inpath = 'otazo_onlineDls_par2h2.mat';
%inpath = 'otazo_onlineDls_par2i.mat';
%inpath = 'otazo_onlineDls_par2j.mat';
%inpath = 'otazo_onlineDls_par2k.mat';
%inpath = 'otazo_onlineDls_par2k2.mat';
%inpath = 'otazo_onlineDls_par2k3.mat';
%inpath = 'otazo_onlineDls_par2l.mat';
%inpath = 'otazo_onlineDls_par2m.mat';
%inpath = 'otazo_onlineDls_par2n.mat';
%inpath = 'otazo_onlineDls_par2o.mat';

%inpath = 'invivo_onlineUnitaryDls_par1.mat';
%inpath = 'invivo_onlineUnitaryDls_par2a.mat';
%inpath = 'invivo_onlineUnitaryDls_par2b.mat';
%inpath = 'invivo_onlineUnitaryDls_par2c.mat';
%inpath = 'invivo_onlineUnitaryDls_par3a.mat';
%inpath = 'invivo_onlineUnitaryDls_par3a2.mat';
%inpath = 'invivo_onlineUnitaryDls_par3b.mat';
%inpath = 'invivo_onlineUnitaryDls_par3c.mat';
%inpath = 'invivo_onlineDls_par1.mat';
%inpath = 'invivo_onlineDls_par1b.mat';
%inpath = 'invivo_onlineDls_par2a.mat';
%inpath = 'invivo_onlineDls_par2b.mat';
%inpath = 'invivo_onlineDls_par2c.mat';
%inpath = 'invivo_onlineDls_par2d.mat';
%inpath = 'invivo_onlineDls_par2e.mat';
%inpath = 'invivo_onlineDls_par2f.mat';
%inpath = 'invivo_onlineDls_par2f2.mat';
%inpath = 'invivo_onlineDls_par2g.mat';
%inpath = 'invivo_onlineDls_par2h.mat';
%inpath = 'invivo_onlineDls_par2h2.mat';
%inpath = 'invivo_onlineDls_par2i.mat';
%inpath = 'invivo_onlineDls_par2j.mat';
%inpath = 'invivo_onlineDls_par2k.mat';
%inpath = 'invivo_onlineDls_par2k2.mat';
%inpath = 'invivo_onlineDls_par2k3.mat';
%inpath = 'invivo_onlineDls_par2l.mat';
%inpath = 'invivo_onlineDls_par2m.mat';
%inpath = 'invivo_onlineDls_par2n.mat';

%inpath = 'pincat_onlineUnitaryDls_par1.mat';
%inpath = 'pincat_onlineUnitaryDls_par2a.mat';
%inpath = 'pincat_onlineUnitaryDls_par2b.mat';
%inpath = 'pincat_onlineUnitaryDls_par2c.mat';
%inpath = 'pincat_onlineUnitaryDls_par3a.mat';
%inpath = 'pincat_onlineUnitaryDls_par3a2.mat';
%inpath = 'pincat_onlineUnitaryDls_par3b.mat';
%inpath = 'pincat_onlineUnitaryDls_par3c.mat';
%inpath = 'pincat_onlineUnitaryDls_par5a.mat';
%inpath = 'pincat_onlineUnitaryDls_par5b.mat';
%inpath = 'pincat_onlineUnitaryDls_par5c.mat';
%inpath = 'pincat_onlineDls_par1.mat';
%inpath = 'pincat_onlineDls_par1b.mat';
%inpath = 'pincat_onlineDls_par2a.mat';
%inpath = 'pincat_onlineDls_par2b.mat';
%inpath = 'pincat_onlineDls_par2c.mat';
%inpath = 'pincat_onlineDls_par2d.mat';
%inpath = 'pincat_onlineDls_par2e.mat';
%inpath = 'pincat_onlineDls_par2f.mat';
%inpath = 'pincat_onlineDls_par2f2.mat';
%inpath = 'pincat_onlineDls_par2g.mat';
%inpath = 'pincat_onlineDls_par2h.mat';
%inpath = 'pincat_onlineDls_par2h2.mat';
%inpath = 'pincat_onlineDls_par2i.mat';
%inpath = 'pincat_onlineDls_par2j.mat';
%inpath = 'pincat_onlineDls_par2k.mat';
%inpath = 'pincat_onlineDls_par2k2.mat';
%inpath = 'pincat_onlineDls_par2k3.mat';
%inpath = 'pincat_onlineDls_par2l.mat';
%inpath = 'pincat_onlineDls_par2m.mat';
%inpath = 'pincat_onlineDls_par2n.mat';

SAVE_FIGURE = true;

% Load data (D, Xhat, err, stats, vars)
load(inpath);
[~, name, ~] = fileparts(inpath);
[vars, opts] = eval(sprintf('%s();',name));
IS_OTAZO = isfield(vars,'ps');

% Optimal parameters
NRMSE = nanmean(err.NRMSE,2); % average over seed
[mnrmse, idx] = min(NRMSE(:));
[ip, is, il, im, ir, ig] = ind2sub(size(NRMSE),idx);

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
cfigure([1190, 343]);
ax = [];

% ps/nLines
subplot(2,5,1);
if IS_OTAZO
    % ps
    if (numel(vars.ps) > 1), ax(end + 1) = gca(); end
    plot(vars.ps,squeeze(NRMSE(:,is,il,im,ir,ig)),'b-o');
    xlabel('ps');
    ylabel('NRMSE');
    title(sprintf('ps = %.3f',vars.ps(ip)));
    axis tight; padAxis();
else
    % nLines
    if (numel(vars.nLines) > 1), ax(end + 1) = gca(); end
    plot(vars.nLines,squeeze(NRMSE(:,is,il,im,ir,ig)),'b-o');
    xlabel('nLines');
    ylabel('NRMSE');
    title(sprintf('nLines = %.0f',vars.nLines(ip)));
    axis tight; padAxis();
end

% lambda
subplot(2,5,2);
if (numel(vars.lambda) > 1), ax(end + 1) = gca(); end
semilogx(vars.lambda,squeeze(NRMSE(ip,is,:,im,ir,ig)),'b-o');
xlabel('lambda');
ylabel('NRMSE');
title(sprintf('lambda = %.3f',vars.lambda(il)));
axis tight; padAxis();

% mu
subplot(2,5,3);
if (numel(vars.mu) > 1), ax(end + 1) = gca(); end
semilogx(vars.mu,squeeze(NRMSE(ip,is,il,:,ir,ig)),'b-o');
xlabel('mu');
ylabel('NRMSE');
title(sprintf('mu = %.3f',vars.mu(im)));
axis tight; padAxis();

% dr
subplot(2,5,4);
if (numel(vars.dr) > 1), ax(end + 1) = gca(); end
plot(vars.dr,squeeze(NRMSE(ip,is,il,im,:,ig)),'b-o');
xlabel('dr');
ylabel('NRMSE');
title(sprintf('dr = %d',vars.dr(ir)));
axis tight; padAxis();

% gamma
subplot(2,5,5);
%if (numel(vars.gamma) > 1), ax(end + 1) = gca(); end
plot(vars.gamma,squeeze(NRMSE(ip,is,il,im,ir,:)),'b-o'); hold on;
plot(vars.gamma2,squeeze(NRMSE(ip,is,il,im,ir,:)),'r-x'); hold on;
xlabel('gamma, gamma2');
ylabel('NRMSE');
title(sprintf('gamma = %.3f, gamma2 = %.3f',vars.gamma(ig),vars.gamma2(ig)));
axis tight; padAxis();

% Optimal cost trajectory
subplot(2,5,6);
if isfield(stats,'ccost')
    % cumulative cost
    optcost = squeeze(stats.ccost(ip,is,il,im,ir,ig,:));
else
    % time-t cost
    optcost = squeeze(stats.cost(ip,is,il,im,ir,ig,:));
end
plotOnlineDlsMetric2(optcost,vars.nItersi,vars.nIters,vars.dt,vars.nReps);
xlabel('time + iteration');
ylabel('cost');
title(sprintf('cost, (T, dt, np) = (%d, %d, %d)',vars.T,vars.dt,vars.np));
axis tight; padAxis();

% Optimal NRMSE trajectory
subplot(2,5,7);
plotOnlineDlsMetric2(squeeze(stats.nrmse(ip,is,il,im,ir,ig,:)),vars.nItersi,vars.nIters,vars.dt,vars.nReps);
xlabel('time + iteration');
ylabel('NRMSE');
title(sprintf('NRMSE = %.5f',mnrmse));
axis tight; padAxis();

% Optimal sparsity trajectory
subplot(2,5,8);
plotOnlineDlsMetric2(squeeze(stats.sparsity(ip,is,il,im,ir,ig,:)),vars.nItersi,vars.nIters,vars.dt,vars.nReps);
xlabel('time + iteration');
ylabel('sparsity (%)');
title('sparsity');
axis tight; padAxis();

% Optimal deltaX trajectory
subplot(2,5,9);
plotOnlineDlsMetric2(squeeze(stats.deltaX(ip,is,il,im,ir,ig,:)),vars.nItersi,vars.nIters,vars.dt,vars.nReps);
xlabel('time + iteration');
ylabel('deltaX');
if isfield(opts,'pgap')
    title(sprintf('deltaX, pgap = [%d, %d, %d]',opts.pgap));
else
    title(sprintf('deltaX, pgap = [%d, %d, %d]',vars.pgap{end}));
end
set(gca,'YScale','log');
axis tight; padAxis();

% Optimal deltaD trajectory
subplot(2,5,10);
plotOnlineDlsMetric2(squeeze(stats.deltaD(ip,is,il,im,ir,ig,:)),vars.nItersi,vars.nIters,vars.dt,vars.nReps);
xlabel('time + iteration');
ylabel('deltaD');
title(sprintf('deltaD, pdim = [%d, %d, %d]',opts.pdim));
set(gca,'YScale','log');
axis tight; padAxis();

% Match axes
if ~isempty(ax), matchAxes(ax,[],'y'); end

% Save figure, if requested
if SAVE_FIGURE
    [~, name, ~] = fileparts(inpath);
    export_fig('-pdf','-transparent',name);
    fprintf('DONE\n');
end
%--------------------------------------------------------------------------

%% [TABLE] *_online*Dls_par*

%{
% Otazo (baseline init)
inpaths = {
%   'otazo_onlineDls_par3a.mat', 'Online DINO-KAT [DCT]';
    'otazo_onlineDls_par3a2.mat','Online DINO-KAT [DCT]';
    'otazo_onlineDls_par3b.mat', 'Online DINO-KAT [oracle]';
    'otazo_onlineDls_par3c.mat', 'Online DINO-KAT [2 passes]';
    'otazo_onlineDls_par3e.mat', 'Online DINO-KAT [1 pass]';

%   'otazo_onlineUnitaryDls_par6a.mat', 'Online DINO-KAT [DCT + unitary]';
%   'otazo_onlineUnitaryDls_par6b.mat', 'Online DINO-KAT [unitary, 1x1]';
%   'otazo_onlineUnitaryDls_par6b2.mat','Online DINO-KAT [unitary, 1x1, 2 passes]';
%   'otazo_onlineUnitaryDls_par6c.mat', 'Online DINO-KAT [DCT + unitary, 1x1]';
%   'otazo_onlineUnitaryDls_par6d.mat', 'Online DINO-KAT [unitary, 1x1, T=10]';
%   'otazo_onlineUnitaryDls_par6e.mat', 'Online DINO-KAT [DCT + unitary]';
%   'otazo_onlineUnitaryDls_par6f.mat', 'Online DINO-KAT [DCT + unitary]';
%   'otazo_onlineUnitaryDls_par6f2.mat','Online DINO-KAT [DCT + unitary]';
%   'otazo_onlineUnitaryDls_par6f3.mat','Online DINO-KAT [DCT + unitary]';
    'otazo_onlineUnitaryDls_par6f4.mat','Online DINO-KAT [unitary]'; % DCT + unitary

    '../ISMRM/batch/batch_otazo_dinokat2.mat' , 'Batch DINO-KAT';
%   '../ISMRM/batch/batch_otazo_dinokat2b.mat', 'Batch DINO-KAT';
%   '../ISMRM/batch/batch_otazo_dinokat3.mat' , 'Batch DINO-KAT';

    '../ISMRM/batch/batch_otazo_fft.mat'     , 'Baseline';
    '../ISMRM/batch/batch_otazo_ktslr.mat'   , 'k-t SLR';
    '../ISMRM/batch/batch_otazo_rpca.mat'    , 'RPCA';
};
SEEDS = 2; SORT_COL = 3; GAIN_ROW = 3;
%}

%{
% Invivo (baseline init)
inpaths = {
    'invivo_onlineDls_par3a.mat', 'Online DINO-KAT [DCT]';
    'invivo_onlineDls_par3b.mat', 'Online DINO-KAT [oracle]';
    'invivo_onlineDls_par3d.mat', 'Online DINO-KAT [2 passes]'; % mem
    'invivo_onlineDls_par3e.mat', 'Online DINO-KAT [1 pass]';
    'invivo_onlineUnitaryDls_par4b.mat', 'Online DINO-KAT [unitary, 2 passes]'; % mem
    'invivo_onlineUnitaryDls_par4c.mat', 'Online DINO-KAT [unitary, 1 pass]';
    '../ISMRM/batch/batch_invivo_dinokat2.mat', 'Batch DINO-KAT 2';           
    '../ISMRM/batch/batch_invivo_fft.mat'     , 'Baseline';
    '../ISMRM/batch/batch_invivo_ktslr.mat'   , 'k-t SLR';
    '../ISMRM/batch/batch_invivo_rpca.mat'    , 'RPCA';
};
SEEDS = 1:3; SORT_COL = 4; GAIN_ROW = 5;
%}

%{
% PINCAT (baseline init)
inpaths = {
    'pincat_onlineDls_par3a.mat', 'Online DINO-KAT [DCT]';
    'pincat_onlineDls_par3b.mat', 'Online DINO-KAT [oracle]';
    'pincat_onlineDls_par3d.mat', 'Online DINO-KAT [2 passes]'; % mem
    'pincat_onlineDls_par3e.mat', 'Online DINO-KAT [1 pass]';

%   'pincat_onlineUnitaryDls_par6a.mat', 'Online DINO-KAT [DCT + unitary]';
    'pincat_onlineUnitaryDls_par6b.mat', 'Online DINO-KAT [unitary]'; % 1x1, 1 pass
%   'pincat_onlineUnitaryDls_par6b2.mat','Online DINO-KAT [unitary, 1x1, 2 passes]';
%   'pincat_onlineUnitaryDls_par6c.mat', 'Online DINO-KAT [DCT + unitary, 1x1]';
%   'pincat_onlineUnitaryDls_par6d.mat', 'Online DINO-KAT [unitary, 1x1, T=10]';
%   'pincat_onlineUnitaryDls_par6f.mat', 'Online DINO-KAT [DCT + unitary]';

%   '../ISMRM/batch/batch_pincat_dinokat2.mat' , 'Batch DINO-KAT';
%   '../ISMRM/batch/batch_pincat_dinokat2b.mat', 'Batch DINO-KAT';
    '../ISMRM/batch/batch_pincat_dinokat3.mat' , 'Batch DINO-KAT';

    '../ISMRM/batch/batch_pincat_fft.mat'     , 'Baseline';
    '../ISMRM/batch/batch_pincat_ktslr.mat'   , 'k-t SLR';
    '../ISMRM/batch/batch_pincat_rpca.mat'    , 'RPCA';
};
SEEDS = 2; SORT_COL = 2; GAIN_ROW = 3;
%}

%{
% Otazo (baseline + batch init)
% seed = [1, 2, 3]
inpaths = {
    'otazo_onlineDls_par3a.mat', 'Online [DCT]';
    'otazo_onlineDls_par3b.mat', 'Online [oracle]';
    'otazo_onlineDls_par3c.mat', 'Online [nReps = 2]';
    'otazo_onlineDls_par3d.mat', 'Online [nReps = 2, mem]';
    'otazo_onlineDls_par3e.mat', 'Online [nReps = 1]';
    'otazo_onlineDls_par3f2.mat','Online [nReps = 1, RPCA init]';
    'otazo_onlineDls_par3g.mat', 'Online [nReps = 1, RPCA init]';
    'otazo_onlineUnitaryDls_par4a.mat', 'Online [unitary, nReps = 2]';
    'otazo_onlineUnitaryDls_par4b.mat', 'Online [unitary, nReps = 2, mem]';
    'otazo_onlineUnitaryDls_par4c.mat', 'Online [unitary, nReps = 1]';
    'otazo_onlineUnitaryDls_par4d.mat', 'Online [unitary, nReps = 1, RPCA init]';
};
%   'otazo_onlineDls_par3f.mat', 'Online [nReps = 1, RPCA init]';
% seed = 1 only
inpaths_batch = {
    '../ISMRM/batch/batch_otazo_dinokat1.mat', 'Batch [RPCA init]';
    '../ISMRM/batch/batch_otazo_dinokat2.mat', 'Batch';
    '../ISMRM/batch/batch_otazo_fft.mat', 'Baseline';
    '../ISMRM/batch/batch_otazo_ktslr.mat', 'k-t SLR';
    '../ISMRM/batch/batch_otazo_rpca.mat', 'RPCA';
};
%    '../ISMRM/batch/batch_otazo_dinokat1b.mat', 'Batch [RPCA init]';
%    '../ISMRM/batch/batch_otazo_dinokat2b.mat', 'Batch';
%    '../ISMRM/batch/batch_otazo_dinokat3.mat',  'Batch';
inpaths = cat(1,inpaths,inpaths_batch);
SEEDS = 1; SORT_COL = 3; GAIN_ROW = 4;
%}

%{
% Invivo (baseline + batch init)
% seed = [1, 2, 3]
inpaths = {'invivo_onlineDls_par3a.mat', 'Online [DCT]';
           'invivo_onlineDls_par3b.mat', 'Online [oracle]';
           'invivo_onlineDls_par3c.mat', 'Online [nReps = 2]';
           'invivo_onlineDls_par3d.mat', 'Online [nReps = 2, mem]';
           'invivo_onlineDls_par3e.mat', 'Online [nReps = 1]';
           'invivo_onlineDls_par3f.mat', 'Online [nReps = 1, k-t SLR init]';
           'invivo_onlineDls_par3g.mat', 'Online [nReps = 2, k-t SLR init]';
           'invivo_onlineUnitaryDls_par4a.mat', 'Online [unitary, nReps = 2]';
           'invivo_onlineUnitaryDls_par4b.mat', 'Online [unitary, nReps = 2, mem]';
           'invivo_onlineUnitaryDls_par4c.mat', 'Online [unitary, nReps = 1]';
           'invivo_onlineUnitaryDls_par4d.mat', 'Online [unitary, nReps = 1, k-t SLR init]';
           'invivo_onlineUnitaryDls_par4e.mat', 'Online [unitary, nReps = 2, k-t SLR init]';
};
% seed = 1 only
inpaths_batch = {
    '../ISMRM/batch/batch_invivo_dinokat1.mat', 'Batch [k-t SLR init]';
    '../ISMRM/batch/batch_invivo_dinokat2.mat', 'Batch';
    '../ISMRM/batch/batch_invivo_fft.mat', 'Baseline';
    '../ISMRM/batch/batch_invivo_ktslr.mat', 'k-t SLR';
    '../ISMRM/batch/batch_invivo_rpca.mat', 'RPCA';
};
%    '../ISMRM/batch/batch_invivo_dinokat1b.mat', 'Batch [k-t SLR init]';
%    '../ISMRM/batch/batch_invivo_dinokat2b.mat', 'Batch';
%    '../ISMRM/batch/batch_invivo_dinokat3.mat',  'Batch';
inpaths = cat(1,inpaths,inpaths_batch);
SEEDS = 1; SORT_COL = 2; GAIN_ROW = 4;
%}

%{
% PINCAT (baseline + batch init)
% seed = [1, 2, 3]
inpaths = {'pincat_onlineDls_par3a.mat', 'Online [DCT]';
           'pincat_onlineDls_par3b.mat', 'Online [oracle]';
           'pincat_onlineDls_par3c.mat', 'Online [nReps = 2]';
           'pincat_onlineDls_par3d.mat', 'Online [nReps = 2, mem]';
           'pincat_onlineDls_par3e.mat', 'Online [nReps = 1]';
           'pincat_onlineDls_par3f2.mat','Online [nReps = 1, k-t SLR init]';
           'pincat_onlineDls_par3g.mat', 'Online [nReps = 2, k-t SLR init]';
           'pincat_onlineUnitaryDls_par4a.mat', 'Online [unitary, nReps = 2]';
           'pincat_onlineUnitaryDls_par4b.mat', 'Online [unitary, nReps = 2, mem]';
           'pincat_onlineUnitaryDls_par4c.mat', 'Online [unitary, nReps = 1]';
           'pincat_onlineUnitaryDls_par4d.mat', 'Online [unitary, nReps = 1, k-t SLR init]';
};
%          'pincat_onlineDls_par3f.mat', 'Online [nReps = 1, k-t SLR init]';
% seed = 1 only
inpaths_batch = {
    '../ISMRM/batch/batch_pincat_dinokat1.mat', 'Batch [k-t SLR init]';
    '../ISMRM/batch/batch_pincat_dinokat2.mat', 'Batch';
    '../ISMRM/batch/batch_pincat_fft.mat', 'Baseline';
    '../ISMRM/batch/batch_pincat_ktslr.mat', 'k-t SLR';
    '../ISMRM/batch/batch_pincat_rpca.mat', 'RPCA';
};
%    '../ISMRM/batch/batch_pincat_dinokat1b.mat', 'Batch [k-t SLR init]';
%    '../ISMRM/batch/batch_pincat_dinokat2b.mat', 'Batch';
%    '../ISMRM/batch/batch_pincat_dinokat3.mat',  'Batch';
inpaths = cat(1,inpaths,inpaths_batch);
SEEDS = 1; SORT_COL = 3; GAIN_ROW = 4;
%}

% Load data
nData = size(inpaths,1);
data = struct();
for k = 1:nData
    datak = load(inpaths{k,1});
    nrmse = datak.stats.nrmse;
    nd = ndims(nrmse);
    if nd == 2
        nd = 3; % handle fft
    end
    args = repmat({':'},1,nd - 1);
    nrmse = squeeze(nrmse(args{:},size(nrmse,nd))); % Last iteration only
    nrmse = nanmean(nrmse(:,SEEDS),2); % Average over trials
    data(k).label = inpaths{k,2};
    data(k).nrmse = 100 * nrmse;
end
IS_OTAZO = isfield(datak.vars,'ps');
if IS_OTAZO
    ps = datak.vars.ps;
else
    nLines = datak.vars.nLines;
end

%--------------------------------------------------------------------------
% Generate table
%--------------------------------------------------------------------------
% Headings
labels = struct();
labels.align = 'center';
if IS_OTAZO
    labelFcn = @(ps) sprintf('%dx',round(1 / ps));
    labels.strs = ['Acceleration', arrayfun(labelFcn,ps,'UniformOutput',false)];
else
    if ~isempty(strfind(lower(inpaths{1,1}),'invivo'))
        % Invivo
        NUM_LINES = [5, 10, 15, 20, 25, 30];
        ACCEL = [23, 12, 8, 6, 5, 4];
    else
        % PINCAT
        NUM_LINES = [5, 10, 15, 20, 25, 30];
        ACCEL = [27, 14, 9, 7, 6, 5];
    end
    labelFcn = @(nLines) sprintf('%dx',ACCEL(NUM_LINES == nLines));
    labels.strs = ['Acceleration', arrayfun(labelFcn,nLines,'UniformOutput',false)];
end

% Row data
rows = struct();
for k = 1:nData
    rows(k).name = data(k).label;
    rows(k).data = data(k).nrmse;
end
[rows.fmt]    = deal('%.2f%%');
[rows.nAlign] = deal('right');
[rows.dAlign] = deal('right');

% Flip invivo-PINCAT data
if ~IS_OTAZO
    labels.strs = [labels.strs{1}, flipud(reshape(labels.strs(2:end),[],1))'];
    for k = 1:nData
        rows(k).data = flipud(rows(k).data(:))';
    end
end

% Print table
printTable(labels,rows,SORT_COL,GAIN_ROW);
%--------------------------------------------------------------------------

%% [PLOT] dictionaries

BASE_PATH = '/Users/Brian/Desktop/Online DINO-KAT/TCI-mri';

%{
% otazo
% accel = [4, 8, 12, 16, 20, 24]
inpath = 'otazo_onlineDls_par3e'; % 1 pass
outpath = 'tci_otazo_dict_learned';
idx = 8;
%}

%{
% pincat
% nlines = [5, 10, 15, 20, 25, 30]
inpath = 'pincat_onlineDls_par3e'; % 1 pass
outpath = 'tci_pincat_dict_learned';
idx = 9;
%}

gap = 1; val = 1;
kk = 1;

% Load learned dictionary
data = load(sprintf('%s/%s/data%d.mat',BASE_PATH,inpath,idx));
Dhat = data.D;

% Construct initial dictionary
D0 = dctmtx(320)';

% Display dictionaries
D0 = reshape(D0,[8, 8, 5, 320]);
D0 = squeeze(D0(:,:,kk,:));
Dhat = reshape(Dhat,[8, 8, 5, 320]);
Dhat = squeeze(Dhat(:,:,kk,:));
C0 = cell(16,20);
Cr = cell(16,20);
Ci = cell(16,20);
for k = 1:320
    C0{k} = formatDict(D0(:,:,k));
    Cr{k} = formatDict(real(Dhat(:,:,k)));
    Ci{k} = formatDict(imag(Dhat(:,:,k)));
end
X0 = cell2mov(C0,gap,val);
Xr = cell2mov(Cr,gap,val);
Xi = cell2mov(Ci,gap,val);

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
cfigure([963, 320]);
axgap = [0.02, 0.02];

ha = tight_subplot(1,3,axgap,axgap,axgap); % learned

%subplot(1,3,1);
axes(ha(1)); 
imshow(X0,[0, 1]);
title('Initial','FontWeight','normal');

%subplot(1,3,2);
axes(ha(2)); 
imshow(Xr,[0, 1]);
title('Learned (real)','FontWeight','normal');

%subplot(1,3,3);
axes(ha(3)); 
imshow(Xi,[0, 1]);
title('Learned (imaginary)','FontWeight','normal');

SetFigFontSize(24);
savePDF(outpath);
%--------------------------------------------------------------------------

%% [PLOT] oracle dictionaries

% otazo
inpath = 'dicts/otazo_Dsoup_par1_data4.mat'; % otazo
%inpath = 'dicts/pincat_Dsoup_par1_data3.mat'; % pincat

gap = 1; val = 1;
kk = 1;

% oracle
data = load(inpath);
Dhat = data.Dhat;

% Display dictionaries
Dhat = reshape(Dhat,[8, 8, 5, 320]);
Dhat = squeeze(Dhat(:,:,kk,:));
ISREAL = isreal(Dhat);
Cr = cell(16,20);
Ci = cell(16,20);
for k = 1:320
    Cr{k} = formatDict(real(Dhat(:,:,k)));
    Ci{k} = formatDict(imag(Dhat(:,:,k)));
end
Xr = cell2mov(Cr,gap,val);
Xi = cell2mov(Ci,gap,val);

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
cfigure([963, 320]);
axgap = [0.02, 0.02];

ha = tight_subplot(1,3,axgap,axgap,axgap); % learned

if ISREAL
    axes(ha(1)); 
    imshow(Xr,[0, 1]);
    title('Oracle','FontWeight','normal');
    
    set(ha(2),'Visible','off');
else
    axes(ha(1)); 
    imshow(Xr,[0, 1]);
    title('Oracle (real)','FontWeight','normal');
    
    axes(ha(2)); 
    imshow(Xi,[0, 1]);
    title('Oracle (imaginary)','FontWeight','normal');
end

set(ha(3),'Visible','off');

SetFigFontSize(14);

%{
savePDF('tci_otazo_dict_oracle');
savePDF('tci_pincat_dict_oracle');
%}
%--------------------------------------------------------------------------

%% 0/6 [PLOT, IMAGE, VIDEOS] recon
%  Otazo: accel = [4x ,  8x, 12x, 16x, 20x, 24x]
% Invivo: accel = [23x, 12x,  8x,  6x,  5x,  4x]
% PINCAT: accel = [27x, 14x,  9x,  7x,  6x,  5x]

% Universal constants
FRAME_RATE = 12;
FONT_SIZE = 14;
CMX = gray(256);
CME = parula(64);
GAP = 1;
VAL = 255;

%{
% Otazo (baseline init), [TCI PAPER]
BASE = '/Users/Brian/Desktop/Online DINO-KAT/TCI-mri';
gtpath = '/Users/Brian/Google Drive/MATLAB/Datasets/MRI/otazo_full.mat';
inpaths = {
%    'otazo_onlineDls_par3b', 'Online (oracle)'; % oracle
%    'otazo_onlineDls_par3c', 'Online DINO-KAT'; % 2 passes
    'otazo_onlineDls_par3c', 'OnAIR-LD'; % 2 passes [TCI FINAL]
%    'otazo_onlineDls_par3e', 'Online DINO-KAT'; % 1 pass
%    'otazo_onlineUnitaryDls_par6f4','Online (unitary)'; % DCT + unitary
    'otazo_onlineDls_par3a2','Online (DCT)'; % DCT
%    'batch_otazo_dinokat2', 'Batch DINO-KAT';
%    'batch_otazo_fft'     , 'Baseline';
    'batch_otazo_ktslr'   , 'k-t SLR';
    'batch_otazo_rpca'    , 'L+S'; % RPCA
};
IDX = 3; % ps = 1 ./ [4, 8, 12, 16, 20, 24], seed = [1, 2, 3]
FRAMES = [4, 11];
ROI = {58,1:128}; % {X,Y}
GAMMAX = 0.8;
GAMMAE = 1.5;
ALPHAE = 0.9;
%}

%{
% Otazo (baseline init) [ASILOMAR PAPER]
BASE = '/Users/Brian/Desktop/Online DINO-KAT/TCI-mri';
gtpath = '/Users/Brian/Google Drive/MATLAB/Datasets/MRI/otazo_full.mat';
inpaths = {
    'otazo_onlineDls_par3e' , 'Online DINO-KAT'; % 1 pass
    'otazo_onlineDls_par3a2', 'Online (DCT)';
    'batch_otazo_rpca'      , 'L+S';
};
IDX = 3; % ps = 1 ./ [4, 8, 12, 16, 20, 24], seed = [1, 2, 3]
FRAMES = [4, 12];
ROI = {58,1:128}; % {X,Y}
GAMMAX = 0.8;
GAMMAE = 1.5;
ALPHAE = 0.9;
%}

%{
% Otazo (baseline init) [ASILOMAR POSTER]
BASE = '/Users/Brian/Desktop/Online DINO-KAT/TCI-mri';
gtpath = '/Users/Brian/Google Drive/MATLAB/Datasets/MRI/otazo_full.mat';
inpaths = {
    'otazo_onlineDls_par3e', 'Online DINO-KAT'; % 1 pass
    'batch_otazo_ktslr'    , 'k-t SLR';
    'batch_otazo_rpca'     , 'RPCA';
};
IDX = 3; % ps = 1 ./ [4, 8, 12, 16, 20, 24], seed = [1, 2, 3]
GAMMAX = 0.8;
GAMMAE = 1.6;
ALPHAE = 0.9;
FRAMES = [1, 6];
%}

%{
% Invivo (baseline init) [UNUSED]
BASE = '/Users/Brian/Desktop/Online DINO-KAT/TCI-mri';
gtpath = '/Users/Brian/Google Drive/MATLAB/Datasets/MRI/invivo_full.mat';
inpaths = {
%    'invivo_onlineDls_par3b', 'Online DINO-KAT [oracle]';
    'invivo_onlineDls_par3d', 'Online DINO-KAT [2 passes]'; % mem
%    'invivo_onlineDls_par3e', 'Online DINO-KAT [1 pass]';
    'invivo_onlineUnitaryDls_par4b', 'Online DINO-KAT [unitary, 2 passes]'; % mem
%    'invivo_onlineUnitaryDls_par4c', 'Online DINO-KAT [unitary, 1 pass]';
%    'invivo_onlineDls_par3a', 'Online DINO-KAT [DCT]';
%    'batch_invivo_dinokat2', 'Batch DINO-KAT 2';           
%    'batch_invivo_fft'     , 'Baseline';
    'batch_invivo_ktslr'   , 'k-t SLR';
    'batch_invivo_rpca'    , 'RPCA';
};
IDX = 3; % nLines = [5, 10, 15, 20, 25, 30], seed = [1, 2, 3]
GAMMAX = 1.0;
GAMMAE = 1.0;
ALPHAE = 1.0;
FRAMES = [XX];
%}

%{
% PINCAT (baseline init), [TCI PAPER]
BASE = '/Users/Brian/Desktop/Online DINO-KAT/TCI-mri';
gtpath = '/Users/Brian/Google Drive/MATLAB/Datasets/MRI/pincat_full.mat';
inpaths = {
%    'pincat_onlineDls_par3b', 'Online (oracle)'; % oracle
%    'pincat_onlineDls_par3d', 'Online DINO-KAT'; % 2 passes, mem
    'pincat_onlineDls_par3d', 'OnAIR-LD'; % 2 passes, mem [TCI FINAL]
%    'pincat_onlineDls_par3e', 'Online DINO-KAT'; % 1 pass
%    'pincat_onlineUnitaryDls_par6b', 'Online (unitary)'; % 1x1, 1 pass, unitary
    'pincat_onlineDls_par3a', 'Online (DCT)'; % DCT
%    'batch_pincat_dinokat3' , 'Batch DINO-KAT';
%    'batch_pincat_fft'     , 'Baseline';
    'batch_pincat_ktslr'   , 'k-t SLR';
    'batch_pincat_rpca'    , 'L+S'; % RPCA
};
IDX = 10; % nLines = [5, 10, 15, 20, 25, 30], seed = [1, 2, 3]
FRAMES = [17, 45];
ROI = {65,1:128}; % {X,Y}
GAMMAX = 1.0;
GAMMAE = 1.2;
ALPHAE = 0.8;
%}

%{
% PINCAT (baseline init) [ASILOMAR PAPER]
BASE = '/Users/Brian/Desktop/Online DINO-KAT/TCI-mri';
gtpath = '/Users/Brian/Google Drive/MATLAB/Datasets/MRI/pincat_full.mat';
inpaths = {
    'pincat_onlineDls_par3e', 'Online DINO-KAT'; % 1 pass
    'pincat_onlineDls_par3a', 'Online (DCT)';
    'batch_pincat_rpca'     , 'L+S';
};
IDX = 10; % nLines = [5, 10, 15, 20, 25, 30], seed = [1, 2, 3]
FRAMES = [17, 44];
ROI = {65,1:128}; % {X,Y}
GAMMAX = 1.0;
GAMMAE = 1.2;
ALPHAE = 0.8;
%}

%{
% PINCAT (baseline init) [ASILOMAR POSTER]
BASE = '/Users/Brian/Desktop/Online DINO-KAT/TCI-mri';
gtpath = '/Users/Brian/Google Drive/MATLAB/Datasets/MRI/pincat_full.mat';
inpaths = {
    'pincat_onlineDls_par3e', 'Online DINO-KAT'; % 1 pass
    'batch_pincat_ktslr'    , 'k-t SLR';
    'batch_pincat_rpca'     , 'RPCA';
};
IDX = 9; % nLines = [5, 10, 15, 20, 25, 30], seed = [1, 2, 3]
GAMMAX = 1.0;
GAMMAE = 1.2;
ALPHAE = 0.8;
FRAMES = [15, 49];
%}

% Load ground truth
datat = load(gtpath);
Xtrue = datat.Xtrue;
trueMax = max(abs(Xtrue(:)));

% Load recons
nData = size(inpaths,1);
Xhat = cell(1,nData);
for k = 1:nData
    datak = load(sprintf('%s/%s/data%d.mat',BASE,inpaths{k,1},IDX));
    Xhat{k} = datak.Xhat;
end

% All data
X = [{Xtrue}, Xhat{:}];
xlabels = [{'Reference'}, inpaths{:,2}];

% Compute per-frame errors
clear err;
for k = 1:nData
    err(k) = computeErrorMetrics(Xhat{k},Xtrue); %#ok
end

%% 1/6 Recon movie

% Format data
CMXg = gammaCorrectColormap(CMX,GAMMAX);
VX = abs(cell2mat(X));
VX = min(VX / trueMax,1);
VX = colors2im(VX,CMXg,[0, 1]);

% Play movie
movie = struct();
movie.video = VX;
movie.Fv = FRAME_RATE;
opts = struct();
opts.xlabels = xlabels;
opts.fontSize = FONT_SIZE;
PlayMovie(movie,opts);

%{
./asilomar_otazo_recons.gif
./asilomar_pincat_recons.gif

./tci_otazo_recons.gif
./tci_pincat_recons.gif
%}

%% 2/6 Recon + error movie

% Compute error videos
E = cell(1,nData + 1);
E{1} = zeros(size(Xtrue));
for k = 1:nData
    E{k + 1} = abs(Xhat{k} - Xtrue);
end

% Format X videos
CMXg = gammaCorrectColormap(CMX,GAMMAX);
VX = abs(cell2mat(X));
VX = min(VX / trueMax,1);
VX = colors2im(VX,CMXg,[0, max(VX(:))]);

% Format E videos
CMEg = gammaCorrectColormap(CME,GAMMAE);
VE = cell2mat(E);
VE = colors2im(VE,CMEg,[0, ALPHAE * max(VE(:))]);

% Play movie
movie = struct();
movie.video = cat(1,VX,VE);
movie.Fv = FRAME_RATE;
opts = struct();
opts.xlabels  = xlabels;
opts.ylabels  = {'Errors','Recons'};
opts.fontSize = FONT_SIZE;
PlayMovie(movie,opts);

%{
./asilomar_otazo_errors.gif
./asilomar_pincat_errors.gif

./tci_otazo_errors.gif
./tci_pincat_errors.gif
%}

%% 2/6 (thesis)

% Compute error videos
E = cell(1,nData + 1);
E{1} = zeros(size(Xtrue));
for k = 1:nData
    E{k + 1} = abs(Xhat{k} - Xtrue);
end

% Format X videos
CMXg = gammaCorrectColormap(CMX,GAMMAX);
VX = cell(size(X));
limx = -inf;
for k = 1:numel(VX)
    VX{k} = min(abs(X{k}) / trueMax,1);
    limx = max(limx,max(VX{k}(:)));
end
for k = 1:numel(VX)
    VX{k} = colors2im(VX{k},CMXg,[0, limx]);
end
VX = cell2mov(VX(1:2),1,255);

% Format E videos
CMEg = gammaCorrectColormap(CME,GAMMAE);
VE = cell(size(E));
limv = -inf;
for k = 1:numel(VE)
    limv = max(limv,max(E{k}(:)));
end
for k = 1:numel(VE)
    VE{k} = colors2im(E{k},CMEg,[0, ALPHAE * limv]);
end
VE = cell2mov(VE(2:end),1,255);

% Save movie frames
XX = cell2mov({VX, VE},1,255);

%{
saveColorFrames(XX,'otazo_errors/frame.png'); % otazo
saveColorFrames(XX,'pincat_errors/frame.png'); % PINCAT
%}

%% 3/6 Single-frame recons

% Find best frames
ourIdx = 1; % our method
comIdx = 2; % competitor method
[~, idx] = sort(err(comIdx).pNRMSE - err(ourIdx).pNRMSE,'descend');
best_frames = idx(1:10) %#ok

% Generate frames
nFrames = numel(FRAMES);
IX = cell(nFrames,nData + 1);
formatFrame = @(X,f) min(abs(X(:,:,f)) / trueMax,1);
for j = 1:nFrames
    IX(j,:) = cellfun(@(Xk) formatFrame(Xk,FRAMES(j)),X,'UniformOutput',false);
end

% Format frames
CMXg = gammaCorrectColormap(CMX,GAMMAX);
IX = cellfun(@(X) colors2im(X,CMXg,[0, 1]),IX,'UniformOutput',false);
IX = cell2mov(IX,GAP,VAL);

% Plot frames
imshow(IX,[0, 1]);
ylabels = arrayfun(@(f) sprintf('Frame %d',f),FRAMES,'UniformOutput',false);
addXLabels(xlabels,[],0.015,'top');
addYLabels(ylabels,[],0.005,'left');
SetFigFontSize(FONT_SIZE);

%{
cfigure([705, 369],gcf()); % 3 algos
savePDF('asilomar_otazo_recons_f4f12');
savePDF('asilomar_pincat_recons_f17f44');

cfigure([818, 369],gcf()); % 4 algos
savePDF('tci_otazo_recons_f4f11');
savePDF('tci_pincat_recons_f17f45');
%}

%% 4/6 Single-frame recons + errors

% Find best frames
ourIdx = 1; % our method
comIdx = 2; % competitor method
[~, idx] = sort(err(comIdx).pNRMSE - err(ourIdx).pNRMSE,'descend');
best_frames = idx(1:10) %#ok

% Generate frames
nFrames = numel(FRAMES);
IX = cell(nFrames,nData + 1);
formatFrame = @(X,f) min(abs(X(:,:,f)) / trueMax,1);
for j = 1:nFrames
    IX(j,:) = cellfun(@(Xk) formatFrame(Xk,FRAMES(j)),X,'UniformOutput',false);
end

% Generate errors
IE = cell(nFrames,nData + 1);
for j = 1:nFrames
    f = FRAMES(j);
    Xtruef = Xtrue(:,:,f);
    IE{j,1} = zeros(size(Xtruef));
    for k = 1:nData
        IE{j,k + 1} = abs(Xhat{k}(:,:,f) - Xtruef);
    end
end

% Extract relevant data
IX = IX(:,1:2); % Truth and online
IE = IE(:,2:end);
XLABELS = [xlabels(1:2), xlabels(2:end)];
YLABELS = arrayfun(@(f) sprintf('Frame %d',f),FRAMES,'UniformOutput',false);

% Format frames
CMXg = gammaCorrectColormap(CMX,GAMMAX);
IX = cellfun(@(X) colors2im(X,CMXg,[0, 1]),IX,'UniformOutput',false);
IX = cell2mov(IX,GAP,VAL);

% Format errors
CMEg = gammaCorrectColormap(CME,GAMMAE);
Elim = [0, ALPHAE * max(vec(cellfun(@(X) max(X(:)),IE)))];
IE = cellfun(@(E) colors2im(E,CMEg,Elim),IE,'UniformOutput',false);
IE = cell2mov(IE,GAP,VAL);

% Plot frames + errors
figure;
imshow(cell2mov({IX, IE},GAP,VAL),[0, 1]);
addXLabels(XLABELS,[],0.015,'top');
addYLabels(YLABELS,[],0.005,'left');

% Add colorbar #1
ax1 = gca();
ax2 = axes('Position',get(ax1,'Position'),'Visible','off');
set(colorbar('Peer',ax1),'Visible','off');
colorbar('Peer',ax2);
colormap(ax2,CMEg);
caxis(ax2,Elim);

%{
% Add colorbar #2
c = colorbar();
colormap(CMEg);
yTicks = get(c,'YTick');
yTicks2 = (Elim(2) / size(CMEg,1)) * yTicks;
yTickLabels = arrayfun(@(t) sprintf('%.2f',t),yTicks2,'UniformOutput',false);
set(c,'YTickLabel',yTickLabels);
%}

% Set font size
SetFigFontSize(FONT_SIZE);

%{
cfigure([886, 369],gcf()); % 3 algos
savePDF('asilomar_otazo_errors_f4f12');
savePDF('asilomar_pincat_errors_f17f44');

cfigure([1019, 369],gcf()); % 4 algos
savePDF('tci_otazo_errors_f4f11');
savePDF('tci_pincat_errors_f17f45');
%}

%% 5/6 Per-frame errors

% Setup figure
CM = linspecer(nData);
cfigure([468, 273]);
nt = size(Xtrue,3);
st = {'-','.-','*-','o-'}; st = repmat(st,1,3);
ms = [6, 12, 6, 6]; ms = repmat(ms,1,3);
lw = [1, 1, 1, 1]; lw = repmat(lw,1,3);

% Sort frames
[~, idx] = sort([err.NRMSE],'descend');

% Plot pNRMSEs
phndl = zeros(1,nData);
for k = 1:nData
    k1 = idx(k);
    k2 = nData + 1 - k;
    phndl(k) = plot(1:nt,100 * err(k1).pNRMSE,st{k2},'Color',CM(k2,:),'Linewidth',lw(k2),'Markersize',ms(k2));
    hold on;
end
xlabel('Frame');
ylabel('NRMSE (%)');
title('Per-frame NRMSEs','FontWeight','Normal');
legend(phndl,inpaths{idx,2},'Location','NE');
axis tight; padAxis();

% Adjust x-ticks
xTicks = get(gca,'XTick'); xTicks(1) = 1;
set(gca,'XTick',xTicks);

% Set font size
SetFigFontSize(FONT_SIZE);

%{
savePDF('asilomar_otazo_pnrmse');
savePDF('asilomar_pincat_pnrmse');

savePDF('tci_otazo_pnrmse');
savePDF('tci_pincat_pnrmse');
%}

%% 6/6 y-t profiles

%{
% To pick a y slice:
imshow(abs(Xtrue(:,:,11)),[0, 1]);
%}

% Construct y-t profiles
nt = size(Xtrue,3);
YT = cell(size(X));
for k = 1:numel(X)
    YT{k} = reshape(X{k}(ROI{2},ROI{1},:),[],nt);
    YT{k} = min(abs(YT{k}) / trueMax,1);
end

% Format frames
CMXg = gammaCorrectColormap(CMX,GAMMAX);
YT = cellfun(@(X) colors2im(X,CMXg,[0, 1]),YT,'UniformOutput',false);
YT = cell2mov(YT,GAP,VAL);

% Plot y-t profiles
figure();
imshow(YT,[0, 1]);
addXLabels(xlabels,[],0.015,'top');
SetFigFontSize(FONT_SIZE);

%{
cfigure([1136, 700],gcf()); % 3 algos
savePDF('asilomar_otazo_ytprofiles');
savePDF('asilomar_pincat_ytprofiles');

cfigure([1136, 700],gcf()); % 4 algos
savePDF('tci_otazo_ytprofiles');
savePDF('tci_pincat_ytprofiles');
%}

%% [TABLE,PLOT] road_interp_par* 

%{
% road396, SNR=inf
inpaths = {'road_interp_par1.mat', 'paint2, 1';
           'road_interp_par2.mat', 'grid2, cubic';
           'road_interp_par3.mat', 'tri3, linear';
           'road_interp_par4.mat', 'tri3, natrual'};
name = 'road_interp_128_psnr.pdf';
%}
%{
% road396, SNR=25
inpaths = {'road_interp_par1c.mat', 'paint2, 1';
           'road_interp_par2c.mat', 'grid2, cubic';
           'road_interp_par3c.mat', 'tri3, linear';
           'road_interp_par4c.mat', 'tri3, natrual'};
name = 'road_snr25_interp_128_psnr.pdf';
%}
%{
% road396, SNR=19
inpaths = {'road_interp_par1c2.mat', 'paint2, 1';
           'road_interp_par2c2.mat', 'grid2, cubic';
           'road_interp_par3c2.mat', 'tri3, linear';
           'road_interp_par4c2.mat', 'tri3, natrual'};
name = 'road_snr19_interp_128_psnr.pdf';
%}
%{
% road100, SNR=inf
inpaths = {'road_interp_par1d.mat', 'paint2, 1';
           'road_interp_par2d.mat', 'grid2, cubic';
           'road_interp_par3d.mat', 'tri3, linear';
           'road_interp_par4d.mat', 'tri3, natrual'};
name = 'road2_interp_128_psnr.pdf';
%}
%{
% road100, SNR=25
inpaths = {'road_interp_par1e.mat', 'paint2, 1';
           'road_interp_par2e.mat', 'grid2, cubic';
           'road_interp_par3e.mat', 'tri3, linear';
           'road_interp_par4e.mat', 'tri3, natrual'};
name = 'road2_snr25_interp_128_psnr.pdf';
%}
%{
% road100, SNR=19
inpaths = {'road_interp_par1e2.mat', 'paint2, 1';
           'road_interp_par2e2.mat', 'grid2, cubic';
           'road_interp_par3e2.mat', 'tri3, linear';
           'road_interp_par4e2.mat', 'tri3, natrual'};
name = 'road2_snr19_interp_128_psnr.pdf';
%}
pSortIdx = 3;
metric = 'pSNR'; % NRMSE, pSNR
SAVE_FIGURE = true;

% Load data
data = struct();
for i = 1:size(inpaths,1)
    datai = load(inpaths{i,1},'vars','err');
    data(i).p = squeeze(datai.vars.p);
    data(i).(metric) = squeeze(datai.err.(metric));
    data(i).label = inpaths{i,2};
end
nData = numel(data);

% Sort by metric
val = @(data) data.(metric)(pSortIdx);
[~, idx] = sort(arrayfun(val,data),'descend');

%--------------------------------------------------------------------------
% Print table
%--------------------------------------------------------------------------
% Headings
labels = struct();
labels.align = 'center';
labelFcn = @(p) sprintf('%d%%',100 * p);
labels.strs  = [metric, arrayfun(labelFcn,data(1).p,'UniformOutput',false)];

% Row data
rows = struct();
for i = 1:nData
    ii = idx(i);
    rows(i).name = data(ii).label;
    rows(i).data = data(ii).(metric);
end
[rows.fmt]    = deal('%.2f');
[rows.nAlign] = deal('right');
[rows.dAlign] = deal('right');

% Print table
printTable(labels,rows);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
cfigure([281, 231]);
cm = linspecer(nData);

% Plot metric
phndl = zeros(1,nData);
pstr = cell(1,nData);
for i = 1:nData
    ii = idx(i);
    phndl(i) = plot(data(ii).p,data(ii).(metric),'-o','Color',cm(i,:));
    pstr{i} = data(ii).label;
    hold on;
end
xlabel('p');
ylabel(metric);
legend(phndl,pstr{:},'Location','SW');
axis tight; padAxis();

% Save figure, if requested
if SAVE_FIGURE
    export_fig('-pdf','-transparent',name);
    fprintf('DONE\n');
end
%--------------------------------------------------------------------------

%% [PLOT] road_rpca_par*

% Knobs
%inpath = 'road_rpca_par1.mat';
%inpath = 'road_rpca_par1b.mat';
%inpath = 'road_rpca_par2.mat';
%inpath = 'road_rpca_par2b.mat';
%inpath = 'road_rpca_par3.mat';
inpath = 'road_rpca_par3b.mat';
%inpath = 'road_rpca_par4.mat';
%inpath = 'road_rpca_par4b.mat';
SAVE_FIGURE = true;

% Load data (Lhat, Shat, err, stats, vars)
load(inpath);
[~, name, ~] = fileparts(inpath);
[vars, ~] = eval(sprintf('%s();',name));

% Optimal parameters
NRMSE = nanmean(err.NRMSE,2); % average over seed
[mnrmse, idx] = min(NRMSE(:));
[ip, is, il, im] = ind2sub(size(NRMSE),idx);

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
cfigure([725, 343]);
ax = [];

% p
subplot(2,3,1);
if (numel(vars.p) > 1), ax(end + 1) = gca(); end
plot(vars.p,squeeze(NRMSE(:,is,il,im)),'b-o');
xlabel('p');
ylabel('NRMSE');
title(sprintf('p = %.3f',vars.p(ip)));
axis tight; padAxis();

% lambda
subplot(2,3,2);
if (numel(vars.lambda) > 1), ax(end + 1) = gca(); end
semilogx(vars.lambda,squeeze(NRMSE(ip,is,:,im)),'b-o');
xlabel('lambda');
ylabel('NRMSE');
title(sprintf('lambda = %.3f',vars.lambda(il)));
axis tight; padAxis();

% mu
subplot(2,3,3);
if (numel(vars.mu) > 1), ax(end + 1) = gca(); end
semilogx(vars.mu,squeeze(NRMSE(ip,is,il,:)),'b-o');
xlabel('mu');
ylabel('NRMSE');
title(sprintf('mu = %.3f',vars.mu(im)));
axis tight; padAxis();

% Optimal cost trajectory
subplot(2,3,4);
plot(1:vars.nIters,squeeze(stats.cost(ip,is,il,im,:)),'b-');
xlabel('iteration');
ylabel('cost');
title('cost');
axis tight; padAxis();

% Optimal NRMSE trajectory
subplot(2,3,5);
plot(1:vars.nIters,squeeze(stats.nrmse(ip,is,il,im,:)),'b-');
xlabel('iteration');
ylabel('NRMSE');
title(sprintf('NRMSE = %.5f',mnrmse));
axis tight; padAxis();

% Optimal delta trajectory
subplot(2,3,6);
plot(1:vars.nIters,squeeze(stats.delta(ip,is,il,im,:)),'b-');
xlabel('iteration');
ylabel('deltaX');
title('delta');
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

%% [PLOT] road_dls_par*

% Knobs
%inpath = 'road_dls_par1a.mat';
%inpath = 'road_dls_par1a2.mat';
SAVE_FIGURE = true;

% Load data (D, Xhat, err, stats, vars)
load(inpath);
[~, name, ~] = fileparts(inpath);
[vars, opts] = eval(sprintf('%s();',name));

% Optimal parameters
NRMSE = nanmean(err.NRMSE,2); % average over seed
[mnrmse, idx] = min(NRMSE(:));
[ip, is, il, im, ir] = ind2sub(size(NRMSE),idx);

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
cfigure([1190, 343]);
ax = [];

% p
subplot(2,5,1);
if (numel(vars.p) > 1), ax(end + 1) = gca(); end
plot(vars.p,squeeze(NRMSE(:,is,il,im,ir)),'b-o');
xlabel('p');
ylabel('NRMSE');
title(sprintf('p = %.3f',vars.p(ip)));
axis tight; padAxis();

% lambda
subplot(2,5,2);
if (numel(vars.lambda) > 1), ax(end + 1) = gca(); end
semilogx(vars.lambda,squeeze(NRMSE(ip,is,:,im,ir)),'b-o');
xlabel('lambda');
ylabel('NRMSE');
title(sprintf('lambda = %.3f',vars.lambda(il)));
axis tight; padAxis();

% mu
subplot(2,5,3);
if (numel(vars.mu) > 1), ax(end + 1) = gca(); end
semilogx(vars.mu,squeeze(NRMSE(ip,is,il,:,ir)),'b-o');
xlabel('mu');
ylabel('NRMSE');
title(sprintf('mu = %.3f',vars.mu(im)));
axis tight; padAxis();

% dr
subplot(2,5,4);
if (numel(vars.dr) > 1), ax(end + 1) = gca(); end
plot(vars.dr,squeeze(NRMSE(ip,is,il,im,:)),'b-o');
xlabel('dr');
ylabel('NRMSE');
title(sprintf('dr = %d',vars.dr(ir)));
axis tight; padAxis();

% Optimal cost trajectory
subplot(2,5,6);
plot(1:vars.nIters,squeeze(stats.cost(ip,is,il,im,ir,:)),'b-');
xlabel('iteration');
ylabel('cost');
title(sprintf('cost, (np) = (%d)',vars.np));
axis tight; padAxis();

% Optimal NRMSE trajectory
subplot(2,5,7);
plot(1:vars.nIters,squeeze(stats.nrmse(ip,is,il,im,ir,:)),'b-');
xlabel('iteration');
ylabel('NRMSE');
title(sprintf('NRMSE = %.5f',mnrmse));
axis tight; padAxis();

% Optimal sparsity trajectory
subplot(2,5,8);
plot(1:vars.nIters,squeeze(stats.sparsity(ip,is,il,im,ir,:)),'b-');
xlabel('iteration');
ylabel('sparsity (%)');
title(sprintf('sparsity, idim = [%d, %d, %d]',vars.idim));
axis tight; padAxis();

% Optimal deltaX trajectory
subplot(2,5,9);
plot(1:vars.nIters,squeeze(stats.deltaX(ip,is,il,im,ir,:)),'b-');
xlabel('iteration');
ylabel('deltaX');
title(sprintf('deltaX, pgap = [%d, %d, %d]',opts.pgap));
set(gca,'YScale','log');
axis tight; padAxis();

% Optimal deltaD trajectory
subplot(2,5,10);
plot(1:vars.nIters,squeeze(stats.deltaD(ip,is,il,im,ir,:)),'b-');
xlabel('iteration');
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

%% [PLOT] road_onlineDls_par*

% Knobs
%inpath = 'road_onlineDls_par1.mat';
%inpath = 'road_onlineDls_par1b.mat';
%inpath = 'road_onlineDls_par2.mat';
%inpath = 'road_onlineDls_par2b.mat';
%inpath = 'road_onlineDls_par2c.mat';
%inpath = 'road_onlineDls_par3.mat';
%inpath = 'road_onlineDls_par4.mat';
%inpath = 'road_onlineDls_par5a.mat';
%inpath = 'road_onlineDls_par5a2.mat';
%inpath = 'road_onlineDls_par5b.mat';
%inpath = 'road_onlineDls_par5b2.mat';
%inpath = 'road_onlineDls_par6.mat';
%inpath = 'road_onlineDls_par7a.mat';
%inpath = 'road_onlineDls_par7a2.mat';
%inpath = 'road_onlineDls_par7b.mat';
%inpath = 'road_onlineDls_par7b2.mat';
%inpath = 'road_onlineDls_par7c.mat';
%inpath = 'road_onlineDls_par8.mat';
%inpath = 'road_onlineDls_par9.mat';
%inpath = 'road_onlineDls_par9b.mat';
%inpath = 'road_onlineDls_par9c.mat';
%inpath = 'road_onlineDls_par9d.mat';
%inpath = 'road_onlineDls_par10.mat';
%inpath = 'road_onlineDls_par10b.mat';
%inpath = 'road_onlineDls_par11.mat';
SAVE_FIGURE = true;

% Load data (D, Xhat, err, stats, vars)
load(inpath);
[~, name, ~] = fileparts(inpath);
[vars, opts] = eval(sprintf('%s();',name));

% Optimal parameters
NRMSE = nanmean(err.NRMSE,2); % average over seed
[mnrmse, idx] = min(NRMSE(:));
[ip, is, il, im, ir, ig] = ind2sub(size(NRMSE),idx);

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
cfigure([1190, 343]);
ax = [];

% p
subplot(2,5,1);
if (numel(vars.p) > 1), ax(end + 1) = gca(); end
plot(vars.p,squeeze(NRMSE(:,is,il,im,ir,ig)),'b-o');
xlabel('p');
ylabel('NRMSE');
title(sprintf('p = %.3f',vars.p(ip)));
axis tight; padAxis();

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
if (numel(vars.gamma) > 1), ax(end + 1) = gca(); end
plot(vars.gamma,squeeze(NRMSE(ip,is,il,im,ir,:)),'b-o'); hold on;
plot(vars.gamma2,squeeze(NRMSE(ip,is,il,im,ir,:)),'r-x'); hold on;
xlabel('gamma, gamma2');
ylabel('NRMSE');
title(sprintf('gamma = %.3f, gamma2 = %.3f',vars.gamma(ig),vars.gamma2(ig)));
axis tight; padAxis();

% Optimal cost trajectory
subplot(2,5,6);
plotOnlineDlsMetric(squeeze(stats.cost(ip,is,il,im,ir,ig,:)),vars.nItersi,vars.nIters,vars.dt,'b-');
xlabel('time + iteration');
ylabel('cost');
title(sprintf('cost, (T, dt, np) = (%d, %d, %d)',vars.T,vars.dt,vars.np));
axis tight; padAxis();

% Optimal NRMSE trajectory
subplot(2,5,7);
plotOnlineDlsMetric(squeeze(stats.nrmse(ip,is,il,im,ir,ig,:)),vars.nItersi,vars.nIters,vars.dt,'b-');
xlabel('time + iteration');
ylabel('NRMSE');
title(sprintf('NRMSE = %.5f',mnrmse));
axis tight; padAxis();

% Optimal sparsity trajectory
subplot(2,5,8);
plotOnlineDlsMetric(squeeze(stats.sparsity(ip,is,il,im,ir,ig,:)),vars.nItersi,vars.nIters,vars.dt,'b-');
xlabel('time + iteration');
ylabel('sparsity (%)');
title(sprintf('sparsity, idim = [%d, %d, %d]',vars.idim));
axis tight; padAxis();

% Optimal deltaX trajectory
subplot(2,5,9);
plotOnlineDlsMetric(squeeze(stats.deltaX(ip,is,il,im,ir,ig,:)),vars.nItersi,vars.nIters,vars.dt,'b-');
xlabel('time + iteration');
ylabel('deltaX');
title(sprintf('deltaX, pgap = [%d, %d, %d]',opts.pgap));
set(gca,'YScale','log');
axis tight; padAxis();

% Optimal deltaD trajectory
subplot(2,5,10);
plotOnlineDlsMetric(squeeze(stats.deltaD(ip,is,il,im,ir,ig,:)),vars.nItersi,vars.nIters,vars.dt,'b-');
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

%% [PLOT] road_onlineDls_atoms_par*

% Knobs
%inpath = 'road_onlineDls_atoms_par1.mat';
%inpath = 'road_onlineDls_atoms_par2.mat';
SAVE_FIGURE = true;

% Load data (D, Xhat, err, stats, vars)
load(inpath);

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
cfigure([367, 190]);

% np
plot(vars.np,err.NRMSE,'b-o');
xlabel('np');
ylabel('NRMSE');
title(sprintf('NRMSE [%.5f]',min(err.NRMSE)));
axis tight; padAxis();

% Save figure, if requested
if SAVE_FIGURE
    [~, name, ~] = fileparts(inpath);
    export_fig('-pdf','-transparent',name);
    fprintf('DONE\n');
end
%--------------------------------------------------------------------------

%% [TABLE,PLOT] p sweeps
% Plots (p,dr) slices of data

% road396, SNR = inf
ipaths = {'road_interp_par2.mat', 'cubic (2D)';
          'road_interp_par4.mat', 'cubic (3D)'};
%rpaths = {'road_rpca_par2.mat', 'RPCA'};
%rpaths = {'road_rpca_par2b.mat', 'RPCA'};
dpaths = {'road_onlineDls_par12a.mat', 'online DINO-KAT';
          'road_onlineDls_par12b.mat', 'online DINO-KAT [fixed D]';
          'road_dls_par2a.mat', 'DINO-KAT';
          'road_dls_par2b.mat', 'DINO-KAT [fixed D]'};
name = 'road396_psweep_SNRinf_psnr.pdf';

%{
% road396, SNR = 19
ipaths = {'road_interp_par2c2.mat', 'cubic (2D)';
          'road_interp_par4c2.mat', 'cubic (3D)'};
%rpaths = {'road_rpca_par4.mat', 'RPCA'};
%rpaths = {'road_rpca_par4b.mat', 'RPCA'};
dpaths = {'road_onlineDls_par12c2.mat', 'online DINO-KAT';
          'road_onlineDls_par12d2.mat', 'online DINO-KAT [fixed D]';
          'road_dls_par2c.mat', 'DINO-KAT';
          'road_dls_par2d.mat', 'DINO-KAT [fixed D]'};
name = 'road396_psweep_SNR19_psnr.pdf';
%}
%{
% road100, SNR = inf
ipaths = {'road_interp_par2d.mat', 'cubic (2D)';
          'road_interp_par4d.mat', 'cubic (3D)'};
dpaths = {'road_onlineDls_par13a.mat', 'online DINO-KAT';
          'road_onlineDls_par13b.mat', 'online DINO-KAT [fixed D]';
          'road_dls_par3a.mat', 'DINO-KAT';
          'road_dls_par3b.mat', 'DINO-KAT [fixed D]'};
name = 'road100_psweep_SNRinf_psnr.pdf';
%}
%{
% road100, SNR = 19
ipaths = {'road_interp_par2e2.mat', 'cubic (2D)';
          'road_interp_par4e2.mat', 'cubic (3D)'};
dpaths = {'road_onlineDls_par13c2.mat', 'online DINO-KAT';
          'road_onlineDls_par13d2.mat', 'online DINO-KAT [fixed D]';
          'road_dls_par3c.mat', 'DINO-KAT';
          'road_dls_par3d.mat', 'DINO-KAT [fixed D]'};
name = 'road100_psweep_SNR19_psnr.pdf';
%}

% road396, SNR = inf [OLD]
%{
ipaths = {'road_interp_par2.mat', 'cubic (2D)';
          'road_interp_par4.mat', 'cubic (3D)'};
%}
%{
% DINO-KAT [8,8,5] (paint2,1)
dpaths = {'road_onlineDls_par5a.mat', 'online DINO-KAT';
          'road_onlineDls_par5b.mat', 'online DINO-KAT [fixed D]';
          'road_dls_par1a.mat'      , 'DINO-KAT'};
name = 'road_psweep_885_psnr.pdf';
%}
%{
% DINO-KAT [8,8,8] (paint2,1)
dpaths = {'road_onlineDls_par7a.mat', 'online DINO-KAT';
          'road_onlineDls_par7b.mat', 'online DINO-KAT [fixed D]'};
name = 'road_psweep_888_psnr.pdf';
%}
%{
% DINO-KAT [8,8,5] (tri3,natural)
dpaths = {'road_onlineDls_par5a2.mat', 'online DINO-KAT';
          'road_onlineDls_par5b2.mat', 'online DINO-KAT [fixed D]',
          'road_dls_par1a2.mat'      , 'DINO-KAT'};
name = 'road_psweep2_885_psnr.pdf';
%}
%{
% DINO-KAT [8,8,8] (tri3,natural)
dpaths = {'road_onlineDls_par7a2.mat', 'online DINO-KAT';
          'road_onlineDls_par7b2.mat', 'online DINO-KAT [fixed D]'};
name = 'road_psweep2_888_psnr.pdf';
%}

pSortIdx = 3;
metric = 'pSNR'; % NRMSE, pSNR
SAVE_FIGURE = true;

% Load data
data = struct();
for i = 1:size(ipaths,1)
    % Interp data
    datai = load(ipaths{i,1},'vars','err');
    data(i).p = squeeze(datai.vars.p);
    data(i).(metric) = squeeze(datai.err.(metric));
    data(i).label = ipaths{i,2};
end
if exist('rpaths','var')
    for j = 1:size(rpaths,1)
        % RPCA data
        dataj = load(rpaths{j,1},'vars','err');
        data(end + 1).p = squeeze(dataj.vars.p); %#ok
        data(end).(metric) = squeeze(dataj.err.(metric));
        data(end).label = rpaths{j,2};
    end
end
for j = 1:size(dpaths,1)
    % DINO-KAT data
    dataj = load(dpaths{j,1},'vars','err');
    errData = dataj.err.(metric); % 7D
    pData = squeeze(dataj.vars.p);
    drData = dataj.vars.dr;
    for k = 1:size(errData,5)
        data(end + 1).p = pData; %#ok
        data(end).(metric) = squeeze(errData(:,:,:,:,k,:,:));
        if isnan(drData(k))
            % No dr
            data(end).label = sprintf('%s',dpaths{j,2});
        else
            % Add dr
            data(end).label = sprintf('%s [dr = %d]',dpaths{j,2},drData(k));
        end
    end
end
nData = numel(data);

% Sort by metric
val = @(data) data.(metric)(pSortIdx);
[~, idx] = sort(arrayfun(val,data),'descend');

%--------------------------------------------------------------------------
% Print table
%--------------------------------------------------------------------------
% Headings
labels = struct();
labels.align = 'center';
labelFcn = @(p) sprintf('%d%%',100 * p);
labels.strs  = [metric, arrayfun(labelFcn,data(1).p,'UniformOutput',false)];

% Row data
rows = struct();
for i = 1:nData
    ii = idx(i);
    rows(i).name = data(ii).label;
    rows(i).data = data(ii).(metric);
end
[rows.fmt]    = deal('%.2f');
[rows.nAlign] = deal('right');
[rows.dAlign] = deal('right');

% Print table
printTable(labels,rows);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
cfigure([492, 231]);
cm = linspecer(nData);

% Plot metric
phndl = zeros(1,nData);
pstr = cell(1,nData);
for i = 1:nData
    ii = idx(i);
    phndl(i) = plot(data(ii).p,data(ii).(metric),'-o','Color',cm(i,:));
    pstr{i} = data(ii).label;
    hold on;
end
xlabel('p');
ylabel(metric);
legend(phndl,pstr{:},'Location','EastOutside');
axis tight; padAxis();

% Save figure, if requested
if SAVE_FIGURE
    export_fig('-pdf','-transparent',name);
    fprintf('DONE\n');
end
%--------------------------------------------------------------------------

%% [VIDEO,PLOT] movies, error maps, per-frame pSNRs

%--------------------------------------------------------------------------
% Knobs
%--------------------------------------------------------------------------
%dr = [1, 5]; p = [0.5]; % for per-frame SNR figures
dr = [1]; p = [0.6, 0.7]; % for single frame figures
%dr = [1]; p = [0.6]; % for movies

BASE = '/Users/Brian/Desktop/Online DINO-KAT/ICIP';
SAVE_FIGURE = true;
FONT_SIZE = 12;
LINEWIDTH = 1.5;

%{
% road396, SNR = inf
%ONLINE_DLS = 'road_onlineDls_par5a';
ONLINE_DLS = 'road_onlineDls_par12a';
FIXED_DLS = 'road_onlineDls_par12b';
%BATCH_DLS = 'road_dls_par2a';
%RPCA = 'road_rpca_par2';
%RPCA = 'road_rpca_par2b';
CUBIC_2D = 'road_interp_par2';
CUBIC_3D = 'road_interp_par4';
BASE_OUT_NAME = 'road396_SNRinf.pdf';
YLIM = [25.3, 30.65]; % HACK
%}
%{
% road396, SNR = 19
ONLINE_DLS = 'road_onlineDls_par12c2';
FIXED_DLS = 'road_onlineDls_par12d2';
%BATCH_DLS = 'road_dls_par2c';
%RPCA = 'road_rpca_par4';
%RPCA = 'road_rpca_par4b';
CUBIC_2D = 'road_interp_par2c2';
CUBIC_3D = 'road_interp_par4c2';
BASE_OUT_NAME = 'road396_SNR19.pdf';
YLIM = [21.2, 26.9];
%}

%{
% road100, SNR = inf
ONLINE_DLS = 'road_onlineDls_par13a';
%BATCH_DLS = 'road_dls_par3a';
CUBIC_2D = 'road_interp_par2d';
CUBIC_3D = 'road_interp_par4d';
BASE_OUT_NAME = 'road100_SNRinf.pdf';
%}
%{
% road100, SNR = 19
ONLINE_DLS = 'road_onlineDls_par13c2';
%BATCH_DLS = 'road_dls_par3c';
CUBIC_2D = 'road_interp_par2e2';
CUBIC_3D = 'road_interp_par4e2';
BASE_OUT_NAME = 'road100_SNR19.pdf';
%}
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Load data
%--------------------------------------------------------------------------
np = numel(p);
nr = numel(dr);
data = cell(np,0);
xlabels = cell(1,0);
ylabels = cell(1,0);
genPath = @(path) sprintf('%s/%s',BASE,path);

% Load vars
[vars, ~] = eval(sprintf('%s();',ONLINE_DLS));
if ~isfield(vars,'SNR'), vars.SNR = inf; end
na = 0;

% Load (Xture, Y)
for i = 1:np
    [Xtruei, Yi, ~] = loadData(vars.inpath,vars.idim,vars.T,vars.dt,p(i),vars.seed,vars.SNR);
    data{i,na + 1} = Xtruei;
    data{i,na + 2} = Yi;
end
xlabels{na + 1} = 'Truth';
xlabels{na + 2} = 'Corrupted';
na = na + 2;

% Load online DINO-KAT
for j = 1:nr
    drIdxj = find(ismember(vars.dr,dr(j)));
    for i = 1:np
        pIdxi = find(ismember(vars.p,p(i)));
        idxi = pIdxi + numel(vars.p) * (drIdxj - 1);
        tmp = load(genPath(sprintf('%s/data%d.mat',ONLINE_DLS,idxi)));
        data{i,na + 1} = tmp.Xhat;
    end
    %xlabels{na + 1} = sprintf('online DINO-KAT [r = %d]',dr(j));
    xlabels{na + 1} = sprintf('Online (r = %d)',dr(j));
    na = na + 1;
end

% Load online DINO-KAT fixed D
if exist('FIXED_DLS','var') && ~isempty(FIXED_DLS)
    for i = 1:np
        pIdxi = find(ismember(vars.p,p(i)));
        idxi = pIdxi;
        tmp = load(genPath(sprintf('%s/data%d.mat',FIXED_DLS,idxi)));
        data{i,na + 1} = tmp.Xhat;
    end
    xlabels{na + 1} = 'Online (DCT D)';
    na = na + 1;
end

% Load DINO-KAT
if exist('BATCH_DLS','var') && ~isempty(BATCH_DLS)
    for j = 1:nr
        drIdxj = find(ismember(vars.dr,dr(j)));
        for i = 1:np
            pIdxi = find(ismember(vars.p,p(i)));
            idxi = pIdxi + numel(vars.p) * (drIdxj - 1);
            tmp = load(genPath(sprintf('%s/data%d.mat',BATCH_DLS,idxi)));
            data{i,na + 1} = tmp.Xhat;
        end
        %xlabels{na + 1} = sprintf('DINO-KAT [r = %d]',dr(j));
        xlabels{na + 1} = sprintf('(r = %d)',dr(j));
        na = na + 1;
    end
end

% Load RPCA
if exist('RPCA','var') && ~isempty(RPCA)
    for i = 1:np
        pIdxi = find(ismember(vars.p,p(i)));
        idxi = pIdxi;
        tmp = load(genPath(sprintf('%s/data%d.mat',RPCA,idxi)));
        data{i,na + 1} = tmp.Lhat + tmp.Shat;
    end
    xlabels{na + 1} = 'RPCA';
    na = na + 1;
end

% Load interp
for i = 1:np
    pIdxi = find(ismember(vars.p,p(i)));
    idxi = pIdxi;
    tmp1 = load(genPath(sprintf('%s/data%d.mat',CUBIC_3D,idxi)));
    tmp2 = load(genPath(sprintf('%s/data%d.mat',CUBIC_2D,idxi)));
    data{i,na + 1} = tmp1.Xhat;
    data{i,na + 2} = tmp2.Xhat;
end
xlabels{na + 1} = 'Interp (3D)';
xlabels{na + 2} = 'Interp (2D)';
na = na + 2;

% y-labels
for i = 1:np
    ylabels{i} = sprintf('%.0f%%',100 * p(i));
end
%--------------------------------------------------------------------------

% Compute errors
pSNR = nan(np,na);
ppSNR = cell(np,na);
for i = 1:np
    for j = 2:na
        errij = computeErrorMetrics(data{i,j},data{i,1});
        pSNR(i,j) = errij.pSNR;
        ppSNR{i,j} = errij.ppSNR;
    end
end
pSNR = pSNR %#ok

%{
%--------------------------------------------------------------------------
% Single-frame images
%--------------------------------------------------------------------------
% Knobs
gap = 1;
val = 1;
topFrameIdx = 1;
pVal = 0.6;
gamma = 1.0;
incl2d = false;

% Find best frames
dc = 1.0 - incl2d;
pIdx = find(ismember(p,pVal));
%[~, order] = sort(ppSNR{pIdx,3} + (ppSNR{pIdx,3} - ppSNR{pIdx,end - dc}),'descend');
[~, order] = sort(ppSNR{pIdx,3} + (ppSNR{pIdx,3} - ppSNR{pIdx,4}),'descend');
TOP_5_FRAMES = order(1:5) %#ok
fIdx = order(topFrameIdx);

% Generate image
C = cellfun(@(D) D(:,:,fIdx),data(pIdx,1:(end - dc)),'UniformOutput',false);
improvf = ppSNR{pIdx,3}(fIdx) - ppSNR{pIdx,4}(fIdx) %#ok
CC = cell2mov(C,gap,val);

% Gamma correction
if gamma ~= 1
    cmg = gammaCorrectColormap(gray(256),gamma);
    CC = colors2im(CC,cmg,[0, 1]);
end

% Plot image
imshow(CC,[0, 1]);
xLabels = xlabels(1:(end - dc));
xLabels{2} = sprintf('Corrupted (%.0f%%)',100 * pVal);
addXLabels(xLabels,[],0.03,'top');
addYLabels({sprintf('Frame %d',fIdx)},[],0.005,'left');
SetFigFontSize(FONT_SIZE);

% Save figure
if SAVE_FIGURE
    [~, name, ~] = fileparts(BASE_OUT_NAME);
    name = sprintf('%s_f%d_p%.0f.pdf',name,fIdx,100 * pVal);
    export_fig('-pdf','-transparent',name);
    fprintf('%s WRITTEN\n',name);
end
%--------------------------------------------------------------------------
%}

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recon + error movies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Knobs
cm = parula(64);
alpha = 0.80;  % Alpha < 1 to clip large values
gamma = 1.20;  % Gamma < 1 to decrease contrast

% Compute error videos
E = cell(np,na - 2);
for i = 1:np
    for j = 3:na
        E{i,j - 2} = abs(data{i,j} - data{i,1});
    end
end

% Format recon videos
M = cell2mat(data);
M = colors2im(M,gray(256),[0, 1]);

% Format error videos
cmg = gammaCorrectColormap(cm,gamma);
E = cell2mat(E);
E = colors2im(E,cmg,[0, alpha * max(E(:))]);

% For ground truth-only
M1 = cell2mat(data(:,1:2));
M1 = colors2im(M1,gray(256),[0, 1]);

video = cat(1,M,cat(2,M1,E));

% Visualize results
movie.video = video;
opts.xlabels  = xlabels;
opts.ylabels  = ylabels(end:-1:1);
opts.fontSize = FONT_SIZE;
opts.mag = 0.5;
PlayMovie(movie,opts);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recon movies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
movie.video = max(0,min(cell2mat(data),1));
opts.xlabels  = xlabels;
opts.ylabels  = ylabels(end:-1:1);
opts.fontSize = FONT_SIZE;
PlayMovie(movie,opts);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

%{
%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
cm = [zeros(2,3); linspecer(na - 2)];
nt = size(data{1},3);

for i = 1:np
    %cfigure([703, 296]);
    %cfigure([519, 274]);
    cfigure([597, 232]);
    pi = 100 * p(i);
    
    % Plot pSNRs
    phndl = nan(1,5);
    for j = na:-1:3
        phndl(j) = plot(1:nt,ppSNR{i,j},'-','Color',cm(j,:),'Linewidth',LINEWIDTH);
        hold on;
    end
    %{
    for j = 3:na
        plot([1, nt],pSNR(i,j) * [1, 1],'--','Color',cm(j,:),'Linewidth',LINEWIDTH);
    end
    %}
    xlabel('Frame');
    ylabel('PSNR (dB)');
    %title(sprintf('Per-frame PSNR [p = %.0f%%]',pi));
    legend(phndl(3:na),xlabels{3:na},'Location','Best');
    axis tight; padAxis();
    SetFigFontSize(FONT_SIZE);
    
    % HACKs
    xTick = get(gca,'XTick'); xTick(1) = 1; set(gca,'XTick',xTick);
    if exist('YLIM','var'), set(gca,'YLim',YLIM); end % HACK
    if exist('XTICK','var'), set(gca,'XTick',XTICK); end % HACK
    
    % Save results
    if SAVE_FIGURE
        [~, baseName, ext] = fileparts(BASE_OUT_NAME);
        name = sprintf('%s_psnrs_p%.0f%s',baseName,pi,ext);
        export_fig('-pdf','-transparent',name);
    end
end
%--------------------------------------------------------------------------
%}

%% [SPOT CHECK] recons

BASE = '/Users/Brian/Desktop/Online DINO-KAT/ICIP';

%{
% road396, SNR = inf
DINOKAT_PATH = 'road_onlineDls_par2c';
CUBIC2D_PATH = 'road_interp_par1'; % paint2,1
CUBIC3D_PATH = 'road_interp_par4'; % tri3,natural
%}

%{
% road396, SNR = 19dB
DINOKAT_PATH = 'road_onlineDls_par9d';
CUBIC2D_PATH = 'road_interp_par1c2'; % paint2,1
CUBIC3D_PATH = 'road_interp_par4c2'; % tri3,natural
%}

%{
% road100, SNR = inf
DINOKAT_PATH = 'road_onlineDls_par10b';
CUBIC2D_PATH = 'road_interp_par1d'; % paint2,1
CUBIC3D_PATH = 'road_interp_par4d'; % tri3,natural
%}

%{
% road100, SNR = 19dB
DINOKAT_PATH = 'road_onlineDls_par11';
CUBIC2D_PATH = 'road_interp_par1e2'; % paint2,1
CUBIC3D_PATH = 'road_interp_par4e2'; % tri3,natural
%}

% online DINO-KAT
datad = load(sprintf('%s.mat',DINOKAT_PATH));
Xhat = datad.Xhat;

% Ground truth
vars = datad.vars;
if ~isfield(vars,'SNR'), vars.SNR = inf; end
[Xtrue, Y, M] = loadData(vars.inpath,vars.idim,vars.T,vars.dt,vars.p,vars.seed,vars.SNR);

% cubic 2D
data2 = load(sprintf('%s/%s/data3.mat',BASE,CUBIC2D_PATH));
Xi2 = data2.Xhat;

% cubic 3D
data3 = load(sprintf('%s/%s/data3.mat',BASE,CUBIC3D_PATH));
Xi3 = data3.Xhat;

% Compute metrics
err_DINO = computeErrorMetrics(Xhat,Xtrue) %#ok
err_CUBIC2D = computeErrorMetrics(Xi2,Xtrue) %#ok
err_CUBIC3D = computeErrorMetrics(Xi3,Xtrue) %#ok

%{
% Spot check 2D interp
Xi = interpVideo(Y,M,'paint2',1);
err_EXTRA = computeErrorMetrics(Xi,Xtrue) %#ok
%}

%% [RUNTIMES]

% road396, SNR = inf, Online (r = 1)
%
% RESULT: ~2.5s per outer iteration
%
run = @road_onlineDls_run;
parFcn = @road_onlineDls_par12a;
idx = 2;

%{
% road396, SNR = inf, Online (DCT D)
%
% RESULT: ~1.6s per outer iteration
%
run = @road_onlineDls_run;
parFcn = @road_onlineDls_par12b;
idx = 2;
%}

%{
% road396, SNR = 19, Online (r = 1)
%
% RESULT: ~2.4s per outer iteration
%
run = @road_onlineDls_run;
parFcn = @road_onlineDls_par12c2;
idx = 3;
%}

%{
% road396, SNR = 19, Online (DCT D)
%
% RESULT: ~1.75s per outer iteration
%
run = @road_onlineDls_run;
parFcn = @road_onlineDls_par12d2;
idx = 3;
%}

% Run results
run(parFcn,idx);

%%
% Compare dls, dls2, and dls3
%
% Brian Moore
% brimoor@umich.edu
%

rng(42);

% Data knobs
path = '/Users/Brian/Archive/MATLAB/video_data/road.mat';
idim = [128, nan, 1];               % Raw video dimensions
pdim = [8, 8, 1];                   % Patch sizes
pgap = [1, 1, 1];                   % Patch strides
PROB = 0.50;                        % Missing data probability

% Algorithm knobs
lambda    = 0.01;                   % Dictionary regularization parameter
mu        = 0.05;                   % Sparsity regularization
D0        = dctmtx(prod(pdim))';    % Initial dictionary
type      = 'hard';                 % Sparsity regularizer
dr        = nan;                    % Atom rank constraint
ddim      = [8, 8];                 % Atom dimensions
fixedD    = false;                  % Fixed dictionary?
flag      = 1;                      % Print stats?
nIters    = 20;                     % # outer DLS iterations
% vvv  dls1 only  vvv
nItersDB1 = 01;                     % # (D,B) updates per outer iteration
% vvv  dls2 only  vvv
nItersDB2 = 20;                     % # (D,B) updates per outer iteration
accelDB2  = true;                   % Nesterov (D,B) updates?
tauDB2    = nan;                    % Relative (D,B) step size 
% vvv  dls3 only  vvv
nItersD3  = 20;                     % # inner D updates
accelD3   = true;                   % Nesterov D updates?
tauD3     = nan;                    % Relative D step size
nItersB3  = 20;                     % # inner B updates
accelB3   = true;                   % Nesterov B updates?
tauB3     = nan;                    % Relative B step size

% Load video
Xtrue = loadVideo(path,idim);
[ny, nx, nt] = size(Xtrue);
[p, m] = size(D0);

% Missing data
M = (rand(ny,nx,nt) > PROB);
Y = Xtrue;
Y(~M) = 0;

% Interpolate
Yi = interpVideo(Y,M);
nrmsei = norm(Yi(:) - Xtrue(:)) / norm(Xtrue(:)) %#ok

% dls
opts1 = struct();
opts1.A        = 1;
opts1.M        = M;
opts1.xdim     = [ny, nx, nt];
opts1.pdim     = pdim;
opts1.pgap     = pgap;
opts1.type     = type;
opts1.dr       = dr;
opts1.ddim     = ddim;
opts1.fixedD   = fixedD;
opts1.nIters   = nIters;
opts1.nItersDB = nItersDB1;
opts1.nItersX  = -1;
opts1.D0       = D0;
opts1.Xtrue    = Xtrue;
opts1.flag     = flag;
[Xhat1, D1, B1, stats1] = dls(Yi,lambda,mu,opts1);

% dls2
opts2 = struct();
opts2.A        = 1;
opts2.M        = M;
opts2.xdim     = [ny, nx, nt];
opts2.pdim     = pdim;
opts2.pgap     = pgap;
opts2.type     = type;
opts2.dr       = dr;
opts2.ddim     = ddim;
opts2.fixedD   = fixedD;
opts2.nIters   = nIters;
opts2.nItersDB = nItersDB2;
opts2.nItersX  = -1;
opts2.D0       = D0;
opts2.Xtrue    = Xtrue;
opts2.accel    = accelDB2;
opts2.tau      = tauDB2;
opts2.flag     = flag;
[Xhat2, D2, B2, stats2] = dls2(Yi,lambda,mu,opts2);

% dls3
opts3 = struct();
opts3.A       = 1;
opts3.M       = M;
opts3.xdim    = [ny, nx, nt];
opts3.pdim    = pdim;
opts3.pgap    = pgap;
opts3.type    = type;
opts3.dr      = dr;
opts3.ddim    = ddim;
opts3.fixedD  = fixedD;
opts3.nIters  = nIters;
opts3.nItersD = nItersD3;
opts3.nItersB = nItersB3;
opts3.nItersX = -1;
opts3.D0      = D0;
opts3.Xtrue   = Xtrue;
opts3.accelD  = accelD3;
opts3.accelB  = accelB3;
opts3.tauD    = tauD3;
opts3.tauB    = tauB3;
opts3.flag    = flag;
[Xhat3, D3, B3, stats3] = dls3(Yi,lambda,mu,opts3);

%{
%--------------------------------------------------------------------------
% Visualize results
%--------------------------------------------------------------------------
movie.video = cell2mat({Xtrue,Y,Yi;Xhat1,Xhat2,Xhat3});
opts.xlabels = {'Truth','Observations','Interpolated';
                'dls','dls2','dls3'};
PlayMovie(movie,opts);
%--------------------------------------------------------------------------
%}

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
cfigure([890, 222]);
cm = linspecer(3);

% cost
subplot(1,3,1);
phndl = zeros(1,3);
phndl(1) = plot(1:nIters,stats1.cost,'-o','Color',cm(1,:)); hold on;
phndl(2) = plot(1:nIters,stats2.cost,'-o','Color',cm(2,:));
phndl(3) = plot(1:nIters,stats3.cost,'-o','Color',cm(3,:));
xlabel('iteration');
title('cost');
legend(phndl,'dls','dls2','dls3');
axis tight; padAxis();

% nrmse
subplot(1,3,2);
plot(1:nIters,stats1.nrmse,'-o','Color',cm(1,:)); hold on;
plot(1:nIters,stats2.nrmse,'-o','Color',cm(2,:));
plot(1:nIters,stats3.nrmse,'-o','Color',cm(3,:));
plot([1, nIters],nrmsei * [1, 1],'k--');
xlabel('iteration');
title('nrmse');
axis tight; padAxis();

% time
subplot(1,3,3);
plot(1:nIters,cumsum(stats1.time),'-o','Color',cm(1,:)); hold on;
plot(1:nIters,cumsum(stats2.time),'-o','Color',cm(2,:));
plot(1:nIters,cumsum(stats3.time),'-o','Color',cm(3,:));
xlabel('iteration');
title('time');
axis tight; padAxis();
%--------------------------------------------------------------------------

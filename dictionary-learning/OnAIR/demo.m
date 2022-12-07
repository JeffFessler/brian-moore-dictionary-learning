%%
% OnAIR demo
%
% Brian Moore
% brimoor@umich.edu
%

rng(42);

% Data knobs
%path = '/Users/Brian/Archive/MATLAB/video_data/road.mat';
path = '/Users/Brian/Archive/MATLAB/video_data/road396.mat';
idim = [128, nan, 200];             % Raw video dimensions
T    = 5;                           % Streaming size
dt   = 1;                           % Temporal patch stride
pdim = [8, 8, T];                   % Patch sizes
pgap = [2, 2, dt];                  % Patch strides
p    = 0.50;                        % Missing data probability

% OnAIR knobs
lambda   = 0.01;                    % Dictionary regularization parameter
mu       = 0.05;                    % Sparsity regularization
gamma    = 0.9;                     % Forgetting factor (online)
D0       = dctmtx(prod(pdim))';     % Initial dictionary
type     = 'hard';                  % Sparsity regularizer
dr       = 5;                       % Atom rank constraint
ddim     = [prod(pdim(1:2)), T];    % Atom dimensions
unitaryD = false;                   % Unitary dictionary?
fixedD   = false;                   % Fixed dictionary?
flag     = 1;                       % Print stats?
nItersDB = 1;                       % # (D,B) updates per outer iteration
nIters1 = 10;                       % # outer iters for first batch
nIters  = 10;                       % # outer iters after first batch

% Load video
Xtrue = loadVideo(path,idim,T,dt);
[ny, nx, nt] = size(Xtrue);

% Missing data
M = (rand(ny,nx,nt) > p);
Y = Xtrue;
Y(~M) = 0;

% Interpolated video
Xi = interpVideo(Y,M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OnAIR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization
Xhat = zeros(ny,nx,nt);
tt = 1:dt:(nt + 1 - T);
ni = numel(tt);
cost = nan(nIters1,ni);

% Data
M1 = M(:,:,1:T);
Xi1 = Xi(:,:,1:T);
Xtrue1 = Xtrue(:,:,1:T);
Y1 = Y(:,:,1:T);

% OnAIR (t = 1)
opts = struct();
opts.A        = 1;
opts.M        = M1;
opts.xdim     = [ny, nx, T];
opts.pdim     = pdim;
opts.pgap     = pgap;
opts.type     = type;
opts.dr       = dr;
opts.ddim     = ddim;
opts.unitaryD = unitaryD;
opts.fixedD   = fixedD;
opts.nIters   = nIters1;
opts.nItersDB = nItersDB;
opts.nItersX  = -1; % Exact updates
opts.X0       = Xi1;
opts.D0       = D0;
opts.Xtrue    = Xtrue1;
opts.accel    = false;
opts.tau      = 1;
opts.flag     = flag;
[Xt, Dt, Bt, params, stats] = OnAIR(Y1,lambda,mu,gamma,opts);
cost(:,1) = stats.cost(:);

% Initialize reconstruction
Xhat(:,:,1:T) = Xt;

% OnAIR (t = 2,3,...)
for i = 2:ni
    % Data
    t = tt(i);
    Mt = M(:,:,t:(t + T - 1));
    Xit = Xi(:,:,t:(t + T - 1));
    Xtruet = Xtrue(:,:,t:(t + T - 1));
    Yt = Y(:,:,t:(t + T - 1));
    
    % Online DLS
    opts.M      = Mt;
    opts.nIters = nIters;
    %opts.X0    = zeroFill(Xhat,Yt,Mt,t,dt);
    opts.X0     = interpFill(Xhat,Xit,t,dt);
    opts.D0     = Dt;
    opts.B0     = Bt;
    opts.params = params;
    opts.Xtrue  = Xtruet;
    [Xt, Dt, Bt, params, stats] = OnAIR(Yt,lambda,mu,gamma,opts);
    cost(1:nIters,i) = stats.cost(:);
    
    % Update reconstruction
    Xhat = updateRecon(Xhat,Xt,gamma,t,dt);  % gamma-weighted average
    %Xhat = updateRecon(Xhat,Xt,1,t,dt);      % unweighted average
    %Xhat = updateRecon(Xhat,Xt,0,t,dt);      % use latest estimate
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualize results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
movie.video = cell2mat({Xtrue,Xi,Xhat});
opts.xlabels = {'Truth','Interpolated','OnAIR'};
PlayMovie(movie,opts);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cm = linspecer(2);
cfigure();

% OnAIR cost
plot(1:nIters1,cost(:,1),'-o','Color',cm(1,:)); hold on;
for i = 2:ni
    xx = tt(i) + (0:(nIters - 1));
    plot(xx,cost(1:nIters,i),'-o','Color',cm(2,:));
end
xlabel('time + iteration');
ylabel('cost');
title('OnAIR');
axis tight; padAxis();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

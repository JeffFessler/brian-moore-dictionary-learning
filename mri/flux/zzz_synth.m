%% Low-rank plus row-sparse model

% Knobs
m = 100;
n = 500;
gamma = 0.05;
sigma = 1;
theta1 = 10;
theta2 = 5;

% Ltrue
u1 = randn(n,1); u1 = u1 / norm(u1);
v1 = randn(m,1); v1 = v1 / norm(v1);
Ltrue = theta1 * (u1 * v1');

% Strue
Strue = (theta2 / sqrt(gamma * n)) * rand(n,m) .* (2 * (rand(n,m) < 0.5) - 1);
Strue(rand(n,1) > gamma,:) = 0;
%{
u2 = rand(n,1); u2(rand(n,1) > gamma) = 0; u2 = u2 / norm(u2);
v2 = randn(m,1); v2 = v2 / norm(v2);
Strue = theta2 * (u2 * v2');
%}

% Noise
Xn = (sigma / sqrt(n)) * randn(n,m);

% Observations
Xtrue = Ltrue + Strue;
Y = Xtrue + Xn;

% Perform optimization
model.r = 1;
model.lambdaL = sigma * sqrt(2);
model.lambdaS = 1.5 * sigma * sqrt(2 / (gamma * n));
opt.tau = 0.75;
opt.tol = -1; % Always run all iterations
opt.Nmax = 100;
opt.methodL = 'OptShrink'; % {'SVT','OptShrink'}
opt.methodS = 'mixed21';
truth.Xtrue = Xtrue;
truth.Ltrue = Ltrue;
truth.Strue = Strue;
[Lhat Shat NRMSE] = LpS_MRI(Y,model,opt,truth);
Xhat = Lhat + Shat;

% Compute errors
Xerr = Xtrue - Xhat;
Lerr = Ltrue - Lhat;
Serr = Strue - Shat;

% Plot NRMSEs
figure;
plotfcn = @plot; % {'semilogy','plot'}
subplot(1,3,1); plotfcn(NRMSE.X,'b-o'); xlabel('Iteration'); title('NRMSE(Xhat)');
subplot(1,3,2); plotfcn(NRMSE.L,'r-o'); xlabel('Iteration'); title('NRMSE(Lhat)');
subplot(1,3,3); plotfcn(NRMSE.S,'g-o'); xlabel('Iteration'); title('NRMSE(Shat)');

% Compute colorbar views
xx = [Xtrue(:)' Xhat(:)' Xerr(:)'];
ll = [Ltrue(:)' Lhat(:)' Lerr(:)'];
ss = [Strue(:)' Shat(:)' Serr(:)'];
viewX = [min(xx) max(xx)];
viewL = [min(ll) max(ll)];
viewS = [min(ss) max(ss)];

% Plot results
figure;
subplot(3,3,1); pcolor(Xtrue); shading flat; colorbar; caxis(viewX); title('Xtrue');
subplot(3,3,2); pcolor(Ltrue); shading flat; colorbar; caxis(viewL); title('Ltrue');
subplot(3,3,3); pcolor(Strue); shading flat; colorbar; caxis(viewS); title('Strue');
subplot(3,3,4); pcolor(Xhat); shading flat; colorbar; caxis(viewX); title('Xhat');
subplot(3,3,5); pcolor(Lhat); shading flat; colorbar; caxis(viewL); title('Lhat');
subplot(3,3,6); pcolor(Shat); shading flat; colorbar; caxis(viewS); title('Shat');
subplot(3,3,7); pcolor(Xerr); shading flat; colorbar; caxis(viewX); title('Xtrue - Xhat');
subplot(3,3,8); pcolor(Lerr); shading flat; colorbar; caxis(viewL); title('Ltrue - Lhat');
subplot(3,3,9); pcolor(Serr); shading flat; colorbar; caxis(viewS); title('Strue - Shat');

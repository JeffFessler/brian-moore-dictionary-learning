%% KT-SLR demo on invivo perfusion data

tic;

% Knobs
%nLines = 20;  % # k-space radial lines to sample per frame
nLines = 15;  % # k-space radial lines to sample per frame

% Load data
load('invivo_perfusion.mat');
[n1, n2, n3] = size(x);

% Generate sampling mask
mask = strucrand(n1,n2,n3,nLines);
mask = fftshift(fftshift(mask,1),2); 

%% Define the forward and backward operators A and Atranspose (At)

S  = find(mask ~= 0);
A  = @(z)  A_fhp3D(z,S);
At = @(z) At_fhp3D(z,S,n1,n2,n3);

%% Define the forward and backward gradient operators

% This function has an option of changing the step sizes of the gradients along the x, y and t dimensions
% For ex: step_size = [1,1,0.303] implies Grad_x(U) = [U(x+1,y,t)-U(x,y,t)]/1,
% Grad_y(U) = [U(x,y+1,t)-U(x,y,t)]/1, Grad_t(U) = [U(x,y,t+1)-U(x,y,t)]/0.303
[D, Dt] = defDDt([1, 1, 0.303]);

%% First guess, direct IFFT

b      = A(x); 
x_init = At(b);

%% Call kt slr using augmented Lagrangian with continuation (refer: S.G.Lingala et al, ISBI 2011)

opts.mu1        = 1e-12; % Regularization parameter for schatten p-norm
opts.mu2        = 4e-10; % Reg. parameter for spatiotemporal TV norm * Note: temporal wt weighted 10 times higher than spatial weight (see DefDDtmod.m)
opts.p          = 0.1;   % The value of p in Schatten p-norm; p=0.1: non convex; p = 1: convex
opts.beta1      = 1e6;   % The continuation parameter for low rank norm; initialize it. 
opts.beta2      = 1e6;   % The continuation parameter for the TV norm; Initialize it. 
opts.beta1rate  = 50;    % The continuation parametr increment for low rank norm
opts.beta2rate  = 25;    % Similar increment for TV norm
opts.outer_iter = 7;     % # outer iterations
opts.inner_iter = 50;    % # inner iterations
[recon, cost, opts] = minSNandTV(A,At,D,Dt,x_init,b,1,opts);

%% Plot results

figure();
colormap(gray);

subplot(3,3,1); imagesc(abs(x(:,:,24)));                title('Fully sampled gold standard, a spatial frame'); 
subplot(3,3,2); imagesc(abs(x_init(:,:,24)));           title('Direct IFFT, a spatial frame'); 
subplot(3,3,3); imagesc(abs(recon(:,:,24)));            title('k-t SLR, a spatial frame'); 
subplot(3,3,4); imagesc(abs(squeeze(x(100,:,:))));      title('Gold standard, image time profile'); 
subplot(3,3,5); imagesc(abs(squeeze(x_init(100,:,:)))); title('Direct IFFT, image time profile'); 
subplot(3,3,6); imagesc(abs(squeeze(recon(100,:,:))));  title('k-t SLR, image time profile'); 
subplot(3,3,8); imagesc(abs(fftshift(mask(:,:,5))));    title('Radial sampling: one frame'); 
subplot(3,3,9); plot(abs(cost),'linewidth',2);          title('Cost v/s iteration');  

tic
clear; clc;close all;

%%  
load aperiodic_pincat.mat; % load the PINCAT phantom
x=new; 

[n1,n2,n3]=size(x);
line = 24; % Number of radial rays per frame

[T3D] = strucrand(n1,n2,n3,line); % Generate the (kx,ky)-t sampling pattern; Each frame has uniformly spaced radial rays, with random rotations across frames

%%

mask = fftshift(fftshift(T3D,1),2);
S=find(mask~=0);


%%  Define the forward and backward operators A and Atranspose (At)
A = @(z)A_fhp3D(z,S); % The forward Fourier sampling operator
At=@(z)At_fhp3D(z,S,n1,n2,n3); % The backward Fourier sampling operator

 
%% Define the forward and backward gradient operators.
% This function has an option of changing the step sizes of the gradients along the x, y and t dimensions
% For ex: step_size = [1,1,0.303] implies Grad_x(U) = [U(x+1,y,t)-U(x,y,t)]/1,
% Grad_y(U) = [U(x,y+1,t)-U(x,y,t)]/1, Grad_t(U) = [U(x,y,t+1)-U(x,y,t)]/0.303
step_size = [1,1,1];
[D,Dt] = defDDt(step_size);

%% First guess, direct IFFT
b = A(x);
x_init = At(b);

%% Call k-t SLR using augmented Lagrangian with continuation (refer: S.G.Lingala et al, ISBI 2011)

mu1 =3e2; % Regularization parameter for the schatten p-norm
mu2 =1e-3; % Regularization parameter for the spatio-temporal TV norm
opts.mu1=mu1; 
opts.mu2 = mu2;
opts.p=0.1; % The value of p in Schatten p-norm; p = 0.1- non-convex, p = 1-convex nuclear norm. 
opts.beta1=1e-9; % The continuation parameter for the Schatten p-norm, initialize it
opts.beta2=1e-2; % The continuation parameter for the TV norm, initialize it
opts.beta1rate = 25; % Continuation parmeter increment for the Schatten p-norm
opts.beta2rate = 25; % Continuation parameter increment for the TV norm
opts.outer_iter =9; % Number of outer iterations
opts.inner_iter = 50; % Number of inner iterations


%% Call k-t SLR
[Recon,cost,opts] = minSNandTV(A,At,D,Dt,x_init,b,T3D,opts);

 %%  Display the results
 
close all; 

figure(1); colormap(gray);
subplot(3,3,1); imagesc(abs(x(:,:,24))); title('Gold standard, a spatial frame'); 
subplot(3,3,2); imagesc(abs(x_init(:,:,24))); title('Direct IFFT, a spatial frame'); 
subplot(3,3,3); imagesc(abs(Recon(:,:,24))); title('k-t SLR, a spatial frame'); 
subplot(3,3,4); imagesc(abs(squeeze(x(60,:,:)))); title('Gold standard, image time profile'); 
subplot(3,3,5); imagesc(abs(squeeze(x_init(60,:,:)))); title('Direct IFFT, image time profile'); 
subplot(3,3,6); imagesc(abs(squeeze(Recon(60,:,:)))); title('k-t SLR, image time profile'); 
subplot(3,3,8); imagesc(abs(T3D(:,:,5))); title('Radial sampling: one frame'); 
subplot(3,3,9); plot(abs(cost),'linewidth',2); title('Cost v/s iteration');  

toc

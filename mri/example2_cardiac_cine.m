%
% L+S reconstruction of undersampled multicoil cardiac cine MRI
%
% Ricardo Otazo (2013)
%

% Load undersampled data 
load('cardiac_cine_R6.mat');
[nx ny nt nc] = size(kdata);
nd = nx * ny;

%--------------------------------------------------------------------------
% L + S reconstruction
%--------------------------------------------------------------------------
E = Emat_xyt(kdata(:,:,:,1) ~= 0,b1);
T = TempFFT(3);
d = kdata;
C = max(svd(reshape(E' * d,[nd nt]))); % scaling constant for Otazo's SVT

% Correlated motion in the background
lambdaL = 0.0025;
lambdaS = 0.00125;

% Stationary background
%lambdaL = 0.01;
%lambdaS = 0.0025;
%r = 1;

tol = 0.0025;
MaxIters = 50;

[L S cost deltaM time its] = lps_ist(E,T,d,C * lambdaL,lambdaS,tol,MaxIters);
%[L S cost deltaM time its] = SLRN_MRI(E,T,d,r,lambdaS,tol,MaxIters);
%[L S cost deltaM time its] = RPCA_DN_MRI(E,T,d,lambda,epsilon,tol,MaxIters,NInnerIters);

M = L + S;
%--------------------------------------------------------------------------

% Extract 4 frames
Md = M(65:192,65:192,2);
Md = cat(2,Md,M(65:192,65:192,8));
Md = cat(2,Md,M(65:192,65:192,14));
Md = cat(2,Md,M(65:192,65:192,20));
Ld = L(65:192,65:192,2);
Ld = cat(2,Ld,L(65:192,65:192,8));
Ld = cat(2,Ld,L(65:192,65:192,14));
Ld = cat(2,Ld,L(65:192,65:192,20));
Sd = S(65:192,65:192,2);
Sd = cat(2,Sd,S(65:192,65:192,8));
Sd = cat(2,Sd,S(65:192,65:192,14));
Sd = cat(2,Sd,S(65:192,65:192,20));

%--------------------------------------------------------------------------
% Plot images
%--------------------------------------------------------------------------
figure;

subplot(3,1,1);
imshow(abs(Md),[0 1]);
ylabel('L + S');

subplot(3,1,2);
imshow(abs(Ld),[0 1]);
ylabel('L');

subplot(3,1,3);
imshow(abs(Sd),[0 1]);
ylabel('S');
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Plot algorithm convergence
%--------------------------------------------------------------------------
figure;

subplot(1,3,1);
plot(1:its,cost,'b-o');
title('Cost Function');
xlabel('Iteration');
set(gca,'XLim',[1 its]);

subplot(1,3,2);
semilogy(1:its,deltaM,'b-o');
title('Convergence Criterion');
xlabel('Iteration');
set(gca,'XLim',[1 its]);

subplot(1,3,3);
semilogy(1:its,cumsum(time),'b-o');
title('Time (Seconds)');
xlabel('Iteration');
set(gca,'XLim',[1 its]);
%--------------------------------------------------------------------------

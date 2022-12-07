%
% L + S reconstruction of dynamic contrast-enhanced abdominal MRI acquired
% with golden-angle radial sampling
% 
% Temporal resolution is flexible and determined by the user in the 
% variable nspokes 
%
% Ricardo Otazo (2013)
%

% Number of spokes to be used per frame (Fibonacci number)
Nspokes = 21;

% Load radial data
load('abdomen_dce_ga.mat');
[nx ny nc] = size(b1); %#ok
nd = nx * ny;
[nr ntviews nc] = size(kdata);

% Number of frames
nt = floor(ntviews / Nspokes);

% Crop the data according to the number of spokes per frame
kdata = kdata(:,1:(nt * Nspokes),:);
k = k(:,1:(nt * Nspokes));
w = w(:,1:(nt * Nspokes));

% Sort the data into a time-series of undersampled images
for ii = 1:nt
    kdatau(:,:,:,ii) = kdata(:,((ii-1) * Nspokes + 1):(ii * Nspokes),:); %#ok
    ku(:,:,ii) = k(:,((ii-1) * Nspokes + 1):(ii * Nspokes)); %#ok
    wu(:,:,ii) = w(:,((ii-1) * Nspokes + 1):(ii * Nspokes)); %#ok
end

%--------------------------------------------------------------------------
% L + S reconstruction
%--------------------------------------------------------------------------
% Multicoil NUFFT operator
E = MCNUFFT(ku,wu,b1);
d = kdatau;
recon_nufft = E' * d;
C = max(svd(reshape(recon_nufft,[nd nt]))); % scaling constant for Otazo's SVT
clear kdata k ku wu w;

% Regularization parameters
lambdaL = 0.025;
lambdaS = 0.5 * max(abs(recon_nufft(:)));

MaxIters = 20;
tol = 0.0025;

[L S cost deltaM time its] = lps_tv(E,d,C * lambdaL,lambdaS,tol,MaxIters);

L = flipdim(L,1);
S = flipdim(S,1);
M = L + S;
%--------------------------------------------------------------------------

% Extract 4 frames
Md = M(65:336,:,1);
Md = cat(2,Md,M(65:336,:,9));
Md = cat(2,Md,M(65:336,:,16));
Md = cat(2,Md,M(65:336,:,25));
Ld = L(65:336,:,1);
Ld = cat(2,Ld,L(65:336,:,9));
Ld = cat(2,Ld,L(65:336,:,16));
Ld = cat(2,Ld,L(65:336,:,25));
Sd = S(65:336,:,1);
Sd = cat(2,Sd,S(65:336,:,9));
Sd = cat(2,Sd,S(65:336,:,16));
Sd = cat(2,Sd,S(65:336,:,25));

%--------------------------------------------------------------------------
% Display results
%--------------------------------------------------------------------------
figure;

subplot(3,1,1);
imshow(abs(Md),[0 5e-4]);
ylabel('L + S');

subplot(3,1,2);
imshow(abs(Ld),[0 5e-4]);
ylabel('L');

subplot(3,1,3);
imshow(abs(Sd),[0 5e-4]);
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

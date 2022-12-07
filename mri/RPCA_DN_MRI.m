function [L S cost deltaM time its NRMSE SNR] = RPCA_DN_MRI(E,T,d,lambda,epsilon,tol,MaxIters,NInnerIters,Mtrue)
% 
% [L S cost deltaM time its] = RPCA_DN_MRI(E,T,d,lambda,epsilon,tol,MaxIters,NInnerIters);
% [L S cost deltaM time its NRMSE SNR] = RPCA_DN_MRI(E,T,d,lambda,epsilon,tol,MaxIters,NInnerIters,Mtrue);
% 
% d: undersampled k-t data (nx,ny,nt,nc)
% E: data acquisition operator
% T: sparsifying transform
% lambda: regularization parameter
% epsilon: denoising parameter
% tol: stopping tolerance
% MaxIters: max number of allowed iterations
% NInnerIters: number of inner PCG iterations to apply during M-update
% [OPTIONAL] Mtrue: true target image
% 
% L + S reconstruction of undersampled dynamic MRI data using a denoising
% formulation of nuclear + ell-1 norm minimization
% 
% Brian Moore (2013)
% 

% Check to see if we should compute NRMSEs and SNRs
if (nargin == 9)
    HaveTrueImages = true;
    Mtrue = abs(Mtrue); % absolute value
    %Mtrue = Mtrue - min(Mtrue(:)); % min value = 0
    %Mtrue = Mtrue / max(Mtrue(:)); % max value = 1
    NRMSEfcn = @(M) 100 * norm(M(:) - Mtrue(:),2) / norm(Mtrue(:),2);
    SNRfcn = @(M) 20 * log10(norm(Mtrue(:),2) / norm(M(:) - Mtrue(:),2));
else
    HaveTrueImages = false;
end
%--------------------------------------------------------------------------
% Parse inputs
%--------------------------------------------------------------------------
[nx ny nt nc] = size(d);
nd = nx * ny;
%rho = nnz(d) / (4 * sum(abs(d(:))));
rho = lambda;

M = reshape(E' * d,[nd nt]);
%M = zeros(nd,nt);
S = zeros(nd,nt);
U1 = zeros(nd,nt);
U2 = zeros(nx,ny,nt,nc);
%--------------------------------------------------------------------------

% L + S iterations
its = 0;
cost = [];
deltaM = [];
time = [];
NRMSE = [];
SNR = [];
fprintf('\n*** RPCA DN MRI ***\n')
while(1)
    % Record iteration time
    tic;
    
    % Update iteration count
	its = its + 1;
    
    % Low-rank update (Singular value thresholding)
    [UL SL VL] = svd(M - S - U1,0);
    sL = soft(diag(SL),1 / rho);
    L = UL * diag(sL) * VL';
    
	% Sparse update (Soft thresholding)
	S = reshape(T' * soft(T * reshape(M - L - U1,[nx ny nt]),lambda / rho),[nd nt]);
    
    % Denoising update (Euclidean ball projection)
    V = EuclidBallProj(E * reshape(M,[nx ny nt]) - U2,d,epsilon);
    
    % M-update (BBGM)
    Mlast = M;
    M = BBGM(L + S + U1,V + U2,E,nx,ny,nt,NInnerIters);
    
    % ADMM dual variable updates
    U1 = U1 + (L + S - M);
    U2 = U2 + (V - E * reshape(M,[nx ny nt]));
    
    % Print algorithm status
	temp = T * reshape(S,[nx ny nt]);
    cost(its,1) = sum(sL) + lambda * norm(temp(:),1); %#ok  cost function
    deltaM(its,1) = norm(M(:) - Mlast(:),2) / norm(Mlast(:),2); %#ok relative change
    time(its,1) = toc; %#ok
    fprintf('Iter: %d, Cost: %f3, deltaM: %f3\n',its,cost(its,1),deltaM(its,1)); 
    
    % Compute NRMSE/SNRs
    if HaveTrueImages
        Mhat = abs(L + S);
        %Mhat = Mhat - min(Mhat(:)); % min value = 0
        %Mhat = Mhat / max(Mhat(:)); % max value = 1
        NRMSE(its,1) = NRMSEfcn(Mhat); %#ok
        SNR(its,1) = SNRfcn(Mhat); %#ok
    end
    
    % Stopping criteria 
    if ((its >= MaxIters) || (deltaM(its,1) < tol))
        break;
    end
end

% Return L and S as images over time
L = reshape(L,[nx ny nt]);
S = reshape(S,[nx ny nt]);

end

% Soft-thresholding
function sY = soft(Y,tau)

sY = max(abs(Y) - tau,0) .* sign(Y);
sY(isnan(sY)) = 0;

end

% Orthogonal projection onto Euclidean ball
function Vebj = EuclidBallProj(V,X,epsilon)

VmX = V - X;
nrm = norm(VmX(:),2);
if (nrm > epsilon)
    Vebj = X + (epsilon / nrm) * VmX;
else
    Vebj = V;
end

end

% BBGM for the M-update
function M = BBGM(B1,B2,E,nx,ny,nt,NInnerIters)

% Parse inputs
nd = nx * ny;

% Gradient
grad = @(M) M + reshape(E' * (E * reshape(M,[nx ny nt]) - B2),[nd nt]) - B1;

% BBGM
M = B1 + reshape(E' * B2,[nd nt]); % initial M
for iter = 1:NInnerIters
    % Gradient
    gM = grad(M);
    
    % Step size
    if (iter == 1)
        alpha = 1;
    else
        MmMlast = M(:) - Mlast(:);
        alpha = norm(MmMlast,2)^2 / dot(MmMlast,gM(:) - gMlast(:));
    end
    
    % Save last state
    if (iter < NInnerIters)
        Mlast = M;
        gMlast = gM;
    end
    
    % Update
    M = M - alpha * gM;
    
    % Display progress
    %fprintf('   Inner Iter: %d, norm(grad): %f3\n',iter,norm(gM(:),2));
end

end

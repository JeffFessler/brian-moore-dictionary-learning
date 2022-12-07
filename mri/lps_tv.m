function [L S cost deltaM time its] = lps_tv(E,d,lambdaL,lambdaS,tol,MaxIters)
% 
% [L S cost deltaM time its] = lps_ist(E,d,lambdaL,lambdaS,tol,MaxIters);
% 
% d: undersampled k-t data (nx,ny,nt,nc)
% E: data acquisition operator
% lambdaL: nuclear-norm weight
% lambdaS: temporal TV weight
% tol: stopping tolerance
% MaxIters: max number of allowed iterations
% 
% L + S reconstruction of undersampled dynamic MRI data using iterative
% soft-thresholding of singular values of L and iterative clipping of
% entries of S (temporal total variation on S)
% 
% Ricardo Otazo (2013)
% 

% Initializations
M = E' * d;
[nx ny nt] = size(M);
M = reshape(M,[nx * ny,nt]);
S = zeros(nx * ny,nt);
z = zeros(nx * ny,nt-1);
its = 0;

% L + S Algorithm iterations
cost = [];
deltaM = [];
time = [];
fprintf('\n*** L + S TV ***\n')
while(1)
    % Record iteration time
    tic;
    
    % Update iteration count
	its = its + 1;
    
    % Save last iterate
    Mlast = M;
    
    % Low-rank update
    [UL SL VL] = svd(M - S,0);
    sL = SoftThresh(diag(SL),lambdaL);
    %sL = SoftThresh(diag(SL),SL(1) * lambdaL); % This is not what your paper said you were doing
    L = UL * diag(sL) * VL';
    
    % Sparse update - TV using clipping
    z = z + 0.25 * diff(M - L,1,2);
    z = sign(z) .* reshape(max(min(abs(z(:)),lambdaS / 2),-lambdaS / 2),[nx * ny,nt - 1]);
    z(isnan(z)) = 0;
    adjDz(:,1) = -z(:,1);
    adjDz(:,2:nt-1) = -diff(z,1,2);
    adjDz(:,nt) = z(:,end); %#ok
    S = M - L - adjDz;
    
    % Data consistency
    resid = E * reshape(L + S,[nx ny nt]) - d;
    M = L + S - reshape(E' * resid,[nx * ny,nt]);
    
    % Print algorithm status
    cost(its,1) = norm(resid(:),2)^2 + lambdaL * sum(sL) + lambdaS * norm(diff(S,1,2),1); %#ok cost function
    deltaM(its,1) = norm(M(:) - Mlast(:),2) / norm(Mlast(:),2); %#ok relative change
    time(its,1) = toc; %#ok
    fprintf('Iter: %d, Cost: %f3, deltaM: %f3\n',its,cost(its,1),deltaM(its,1)); 
    
    % Stopping criteria 
    if ((its >= MaxIters) || (deltaM(its,1) < tol))
        break;
    end
end

% Return L and S as images over time
L = reshape(L,[nx ny nt]);
S = reshape(S,[nx ny nt]);

end

% Soft-thresholding function
function y = SoftThresh(x,p)

y = max(abs(x) - p,0) .* sign(x);
y(isnan(y)) = 0;

end

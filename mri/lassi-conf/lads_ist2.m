function [L,S,D,B,ite,time,cost,mse,delta,sparsity,amse] = lads_ist2(param)
% Syntax:   [L,S,D,B,ite,time,cost,mse,delta,sparsity,amse] = lads_ist2(param);

% Check if we're using OptShrink
USE_OPTSHRINK = isfield(param,'r');

% Initial dictionary
D = param.D;
[~, K] = size(D);

% Patch sizes
lx = param.patchsize(1);
ly = param.patchsize(2);
lt = param.patchsize(3);

% Patch overlap strides
slidingDis   = param.slidedis(1);
tSlideingDis = param.slidedis(2);

% Atomic rank
p = param.p;

% Weights
l2 = param.lambda_B;
ls = param.lambda_S;

% Initialize L and S
L = param.L0;
S = param.S0;
[nx, ny, nt] = size(L);
L = reshape(L,[nx * ny,nt]);
S = reshape(S,[nx * ny,nt]);

% Initialize M
resk = param.E * reshape(L + S,[nx,ny,nt]) - param.d;    
M    = L + S - reshape(param.E' * resk,[nx * ny,nt]);

% Initialize blocks
[blocks, ~, rows, cols, temps] = my_im2col_3D(reshape(S,nx,ny,nt),[lx,ly,lt],slidingDis,tSlideingDis);

% Initialize B
%B = sparse(zeros(K,size(blocks,2)));
B  = zeros(K,size(blocks,2));

% Iterations
if USE_OPTSHRINK
    disp('***** L + ADS (OptShrink) *****')
else
    disp('***** L + ADS (SVT) *****')
end
ite      = 0;
time     = nan(param.nite,1);
cost     = nan(param.nite,1);
mse      = nan(param.nite,1);
delta    = nan(param.nite,1);
sparsity = nan(param.nite,1);
amse     = nan(param.nite,1);
while true
    itimer = tic;
	ite = ite + 1;
    
    % (D,B) updates
    [D,B,~] = DLSVDKronL0v(blocks,D,B,l2(ite),lx * ly,lt,p,param.niteDL,0);
    ZZ      = D * B;
    
    % (L,S,M) updates
	M0 = M;
    for pt = 1:param.niteLS
        % Last iterates
        Ltild = M - S;
        Stild = M - L;
        
        % L update
        if USE_OPTSHRINK
            % OptShrink
            L = OptShrink(Ltild,param.r);
        else
            % SVT
            [Ut, St, Vt] = svd(Ltild,0);
            St = diag(SoftThresh(diag(St),param.lambda_L));
            L  = Ut * St * Vt';
        end
        
        % S update
        if pt == 1
            IMO = zeros(nx,ny,nt); 
            if ite == 1
                WO = zeros(nx,ny,nt); % The diagonal inversion matrix in the S update
            end
            Nr = 30000;
            % Summing patch approximations at 3D locations
            for jj = 1:Nr:size(blocks,2)
                jumpSize = min(jj + Nr - 1,size(blocks,2));
                ZZ2 = ZZ(:,jj:jumpSize);
                for ii = jj:jumpSize
                    col   = cols(ii);
                    row   = rows(ii);
                    temp  = temps(ii);
                    block = reshape(ZZ2(:,ii - jj + 1),[lx,ly,lt]);    
                    IMO(row:(row + lx - 1),col:(col + ly - 1),temp:(temp + lt - 1)) = IMO(row:(row + lx - 1),col:(col + ly - 1),temp:(temp + lt - 1)) + block;
                    if ite == 1
                     WO(row:(row + lx - 1),col:(col + ly - 1),temp:(temp + lt - 1)) = WO(row:(row + lx - 1),col:(col + ly - 1),temp:(temp + lt - 1)) + ones(lx,ly,lt);
                    end
                end
            end
            if ite == 1
                WO = ones(nx,ny,nt) + (2 * ls) * WO;
            end
        end
        S = (Stild + (2 * ls) * reshape(IMO,nx * ny,nt)) ./ reshape(WO,nx * ny,nt);
        
        % M update
        resk = param.E * reshape(L + S,[nx,ny,nt]) - param.d;    
        M    = L + S - reshape(param.E' * resk,[nx * ny,nt]);
    end
	
    % Extract patches
    [blocks, ~, rows, cols, temps] = my_im2col_3D(reshape(S,nx,ny,nt),[lx,ly,lt],slidingDis,tSlideingDis);
    
    % Save stats
    time(ite)     = toc(itimer);
    if ~USE_OPTSHRINK
        cost(ite) = 0.5 * norm(resk(:),2)^2 + param.lambda_L * sum(diag(St)) + param.lambda_S * (norm(blocks(:) - ZZ(:))^2 + l2(end)^2 * nnz(B));
    end
    mse(ite)      = norm(param.Lfull(:) - (L(:) + S(:)));
    delta(ite)    = norm(M(:) - M0(:)) / norm(M0(:));
    sparsity(ite) = 100 * (nnz(abs(B)) / numel(B));
    amse(ite)     = norm(abs(param.Lfull(:)) - abs(L(:) + S(:)));
    
    % Print stats
	fprintf(' ite: %03d, time: %.2f, cost: %.3f, mse: %.3f, amse: %.3f, delta: %.2e, sparsity: %.2f%%\n',ite,time(ite),cost(ite),mse(ite),amse(ite),delta(ite),sparsity(ite));
    
    % Check stopping criteria
    if (ite >= param.nite) || (delta(ite) < param.tol), break; end % consider dropping param.tol and using only param.nite
end
L = reshape(L,nx,ny,nt);
S = reshape(S,nx,ny,nt);

time     = time(1:ite);
cost     = cost(1:ite);
mse      = mse(1:ite);
delta    = delta(1:ite);
sparsity = sparsity(1:ite);
amse     = amse(1:ite);

% Soft thresholding
function y = SoftThresh(x,lambda)
y = max(abs(x) - lambda,0) .* sign(x);
y(isnan(y)) = 0;

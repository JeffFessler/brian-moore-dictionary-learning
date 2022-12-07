function [D,X,obj]= DLSVDKronL0v(Y,D,X,l2,na,nb,p,numiter,ct)
%Y: training data matrix
%D: initial dictionary (with kronecker structure)
%X: initial sparse code matrix
%l2: sparsity penalty weight
%na: number of rows in reshaped dictionary atom (e.g., 64)
%nb: number of columns in reshaped dictionary atom (e.g., 8 or 12)
%p: maximum allowed rank per atom (p < min(na,nb))
%numiter: Number of iterations of block coordinate descent minimization
%ct: if set to 1, the objective value is computed in each iteration

[n, K] = size(D);
Z = [1; zeros(n - 1,1)]; %If there is a zero row in X, the corresponding dictionary atom update solution is non-unique and is set to Z.

sd = 1:K;  %order of sequencing over dictionary atoms in block coordinate descent
%sd = randperm(K);

obj = zeros(1,numiter);

for it = 1:numiter
    % Dictionary update and sparse coding
    itimer = tic;
    
    %sd = randperm(K);
    ZP = D' * Y;
    [colx, rowx, vx] = find(X');
    
    for j = sd  
        ind2 = colx(rowx == j);
        gamm = vx(rowx == j)';
        Dj2  = D(:,j);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % brimoor: leave X as full matrix
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %sparse coding 
        %h = -(Dj2' * D) * sparse(X);
        h =  -(Dj2' * D) * X;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        h(ind2) = h(ind2) + gamm;
        h = ZP(j,:) + h;   
        h = h .* (abs(h) >= l2); %hard-thresholding to find jth row of X (l0 norm)
        %ll = l2 / 2; h = ((h >= ll) .* (h - ll)) + ((h <= (-ll)) .* (h + ll)); % soft thresholding to find jth row of X (uncomment for l1 norm)
        ind = find(h);
        
        %atom update
        if ~any(h)
            Dj2 = Z;
        else
            [ind3, inna, ~] = intersect(ind2,ind);
            Dj2 = Y(:,ind) * h(ind)' - D * (X(:,ind) * h(ind)') + Dj2 * (gamm(:,inna) * h(ind3)');
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % brimoor: Skip SVD if possible. Also, svd() is faster than
            %          svds() for small matrices
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if p < min(na,nb)
                %[b,S,a] = svds(reshape(Dj2,na,nb),p);
                [bb, SS, aa] = svd(reshape(Dj2,na,nb),'econ');
                b = bb(:,1:p);
                S = SS(1:p,1:p);
                a = aa(:,1:p);
                
                %Dj2 = kron(a,b); %update jth column of D
                Dj2  = reshape(b * S * a',n,1);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            Dj2 = Dj2 / norm(Dj2,2);
        end
        
        % Record new values
        X(j,ind2) = 0;
        X(j,ind)  = h(ind);
        D(:,j)    = Dj2;
    end
    
    if ct == 1
        obj(it) = norm(Y - D * X,'fro')^2 + l2^2 * nnz(X);
    elseif ct == 2
        fprintf('  DLSVDKronL0v[%.2fs]\n',toc(itimer));
    end
end

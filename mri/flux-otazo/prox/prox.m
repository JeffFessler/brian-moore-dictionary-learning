function Yhat = prox(Y,str,params)
% Syntax:   Yhat = prox(Y,str,params);
%
% Apply proximal operator
%

% Apply proximal operator
switch lower(str)
    case 'soft'
        % Apply soft-thresholding
        Yhat = soft(Y,params{:});
    case 'softp'
        % Apply ell-p shrinkage
        Yhat = pShrink(Y,params{:});
    case 'mixed21'
        % Apply mixed 2-1 norm shrinkage
        Yhat = mixed21(Y,params{:});
    case 'svt'
        % Apply singular value thresholding
        Yhat = SVT(Y,params{:});
    case 'svt_llr'
        % Apply locally low-rank (LLR) singular value thresholding
        mask = params{1};
        params = {params{3:end} params{2}};
        if isempty(mask)
            % No mask
            Yhat = SVT_LLR(Y,params{:});
        else
            % Embed/remask data
            Yhat = masker(SVT_LLR(embed(Y,mask),params{:}),mask);
        end
    case 'svtp'
        % Apply ell-p singular value shrinkage
        Yhat = SVTp(Y,params{:});
    case 'svtp_llr'
        % Apply locally low-rank (LLR) ell-p singular value shrinkage
        mask = params{1};
        params = {params{3:end} params{2}};
        if isempty(mask)
            % No mask
            Yhat = SVTp_LLR(Y,params{:});
        else
            % Embed/remask data
            Yhat = masker(SVTp_LLR(embed(Y,mask),params{:}),mask);
        end
    case 'optshrink'
        % Apply OptShrink
        Yhat = OptShrink(Y,params{:});
    case 'optshrink_llr'
        % Apply locally low-rank (LLR) OptShrink
        mask = params{1};
        params = {params{3:end} params{2}};
        if isempty(mask)
            % No mask
            Yhat = OptShrink_LLR(Y,params{:});
        else
            % Embed/remask data
            Yhat = masker(OptShrink_LLR(embed(Y,mask),params{:}),mask);
        end
    case 'pca'
        % Apply PCA
        Yhat = PCA(Y,params{:});
    case 'pca_llr'
        % Apply locally low-rank (LLR) PCA
        mask = params{1};
        params = {params{3:end} params{2}};
        if isempty(mask)
            % No mask
            Yhat = PCA_LLR(Y,params{:});
        else
            % Embed/remask data
            Yhat = masker(PCA_LLR(embed(Y,mask),params{:}),mask);
        end
    case 'fro'
        % Apply frobenius shrinkage
        Yhat = FroShrink(Y,params{:});
    case 'ball2'
        % Apply Euclidean ball projection
        Yhat = EuclideanBallProject(Y,params{:});
end

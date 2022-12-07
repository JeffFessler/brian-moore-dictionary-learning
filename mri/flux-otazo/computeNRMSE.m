function nrmse = computeNRMSE(Xhat,Xtrue,metric,ROIs)
% Syntax:   nrmse = computeNRMSE(Xhat,Xtrue);
%           nrmse = computeNRMSE(Xhat,Xtrue,metric,ROIs);

% Parse inputs
if ~exist('metric','var') || isempty(metric)
    % Default = Whole image
    metric = 0;
end

% Compute NRMSE
nf = size(Xtrue,3);
NRMSEfcn = @(Xhat,Xtrue) norm(Xhat(:) - Xtrue(:)) / norm(Xtrue(:));
switch metric
    case 1
        % Heart-only
        M1 = repmat(ROIs{1},[1 1 nf]);
        nrmse = NRMSEfcn(Xhat(M1),Xtrue(M1));
    case 2
        % Body-only
        M2 = repmat(ROIs{2},[1 1 nf]);
        nrmse = NRMSEfcn(Xhat(M2),Xtrue(M2));
    otherwise
        % Whole image
        nrmse = NRMSEfcn(Xhat,Xtrue);
end

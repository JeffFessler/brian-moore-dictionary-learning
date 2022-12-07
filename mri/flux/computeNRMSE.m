function nrmse = computeNRMSE(Xhat,Xtrue,metric,ROIs)
% Syntax:   nrmse = computeNRMSE(Xhat,Xtrue);
%           nrmse = computeNRMSE(Xhat,Xtrue,metric,ROIs);

% Parse inputs
if ~exist('metric','var') || isempty(metric)
    % Default = All data
    metric = 1;
end

% Compute NRMSE
nf = size(Xtrue,3);
NRMSEfcn = @(Xhat,Xtrue) norm(Xhat(:) - Xtrue(:)) / norm(Xtrue(:));
switch metric
    case 1
        % X
        nrmse = NRMSEfcn(Xhat,Xtrue);
    case 2
        % Lesion #1
        M1 = repmat(ROIs{1},[1 1 nf]);
        nrmse = NRMSEfcn(Xhat(M1),Xtrue(M1));
    case 3
        % Lesion #2
        M2 = repmat(ROIs{2},[1 1 nf]);
        nrmse = NRMSEfcn(Xhat(M2),Xtrue(M2));
    case 4
        % Lesion #3
        M3 = repmat(ROIs{3},[1 1 nf]);
        nrmse = NRMSEfcn(Xhat(M3),Xtrue(M3));
    case 5
        % All lesions
        M1 = repmat(ROIs{1},[1 1 nf]);
        M2 = repmat(ROIs{2},[1 1 nf]);
        M3 = repmat(ROIs{3},[1 1 nf]);
        M = (M1 | M2 | M3);
        nrmse = NRMSEfcn(Xhat(M),Xtrue(M));
    otherwise
        % Error
        error('Unsupported metric');
end

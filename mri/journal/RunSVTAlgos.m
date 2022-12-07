function RunSVTAlgos(nc,nt,SNR,normA,T,basepath,algoStr,minLambdaL,maxLambdaL,NlambdaL,minLambdaS,maxLambdaS,NlambdaS,nb)
% Syntax:   RunSVTAlgos(nc,nt,SNR,normA,T,basepath,algoStr,minLambdaL,maxLambdaL,NlambdaL,minLambdaS,maxLambdaS,NlambdaS,nb);

% Load Fessler phantom
[Y A Xtrue Xtrue480 Xfft mask ROIs nd splineInterp] = GenerateFesslerPhantom(nc,nt,SNR);
Xtrue = masker(Xtrue,mask);
X0 = masker(Xfft,mask);

% Save nt-frame truth
truth = struct('Xtrue',Xtrue);

% NRMSE function handle
NRMSEfcn = @(Xhat,Xtrue) norm(Xhat(:) - Xtrue(:)) / norm(Xtrue(:));

% Constants
CL = svd(Xtrue); CL = CL(2); % lambdaL rel. to 2nd sing. val
CS = max(abs(Xtrue(:)));

% Initializations
lambdaL = CL * logspace(log10(minLambdaL),log10(maxLambdaL),NlambdaL);
lambdaS = CS * logspace(log10(minLambdaS),log10(maxLambdaS),NlambdaS);
recon = repmat(struct('X',[],'L',[],'S',[],'NRMSE',[]),[NlambdaL NlambdaS]);
NRMSE = repmat(struct('X',[],'ROI',[],'ROI1',[],'ROI2',[],'ROI3',[]),[NlambdaL NlambdaS]);
if exist('nb','var')
    [nx ny] = size(mask);
    Nb = floor(nx / nb(1)) * floor(ny / nb(2));
    lambdaL = lambdaL / Nb;
end

% Construct path
if exist('nb','var')
    % Append locally low-rank (LLR) block size
    basepath = sprintf('%s_nb_%i',basepath,nb);
end
basepath = sprintf('%s_nc_%i_SNR_%i_nt_%i',basepath,nc,SNR,nt);

% Perform reconstructions
for j = 1:NlambdaS
    % SVT formulations
    for i = 1:NlambdaL
        % Perform reconstruction
        rtimer = tic;
        method = struct();
        method.normA = normA;
        method.str = 'SVT';
        method.lambdaL = lambdaL(i);
        method.lambdaS = lambdaS(j);
        if exist('nb','var')
            method.nb = nb;
        end
        method.X0 = X0;
        [NRMSEij reconij] = RunLSAlgo(algoStr,Y,A,T,nd,1,nt,method,truth);
        reconij.NRMSE = NRMSEij.X;
        recon(i,j) = reconij;
        
        % Compute NRMSEs
        Xhat = splineInterp(abs(embed(reconij.X,mask)));
        NRMSE(i,j).X = NRMSEfcn(Xhat,Xtrue480);
        M1 = repmat(ROIs{1},[1 1 480]);
        NRMSE(i,j).ROI1 = NRMSEfcn(Xhat(M1),Xtrue480(M1));
        M2 = repmat(ROIs{2},[1 1 480]);
        NRMSE(i,j).ROI2 = NRMSEfcn(Xhat(M2),Xtrue480(M2));
        M3 = repmat(ROIs{3},[1 1 480]);
        NRMSE(i,j).ROI3 = NRMSEfcn(Xhat(M3),Xtrue480(M3));
        M = (M1 | M2 | M3);
        NRMSE(i,j).ROI = NRMSEfcn(Xhat(M),Xtrue480(M));
        
        % Display progress
        fprintf('*** SVT (%i,%i)/(%i,%i) complete - Time = %.2fs ***\n',j,i,NlambdaS,NlambdaL,toc(rtimer));
    end
end

% Save NRMSEs
NRMSEpath = sprintf('%s_NRMSE.mat',basepath);
if exist('nb','var')
    save(NRMSEpath,'NRMSE','lambdaL','lambdaS','algoStr','SNR','nt','nc','nb','-v7.3');
else
    save(NRMSEpath,'NRMSE','lambdaL','lambdaS','algoStr','SNR','nt','nc','-v7.3');
end

% Save reconstructions
stimer = tic;
reconpath = sprintf('%s_recon.mat',basepath);
if exist('nb','var')
    save(reconpath,'recon','lambdaL','lambdaS','algoStr','SNR','nt','nc','mask','splineInterp','X0','Xtrue','ROIs','nb','-v7.3');
else
    save(reconpath,'recon','lambdaL','lambdaS','algoStr','SNR','nc','nt','mask','splineInterp','X0','Xtrue','ROIs','-v7.3');
end
D = dir(reconpath);
fprintf('Data saved - Size = %.2fGB - Time = %.2fs\n',D.bytes / 1e9,toc(stimer));

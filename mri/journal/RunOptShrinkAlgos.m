function RunOptShrinkAlgos(nc,nt,SNR,normA,T,basepath,algoStr,Nr,minLambdaS,maxLambdaS,NlambdaS,nb)
% Syntax:   RunOptShrinkAlgos(nc,nt,SNR,normA,T,basepath,algoStr,Nr,minLambdaS,maxLambdaS,NlambdaS,nb);

% Load Fessler phantom
[Y A Xtrue Xtrue480 Xfft mask ROIs nd splineInterp] = GenerateFesslerPhantom(nc,nt,SNR);
Xtrue = masker(Xtrue,mask);
X0 = masker(Xfft,mask);

% Save nt-frame truth
truth = struct('Xtrue',Xtrue);

% NRMSE function handle
NRMSEfcn = @(Xhat,Xtrue) norm(Xhat(:) - Xtrue(:)) / norm(Xtrue(:));

% Constants
CS = max(abs(Xtrue(:)));

% Initializations
r = 1:Nr;
lambdaS = CS * logspace(log10(minLambdaS),log10(maxLambdaS),NlambdaS);
recon = repmat(struct('X',[],'L',[],'S',[],'NRMSE',[]),[Nr NlambdaS]);
NRMSE = repmat(struct('X',[],'ROI',[],'ROI1',[],'ROI2',[],'ROI3',[]),[Nr NlambdaS]);

% Construct path
if exist('nb','var')
    % Append locally low-rank (LLR) block size
    basepath = sprintf('%s_nb_%i',basepath,nb);
end
basepath = sprintf('%s_nc_%i_SNR_%i_nt_%i',basepath,nc,SNR,nt);

% Perform reconstructions
for j = 1:NlambdaS
    % OptShrink formulations
    for i = 1:Nr
        if (r(i) > min(nd,nt))
            continue;
        end
        
        % Perform reconstruction
        rtimer = tic;
        method = struct();
        method.normA = normA;
        method.str = 'OptShrink';
        method.r = r(i);
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
        fprintf('*** OptShrink (%i,%i)/(%i,%i) complete - Time = %.2fs ***\n',j,i,NlambdaS,Nr,toc(rtimer));
    end
end

% Save NRMSEs
NRMSEpath = sprintf('%s_NRMSE.mat',basepath);
if exist('nb','var')
    save(NRMSEpath,'NRMSE','r','lambdaS','algoStr','SNR','nt','nc','nb','-v7.3');
else
    save(NRMSEpath,'NRMSE','r','lambdaS','algoStr','SNR','nt','nc','-v7.3');
end

% Save reconstructions
stimer = tic;
reconpath = sprintf('%s_recon.mat',basepath);
if exist('nb','var')
    save(reconpath,'recon','r','lambdaS','algoStr','SNR','nt','nc','mask','splineInterp','X0','Xtrue','ROIs','nb','-v7.3');
else
    save(reconpath,'recon','r','lambdaS','algoStr','SNR','nt','nc','mask','splineInterp','X0','Xtrue','ROIs','-v7.3');
end
D = dir(reconpath);
fprintf('Data saved - Size = %.2fGB - Time = %.2fs\n',D.bytes / 1e9,toc(stimer));

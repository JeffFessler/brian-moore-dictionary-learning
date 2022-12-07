function RunFFTAlgos(nc,nt,SNR,basepath)
% Syntax:   RunFFTAlgos(nc,nt,SNR,basepath);

% Intialize data storage
Nt = length(nt);
recon = repmat(struct('X',[]),[1 Nt]);
NRMSE = repmat(struct('X',[],'ROI',[],'ROI1',[],'ROI2',[],'ROI3',[]),[1 Nt]);

% Construct path
basepath = sprintf('%s_nc_%i_SNR_%i',basepath,nc,SNR);

% NRMSE function handle
NRMSEfcn = @(Xhat,Xtrue) norm(Xhat(:) - Xtrue(:)) / norm(Xtrue(:));

% Loop over # frames
for i = 1:Nt
    % Load Fessler phantom + FFT-based solution
    rtimer = tic;
    [~,~,~,Xtrue480,Xfft,mask,ROIs,~,splineInterp] = GenerateFesslerPhantom(nc,nt(i),SNR);
    recon(i).X = masker(Xfft,mask);
    
    % Compute NRMSEs
    Xfft = splineInterp(abs(Xfft));
    NRMSE(i).X = NRMSEfcn(Xfft,Xtrue480);
    M1 = repmat(ROIs{1},[1 1 480]);
    NRMSE(i).ROI1 = NRMSEfcn(Xfft(M1),Xtrue480(M1));
    M2 = repmat(ROIs{2},[1 1 480]);
    NRMSE(i).ROI2 = NRMSEfcn(Xfft(M2),Xtrue480(M2));
    M3 = repmat(ROIs{3},[1 1 480]);
    NRMSE(i).ROI3 = NRMSEfcn(Xfft(M3),Xtrue480(M3));
    M = (M1 | M2 | M3);
    NRMSE(i).ROI = NRMSEfcn(Xfft(M),Xtrue480(M));
    
    % Display progress
    fprintf('*** Frame %i/%i complete - Time = %.2fs ***\n',i,Nt,toc(rtimer));
end

% Save NRMSEs
NRMSEpath = sprintf('%s_NRMSE.mat',basepath);
save(NRMSEpath,'NRMSE','SNR','nt','nc','-v7.3');

% Save reconstructions
stimer = tic;
reconpath = sprintf('%s_recon.mat',basepath);
save(reconpath,'recon','SNR','nt','nc','mask','splineInterp','ROIs','-v7.3');
D = dir(reconpath);
fprintf('Data saved - Size = %.2fGB - Time = %.2fs\n',D.bytes / 1e9,toc(stimer));

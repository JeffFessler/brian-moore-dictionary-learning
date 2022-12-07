function interp_agg(parFcn)
% Syntax: interp_agg(parFcn);

% Get parameters
[vars, ~] = parFcn();

% Initialize storage
nP = numel(vars.p);
nS = numel(vars.seed);
nt = vars.idim(3);
stats.time = nan(nP,nS,1);
err.NRMSE  = nan(nP,nS,1);
err.pSNR   = nan(nP,nS,1);
err.pNRMSE = nan(nP,nS,nt);
err.ppSNR  = nan(nP,nS,nt);

% Load data
files = findMatchingFiles(vars.rawpath);
nFiles = numel(files);
for i = 1:nFiles
    itimer = tic;
    
    % Extract index
    [~, namei, ~] = fileparts(files{i});
    idx = str2double(regexp(namei,'\d+','match'));
    [ii, jj] = ind2sub([nP, nS],idx);
    
    % Save data
    datai = load(files{i},'stats','err');
    stats.time(ii,jj,:) = datai.stats.time;
    err.NRMSE(ii,jj,:)  = datai.err.NRMSE;
    err.pSNR(ii,jj,:)   = datai.err.pSNR;
    err.pNRMSE(ii,jj,:) = datai.err.pNRMSE;
    err.ppSNR(ii,jj,:)  = datai.err.ppSNR;
    
    % Display progress
    fprintf('File %i/%i processed [Time = %.2fs]\n',i,nFiles,toc(itimer));
end

%{
% Load min-NRMSE reconstruction
[~, im] = nanmin(err.NRMSE(:));
datam   = load(files{im},'Xhat');
Xhat    = datam.Xhat; %#ok

% Save results
save(vars.outpath,'vars','stats','err','Xhat');
%}
save(vars.outpath,'vars','stats','err');
fprintf('DONE\n');

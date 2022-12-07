function irtsetup(irtdir)
% Syntax:   irtsetup(irtdir);

% Handle quick return
if isempty(irtdir), return; end

% Parse input path
irtdir = regexprep(irtdir,'\','/');
if irtdir(end) == '/'
    irtdir = irtdir(1:(end - 1));
end

% Add relevant irt directories to path
addpath([irtdir '/align']);             % image registration
addpath([irtdir '/align/mex']);         % image registration mex files
addpath([irtdir '/blob']);              % blob (KB) basis
addpath([irtdir '/ct']);                % x-ray CT (polyenergetic) recon
addpath([irtdir '/data']);              % example data
addpath([irtdir '/emission']);          % emission image reconstruction
addpath([irtdir '/example']);           % example applications
addpath([irtdir '/fbp']);               % FBP (filtered backprojection) code
addpath([irtdir '/general']);           % generic image reconstruction
addpath([irtdir '/graph']);             % graphics routines
addpath([irtdir '/mri']);               % MRI reconstruction
addpath([irtdir '/mri-rf/yip-spsp']);	% MRI RF pulse design
%addpath([irtdir '/mri/recon']);        % MRI reconstruction - old
addpath([irtdir '/nufft']);             % nonuniform FFT (for a fast projector)
addpath([irtdir '/nufft/table']);       % mex files for NUFFT
addpath([irtdir '/penalty']);           % regularizing penalty functions
addpath([irtdir '/systems']);           % system "matrices"
addpath([irtdir '/systems/tests']);     % tests of systems
addpath([irtdir '/transmission']);      % transmission image reconstruction
addpath([irtdir '/utilities']);         % various utility functions
addpath([irtdir '/wls']);               % weighted least-squares (WLS) estimates
addpath([irtdir '/mex/v7']);            % mex files

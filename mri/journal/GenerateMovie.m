function GenerateMovie(X,outpath,fps,res,lims)
%--------------------------------------------------------------------------
% Syntax:       GenerateMovie(X,outpath,fps);
%               GenerateMovie(X,outpath,fps,res);
%               GenerateMovie(X,outpath,fps,res,lims);
%               
% Inputs:       X is an (ny x nx x nt) cube containing nt video frames
%               
%               outpath is the desired output filename/path (without
%               extension)
%               
%               fps is the desired frames per second for the movie
%               
%               [OPTIONAL] lims = [minlim maxlim] contains the desired
%               limits for image display. The default value is
%               lims = [min(X(:)) max(X(:))]
%               
% Description:  This function generates an .avi video of the input data X
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         August 25, 2014
%--------------------------------------------------------------------------

% Video formats
%format = 'Motion JPEG AVI';
%format = 'Uncompressed AVI';
%format = 'Archival';
%format = 'Motion JPEG 2000';
format = 'MPEG-4';

% Parse inputs
if ~exist('res','var') || isempty(res)
    % Default resolution
    res = size(X,1);
end
[ny nx] = size(imresize(X(:,:,1),[res nan]));
nt = size(X,3);
if ~exist('lims','var') || isempty(lims)
    % Default display limits
    lims = [min(X(:)) max(X(:))];
end

% Create video object (could also use avifile() function)
vidObj = VideoWriter(outpath,format);
vidObj.FrameRate = fps;

% Open the .avi file
vidObj.open();

% Set up a nice figure for recording
scrsz = get(0,'ScreenSize');
xyc = 0.5 * scrsz(3:4);
dim = [nx ny];
f = figure('Name','', ...
           'MenuBar','None', ...
           'NumberTitle','off', ...
           'Position',[(xyc - 0.5 * dim) dim], ...
           'Resize','off');
ax = axes('Position',[0 0 1 1]);

% Record the movie
for k = 1:nt
    % Display the kth frame
    Xk = imresize(X(:,:,k),[res nan],'nearest');
    imshow(Xk,lims,'Border','tight','Parent',ax);
    set(ax,'Position',[0 0 1 1]);
    axis(ax,'off');
    drawnow; pause(0.01);
    
    % Write frame to file
    frame = getframe(f);
    vidObj.writeVideo(frame);
end
close(f);

% Close the .avi file
vidObj.close();

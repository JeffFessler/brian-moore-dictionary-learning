function GenerateMovie(X,outpath,fps,lims)
%--------------------------------------------------------------------------
% Syntax:       GenerateMovie(X,outpath,fps);
%               GenerateMovie(X,outpath,fps,lims);
%               
% Inputs:       X is an (nx x ny x nt) cube containing nt video frames
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
% Date:         March 21, 2014
%--------------------------------------------------------------------------

% Knobs
%format = 'Motion JPEG AVI';
format = 'Uncompressed AVI';
%format = 'Archival';
%format = 'Motion JPEG 2000';
%format = 'MPEG-4';

% Parse inputs
if (nargin < 4)
    lims = [min(X(:)) max(X(:))];
end
[nx ny nt] = size(X);

% Create video object (could also use avifile() function)
vidObj = VideoWriter([outpath '.avi'],format);
vidObj.FrameRate = fps;

% Open the .avi file
vidObj.open();

% Set up a nice figure for recording
scrsz = get(0,'ScreenSize');
x0 = scrsz(3);
y0 = scrsz(4);
f = figure;
ax = gca;
set(f,'Name','');
set(f,'MenuBar','None');
set(f,'NumberTitle','off');
set(f,'Position',[round((x0 - ny) / 2) round((y0 - nx) / 2) ny nx]);
set(f,'Renderer','zbuffer');
set(f,'Color','w');
set(f,'Resize','off');
set(ax,'Position',[0 0 1 1]);
axis(ax,'off');

% Record the movie
for k = 1:nt
    % Display the kth frame
    imshow(X(:,:,k),lims,'Border','tight','Parent',ax);
    drawnow;
    
    % Write frame to file
    frame = getframe(f);
    vidObj.writeVideo(frame);
end
close(f);

% Close the .avi file
vidObj.close();

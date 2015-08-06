% plot 1d fd model data
function FD1Dfnc(frames, path , filename)
close all;

msg = nargchk(1,3,nargin);
if(nargin < 3)
    filename = 'u.dat.';
end
if(nargin < 2)
    path = '/Users/davidson/XCODE-Projects/FD1D/FD1D/solutions/';
end

% get total number of files and data
data = textread(strcat(path,'file_count.txt'));
files = data(1)-1;
if frames > files
    frames = files;
end
C = data(2);
div = data(3);

dt = C/(div);

% new figure
figure(1)
file = sprintf('%s%s%0.3d',path,filename,0);
data = textread(file);
h = plot(data(:,1), data(:,2));
tit = sprintf('time: %2.4f\nCFL: %2.3f',0,C);
ht = title(tit);
set(ht,'fontsize',20);

if(frames == 0), return; end
    
if frames == 1; frames = 2; end
del_files = floor(files/frames);
for i = del_files:del_files:files
    pause(0.005);
    file = sprintf('%s%s%0.3d',path,filename,i);
    data = textread(file);
    % plot the result
    set(h,'ydata',data(:,2));
    tit = sprintf('time: %2.4f\nCFL: %2.3f',dt*i,C);
    set(ht,'string',tit);
    axis ([0 1 -0.1 0.1]);
    grid on
%     if(i == idx || i == idx+1) pause(2); end
end

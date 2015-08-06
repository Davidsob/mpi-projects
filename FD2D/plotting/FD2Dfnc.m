% plot 1d fd model data
function FD2Dfnc(frames, path , filename)
close all;

msg = nargchk(1,3,nargin);
if(nargin < 3)
    filename = 'u.dat.';
end
if(nargin < 2)
    path = '/Users/davidson/XCODE-Projects/FD1D/FD1D/solutions/';
end

% get total number of files and data
dat = textread(strcat(path,'file_count.txt'));
files = dat(1)-1;
C = dat(2);
divx = dat(3);
divy = dat(4);
nx = divx+1;
ny = divy+1;
hx = 1/divx;
hy = 1/divy;
dt = C*min(hx,hy);

x = linspace(0,1,divx+1);
y = linspace(0,1,divy+1);
[X,Y] = meshgrid(x,y);
% new figure
figure(1)
dat = sprintf('%s%s%0.3d',path,filename,0);
d = textread(dat);
Z = reshape(d(:,3),nx,ny);
h = surf(X,Y,Z);
tit = sprintf('time: %2.4f\nCFL: %2.3f',0,C);
ht = title(tit);
set(ht,'fontsize',20);
mxz = max(max(Z));
axis ([0 1 0 1 -2*mxz 2*mxz]);
daspect([1 1 1])
if frames == 1; frames = 2; end
del_files = round(files/(frames-1));
for i = (1+del_files):del_files:files
    pause(0.1);
    file = sprintf('%s%s%0.3d',path,filename,i);
    d = textread(file);
    Z = reshape(d(:,3),nx,ny);
    % plot the result
    set(h,'zdata',Z);
    tit = sprintf('time: %2.4f\nCFL: %2.3f',dt*i,C);
    set(ht,'string',tit);
    axis ([0 1 0 1 -2*mxz 2*mxz]);
    grid on
end

end



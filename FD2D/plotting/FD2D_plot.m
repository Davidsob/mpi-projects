% plot 1d fd model data
close all; clear all; clc;
speed = 1;
filename = 'u.dat.';

path = '/Users/davidson/XCODE-Projects/FD1D/FD1D/solutions/';


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

time = textread(strcat(path,'time.dat'));

x = linspace(0,1,divx+1);
y = linspace(0,1,divy+1);
[X,Y] = meshgrid(x,y);
% new figure
figure(1)
dat = sprintf('%s%s%0.3d',path,filename,0);
d = textread(dat);
Z = reshape(d(:,3),nx,ny);
h = surf(X,Y,Z);
tit = sprintf('time: %2.4f\nCFL: %2.3f',time(1),C);
ht = title(tit);
set(ht,'fontsize',20);
daspect([1 1 1])
% t_c = 0.5;
% idx = round(t_c/dt);
for i = 1:speed:files
    pause(0.005);
    dat = sprintf('%s%s%0.3d',path,filename,i);
    d = textread(dat);
    Z = reshape(d(:,3),nx,ny);
    % plot the result
    set(h,'zdata',Z);
    tit = sprintf('time: %2.4f\nCFL: %2.3f',dt*i,C);
    set(ht,'string',tit);
    axis ([0 1 0 1 -5 5]);
    grid on
    %     if(i == idx || i == idx+1) pause(2); end
end




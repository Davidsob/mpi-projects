% plot a partion of a grid (color, or passing topology)
close all; clear all; clc;

%% open file
path = '/Users/LNLB/mpi-projects/FD2D/';
file = 'grids/graph.dat';

% format ix, jy, nx, ny, color,
partition = textread(strcat(path,file));

% get number of files
file = 'solutions/model_info.dat';
dat = textread(strcat(path,file));
files = dat(1);
parts = dat(2);
C = dat(3);

filename = 'solutions/u.dat.';
fno = 0;
data = {parts,1};
for i = 0:parts-1
    file = sprintf('%s%s%0.3d.proc_%d',path,filename,fno,i);
    data{i+1} = textread(file);
end

time = textread(strcat(path,'solutions/time.dat'));

%% plot the data using a unit step in x, y, (z)
movie = 0;
if(movie)
    writer = avifile('heat.avi');
    writer.Quality = 100;
end

color_by_part = 0;
figure(1)

hndls = {parts,1};

data_idx = 3;

for i = 1:1:parts
    nr = partition(i,3)+1;
    nc = partition(i,4)+1;
    dat = data{i};
    X = reshape(dat(:,1),nc,nr);
    Y = reshape(dat(:,2),nc,nr);
    Z = reshape(dat(:,data_idx),nc,nr);
    
    if color_by_part ~= 1
        hndls{i} = surf(X,Y,Z);
    else
        hndls{i} = surf(X,Y,Z,ones(nc,nr)*i);
    end
    
    hold on;
end

% daspect([1 1 1]);
tit = sprintf('2D-Wave: C = %1.2g, t = %1.3f[s]',C,time(1));
ht = title(tit);
set(ht,'fontsize',20);
if(movie)
    frame = getframe;
    writer = addframe(writer,frame);
end
frames = 1;
if(frames == 0), return; end

%% get extents of surface plot
mins = [1e8, 1e8, 1e8];
maxs = [-1e8, -1e8, -1e8];
for i = 1:parts
    mins_i = [min(min(get(hndls{i},'xdata'))), ...
        min(min(get(hndls{i},'ydata'))), ...
        min(min(get(hndls{i},'zdata')))];
    
    maxs_i = [max(max(get(hndls{i},'xdata'))), ...
        max(max(get(hndls{i},'ydata'))), ...
        max(max(get(hndls{i},'zdata')))];
    
    flag = mins_i < mins;
    mins(flag) = mins_i(flag);
    
    
    flag = maxs_i > maxs;
    maxs(flag) = maxs_i(flag);
    
end

extents = ones(1,6);
extents(1:2:end) = -max(abs(mins),abs(maxs));
extents(2:2:end) = max(abs(mins),abs(maxs));

if(extents(5) == extents(6))
    extents(5) = extents(5)-1;
    extents(6) = extents(6)+1;
end

axis([0 1 0 1 extents(5:6)]);
% axis([0 1 0 1 0 1e6]);
%% make a movie!
for i = 2:length(time)
    pause(0.1);
    for j = 1:parts
        file = sprintf('%s%s%0.3d.proc_%d',path,filename,i-1,j-1);
        nr = partition(j,3)+1;
        nc = partition(j,4)+1;
        dat = textread(file);
        Z = reshape(dat(:,data_idx),nc,nr);
        set(hndls{j},'zdata',Z);
    end
    % plot the result
    tit = sprintf('2D-Wave: C = %1.2g, t = %1.3f[s]',C,time(i));
    set(ht,'string',tit);
    grid on
    
    if(movie)
        frame = getframe;
        writer = addframe(writer,frame);
    end
end

if(movie)
writer = close(writer);
end


% plot a partion of a grid (color, or passing topology)
close all; clear all; clc;

%% open file
path = '/Users/davidson/MPI-Tutorials/FD2D/grids/';
file = 'graph.dat';

% format ix, jy, nx, ny, color,
data = textread(strcat(path,file));

%% plot the data using a unit step in x, y, (z)

figure(1)
for i = 1:size(data,1);
   idx = data(i,1);
   jdx = data(i,2);
   nx = data(i,3);
   ny = data(i,4);
   c = data(i,5);

   plotGridPart(idx,jdx,nx,ny,c);
    
end

daspect([1 1 1]);

figure(2)
for i = 1:size(data,1);
   idx = data(i,1);
   jdx = data(i,2);
   nx = data(i,3);
   ny = data(i,4);
   c = data(i,5);

   plotGridPart(idx,jdx,nx,ny,i);
    
end

daspect([1 1 1]);
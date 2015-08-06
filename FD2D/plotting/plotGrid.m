function plotGridPart(x,y,nx,ny,color)

C = ones(1,4)*color;
for i = 1:nx
    for j = 1:ny
        X = [x + (i-1), x + i, x + i, x + (i-1)];
        Y = [y + (j-1), y + (j-1), y + j, y + j];
        fill(X,Y,color);
    end
end
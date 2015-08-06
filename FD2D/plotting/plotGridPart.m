function plotGridPart(x,y,nx,ny,color)

C = ones(1,4)*color;
for i = 0:nx-1
    for j = 0:ny-1
        X = [x + i, x + i+1, x + i+1, x + i];
        Y = [y + j, y + j, y + j+1, y + j+1];
        fill(Y,X,color);
        hold on;
    end
end
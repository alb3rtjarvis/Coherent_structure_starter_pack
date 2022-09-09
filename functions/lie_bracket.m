function [lb] = lie_bracket(x,y,X,Y,cx,cy)
% - Compute Lie bracket of two vector fields with scalings cx,cy - %
nx = length(x);
ny = length(y);

if ~exist('cx','var')
   cx = ones(ny,nx);
end
if ~exist('cy','var')
    cy = ones(ny,nx);
end

for i = 2:ny-1
    for j = 2:nx-1
        Jx = [(cx(i,j+1)*X(i,j+1,1)-cx(i,j-1)*X(i,j-1,1))/(x(j+1)-x(j-1)),
            (cx(i+1,j)*X(i+1,j,1)-cx(i-1,j)*X(i-1,j,1))/(y(i+1)-y(i-1));
              (cx(i,j+1)*X(i,j+1,2)-cx(i,j-1)*X(i,j-1,2))/(x(j+1)-x(j-1)),
              (cx(i+1,j)*X(i+1,j,2)-cx(i-1,j)*X(i-1,j,2))/(y(i+1)-y(i-1))]; 
        Jy = [(cy(i,j+1)*Y(i,j+1,1)-cy(i,j-1)*Y(i,j-1,1))/(x(j+1)-x(j-1)),
            (cy(i+1,j)*Y(i+1,j,1)-cy(i-1,j)*Y(i-1,j,1))/(y(i+1)-y(i-1));
              (cy(i,j+1)*Y(i,j+1,2)-cy(i,j-1)*Y(i,j-1,2))/(x(j+1)-x(j-1)),
              (cy(i+1,j)*Y(i+1,j,2)-cy(i-1,j)*Y(i-1,j,2))/(y(i+1)-y(i-1))];
          
        lb(i,j,:) = Jy*(cx(i,j)*[X(i,j,1); X(i,j,2)]) - Jx*(cy(i,j)*[Y(i,j,1); Y(i,j,2)]);    
    end
end
end


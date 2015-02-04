% Filename eel2D_straightswimmer.m
% Created By Amneet Bhalla

clear all;
clc;


%Background Mesh Properties
Ny = 16*4*4*4; Ly = 4;
Nx = 32*4*4*4; Lx = 8;
dy = Ly/Ny; dx = Lx/Nx;

L  = 1;                 %length of the body
Wh = 0.04*L;            %width of the head.
Xh = 0.04;   

% No of points on backbone.
BodyNx = ceil(L/dx);
HeadNx = ceil(0.04/dx);
Coord{BodyNx,2} = [];
NumMatPoints    =  0;

angleRotation  = 0;
RotationMatrix = [ cos(angleRotation) -sin(angleRotation);
                   sin(angleRotation)  cos(angleRotation)];

%store from nose to head of the body.
for i = 1:HeadNx
     
    x                 = (i-1)*dx;
    y                 = 0.125* ((x + 0.03125)/1.03125)*sin(2*pi*x);
    height            = sqrt(2*Wh*x - x^2);
    NumPointsInHeight = ceil(height/dy);
    
    ycoord_up = []; ycoord_down = []; xcoord_up = []; xcoord_down = [];
    for j = 1:NumPointsInHeight
       
        xshifted        = x - 0.5;
        yshifted        = y;
        RotatedCoord    = RotationMatrix*[xshifted;yshifted];
        
        
        xcoord_up(j)   = RotatedCoord(1) -(j-1)*dy*sin(angleRotation);
        xcoord_down(j) = RotatedCoord(1) + j*dy*sin(angleRotation);
        ycoord_up(j)   = RotatedCoord(2) + (j-1)*dy*cos(angleRotation);
        ycoord_down(j) = RotatedCoord(2) -j*dy*cos(angleRotation);
        NumMatPoints   = NumMatPoints+2;
    end
    
    Coord{i,1} = cat(2,xcoord_up,xcoord_down);
    Coord{i,2} = cat(2,ycoord_up,ycoord_down);
    
end


%store from head to tail of the body.
for i = HeadNx+1 : BodyNx
    
    x     = (i-1)*dx;
    y     = 0.125* ((x + 0.03125)/1.03125)*sin(2*pi*x);
    height            = Wh*(L - x)/(L - Xh);
    NumPointsInHeight = ceil(height/dy);
    
    ycoord_up = []; ycoord_down = []; xcoord_up = []; xcoord_down = [];
    for j = 1:NumPointsInHeight
        
        xshifted        = x - 0.5;
        yshifted        = y;
        RotatedCoord    = RotationMatrix*[xshifted;yshifted];
        
        xcoord_up(j)   = RotatedCoord(1) -(j-1)*dy*sin(angleRotation);
        xcoord_down(j) = RotatedCoord(1) + j*dy*sin(angleRotation);
        ycoord_up(j)   = RotatedCoord(2) + (j-1)*dy*cos(angleRotation);
        ycoord_down(j) = RotatedCoord(2) -j*dy*cos(angleRotation);
        NumMatPoints = NumMatPoints+2;
    end
    
    Coord{i,1} = cat(2,xcoord_up,xcoord_down);
    Coord{i,2} = cat(2,ycoord_up,ycoord_down);
    
end
    

% plot the fish
for i = 1:size(Coord,1)
    
    plot(Coord{i,1}(:),Coord{i,2}(:),'.');
    hold on;
end
    

% write the coordinates into txt file.
fid = fopen('./eel2d.vertex','wt');
fprintf(fid,'%d\n',NumMatPoints);

for i = 1:size(Coord,1)
    
    for j = 1: length( Coord{i,1}  )
        
         fprintf(fid,'%f\t%f\n', Coord{i,1}(j), Coord{i,2}(j) );
    end

end

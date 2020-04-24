%UNDISTORT unwrap part of the image onto a plane perpendicular to the
%camera axis
%   B = UNDISTORT(OCAM_MODEL, A, FC, DISPLAY)
%   A is the input image
%   FC is a factor proportional to the distance of the camera to the plane;
%   start with FC=5 and then tune the parameter to change the result.
%   DISPLAY visualizes the output image if set to 1; its default value is
%   0.
%   B is the final image
%   Note, this function uses nearest neighbour interpolation to unwrap the
%   image point. Better undistortion methods can be implemented using
%   bilinear or bicub interpolation.
%   Note, if you want to change the size of the final image, change Nwidth
%   and Nheight
%   Author: Davide Scaramuzza, 2009

function Nimg = undistort( ocam_model, img , display, num_ims)

% Parameters of the new image
Nwidth = size(img,1).*2; %size of the final image
Nheight = size(img,2);

if ~isfield(ocam_model,'pol') 
    width = ocam_model.width;
    height = ocam_model.height;
    %The ocam_model does not contain the inverse polynomial pol
    ocam_model.pol = findinvpoly(ocam_model.ss,sqrt((width/2)^2+(height/2)^2));
end

if nargin < 3
    fc = 5;%distance of the plane from the camera, change this parameter to zoom-in or out
    display = 0;
end
    
 %if length(size(img)) == 3
    if num_ims == 1
        Nimg = zeros(Nheight, Nwidth, 3);
    else
        Nimg = zeros(Nheight, Nwidth, 3, num_ims);
    end
% else
%     Nimg = zeros(Nheight, Nwidth);
% end

[i,j] = meshgrid(1:Nheight,1:Nwidth);
theta = 2.*pi.*(j./Nwidth)+pi;
phi = (pi./2)-((i./Nheight).*pi);
Nx = -sin(phi);
Ny = cos(phi).*sin(theta);
Nz = -cos(phi).*cos(theta);
M = [Nx(:)';Ny(:)';Nz(:)'];

m = world2cam_fast( M , ocam_model );

if num_ims == 1
    if length(size(img)) == 2
        I(:,:,1) = img;
        I(:,:,2) = img;
        I(:,:,3) = img;
    elseif length(size(img)) == 3
        I(:,:,1) = img(:,:,1);
        I(:,:,2) = img(:,:,2);
        I(:,:,3) = img(:,:,3);
    end
else
    I(:,:,1,:) = img(:,:,1,:);
    I(:,:,2,:) = img(:,:,2,:);
    I(:,:,3,:) = img(:,:,3,:);
end

[r,g,b] = get_color_from_imagepoints( I, m', num_ims);
if num_ims == 1
    Nimgr = reshape(r,Nwidth,Nheight)';
    Nimgg = reshape(g,Nwidth,Nheight)';
    Nimgb = reshape(b,Nwidth,Nheight)';
    Nimg = uint16(Nimg);
    Nimg(:,:,1) = Nimgr;
    Nimg(:,:,2) = Nimgg;
    Nimg(:,:,3) = Nimgb;
else
    Nimgr = rot90(reshape(r,Nwidth,Nheight,num_ims),3);
    Nimgg = rot90(reshape(g,Nwidth,Nheight,num_ims),3);
    Nimgb = rot90(reshape(b,Nwidth,Nheight,num_ims),3);
    Nimg = uint16(Nimg);
    Nimg(:,:,1,:) = Nimgr;
    Nimg(:,:,2,:) = Nimgg;
    Nimg(:,:,3,:) = Nimgb;
end
    

if display
    figure(1); imshow(Nimg);
end

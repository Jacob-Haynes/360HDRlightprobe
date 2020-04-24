function result_im = merge360(Im, cal_lower, cal_upper, cal_R, cal_col, num_ims)
%% split
Im1 = Im(1:cal_lower.height,:,:); %lower
Im2 = Im(cal_lower.height+1:end,:,:); %upper

%% undistort fisheye
UndistIm1 = undistort(cal_lower, Im1 , 0, num_ims);
UndistIm2 = undistort(cal_upper, Im2, 0, num_ims);
UndistIm2 = flip(UndistIm2,2);

%% rotate 2nd image to align with first
rotated_im = UndistIm1;
pix_perTheta = size(UndistIm1,2)/360;
pix_perPhi = size(UndistIm2,1)/180;
%mesh
[i,j] = meshgrid(1:size(UndistIm2,1),1:size(UndistIm2,2));
%center
ci = i'-(size(UndistIm2,1)/2);
cj = j'-(size(UndistIm2,2)/2);
%radians
theta = deg2rad(cj./pix_perTheta);
phi = deg2rad(ci./pix_perTheta);
%cartesian
[x,y,z]=sph2cart(theta,phi,1);
%rotate
temp = [x(:),y(:),z(:)]*cal_R.';
sz = size(x);
xrot=reshape(temp(:,1),sz);
yrot=reshape(temp(:,2),sz);
zrot=reshape(temp(:,3),sz);
%spherical
[theta,phi,~] = cart2sph(xrot,yrot,zrot);
%new index
Nj = rad2deg(theta)*pix_perTheta + (size(UndistIm2,2)/2);
Ni = rad2deg(phi)*pix_perPhi + (size(UndistIm2,1)/2);
%interp
rotated_im(:,:,1) = interp2(double(UndistIm2(:,:,1)),Nj,Ni);
rotated_im(:,:,2) = interp2(double(UndistIm2(:,:,2)),Nj,Ni);
rotated_im(:,:,3) = interp2(double(UndistIm2(:,:,3)),Nj,Ni);

%% colour correction
lab = rgb2lab(rotated_im);
cols = reshape(lab, [],3);
new_cols = double(cols)*cal_col;
col_cor_im_lab = reshape(new_cols, size(rotated_im,1),size(rotated_im,2),3);
col_cor_im = lab2rgb(col_cor_im_lab);

%convert Im1 to lab and back to rgb to scale it right as it scales it down
%to decimal rather than thousands.
UndistIm1_lab = rgb2lab(UndistIm1);
UndistIm1_rgb = lab2rgb(UndistIm1_lab);
%% blend the two images
% create mask for Im1 - blur out everything 0-90 270-360
mask1 = zeros(size(UndistIm1_rgb));
x1 = 90*pix_perTheta;
x2 = 270*pix_perTheta;
grad1 = linspace(0,1,300);
A = repmat(grad1,size(mask1,1),1,3);
B = flip(A,2);
mask1(:,x1:x2,:) = 1;
mask1(:,(x1-150):(x1+149),:)=A;
mask1(:,(x2-150):(x2+149),:)=B;
%aply mask
dubIm1 = double(UndistIm1_rgb);
maskedIm1 = dubIm1.*mask1;
% create mask for Im2 - blur out center 90-270
mask2 = ones(size(col_cor_im))-mask1;
%aply mask
dubIm2 = double(col_cor_im);
maskedIm2 = dubIm2.*mask2;
%% result fuse
result_im = maskedIm1+maskedIm2;
end
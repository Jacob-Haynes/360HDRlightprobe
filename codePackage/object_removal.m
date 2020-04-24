function result = object_removal(lower, upper, R)
%% properties
clc
fprintf('Loading images... \n');
%load min of 3 images
num_ims = 3; %default
prompt = 'How many images are you using? (default = 3)';
str = input(prompt);
if isempty(str)
    num_ims = 3;
else
    num_ims = str;
end
for i = 1:num_ims
    prompt = 'Image name: ';
    filename = input(prompt, 's');
    im(:,:,:,i) = imread(filename);
end
%% split images
ImLower = im(1:lower.height,:,:,:); %lower
ImUpper = im(lower.height+1:end,:,:,:); %upper
%% undistort fisheye
fprintf('Undistorting images... \n');
UndistIm1 = undistort(lower, ImLower , 0, num_ims);
UndistIm2 = undistort(upper, ImUpper, 0, num_ims);
UndistIm2 = flip(UndistIm2,2);
imageSize = size(UndistIm1);
clear ImLower ImUpper
%% Colour calibration
% match points in each image
for i = 1:num_ims
%     fprintf('Correcting colours for image %i... \n', i)
%     gray1 = histeq(imadjust(im2double(rgb2gray(UndistIm1(:,:,:,i)))));
%     gray2 = histeq(imadjust(im2double(rgb2gray(UndistIm2(:,:,:,i)))));
%     points1 = detectSURFFeatures(gray1);%, 'ROI', [1,size(gray1,1)./6, size(gray1,2), size(gray1,1)./3.*2]);
%     points2 = detectSURFFeatures(gray2);%, 'ROI', [1,size(gray2,1)./6, size(gray2,2), size(gray2,1)./3.*2]);
%     [feat1, points1] = extractFeatures(gray1, points1);
%     [feat2, points2] = extractFeatures(gray2, points2);
%     indexPairs = matchFeatures(feat1, feat2, 'Unique', true);
%     matched1 = points1(indexPairs(:,1),:).Location;
%     matched2 = points2(indexPairs(:,2),:).Location;
%     inlirs = [1:length(matched1)];
%     col(:,:,i) = calibrate_col(UndistIm1, UndistIm2, inlirs, matched1, matched2);
    col(:,:,i) = eye(3);
end
clear UndistIm1 UndistIm2 inlirs gray1 gray2 points1 points2 feat1 feat2 indexPairs matched1 matched2
%% Merge
R = [-0.99591,-0.089274,-0.01383;0.089566,-0.99574,-0.022112;-0.011797,-0.02326,0.99966];
 %debugging
for i = 1:num_ims
    fprintf('Merging Image %i... \n', i);
    merged(:,:,:,i) = merge360(im(:,:,:,i),lower,upper,R,col(:,:,i),1);
end
%% Align images to image 1
    gray1 = histeq(imadjust(im2double(rgb2gray(merged(:,:,:,1)))));
    points1 = detectSURFFeatures(gray1);%, 'ROI', [1,size(gray1,1)./6, size(gray1,2), size(gray1,1)./3.*2]);
    [feat1, points1] = extractFeatures(gray1, points1);
    rotated_im = merged(:,:,:,:);
for i = 2:num_ims
    fprintf('Algining image %i to 1... \n',i);
    %match feats
    gray2 = histeq(imadjust(im2double(rgb2gray(merged(:,:,:,i)))));
    points2 = detectSURFFeatures(gray2);%, 'ROI', [1,size(gray2,1)./6, size(gray2,2), size(gray2,1)./3.*2]);
    [feat2, points2] = extractFeatures(gray2, points2);
    indexPairs = matchFeatures(feat1, feat2, 'Unique', true);
    matched1 = points1(indexPairs(:,1),:).Location;
    matched2 = points2(indexPairs(:,2),:).Location;
    %convert 3d
    pix_perTheta = imageSize(2)/360;
    pix_perPhi = imageSize(1)/180;
    %add r
    pad1 = ones(length(matched1),1);
    matched1_3D_flat = [matched1,pad1];
    matched2_3D_flat = [matched2,pad1];
    % center
    center(1) = size(gray1,2)./2;
    center(2) = size(gray1,1)./2;
    center(3) = 0;
    c1 = matched1_3D_flat - center;
    c2 = matched2_3D_flat - center;
    %convert to rad
    d1 = [c1(:,2)./pix_perPhi,c1(:,1)./pix_perTheta,c1(:,3)];
    d2 = [c2(:,2)./pix_perPhi,c2(:,1)./pix_perTheta,c2(:,3)];
    angle1 = [deg2rad(d1(:,1:2)),d1(:,3)];
    angle2 = [deg2rad(d2(:,1:2)),d2(:,3)];
    %sph 2 cart
    [x1,y1,z1] = sph2cart(angle1(:,2), angle1(:,1),angle1(:,3));
    [x2,y2,z2] = sph2cart(angle2(:,2),angle2(:,1),angle2(:,3));
    matched1_3D = [x1,y1,z1];
    matched2_3D = [x2,y2,z2];
    %ransac
    [A(:,:,i), ~] = RANSAC(matched1_3D, matched2_3D, 'pos');
    % rotate 2nd image to align with first
    pix_perTheta = size(merged,2)/360;
    pix_perPhi = size(merged,1)/180;
    %mesh
    [a,b] = meshgrid(1:size(merged,1),1:size(merged,2));
    %center
    ci = a'-(size(merged,1)/2);
    cj = b'-(size(merged,2)/2);
    %radians
    theta = deg2rad(cj./pix_perTheta);
    phi = deg2rad(ci./pix_perTheta);
    %cartesian
    [x,y,z]=sph2cart(theta,phi,1);
    %rotate
    temp = [x(:),y(:),z(:)]*A(:,:,i).';
    sz = size(x);
    xrot=reshape(temp(:,1),sz);
    yrot=reshape(temp(:,2),sz);
    zrot=reshape(temp(:,3),sz);
    %spherical
    [theta,phi,~] = cart2sph(xrot,yrot,zrot);
    %new index
    Nj = rad2deg(theta)*pix_perTheta + (size(merged,2)/2);
    Ni = rad2deg(phi)*pix_perPhi + (size(merged,1)/2);
    %interp
    rotated_im(:,:,1,i) = interp2(double(merged(:,:,1,i)),Nj,Ni);
    rotated_im(:,:,2,i) = interp2(double(merged(:,:,2,i)),Nj,Ni);
    rotated_im(:,:,3,i) = interp2(double(merged(:,:,3,i)),Nj,Ni);
end
%% combine the images taking the median pixel value
result = median(rotated_im, 4);
end

function calibrate360
%% properties
clear
clc
%defaults
fileName = 'Cal_A.tiff';
autofeatures = 1;
surf = 1;
properties_set = 0;
%get preferences
while properties_set == 0
    fprintf('Current properties: \n Calibration file = "%s" \n auto feature selection = "%i" \n surf features = "%i" \n', fileName,autofeatures,surf);
    prompt = 'Do you want to change properties? Y/N [default = N]: ';
    str = input(prompt,'s');
    if isempty(str)
        str = 'N';
    end
    if str == 'N'
        properties_set = 1;
    else
    prompt = 'Calibration file: ';
    fileName = input(prompt,'s');
    prompt = 'Auto feature selection: ';
    autofeatures = input(prompt);
    prompt = 'Surf features: ';
    surf = input(prompt);
    end
end
%% Get camera calibration paramaters
%grid squares x=7 y=4 x=32mm y=29mm,  10 images
% from ocam_calib
% polynomial coefficients ss, invpol
% centers xc yc
% affine parameters "c", "d", "e"
% lower 5 15
lower.ss = [ -9.343383e+02 0.000000e+00 2.546344e-04 -2.450943e-08 9.461811e-11  ];
lower.invpol = [ 1481.140014 812.598208 -58.068926 152.720328 50.274832 -4.348677 50.556259 -24.269045 -10.887808 52.821962 3.440259 -29.792634 -2.530701 8.810030 2.566491 ];
lower.xc = 1725.904819;
lower.yc = 1728.486293;
lower.c = 0.993223;
lower.d = 0.000077;
lower.e = -0.000119;
lower.height = 3456;
lower.width = 3456;
% upper 5 16
upper.ss = [-8.706589e+02 0.000000e+00 3.503911e-04 -1.414086e-07 1.222870e-10 ];
upper.invpol = [ 1464.982183 885.649072 -82.652189 119.590497 99.176419 -53.771453 58.137939 43.968097 -84.205343 7.232620 77.777963 -6.422891 -38.232782 -2.010324 9.544693 2.576489 ];
upper.xc = 1725.066534;
upper.yc = 1725.747413;
upper.c = 0.991498;
upper.d = 0.001621;
upper.e = 0.000991;
upper.height = 3456;
upper.width = 3456;
%FOV
f = 210;
%% load
% image file
combineIm = imread(fileName);
%automatic features or manual
auto = autofeatures;
%% Split the two images
Im1 = combineIm(1:lower.height,:,:); %lower
Im2 = combineIm(lower.height+1:end,:,:); %upper

%% undistort fisheye
fprintf('Undistorting image... \n');
UndistIm1 = undistort(lower, Im1 , 0,1);
UndistIm2 = undistort(upper, Im2, 0,1);
UndistIm2 = flip(UndistIm2,2);
imageSize = size(UndistIm1);
clear Im1 Im2
%% Feature match
fprintf('Finding features... \n');
gray1 = histeq(imadjust(im2double(rgb2gray(UndistIm1))));
gray2 = histeq(imadjust(im2double(rgb2gray(UndistIm2))));
% SURF features
if surf == 1
    %extract features in the relevent areas of the image
    %swap lower L and R as its mirrored...
    points1 = detectSURFFeatures(gray1);%, 'ROI', [1,size(gray1,1)./6, size(gray1,2), size(gray1,1)./3.*2]);
    points2 = detectSURFFeatures(gray2);%, 'ROI', [1,size(gray2,1)./6, size(gray2,2), size(gray2,1)./3.*2]);
    [feat1, points1] = extractFeatures(gray1, points1);
    [feat2, points2] = extractFeatures(gray2, points2);
    %match features
    indexPairs = matchFeatures(feat1, feat2, 'Unique', true);
    matched1 = points1(indexPairs(:,1),:).Location;
    matched2 = points2(indexPairs(:,2),:).Location;
    %manual select
    if auto == 0
        good_points = 0;
        while good_points == 0
            chosenIdx = select_matched_points(matched1,matched2,gray1,gray2);
            showMatchedFeatures(gray2,gray1,matched2(chosenIdx),matched1(chosenIdx),'montage')
            prompt = 'Do you want to try again? Y/N [default = N]: ';
            str = input(prompt,'s');
            if isempty(str)
                str = 'N';
            end
            if str == 'N'
                good_points = 1;
            end
        end
        matched1 = matched1(chosenIdx);
        matched2 = matched2(chosenIdx);
    end
else
   fprintf('Select points in image one. Remember order');
   matched1 = getpoints(gray1);
   fprintf('Select matching points in the same order');
   matched2 = getpoints(gray2);
end

clear points1 points2 feat1 feat2 indexPairs


%% convert to 3D location
fprintf('Converting images to spherical 3D... \n');
pix_perTheta = imageSize(2)/360;
pix_perPhi = imageSize(1)/180;
%add r
pad1 = ones(length(matched1),1);
matched1_3D_flat = [matched1,pad1];
matched2_3D_flat = [matched2,pad1];
% center
center(1) = size(UndistIm2,2)./2;
center(2) = size(UndistIm2,1)./2;
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

clear matched1_3D_flat matched2_3D_flat center c1 c2 d1 d2 angle1 angle2 x1 x2 y1 y2 z1 z2 gray1 gray2

%% RANSAC rotation matrix
fprintf('Performing RANSAC for rotation matrix... \n');
[R, inlirs] = RANSAC(matched1_3D, matched2_3D, 'pos');
clear matched1_3D matched2_3D

%% colour calibration
fprintf('Performing colour calibration... \n');
col = calibrate_col(UndistIm1, UndistIm2, inlirs, matched1, matched2);

%% check with image
fprintf('Merging Images and correcting colours... \n');
result = merge360(combineIm,lower,upper,R,col,1);
fprintf('DONE.\n R = \n');
disp(R)
fprintf('Colour calibration = \n');
disp(col)
figure(1)
imshow(result);
figure(2)
imshow(histeq(result));

%% save calibration values
cal_ok = 0;
while cal_ok == 0
    prompt = 'Was the callibration successful? Y/N [default = N]: ';
    str = input(prompt,'s');
    if str == 'Y'
        cal_ok = 1;
    else
        calibrate360
    end
end
if cal_ok ==1
    csvwrite('R_calibration.csv', R);
end
end

%part 2 use 3 images with me in different palces
% form the combined image for all of them
% when brightening post alignment only use features that are super close
% after alignment to not use matches from me
% could use the matches from me to define a region to be removed via
% gradient mask or use the median value from the 3 images (post
% normalising)

%try using a 'calibration' image to get R then use the same R for rotating
%a test image
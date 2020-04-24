function lightprobe
clear
clc
%defaults
fileName = 'Cal_A.tiff';
R = csvread('R_calibration.csv');
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

%check file
properties_set = 0;
while properties_set == 0
    fprintf('Current image to be processed: \n Calibration file = "%s" \n ', fileName,autofeatures,surf);
    prompt = 'Do you want to change file? Y/N [default = N]: ';
    str = input(prompt,'s');
    if isempty(str)
        str = 'N';
    end
    if str == 'N'
        properties_set = 1;
    else
    prompt = 'image file name: ';
    fileName = input(prompt,'s');
    end
end

%% load image
fprintf('Loading image... \n');
Im = imread(fileName);

%% colour match
fprintf('Matching colours... \n');
% split
Im1 = Im(1:lower.height,:,:); %lower
Im2 = Im(lower.height+1:end,:,:); %upper
% undistort fisheye
UndistIm1 = undistort(lower, Im1 , 0);
UndistIm2 = undistort(upper, Im2, 0);
UndistIm2 = flip(UndistIm2,2);
% features
points1 = detectSURFFeatures(gray1);
points2 = detectSURFFeatures(gray2);
[feat1, points1] = extractFeatures(gray1, points1);
[feat2, points2] = extractFeatures(gray2, points2);
%match features
indexPairs = matchFeatures(feat1, feat2, 'Unique', true);
matched1 = points1(indexPairs(:,1),:).Location;
matched2 = points2(indexPairs(:,2),:).Location;
% need some way to get inlirs... perhaps shortest distance on the rotated
col = calibrate_col(UndistIm1, UndistIm2, inlirs, matched1, matched2)
%% merge
fprintf('Merging... \n');
result_im = merge360(Im, lower, upper, R, col);

imshow(result_im);

end
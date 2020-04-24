function col = calibrate_col(UndistIm1, UndistIm2, inlirs, matched1, matched2)
%% Get colour vals of matched positions
Lab_1 = rgb2lab(UndistIm1);
Lab_2 = rgb2lab(UndistIm2);
col_1(:,1) = interp2(double(Lab_1(:,:,1)),matched1(inlirs,1),matched1(inlirs,2));
col_1(:,2) = interp2(double(Lab_1(:,:,2)),matched1(inlirs,1),matched1(inlirs,2));
col_1(:,3) = interp2(double(Lab_1(:,:,3)),matched1(inlirs,1),matched1(inlirs,2));

col_2(:,1) = interp2(double(Lab_2(:,:,1)),matched2(inlirs,1),matched2(inlirs,2));
col_2(:,2) = interp2(double(Lab_2(:,:,2)),matched2(inlirs,1),matched2(inlirs,2));
col_2(:,3) = interp2(double(Lab_2(:,:,3)),matched2(inlirs,1),matched2(inlirs,2));

%% RANSAC
[col, inlirs] = RANSAC(col_1, col_2, 'col');

end
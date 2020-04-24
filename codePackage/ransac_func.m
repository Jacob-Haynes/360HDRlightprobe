function [R inindex] = ransac_func(x, y, n, k, t, d, findH, findDist)
%     x,y  observed matched feature points
%     n - min num of points required
%     k - max num of iterations
%     t - threshold value for fitting
%     d - ratio of inliers
num_points=size(x,1);
inlirthresh = round(d*num_points);
inlir_num = zeros(1,k);
Rlib = cell(1,k);

for i = 1:k
    % randomly select sample points
    random_pairs = randperm(size(x,1),n);
    %     sample = data(random_pairs, :);
    R1 = findH(x(random_pairs,:), y(random_pairs,:));
    
    % count inliers, if more than inlier thresh, refit
    Dist = findDist(R1, x, y);
    inlir = find(Dist < t);
    inlir_num(i) = length(inlir);
    if inlir_num(i) > inlirthresh
        Rlib{i} = findH(x(inlir,:),y(inlir,:));
    end
end
%choose the set with most inliers
[~,index] = max(inlir_num);
R = Rlib{index};
empt = isempty(R);
if empt == 1
    R = eye(3);
end
Dist = findDist(R, x, y);
inindex = find(Dist < t);
end
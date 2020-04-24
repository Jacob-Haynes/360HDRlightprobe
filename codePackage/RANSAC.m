function [R ptindx] = RANSAC(points_in, points_out, type)
%find homo between matched points using ransac
%points=[x1,y1;x2,y2; etc]
%ptindx = inlier index
if type == 'pos'
    n = 7; % min num of points required %think 3 for rotation?
    k = 7000; % number of iterations
    t = 0.3; % inlier threshold %how many degrees is the threshold
    d = 0.2; % ratio of inliers
    [R ptindx] = ransac_func(points_in, points_out, n, k, t, d, @homo_estimate, @dist_func);
end
if type == 'col'
    n = 7;
    k = 7000;
    t = 2.3;
    d = 0.3;
    [R ptindx] = ransac_func(points_in, points_out, n, k, t, d, @col_estimate, @col_dif);
end
end

function dist = dist_func(R, points_in, points_out)
points_temp = points_out*R;
dist = zeros(length(points_temp),1);
for i = 1:length(points_temp)
    dist(i) = atan2(norm(cross(points_in(i,:), points_temp(i,:))),dot(points_in(i,:),points_temp(i,:)));
end
end

function R = homo_estimate(points_in, points_out)
    [ctform,~,~] = absor(points_in', points_out','doTrans', 0);
    R = ctform.R;
end

function R = col_estimate(points_in, points_out)
    [ctform,~,~] = absor(points_in',points_out','doTrans',0);
    R = ctform.R;
end

function dist = col_dif(R, points_in, points_out)
points_temp = points_out*R;
deltaE = sqrt((points_temp(:,1)-points_in(:,1)).^2 + (points_temp(:,2)-points_in(:,2)).^2 + (points_temp(:,3)-points_in(:,3)).^2);
dist = deltaE;
end
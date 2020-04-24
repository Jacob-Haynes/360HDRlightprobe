function chosenIdx = select_matched_points(matched1,matched2, image1, image2)
showMatchedFeatures(image2,image1,matched2,matched1,'montage')
hold on
fprintf('click the matches you want to keep on the LEFT image. Press RETURN when done. \n');
[mx,my] = ginput;
selected = zeros(length(mx),2);
for i = 1:length(mx)
    %compute Euclidean distances:
    distances = sqrt(sum(bsxfun(@minus, matched2, [mx(i),my(i)]).^2,2));
    %find the smallest distance and use that as an index into B:
    closest = matched2((distances==min(distances)),:);
    % Give "true" if the element in "a" is a member of "b".
    logical = ismember(matched2,closest);
    % Extract the elements of a at those indexes.
    indexes = find(logical);
    selected(i,:) = indexes';
end
chosenIdx = unique(selected,'rows');
hold off
end
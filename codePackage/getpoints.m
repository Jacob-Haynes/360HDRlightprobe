function points = getpoints(im)
figure(2)
imshow(im);
hold on
fprintf('Click to select points, press RETURN when finished, BACKSPACE removes the previous point.\n');
[xi,yi]=getpts;
points = [xi,yi];
hold off
end

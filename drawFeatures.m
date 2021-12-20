function [] = drawFeatures( img, loc )
% Function: 画出特征点
figure;
imshow(img);
hold on;
plot(loc(:,2),loc(:,1),'+g');
end
function [] = drawMatched( matched, img1, img2, loc1, loc2 )
% Function: 画匹配上的特征点
img3 = appendimages(img1,img2);

% 特征点连线
figure('Position', [100 100 size(img3,2) size(img3,1)]);
colormap gray;
imagesc(img3);
hold on;
cols1 = size(img1,2);
n = size(matched,2);
colors = ['c','m','y'];
colors_n = length(colors);
for i = 1: n
  if (matched(i) > 0)
    color = colors(randi(colors_n));
    line([loc1(i,2) loc2(matched(i),2)+cols1], ...
         [loc1(i,1) loc2(matched(i),1)], 'Color', color);
  end
end
hold off;
num = sum(matched > 0);
fprintf('Found %d matches.\n', num);

end

function im = appendimages(image1, image2)
% im = appendimages(image1, image2)
rows1 = size(image1,1);
rows2 = size(image2,1);

if (rows1 < rows2)
     image1(rows2,1) = 0;
else
     image2(rows1,1) = 0;
end

im = [image1 image2];
end


function [ matched ] = match( des1, des2 )
% Function: 匹配关键点并返回匹配的两个关键点的索引

% 匹配的特征向量从最近邻到第二近邻的角度小于 distRatio？
distRatio = 0.6;

des2t = des2';
n = size(des1,1);
matched = zeros(1,n);
for i = 1 : n
   dotprods = des1(i,:) * des2t;
   [values,index] = sort(acos(dotprods));
   if (values(1) < distRatio * values(2))
      matched(i) = index(1);
   else
      matched(i) = 0;
   end
end

end
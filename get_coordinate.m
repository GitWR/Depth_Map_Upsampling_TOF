function [zuo_b_c,zuo_b_l] = get_coordinate(a,b)
% 获取输入图片的坐标，以用于计算term2的K
zuo_b_c = cell(a,b); % 每个cell单元里面存储的都是一个列向量-->[x,y]'
zuo_b_l = cell(a,b); 

for i = 1 : a
    for j = 1 : b
        zuo_b_c{i,j} = [i,j];
        zuo_b_l{i,j} = [i,j]';
    end
end
end


function [zuo_b_c,zuo_b_l] = get_coordinate(a,b)
% ��ȡ����ͼƬ�����꣬�����ڼ���term2��K
zuo_b_c = cell(a,b); % ÿ��cell��Ԫ����洢�Ķ���һ��������-->[x,y]'
zuo_b_l = cell(a,b); 

for i = 1 : a
    for j = 1 : b
        zuo_b_c{i,j} = [i,j];
        zuo_b_l{i,j} = [i,j]';
    end
end
end


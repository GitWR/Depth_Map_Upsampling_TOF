function [x] = CG_gradient(A,b,x_0,error)
% CG_GRADIENT 
% A 表示系数矩阵
% b 表示值矩阵
% x_0 表示目标图片的初始值
% error 表示迭代停止条件
n = size(A,1); % 最大迭代次数
x = x_0; % 初始化
r = b - A * x; % 初始搜寻方向,r0
d = r; % 后续更新的搜索方向,d0
for i = 0 : (n-1) % 这样可以节省内存
    yita = (r' * r) / (d' * A * d); % 学习率，yita0
    x = x + yita * d; % 更新x,x1 = x0 + yita0 * d0
    r1 = b - A * x; % 更新r, actually is r1 = b - A * x1
    if (norm(r1) <= error || (i == n-1))
        x_h = x;
        break;
    end
    beta = norm(r1)^2 / norm(r)^2; % beta0
    d = r1 + beta * d; % 更新d,d1
    r = r1; % 更新r,r1
    fprintf('第%d次迭代+误差为%6f\n',i,norm(r1));
    plot(norm(r1),i);
    hold on
end
x = x_h;
end


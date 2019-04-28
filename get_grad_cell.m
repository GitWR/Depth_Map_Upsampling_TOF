function [ti_d] = get_grad_cell(I_grad_x,I_grad_y,a,b)
% 用cell把矩阵每个点的梯度组合成 3\times2 的矩阵进行保存
ti_d = cell(a,b);
temp = zeros(3,2); % 每个RGB的像素点的梯度为一个3\times2 的矩阵
for i = 1: a
    for j = 1 : b
        r_x = I_grad_x(:,:,1);
        g_x = I_grad_x(:,:,2);
        b_x = I_grad_x(:,:,3);
        r_y = I_grad_y(:,:,1);
        g_y = I_grad_y(:,:,2);
        b_y = I_grad_y(:,:,3);
        temp = [r_x(i,j),r_y(i,j);g_x(i,j),g_y(i,j);b_x(i,j),b_y(i,j)];
        zegma_pixel = temp' * temp;
        ti_d{i,j} = zegma_pixel;
    end
end
end


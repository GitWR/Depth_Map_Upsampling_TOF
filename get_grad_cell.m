function [ti_d] = get_grad_cell(I_grad_x,I_grad_y,a,b)
% ��cell�Ѿ���ÿ������ݶ���ϳ� 3\times2 �ľ�����б���
ti_d = cell(a,b);
temp = zeros(3,2); % ÿ��RGB�����ص���ݶ�Ϊһ��3\times2 �ľ���
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


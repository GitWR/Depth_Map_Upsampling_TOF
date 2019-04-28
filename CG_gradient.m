function [x] = CG_gradient(A,b,x_0,error)
% CG_GRADIENT 
% A ��ʾϵ������
% b ��ʾֵ����
% x_0 ��ʾĿ��ͼƬ�ĳ�ʼֵ
% error ��ʾ����ֹͣ����
n = size(A,1); % ����������
x = x_0; % ��ʼ��
r = b - A * x; % ��ʼ��Ѱ����,r0
d = r; % �������µ���������,d0
for i = 0 : (n-1) % �������Խ�ʡ�ڴ�
    yita = (r' * r) / (d' * A * d); % ѧϰ�ʣ�yita0
    x = x + yita * d; % ����x,x1 = x0 + yita0 * d0
    r1 = b - A * x; % ����r, actually is r1 = b - A * x1
    if (norm(r1) <= error || (i == n-1))
        x_h = x;
        break;
    end
    beta = norm(r1)^2 / norm(r)^2; % beta0
    d = r1 + beta * d; % ����d,d1
    r = r1; % ����r,r1
    fprintf('��%d�ε���+���Ϊ%6f\n',i,norm(r1));
    plot(norm(r1),i);
    hold on
end
x = x_h;
end


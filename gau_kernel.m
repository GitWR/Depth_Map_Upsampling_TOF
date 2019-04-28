function [kernel] = gau_kernel(f)
kernel = zeros(2 * f + 1,2 * f + 1); % ��fΪ�뾶 
%for d = 1:f
d = f;
  %value = 1 / (2 * d + 1)^2 ; 
  delta = 1; % Ĭ������Ϊ1
  for i = -d : d
     for j = -d : d
        kernel(f-i+1,f-j+1) = (1/(2 * pi)) * (exp(-((i)^2+(j)^2))); % kernel(f+1-i,f+1-j) + value
     end
  end
%end
sum_weight = sum(sum(kernel)); % �ܺ�
kernel = kernel ./ sum_weight; % ��һ��
end


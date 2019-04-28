function [kernel] = gau_kernel(f)
kernel = zeros(2 * f + 1,2 * f + 1); % 以f为半径 
%for d = 1:f
d = f;
  %value = 1 / (2 * d + 1)^2 ; 
  delta = 1; % 默认设置为1
  for i = -d : d
     for j = -d : d
        kernel(f-i+1,f-j+1) = (1/(2 * pi)) * (exp(-((i)^2+(j)^2))); % kernel(f+1-i,f+1-j) + value
     end
  end
%end
sum_weight = sum(sum(kernel)); % 总和
kernel = kernel ./ sum_weight; % 归一化
end


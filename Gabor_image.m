function [Transition] = Gabor_image(A)

% GABOR_IMAGE 
% 对图像进行Gabor滤波
% A 表示给定的color图像
Transition = zeros(size(A,1),size(A,2)); % 输出的每一个map的大小
nscale = 5; % 尺度
norient = 8; % 方向
minWaveLength = 3;
mult = 1.7; % 连续滤波器之间的比例因子
sigmaOnf = 0.65; % standard  derivation
dThetaOnSigma = 1.5;
feedback = 0;
A_gray = double(rgb2gray(uint8(A))); 
A_filtered = Log_gabor(A_gray,nscale, norient, minWaveLength, mult,sigmaOnf,dThetaOnSigma,feedback);
count = 0;
for i = 1 : 5
    for j = 1 : 8
        count = count + 1;
        Transition = Transition + real(A_filtered{i,j}); 
        %figure;
        %imshow(double(Transition));
    end
end

end


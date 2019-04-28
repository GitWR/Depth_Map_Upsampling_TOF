function [Data] = data_process(input_depmap,input_colmap,input_gth)
% 对数据进行裁剪
% 大小为200 * 200
bicu_d = imresize(input_depmap,[1000,1000],'bicubic');
enlarge_depmap = imresize(input_depmap,[1000,1000],'nearest');
input_colmap = imresize(input_colmap,[1000,1000],'nearest');
[a,b] = size(enlarge_depmap);
P_dp = enlarge_depmap(a-600:end-401,200:399);

tic
Edg_color = double(real(Gabor_image(input_colmap)));
toc
disp('Geting the edge information')
 
c1 = input_colmap(:,:,1);
c2 = input_colmap(:,:,2);
c3 = input_colmap(:,:,3);
P_c1 = c1(a-600:end-401,200:399); 
P_c2 = c2(a-600:end-401,200:399);
P_c3 = c3(a-600:end-401,200:399);
P_c(:,:,1) = P_c1;
P_c(:,:,2) = P_c2;
P_c(:,:,3) = P_c3;
P_gt = input_gth(a-600:end-401,200:399);
P_bicu = bicu_d(a-600:end-401,200:399);
P_edg = Edg_color(a-600:end-401,200:399); 
Data.d = P_dp; % 原始上采样后裁剪的深度图
Data.c = P_c; % 裁剪的RGB图
Data.g = P_gt; % 裁剪的gth
Data.b = P_bicu; % 裁剪的depth_guide
Data.e = P_edg; % 裁剪的边界响应
save Patch Data
end


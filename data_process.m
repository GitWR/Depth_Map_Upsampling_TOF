function [Data] = data_process(input_depmap,input_colmap,input_gth)
% �����ݽ��вü�
% ��СΪ200 * 200
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
Data.d = P_dp; % ԭʼ�ϲ�����ü������ͼ
Data.c = P_c; % �ü���RGBͼ
Data.g = P_gt; % �ü���gth
Data.b = P_bicu; % �ü���depth_guide
Data.e = P_edg; % �ü��ı߽���Ӧ
save Patch Data
end


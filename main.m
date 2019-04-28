clc
clear all,
close all,

% 输入图片
% input_depmap = double(imread('data/moebius-depth.png'));
% input_colmap = double(imread('data/moebius-color.png'));
% input_gth = double(imread('data/mobius_big.png'));
% input_gth = imresize(input_gth,[1110,1390]);

% Data = data_process(input_depmap,input_colmap,input_gth); % 裁剪数据
load ('ToFSimulated_2')
% 为了匹配维数的整数关系
% G_ori = DepthSample; %double(Data.d); % lr depth map
% G_ori = G_ori(500:699,500:699);

I = Color; % double(Data.c); % corresponded color map
% I = I(500:699,500:699,:);

% gabor_maps = gen_Gabor_maps(I,5,8,11,11);
% temp = zeros(200,200);
% for c = 1 : 40
%    temp = temp + real(gabor_maps(:,:,c));
% end

G_bcubic = DepthGuide; % double(Data.b); % bicubic 插值结果
G_bcubic = G_bcubic(520:719,120:319); % 500:699,500:699

G_gth = DepthGT; %double(Data.g); % 裁剪后的gth2
G_ori = imresize(G_gth,1/2);
up_sample = zeros(size(G_gth,1),size(G_gth,2));
up_sample(1:2:end,1:2:end) = G_ori;
G_ori = up_sample(520:719,120:319);

G_gth = G_gth(520:719,120:319);

G_edg = Gabor_image(I); % 裁剪后的边界响应
% level = graythresh(G_edg/sum(G_edg(:)));
% bw = im2bw(G_edg/sum(G_edg(:)),level);
% imshow(bw)
%% 相关参数的初始化
win = 4; % 搜索窗口的半径大小 9 * 9
ser = 3; % 邻域窗口的大小 7 * 7
lamda_N = 0.0000008; % 按照文中给的结果0.006 ,这一项不能太大
lamda_s = 0.09; % 同上,0.008 , 5的话很差 0.09
Patch_size = 200;

%% Optimization
D = MRF_D(G_ori,G_bcubic,I,G_edg,Patch_size,lamda_N,lamda_s,win,ser); % 计算结果矩阵D，利用MRF_CG算法. 利用的是Obj中的第一和第三项
hr_depth = reshape(D, Patch_size, Patch_size);

% tic
% temp1 = bsxfun(@minus,hr_depth,G_gth); % 偏差                                                                                                                                                                                            项
% temp2 = bsxfun(@times,temp1,temp1); % 偏差的平方
% temp3 = sum(sum(temp2)); % 平方和
% RMES = sqrt(temp3 / Patch_size.^2);
% toc
% disp('computing the RMES')
% fprintf('均方根误差为%f',RMES);
[m,n] = size(G_gth);
Mask = zeros(m, n);
Mask(G_gth>0) = 1;
G_gth(Mask<1) = 0;
Diff = abs(G_gth - hr_depth).* Mask;
MAD = sum(Diff(:))/sum(Mask(:)); % 每个像素点的平均值误差值
RSME = sqrt(sum(Diff(:).^2)  / sum(Mask(:)));
fprintf('平均误差为%f\n',MAD);
fprintf('平均误差为%f\n',RSME);

imshow(uint8(G_gth)),title('groundtruth');
figure
imshow(uint8(hr_depth)),title('recovery');

%% 实现Segmentation [17]
% ratio = 0.5;
% kernelsize = 2;
% Iseg = vl_quickseg(I, ratio, kernelsize, maxdist)-->maxdist = [10,20];







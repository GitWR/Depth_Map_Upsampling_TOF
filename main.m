clc
clear all,
close all,

% ����ͼƬ
% input_depmap = double(imread('data/moebius-depth.png'));
% input_colmap = double(imread('data/moebius-color.png'));
% input_gth = double(imread('data/mobius_big.png'));
% input_gth = imresize(input_gth,[1110,1390]);

% Data = data_process(input_depmap,input_colmap,input_gth); % �ü�����
load ('ToFSimulated_2')
% Ϊ��ƥ��ά����������ϵ
% G_ori = DepthSample; %double(Data.d); % lr depth map
% G_ori = G_ori(500:699,500:699);

I = Color; % double(Data.c); % corresponded color map
% I = I(500:699,500:699,:);

% gabor_maps = gen_Gabor_maps(I,5,8,11,11);
% temp = zeros(200,200);
% for c = 1 : 40
%    temp = temp + real(gabor_maps(:,:,c));
% end

G_bcubic = DepthGuide; % double(Data.b); % bicubic ��ֵ���
G_bcubic = G_bcubic(520:719,120:319); % 500:699,500:699

G_gth = DepthGT; %double(Data.g); % �ü����gth2
G_ori = imresize(G_gth,1/2);
up_sample = zeros(size(G_gth,1),size(G_gth,2));
up_sample(1:2:end,1:2:end) = G_ori;
G_ori = up_sample(520:719,120:319);

G_gth = G_gth(520:719,120:319);

G_edg = Gabor_image(I); % �ü���ı߽���Ӧ
% level = graythresh(G_edg/sum(G_edg(:)));
% bw = im2bw(G_edg/sum(G_edg(:)),level);
% imshow(bw)
%% ��ز����ĳ�ʼ��
win = 4; % �������ڵİ뾶��С 9 * 9
ser = 3; % ���򴰿ڵĴ�С 7 * 7
lamda_N = 0.0000008; % �������и��Ľ��0.006 ,��һ���̫��
lamda_s = 0.09; % ͬ��,0.008 , 5�Ļ��ܲ� 0.09
Patch_size = 200;

%% Optimization
D = MRF_D(G_ori,G_bcubic,I,G_edg,Patch_size,lamda_N,lamda_s,win,ser); % ����������D������MRF_CG�㷨. ���õ���Obj�еĵ�һ�͵�����
hr_depth = reshape(D, Patch_size, Patch_size);

% tic
% temp1 = bsxfun(@minus,hr_depth,G_gth); % ƫ��                                                                                                                                                                                            ��
% temp2 = bsxfun(@times,temp1,temp1); % ƫ���ƽ��
% temp3 = sum(sum(temp2)); % ƽ����
% RMES = sqrt(temp3 / Patch_size.^2);
% toc
% disp('computing the RMES')
% fprintf('���������Ϊ%f',RMES);
[m,n] = size(G_gth);
Mask = zeros(m, n);
Mask(G_gth>0) = 1;
G_gth(Mask<1) = 0;
Diff = abs(G_gth - hr_depth).* Mask;
MAD = sum(Diff(:))/sum(Mask(:)); % ÿ�����ص��ƽ��ֵ���ֵ
RSME = sqrt(sum(Diff(:).^2)  / sum(Mask(:)));
fprintf('ƽ�����Ϊ%f\n',MAD);
fprintf('ƽ�����Ϊ%f\n',RSME);

imshow(uint8(G_gth)),title('groundtruth');
figure
imshow(uint8(hr_depth)),title('recovery');

%% ʵ��Segmentation [17]
% ratio = 0.5;
% kernelsize = 2;
% Iseg = vl_quickseg(I, ratio, kernelsize, maxdist)-->maxdist = [10,20];







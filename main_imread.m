clc
clear all,
close all,

% ����ͼƬ
% input_depmap = double(imread('data/moebius-depth.png'));
Color = double(imread('data/moebius-color.png'));
Color = imresize(Color,[1088,1376]);
G_gth = double(imread('data/mobius_big.png'));
G_ori = imresize(G_gth,1/2);
DepthGuide = imresize(G_ori,[1088,1376],'bicubic');
up_sample = zeros(size(G_gth,1),size(G_gth,2));
up_sample(1:2:end,1:2:end) = G_ori;
G_ori = up_sample(700:899,550:749);
% input_gth = imresize(input_gth,[1110,1390]);

% Data = data_process(input_depmap,input_colmap,input_gth); % �ü�����
% load ('ToFSimulated_2')
% Ϊ��ƥ��ά����������ϵ
% G_ori = DepthSample; %double(Data.d); % lr depth map
% G_ori = G_ori(500:699,500:699);

I = Color; % double(Data.c); % corresponded color map
I = I(700:899,550:749,:);

G_bcubic = DepthGuide; % double(Data.b); % bicubic ��ֵ���
G_bcubic = G_bcubic(700:899,550:749);

% G_gth = DepthGT; %double(Data.g); % �ü����gth2
% G_ori = imresize(G_gth,1/16);
% up_sample = zeros(size(G_gth,1),size(G_gth,2)); % load����ʱ���������ΪĬ���Ѿ��������29-33
% up_sample(1:16:end,1:16:end) = G_ori;
% G_ori = up_sample(500:699,500:699);

G_gth = G_gth(700:899,550:749);

%G_edg = Gabor_image(I); % �ü���ı߽���Ӧ
G_edg = gen_Gabor_maps(I,5,8,10,10);
imshow(G_edg,[]);
colormap
colormap jet
%G_edg = (G_edg - min(G_edg(:))) / (max(G_edg(:))-min(G_edg(:)));
%% ��ز����ĳ�ʼ��
win = 4; % �������ڵİ뾶��С 9 * 9
ser = 3; % ���򴰿ڵĴ�С 7 * 7
lamda_N = 0.0000008; % �������и��Ľ��0.006 ,��һ���̫��
lamda_s = 0.09; % ͬ��,0.008 , 5�Ļ��ܲ�
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
%subplot(1,3,2),imshow(uint8(I)),title('RGB');
imshow(uint8(hr_depth)),title('recovery');

%% ʵ��Segmentation [17]
% ratio = 0.5;
% kernelsize = 2;
% Iseg = vl_quickseg(I, ratio, kernelsize, maxdist)-->maxdist = [10,20];







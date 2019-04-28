function [x_h] = MRF_D(G_ori,G_bcubic,I,G_edg,Patch_size,lamda_N,lamda_s,win,ser)
% goal: Ax = b 
% author: Rui Wang
% date: 2019.03.12 
% referenced: an application of MRF

tic
[a,b] = size(G_ori);
index = zeros(a,b);
index(G_ori > 0) = 1;
M = spdiags(index(:),0,a*b,a*b);
B = G_ori(:);
% G_column = G_ori(:); %����������,ֱ�Ӱ�����ȥ��
% G_column_index = find(G_column~=0); % �ҵ�����������                                                                                                                                                                                                                                                                              
% G_column_valid = G_column(G_column_index); % ȡ����ֵ�õ�
% M =  sparse([],[],[],length(G_column_valid),length(G_column),0); % zeros(length(G_column_valid),length(G_column));����һ��ϡ��Ȩ�ؾ��󣬳�ֵΪ��ϡ�����
% for k = 1 : length(G_column_valid)
%   M(k,G_column_index(k)) = 1; % ������ֵ����ֵ�ò��֣�data term
% end
% M = sparse(M);
T_1 = M; % M' * M; % ��һ���ֵ-->D'M'MD

toc
disp('Building the Laplace Weighting for First term')

tic
YUV = double(rgb_2_yuv(I));%color image
toc
disp('Color space converted into YUV ')

I_pad = padarray(I, [win, win], 'symmetric');
YUV_pad = padarray(YUV, [ser, ser], 'symmetric');

tic
K = NLM_construct_11(I,I_pad, G_bcubic,Patch_size,win,ser); % 40000 * 40000
toc
disp('Building the Laplace weighting for NLM,a simplified version')

tic
W = Smoothness_build(G_bcubic,YUV_pad,G_edg,Patch_size,win,ser);
toc
disp('Builded NLM matrix');

% II = sparse(eye(Patch_size.^2)); % ��ʱ̫���޷����
% T_1 = M' * M;
T_2 = lamda_s * (W' * W); 
T_3 = lamda_N * (K * K');
x_h = (T_1 + T_2 + T_3)\B; % solve sparse linear system
% tic
% T = T_1 + T_2 + T_3;
% x_0 = B; % Ŀ��HR depthmap�ĳ�ʼֵ����ΪG��ֵ��ԭlr depmap ��ֵ���ֵ��
% error = 1.0e-5;
% tic
% x_h = CG_gradient(T,B,x_0,error); % ���ù����ݶȷ����Ŀ��
% toc
% disp('Optimizing using CG')

clear W
clear K
end


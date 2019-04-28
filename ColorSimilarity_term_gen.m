%% ColorSimilarity_term_gen.m
% This function generates the Color similarity term
% this work was originally proposed by J. Park on iccv2011
% "High Quality Depth Map Upsampling for 3D-ToF Cameras"
%
% Author: sgzthu @ Deptrum
% Last Modified: 2019.03.19
% 
function W_pq = ColorSimilarity_term_gen(im,para)
%%
[M N ~] = size(im);
Num = M*N;
im = double(im);
rgb2yuv_coef = [0.299 -0.14713 0.615;...
                0.587 -0.28886 -0.51499;...
                0.114 0.436 -0.10001];
yuv = zeros(size(im));
yuv(:,:,1) = rgb2yuv_coef(1,1).*im(:,:,1)+rgb2yuv_coef(2,1).*im(:,:,2)+rgb2yuv_coef(3,1).*im(:,:,3);
yuv(:,:,2) = rgb2yuv_coef(1,2).*im(:,:,1)+rgb2yuv_coef(2,2).*im(:,:,2)+rgb2yuv_coef(3,2).*im(:,:,3);
yuv(:,:,3) = rgb2yuv_coef(1,3).*im(:,:,1)+rgb2yuv_coef(2,3).*im(:,:,2)+rgb2yuv_coef(3,3).*im(:,:,3);

%% calculate the matrix W_pq
% first order neighbor
neighbor = [-1 1 0 0 0;...
            0  0 0 -1 1];    
sigma = para.sigma_Color;     
index_Matrix = reshape([1:Num],M,N);
W_pq = sparse(Num,Num);
for xy = neighbor
    x = xy(1); y = xy(2);
    % calculate the w_pq for all [p, q=p+(x,y)] pairs
    left = max(1+x,1);
    right = min(M+x,M);
    up = max(1+y,1);
    down = min(N+y,N);
    CD = zeros(size(yuv));
    CD(left-x:right-x,up-y:down-y,:) = yuv(left-x:right-x,up-y:down-y,:)-yuv(left:right,up:down,:);
    CD = sum(CD.^2,3);
    W_q = exp(-CD.^2./sigma^2);
    % calculate the indes of [p,q] pairs
    q_index = index_Matrix(left:right,up:down);
    p_index = index_Matrix(left-x:right-x,up-y:down-y);
    % assign the valid k_pq value to k_pq matrix
    W_q = W_q(left-x:right-x,up-y:down-y);
    W_pq_temp = sparse(p_index(:),q_index(:),W_q(:),Num,Num);
    W_pq = W_pq+W_pq_temp;
end

end
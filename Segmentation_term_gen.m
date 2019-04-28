%% Segmentation_term_gen.m
% This function generates the segmentation term
% this work was originally proposed by J. Park on iccv2011
% "High Quality Depth Map Upsampling for 3D-ToF Cameras"
%
% Author: sgzthu @ Deptrum
% Last Modified: 2019.03.19
% 
function W_pq = Segmentation_term_gen(seg,para)
%%
[M N ~] = size(seg);
Num = M*N;

%% calculate the matrix W_pq
% first order neighbor
neighbor = [-1 1 0 0 0;...
            0  0 0 -1 1]; 
index_Matrix = reshape([1:Num],M,N);
W_pq = sparse(Num,Num);
for xy = neighbor
    x = xy(1); y = xy(2);
    % calculate the w_pq for all [p, q=p+(x,y)] pairs
    left = max(1+x,1);
    right = min(M+x,M);
    up = max(1+y,1);
    down = min(N+y,N);
    W_q = zeros(M,N);
    W_q(left-x:right-x,up-y:down-y) = double(seg(left-x:right-x,up-y:down-y)==seg(left:right,up:down));
    W_q(left-x:right-x,up-y:down-y) = W_q(left-x:right-x,up-y:down-y)*(1-para.seg_penalty)+para.seg_penalty;
    % calculate the indes of [p,q] pairs
    q_index = index_Matrix(left:right,up:down);
    p_index = index_Matrix(left-x:right-x,up-y:down-y);
    % assign the valid k_pq value to k_pq matrix
    W_q = W_q(left-x:right-x,up-y:down-y);
    W_pq_temp = sparse(p_index(:),q_index(:),W_q(:),Num,Num);
    W_pq = W_pq+W_pq_temp;
end

end
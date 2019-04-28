%% NLS_term_gen.m
% This function generates the non local structure term
% this work was originally proposed by J. Park on iccv2011
% "High Quality Depth Map Upsampling for 3D-ToF Cameras"
%
% Author: sgzthu @ Deptrum
% Last Modified: 2019.03.19
% 
function K_pq = NLS_term_gen(im,para)
%%
% work on illuminance only
if size(im,3)>1
    im = rgb2ycbcr(im);
    im = im(:, :, 1);
end

[M N] = size(im);
Num = M*N;

% calculate matrix: Nabla_p*Nabla_p.'
Nabla_x = zeros(M,N);
Nabla_y = zeros(M,N);
Nabla_x(2:end,:) = im(2:end,:)-im(1:end-1,:);
Nabla_y(:,2:end) = im(:,2:end)-im(:,1:end-1);
Nabla_matrix = zeros([M,N,3]);
Nabla_matrix(:,:,1)=Nabla_x.*Nabla_x;       % 2*2对称矩阵，存三个元素就ok
Nabla_matrix(:,:,2)=Nabla_x.*Nabla_y;
Nabla_matrix(:,:,3)=Nabla_y.*Nabla_y;
% calculate Sigma_p_inverse
Sigma_matrix = zeros([M,N,3]);
Ones_matrix = ones([M,N,3]);
A_norm = zeros([M,N,3]);
Q = para.NLS_window_size;
for x = -Q:Q
    for y = -Q:Q
        % range of p':(left:right,up:down);
        % range of p :(left-x:right-x,up-y:down-y)
        left = max(1+x,1);
        right = min(M+x,M);
        up = max(1+y,1);
        down = min(N+y,N);
        Sigma_matrix(left-x:right-x,up-y:down-y,:)=Sigma_matrix(left-x:right-x,up-y:down-y,:)...
                                          +Nabla_matrix(left:right,up:down,:);
        A_norm(left-x:right-x,up-y:down-y,:)=A_norm(left-x:right-x,up-y:down-y,:)...
                                          +Ones_matrix(left:right,up:down,:);
    end
end
Sigma_matrix = Sigma_matrix./A_norm;

% %% 求逆；2*2对称矩阵求逆，也是对称2x2，存三个元素就ok
% Sigma_denominator = Sigma_matrix(:,:,1).*Sigma_matrix(:,:,3)-Sigma_matrix(:,:,2).^2;
% Sigma_inverse_matrix = zeros([M,N,3]);
% Sigma_inverse_matrix(:,:,1)=Sigma_matrix(:,:,3)./Sigma_denominator;
% Sigma_inverse_matrix(:,:,2)=-1*Sigma_matrix(:,:,2)./Sigma_denominator;
% Sigma_inverse_matrix(:,:,3)=Sigma_matrix(:,:,1)./Sigma_denominator;

% 不求逆
Sigma_inverse_matrix = Sigma_matrix;

% calculate the matrix pq
K_pq = sparse(Num,Num);
index_Matrix = reshape([1:Num],M,N);
sigma = para.sigma_NLS;
for x = -Q:Q
    for y = -Q:Q
        % calculate the k_pq for all [p, q=p+(x,y)] pairs
        left = max(1+x,1);
        right = min(M+x,M);
        up = max(1+y,1);
        down = min(N+y,N);
        ep = zeros(M,N);
        ep(left-x:right-x,up-y:down-y) = exp((2*Sigma_inverse_matrix(left-x:right-x,up-y:down-y,2).*x.*y...
             -Sigma_inverse_matrix(left-x:right-x,up-y:down-y,1).*x^2-Sigma_inverse_matrix(left-x:right-x,up-y:down-y,3).*y^2)./sigma^2);
        eq = zeros(M,N);
        eq(left:right,up:down) = exp((2*Sigma_inverse_matrix(left:right,up:down,2).*x.*y...
             -Sigma_inverse_matrix(left:right,up:down,1).*x^2-Sigma_inverse_matrix(left:right,up:down,3).*y^2)./sigma^2);
        K_q = 0.5*(ep+eq);
%         K_q(isnan(K_q)) = 0;      % 用逆的话需要判断奇异情况
        % calculate the indes of [p,q] pairs
        q_index = index_Matrix(left:right,up:down);
        p_index = index_Matrix(left-x:right-x,up-y:down-y);
        % assign the valid k_pq value to k_pq matrix
        K_q = K_q(left-x:right-x,up-y:down-y);
        K_pq_temp = sparse(p_index(:),q_index(:),K_q(:),Num,Num);
        K_pq = K_pq+K_pq_temp;
    end
end

% normalize each row (sum of coefficients for each p)
norm = 1./sum(K_pq.');  
[indx,indy,value] = find(K_pq); % 找出所有非零值坐标
norm = sparse(indx,indy,norm(indx),Num,Num); 
K_pq = K_pq .* norm;

% 小于阈值的抛弃掉，提高稀疏程度
[indx,indy,value] = find(K_pq);
value = value.*(value>para.NLS_th);           
K_pq = sparse(indx,indy,value,Num,Num);  

end
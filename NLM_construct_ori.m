function K = NLM_construct_ori(I,I_pad, G_bicubic, Patch_size,win,ser)

% non-local term;
sigma_r = 0.02; % 3
w = win; % Search window radius
sigma_nlm = 0.5; % 0.02

col_map = I_pad;
K = zeros(Patch_size^2); % 4e4 \times 4e4

% tic
% [I_grad_x,I_grad_y] = gradient(I); % ��colorͼ��ÿ�����ص���ݶ�ֵ
% ti_d = get_grad_cell(I_grad_x,I_grad_y,a,b);
% toc
% disp('Computing gradient')

% tic
% [zuo_b_c,zuo_b_l] = get_coordinate(a,b);
% toc
% disp('Geting the coordinate')
tic
[y,x] = meshgrid(1:200,1:200);
toc
disp('generate the coordinate')

for i = 1 : Patch_size
      for j = 1 : Patch_size
          
        NLM_win = zeros(2 * w + 1); % �������ڵĴ�С
        H = zeros(Patch_size + 2 * w); % 68 * 68����һ���������ڴ����Ͻǵ�һ���㿪ʼ����СΪ9 * 9
        Gau_kernel = zeros(2*w+1); % 9 * 9Ϊ������Χ������������ȡ�ã�����Ϊ5 * 5
        j_min = max(j-w,1); % initial-->1
        j_max = min(j+w,Patch_size); % initial-->5
        i_min = max(i-w,1); % 1
        i_max = min(i+w,Patch_size); % 5
        tic
        
        zegma = G_bicubic(j_min:j_max,i_min:i_max); % ��������ڵ��ݶȾ���ĺͣ���ʽ��8��
%         t = zeros(2,2,size(zegma,1)*size(zegma,2));
%         t = cat(3,zegma{:}); % ��cellת����3ά����,����ת��
%         t_sum = double(sum(t,3) ./ (size(zegma,1)*size(zegma,2))); % ���յ�����ά�����
%         t_sum_3d = repmat(t_sum,[1,1,size(zegma,1)*size(zegma,2)]); % 2dת3d
%         toc
%         disp('Convert the zegma')
%         tic
%         zuo_b_patch_c = zuo_b_c(j_min:j_max,i_min:i_max);
%         zuo_b_3d_c = zeros(1,2,size(zegma,1)*size(zegma,2));
%         zuo_b_3d_c = cat(3,zuo_b_patch_c{:});
%         toc
%         disp('Convert the x-axis coordiante')
%         
%         tic
%         zuo_b_patch_l = zuo_b_l(j_min:j_max,i_min:i_max);
%         zuo_b_3d_l = zeros(1,2,size(zegma,1)*size(zegma,2));
%         zuo_b_3d_l = cat(3,zuo_b_patch_l{:});
%         toc
%         disp('Convert the y-axis coordiante')
%         
        tic
        diff_R = abs(col_map(j:j+2*w,i:i+2*w,1)-col_map(j+w,i+w,1)); 
        diff_G = abs(col_map(j:j+2*w,i:i+2*w,2)-col_map(j+w,i+w,2));
        diff_B = abs(col_map(j:j+2*w,i:i+2*w,3)-col_map(j+w,i+w,3));
        color_term = exp(-(diff_R.^2+diff_G.^2+diff_B.^2) / (2 * sigma_nlm^2));
        toc
        disp('computing the color term')
%         
%         tic
%         k_pq = zeros(1,1,size(zegma,1) * size(zegma,2));
%         for k = 1 : size(zegma,1) * size(zegma,2)
%             k_pq(:,:,k) = (zuo_b_3d_c(:,:,k)-zuo_b_c{j ,i}) * (zuo_b_3d_l(:,:,k)-zuo_b_l{j,i}); % * t_sum_3d(:,:,k) * 
%         end
%         k_pq_matrix = reshape(k_pq,j_max-j_min+1,i_max-i_min+1); % ��ص�ǰ�����򴰿ڴ�С,���˱߽�Ϊ5 * 5�⣬���඼��9 * 9
%         toc
        tic
        diff_x = abs(x(j_min:j_max,i_min:i_max)-x(j,i)); 
        diff_y = abs(y(j_min:j_max,i_min:i_max)-y(j,i));
        toc
        disp('computing the spatial distance term')
%       k_pq_matrix_c = bsxfun(@minus,k_pq_matrix,mean(k_pq_matrix));
%       k_pq_matrix = bsxfun(@rdivide,k_pq_matrix_c,std(k_pq_matrix));
        Gau_kernel((j_min:j_max)-j+w+1,(i_min:i_max)-i+w+1) = exp(-(diff_x.^2 + diff_y.^2)/(2 * sigma_r^2));
        bi_filter = Gau_kernel .* color_term; 
        NLM_win = bi_filter;
        NLM_win(w+1,w+1) = 1; % �������ֵΪ0,���Կ���Gau_kernelֻ�к���5 * 5�Ĳ�����ֵ
        H(j:j+2*w,i:i+2*w) = NLM_win; % ��Ϊ HΪ68 * 68��
        mask = H(w+1:end-w,w+1:end-w); % ȡ��һ��block
        K((i-1) * Patch_size + j,:) = mask(:)'; % ������ˣ������ھ��
      end
end
K = sparse(K);
end


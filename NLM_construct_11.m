function K = NLM_construct_11(I,I_pad,G_bicubic, Patch_size,win,ser)

% non-local term;
[a,b] = size(G_bicubic);
w = win; % Search window radius
f = ser; % \mathcal{A} = 5 * 5,��������Ϊ202 * 202
sigma_nlm = 0.05;

K = zeros(Patch_size^2); % 4e4 \times 4e4
% K = sparse(Patch_size^2,Patch_size^2); % ����һ��ϡ��Ȩ�ؾ��󣬳�ֵΪ��ϡ�����,��Ȼ�Ļ��ڴ�᲻��
% K = [];
tic
[I_grad_x,I_grad_y] = gradient(I); % ��colorͼ��ÿ�����ص���ݶ�ֵ
ti_d = get_grad_cell(I_grad_x,I_grad_y,a,b);
toc
disp('Computing gradient')

tic
[zuo_b_c,zuo_b_l] = get_coordinate(a,b);
toc
disp('Geting the coordinate')

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
        zegma = ti_d(j_min:j_max,i_min:i_max); % ��������ڵ��ݶȾ���ĺͣ���ʽ��8��
        t = zeros(2,2,size(zegma,1)*size(zegma,2));
        t = cat(3,zegma{:}); % ��cellת����3ά����,����ת��
        t_sum = double(sum(t,3) ./ (size(zegma,1)*size(zegma,2))); % ���յ�����ά�����
        t_sum_3d = repmat(t_sum,[1,1,size(zegma,1)*size(zegma,2)]); % 2dת3d
        toc
        disp('Convert the zegma')
        tic
        zuo_b_patch_c = zuo_b_c(j_min:j_max,i_min:i_max);
        zuo_b_3d_c = zeros(1,2,size(zegma,1)*size(zegma,2));
        zuo_b_3d_c = cat(3,zuo_b_patch_c{:});
        toc
        disp('Convert the x-axis coordiante')
        tic
        zuo_b_patch_l = zuo_b_l(j_min:j_max,i_min:i_max);
        zuo_b_3d_l = zeros(1,2,size(zegma,1)*size(zegma,2));
        zuo_b_3d_l = cat(3,zuo_b_patch_l{:});
        toc
        disp('Convert the y-axis coordiante')
        k_pq = zeros(1,1,size(zegma,1) * size(zegma,2));
        for k = 1 : size(zegma,1) * size(zegma,2)
            k_pq(:,:,k) = double((zuo_b_3d_c(:,:,k)-zuo_b_c{j,i}) * t_sum_3d(:,:,k) * (zuo_b_3d_l(:,:,k)-zuo_b_l{j,i}));
        end
        k_pq_matrix = reshape(k_pq,j_max-j_min+1,i_max-i_min+1); % ��ص�ǰ�����򴰿ڴ�С,���˱߽�Ϊ5 * 5�⣬���඼��9 * 9
        Gau_kernel((j_min:j_max)-j+w+1,(i_min:i_max)-i+w+1) = exp(-k_pq_matrix);
        NLM_win = Gau_kernel;
        NLM_win(w+1,w+1) = 1; % �������ֵΪ0,���Կ���Gau_kernelֻ�к���5 * 5�Ĳ�����ֵ
        H(j:j+2*w,i:i+2*w) = NLM_win; % ��ΪHΪ68 * 68��
        mask = H(w+1:end-w,w+1:end-w); % ȡ��һ��block
        K((i-1) * Patch_size + j,:) = mask(:); % ������ˣ������ھ��
        
    end
end
K = sparse(K);
end


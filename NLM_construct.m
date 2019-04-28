 
function K = NLM_construct(I,G_bicubic, Patch_size,win,ser)

% non-local term;
col_map = I;
[a,b] = size(G_bicubic);
w = win; % Search window radius
f = ser; % \mathcal{A} = 5 * 5,��������Ϊ202 * 202
sigma_nlm = 0.02;

K = zeros(Patch_size^2); % 4e4 \times 4e4
% K = sparse([],[],[],Patch_size^2,Patch_size^2,0); % ����һ��ϡ��Ȩ�ؾ��󣬳�ֵΪ��ϡ�����,��Ȼ�Ļ��ڴ�᲻��                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             

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
        Gau_kernel = zeros(2 * w+1); % 9 * 9Ϊ������Χ������������ȡ�ã�����Ϊ5 * 5
        j_min = max(j-w,1); % initial-->1
        j_max = min(j+w,Patch_size); % initial-->5
        i_min = max(i-w,1); % 1
        i_max = min(i+w,Patch_size); % 5
                                                         
        diff_R = abs(col_map(j_min:j_max,i_min:i_max,1)-col_map(j,i,1)); 
        diff_G = abs(col_map(j_min:j_max,i_min:i_max,2)-col_map(j,i,2));
        diff_B = abs(col_map(j_min:j_max,i_min:i_max,3)-col_map(j,i,3));
        Gau_kernel((j_min:j_max)-j+w+1,(i_min:i_max)-i+w+1) = ...
        exp(-(diff_R.^2+diff_G.^2+diff_B.^2) / (2 * sigma_nlm^2));
        NLM_win = Gau_kernel;
        NLM_win = NLM_win;
        NLM_win(w+1,w+1) = 1; % �������ֵΪ0,���Կ���Gau_kernelֻ�к���5 * 5�Ĳ�����ֵ
        H(j:j+2*w,i:i+2*w) = NLM_win; % ��ΪHΪ68 * 68��
        mask = H(w+1:end-w,w+1:end-w); % ȡ��һ��block
        K((i-1) * Patch_size + j,:) = mask(:); % ������ˣ������ھ��
    end
end
K = sparse(K); %spdiags(K,0,a*b,a*b);
end







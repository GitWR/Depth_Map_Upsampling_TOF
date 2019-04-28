 
function K = NLM_construct(I,G_bicubic, Patch_size,win,ser)

% non-local term;
col_map = I;
[a,b] = size(G_bicubic);
w = win; % Search window radius
f = ser; % \mathcal{A} = 5 * 5,下面扩充为202 * 202
sigma_nlm = 0.02;

K = zeros(Patch_size^2); % 4e4 \times 4e4
% K = sparse([],[],[],Patch_size^2,Patch_size^2,0); % 构造一个稀疏权重矩阵，初值为零稀疏矩阵,不然的话内存会不足                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             

tic
[I_grad_x,I_grad_y] = gradient(I); % 求color图中每个像素点的梯度值
ti_d = get_grad_cell(I_grad_x,I_grad_y,a,b);
toc
disp('Computing gradient')

tic
[zuo_b_c,zuo_b_l] = get_coordinate(a,b);
toc
disp('Geting the coordinate')

for i = 1 : Patch_size
    for j = 1 : Patch_size
        NLM_win = zeros(2 * w + 1); % 搜索窗口的大小
        H = zeros(Patch_size + 2 * w); % 68 * 68，第一个搜索窗口从左上角第一个点开始，大小为9 * 9
        Gau_kernel = zeros(2 * w+1); % 9 * 9为搜索范围，领域在里面取得，领域为5 * 5
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
        NLM_win(w+1,w+1) = 1; % 和自身的值为0,可以看到Gau_kernel只有后面5 * 5的部分有值
        H(j:j+2*w,i:i+2*w) = NLM_win; % 因为H为68 * 68的
        mask = H(w+1:end-w,w+1:end-w); % 取出一个block
        K((i-1) * Patch_size + j,:) = mask(:); % 交错相乘，类似于卷积
    end
end
K = sparse(K); %spdiags(K,0,a*b,a*b);
end







 
function W = Smoothness_build(G_bcubic,YUV_pad,G_edge,Patch_size,win,ser)

% smoothness term
w = win; % Search window radius
f = ser;
W = zeros(Patch_size^2);
% W = sparse(Patch_size^2,Patch_size^2); % ����һ��ϡ��Ȩ�ؾ��󣬳�ֵΪ��ϡ�����?,��Ȼ�Ļ��ڴ�᲻��?
% W = [];
col_map = YUV_pad; %padarray(YUV, [f,f], 'symmetric');  % ��Ե���ۣ����������������?
sigma_I  = 0.2; % YUV weighted  0.5
sigma_g = 4; % depth weighted
sigma_p = 11; %color weighted

for i = 1 : Patch_size
    for j = 1 : Patch_size 
        
        w_pq = zeros(2 * w + 1); % �������ڵĴ�С 7 * 7
        H = zeros(Patch_size + 2 * w); % 66 * 66����һ���������ڴ����Ͻǵ�һ���㿪ʼ����СΪ9 * 9
        w_d = zeros(2 * w + 1); % 9 * 9Ϊ������?
        w_e = zeros(2 * w + 1);
        w_s = zeros(2 * w + 1);
        
        j_min = max(j - w,1); % initial-->1
        j_max = min(j + w,Patch_size); % initial-->5
        i_min = max(i-w,1); % 1
        i_max = min(i + w,Patch_size); % 5
        
        P_x = col_map(j:j + 2 * f,i:i + 2 * f,:); % ��xΪ���ĵ�һ����2f+1��ֱ����С������
        
        tic
        diff_R = abs(col_map(j:j + 2 * f,i:i + 2 * f,1)-col_map(j + f,i + f,1)); 
        diff_G = abs(col_map(j:j + 2 * f,i:i + 2 * f,2)-col_map(j + f,i + f,2));
        diff_B = abs(col_map(j:j + 2 * f,i:i + 2 * f,3)-col_map(j + f,i + f,3));
%       diff_R = abs(col_map(j_min:j_max,i_min:i_max,1)-col_map(j,i,1)); 
%       diff_G = abs(col_map(j_min:j_max,i_min:i_max,2)-col_map(j,i,2));
%       diff_B = abs(col_map(j_min:j_max,i_min:i_max,3)-col_map(j,i,3));
        w_c =exp(-(diff_R.^2 + diff_G.^2 + diff_B.^2) / (2 * sigma_I^2)); % ��ʽ��9��
        toc
        disp('Compute the YUV weighting term')
        
        tic
        for k = j_min + f:j_max+f
            for l = i_min + f:i_max + f
                if ~(k == (j + f) && l == (i + f)) % λ�ò�Ҫ�غ�
                    P_y = col_map(k - f:k + f,l - f:l + f,:);
                    diff = sum((P_x - P_y),3); 
                    k_pq = exp(-sum(sum (w_c .* (diff.^2)))/(2 * sigma_p^2));
                    w_pq(k - f - j + w + 1,l - f - i + w + 1) = k_pq;
                end
            end
        end
        toc
        disp('inner color term')
        
        tic
        w_d((j_min:j_max) - j + w + 1,(i_min:i_max) - i + w + 1)=...
        exp(-(abs(G_bcubic(j_min:j_max,i_min:i_max)-G_bcubic(j,i))).^2 / (2 * sigma_g^2)); % ��ʽ12
        toc
        disp('Compute the Depth weighting term')
        
        tic
        w_e((j_min:j_max) - j + w + 1,(i_min:i_max) - i + w + 1) = 1 ./ (sqrt(G_edge(j,i)^2 + G_edge(j_min:j_max,i_min:i_max).^2) + 1);
        toc
        disp('Compute the edge weighting term')
        
        w_pq = w_pq .* w_d;
        w_pq = -w_pq / sum(w_pq(:)); % normalization
        w_pq(w+1,w+1) = 1; % ��������ֵΪ0,���Կ���Gau_kernelֻ�к���5 * 5�Ĳ�����ֵ
        H(j:j + 2 * w,i:i + 2 * w) = w_pq; % ��ΪHΪ68 * 68��
        mask = H(w + 1:end - w,w + 1 :end - w); % ȡ��һ��patch
        W((i-1) * Patch_size + j,:) = mask(:); % ������ˣ������ھ���?
    end
end
W = sparse(W);
end








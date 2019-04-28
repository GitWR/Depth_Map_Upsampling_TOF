 
function NLM = NLM_construct(colormap_temp,guidedmap_temp, block_size,halfwindow_size,halfsearch_size)

%non-local term
w=halfwindow_size;%half size of the search window
f=halfsearch_size;%block size to diff, Neighbor = 7 *7,下面扩充为66 * 66

sigma1=4; %depth weighted
sigma_c=0.35;
sigma2=sigma_c*1.732*11; %color weighted
% sigma3=3; %gausian weighted 
sigma4=0.1;%bi weighted
% [Dx,Dy]= gradient(guidedmap_temp);
% D = sqrt(Dx.^2+Dy.^2);
% norm_D = (D-min(D(:)))/(max(D(:))-min(D(:)));
% sigma1 = 3*exp(-log(6)*norm_D);
NLM=zeros(block_size^2); %3600 * 3600
image = padarray(colormap_temp, [f,f], 'symmetric');  % 边缘对折
% [X,Y] = meshgrid(-f:f,-f:f);
%G_kernel = exp(-(X.^2+Y.^2)/(2*sigma3^2));%gaussian kernal

       for t=1:block_size
            for s=1:block_size
                
                 NLM_temp=zeros(2*w+1); % 搜索窗口的大小 
                 H=zeros(block_size+2*w);% 68 * 68，第一个搜索域从左上角第一个点开始，大小为9 * 9
                 guided_kernel=zeros(2*w+1); % 9 * 9为搜索范围，领域在里面取得，领域为5 * 5
                 sMin = max(s-w,1); % initial 1
                 sMax = min(s+w,block_size); % initial 5
                 tMin = max(t-w,1); % 1
                 tMax = min(t+w,block_size); % 5
                 sigma_P=image(s:s+2*f,t:t+2*f,:); % patch 大小 7 * 7 （文中的w * w）
                 guided_kernel((sMin:sMax)-s+w+1,(tMin:tMax)-t+w+1)=...
                 exp(-(abs(guidedmap_temp(sMin:sMax,tMin:tMax)-guidedmap_temp(s,t))).^2/(2*sigma1^2)); % a_{xy}^{\hat{D}}.(7)
                 dR = abs(sigma_P(:,:,1)-image(s+f,t+f,1)); % 邻域窗口的中心点的值 R
                 dG = abs(sigma_P(:,:,2)-image(s+f,t+f,2));
                 dB = abs(sigma_P(:,:,3)-image(s+f,t+f,3));
                 color_kernel = exp(-(dR.^2+dG.^2+dB.^2)/(2*3*sigma4^2)); % 公式（9）中的第二部分
                 for k=sMin+f:sMax+f
                    for l=tMin+f:tMax+f
                        if ~(k==(s+f) && l==(t+f)) % w * w 的中心点，也在（x,x）处
                          sigma_Q=image(k-f:k+f,l-f:l+f,:); % 与RGB图像的映射关系
                          diff=sum((sigma_P-sigma_Q),3); % 差异
                          kpq=exp(-sum(sum (color_kernel.*(diff.^2)))/(2*3*sigma2^2)); % 公式8，因为他们都是在w * w的patch里面计算得到的                             
                          NLM_temp(k-f-s+w+1,l-f-t+w+1)=-kpq; % 公式8 a_{xy}^{\hat{D}}
                        end
                     end
                 end
                  NLM_temp=NLM_temp.*guided_kernel; % 可以看到guided_kernel只有后面5 * 5的部分有值
                  NLM_temp=-NLM_temp/sum(NLM_temp(:)); % 一个搜索窗口内的结果 公式6
                  NLM_temp(w+1,w+1)=1; % 和自身的值为0
                  H(s:s+2*w,t:t+2*w)=NLM_temp; % 因为H为68 * 68的
                  mask=H(w+1:end-w,w+1:end-w); % 取出一个block
                  NLM((t-1)*block_size+s,:)=mask(:); % 交错相乘，类似于卷积
             end
       end
end








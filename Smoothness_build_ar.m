 
function NLM = NLM_construct(colormap_temp,guidedmap_temp, block_size,halfwindow_size,halfsearch_size)

%non-local term
w=halfwindow_size;%half size of the search window
f=halfsearch_size;%block size to diff, Neighbor = 7 *7,��������Ϊ66 * 66

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
image = padarray(colormap_temp, [f,f], 'symmetric');  % ��Ե����
% [X,Y] = meshgrid(-f:f,-f:f);
%G_kernel = exp(-(X.^2+Y.^2)/(2*sigma3^2));%gaussian kernal

       for t=1:block_size
            for s=1:block_size
                
                 NLM_temp=zeros(2*w+1); % �������ڵĴ�С 
                 H=zeros(block_size+2*w);% 68 * 68����һ������������Ͻǵ�һ���㿪ʼ����СΪ9 * 9
                 guided_kernel=zeros(2*w+1); % 9 * 9Ϊ������Χ������������ȡ�ã�����Ϊ5 * 5
                 sMin = max(s-w,1); % initial 1
                 sMax = min(s+w,block_size); % initial 5
                 tMin = max(t-w,1); % 1
                 tMax = min(t+w,block_size); % 5
                 sigma_P=image(s:s+2*f,t:t+2*f,:); % patch ��С 7 * 7 �����е�w * w��
                 guided_kernel((sMin:sMax)-s+w+1,(tMin:tMax)-t+w+1)=...
                 exp(-(abs(guidedmap_temp(sMin:sMax,tMin:tMax)-guidedmap_temp(s,t))).^2/(2*sigma1^2)); % a_{xy}^{\hat{D}}.(7)
                 dR = abs(sigma_P(:,:,1)-image(s+f,t+f,1)); % ���򴰿ڵ����ĵ��ֵ R
                 dG = abs(sigma_P(:,:,2)-image(s+f,t+f,2));
                 dB = abs(sigma_P(:,:,3)-image(s+f,t+f,3));
                 color_kernel = exp(-(dR.^2+dG.^2+dB.^2)/(2*3*sigma4^2)); % ��ʽ��9���еĵڶ�����
                 for k=sMin+f:sMax+f
                    for l=tMin+f:tMax+f
                        if ~(k==(s+f) && l==(t+f)) % w * w �����ĵ㣬Ҳ�ڣ�x,x����
                          sigma_Q=image(k-f:k+f,l-f:l+f,:); % ��RGBͼ���ӳ���ϵ
                          diff=sum((sigma_P-sigma_Q),3); % ����
                          kpq=exp(-sum(sum (color_kernel.*(diff.^2)))/(2*3*sigma2^2)); % ��ʽ8����Ϊ���Ƕ�����w * w��patch�������õ���                             
                          NLM_temp(k-f-s+w+1,l-f-t+w+1)=-kpq; % ��ʽ8 a_{xy}^{\hat{D}}
                        end
                     end
                 end
                  NLM_temp=NLM_temp.*guided_kernel; % ���Կ���guided_kernelֻ�к���5 * 5�Ĳ�����ֵ
                  NLM_temp=-NLM_temp/sum(NLM_temp(:)); % һ�����������ڵĽ�� ��ʽ6
                  NLM_temp(w+1,w+1)=1; % �������ֵΪ0
                  H(s:s+2*w,t:t+2*w)=NLM_temp; % ��ΪHΪ68 * 68��
                  mask=H(w+1:end-w,w+1:end-w); % ȡ��һ��block
                  NLM((t-1)*block_size+s,:)=mask(:); % ������ˣ������ھ��
             end
       end
end








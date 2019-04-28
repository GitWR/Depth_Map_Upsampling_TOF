function featureVector = gen_Gabor_maps(color,u,v,m,n)


% Inputs:
%       u	:	No. of scales (usually set to 5)
%       v	:	No. of orientations (usually set to 8)
%       m	:	No. of rows in a 2-D Gabor filter (an odd integer number usually set to 39)
%       n	:	No. of columns in a 2-D Gabor filter (an odd integer number usually set to 39)
%       m and n is actually the size of each filter
% Output:
%       gaborArray: A u by v cell arrray, element of which are m by n
%       matries; each matrix being a 2-D Gabor filter

img = double(rgb2gray(uint8(color)));

gaborArray = cell(u,v);
fmax = 0.25;
gama = sqrt(2);
eta = sqrt(2);

for i = 1:u
    
    fu = fmax/((sqrt(2))^(i-1));
    alpha = fu/gama;
    beta = fu/eta;
    
    % 尺度和方向的选择，通过角度
    for j = 1:v
        tetav = ((j-1)/8)*pi;
        gFilter = zeros(m,n);  
        for x = 1:m
            for y = 1:n
                xprime = (x-((m+1)/2))*cos(tetav)+(y-((n+1)/2))*sin(tetav);
                yprime = -(x-((m+1)/2))*sin(tetav)+(y-((n+1)/2))*cos(tetav);
                gFilter(x,y) = (fu^2/(pi*gama*eta))*exp(-((alpha^2)*(xprime^2)+(beta^2)*(yprime^2)))*exp(1i*2*pi*fu*xprime);
            end
        end
        gaborArray{i,j} = gFilter;
    end
end

gaborResult = cell(u,v);
for i = 1:u
    for j = 1:v
        gaborResult{i,j} = conv2(img,gaborArray{i,j},'same'); % 利用得到的多尺度多方向的卷积模板对图像进行卷积操作
        % J{u,v} = filter2(G{u,v},I);
    end
end

[n1,m1] = size(img);
%featureVector = zeros(n1,m1,u*v);
featureVector = zeros(n1,m1);
c = 0;

for i = 1:u
    for j = 1:v
        c = c+1;
        gaborAbs = abs(gaborResult{i,j}); 
        figure;
        imshow(gaborAbs);
        featureVector = featureVector + gaborAbs; % 得到uv个和原图大小一样的feature maps
    end
end

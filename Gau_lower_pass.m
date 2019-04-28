clc
clear all
close all
IA = double(imread('art-depth.png'));
[f1,f2] = freqspace(size(IA),'meshgrid');
D = 100 / size(IA,1);
r = f1.^2 + f2.^2;
Hd = ones(size(IA));

for i = 1 : size(IA,1)
    for j = 1 : size(IA,2)
        t = r(i,j) / (D * D);
        Hd(i,j) = exp(-t);
    end
end

Y = fft2(double(IA));
Y = fftshift(Y);
Ya = Y .* Hd;
Ya = ifftshift(Ya);
Ia = real(ifft2(Ya));

imshow(uint8(IA)); title('Ô­Í¼Ïñ')
figure
imshow(uint8(Ia)); title('µÍÍ¨Í¼Ïñ')
figure
imshow(uint8(IA - Ia),[]); title('²ÐÓàÍ¼')

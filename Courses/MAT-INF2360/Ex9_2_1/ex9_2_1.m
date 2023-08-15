clear all;
close all;
clc;

X = double(imread('lena.png','png'));

newvals = X(:,:,1) + X(:,:,2) + X(:,:,3);
Z = newvals/max(max(newvals))*255;
[m n]= size(Z);

for i = 1:m
    for j = 1:n
        if (Z(i,j) < 127.5)
            Z(i,j) = 0.;
        elseif (Z(i,j) >= 127.5)
            Z(i,j) = 255.;
        end
    end
end

imshow(uint8(Z));
%imwrite(uint8(Z),'lena_blackwhite.png','png');
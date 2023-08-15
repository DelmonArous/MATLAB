clear all;
close all;
clc;

X = double(imread('lena.png','png'));

for epsilon = [0.1 0.01 0.001]
    Z = contrastadjust(X,epsilon);
    figure();
    imshow(uint8(Z));
end
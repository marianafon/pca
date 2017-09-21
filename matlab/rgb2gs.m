clear; clc;
% lê a imagem colorida 
A = imread('../source/2000.jpg');
% passa pra grayscale
B = rgb2gray(A);
% exporta a imagem em grayscale
imwrite(B,'../source/2000G.jpg');

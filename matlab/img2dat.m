clear; clc;
% lê a imagem colorida 
A = imread('../source/2000.jpg');
% passa pra grayscale
B = rgb2gray(A);
% exporta a matriz B para um csv
csvwrite('../datasets/2000.dat', B);